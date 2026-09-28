"""Same-native finite hydrogen inventory and reciprocal bound-free exchange.

Counterexample candidate: thermodynamic reverse rates, finite native levels,
and a specified diluted Planck bath. This is not a complete plasma opacity.
"""
from pathlib import Path
import ctypes
import json
import signal
import sys
import time
import numpy as np
from scipy.interpolate import RectBivariateSpline
import def_native_inventory_flow as old

OUT=old.OUT.parent/'def-native-hydrogen-exchange'
C=old.C
K=1.380649e-16
H=1.9864458571489287e-16/C
NA=6.02214076e23
write=old.write
ATOMIC=Path('/home/lpaiu/work/direct-eos-gr33/photon-hhe-coupled/atomic.so')
LEVELS=Path('/home/lpaiu/work/direct-eos-gr33/photon-eos-levels-repaired/levels.so')


class Native(old.Native):
    def __init__(self,cap=4000):
        super().__init__(cap);self.y0=float(self.target[0,0]/self.target[0,:2].sum())
        self.epsH=float(self.target[0,:2].sum());self.nH=self.epsH*self.fan.cx*NA
        array=np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS')
        self.levels=ctypes.CDLL(str(LEVELS)).photon_levels
        self.levels.argtypes=[ctypes.c_int,ctypes.c_double,array,array];self.levels.restype=None
        self.cross=ctypes.CDLL(str(ATOMIC)).hydrogen_cross
        self.cross.argtypes=[ctypes.c_int]*3+[ctypes.c_double,array,array,array];self.cross.restype=None
        const=np.loadtxt(old.OUT.parent/'def-photon-eos-populations/levels-repaired/constants.txt')
        self.binding=109737.3156816/(1+const[2]/const[3])/np.arange(1,11)**2*1.9864458571489287e-16
        self.binding[0]=self.energy[0]*K
        self.oldground=float(np.loadtxt(old.OUT.parent/'def-photon-atomic-rates/fort.93')[0,3])

    def state(self,x,lt,y):
        assert 0<y<1
        d=self.prefix;j=int(np.argmin(abs(np.log(d['T'])-lt)));last=np.log(d['T'][j])
        dr=x-d['log_density_ratio'][j];dt=lt-last;inv=np.exp(-lt)-np.exp(-last)
        fields=d['fields'][j].copy();fields[:316]+=np.where(self.active,self.energy*inv+self.charges*(dr-1.5*dt),0)
        fields[316]+=-self.diss*inv-dr+1.5*dt;fields[317]+=self.molion*inv+dr-1.5*dt
        fields[0]+=np.log((1-y)/y)-np.log((1-self.y0)/self.y0)
        target=self.target.copy();target[0,:2]=self.epsH*np.array([y,1-y])
        a,_,err=self.ion.constrain(self.lr+x,lt,target,fields,target_molecules=self.base['molecular_H_fractions'],tolerance=1e-12)
        for _ in range(5):
            hy=float(a['number_fractions'][0,0]/a['number_fractions'][0,:2].sum())
            if abs(hy/y-1)<1e-8:break
            self.ion.fields[0]+=np.log(hy/y)+np.log((1-y)/(1-hy))
            a,_,err=self.ion.constrain(self.lr+x,lt,target,self.ion.fields.copy(),target_molecules=self.base['molecular_H_fractions'],tolerance=1e-12)
        else:raise AssertionError(('Relative neutral H constraint',x,lt,y,hy))
        lib=self.ion.gas.gas_lib
        def arr(name,n,dtype=ctypes.c_double):
            return np.ctypeslib.as_array((dtype*n).in_dll(lib,'__mod_excitation_block_MOD_'+name)).copy()
        count=int(arr('extrace_count',1,ctypes.c_int)[0]);ids=arr('extrace_ids',636,ctypes.c_int).reshape(318,2)[:count]
        row=np.flatnonzero(np.all(ids==[1,0],axis=1));assert len(row)==1
        L=float((arr('extrace_value',318*6).reshape(318,6)[:count,3]/arr('extrace_scale',318)[:count])[row[0]])
        terms=np.zeros(30);self.levels(1,float(lt),arr('x',5),terms);tail=terms.reshape(3,10)[0]
        terms=tail-np.r_[tail[1:],0.];assert np.all(terms>=0) and L>=0
        fraction=np.zeros(10);fraction[0]=np.exp(-L)
        if L>0:fraction[1:]=-np.expm1(-L)*terms[1:]/sum(terms[1:])
        assert abs(sum(fraction)-1)<1e-13
        return dict(raw=a['eos'],fraction=fraction,affinity=float(self.ion.fields[0]),population_error=err,x=x,lt=lt,y=y,L=L)

    def rates(self,state,Trad,order=16):
        """Number/s and erg/s per atomic H nucleus, separated by photon term.

        A, Rsp, Rstim use the *current* EOS affinity, not a frozen LTE ratio.
        The inverse is a thermodynamic closure, not a microscopic Milne proof.
        """
        T=np.exp(state['lt']);gx,gw=np.polynomial.legendre.leggauss(order)
        edges=np.array([0.,.1,1.,4.,12.,30.,60.]);du=np.diff(edges)/2
        z=((edges[:-1]+edges[1:])[:,None]/2+du[:,None]*gx).ravel();w=(du[:,None]*gw).ravel()
        result=np.zeros((3,2))
        for n,f in enumerate(state['fraction'],1):
            if f==0:continue
            b=self.binding[n-1];logpop=np.log(state['y'])+np.log(f)
            for which,temp in [('abs',Trad),('emit',T)]:
                energy=b+K*temp*z;nu=energy/H
                native_nu=np.ascontiguousarray((self.oldground/n**2+K*temp*z)/float(np.float32(6.6256e-27)))
                bf=np.zeros_like(nu);ff=np.zeros_like(nu);self.cross(len(nu),1,n,float(T),native_nu,bf,ff)
                assert np.all(bf>=0)
                modes=8*np.pi*nu**2/C**2*bf*(K*temp/H)*w
                occupation=np.exp(-energy/(K*Trad))/(-np.expm1(-energy/(K*Trad)))
                if which=='abs':a=np.exp(logpop)*occupation*modes;result[0]+=[a.sum(),a@energy]
                else:
                    reverse=np.exp(logpop+state['affinity']-b/(K*T)-z)*modes
                    result[1]+=[reverse.sum(),reverse@energy]
                    stim=reverse*occupation;result[2]+=[stim.sum(),stim@energy]
        assert np.all(np.isfinite(result)) and np.all(result>0)
        return result


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='bdbbfe161',
        claim='Evolve the retained H bound-free reaction, its native chemical energy and its opposite photon energy/momentum transfer in the actual conservative spherical outflow.',
        decision='Compare reactive and frozen-inventory histories and their direct retarded charge; stop if source consistency, native support or runtime fails.',
        model='Original native EOS and finite H levels; all other elemental ionic inventories frozen to the surface. Translated SYNspec H cross sections; reciprocal reverse rates use current native ionic affinity. Specified outgoing diluted Planck bath matched to saved bolometric luminosity, initially present. This is a thermodynamic inverse-rate model, not a microscopic Milne derivation or full photon transport.',
        budget=dict(native_calls=4000,native_seconds=120,flow_seconds=180,readout_seconds=45,CPU_threads=1,memory_GB=2,pilot_cells=224,production_cells=[448,896]),
        support=dict(log_density=[-18,.02],temperature_K=[240,20000],neutral_fraction=[1e-9,.001],horizon_seconds=.0034344311179287023,domain_m=[-200,1200]),
        gates=dict(native_constitutive=.002,rate_interpolation=.002,quadrature=1e-5,equilibrium_balance=1e-9,primitive_energy=1e-8,baryon=1e-10,energy=1e-8,species=1e-9,mass_refinement=.02,trace_refinement=.02,wave_refinement=.02,optically_thin_absorption=.001),
        stop='No automatic support, grid or horizon increase. Preserve failures. Reassess missing physics rather than declaring a final charge. Do not substitute a fitted relaxation time or count rate-only diagnostics as coupled evolution.',
        bindings={str(p):old.cold.sha(p) for p in [Path(__file__),Path(old.__file__),old.OUT/'dilute-eos.npz',old.OUT/'dilute/cells-896.npz',ATOMIC,LEVELS,old.cold.CACHE/'libfree_eos_native_cold_stable.so']}))


def probe():
    assert not (OUT/'probe.json').exists();start=time.monotonic();signal.alarm(30);native=Native(cap=120)
    m=old.DiluteFlow(224).base;Trad=(float(m.env['Linfinity'])/(4*np.pi*m.RJ**2*m.a0**2*5.670374419e-5))**.25
    rows=[];states=[]
    for x,T,y in [(0.,native.fan.T,native.y0),(-2.,6000.,native.y0),(-6.,1500.,native.y0),(0.,native.fan.T,1e-9),(0.,native.fan.T,.001)]:
        a=native.state(x,np.log(T),y);r=native.rates(a,Trad);rr=native.rates(a,Trad,32)
        qerr=float(np.max(abs(r/rr-1)));assert qerr<1e-5
        # Detailed balance belongs to unbiased LTE: constrain this state's
        # *actual* zero-affinity snapshot, rather than setting a frozen ratio.
        eq=native.ion.snapshot(native.lr+x,np.log(T),np.zeros(318));eqy=eq['number_fractions'][0,0]/eq['number_fractions'][0,:2].sum()
        balanced=dict(a,y=eqy,affinity=0.);bb=native.rates(balanced,T,32)
        balance=float(np.max(abs(bb[0]/(bb[1]+bb[2])-1)));assert balance<1e-9
        rate=.5*r[0]-r[1]-.5*r[2]
        rows.append(dict(x=x,T=T,y=y,affinity=a['affinity'],number_per_H=rate[0],energy_erg_per_H=rate[1],photo_per_neutral=r[0,0]/y*.5,recombination_per_ion=(r[1,0]+.5*r[2,0])/(1-y),quadrature=qerr,balance=balance))
        states.append(a)
    np.savez_compressed(OUT/'probe-states.npz',**{k:np.array([s[k] for s in states]) for k in states[0]})
    write(OUT/'probe.json',dict(classification='Counterexample candidate',passed=True,Trad_K=Trad,y0=native.y0,nH_per_g=native.nH,rows=rows,EOS_calls=native.ion.calls,seconds=time.monotonic()-start,source_sha256=old.cold.sha(__file__)))
    signal.alarm(0);print((OUT/'probe.json').read_text(),flush=True)


def bank():
    assert not (OUT/'bank.npz').exists();start=time.monotonic();signal.alarm(115)
    probe=json.loads((OUT/'probe.json').read_text());native=Native(cap=3980);Trad=probe['Trad_K']
    # The probe measures ~50 ionizations per neutral in this horizon. At the
    # existing dilute cutoff the physical neutral fraction can fall below1e-9.
    # Reassess that composition support BEFORE evolving, without extra paths.
    write(OUT/'bank-plan.json',dict(classification='Counterexample candidate',
        neutral_fraction=[1e-16,.001],density_nodes=17,temperature_nodes=25,
        reason='Measured photo-rate times horizon is about50. The originally proposed1e-9 neutral floor excludes finite-ionization states in the already accepted dilute flow. Use native constraints down to1e-16 and explicitly require relative H abundance accuracy1e-8; do not clip evolved abundances.',
        reuse='Retain the entire accepted153-state fixed-chemistry table as independent interior-composition controls. Native cross sections and level routines are reused unchanged.',
        remaining_calls=3980,remaining_seconds=115,source_sha256=old.cold.sha(__file__)))
    saved=np.load(old.OUT/'dilute-eos.npz');x=saved['x'];lt=np.linspace(saved['lt'][0],saved['lt'][-1],25);ys=np.array([1e-16,.001])
    raw=np.zeros((2,len(x),len(lt),21));rates=np.zeros((2,len(x),len(lt),3,2));done=np.zeros(raw.shape[:3],bool)
    try:
        for iy,y in enumerate(ys):
            for i,xx in enumerate(x):
                for j,t in enumerate(lt):
                    a=native.state(float(xx),float(t),float(y));raw[iy,i,j]=a['raw'];r=native.rates(a,Trad)
                    r[0]/=y;r[1:]/=(1-y);rates[iy,i,j]=r;done[iy,i,j]=True
                np.savez_compressed(OUT/'bank-progress.npz',x=x,lt=lt,ys=ys,raw=raw,rates=rates,done=done,EOS_calls=native.ion.calls)
        np.savez_compressed(OUT/'bank.npz',x=x,lt=lt,ys=ys,raw=raw,rates=rates,rho0=native.fan.rho,cx=native.fan.cx,sunit=saved['sunit'],s0=saved['s0'],y0=native.y0,nH=native.nH,Trad=Trad)
        eos=EOS();checks=[]
        for xx,T,y in [(-17.,310.,2e-14),(-14.5,700.,1e-10),(-11.,1900.,1e-8),(-7.,5700.,native.y0),(-3.1,11200.,.0003),(-.15,13380.,3e-6),(-.01,18200.,.0007)]:
            t=np.log(T);a=native.state(xx,t,y);r=native.rates(a,Trad,32);eos.y=np.array([y]);p,u,g,_,_,cv,_=eos.evaluate(np.array([np.exp(xx)]),np.array([t]));est=eos.reactions(np.array([np.exp(xx)]),np.array([t]))[0]
            errors=[abs(p[0]*eos.rho0*C*C/a['raw'][1]-1),abs(u[0]*C*C/a['raw'][2]-1)]
            checks.append(dict(x=xx,T=T,y=y,constitutive=list(map(float,errors)),rate=float(np.max(abs(est/r-1)))))
        # Reuse every original fixed-composition native point, including cold
        # and dilute states, as an independent y-interpolation control.
        xx,tt=np.meshgrid(x,saved['lt'],indexing='ij');eos.y=np.full(xx.size,native.y0);p,u,*_=eos.evaluate(np.exp(xx.ravel()),tt.ravel())
        cached=max(np.max(abs(p.reshape(xx.shape)*eos.rho0*C*C/saved['raw'][:,:,1]-1)),np.max(abs(u.reshape(xx.shape)*C*C/saved['raw'][:,:,2]-1)))
        passed=bool(cached<.002 and max(max(a['constitutive']) for a in checks)<.002 and max(a['rate'] for a in checks)<.002)
        write(OUT/'bank.json',dict(classification='Counterexample candidate',passed=passed,checks=checks,cached153_constitutive_relative=float(cached),EOS_calls=native.ion.calls,seconds=time.monotonic()-start,source_sha256=old.cold.sha(__file__)))
        print((OUT/'bank.json').read_text(),flush=True);assert passed
    except Exception as exc:
        write(OUT/'bank-failure.json',dict(error=repr(exc),completed=int(done.sum()),EOS_calls=native.ion.calls,seconds=time.monotonic()-start));raise
    finally:
        np.savez_compressed(OUT/'bank-native-states.npz',**{k:np.array([s[k] for s in native.ion.states]) for k in native.ion.states[0]})
    signal.alarm(0)


class EOS:
    def __init__(self):
        from types import SimpleNamespace
        self.d=d=np.load(OUT/('repaired-bank.npz' if (OUT/'repaired-bank.npz').exists() else 'bank.npz'));self.x=d['x'];self.lt=d['lt'];self.ys=d['ys'];self.y0=float(d['y0']);self.y=self.y0
        self.rho0=float(d['rho0']);self.cx=float(d['cx']);self.sunit=float(d['sunit']);self.floor=np.exp(self.x[0]);self.top=np.exp(self.x[-1]);self.nH=float(d['nH']);self.fan=SimpleNamespace(calls=0)
        self.f=[];self.rf=[];self.R=[];self.u0=[]
        rho=self.rho0*np.exp(self.x[:,None]);T=np.exp(self.lt[None,:])
        for a,rate in zip(d['raw'],d['rates']):
            R=float(a[-1,-1,1]/a[-1,-1,0]/T[0,-1]);u0=float(a[-1,-1,2]-1.5*R*T[0,-1]);self.R.append(R);self.u0.append(u0)
            self.f.append([RectBivariateSpline(self.x,self.lt,v) for v in [np.log(a[:,:,1]/rho/T),a[:,:,2]-u0-1.5*R*T,a[:,:,3],a[:,:,13]/a[:,:,0]]])
            self.rf.append([RectBivariateSpline(self.x,self.lt,np.log(rate[:,:,i,j])) for i in range(3) for j in range(2)])

    def limits(self,rho):return np.full_like(rho,self.lt[0]),np.full_like(rho,self.lt[-1])

    def weights(self,rho,lt):
        active=rho>=self.floor;y=np.broadcast_to(self.y,rho.shape)
        assert np.all((y[active]>=self.ys[0]*(1-1e-6))&(y[active]<=self.ys[-1]*(1+1e-9))),('Hydrogen inventory support',float(min(y[active],default=self.y0)),float(max(y[active],default=self.y0)))
        assert max(rho)<=self.top*(1+1e-9),'Hydrogen density support'
        assert np.all((lt[active]>=self.lt[0]-1e-12)&(lt[active]<=self.lt[-1]+1e-12)),'Hydrogen temperature support'
        return np.log(np.maximum(rho,self.floor)),np.clip(lt,*self.lt[[0,-1]]),(y-self.ys[0])/np.diff(self.ys)[0],active,y

    def evaluate(self,rho,lt):
        x,t,w,active,y=self.weights(rho,lt);T=np.exp(t);rows=[]
        for i,(lp,uf,sf,kf) in enumerate(self.f):
            p=np.exp(lp.ev(x,t))*self.rho0*rho*T;u=self.u0[i]+1.5*self.R[i]*T+uf.ev(x,t);cv=1.5*self.R[i]*T+uf.ev(x,t,dy=1)
            rows.append(np.array([p,u,cv,p*(1+lp.ev(x,t,dx=1)),p*(1+lp.ev(x,t,dy=1)),uf.ev(x,t,dx=1),sf.ev(x,t),kf.ev(x,t)]))
        p,u,cv,pr,pt,ur,s,kap=rows[0]*(1-w)+rows[1]*w
        gamma=np.divide(pr+pt*(p/(self.rho0*np.maximum(rho,self.floor))-ur)/cv,p,out=np.full_like(p,5/3),where=active)
        assert np.all(cv[active]>0) and np.all(gamma[active]>1)
        return p/(self.rho0*C*C)*active,u/C**2*active,gamma,T,kap*6.6524587321e-25/1.66053906660e-24*active,cv/C**2,s

    def __call__(self,rho,lt):return self.evaluate(rho,lt)[:5]

    def reactions(self,rho,lt):
        x,t,w,active,y=self.weights(rho,lt)
        a=np.array([np.exp(f.ev(x,t)) for f in self.rf[0]]).T.reshape((-1,3,2))
        b=np.array([np.exp(f.ev(x,t)) for f in self.rf[1]]).T.reshape((-1,3,2))
        r=a*(1-w[:,None,None])+b*w[:,None,None];r[:,0]*=y[:,None];r[:,1:]*=(1-y[:,None,None]);r*=active[:,None,None]
        return r


def warm_repair():
    assert not (OUT/'repaired-bank.npz').exists();start=time.monotonic();signal.alarm(70)
    failed=json.loads((OUT/'bank.json').read_text());assert not failed['passed']
    write(OUT/'warm-interpolation-reassessment.json',dict(classification='Counterexample candidate',
        failure='Native constitutive and original153-state checks pass, but the independent18200K photo/reverse rate interpolation error0.006699 exceeds the unchanged0.002 gate. The failed bank verdict is preserved.',
        repair='Retain all850 completed native states. Add exactly the two midpoints of the warmest two temperature intervals,68 states, to resolve the strongly varying finite-level rates. Keep physical temperature/density/composition domains, production paths and all gates unchanged.',
        remaining_native_calls=880,remaining_native_seconds=70,new_native_states=68,stop='No further table refinement if the original independent controls still fail.',source_sha256=old.cold.sha(__file__)))
    native=Native(cap=880);d=dict(np.load(OUT/'bank.npz'));newt=(d['lt'][-3:-1]+d['lt'][-2:])/2
    raw=np.zeros((2,len(d['x']),2,21));rates=np.zeros((2,len(d['x']),2,3,2))
    for k,y in enumerate(d['ys']):
        for i,x in enumerate(d['x']):
            for j,t in enumerate(newt):
                a=native.state(float(x),float(t),float(y));raw[k,i,j]=a['raw'];r=native.rates(a,float(d['Trad']));r[0]/=y;r[1:]/=(1-y);rates[k,i,j]=r
    order=np.argsort(np.r_[d['lt'],newt]);d['lt']=np.r_[d['lt'],newt][order]
    d['raw']=np.concatenate([d['raw'],raw],axis=2)[:,:,order];d['rates']=np.concatenate([d['rates'],rates],axis=2)[:,:,order]
    np.savez_compressed(OUT/'repaired-bank.npz',**d);eos=EOS();checks=[]
    for row in failed['checks']:
        x,T,y=row['x'],row['T'],row['y'];a=native.state(x,np.log(T),y);r=native.rates(a,float(d['Trad']),32);eos.y=np.array([y])
        p,u,*_=eos(np.array([np.exp(x)]),np.array([np.log(T)]));estimate=eos.reactions(np.array([np.exp(x)]),np.array([np.log(T)]))[0]
        checks.append(dict(x=x,T=T,y=y,constitutive=float(max(abs(p[0]*eos.rho0*C*C/a['raw'][1]-1),abs(u[0]*C*C/a['raw'][2]-1))),rate=float(np.max(abs(estimate/r-1)))))
    passed=bool(max(a['rate'] for a in checks)<.002 and max(a['constitutive'] for a in checks)<.002)
    write(OUT/'repaired-bank.json',dict(classification='Counterexample candidate',passed=passed,checks=checks,EOS_calls=native.ion.calls,seconds=time.monotonic()-start,source_sha256=old.cold.sha(__file__)))
    np.savez_compressed(OUT/'repair-native-states.npz',**{k:np.array([s[k] for s in native.ion.states]) for k in native.ion.states[0]})
    signal.alarm(0);print((OUT/'repaired-bank.json').read_text(),flush=True);assert passed


if __name__=='__main__':globals()[sys.argv[1]]()
