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


if __name__=='__main__':globals()[sys.argv[1]]()
