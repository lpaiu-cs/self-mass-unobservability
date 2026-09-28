"""Counterexample candidate: a bounded, resolved native spherical pulse experiment.

Generalizes the frozen Phase39 first step to arbitrary previous scalar states.
Backward Euler matter and a symmetric scalar discrete gradient: globally first
order is expected. The controlled boundary is not an astrophysical companion.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from types import FunctionType, SimpleNamespace
import argparse
import inspect
import json
import time
import numpy as np
import sympy as sp
from scipy.optimize import brentq
from scipy.linalg import solve_banded
import def_spherical_regular as prior

old,m,e,ld=prior.old,prior.m,prior.e,prior.ld
ROOT=prior.ROOT;OUT=e.g.OUT/'def-resolved-scalar-pulse'
ATOL,SCALE=prior.ATOL,prior.SCALE


def coarse_geometry(n):
    data=np.load(old.imported.initial.OUT/'initial.npz')
    # Snap uniform-radius faces to old faces; sums preserve actual inventories.
    cuts=np.searchsorted(data['rf'],np.linspace(0,float(data['rf'][-1]),n+1))
    cuts[0]=0;cuts[-1]=len(data['r']);assert np.all(np.diff(cuts)>0)
    rf=data['rf'][cuts];V=4*np.pi/3*np.diff(rf**3)
    B=np.array([data['baryons'][a:b].sum() for a,b in zip(cuts[:-1],cuts[1:])],dtype=ld)
    U=np.array([data['energy'][a:b].sum() for a,b in zip(cuts[:-1],cuts[1:])],dtype=ld)
    X=np.array([np.sum(data['baryons'][a:b,None]*data['base'][a:b,5:],axis=0)/B[i]
        for i,(a,b) in enumerate(zip(cuts[:-1],cuts[1:]))],dtype=ld)
    rho=np.exp(data['base'][:,0]);rest=(data['base'][:,5:]/e.g.c.A)@e.g.c.W*e.C**2
    rest_new=(X/e.g.c.A)@e.g.c.W*e.C**2
    rest_energy=np.array([np.sum(rho[a:b]*rest[a:b]*data['volume'][a:b]) for a,b in zip(cuts[:-1],cuts[1:])])
    aeff=B*rest_new/rest_energy
    mf=np.r_[ld(0),np.cumsum(e.GRAV*U)]
    r=[]
    for i in range(n):
        target=1-1/aeff[i]**2
        def equation(x):
            rr=ld(x);mass=mf[i]+e.GRAV*U[i]*(rr**3-rf[i]**3)/(rf[i+1]**3-rf[i]**3)
            return float(2*mass/rr-target)
        grid=np.linspace(float(max(rf[i],rf[i+1]*ld('1e-8'))),float(rf[i+1]),65)
        brackets=[(a,b) for a,b in zip(grid[:-1],grid[1:]) if equation(a)*equation(b)<=0]
        assert brackets,('No moment-compatible radius',i)
        roots=[brentq(equation,a,b,xtol=1e-7,rtol=1e-14) for a,b in brackets]
        r.append(min(roots,key=lambda x:abs(x-float((rf[i]+rf[i+1])/2))))
    r=np.asarray(r,dtype=ld)
    frac=(r**3-rf[:-1]**3)/np.diff(rf**3)
    a=1/np.sqrt(1-2*(mf[:-1]+e.GRAV*U*frac)/r)
    lr=np.log(B/(a*V));target_u=U/(np.exp(lr)*V)-rest_new
    assert np.all(target_u>0)
    guess=np.array([np.sum(data['baryons'][a:b]*data['base'][a:b,1])/B[i] for i,(a,b) in enumerate(zip(cuts[:-1],cuts[1:]))])
    fine_a=data['baryons']/(rho*data['volume'])
    heat=np.array([[np.sum(fine_a[a:b]*data['qscale'][a:b]*data['base'][a:b,j]*data['volume'][a:b])
        for j in [3,4]] for a,b in zip(cuts[:-1],cuts[1:])])/(a*V)[:,None]
    return dict(rf=rf,r=r,volume=V,indices=np.arange(n),baryons=B,energy=U,X=X,
        logrho=lr,guess_logT=guess,target_u=target_u,heat=heat,
        original_baryons=data['original_baryons'],original_mass=data['original_mass'],cuts=cuts)


def metric(star,dEm,R,Pi,Phi):
    z=prior.metric(star,dEm,R,Pi,Phi)
    if not hasattr(star,'previous'):return z
    p=star.previous;f=star.fraction
    da=z['da']-p['da'];x=da/p['a'];ratio=np.empty_like(x);small=abs(x)<ld('1e-6')
    ratio[small]=1-x[small]/2+x[small]**2/3-x[small]**3/4
    ratio[~small]=np.log1p(x[~small])/x[~small]
    Rstar=(p['a']*p['R']+z['a']*R)/2*ratio/p['a']
    K=star.scalar_volume*(Pi**2+p['Pi']**2)/(16*np.pi)
    face=star.face_weight*(Phi**2+p['Phi']**2)/(16*np.pi)
    L,Rg=face[:-1]/2,face[1:]/2;Rg[-1]*=2
    C=1-2*K*(1-f)/star.r-2*L/np.where(star.rf[:-1]>0,star.rf[:-1],1)
    D=1+2*K*f/star.r+2*Rg/star.rf[1:]
    chi=4/(np.sqrt(p['b'])+np.sqrt(z['b']))**2
    work=e.GRAV*star.volume*(star.reference['E']+(dEm+p['dEm'])/2+Rstar)*chi/star.r
    C-=work*(1-f);D+=work*f;assert min(C.min(),D.min())>0
    lam=np.ones(star.n+1,dtype=ld);H=np.empty(star.n,dtype=ld)
    for i in range(star.n-1,-1,-1):H[i]=lam[i+1]/D[i];lam[i]=H[i]*C[i]
    Hgrad=np.r_[H[0],(H[:-1]+H[1:])/2,H[-1]]
    z.update(H=H,Hface=lam,Hgrad=Hgrad,N=H/z['a'],Rstar=Rstar,pair_da=da,
        k=H*(p['b']+z['b'])/2,kf=Hgrad*(p['bf']+z['bf'])/2)
    return z


def wave(star,z):
    p=star.previous;k=z['k'];h=star.h
    coeff=h*h*e.C**2*k/4
    source=z['H']*e.GRAV*star.volume/(star.scalar_volume/(4*np.pi))*star.beta/(p['a']+z['a'])
    conductance=star.scalar_area*z['kf']/star.distance
    factor=coeff/star.scalar_volume
    lower=-factor[1:]*conductance[1:-1];upper=-factor[:-1]*conductance[1:-1]
    diagonal=1+factor*(conductance[:-1]+conductance[1:])-2*coeff*source*z['a']*z['trace']
    rhs=p['psi']+h*e.C*k*p['Pi']/2+coeff*source*(p['a']*p['trace']-z['a']*z['trace'])*p['psi']
    rhs[-1]+=factor[-1]*conductance[-1]*(star.amplitude+star.previous_boundary)/2
    band=np.zeros((3,star.n));band[1]=diagonal;band[0,1:]=upper;band[2,:-1]=lower
    mid=solve_banded((1,1),band,np.asarray(rhs,float)).astype(ld)
    for _ in range(3):
        defect=rhs-diagonal*mid;defect[1:]-=lower*mid[:-1];defect[:-1]-=upper*mid[1:]
        mid+=solve_banded((1,1),band,np.asarray(defect,float)).astype(ld)
    psi=2*mid-p['psi'];Pi=2*(psi-p['psi'])/(h*e.C*k)-p['Pi']
    Phi=np.r_[ld(0),np.diff(psi)/np.diff(star.r),(star.amplitude-psi[-1])/star.distance[-1]]
    defect=diagonal*mid-rhs;defect[1:]+=lower*mid[:-1];defect[:-1]+=upper*mid[1:]
    return psi,Pi,Phi,float(abs(defect).max()/max(star.peak,ld('1e-100')))


class PulseStar(prior.RegularStar):
    def native_aux(self,delta,psi):
        return super().native_aux(delta,psi if self.beta else np.zeros_like(psi))

    def fluid(self,delta,psi,raw):
        z=super().fluid(delta,psi if self.beta else np.zeros_like(psi),raw)
        z['psi']=psi
        return z

    def evaluate(self,delta):
        assert abs(delta[:,2]).max()<.2
        X=self.base[:,5:]+delta[:,5:];assert X.min()>=-1e-15 and abs(X.sum(1)-1).max()<1e-10
        if hasattr(self,'linearization'):
            anchor=self.linearization;z=self.fluid(delta,anchor['psi'],self.tangent_aux(delta))
            z.update({k:anchor[k] for k in self.metric_keys+['Pi','Phi','scalar_residual']})
        else:
            psi,Pi,Phi=[v.copy() for v in self.field_guess]
            for _ in range(12):
                z=self.fluid(delta,psi,self.native_aux(delta,psi));z.update(metric(self,z['dEm'],z['R'],Pi,Phi))
                psin,pin,phin,error=wave(self,z)
                change=abs(psin-psi).max()/max(self.peak,ld('1e-100'))
                psi,Pi,Phi=psin,pin,phin
                if change<ld('2e-19'):break
            else:raise RuntimeError(('Scalar fixed point',float(change)))
            z=self.fluid(delta,psi,self.native_aux(delta,psi));z.update(metric(self,z['dEm'],z['R'],Pi,Phi))
            pn,qn,fn,error=wave(self,z)
            error=max(error,float(abs(pn-psi).max()/max(self.peak,ld('1e-100'))))
            assert error<2e-17,error
            z.update(Pi=Pi,Phi=Phi,scalar_residual=error);self.field_guess=psi,Pi,Phi
        p=self.previous
        z['gstar']=self.beta*(p['a']*p['psi']*p['trace']+z['a']*z['psi']*z['trace'])/(p['a']+z['a'])
        self.finish_moments(delta,z,z['dEm'])
        if hasattr(self,'linearization'):m.finish_fluxes(self,z,self.faces(anchor['N']*anchor['v']/anchor['a'],odd=True))
        z['at']=z['pair_da']/self.h;self.current=z
        return z


def residual(star,delta,previous,older,h,coefficients):
    value,z=old.matter_residual(star,delta,previous,older,h,coefficients);p=previous[1]
    balance=z['dEm']-p['dEm']+(star.reference['E']+(z['dEm']+p['dEm'])/2+z['Rstar'])*z['pair_da']/((p['a']+z['a'])/2)
    balance+=z['gstar']*(z['psi']-p['psi'])+h*e.C*star.divergence(z['Hface']*z['fluxes'][1])/z['H']
    value[:,1]=balance/star.heat0
    force=-z['at']*z['S']-e.C*(z['N']*z['nur']*z['E']-star.equilibrium_gravity)
    force+=e.C*(z['N']*z['P']-star.equilibrium_pressure)*star.area_difference_over_volume
    force+=e.C*z['N']*star.beta*z['psi']*z['trace']*(z['Phi'][:-1]+z['Phi'][1:])/2
    value[:,2]=(z['dU'][:,2]-p['dU'][:,2]+h*e.C*star.divergence(z['fluxes'][2]-star.equilibrium_pressure_flux)-h*force)/star.momentum_scale
    return value,z


class Tangent(PulseStar):
    evaluate=prior.Tangent.evaluate


# prior.Tangent uses zero-argument super bound to its original class; use the
# same tiny adapter with this class rather than mutating any frozen module.
def tangent_evaluate(self,delta):
    raw=self.tangent_aux(delta);projected=self.projected_delta;self.tangent_aux=lambda _:raw
    try:return PulseStar.evaluate(self,projected)
    finally:del self.tangent_aux
Tangent.evaluate=tangent_evaluate
tangent=FunctionType(old.tangent.__code__,dict(vars(old),Tangent=Tangent))
jacobian=FunctionType(m.prior.jacobian.__code__,dict(vars(m.prior),residual=residual))
stage=FunctionType(prior.stage.__code__,dict(prior.stage.__globals__,residual=residual,tangent=tangent,jacobian=jacobian))


def initialize(pool,beta,peak,h):
    saved=dict(np.load(OUT/'initial.npz'))
    star=old.imported.initial.original.make_star(old.imported.cached.CachedMaterialPool(pool),saved)
    fn=m.prior.initialize
    star=FunctionType(fn.__code__,dict(fn.__globals__,parent=SimpleNamespace(initialize=lambda _:star)))(pool)
    star.__class__=PulseStar;star.h=h;star.beta=ld(beta);star.peak=ld(peak);star.amplitude=ld(0);star.previous_boundary=ld(0)
    star.area=4*np.pi*star.rf**2;star.scalar_volume=4*np.pi*star.r**2*np.diff(star.rf)
    star.scalar_area=4*np.pi*np.r_[ld(0),star.r[:-1]*star.r[1:],star.r[-1]*star.rf[-1]]
    star.distance=np.r_[ld(1),np.diff(star.r),star.rf[-1]-star.r[-1]];star.face_weight=star.scalar_area*star.distance
    ref=star.reference;star.b0=1/ref['a']**2;star.bf0=np.r_[ld(1),1-2*ref['mf'][1:]/star.rf[1:]]
    zeros=np.zeros(star.n,dtype=ld);faces=np.zeros(star.n+1,dtype=ld)
    met=prior.metric(star,zeros,ref['R'],zeros,faces);ref.update(met);m.attach_equilibrium(star)
    z=star.fluid(np.zeros_like(star.base),zeros,star.initial_aux.copy());z.update(met,Pi=zeros.copy(),Phi=faces,scalar_residual=0.,pair_da=zeros.copy())
    star.finish_moments(np.zeros_like(star.base),z,zeros);z['at']=zeros.copy()
    star.previous=z;star.field_guess=zeros.copy(),zeros.copy(),faces.copy()
    star.metric_keys=list(met)+['pair_da'];star.initial_guess=np.zeros_like(star.base)
    return star,z


def symbolic():
    p0,p1,a0,a1,t0,t1,beta,x=sp.symbols('p0 p1 a0 a1 t0 t1 beta x')
    g=beta*(a0*p0*t0+a1*p1*t1)/(a0+a1)
    assert sp.simplify(g.subs(p1,2*x-p0)-beta*(2*a1*t1*x+(a0*t0-a1*t1)*p0)/(a0+a1))==0
    U,B,C,R=sp.symbols('U B C R',positive=True)
    ae=B*C/R;rho=B/ae
    assert sp.simplify(U/rho-C-ae*(U-R)/B)==0
    # sin^8 has the first seven derivatives zero at a pulse endpoint.
    q=sp.Symbol('q');pulse=sp.sin(sp.pi*q)**8
    assert all(sp.diff(pulse,q,k).subs(q,0)==0 for k in range(8))
    return dict(classification='Proven',passed=True,scope='General affine reciprocal source, positive coarse thermal energy under rest-energy metric matching, and vanishing finite drive jet through order seven at cutoff. No continuum or physical signal claim.')


def prepare():
    assert not OUT.exists();geometry=coarse_geometry(24);OUT.mkdir()
    np.savez_compressed(OUT/'geometry.npz',**geometry)
    files=[Path(__file__),Path(prior.__file__),Path(old.__file__),Path(m.__file__),Path(m.prior.__file__),Path(m.method.__file__),Path(old.two.__file__),
        old.imported.initial.OUT/'initial.npz',old.imported.initial.OUT/'initial-manifest.json',OUT/'geometry.npz']
    plan=dict(classification='Counterexample candidate',bindings={p.relative_to(ROOT).as_posix():e.digest(p) for p in files},cells=24,
        pulse_duration_crossing_times=2,total_duration_crossing_times=4,amplitude='0.0002',beta=-4,
        time_steps=[48,96,192],paths=['driven','undriven','decoupled'],
        drive='Prescribed scalar boundary sin(pi*t/T)^8 for 0<t<T, exactly zero after T; initial scalar zero. Finite-energy boundary-controlled scattering experiment, not a companion waveform or measured amplitude.',
        primary_readout='Baryon-weighted RMS Jordan log-temperature difference. After cutoff compare driven-minus-undriven with beta=0 scalar-stress-only control. No scalar charge or free-fall inference from temperature alone.',
        comparator='Same time-dependent undriven background plus any fixed polynomial in boundary drive and its first four time derivatives, with zero increment when the drive jet is zero. Tail readout tests this restricted class only; propagating scalar and acoustic memory are not automatically a thermal relaxation mode.',
        gates=dict(native=ATOL.astype(float).tolist(),scalar=2e-17,identity=2e-15,
            minimum_time_order=.8,maximum_relative_time_difference=.1,minimum_control_separation_over_time_error=10.,minimum_Jordan_temperature_rms=1e-7),
        symbolic=symbolic(),budget=dict(initial_and_pilot_timeout_seconds=120,workers=4,blas_threads=1,
            maximum_production_timeout_seconds=900,production_requires_measured_pilot=True,maximum_path_attempts=9,automatic_expansion=False),
        scope='New 24-cell conservative discretization of saved baryons, species and material energy. First-order finite-time response at fixed spatial grid, not spatial convergence, physical EOS certification, orbital matching or observational closure.')
    e.write(OUT/'plan.json',plan);print('PREPARED',json.dumps(plan),flush=True)


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,sha in p['bindings'].items():assert e.digest(ROOT/rel)==sha,rel
    return p


def build_initial(pool):
    p=bindings();g=dict(np.load(OUT/'geometry.npz'));began=time.monotonic()
    rows=list(pool.map(old.imported.initial.original.restrict_cell,zip(range(p['cells']),g['logrho'],g['guess_logT'],g['target_u'],g['X']),chunksize=1))
    aux=np.asarray([r[1] for r in rows],dtype=ld);rho=np.exp(g['logrho']);rest=(g['X']/e.g.c.A)@e.g.c.W*e.C**2
    qscale=rho*(rest+aux[:,2])+aux[:,1]
    base=np.column_stack([g['logrho'],[r[0] for r in rows],np.zeros(p['cells']),g['heat']/qscale[:,None],g['X']]).astype(ld)
    g.update(aux=aux,base=base,qscale=qscale);np.savez_compressed(OUT/'initial.npz',**g)
    star,z=initialize(pool,-4,ld(p['amplitude']),ld(1))
    B=z['B']*star.volume;U=z['E']*star.volume
    report=dict(classification='Counterexample candidate',cells=star.n,native_calls=sum(r[3] for r in rows),seconds=time.monotonic()-began,
        baryon_relative_defect=float(abs(B/g['baryons']-1).max()),energy_relative_defect=float(abs(U/g['energy']-1).max()),
        maximum_native_thermal_residual=max(abs(r[2]) for r in rows),radius_cm=float(star.rf[-1]),crossing_seconds=float(star.rf[-1]/e.C))
    report['passed']=bool(report['baryon_relative_defect']<1e-11 and report['energy_relative_defect']<1e-11 and report['maximum_native_thermal_residual']<1e-8)
    e.write(OUT/'initial-result.json',report);assert report['passed'],report
    e.write(OUT/'initial-manifest.json',dict(sha256={q.relative_to(ROOT).as_posix():e.digest(q) for q in [OUT/'initial.npz',OUT/'initial-result.json']}))
    print('INITIAL',json.dumps(report),flush=True)


def run_path(pool,label,steps,limit=None):
    p=bindings();folder=OUT/f'{label}-{steps}';assert not folder.exists();folder.mkdir()
    R=ld(json.loads((OUT/'initial-result.json').read_text())['radius_cm']);tc=R/e.C;h=4*tc/steps
    peak=ld(0) if label=='undriven' else ld(p['amplitude']);beta=0 if label=='decoupled' else -4
    star,z=initialize(pool,beta,peak,h);delta=np.zeros_like(star.base);start=time.monotonic()
    history=[];budgets=[];records=[]
    try:
        for j in range(1,(limit or steps)+1):
            previous=(delta,z);star.previous=z;star.initial_guess=delta.copy()
            star.previous_boundary=star.amplitude
            t=ld(j)*h;star.amplitude=peak*np.sin(ld(str(np.pi))*t/(2*tc))**8 if j<steps//2 else ld(0)
            def log(row):
                row['step']=j;records.append(row)
                with (folder/'iterations.jsonl').open('a') as stream:stream.write(json.dumps(row)+'\n')
            delta,z,value=stage(star,previous,log)
            before=previous[1]
            work=z['Hgrad'][-1]*(z['bf'][-1]+before['bf'][-1])/2*star.scalar_area[-1]/(4*np.pi*e.GRAV)*(z['Phi'][-1]+before['Phi'][-1])/2*(star.amplitude-star.previous_boundary)
            dm=(z['dmf'][-1]-before['dmf'][-1])/e.GRAV
            res=np.sum(z['H']*star.volume*star.heat0*value[:,1]);defect=dm-work-res
            scale=np.sum(star.volume*(abs(z['Ephi'])+abs(before['Ephi'])+abs(z['dEm']-before['dEm'])))+abs(work)
            identity=float(abs(defect)/max(scale,ld('1e-100')))
            assert identity<2e-15,(j,identity)
            history.append(dict(t=float(t),boundary=float(star.amplitude),native_norm=float(np.max(abs(value)/ATOL)),scalar=z['scalar_residual'],identity=identity))
            budgets.append([float(work),float(dm),float(res),float(defect)])
            np.savez_compressed(folder/f'step-{j:03d}.npz',delta=delta,residual=value,**{k:v for k,v in z.items() if isinstance(v,np.ndarray)})
            print('PULSE',label,steps,j,history[-1]['native_norm'],flush=True)
        result=dict(classification='Counterexample candidate',completed_steps=len(history),steps=steps,seconds=time.monotonic()-start,
            native_calls=star.pool.evaluations,history=history,budgets=budgets)
        e.write(folder/'result.json',result);return result
    except Exception as error:
        if hasattr(star,'last_state'):np.savez_compressed(folder/'failed-iterate.npz',delta=star.last_delta,**{k:v for k,v in star.last_state.items() if isinstance(v,np.ndarray)})
        e.write(folder/'failure.json',dict(error=repr(error),seconds=time.monotonic()-start,step=j,native_calls=star.pool.evaluations));raise


def pilot():
    p=bindings();assert not (OUT/'pilot.json').exists()
    with ProcessPoolExecutor(max_workers=4,initializer=old.imported.initial.original.worker_init) as pool:
        if not (OUT/'initial.npz').exists():build_initial(pool)
        result=run_path(pool,'driven',48,limit=4)
    e.write(OUT/'pilot.json',dict(classification='Counterexample candidate',passed=True,**{k:result[k] for k in ['seconds','native_calls','completed_steps']},
        scope='Timing only; early small pulse samples do not establish response amplitude. Production must account for later stronger field separately.'))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','pilot'])
    globals()[parser.parse_args().action]()
