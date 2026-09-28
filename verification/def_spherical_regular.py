"""Counterexample candidate: regular scalar plus native spherical matter, one step.

The mass adjoint, finite material work and shared transport faces are one
discretization. No division by scalar momentum or compensating matter heating.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from types import FunctionType, SimpleNamespace
import argparse
import inspect
import json
import time

import numpy as np
import def_spherical_coupled as old
import def_heat_metric_secant as secant

m,e,ld=old.m,old.e,old.ld
ATOL,SCALE=old.ATOL,old.SCALE
ROOT=e.ROOT
OUT=e.g.OUT/'def-spherical-regular'


def metric(star,dEm,R,Pi,Phi):
    ref=star.reference; f=star.fraction; gV=e.GRAV*star.volume
    K=star.scalar_volume*Pi**2/(8*np.pi)
    face=star.face_weight*Phi**2/(8*np.pi)
    L,Rg=face[:-1]/2,face[1:]/2; Rg[-1]*=2
    left=2*L/np.where(star.rf[:-1]>0,star.rf[:-1],1)
    C=1-2*K*(1-f)/star.r-left
    D=1+2*K*f/star.r+2*Rg/star.rf[1:]
    source=gV*dEm+K*star.b0+L*star.bf0[:-1]+Rg*star.bf0[1:]
    dmf=np.zeros(star.n+1,dtype=ld)
    for i in range(star.n): dmf[i+1]=(C[i]*dmf[i]+source[i])/D[i]
    dm=(1-f)*dmf[:-1]+f*dmf[1:]
    db=-2*dm/star.r; dbf=np.r_[ld(0),-2*dmf[1:]/star.rf[1:]]
    b,bf=star.b0+db,star.bf0+dbf
    assert min(b.min(),bf.min())>0
    root=np.sqrt(1+db/star.b0)
    da=ref['a']*(-db/star.b0)/(root*(1+root)); a=ref['a']+da
    x=da/ref['a']; ratio=np.empty_like(x); small=abs(x)<ld('1e-6')
    ratio[small]=1-x[small]/2+x[small]**2/3-x[small]**3/4
    ratio[~small]=np.log1p(x[~small])/x[~small]
    Rstar=(ref['a']*ref['R']+a*R)/2*ratio/ref['a']
    chi=4/(np.sqrt(star.b0)+np.sqrt(b))**2
    work=gV*(ref['E']+dEm/2+Rstar)*chi/star.r
    Cs=(1+C)/2-work*(1-f); Ds=(1+D)/2+work*f
    assert min(Cs.min(),Ds.min())>0
    lam=np.ones(star.n+1,dtype=ld); H=np.empty(star.n,dtype=ld)
    for i in range(star.n-1,-1,-1):
        H[i]=lam[i+1]/Ds[i]; lam[i]=H[i]*Cs[i]
    Hgrad=np.r_[H[0],(H[:-1]+H[1:])/2,H[-1]]
    Ephi=(K*b+L*bf[:-1]+Rg*bf[1:])/gV
    Et=ref['E']+dEm+Ephi; Rt=R+Ephi; mass=ref['m']+dm
    return dict(a=a,da=da,b=b,bf=bf,db=db,dbf=dbf,m=mass,mf=ref['mf']+dmf,
        dmf=dmf,N=H/a,H=H,Hface=lam,Hgrad=Hgrad,Rstar=Rstar,
        k=H*(star.b0+b)/2,kf=Hgrad*(star.bf0+bf)/2,
        Etot=Et,Ephi=Ephi,dEt=dEm+Ephi,
        nur=a*a*(mass/star.r**2+4*np.pi*e.GRAV*star.r*Rt),
        ar=a*a*(4*np.pi*e.GRAV*star.r*Et-mass/star.r**2))


def scalar_coefficients(star,z,trace):
    k,kf=z['k'],z['kf']
    force=star.h**2*e.C**2*k/2*z['H']*e.GRAV*star.volume/(star.scalar_volume/(4*np.pi))
    force*=old.BETA*z['a']*trace/(star.reference['a']+z['a'])
    conductance=star.scalar_area*kf/star.distance
    factor=star.h**2*e.C**2*k/(4*star.scalar_volume)
    lower=-factor[1:]*conductance[1:-1]; upper=-factor[:-1]*conductance[1:-1]
    diagonal=1+factor*(conductance[:-1]+conductance[1:])-force
    rhs=np.zeros(star.n,dtype=ld);rhs[-1]=factor[-1]*conductance[-1]*star.amplitude/2
    return lower,diagonal,upper,rhs,k


# Same tridiagonal solve and long-double residual refinement, new coefficients.
wave=FunctionType(old.wave.__code__,dict(vars(old),scalar_coefficients=scalar_coefficients))


class RegularStar(old.CoupledStar):
    def evaluate(self,delta):
        assert abs(delta[:,2]).max()<.2
        X=self.base[:,5:]+delta[:,5:]
        assert X.min()>=-1e-15 and abs(X.sum(1)-1).max()<1e-10
        if hasattr(self,'linearization'):
            anchor=self.linearization
            z=self.fluid(delta,anchor['psi'],self.tangent_aux(delta))
            z.update({k:anchor[k] for k in self.metric_keys+['Pi','Phi','scalar_residual']})
        else:
            psi,Pi,Phi=[v.copy() for v in self.field_guess]
            for _ in range(12):
                z=self.fluid(delta,psi,self.native_aux(delta,psi))
                z.update(metric(self,z['dEm'],z['R'],Pi,Phi))
                psin,pin,phin=wave(self,z,z['trace'])
                error=abs(psin-psi).max()/max(abs(self.amplitude),ld('1e-100'))
                psi,Pi,Phi=psin,pin,phin
                if error<ld('2e-19'):break
            else:raise RuntimeError(('Scalar fixed point',float(error)))
            z=self.fluid(delta,psi,self.native_aux(delta,psi))
            z.update(metric(self,z['dEm'],z['R'],Pi,Phi))
            lo,di,up,rhs,k=scalar_coefficients(self,z,z['trace'])
            defect=di*psi/2-rhs;defect[1:]+=lo*psi[:-1]/2;defect[:-1]+=up*psi[1:]/2
            error=max(abs(defect).max(),abs(psi-self.h*e.C*k*Pi/2).max())/max(abs(self.amplitude),ld('1e-100'))
            assert error<ld('2e-17'),float(error)
            z.update(Pi=Pi,Phi=Phi,scalar_residual=float(error))
            self.field_guess=psi,Pi,Phi
        z['gstar']=old.BETA*z['psi']*z['a']*z['trace']/(self.reference['a']+z['a'])
        z['scalar_work']=np.zeros(self.n,dtype=ld)
        self.finish_moments(delta,z,z['dEm'])
        if hasattr(self,'linearization'):
            m.finish_fluxes(self,z,self.faces(anchor['N']*anchor['v']/anchor['a'],odd=True))
        z['at']=z['da']/self.h
        self.current=z
        return z


def residual(star,delta,previous,older,h,coefficients):
    value,z=old.residual(star,delta,previous,older,h,coefficients)
    assert np.all(previous[0]==0) and np.all(older[0]==0), 'First step only'
    ref=star.reference;abar=(ref['a']+z['a'])/2
    balance=z['dEm']+(ref['E']+z['dEm']/2+z['Rstar'])*z['da']/abar+z['gstar']*z['psi']
    balance+=h*e.C*star.divergence(z['Hface']*z['fluxes'][1])/z['H']
    value[:,1]=balance/star.heat0
    return value,z


class Tangent(RegularStar):
    def evaluate(self,delta):
        raw=self.tangent_aux(delta);projected=self.projected_delta
        self.tangent_aux=lambda _:raw
        try:return super().evaluate(projected)
        finally:del self.tangent_aux


tangent=FunctionType(old.tangent.__code__,dict(vars(old),Tangent=Tangent))
jacobian=FunctionType(m.prior.jacobian.__code__,dict(vars(m.prior),residual=residual))
# Reuse the native Newton budget and line search; the previous Cauchy state
# remains zero even when the old-model endpoint is an initial guess.
source=inspect.getsource(old.stage).replace('delta = previous[0].copy()', 'delta = star.initial_guess.copy()')
source=source.replace('    star.composition_data = composition_response(star,delta,z)',
    '    if np.max(abs(value)/ATOL)>1:\n        star.composition_data = composition_response(star,delta,z)')
exec(compile(source,__file__,'exec'),dict(vars(old),residual=residual,tangent=tangent,jacobian=jacobian),ns:={})
stage=ns['stage']


def initialize(pool,amplitude,h,cells):
    if cells==5735:
        star=old.imported.initialize(pool)
    else:
        saved=dict(np.load(old.imported.initial.OUT/'initial.npz'))
        for key in ['r','volume','indices','base','qscale','aux']:saved[key]=saved[key][:cells]
        saved['rf']=saved['rf'][:cells+1]
        star=old.imported.initial.original.make_star(old.imported.cached.CachedMaterialPool(pool),saved)
        fn=m.prior.initialize
        star=FunctionType(fn.__code__,dict(fn.__globals__,parent=SimpleNamespace(initialize=lambda _:star)))(pool)
    star.__class__=RegularStar;star.h,star.amplitude=h,amplitude
    star.area=4*np.pi*star.rf**2
    star.scalar_volume=4*np.pi*star.r**2*np.diff(star.rf)
    star.scalar_area=4*np.pi*np.r_[ld(0),star.r[:-1]*star.r[1:],star.r[-1]*star.rf[-1]]
    star.distance=np.r_[ld(1),np.diff(star.r),star.rf[-1]-star.r[-1]]
    star.face_weight=star.scalar_area*star.distance
    ref=star.reference;star.b0=1/ref['a']**2
    star.bf0=np.r_[ld(1),1-2*ref['mf'][1:]/star.rf[1:]]
    zero=np.zeros(cells,dtype=ld);faces=np.zeros(cells+1,dtype=ld)
    met=metric(star,zero,ref['R'],zero,faces);star.metric_keys=list(met);ref.update(met)
    star.field_guess=zero.copy(),zero.copy(),faces.copy()
    star.initial_guess=np.zeros_like(star.base)
    m.attach_equilibrium(star)
    return star


def assess(star,delta,z,value):
    ref=star.reference
    boundary=z['Hgrad'][-1]*(star.bf0[-1]+z['bf'][-1])/2*star.scalar_area[-1]/(4*np.pi)*z['Phi'][-1]/2*star.amplitude/e.GRAV
    material_flux=star.h*e.C*star.area*z['Hface']*z['fluxes'][1]
    total=z['dmf'][-1]/e.GRAV
    defect=total-boundary+material_flux[-1]-material_flux[0]
    residual_energy=np.sum(z['H']*star.volume*star.heat0*value[:,1])
    scalar_scale=np.sum(star.volume*z['Ephi'])+abs(boundary)
    identity_scale=max(scalar_scale+np.sum(abs(z['H']*star.volume*z['dEm'])),ld('1e-100'))
    identity_error=abs(defect-residual_energy)/identity_scale
    allowance=np.sum(z['H']*star.volume*star.heat0)*ATOL[1]+ld('2e-15')*scalar_scale
    # The original baryon/isotope tolerances are reused; old energy flux is
    # replaced by the explicitly registered material work equation above.
    baryon=star.h*e.C*star.area*z['fluxes'][0]
    berror=(z['dU'][:,0]+np.diff(baryon)/star.volume)/star.B0
    inventory=np.sum(z['dBX']*star.volume[:,None],axis=0)/np.sum(star.B0*star.volume)
    cone_source=inspect.getsource(old.two.cones).replace('def cones(', 'def frame_cones(').replace('e.TAU',"z['tau_cond']")
    namespace=dict(vars(old.two));exec(compile(cone_source,__file__,'exec'),namespace)
    cone=namespace['frame_cones'](z)
    result=dict(classification='Counterexample candidate',cells=star.n,
        native_residual_norm=float(np.max(abs(value)/ATOL)),scalar_residual=z['scalar_residual'],
        mass_boundary_identity_relative_error=float(identity_error),energy_allowance_ratio=float(abs(defect)/allowance),
        boundary_scalar_work_erg=float(boundary),mass_energy_increment_erg=float(total),energy_defect_erg=float(defect),
        residual_energy_erg=float(residual_energy),scalar_energy_erg=float(np.sum(star.volume*z['Ephi'])),
        maximum_baryon_relative_defect=float(abs(berror).max()),isotope_inventory_defect=float(abs(inventory).max()),
        cone=cone,max_phi=float(abs(z['psi']).max()),maximum_primitive_changes=np.max(abs(delta[:,:5]),axis=0).astype(float).tolist(),
        native_calls=star.pool.evaluations,minimum_metric_b=float(z['b'].min()),minimum_adjoint=float(z['H'].min()))
    result['passed']=bool(result['native_residual_norm']<=1 and z['scalar_residual']<2e-17
        and identity_error<ld('2e-15') and abs(defect)<=allowance and abs(berror).max()<ld('1e-9')
        and abs(inventory).max()<ld('1e-9') and cone['sampled_cone_inside_light_cone'])
    return result


def prepare(cells):
    folder=OUT/f'cells-{cells}';assert not folder.exists()
    proofs=secant.check()
    files=[Path(__file__),Path(old.__file__),Path(secant.__file__),Path(old.imported.__file__),
        Path(m.__file__),Path(m.prior.__file__),Path(m.method.__file__),Path(old.two.__file__),
        old.imported.initial.OUT/'initial-manifest.json',old.imported.initial.OUT/'initial.npz']
    budget=dict(workers=2 if cells<5735 else 8,blas_threads=1,hard_timeout_seconds=90 if cells<5735 else 900,
        maximum_paired_runs=1,maximum_iterations_per_path=24,automatic_expansion=False)
    if cells==5735:
        pilot=json.loads((OUT/'cells-16/result.json').read_text());assert pilot['passed']
        budget.update(pilot_seconds=pilot['seconds'],pilot_native_calls=sum(p['native_calls'] for p in pilot['paths']),
            estimated_wall_seconds=[150,850],estimate_basis='Measured 16-cell assembly plus prior 5735-cell native first steps 142-330 seconds per path. New full-grid scans and new lapse convergence unmeasured; 900-second hard cap. Old endpoints only initial guesses; exact raw cache reuse.')
        files += [OUT/'cells-16/result.json',old.OUT/'path-0/endpoint.npz',e.g.OUT/'def-spherical-balanced-wave/endpoint.npz']
    folder.mkdir(parents=True)
    e.write(folder/'plan.json',dict(classification='Counterexample candidate',cells=cells,h_seconds=str(e.TAU/128),
        amplitudes=['0','0.000001'],bindings={p.relative_to(ROOT).as_posix():e.digest(p) for p in files},budget=budget,
        claim='One native matter, two heat currents, 26-species transport and regular scalar spherical first step with radial metric constraint and a shared mass adjoint.',
        decision='Accept this finite assembly only if native residual, scalar equation, weighted mass-boundary identity, baryons, species and local rest cone all pass; otherwise retain failure and stop before larger runs.',
        tolerances=dict(native=ATOL.astype(float).tolist(),scalar_equation=2e-17,mass_boundary_identity=2e-15,
            energy_allowance='sum(H V heat0)*1e-9 + 2e-15*(scalar_energy+abs(boundary_scalar_work))',
            baryons=1e-9,species=1e-9),symbolic=proofs,
        domain='Full original star' if cells==5735 else 'Inner 16 original cells with a new reflecting matter wall and lapse anchor at the truncated outer radius; assembly control only.',
        limitations=['First step only, no time convergence.','Radial mass constraint and discrete total energy only, no full Einstein constraint convergence.',
            'Scalar boundary ramp is a numerical control, not a companion or physical perturbation calibration.','No physical EOS certification or observational inference.']))
    print('PREPARED',folder,flush=True)


def run(cells):
    folder=OUT/f'cells-{cells}';plan=json.loads((folder/'plan.json').read_text())
    for rel,sha in plan['bindings'].items():assert e.digest(ROOT/rel)==sha,rel
    assert not (folder/'result.json').exists() and not (folder/'failure.json').exists()
    began=time.monotonic();results=[];endpoints=[]
    with ProcessPoolExecutor(max_workers=plan['budget']['workers'],initializer=old.imported.initial.original.worker_init) as pool:
        for index,amp in enumerate(plan['amplitudes']):
            path=folder/f'path-{index}';path.mkdir();rows=[]
            star=initialize(pool,ld(0),ld(plan['h_seconds']),cells)
            zero=np.zeros_like(star.base);z0=star.evaluate(zero);star.amplitude=ld(amp)
            if cells==5735:
                warm=old.OUT/'path-0/endpoint.npz' if index==0 else e.g.OUT/'def-spherical-balanced-wave/endpoint.npz'
                saved=np.load(warm);star.initial_guess=saved['delta'].copy()
                y=star.base+saved['delta'];logA=-2*saved['psi']**2
                star.material_cache.update({e.material_key(row):raw for row,raw in zip(zip(y[:,0]-3*logA,y[:,1]-logA,y[:,5:]),saved['raw'])})
            def log(row):
                rows.append(row)
                with (path/'iterations.jsonl').open('a') as stream:stream.write(json.dumps(row)+'\n')
                print('REGULAR',cells,index,row['iteration'],row['residual_norm'],flush=True)
            try:
                delta,z,value=stage(star,(zero,z0),log)
                np.savez_compressed(path/'endpoint.npz',delta=delta,residual=value,**{k:v for k,v in z.items() if isinstance(v,np.ndarray)})
                result=assess(star,delta,z,value);result.update(iterations=len(rows),elapsed_seconds=time.monotonic()-began)
                e.write(path/'result.json',result);assert result['passed'],result
                results.append(result);endpoints.append(delta)
            except Exception as error:
                if hasattr(star,'last_state'):
                    np.savez_compressed(path/'failed-iterate.npz',delta=star.last_delta,**{k:v for k,v in star.last_state.items() if isinstance(v,np.ndarray)})
                e.write(folder/'failure.json',dict(classification='Counterexample candidate',path=index,error=repr(error),seconds=time.monotonic()-began))
                raise
    result=dict(classification='Counterexample candidate',passed=True,paths=results,seconds=time.monotonic()-began,
        maximum_driven_minus_undriven=np.max(abs(endpoints[1]-endpoints[0]),axis=0).astype(float).tolist(),
        full_duration_evolution=False,time_convergence_verified=False,observational_closure=False)
    e.write(folder/'result.json',result)
    e.write(folder/'manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():e.digest(p) for p in folder.rglob('*') if p.is_file() and p.name!='manifest.json'}))
    print('PASS',cells,json.dumps(result),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','run'])
    parser.add_argument('--cells',type=int,choices=[16,5735],default=16)
    args=parser.parse_args();globals()[args.action](args.cells)
