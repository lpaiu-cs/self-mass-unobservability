"""Counterexample candidate: native molecular matter + spherical DEF wave.

One backward-Euler matter / midpoint scalar step, with reciprocal discrete
work. The finite boundary wave is a control, not a companion-matched drive.
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
from scipy.linalg import solve_banded
from scipy.sparse.linalg import splu

import gr_molecular_compatible_step as imported
import gr_two_carrier_evolution as two

m, e, ld = imported.m, imported.e, imported.ld
OUT = e.g.OUT/'def-spherical-coupled'
BETA = ld(-4)
ATOL, SCALE = m.ATOL, m.SCALE


def allocation(face):
    result = (face[:-1]+face[1:])/2
    result[-1] += face[-1]/2
    return result


def symbolic():
    u0,u1,b0,b1 = sp.symbols('u0 u1 b0 b1')
    assert sp.expand(u1*u1*b1-u0*u0*b0
        -(b1+b0)/2*(u1-u0)*(u1+u0)
        -(u1*u1+u0*u0)/2*(b1-b0)) == 0
    p,q,b,bf,H,Hf = sp.symbols('p q b bf H Hf')
    # One interior face's kinetic + half gradient energy + outward flux.
    face = p*b*Hf*bf + bf*(H*q-H*p)/2-bf*(H*p+H*q)/2
    assert sp.expand(face-p*bf*H*(b-1)-p*b*bf*(Hf-H)) == 0
    # General cell/neighbor k values, retaining the cell's b in k=H*b.
    k,kq = sp.symbols('k kq')
    face = p*b*Hf*bf+bf*(kq*q-k*p)/2-bf*(k*p+kq*q)/2
    assert sp.expand(face.subs(k,H*b)-p*b*bf*(Hf-H)) == 0
    return dict(classification='Proven',passed=True,
        scope='Finite product identity and shared-face scalar work identity for the declared staggered energy quadrature; no continuum error bound.')


def selfcheck():
    symbolic()
    s = SimpleNamespace(n=4,r=np.array([1,3,5,7],dtype=ld)*ld('1e6'),
        rf=np.array([0,2,4,6,8],dtype=ld)*ld('1e6'),h=ld('1e-5'),amplitude=ld('1e-6'))
    s.volume=4*np.pi/3*np.diff(s.rf**3);s.area=4*np.pi*s.rf**2
    s.partial_volume=4*np.pi/3*(s.r**3-s.rf[:-1]**3)
    s.distance=np.r_[ld(1),np.diff(s.r),s.rf[-1]-s.r[-1]];s.face_weight=s.area*s.distance
    E=np.full(s.n,ld('1e25'));R=E/100
    mf=np.r_[ld(0),np.cumsum(e.GRAV*E*s.volume)];mass=mf[:-1]+e.GRAV*E*s.partial_volume
    s.b0=1-2*mass/s.r;s.bf0=np.r_[ld(1),1-2*mf[1:]/s.rf[1:]]
    s.reference=dict(E=E,m=mass,mf=mf,a=1/np.sqrt(s.b0))
    s.divergence=lambda f:np.diff(s.area*f)/s.volume
    zero=np.zeros(s.n,dtype=ld);zf=np.zeros(s.n+1,dtype=ld)
    old=metric(s,zero,R,zero,zf);s.scalar_reference=dict(old,trace=-E+3*R)
    z=old
    for _ in range(12):
        psi,Pi,Phi=wave(s,z,s.scalar_reference['trace'])
        z=metric(s,zero,R,Pi,Phi)
    work=scalar_work(s,z,s.scalar_reference['trace'],psi,Pi,Phi)
    assert work['scalar_residual']<2e-17,work['scalar_residual']
    assert work['scalar_relative_balance']<2e-15,work['scalar_relative_balance']
    assert np.max(abs(z['dEt']-z['Ephi'])/(1+z['Ephi']))<ld('1e-17')
    return dict(classification='Counterexample candidate',passed=True,
        scalar_residual=work['scalar_residual'],scalar_relative_balance=work['scalar_relative_balance'])


def metric(star, dEm, Rm, Pi, Phi):
    """Solve the linear-in-mass scalar energy constraint by a radial scan."""
    zero = star.reference
    kin = Pi*Pi/(8*np.pi*e.GRAV)
    grad = star.face_weight*Phi*Phi/(8*np.pi*e.GRAV)
    left, right = grad[:-1]/(2*star.volume), grad[1:]/(2*star.volume)
    right[-1] *= 2
    dEt = np.empty(star.n,dtype=ld); dmf = np.zeros(star.n+1,dtype=ld)
    for i in range(star.n):
        incoming = dmf[i]
        bl = star.bf0[i] if i == 0 else star.bf0[i]-2*incoming/star.rf[i]
        numerator = (dEm[i]+kin[i]*(star.b0[i]-2*incoming/star.r[i])
            +left[i]*bl+right[i]*(star.bf0[i+1]-2*incoming/star.rf[i+1]))
        denominator = 1+2*e.GRAV*(kin[i]*star.partial_volume[i]/star.r[i]
            +right[i]*star.volume[i]/star.rf[i+1])
        dEt[i] = numerator/denominator
        dmf[i+1] = incoming+e.GRAV*star.volume[i]*dEt[i]
    dm = dmf[:-1]+e.GRAV*star.partial_volume*dEt
    db = -2*dm/star.r
    dbf = np.r_[ld(0),-2*dmf[1:]/star.rf[1:]]
    b,bf = star.b0+db,star.bf0+dbf
    assert np.min(b)>0 and np.min(bf)>0
    a = 1/np.sqrt(b)
    factor = np.sqrt(1+db/star.b0)
    da = zero['a']*(-db/star.b0)/(factor*(1+factor))
    Ephi = kin*b+allocation(grad*bf)/star.volume
    mass,mf = zero['m']+dm,zero['mf']+dmf
    Et,Rt = zero['E']+dEm+Ephi,Rm+Ephi
    Hprime = 4*np.pi*e.GRAV*star.r*a*a*(Et+Rt)
    logHf = np.r_[-np.cumsum((Hprime*np.diff(star.rf))[::-1])[::-1],ld(0)]
    logH = logHf[1:]-Hprime*(star.rf[1:]-star.r)
    N = np.exp(logH)/a
    return dict(a=a,da=da,b=b,bf=bf,db=db,dbf=dbf,m=mass,mf=mf,N=N,
        logH=logH,logHf=logHf,Etot=Et,Ephi=Ephi,dEt=dEt,
        nur=a*a*(mass/star.r**2+4*np.pi*e.GRAV*star.r*Rt),
        ar=a*a*(4*np.pi*e.GRAV*star.r*Et-mass/star.r**2))


def scalar_coefficients(star,z,trace):
    old = star.scalar_reference
    b,bf = (z['b']+old['b'])/2,(z['bf']+old['bf'])/2
    logH,logHf = (z['logH']+old['logH'])/2,(z['logHf']+old['logHf'])/2
    H,Hf = np.exp(logH),np.exp(logHf)
    k,kf = H*b,Hf*bf
    force = star.h**2*np.pi*e.GRAV*e.C**2*k*H*BETA*(trace+old['trace'])/2
    conductance = star.area*kf/star.distance
    conductance[0] = 0
    factor = star.h**2*e.C**2*k/(4*star.volume)
    lower = -factor[1:]*conductance[1:-1]
    upper = -factor[:-1]*conductance[1:-1]
    diagonal = 1+factor*(conductance[:-1]+conductance[1:])-force
    rhs = np.zeros(star.n,dtype=ld)
    rhs[-1] = factor[-1]*conductance[-1]*star.amplitude/2
    return lower,diagonal,upper,rhs,k,kf,b,bf,logH,logHf


def wave(star,z,trace):
    lower,diagonal,upper,rhs,k,*_ = scalar_coefficients(star,z,trace)
    band = np.zeros((3,star.n));band[1]=diagonal
    band[0,1:],band[2,:-1] = upper,lower
    mid = solve_banded((1,1),band,np.asarray(rhs,float)).astype(ld)
    # Binary64 factorization, extended-precision residual refinement.
    for _ in range(3):
        defect = rhs-diagonal*mid
        defect[1:] -= lower*mid[:-1]; defect[:-1] -= upper*mid[1:]
        mid += solve_banded((1,1),band,np.asarray(defect,float)).astype(ld)
    psi = 2*mid
    Pi = 4*mid/(star.h*e.C*k)
    Phi = np.r_[ld(0),np.diff(psi)/np.diff(star.r),(star.amplitude-psi[-1])/star.distance[-1]]
    return psi,Pi,Phi


def scalar_work(star,z,trace,psi,Pi,Phi):
    coeff = scalar_coefficients(star,z,trace)
    lower,diagonal,upper,rhs,k,kf,b,bf,logH,logHf = coeff
    mid,pm,fm = psi/2,Pi/2,Phi/2
    defect = diagonal*mid-rhs
    defect[1:] += lower*mid[:-1];defect[:-1] += upper*mid[1:]
    scale = max(abs(star.amplitude),ld('1e-100'))
    wave_error = np.max(abs(defect))/scale
    definition_error = np.max(abs(psi-star.h*e.C*k*pm))/scale
    velocity = k*pm
    face_velocity = np.r_[ld(0),(velocity[:-1]+velocity[1:])/2,star.amplitude/(star.h*e.C)]
    flux = -bf*fm*face_velocity/(4*np.pi*e.GRAV)
    direct = BETA*mid*(trace+star.scalar_reference['trace'])/2*psi/star.h
    temporal = (Pi*Pi/2*z['db']+allocation(star.face_weight*Phi*Phi/2*z['dbf'])/star.volume)/(8*np.pi*e.GRAV*star.h)
    # Stable differences of H=N*a avoid subtracting nearly equal unit values.
    dHleft = np.exp(logH)*np.expm1(logHf[:-1]-logH)
    dHright = np.exp(logH)*np.expm1(logHf[1:]-logH)
    spatial = e.C*pm*b/(4*np.pi*e.GRAV*star.volume)*(star.area[1:]*fm[1:]*bf[1:]*dHright
        -star.area[:-1]*fm[:-1]*bf[:-1]*dHleft)
    work = direct+temporal+spatial
    balance = z['Ephi']+star.h*e.C*star.divergence(flux)-star.h*work
    denom = np.maximum(z['Ephi']+abs(star.h*e.C*star.divergence(flux))+abs(star.h*work),ld('1e-100'))
    return dict(scalar_flux=flux,scalar_work=work,direct_work=direct,geometric_work=temporal+spatial,
        scalar_balance=balance,scalar_relative_balance=float(np.max(abs(balance)/denom)),
        scalar_residual=float(max(wave_error,definition_error)))


class CoupledStar(m.CompatibleStar):
    def native_aux(self,delta,psi):
        y = self.base+delta; logA = -2*psi*psi
        rows = list(zip(y[:,0]-3*logA,y[:,1]-logA,y[:,5:]))
        keys = [e.material_key(row) for row in rows]
        missing = {key:row for key,row in zip(keys,rows) if key not in self.material_cache}
        self.material_cache.update(zip(missing,self.pool.map(two.material,missing.values(),chunksize=8)))
        raw = np.asarray([self.material_cache[key] for key in keys],dtype=ld)
        if len(self.material_cache)>8*self.n:
            self.material_cache = dict(zip(keys,raw))
        return raw

    def fluid(self,delta,psi,raw):
        zero = self.reference
        logA = -2*psi*psi; A = np.exp(logA)
        aux = raw.copy()
        aux[:,1] *= np.exp(4*logA)
        aux[:,[2,9,10]] *= A[:,None]
        aux[:,[21,24,27]] *= np.exp(-2*logA)[:,None]
        rho,T = zero['rho']*np.exp(delta[:,0]),zero['T']*np.exp(delta[:,1])
        drestJ = (delta[:,5:]/e.g.c.A)@e.g.c.W*e.C**2
        drest = drestJ+np.expm1(logA)*(zero['rest']+drestJ)
        rest,u,P = zero['rest']+drest,aux[:,2],aux[:,1]
        eps0 = zero['rho']*(zero['rest']+zero['u'])
        deps = zero['rho']*((zero['rest']+zero['u'])*np.expm1(delta[:,0])
            +np.exp(delta[:,0])*(drest+u-zero['u']))
        v = delta[:,2]; Q = (self.base[:,3]+delta[:,3])*self.qscale
        W = 1/np.sqrt(1-v*v)
        dE = (deps+(P+eps0)*v*v+2*Q*v)/(1-v*v)
        eps = rho*(rest+u);w=eps+P
        S = (w*v+Q*(1+v*v))*W*W
        R = (eps*v*v+P+2*Q*v)*W*W
        Qrad = (self.base[:,4]+delta[:,4])*self.qscale
        transport = 16*ld('5.670400e-5')*T**3/(3*rho)
        Krad,Kcond = transport/aux[:,24],transport/aux[:,27]
        return dict(aux=aux,raw=raw,rho=rho,T=T,P=P,u=u,rest=rest,w=w,v=v,Q=Q,W=W,D=rho*W,
            E=eps0+dE,dEm=dE,S=S,R=R,Qrad=Qrad,Qcond=Q-Qrad,Krad=Krad,Kcond=Kcond,K=Krad+Kcond,
            tau_rad=1/(e.C*rho*aux[:,24]),tau_cond=e.TAU/A,trace=-eps+3*P,psi=psi,logA=logA,
            beta_delta=np.column_stack([-logA-np.log(Kcond/zero['Kcond'])-2*delta[:,1],-5*delta[:,1]]))

    def evaluate(self,delta):
        assert np.max(abs(delta[:,2]))<.2
        X = self.base[:,5:]+delta[:,5:]
        assert X.min()>=-1e-15 and np.max(abs(X.sum(1)-1))<1e-10
        if hasattr(self,'linearization'):
            z = self.fluid(delta,self.linearization['psi'],self.tangent_aux(delta))
            for key in self.metric_keys:
                z[key] = self.linearization[key]
            for key in ['scalar_flux','scalar_work','direct_work','geometric_work','scalar_balance',
                        'scalar_relative_balance','scalar_residual','Pi','Phi']:
                z[key] = self.linearization[key]
        else:
            psi,Pi,Phi = [v.copy() for v in self.field_guess]
            for iteration in range(12):
                z = self.fluid(delta,psi,self.native_aux(delta,psi))
                z.update(metric(self,z['dEm'],z['R'],Pi,Phi))
                psin,pin,phin = wave(self,z,z['trace'])
                error = np.max(abs(psin-psi))/max(abs(self.amplitude),ld('1e-100'))
                psi,Pi,Phi = psin,pin,phin
                if error<ld('2e-19'):
                    break
            else:
                raise RuntimeError(('Scalar-metric fixed point',float(error)))
            z = self.fluid(delta,psi,self.native_aux(delta,psi))
            z.update(metric(self,z['dEm'],z['R'],Pi,Phi))
            z.update(scalar_work(self,z,z['trace'],psi,Pi,Phi),Pi=Pi,Phi=Phi)
            assert z['scalar_residual']<2e-17,z['scalar_residual']
            self.field_guess = psi,Pi,Phi
        self.finish_moments(delta,z,z['dEm'])
        if hasattr(self,'linearization'):
            anchor=self.linearization
            m.finish_fluxes(self,z,self.faces(anchor['N']*anchor['v']/anchor['a'],odd=True))
        mass_rate = -e.GRAV*e.C*self.area*(z['fluxes'][1]+z['scalar_flux'])
        mdot = mass_rate[:-1]+self.fraction*np.diff(mass_rate)
        z['at'] = z['a']**3*mdot/self.r
        self.current = z
        return z

    def gradient(self,value,odd=False):
        z = getattr(self,'current',None)
        if z is not None and (value is z['Q'] or value is z['Qrad']):
            return z['a']/z['N']*self.divergence(self.faces(z['N']*value/z['a'],odd=True))-value*(z['nur']-z['ar']+2/self.r)
        return e.Star.gradient(self,value,odd)

    def tangent_aux(self,delta):
        delta = delta.copy()
        last,older,h,coefficients = self.composition_context
        projected = m.species(self,delta,last,older,h,coefficients,self.linearization)
        delta[:,5:] = self.anchor[:,5:]+projected-self.projected_anchor
        self.projected_delta = delta
        change = delta-self.anchor
        raw = self.linearization['raw'].copy()
        factors = np.einsum('nki,ni->nk',self.composition_data['inverse'],change[:,5:])
        extra = np.einsum('nk,nkj->nj',factors,self.composition_data['coefficients'])
        raw[:,1] *= np.exp(raw[:,5]*change[:,0]+raw[:,6]*change[:,1]+extra[:,0])
        raw[:,2] += raw[:,9]*change[:,0]+raw[:,10]*change[:,1]+extra[:,1]
        for q,j in [(2,24),(3,27)]:
            raw[:,j] *= np.exp(raw[:,j+1]*change[:,0]+raw[:,j+2]*change[:,1]+extra[:,q])
        raw[:,21] = raw[:,24]*raw[:,27]/(raw[:,24]+raw[:,27])
        return raw


# Reuse all 31 finite matter equations. Only the actual proper conduction time
# changes under the frame map; scalar exchange is added explicitly below.
_source = inspect.getsource(m.prior.residual).replace('def residual(', 'def matter_residual(')
_source = _source.replace("z['Kcond'],e.TAU,0", "z['Kcond'],z['tau_cond'],0")
exec(compile(_source,__file__,'exec'),dict(vars(m.prior)),_namespace := {})
matter_residual = _namespace['matter_residual']


def residual(star,delta,previous,older,h,coefficients):
    value,z = matter_residual(star,delta,previous,older,h,coefficients)
    value[:,1] += h*z['scalar_work']/star.heat0
    force = -z['at']*z['S']-e.C*(z['N']*z['nur']*z['E']-star.equilibrium_gravity)
    force += e.C*(z['N']*z['P']-star.equilibrium_pressure)*star.area_difference_over_volume
    phi_gradient = (z['Phi'][:-1]+z['Phi'][1:])/2
    force += e.C*z['N']*BETA*z['psi']*z['trace']*phi_gradient
    value[:,2] = (z['dU'][:,2]+h*e.C*star.divergence(z['fluxes'][2]-star.equilibrium_pressure_flux)-h*force)/star.momentum_scale
    return value,z


class Tangent(CoupledStar):
    def evaluate(self,delta):
        raw = self.tangent_aux(delta)
        projected = self.projected_delta
        # CoupledStar consumes this precomputed raw only at its tangent branch.
        self.tangent_aux = lambda _: raw
        try:
            return super().evaluate(projected)
        finally:
            del self.tangent_aux


def tangent(star,delta,z):
    model = Tangent.__new__(Tangent);model.__dict__=dict(star.__dict__)
    model.anchor,model.linearization = delta.copy(),z
    model.projected_anchor = m.species(star,delta,*star.composition_context,z)
    return model


jacobian = FunctionType(m.prior.jacobian.__code__,dict(vars(m.prior),residual=residual))


def initialize(pool,amplitude,h):
    star = imported.initialize(pool)
    star.__class__ = CoupledStar
    star.h,star.amplitude = h,amplitude
    star.area = 4*np.pi*star.rf**2
    star.distance = np.r_[ld(1),np.diff(star.r),star.rf[-1]-star.r[-1]]
    star.face_weight = star.area*star.distance
    zero = star.reference
    star.b0 = 1/zero['a']**2
    star.bf0 = np.r_[ld(1),1-2*zero['mf'][1:]/star.rf[1:]]
    field = np.zeros(star.n,dtype=ld);faces=np.zeros(star.n+1,dtype=ld)
    met = metric(star,field,zero['R'],field,faces)
    star.metric_keys = list(met)
    zero.update(met)
    star.scalar_reference = dict(met,trace=-zero['rho']*(zero['rest']+zero['u'])+3*zero['P'])
    star.field_guess = field.copy(),field.copy(),faces.copy()
    m.attach_equilibrium(star)
    return star


def composition_response(star,delta,z):
    # Reuse the native two-neighbor probes in Jordan variables.
    view = m.CompatibleStar.__new__(m.CompatibleStar);view.__dict__=dict(star.__dict__)
    view.base = star.base.copy()
    view.base[:,0] -= 3*z['logA'];view.base[:,1] -= z['logA']
    return m.composition_response(view,delta,dict(z,aux=z['raw']))


def stage(star,previous,log):
    delta = previous[0].copy();h=star.h;coefficients=m.weights(h,None)
    value,z = residual(star,delta,previous,previous,h,coefficients)
    star.composition_context = previous,previous,h,coefficients
    star.composition_data = composition_response(star,delta,z)
    for iteration in range(24):
        norm = float(np.max(abs(value)/ATOL))
        log(dict(iteration=iteration,residual_norm=norm,scalar_residual=z['scalar_residual'],
            maximum_absolute_residual=np.max(abs(value),axis=0).astype(float).tolist()))
        star.last_delta,star.last_state = delta.copy(),z
        if norm<=1:
            return delta,z,value
        if iteration==23:
            break
        factor = splu(jacobian(tangent(star,delta,z),delta,previous,previous,h,coefficients))
        correction = -factor.solve(np.asarray(value[:,:5]/SCALE,float).ravel()).reshape(star.n,5).astype(ld)*SCALE
        fraction = min(1.,.1/max(float(abs(correction[:,:2]).max()),1e-300),.01/max(float(abs(correction[:,2]).max()),1e-300))
        for _ in range(8):
            candidate = delta.copy();candidate[:,:5] += fraction*correction
            candidate[:,5:] = m.species(star,candidate,previous,previous,h,coefficients,z)
            trial,state = residual(star,candidate,previous,previous,h,coefficients)
            score = float(np.max(abs(trial)/ATOL))
            if score<=1 or score<norm*(1-1e-4*fraction):
                delta,value,z = candidate,trial,state
                break
            fraction /= 2
        else:
            raise RuntimeError(('Native joint line search failed',norm,score))
    raise RuntimeError(('Native joint iteration budget',norm))


def prepare():
    assert not OUT.exists()
    check=selfcheck()
    files = [Path(__file__),Path(imported.__file__),Path(imported.initial.__file__),Path(m.__file__),
        Path(m.prior.__file__),Path(m.method.__file__),Path(two.__file__),
        imported.initial.OUT/'initial-manifest.json',imported.initial.OUT/'initial.npz',imported.OUT/'manifest.json']
    template = json.loads((imported.OUT/'plan.json').read_text())
    OUT.mkdir()
    e.write(OUT/'plan.json',dict(classification='Counterexample candidate',
        bindings={p.relative_to(e.ROOT).as_posix():e.digest(p) for p in files},
        cells=5735,h_seconds=str(e.TAU/128),boundary_amplitudes=['0','0.000001'],
        method='Same molecular initial Cauchy state, phi=Pi=Phi=0. One BE native matter step coupled to a midpoint massless DEF scalar step. Scalar and radial Einstein constraints solved inside every native matter residual. Staggered scalar gradient energy and reciprocal finite work. New piecewise-constant quadrature of log(N*a); fixed zero-scalar hydrostatic momentum projection. Reflecting matter wall, prescribed linear Dirichlet scalar ramp at the wall.',
        beta=-4,frame='A=exp(-2*phi^2); native Jordan EOS and opacity; Einstein rho,T,P,u,kappa,tau transformations.',
        decision='Can the full 5735-cell native matter/heat/26-species and scalar wave evolve simultaneously, pass original matter gates and reciprocal field budgets, and produce a driven-minus-undriven endpoint?',
        nonlinear_absolute_tolerances=ATOL.astype(float).tolist(),maximum_iterations=24,
        scalar_residual_tolerance=2e-17,scalar_relative_budget_tolerance=2e-15,
        finite_conservation_gates=template['finite_conservation_gates'],
        budget=dict(workers=8,blas_threads=1,hard_timeout_seconds=1200,maximum_paired_runs=1,
            expected_wall_seconds=[350,1000],basis='Previous native molecular single step 156.5s, including 12905 material calls; two paths plus nested scalar/metric scans. New scalar overhead unmeasured, capped by 1200s. No automatic amplitude/grid/time expansion.'),
        symbolic=symbolic(),finite_operator_check=check,limitations=['Boundary wave is a numerical control, not companion-matched orbital forcing.',
            'Single step; no time convergence or continuum error certificate.',
            'Discrete reciprocal energy work is a spatial/time discretization; no independent Einstein angular/momentum constraint convergence proof.',
            'No physical EOS/atmosphere or observational closure.']))
    print('PREPARED bounded simultaneous paired spherical evolution',flush=True)


def run():
    plan = json.loads((OUT/'plan.json').read_text())
    for rel,sha in plan['bindings'].items():assert e.digest(e.ROOT/rel)==sha,rel
    assert not (OUT/'result.json').exists() and not (OUT/'failure.json').exists()
    start=time.monotonic();results=[];endpoints=[]
    with ProcessPoolExecutor(max_workers=8,initializer=imported.initial.original.worker_init) as pool:
        for index,amplitude in enumerate(plan['boundary_amplitudes']):
            folder=OUT/f'path-{index}';folder.mkdir()
            star=initialize(pool,ld(amplitude),ld(plan['h_seconds']))
            zero=np.zeros_like(star.base)
            # The initial state has the initial boundary value zero.
            star.amplitude=ld(0);z0=star.evaluate(zero);star.amplitude=ld(amplitude)
            records=[]
            def log(row):
                records.append(row)
                with (folder/'iterations.jsonl').open('a') as stream:stream.write(json.dumps(row)+'\n')
                print('JOINT',index,row['iteration'],row['residual_norm'],row['scalar_residual'],flush=True)
            try:
                delta,z,value=stage(star,(zero,z0),log)
                # Reuse the actual local rest principal cone with transformed tau.
                cone_source=inspect.getsource(two.cones).replace('def cones(', 'def frame_cones(').replace('e.TAU',"z['tau_cond']")
                ns=dict(vars(two));exec(compile(cone_source,__file__,'exec'),ns)
                cone=ns['frame_cones'](z)
                energy=star.h*e.C*star.area*z['fluxes'][1]
                baryon=star.h*e.C*star.area*z['fluxes'][0]
                corrected=dict(z,dU=z['dU'].copy())
                corrected['dU'][:,1] += star.h*z['scalar_work']
                budget,ed,bd=m.budget(star,corrected,energy,baryon)
                assert m.budget_passed(budget,plan),(budget,'matter conservation')
                assert cone['sampled_cone_inside_light_cone'],cone
                assert z['scalar_relative_balance']<plan['scalar_relative_budget_tolerance'],z['scalar_relative_balance']
                assert z['scalar_residual']<plan['scalar_residual_tolerance']
                np.savez_compressed(folder/'endpoint.npz',delta=delta,residual=value,**{k:v for k,v in z.items() if isinstance(v,np.ndarray)})
                total=z['dEm']+z['Ephi']+star.h*e.C*star.divergence(z['fluxes'][1]+z['scalar_flux'])
                allowance=star.heat0*ld('1e-9')+ld('2e-15')*(z['Ephi']+abs(star.h*e.C*star.divergence(z['scalar_flux'])))
                assert np.max(abs(total)/allowance)<=1,float(np.max(abs(total)/allowance))
                result=dict(classification='Counterexample candidate',passed=True,amplitude=amplitude,iterations=len(records),
                    native_residual_norm=float(np.max(abs(value)/ATOL)),scalar_residual=z['scalar_residual'],
                    scalar_relative_balance=z['scalar_relative_balance'],total_energy_allowance_ratio=float(np.max(abs(total)/allowance)),
                    budget=budget,cone=cone,native_calls=star.pool.evaluations,
                    scalar_energy_erg=float(np.sum(z['Ephi']*star.volume)),
                    direct_matter_exchange_erg=float(-star.h*np.sum(z['direct_work']*star.volume)),
                    geometric_matter_exchange_erg=float(-star.h*np.sum(z['geometric_work']*star.volume)),
                    max_phi=float(np.max(abs(z['psi']))),elapsed_seconds=time.monotonic()-start)
                e.write(folder/'result.json',result);results.append(result);endpoints.append(delta)
            except Exception as error:
                if hasattr(star,'last_state'):
                    np.savez_compressed(folder/'failed-iterate.npz',delta=star.last_delta,**{k:v for k,v in star.last_state.items() if isinstance(v,np.ndarray)})
                e.write(OUT/'failure.json',dict(classification='Counterexample candidate',path=index,error=repr(error),seconds=time.monotonic()-start))
                raise
    difference=endpoints[1]-endpoints[0]
    result=dict(classification='Counterexample candidate',passed=True,paths=results,
        maximum_driven_minus_undriven=np.max(abs(difference),axis=0).astype(float).tolist(),
        seconds=time.monotonic()-start,simultaneous_native_scalar_matter_metric=True,
        companion_matched=False,time_convergence_verified=False,observational_closure=False)
    e.write(OUT/'result.json',result)
    e.write(OUT/'manifest.json',dict(sha256={p.relative_to(e.ROOT).as_posix():e.digest(p) for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json'}))
    print('PASS simultaneous paired spherical evolution',json.dumps(result),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['symbolic','selfcheck','prepare','run'])
    print(globals()[parser.parse_args().action]())
