"""Exterior-matched charge and a fixed-inventory thermal derivative.

Counterexample candidate: constrained snapshots, not hydrostatic equilibria or
an executed orbit. The leading derivative keeps every shell baryon and radius
fixed. Scalar-vacuum matching reuses the tested Just map without a far cutoff.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import argparse
import json
import time
import mpmath as mp
import numpy as np
import sympy as sp
import def_resolved_scalar_pulse as s
import gr_scalar_nonlinear_exterior as exterior

e,ld=s.e,s.ld
OUT=e.g.OUT/'def-normalized-charge'


def tri(diagonal,off,rhs):
    d=diagonal.copy();b=rhs.copy()
    for i in range(1,len(d)):
        factor=off[i-1]/d[i-1];d[i]-=factor*off[i-1];b[i]-=factor*b[i-1]
    b[-1]/=d[-1]
    for i in range(len(d)-2,-1,-1):b[i]=(b[i]-off[i]*b[i+1])/d[i]
    return b


def lapse(star,E,P,b,bf,Phi):
    faces=star.face_weight*Phi**2/(8*np.pi)
    left=faces[:-1]/2;right=faces[1:]/2;right[-1]*=2
    work=e.GRAV*star.volume*(E+P)/(b*star.r)
    C=1-2*left/np.where(star.rf[:-1]>0,star.rf[:-1],1)-work*(1-star.fraction)
    D=1+2*right/star.rf[1:]+work*star.fraction
    H=np.empty_like(E);lam=E[0]*0+1
    for i in range(star.n-1,-1,-1):H[i]=lam/D[i];lam=H[i]*C[i]
    return H


def scalar(star,E,P,bf,H,mu,G):
    R=star.rf[-1];hg=(H[:-1]+H[1:])/2
    c=star.scalar_area[1:-1]/(4*np.pi*R)*hg*bf[1:-1]/np.diff(star.r)
    ce=bf[-1]/(star.distance[-1]/(H[-1]*star.r[-1])+G)
    potential=e.GRAV*star.volume/R*H*star.beta*(-E+3*P)
    diagonal=np.r_[c,ce]+np.r_[potential[0]*0,c]-potential
    w=tri(diagonal,-c,potential)  # psi=1+w avoids subtraction of nearly equal fields.
    flux=-np.sum(potential*(1+w))
    residual=diagonal*w-potential;residual[1:]-=c*w[:-1];residual[:-1]-=c*w[1:]
    return w,flux,diagonal,c,residual


def snapshot(pool,phi0,theta,beta=-4):
    star,_=s.initialize(pool,beta,ld('.001'),ld(1))
    phi0=ld(phi0);theta=np.asarray(theta,dtype=ld)
    delta=np.zeros_like(star.base);w=np.zeros(star.n,dtype=ld)
    Phi=np.zeros(star.n+1,dtype=ld);q=ld(0);rows=[]
    for iteration in range(30):
        phi=phi0*(1+w);la=star.beta*phi**2/2
        delta[:,1]=theta+la
        raw=star.native_aux(delta,phi);z=star.fluid(delta,phi,raw)
        met=s.prior.metric(star,z['dEm'],z['P'],np.zeros(star.n,dtype=ld),Phi)
        H=lapse(star,z['E'],z['P'],met['b'],met['bf'],Phi)
        mu=met['mf'][-1]/star.rf[-1]
        ext=exterior.exact(mp.mpf(str(mu)),mp.mpf(str(q)))
        G=ld(str(ext[2]))
        wn,flux,diagonal,c,residual=scalar(star,z['E'],z['P'],met['bf'],H,mu,G)
        qn=phi0*flux/met['bf'][-1]
        Phin=np.r_[ld(0),phi0*np.diff(wn)/np.diff(star.r),qn/(H[-1]*star.r[-1])]
        lr=-np.log1p(met['da']/star.reference['a'])
        error=max(float(abs(wn-w).max()),float(abs(lr-delta[:,0]).max()))
        rows.append(dict(iteration=iteration,change=error,native_calls=star.pool.evaluations))
        if error<2e-17:break
        delta[:,0]=lr;w,Phi,q=wn,Phin,qn
    else:raise RuntimeError(('Snapshot iteration cap',rows))
    ext=exterior.exact(mp.mpf(str(mu)),mp.mpf(str(qn)))
    mass=mu+qn*qn*ld(str(ext[0]))
    alpha=flux/met['bf'][-1]*ld(str(ext[1]))/mass
    surface=phi0*(1+wn[-1])+Phin[-1]*star.distance[-1]
    match=abs(surface+qn*ld(str(ext[2]))-phi0)/max(abs(phi0),ld(1e-30))
    baryon=abs(np.exp(delta[:,0])*met['a']/star.reference['a']-1).max()
    field_residual=abs(residual).max()/max(abs(flux),ld(1e-30))
    assert max(float(match),float(baryon),float(field_residual))<2e-13
    z.update(met,H=H,normalized_field=1+wn,normalized_field_increment=wn,
        Phi=Phin,delta=delta,theta=theta)
    return star,z,dict(alpha_over_phi_infinity=float(alpha),phi_infinity=float(phi0),
        surface_flux_q=float(qn),surface_phi=float(surface),mass_geom_cm=float(mass*star.rf[-1]),
        exterior_match_residual=float(match),baryon_relative_residual=float(baryon),
        scalar_relative_residual=float(field_residual),native_calls=star.pool.evaluations,iterations=rows)


def leading(star,base,theta):
    """Analytic native EOS differential and the constrained mass differential."""
    dtype=np.result_type(theta.dtype,base['E'].dtype)
    er=base['E']+base['rho']*base['raw'][:,9];et=base['rho']*base['raw'][:,10]
    dmf=np.zeros(star.n+1,dtype=dtype);de=np.empty(star.n,dtype=dtype)
    for i in range(star.n):
        k=base['a'][i]**2/star.r[i];f=star.fraction[i];gv=e.GRAV*star.volume[i]
        dmf[i+1]=((1-gv*er[i]*k*(1-f))*dmf[i]+gv*et[i]*theta[i])/(1+gv*er[i]*k*f)
        de[i]=(dmf[i+1]-dmf[i])/gv
    dr=-base['a']**2*((1-star.fraction)*dmf[:-1]+star.fraction*dmf[1:])/star.r
    dp=base['P']*(base['raw'][:,5]*dr+base['raw'][:,6]*theta)
    E,P=base['E']+de,base['P']+dp
    b=base['b']-2*((1-star.fraction)*dmf[:-1]+star.fraction*dmf[1:])/star.r
    bf=base['bf']-np.r_[0,2*dmf[1:]/star.rf[1:]]
    H=lapse(star,E,P,b,bf,np.zeros(star.n+1,dtype=dtype))
    mu=(base['mf'][-1]+dmf[-1])/star.rf[-1]
    G=-bf[-1]*np.log1p(-2*mu)/(2*mu)
    _,flux,*_=scalar(star,E,P,bf,H,mu,G)
    return flux/mu


def symbolic():
    B,a,V,rho,U,R,C=sp.symbols('B a V rho U R C',positive=True)
    mass,source,dm,ds=sp.symbols('mass source dm ds')
    h=sp.Symbol('h')
    assert sp.diff((source+h*ds)/(mass+h*dm),h).subs(h,0)==ds/mass-source*dm/mass**2
    assert sp.diff(sp.log(B/(a*V)),a)==-1/a
    x=sp.Symbol('x');M=sp.Matrix([[2,-1],[-1,2]]);p=sp.Matrix([x,2*x])
    K=M-sp.diag(*p);w=K.inv()*p
    assert sp.simplify(K*(sp.ones(2,1)+w)-(M*sp.ones(2,1)))==sp.zeros(2,1)
    return dict(classification='Proven',passed=True,
        scope='Fixed shell B,V implies dlnrho=-dlna; full charge derivative includes mass normalization. Solving for psi-1 is algebraically equivalent. No physical coefficient-error bound.')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    files=[Path(__file__),Path(s.__file__),Path(s.prior.__file__),Path(exterior.__file__),
        s.OUT/'initial.npz',s.OUT/'initial-manifest.json',s.OUT/'companion-benchmark.json',
        e.g.OUT/'gr-resolved-response-milestone-manifest.json',exterior.OUT/'manifest.json']
    e.write(OUT/'plan.json',dict(classification='Counterexample candidate',
        bindings={p.relative_to(s.ROOT).as_posix():e.digest(p) for p in files},
        symbolic=symbolic(),phi0=.001,thermal_steps=[.001,.0003,.0001],
        directions=['uniform','sign_of_leading_charge_gradient'],
        gates=dict(snapshot_residual=2e-13,leading_derivative_relative_agreement=.005,
            derivative_last_two_relative_difference=.005,complex_step_relative_agreement=1e-10),
        budget=dict(workers=4,blas_threads=1,hard_timeout_seconds=180,maximum_snapshot_iterations=30,
            maximum_snapshots=15,automatic_expansion=False),
        claim='Construct finite nonzero scalar background snapshots on the Phase40 full-radius grid with fixed shell baryons, composition, radii and declared Jordan temperature. Compute normalized exterior charge and the leading fixed-inventory thermal derivative including ADM mass normalization. Compare two directions against native finite-background central differences.',
        physical_scope='Constraint-consistent quasistatic scalar readout on constrained Cauchy snapshots, not hydrostatic equilibrium or physical atmosphere. No time evolution, scalar-frequency response, fluid resolvent, orbital stationarity or observed signal. Leading derivative is phi0->0; finite phi0 native comparisons do not certify its remainder.',
        boundary='Just vacuum map has no distant-radius cutoff; the finite-pressure truncated matter surface is a declared wall, not a physical vacuum photosphere.',
        decisions='Use the thermal readout norm to state the minimum normalized fluid amplification needed for the existing 1e-9 target. Without a fluid resolvent bound this is not a no-go. Do not run a long orbit.'))


def run():
    started=time.monotonic();mp.mp.dps=70
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,h in plan['bindings'].items():assert e.digest(s.ROOT/rel)==h,rel
    assert not (OUT/'result.json').exists();results=[]
    with ProcessPoolExecutor(max_workers=4,initializer=s.old.imported.initial.original.worker_init) as pool:
        star,base,zero=snapshot(pool,0,np.zeros(24));results.append(dict(name='zero_background',**zero))
        np.savez_compressed(OUT/'zero.npz',**{k:v for k,v in base.items() if isinstance(v,np.ndarray)})
        gradients=[]
        for step in [1e-20,1e-25]:
            h=[]
            for i in range(star.n):
                t=np.zeros(star.n,dtype=np.clongdouble);t[i]=1j*step
                h.append(float(leading(star,base,t).imag/step))
            gradients.append(np.asarray(h))
        g=gradients[-1];norm=float(abs(g).sum())
        assert np.max(abs(gradients[0]-g))<1e-10*norm
        np.savez_compressed(OUT/'gradient.npz',gradient=g,steps=np.array([1e-20,1e-25]),gradients=gradients)
        _,finite,background=snapshot(pool,.001,np.zeros(24));results.append(dict(name='finite_background',**background))
        np.savez_compressed(OUT/'finite.npz',**{k:v for k,v in finite.items() if isinstance(v,np.ndarray)})
        comparisons=[]
        for name,direction in [('uniform',np.ones(24)),('sign',np.sign(g))]:
            expected=float(g@direction);values=[]
            for step in plan['thermal_steps']:
                pair=[]
                for sign in [-1,1]:
                    _,state,row=snapshot(pool,.001,sign*step*direction)
                    label=f'{name}-{step}-{sign}';results.append(dict(name=label,**row));pair.append(row['alpha_over_phi_infinity'])
                    np.savez_compressed(OUT/f'{label}.npz',**{k:v for k,v in state.items() if isinstance(v,np.ndarray)})
                values.append((pair[1]-pair[0])/(2*step))
                print('NATIVE THERMAL CHARGE',name,step,values[-1],flush=True)
            comparisons.append(dict(direction=name,leading_derivative=expected,native_central_differences=values,
                relative_leading_difference=abs(values[-1]-expected)/max(abs(expected),1e-30),
                last_two_relative_difference=abs(values[-1]-values[-2])/max(abs(expected),1e-30)))
        _,_,control=snapshot(pool,.001,np.zeros(24),beta=0);results.append(dict(name='beta_zero',**control))
        assert control['alpha_over_phi_infinity']==0
    benchmark=json.loads((s.OUT/'companion-benchmark.json').read_text())
    target=benchmark['numerical_target']['charge_over_phi0_residual'];du=benchmark['leading_drive']['maximum_delta_u']
    passed=all(r['relative_leading_difference']<.005 and r['last_two_relative_difference']<.005 for r in comparisons)
    result=dict(classification='Counterexample candidate',passed=passed,rows=results,comparisons=comparisons,
        leading_thermal_gradient_l1=norm,required_temperature_linf_for_target=target/norm,
        required_normalized_fluid_gain=target/(.001**2*du*norm),
        leading_derivative_certificate=False,fluid_gain_bounded=False,physical_benchmark_closed=False,
        finite_background_constructed=True,hydrostatic_equilibrium_solved=False,
        seconds=time.monotonic()-started,total_native_calls=sum(r['native_calls'] for r in results),symbolic=symbolic())
    e.write(OUT/'result.json',result)
    e.write(OUT/'manifest.json',dict(sha256={p.relative_to(s.ROOT).as_posix():e.digest(p) for p in OUT.iterdir() if p.is_file() and p.name not in ['manifest.json','run.log']}))
    print('FINAL',json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run'])
    globals()[p.parse_args().action]()
