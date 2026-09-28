"""Counterexample candidate: solve switching face velocities together.

Eliminate the five-field sparse linear system and solve the small piecewise
linear donor problem. This constructs a guess only; native gates decide it.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json

import numpy as np
from scipy.optimize import root
import gr_step42_analysis as test

full,m,e,ld=test.full,test.m,test.e,test.ld
OUT=m.BASE/'step42-complementarity'


def correction(star,delta,z,previous,older,h,coefficients):
    value,z=m.residual(star,delta,previous,older,h,coefficients)
    matrix=m.jacobian(m.tangent(star,delta,z),delta,previous,older,h,coefficients)
    factor=m.splu(matrix)
    x0=-factor.solve(np.asarray(value[:,:5]/m.SCALE,float).ravel())
    ratio=np.asarray(z['N']/z['a'],float)
    face0=star.faces(z['N']*z['v']/z['a'],odd=True).astype(float)
    velocity_change=ratio*x0.reshape(star.n,5)[:,2]*float(m.SCALE[2])
    face_trial=face0+star.faces(velocity_change,odd=True).astype(float)
    selected=np.flatnonzero((face0>=0)!=(face_trial>=0))
    reports=[]
    for expansion in range(4):
        assert len(selected)>0 and np.all((selected>0)&(selected<star.n))
        # For either anchor sign s, native minus fixed-donor flux is
        # (U_left-U_right)*max(-s*v_face,0). Each face touches two cells.
        U=np.zeros((5*star.n,len(selected)))
        for j,face in enumerate(selected):
            scale=float(h*e.C*4*np.pi*star.rf[face]**2)
            for field,(quantity,normalization) in enumerate([(z['B'],star.B0),(z['E'],star.heat0),(z['AS'],star.momentum_scale)]):
                jump=float(quantity[face-1]-quantity[face])
                U[5*(face-1)+field,j]=scale*jump/float(star.volume[face-1]*normalization[face-1]*m.SCALE[field])
                U[5*face+field,j]=-scale*jump/float(star.volume[face]*normalization[face]*m.SCALE[field])
        inverse=factor.solve(U)
        velocity=inverse.reshape(star.n,5,-1)[:,2,:]*float(m.SCALE[2])*ratio[:,None]
        t=np.asarray((star.rf[1:-1]-star.r[:-1])/np.diff(star.r),float)
        projected=(velocity[:-1]+t[:,None]*(velocity[1:]-velocity[:-1]))[selected-1]
        anchor_sign=np.where(face0[selected]>=0,1.,-1.)
        scales=np.maximum(np.maximum(abs(face0[selected]),abs(face_trial[selected]-face0[selected])),1e-30)
        K=projected*scales[None,:]/scales[:,None]
        a=face_trial[selected]/scales
        def fun(w):
            return w-a+K@np.maximum(-anchor_sign*w,0)
        def jac(w):
            return np.eye(len(w))-K*(anchor_sign*((-anchor_sign*w)>0))[None,:]
        answer=root(fun,face0[selected]/scales,jac=jac,method='hybr',options={'xtol':1e-10,'maxfev':2000})
        error=float(np.max(abs(fun(answer.x))))
        x=x0-inverse@(scales*np.maximum(-anchor_sign*answer.x,0))
        candidate=delta.copy()
        candidate[:,:5]+=x.reshape(star.n,5).astype(ld)*m.SCALE
        candidate[:,5:]=m.species(star,candidate,previous,older,h,coefficients,z)
        native,state=m.residual(star,candidate,previous,older,h,coefficients)
        donor=star.faces(state['N']*state['v']/state['a'],odd=True)
        changed=np.flatnonzero((donor>=0)!=(face0>=0))
        missing=np.setdiff1d(changed,selected)
        score=float(np.max(abs(native)/m.ATOL))
        reports.append(dict(expansion=expansion,faces=len(selected),reduced_residual=error,
            solver_success=bool(answer.success),evaluations=answer.nfev,native_score=score,new_faces=len(missing)))
        print('SIMULTANEOUS DONOR SOLVE',json.dumps(reports[-1]),flush=True)
        if len(missing):
            selected=np.union1d(selected,missing)
            continue
        return candidate,native,state,reports
    raise RuntimeError(('Unresolved donor set',reports))


def check():
    assert not OUT.exists()
    # Exact scalar upwind identity used in the reduced equation.
    for anchor in (-1,1):
        for velocity in (-2.,0.,3.):
            native=velocity*(5 if velocity>=0 else 2)
            fixed=velocity*(5 if anchor>0 else 2)
            assert native-fixed==(5-2)*max(-anchor*velocity,0)
    plan=full.bindings();OUT.mkdir()
    e.write(OUT/'plan.json',dict(classification='Counterexample candidate',source_sha256=e.digest(Path(__file__)),
        failed_plan_sha256=e.digest(full.OUT/'plan.json'),rule='Solve the piecewise linear donor correction jointly; native equations and tolerances decide acceptance.'))
    with ProcessPoolExecutor(max_workers=15,initializer=e.worker_init) as pool:
        star=m.initialize(pool);times=full.time_nodes(plan,1)
        older,previous=[full.restored(star,full.OUT/'path-1',k,times[k]) for k in (40,41)]
        delta,z=test.failed_state(star)
        h=times[42]-times[41];coefficients=m.weights(h,times[41]-times[40])
        with np.load(test.OUT/'lagged_composition.npz') as cp:
            star.composition_data=dict(coefficients=cp['composition_coefficients'],inverse=cp['composition_inverse'])
        star.composition_context=(previous,older,h,coefficients)
        rows=[]
        for iteration in range(14,24):
            delta,value,z,trials=correction(star,delta,z,previous,older,h,coefficients)
            score=float(np.max(abs(value)/m.ATOL))
            rows.append(dict(iteration=iteration,residual_norm=score,maximum_absolute_residual=np.max(abs(value),axis=0).astype(float).tolist(),donor_solve=trials))
            e.write(OUT/'iterations.json',rows)
            if score<=1:
                break
        if score>1:
            e.write(OUT/'failure.json',dict(classification='Counterexample candidate',native_score=score))
            raise AssertionError(score)
        with np.load(full.OUT/'path-1/step-0041.npz') as last,np.load(full.OUT/'path-1/step-0040.npz') as before:
            c0,c1,c2=coefficients;area=4*np.pi*star.rf**2
            energy=(-c1*last['integrated_energy_flux']-c2*before['integrated_energy_flux']+h*e.C*area*z['fluxes'][1])/c0
            baryon=(-c1*last['integrated_baryon_flux']-c2*before['integrated_baryon_flux']+h*e.C*area*z['fluxes'][0])/c0
        budget,ed,bd=m.budget(star,z,energy,baryon);cone=m.parent.cones(z)
        assert m.budget_passed(budget,plan) and cone['sampled_cone_inside_light_cone']
        m.save_state(OUT/'step-0042.npz',delta,z,times[42],energy,baryon,ed,bd)
        e.write(OUT/'result.json',dict(classification='Counterexample candidate',passed=True,native_score=score,iterations=rows,budget=budget,cone=cone,full_duration_completed=False))
        print('NATIVE COMPLEMENTARITY STEP PASSED',score,flush=True)
    e.write(OUT/'manifest.json',dict(sha256={p.relative_to(e.ROOT).as_posix():e.digest(p) for p in [Path(__file__),test.OUT/'manifest.json',*OUT.iterdir()] if p.is_file()}))


if __name__=='__main__':
    check()
