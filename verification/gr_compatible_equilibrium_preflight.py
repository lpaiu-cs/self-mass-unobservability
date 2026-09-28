"""Counterexample candidate: mandatory preflight for the redesigned GR run.

Known-failure control, acoustic power, TOV consistency, native nonzero-drive
checks and the fixed two-pulse time selection precede every full-grid pilot.
The original failure and its sources are never overwritten.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from types import SimpleNamespace
import json
import os

import numpy as np
from scipy.linalg import eig, expm, solve, lu_factor, lu_solve
import sympy as sp

import gr_compatible_equilibrium_evolution as model
import gr_time_convergence_cause as control
import gr_acoustic_adjoint_cause as acoustic

ld = model.ld


def register():
    assert not model.BASE.exists()
    model.BASE.mkdir()
    model.e.write(model.BASE/'design.json', dict(classification='Counterexample candidate',
        spatial_discretization='Dual pressure interpolation plus minus one-half the maximum adjacent native acoustic impedance times the velocity jump. The coefficient is fixed, not fitted to measured growth. Unchanged shared energy/baryon/isotope fluxes. Balance only the stationary hydrostatic momentum quadrature using a fixed reference, with equilibrium flux/source differences evaluated before divergence.',
        equilibrium_scope='The imported conserved TOV projection with v=Q=0 is the hydrostatic reference. It is not a thermal equilibrium. Do not subtract any heat, species, energy or baryon residual. Preserve the original nonzero Q and full heat laws in the actual evolution.',
        reference='https://arxiv.org/abs/2108.02960, section 3.2; equilibrium subtraction principle only, not an implementation of that paper or a proof for this two-carrier model.',
        gates=dict(acoustic_maximum_real_per_second=1e-10, acoustic_adjoint_relative=1e-12,
            native_core_growth_over_full_duration=.01, manufactured_TOV_minimum_order=.9,
            derivative_relative=.005, pulse_pair_order=1.5, finest_pulse_relative_error=.01,
            pulse_reference_crosscheck_relative=1e-8, preserved_native_rows=[0,1,3,4]+list(range(5,31))),
        fixed_pulses=['Adiabatic density Gaussian with width r[24], maximum log-density 1e-6.',
                      'Regular radial velocity Gaussian with width r[24], coefficient v/c=1e-8.'],
        selection='Choose the first passing three-level triple among 1,2,4,8,16,32,64,128 for BOTH fixed nonzero pulses. Compare all three fluid fields over all 39 common times and central 16 cells. Stop if no triple passes; never retune thresholds after observing production.',
        native_pilot='Full 5735 cells and 31 native equations, duration 4*tau_cond, 8/16/32 BDF steps; every path from the same original conserved initial state. Apply the unchanged original five-field time gate and all residual/conservation/cone gates.',
        production='Only after preflight and native pilot pass: full original duration with the selected multiplier and 1/2/4 refinements. No old accepted time-state import. Report all common times, regions and endpoints; retain failures.',
        freeze='Engineering implementation is frozen by preflight.json source hashes after passing these fixed gates; production preparation rechecks them.',
        predecessor='v1 dual interpolation alone failed both frozen-composition and advected-composition growth controls. This v2 adds consistent acoustic velocity-jump dissipation; all numerical acceptance thresholds are unchanged.',
        additional_positive_control='Periodic smooth linear acoustic wave, one quarter crossing, 32/64/128/256 cells; compare the exact semidiscrete propagator with continuum translation and require order at least 0.9. The physical pulse must converge as acoustic viscosity vanishes with cell width.',
        linear_control='Native rho/T derivatives and native EOS responses to all 26 advected isotopes through their linearized two-neighbor image; heat perturbations held zero, full initial metric reconstructed. Verify this image against the native isotope-minus-baryon equations. Compare sequential and whole-interval matrix exponentials.',
        original_failed_verdict_sha256=model.e.digest(model.SOURCE/'time-refinement.json')))


def acoustic_checks(star, cells):
    D,G,dual,volume,error=acoustic.operators(star,cells)
    speed=1.16e8
    aux=np.zeros((star.n,30),dtype=ld)
    aux[:,5]=aux[:,10]=1
    aux[:,9]=(speed/model.e.C)**2
    z=dict(rho=np.ones(star.n,dtype=ld),w=np.ones(star.n,dtype=ld),
           N=np.ones(star.n,dtype=ld),P=aux[:,9],aux=aux)
    damping=np.empty_like(D)
    for j in range(cells):
        z['v']=np.zeros(star.n,dtype=ld)
        z['v'][j]=1
        damping[:,j]=-model.e.C*star.divergence(model.acoustic_stress(star,z))[:cells]
    matrices=[np.block([[np.zeros_like(D),-speed*D],[-speed*p,viscosity]])
              for p,viscosity in [(G,np.zeros_like(D)),(dual,np.zeros_like(D)),(dual,damping)]]
    return [eig(a,right=False) for a in matrices],error


def dissipative_symbolic():
    t,p,q,u,v,Z=sp.symbols('t p q u v Z',real=True)
    pressure=t*p+(1-t)*q-Z*(v-u)/2
    work=(p-q)*((1-t)*u+t*v)+u*(pressure-p)+v*(q-pressure)
    assert sp.expand(work-Z*(v-u)**2/2)==0
    return dict(classification='Proven',face_energy_rate='-A*Z*(v_R-v_L)^2/2 <= 0 for Z >= 0',
        scope='Constant-medium linear acoustic energy only. The velocity-jump stress is O(cell width) for smooth fields; no nonlinear GR entropy theorem.')


def acoustic_wave():
    errors=[]
    for n in [32,64,128,256]:
        dx=1/n
        k=2*np.pi
        off=-1j*np.sin(k*dx)/dx
        matrix=np.array([[0,off],[off,-2*np.sin(k*dx/2)**2/dx]])
        actual=expm(.25*matrix)@np.ones(2)
        errors.append(float(np.max(abs(actual-np.exp(-.25j*k)))))
    return dict(classification='Counterexample candidate',cells=[32,64,128,256],
                maximum_field_errors=errors,orders=np.log2(np.array(errors[:-1])/errors[1:]).tolist())


def manufactured_tov():
    errors = []
    epsilon = ld('.3')/(4*np.pi)
    compactness = ld('.2')
    surface = np.sqrt(1-compactness)
    for cells in [32,64,128,256]:
        rf = np.linspace(ld(0),ld(1),cells+1)**ld('1.3')
        r = ((rf[:-1]**3+rf[1:]**3)/2)**(ld(1)/3)
        volume = 4*np.pi/3*np.diff(rf**3)
        area = 4*np.pi*rf**2
        q = np.sqrt(1-compactness*r*r)
        N = (3*surface-q)/2
        P = epsilon*(q-surface)/(3*surface-q)
        nur = compactness*r/(q*(3*surface-q))
        faces = model.pressure_faces(SimpleNamespace(r=r,rf=rf),N*P)
        defect = np.diff(area*faces)/volume-N*P*np.diff(area)/volume+N*nur*epsilon
        errors.append(dict(cells=cells,volume_l1=float(np.sum(abs(defect)*volume)),
                           interior_linf=float(np.max(abs(defect[:-2])))))
    values=np.array([[row['volume_l1'],row['interior_linf']] for row in errors])
    orders=np.log2(values[:-1]/values[1:])
    return dict(classification='Counterexample candidate',rows=errors,orders=orders.tolist(),
        scope='Analytic constant-density Schwarzschild interior, smooth nonuniform radial meshes. The reflecting cell-pressure outer closure is checked in volume L1; interior L-infinity excludes its last two cells. This is consistency, not a global truncation-error bound.')


def native_checks(star):
    model.e.worker_init()
    zero=np.zeros_like(star.base)
    z0=star.evaluate(zero)
    null=zero.copy()
    null[:,3:5]=-star.base[:,3:5]
    zn=star.evaluate(null)
    hydro,_=model.residual(star,null,(null,zn),(null,zn),ld(1),model.weights(ld(1),None))
    assert np.max(abs(hydro[:,:3])) < ld('1e-25')
    same=[0,1,3,4]+list(range(5,31))
    tests=[]
    for cell,column,amplitude in [(0,1,ld('1e-6')),(3142,2,ld('1e-8')),(5704,3,ld('1e-13'))]:
        delta=zero.copy()
        delta[cell,column]=amplitude
        y=star.base+delta
        row=(y[cell,0],y[cell,1],y[cell,5:])
        star.material_cache[model.e.material_key(row)]=model.parent.parent.material(row)
        corrected,z=model.residual(star,delta,(zero,z0),(zero,z0),model.e.TAU,model.weights(model.e.TAU,None))
        ordinary,_=model.prior.residual(star,delta,(zero,z0),(zero,z0),model.e.TAU,model.weights(model.e.TAU,None))
        assert np.array_equal(corrected[:,same],ordinary[:,same])
        assert np.max(abs(corrected[:,2])) > 1e-15
        tests.append(dict(cell=cell,column=column,native_other_rows_exact=True,
            maximum_momentum_residual=float(np.max(abs(corrected[:,2])))))
    actual,_=model.residual(star,zero,(zero,z0),(zero,z0),model.e.TAU,model.weights(model.e.TAU,None))
    assert np.max(abs(actual[:,3:5])) > 1e-17
    return dict(classification='Counterexample candidate',hydrostatic_zero_heat_null_maximum=float(np.max(abs(hydro[:,:3]))),
        probes=tests,actual_thermal_drive_not_removed=True,maximum_initial_heat_residual=float(np.max(abs(actual[:,3:5]))))


class FrozenCompatible(model.CompatibleStar):
    def evaluate(self,delta):
        aux=self.initial_aux.copy()
        extra=np.zeros((self.n,4),dtype=ld)
        extra[:control.NC]=np.einsum('nk,nkj->nj',self.composition_coordinates,self.control_composition)
        aux[:,1]*=np.exp(aux[:,5]*delta[:,0]+aux[:,6]*delta[:,1]+extra[:,0])
        aux[:,2]+=aux[:,9]*delta[:,0]+aux[:,10]*delta[:,1]+extra[:,1]
        for q,j in [(2,24),(3,27)]:
            aux[:,j]*=np.exp(aux[:,j+1]*delta[:,0]+aux[:,j+2]*delta[:,1]+extra[:,q])
        aux[:,21]=aux[:,24]*aux[:,27]/(aux[:,24]+aux[:,27])
        y=self.base+delta
        self.material_cache={model.e.material_key(row):a for row,a in zip(zip(y[:,0],y[:,1],y[:,5:]),aux)}
        return super().evaluate(delta)


def evaluate_linear(x):
    star,n=control.STAR,control.NC
    x=np.asarray(x,dtype=ld).reshape(n,4)
    delta=control.ZERO.copy()
    delta[:n,:3]=x[:,:3]*control.SCALE
    displacement=np.zeros(star.n,dtype=ld)
    displacement[:n]=x[:,3]*control.SCALE[2]
    speed=star.faces(control.Z0['N']*displacement/control.Z0['a'],odd=True)
    area=4*np.pi*star.rf**2
    B=star.B0
    # This displacement integrates the scaled cell velocity. The two factors
    # below are the exact linearized, symmetric-donor species/baryon flux
    # difference. All 26 isotope perturbations are retained in their image.
    left=model.e.C*area[:n]*np.r_[B[0],B[:n-1]]*speed[:n]/(2*star.volume[:n]*B[:n])
    right=-model.e.C*area[1:n+1]*B[1:n+1]*speed[1:n+1]/(2*star.volume[:n]*B[:n])
    star.composition_coordinates=np.column_stack([left,right])
    delta[:n,5:]=np.einsum('nk,nki->ni',star.composition_coordinates,star.control_basis)
    z=star.evaluate(delta)
    spatial,_=model.residual(star,delta,(delta,z),(delta,z),ld(1),model.weights(ld(1),None))
    mass,_=model.residual(star,delta,(control.ZERO,control.Z0),(control.ZERO,control.Z0),ld('1e-40'),model.weights(ld(1),None))
    return (np.column_stack([control.rows(spatial).reshape(n,3),-x[:,2]]).ravel(),
            np.column_stack([control.rows(mass).reshape(n,3),x[:,3]]).ravel())


def column(j):
    x=np.zeros(4*control.NC,dtype=ld)
    x[j]=control.EPS
    p,mp=evaluate_linear(x)
    q,mq=evaluate_linear(-x)
    return j,np.asarray((p-q)/(2*control.EPS),float),np.asarray((mp-mq)/(2*control.EPS),float)


def build_operator(cells,workers):
    control.STAR=model.initialize(None)
    control.ZERO=np.zeros_like(control.STAR.base)
    control.Z0=control.STAR.evaluate(control.ZERO)
    control.NC=cells
    star=control.STAR
    X=star.base[:cells+1,5:]
    star.control_basis=np.stack([np.vstack([X[:1],X[:-1]])[:cells]-X[:cells],X[1:]-X[:cells]],axis=1)
    with ProcessPoolExecutor(max_workers=workers,initializer=model.e.worker_init) as pool:
        view=SimpleNamespace(base=star.base[:cells+1],n=cells+1,pool=pool)
        data=model.composition_response(view,control.ZERO[:cells+1],dict(aux=control.Z0['aux'][:cells+1]))
    star.control_composition=data['coefficients'][:cells]
    star.composition_coordinates=np.zeros((cells,2),dtype=ld)
    star.__class__=FrozenCompatible
    control.Z0=star.evaluate(control.ZERO)
    probe=np.zeros((cells,4),dtype=ld)
    probe[:,3]=ld('.001')*np.sin(np.arange(cells)+ld('.3'))
    evaluate_linear(probe)
    image=np.einsum('nk,nki->ni',star.composition_coordinates,star.control_basis)
    star.composition_coordinates=np.zeros((cells,2),dtype=ld)
    values=[]
    for sign in [1,-1]:
        delta=control.ZERO.copy()
        delta[:cells,2]=sign*probe[:,3]*control.SCALE[2]
        z=star.evaluate(delta)
        value,_=model.residual(star,delta,(delta,z),(delta,z),ld(1),model.weights(ld(1),None))
        values.append(value[:cells,5:]-star.base[:cells,5:]*value[:cells,0,None])
    native_image=-(values[0]-values[1])/2
    image_error=float(np.max(abs(image-native_image))/np.max(abs(native_image)))
    assert image_error < .005,image_error
    count=4*cells
    K,M=np.empty((count,count)),np.empty((count,count))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for j,k,m in pool.map(column,range(count),chunksize=1):
            K[:,j],M[:,j]=k,m
            if (j+1)%96==0: print('PREFLIGHT OPERATOR',cells,j+1,count,flush=True)
    direction=np.sin(np.arange(count)+ld('.3'))*ld('1e-4')
    p,mp=evaluate_linear(direction)
    q,mq=evaluate_linear(-direction)
    checks=[]
    for matrix,actual in [(K,(p-q)/2),(M,(mp-mq)/2)]:
        actual=actual.reshape(cells,4)
        predicted=(matrix@np.asarray(direction,float)).reshape(cells,4)
        checks.append((np.max(abs(predicted-actual),axis=0)/np.maximum(np.max(abs(actual),axis=0),ld('1e-30'))).astype(float).tolist())
    return M,K,checks,image_error


def integrate_linear(M,K,initial,plan,refinement):
    times=model.parent.wall.prior.time_nodes(plan,refinement)
    previous,older=initial.copy(),initial.copy()
    cache,history={},[]
    for step in range(1,len(times)):
        h=times[step]-times[step-1]
        c0,c1,c2=map(float,model.weights(h,None if step==1 else times[step-1]-times[step-2]))
        key=c0,float(h)
        if key not in cache:
            matrix=c0*M+float(h)*K
            scale=np.max(abs(matrix),axis=1)
            matrix=matrix/scale[:,None]
            cache[key]=lu_factor(matrix),scale,matrix
        factor,scale,matrix=cache[key]
        rhs=(-c1*(M@previous)-c2*(M@older))/scale
        current=lu_solve(factor,rhs)
        assert np.max(abs(matrix@current-rhs))/(1+np.max(abs(rhs)))<1e-10
        older,previous=previous,current
        if step%refinement==0:
            history.append(current.reshape(control.NC,4)[:16,:3]*np.asarray(control.SCALE,float))
    return np.asarray(history)


def pulse_selection(M,K,plan,gate):
    star,z=control.STAR,control.Z0
    count=control.NC
    A=solve(M,-K)
    values=eig(A,right=False)
    duration=plan['duration_seconds']
    assert max(values.real)*duration <= gate['native_core_growth_over_full_duration'],max(values.real)
    times=np.asarray(model.parent.wall.prior.time_nodes(plan,1)[1:],float)
    reference_matrix=expm(duration*A)
    propagators={}
    beta=(z['P']/z['rho']-z['aux'][:,9])/z['aux'][:,10]
    width=star.r[24]
    radial=np.asarray(star.r[:count]/width,float)
    profile=np.exp(-radial*radial)
    references,all_histories,records={}, {}, []
    refinements=[1,2,4,8,16,32,64,128]
    for pulse in [0,1]:
        perturbation=np.zeros((count,3))
        if pulse==0:
            perturbation[:,0]=1e-6*profile
            perturbation[:,1]=np.asarray(beta[:count],float)*perturbation[:,0]
        else:
            perturbation[:,2]=1e-8*radial*profile
        initial=np.column_stack([perturbation/np.asarray(control.SCALE,float),np.zeros(count)]).ravel()
        state=initial.copy()
        reference=[]
        previous_time=0.
        for t in times:
            h=float(t-previous_time)
            if h not in propagators:propagators[h]=expm(h*A)
            state=propagators[h]@state
            reference.append(state.reshape(count,4)[:16,:3]*np.asarray(control.SCALE,float))
            previous_time=t
        reference=np.asarray(reference)
        exact=(reference_matrix@initial).reshape(count,4)[:16,:3]*np.asarray(control.SCALE,float)
        normal=np.max(abs(reference),axis=(0,1))
        assert np.all(normal>0)
        cross=np.max(abs(reference[-1]-exact),axis=0)/normal
        assert max(cross)<gate['pulse_reference_crosscheck_relative'],cross
        histories={}
        for refinement in refinements:
            histories[refinement]=integrate_linear(M,K,initial,plan,refinement)
        errors=np.array([np.max(abs(histories[r]-reference),axis=(0,1))/normal for r in refinements])
        pairs=np.array([np.max(abs(histories[b]-histories[a]),axis=(0,1)) for a,b in zip(refinements[:-1],refinements[1:])])
        orders=np.log2(pairs[:-1]/pairs[1:])
        records.append(dict(pulse=pulse,relative_time_errors=errors.tolist(),pair_orders=orders.tolist(),reference_crosscheck=cross.tolist()))
        references[f'pulse_{pulse}']=reference
        all_histories.update({f'pulse_{pulse}_r{r}':h for r,h in histories.items()})
    passed=[]
    for j,r in enumerate(refinements[:-2]):
        if all(min(record['pair_orders'][j])>=gate['pulse_pair_order'] and max(record['relative_time_errors'][j+2])<=gate['finest_pulse_relative_error'] for record in records):
            passed.append(r)
    assert passed,('No predeclared time triple passed',records)
    return dict(classification='Counterexample candidate',maximum_core_mode_real=float(max(values.real)),
        refinements=refinements,pulses=records,passing_coarse_multipliers=passed,selected_time_multiplier=passed[0]),dict(**references,**all_histories)


def check(workers=15):
    design=json.loads((model.BASE/'design.json').read_text())
    gate=design['gates']
    assert not (model.BASE/'preflight.json').exists()
    assert model.e.digest(model.SOURCE/'time-refinement.json')==design['original_failed_verdict_sha256']
    assert model.prior.symbolic()['passed']
    assert acoustic.symbolic()['dual_face_work']==0
    dissipation=dissipative_symbolic()
    wave=acoustic_wave()
    assert min(wave['orders'])>=gate['manufactured_TOV_minimum_order'],wave
    star=model.initialize(None)
    known=[]
    for cells in [64,128,256]:
        # Verify the actual candidate pressure routine, not only the old toy.
        test=np.sin(np.arange(star.n,dtype=ld))
        t=(star.rf[1:-1]-star.r[:-1])/np.diff(star.r)
        assert np.array_equal(model.pressure_faces(star,test)[1:-1],test[:-1]+(1-t)*(test[1:]-test[:-1]))
        (original,dual,corrected),error=acoustic_checks(star,cells)
        assert max(original.real)>10 and max(abs(dual.real))<gate['acoustic_maximum_real_per_second']
        assert max(corrected.real)<=gate['acoustic_maximum_real_per_second']
        assert error<gate['acoustic_adjoint_relative']
        known.append(dict(cells=cells,old_growth=float(max(original.real)),dual_growth=float(max(abs(dual.real))),new_growth=float(max(corrected.real)),adjoint_relative=error))
    tov=manufactured_tov()
    assert np.min(tov['orders'])>=gate['manufactured_TOV_minimum_order'],tov
    native=native_checks(star)
    print('PREFLIGHT NATIVE',json.dumps(native),flush=True)
    M,K,derivative,image_error=build_operator(128,workers)
    assert np.max(derivative)<gate['derivative_relative'],derivative
    template=json.loads((model.SOURCE/'plan.json').read_text())
    np.savez_compressed(model.BASE/'preflight-operator.npz',M=M,K=K)
    pulses,histories=pulse_selection(M,K,template,gate)
    np.savez_compressed(model.BASE/'preflight-arrays.npz',M=M,K=K,**histories)
    result=dict(classification='Counterexample candidate',passed=True,
        implementation_sha256=model.e.digest(Path(model.__file__)),preflight_source_sha256=model.e.digest(Path(__file__)),
        design_sha256=model.e.digest(model.BASE/'design.json'),symbolic=acoustic.symbolic(),negative_and_positive_controls=known,
        manufactured_tov=tov,native_checks=native,derivative_checks=derivative,pulse_time_selection=pulses,
        acoustic_dissipation=dissipation,continuum_acoustic_wave=wave,native_species_image_relative_error=image_error,
        selected_time_multiplier=pulses['selected_time_multiplier'],full_native_pilot_passed=False,
        full_GR_time_convergence_fixed=False,physical_EOS_certified=False)
    model.e.write(model.BASE/'preflight.json',result)
    files=[Path(__file__),Path(model.__file__),model.BASE/'design.json',model.BASE/'preflight.json',model.BASE/'preflight-arrays.npz',model.BASE/'preflight-operator.npz']
    model.e.write(model.BASE/'preflight-manifest.json',dict(sha256={p.relative_to(model.e.ROOT).as_posix():model.e.digest(p) for p in files}))
    print('GR REDESIGN PREFLIGHT',json.dumps(result),flush=True)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('command',choices=['register','check'])
    parser.add_argument('--workers',type=int,default=15)
    args=parser.parse_args()
    assert set(os.sched_getaffinity(0))<=set(range(16)) and 1<=args.workers<=16
    register() if args.command=='register' else check(args.workers)
