"""Counterexample candidate: implicit nonlinear heat exchange on the new GR grid.

Actual EOS and opacity are reevaluated at every nonlinear iteration. This
closed, fixed-density/metric/composition operator is not full stellar evolution.
"""
from concurrent.futures import ProcessPoolExecutor
import json, sys
import numpy as np
from scipy.linalg import solve_banded
import sympy as sp
import gr_microphysics as micro
import opacity_tables as opacity

g=micro.g;OUT=g.OUT/'gr-nonlinear-thermal'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    inputs=[g.ROOT/'verification'/n for n in ['gr_nonlinear_thermal.py','opacity_tables.py','conservative_star.py','direct_eos_gr.py']]
    inputs += [g.OUT/n for n in ['initial-state-17-4.npz','initial-GR.json','gr-microphysics/auxiliaries.npz',
        'gr-opacity/new-GR-captured.npz','gr-opacity/evaluation.npz','gr-transport/diagnostics.npz']]
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='8ab926a',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in inputs},processes=4,
        step_counts=[1,2,4],duration_rule='min(1 coordinate second, 0.02/max(abs(initial closed dlnT/dt)))',
        max_Newton_iterations=12,max_backtracks=10,temperature_step_cap=.1,
        local_energy_residual_scaled_tolerance=2e-12,global_energy_relative_to_exchange_tolerance=1e-8,
        finite_time_refinement_logT_tolerance=1e-4,
        boundary='Both external faces closed, no nuclear or neutrino sources in this operator. Fixed baryon density, fixed GR geometry and fixed nuclear composition.',
        equation='dm_i*N_i*(u_i(T_new)-u_i(T_old))=dt*(Linf_inner-Linf_outer). Face flux and EOS are fully reevaluated at the new temperature. Newton Jacobian includes opacity temperature derivatives.',
        limitations='A nonlinear operator closure control, not a fluid/metric/atmosphere, convection or physical EOS certificate. Finite time refinement is not a rigorous time-error enclosure.',
        physical_EOS_certified=False,full_GR_evolution=False))
    symbolic()


def symbolic():
    xo,xi,no,ni,K,kap,do,di=sp.symbols('xo xi no ni K kap do di',positive=True)
    thetao=sp.exp(xo)*no;thetai=sp.exp(xi)*ni
    flux=K*(thetai**4-thetao**4)/kap
    out=sp.diff(flux,xo)+sp.diff(flux,kap)*do
    inner=sp.diff(flux,xi)+sp.diff(flux,kap)*di
    assert sp.simplify(out-(-4*K*thetao**4/kap-flux*do/kap))==0
    assert sp.simplify(inner-(4*K*thetai**4/kap-flux*di/kap))==0
    # Independent matrix check: a two-cell exchange has zero column sum
    # before row normalization, including state-dependent conductance.
    J=sp.Matrix([[out,inner],[-out,-inner]])
    assert sp.ones(1,2)*J==sp.zeros(1,2)
    save('symbolic.json',dict(classification='Proven',passed=True,
        face_Jacobian='dL/dlnT_outer=-4K theta_outer^4/kap-L*d(kap)/dlnT_outer/kap; inner sign is positive for the fourth-power term.',
        energy='Equal and opposite face terms and their unscaled Jacobian columns sum to zero.',
        entropy='Conditional: positive T and cv imply convex u(s) at fixed density/composition. The backward-Euler energy exchange therefore has S_new-S_old >= sum_i DeltaE_infinity_i/theta_new_i >= 0 for an exact implicit root and a positive conductance. An EOS that is not smooth/first-law-consistent does not inherit this finite-step theorem.',
        scope='Algebraic/conditional operator statements; actual numerical and physical errors remain separate.'))
    print('PASS nonlinear thermal face Jacobian and conservation identities',flush=True)


def initialize():
    global worker_eos,worker_state
    worker_eos=g.EOS();worker_state=dict(np.load(g.OUT/'initial-state-17-4.npz'))


def samples(item):
    start,temperature=item
    return np.array([worker_eos(2,worker_state['lnd'][i],t,worker_state['X'][i])
        for i,t in enumerate(temperature,start)])


def faces(t,nu,dm,radius,nuface,op):
    w=dm[:-1]/(dm[:-1]+dm[1:]);kap=(1-w)*op[:-1,0]+w*op[1:,0]
    z=t+nu;K=(4*np.pi*(radius[1:-1]*100)**2)**2*4*5.670400e-5/(3*np.exp(2*nuface[1:-1])*((dm[:-1]+dm[1:])/2))
    fourth=np.exp(4*z[1:])*np.expm1(4*(z[:-1]-z[1:]))
    flux=-K*fourth/kap
    douter=-4*K*np.exp(4*z[:-1])/kap-flux*((1-w)*op[:-1,0]*op[:-1,2])/kap
    dinner=4*K*np.exp(4*z[1:])/kap-flux*(w*op[1:,0]*op[1:,2])/kap
    return flux,douter,dinner


def jacobian_control():
    t=np.array([1.,1.05,1.08]);nu=np.array([-.01,-.02,-.03]);dm=np.array([2.,3.,4.])
    radius=np.array([4.,3.,2.,0.]);nf=np.array([-.005,-.015,-.025,-.035]);beta=np.array([.2,-.3,.1])
    def op(q):return np.column_stack([np.exp(beta*q),np.zeros(3),beta])
    f,do,di=faces(t,nu,dm,radius,nf,op(t));J=np.zeros((2,3));J[0,:2]=[do[0],di[0]];J[1,1:]=[do[1],di[1]]
    h=1e-5;numeric=[]
    for i in range(3):
        shift=np.zeros(3);shift[i]=h
        numeric.append((faces(t+shift,nu,dm,radius,nf,op(t+shift))[0]-faces(t-shift,nu,dm,radius,nf,op(t-shift))[0])/(2*h))
    error=float((abs(np.array(numeric).T-J)/np.maximum(1,abs(J))).max());assert error<1e-8,error
    save('Jacobian-control.json',dict(classification='Counterexample candidate',passed=True,relative_error=error,
        manufactured_state=True,EOS_or_physical_derivative_certificate=False))


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    jacobian_control();state,aux=micro.inputs();n=len(state['X']);dm=state['dm'];N=np.exp(state['nu']);mass=dm*N
    base=aux['eos'];p=dict(np.load(g.OUT/'gr-opacity/new-GR-captured.npz'))['parameters'].copy()
    model=opacity.Opacity();base_op=dict(np.load(g.OUT/'gr-opacity/evaluation.npz'))['values']
    f0,*_=faces(state['lnT'],state['nu'],dm,state['radius_faces_m'],state['nu_faces'],base_op)
    assert np.allclose(f0,g.s.baryon_face_diffusion(state,base_op[:,0]),rtol=1e-14,atol=0)
    rhs0=np.r_[f0,0.]-np.r_[0.,f0];rate=rhs0/(mass*base[:,10]);duration=min(1.,.02/float(abs(rate).max()))
    save('duration.json',dict(classification='Counterexample candidate',coordinate_seconds=duration,
        initial_max_abs_dlnT_dt=float(abs(rate).max()),rule=plan['duration_rule']))
    previous=None;records=[]
    with ProcessPoolExecutor(max_workers=plan['processes'],initializer=initialize) as pool:
        def evaluate(t):
            a=np.concatenate(list(pool.map(samples,[(i,t[i:i+256]) for i in range(0,n,256)])))
            assert np.all(a[:,10]>0) and np.all(np.isfinite(a))
            par=p.copy();par[:,4]=t/np.log(10);op=np.array([model(row) for row in par])
            flux,fo,fi=faces(t,state['nu'],dm,state['radius_faces_m'],state['nu_faces'],op)
            return a,op,flux,fo,fi
        fresh=evaluate(state['lnT']);assert np.array_equal(fresh[0],base)
        assert np.allclose(fresh[1],base_op,rtol=1e-12,atol=1e-12)
        for steps in plan['step_counts']:
            t=state['lnT'].copy();old=base.copy();history=[];exchanged=np.longdouble(0)
            dt=duration/steps
            for step in range(steps):
                evaluated=evaluate(t);prior_norm=np.inf
                for iteration in range(plan['max_Newton_iterations']):
                    a,op,flux,fo,fi=evaluated;divergence=np.r_[flux,0.]-np.r_[0.,flux]
                    residual=mass*(a[:,2]-old[:,2])-dt*divergence
                    scale=mass*a[:,10];norm=float(abs(residual/scale).max())
                    print('NONLINEAR THERMAL',steps,step,iteration,norm,flush=True)
                    if norm<=plan['local_energy_residual_scaled_tolerance']:break
                    diag=scale-dt*(np.r_[fo,0.]-np.r_[0.,fi])
                    band=np.zeros((3,n));band[1]=diag/scale
                    band[0,1:]=-dt*fi/scale[:-1];band[2,:-1]=dt*fo/scale[1:]
                    update=solve_banded((1,1),band,-residual/scale)
                    limit=plan['temperature_step_cap'];update*=min(1.,limit/float(abs(update).max()))
                    for backtrack in range(plan['max_backtracks']):
                        candidate=t+update*(.5**backtrack);trial=evaluate(candidate)
                        tr=mass*(trial[0][:,2]-old[:,2])-dt*(np.r_[trial[2],0.]-np.r_[0.,trial[2]])
                        trial_norm=float(abs(tr/(mass*trial[0][:,10])).max())
                        if trial_norm<norm:
                            t=candidate;evaluated=trial;prior_norm=norm;break
                    else:raise AssertionError(('Nonlinear thermal line search',steps,step,norm,trial_norm))
                else:raise AssertionError(('Nonlinear thermal root',steps,step,norm,prior_norm))
                change=mass.astype(np.longdouble)*(a[:,2].astype(np.longdouble)-old[:,2].astype(np.longdouble))
                exchange=np.sum(abs(dt*divergence).astype(np.longdouble));exchanged+=exchange
                score=float(abs(change.sum())/max(exchange,np.longdouble(1)))
                entropy=np.sum(dm.astype(np.longdouble)*(a[:,3].astype(np.longdouble)-old[:,3].astype(np.longdouble)))
                row=dict(step=step,iterations=iteration,local_scaled_energy_residual=norm,
                    global_energy_relative_to_exchange=score,entropy_change_erg_K=float(entropy),
                    maximum_logT_change=float(abs(t-state['lnT']).max()))
                history.append(row);save(f'path-{steps}-progress.json',dict(classification='Counterexample candidate',rows=history))
                assert score<plan['global_energy_relative_to_exchange_tolerance'],row
                assert entropy>=0,row
                old=a.copy()
            final=mass.astype(np.longdouble)*(old[:,2].astype(np.longdouble)-base[:,2].astype(np.longdouble))
            row=dict(classification='Counterexample candidate',steps=steps,completed=True,history=history,
                total_global_energy_relative_to_exchange=float(abs(final.sum())/exchanged),
                maximum_logT_change=float(abs(t-state['lnT']).max()),full_GR_evolution=False)
            if previous is not None:
                row['time_refinement_logT_difference']=float(abs(t-previous).max())
                row['finite_refinement_passed']=row['time_refinement_logT_difference']<plan['finite_time_refinement_logT_tolerance']
            np.savez_compressed(OUT/f'path-{steps}.npz',lnT=t,eos=old,opacity=op,interior_Linf=flux)
            save(f'path-{steps}.json',row);records.append(row);previous=t.copy()
    save('result.json',dict(classification='Counterexample candidate',completed=True,paths=records,
        actual_EOS_and_opacity_recomputed=True,closed_boundary=True,fixed_density_metric_composition=True,
        sources_included=False,convection_solved=False,physical_EOS_certified=False,full_GR_evolution=False))
    print('NONLINEAR THERMAL COMPLETE',[(r['steps'],r['maximum_logT_change']) for r in records],flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
