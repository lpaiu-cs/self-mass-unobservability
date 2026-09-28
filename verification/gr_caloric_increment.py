"""Counterexample candidate: resolve heat increments before adding large offsets.

The fixed-density caloric chart integrates the actual EOS cv*T along log T.
Its quadrature/native link is tested separately from discrete conservation.
"""
from concurrent.futures import ProcessPoolExecutor
import json,sys
import numpy as np
from scipy.linalg import solve_banded
import sympy as sp
import gr_nonlinear_thermal as old

g=old.g;OUT=g.OUT/'gr-caloric-increment'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    diagnosis=json.loads((g.OUT/'gr-thermal-endpoint-replay/diagnosis.json').read_text())
    assert diagnosis['failure_reproduced'] and not diagnosis['original_global_energy_gate_passed']
    plan=json.loads((old.OUT/'plan.json').read_text())
    plan.update(checkpoint='9cf504d',processes=8,
        caloric_coordinate='U(delta)-U(0)=integral_0^delta [cv*T](lnT_old+z) dz. Store delta separately from the absolute binary64 log temperature. Four-point Gauss integration multiplies delta after summing capacities. No large absolute energy subtraction enters the residual.',
        Newton_method='Use endpoint cv*T as the caloric quasi-Newton diagonal and the original analytic face/opacity derivatives. Require both original local residual and original global energy criteria before accepting a step.',
        quadrature_nodes=4,independent_quadrature_nodes=2,finite_quadrature_difference_tolerance=1e-8,
        native_link='Report the original absolute-energy endpoint difference against the caloric increment with an explicit endpoint ulp budget. It is not used to silently certify the new coordinate or its continuous error.',
        limits='Fixed density/metric/composition closed operator. The integral chart equals the native energy only under first-law differentiability along the path. Native branch jumps, quadrature remainder, EOS primitive errors and physical calibration remain separate.',
        diagnostic_binding={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in [
            g.ROOT/'verification/gr_caloric_increment.py',g.OUT/'gr-thermal-endpoint-replay/manifest.json',
            g.OUT/'gr-nonlinear-thermal/failure-manifest.json']})
    save('plan.json',plan);symbolic()


def symbolic():
    h=sp.symbols('h',real=True);z=sp.symbols('z',real=True);c=sp.symbols('c0:4',real=True)
    P=sum(c[i]*z**i for i in range(4));nodes=[(1-1/sp.sqrt(3))/2,(1+1/sp.sqrt(3))/2]
    gauss=h*sum(P.subs(z,h*q) for q in nodes)/2
    assert sp.simplify(gauss-sp.integrate(P,(z,0,h)))==0
    t=np.array([1.,1.05,1.08]);nu=np.array([-.01,-.02,-.03]);dm=np.array([2.,3.,4.])
    radius=np.array([4.,3.,2.,0.]);nf=np.array([-.005,-.015,-.025,-.035])
    op=np.column_stack([np.ones(3),np.zeros(3),np.zeros(3)])
    ordinary=old.faces(t,nu,dm,radius,nf,op);shifted=fluxes(t,np.zeros(3),nu,dm,radius,nf,op)
    assert all(np.allclose(a,b,rtol=1e-14,atol=0) for a,b in zip(ordinary,shifted))
    delta=np.array([1e-18,0.,0.]);constant=np.full(3,16.)
    assert np.array_equal(constant+delta,constant)
    assert fluxes(constant,delta,np.zeros(3),dm,radius,nf,op)[0][0]!=0
    assert np.all(old.faces(constant+delta,np.zeros(3),dm,radius,nf,op)[0]==0)
    n=4;constant=sp.factorial(n)**4/((2*n+1)*sp.factorial(2*n)**3)
    save('symbolic.json',dict(classification='Proven',passed=True,
        caloric_identity='At fixed density and nuclear composition, if du/dlnT=cv*T along the path, Delta u=delta*integral_0^1 [cv*T](lnT_old+delta*x) dx. An arbitrary constant energy offset cancels exactly from this representation.',
        Gaussian_polynomial_control='Two-point Gauss integrates every cubic capacity exactly, including negative delta.',
        compensated_face_control='Zero shifts reproduce the original face value/Jacobian. A 1e-18 log-temperature difference at absolute logT=16 remains a nonzero flux in the separate-shift formula, while adding it to binary64 absolute temperatures erases it.',
        conditional_four_point_remainder=f'If the capacity has an eighth log-temperature derivative bounded by M8 on the segment, the exact four-point Gaussian remainder is at most {constant}*abs(delta)^9*M8. Node/evaluation rounding is additional.',
        physical_EOS_or_native_derivative_certificate=False))


def initialize():
    global eos,state,nodes4,weights4,nodes2,weights2
    eos=g.EOS();state=dict(np.load(g.OUT/'initial-state-17-4.npz'))
    nodes4,weights4=np.polynomial.legendre.leggauss(4);nodes4=(nodes4+1)/2;weights4=weights4/2
    nodes2,weights2=np.polynomial.legendre.leggauss(2);nodes2=(nodes2+1)/2;weights2=weights2/2


def samples(item):
    start,shift,delta,number=item;nodes,weights=(nodes4,weights4) if number==4 else (nodes2,weights2)
    rows=[];integrals=[];entropies=[]
    for j,(s,d) in enumerate(zip(shift,delta),start):
        def call(offset):return eos(2,state['lnd'][j],float(np.longdouble(state['lnT'][j])+offset),state['X'][j])
        if d==0:
            integrals.append(np.longdouble(0));entropies.append(np.longdouble(0))
            if number==4:rows.append(call(s))
            continue
        offsets=s+np.longdouble(d)*nodes
        capacities=np.array([call(offset)[10] for offset in offsets])
        assert np.all(capacities>0)
        integrals.append(np.longdouble(d)*np.sum(capacities.astype(np.longdouble)*weights.astype(np.longdouble)))
        temperatures=np.exp(np.longdouble(state['lnT'][j])+offsets)
        entropies.append(np.longdouble(d)*np.sum(capacities.astype(np.longdouble)*weights.astype(np.longdouble)/temperatures))
        if number==4:rows.append(call(s+np.longdouble(d)))
    return np.array(integrals),np.array(rows),np.array(entropies)


def fluxes(t0,shift,nu,dm,radius,nuface,op):
    # Preserve tiny shifts in the temperature difference even when adding
    # them to the large absolute log temperature would round to zero.
    w=dm[:-1]/(dm[:-1]+dm[1:]);kap=(1-w)*op[:-1,0]+w*op[1:,0]
    K=(4*np.pi*(radius[1:-1]*100)**2)**2*4*5.670400e-5/(3*np.exp(2*nuface[1:-1])*((dm[:-1]+dm[1:])/2))
    z=t0+nu;s=shift.astype(np.longdouble)
    fourth=np.exp(4*z[1:].astype(np.longdouble))*np.exp(4*s[1:])*np.expm1(4*((z[:-1]-z[1:]).astype(np.longdouble)+s[:-1]-s[1:]))
    L=-K.astype(np.longdouble)*fourth/kap
    do=-4*K*np.exp(4*z[:-1])*np.exp(4*s[:-1])/kap-L*((1-w)*op[:-1,0]*op[:-1,2])/kap
    di=4*K*np.exp(4*z[1:])*np.exp(4*s[1:])/kap-L*(w*op[1:,0]*op[1:,2])/kap
    return L,do,di


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in {**plan['bindings'],**plan['diagnostic_binding']}.items():assert g.c.sha(g.ROOT/rel)==digest,rel
    state,aux=old.micro.inputs();base=aux['eos'];dm=state['dm'];mass=dm*np.exp(state['nu']);n=len(dm)
    p=dict(np.load(g.OUT/'gr-opacity/new-GR-captured.npz'))['parameters'].copy();model=old.opacity.Opacity()
    duration=json.loads((g.OUT/'gr-nonlinear-thermal/duration.json').read_text())['coordinate_seconds']
    records=[];previous=None;ld=np.longdouble
    with ProcessPoolExecutor(max_workers=plan['processes'],initializer=initialize) as pool:
        def evaluate(shift,delta,number=4):
            pieces=list(pool.map(samples,[(i,shift[i:i+128],delta[i:i+128],number) for i in range(0,n,128)]))
            integral=np.concatenate([r[0] for r in pieces])
            if number==2:return integral
            a=np.concatenate([r[1] for r in pieces]);par=p.copy()
            par[:,4]=np.asarray(state['lnT'].astype(ld)+shift+delta,float)/np.log(10)
            op=np.array([model(row) for row in par]);L,fo,fi=fluxes(state['lnT'],shift+delta,state['nu'],dm,state['radius_faces_m'],state['nu_faces'],op)
            return integral,a,op,L,fo,fi,np.concatenate([r[2] for r in pieces])
        for steps in plan['step_counts']:
            shift=np.zeros(n,dtype=ld);old_a=base.copy();total_energy=np.zeros(n,dtype=ld);exchanged=ld(0);history=[];dt=duration/steps
            for step in range(steps):
                delta=np.zeros(n);evaluated=evaluate(shift,delta)
                for iteration in range(plan['max_Newton_iterations']):
                    U,a,op,L,fo,fi,S=evaluated;div=np.r_[L,ld(0)]-np.r_[ld(0),L]
                    energy=mass.astype(ld)*U;target=ld(dt)*div;residual=energy-target
                    scale=mass.astype(ld)*a[:,10];norm=float(abs(residual/scale).max())
                    exchange=abs(target).sum();global_error=float(abs(residual.sum())/max(ld(1),exchange))
                    print('CALORIC INCREMENT',steps,step,iteration,norm,global_error,flush=True)
                    if norm<=plan['local_energy_residual_scaled_tolerance'] and global_error<=plan['global_energy_relative_to_exchange_tolerance']:break
                    band=np.zeros((3,n));band[1]=np.asarray((scale-dt*(np.r_[fo,ld(0)]-np.r_[ld(0),fi]))/scale,float)
                    band[0,1:]=np.asarray(-dt*fi/scale[:-1],float);band[2,:-1]=np.asarray(dt*fo/scale[1:],float)
                    update=solve_banded((1,1),band,np.asarray(-residual/scale,float))
                    update*=min(1.,plan['temperature_step_cap']/max(float(abs(update).max()),1e-300))
                    # The merit includes the global gate, so convergence of
                    # large local energy scales cannot mask a lost increment.
                    merit=max(norm/plan['local_energy_residual_scaled_tolerance'],global_error/plan['global_energy_relative_to_exchange_tolerance'])
                    for backtrack in range(plan['max_backtracks']):
                        proposed=delta+update*(.5**backtrack);trial=evaluate(shift,proposed)
                        tr=mass.astype(ld)*trial[0]-ld(dt)*(np.r_[trial[3],ld(0)]-np.r_[ld(0),trial[3]])
                        tn=float(abs(tr/(mass.astype(ld)*trial[1][:,10])).max())
                        tg=float(abs(tr.sum())/max(ld(1),abs(ld(dt)*(np.r_[trial[3],ld(0)]-np.r_[ld(0),trial[3]])).sum()))
                        if max(tn/plan['local_energy_residual_scaled_tolerance'],tg/plan['global_energy_relative_to_exchange_tolerance'])<merit:
                            delta=proposed;evaluated=trial;break
                    else:raise AssertionError(('Caloric line search',steps,step,norm,global_error,tn,tg))
                else:raise AssertionError(('Caloric root',steps,step,norm,global_error))
                lower=evaluate(shift,delta,2);qscore=float(np.sum(mass.astype(ld)*abs(U-lower))/max(exchange,ld(1)))
                native=a[:,2].astype(ld)-old_a[:,2].astype(ld)
                budget=32*(np.spacing(abs(a[:,2]))+np.spacing(abs(old_a[:,2]))).astype(ld)
                unresolved=abs(native-U)<=budget
                row=dict(step=step,iterations=iteration,local_scaled_energy_residual=norm,
                    global_energy_relative_to_exchange=global_error,finite_two_four_point_difference_relative_to_exchange=qscore,
                    finite_quadrature_passed=qscore<plan['finite_quadrature_difference_tolerance'],
                    caloric_chart_entropy_change_erg_K=float(dm.astype(ld)@S),
                    native_endpoint_difference_within_32ulp_cells=int(unresolved.sum()),
                    native_endpoint_maximum_difference_outside_32ulp=float(np.maximum(abs(native-U)-budget,0).max()))
                np.savez_compressed(OUT/f'endpoint-{steps}-{step}.npz',old_shift=shift,step_shift=delta,caloric_increment=U,
                    two_point_increment=lower,entropy_increment=S,eos=a,opacity=op,interior_Linf=L,energy_residual=residual)
                history.append(row);save(f'path-{steps}-progress.json',dict(classification='Counterexample candidate',rows=history))
                assert row['finite_quadrature_passed'],row
                assert row['caloric_chart_entropy_change_erg_K']>=0,row
                shift+=delta.astype(ld);total_energy+=energy;exchanged+=exchange;old_a=a.copy()
            record=dict(classification='Counterexample candidate',steps=steps,completed=True,history=history,
                global_increment_energy_relative_to_exchange=float(abs(total_energy.sum())/exchanged),maximum_logT_change=float(abs(shift).max()),
                physical_EOS_certified=False,continuous_caloric_error_certified=False,full_GR_evolution=False)
            if previous is not None:
                error=float(abs(shift-previous).max());record.update(time_refinement_logT_difference=error,
                    finite_refinement_passed=error<plan['finite_time_refinement_logT_tolerance'])
            np.savez_compressed(OUT/f'path-{steps}.npz',lnT_base=state['lnT'],lnT_shift=shift,eos=old_a,total_energy_increment=total_energy)
            save(f'path-{steps}.json',record);records.append(record);previous=shift.copy()
    save('result.json',dict(classification='Counterexample candidate',completed=True,paths=records,
        original_endpoint_failure_preserved=True,physical_EOS_certified=False,continuous_caloric_error_certified=False,
        fixed_density_metric_composition=True,full_GR_evolution=False))


if __name__=='__main__':globals()[sys.argv[1]]()
