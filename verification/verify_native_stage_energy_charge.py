"""Apply captured stage histories to the actual retarded charge and its bound."""
from pathlib import Path
import json,signal,time
import numpy as np
import mpmath as mp
import sympy as sp
import def_native_stage_energy_history as run
import verify_native_anisotropic_gr as independent
import def_native_global_scalar_closure as exterior

OUT=run.OUT/'gr';flow=run.flow;gr=run.gr;C=run.C;G=flow.G;LD=np.longdouble;write=run.write;sha=run.sha


def read(model,d,order):
    model.setup(d,order)
    # Only the requested charge is needed. Reuse characteristic integration,
    # with one outer observer; enclose the global potential analytically.
    model.tx=model.tx[-1:];model.distance=abs(model.tx[:,None]-model.x[None,:])/C;model.sign=np.sign(model.tx[:,None]-model.x[None,:])
    u=model.propagate(model.source)[0][:,0];direct=model.propagate(model.direct)[0][:,0]
    return dict(t=model.t,free_scalar=-u/float(d['M_cm']),direct_scalar=-direct/float(d['M_cm']))


def enclosure(model,d,history):
    previous=json.loads((exterior.OUT/'bound.json').read_text())
    # Reuse only coefficient geometry for the IDENTICAL current background;
    # source norms, emitted energy, port and floor terms are rebuilt here.
    p=json.loads((exterior.OUT/'bound-plan.json').read_text())
    bgfile=str(flow.INPUT/'balanced-20.npz');assert sha(bgfile)==p['bindings'][bgfile]
    model.setup(d,8);iv=mp.iv;iv.dps=40;I=iv.mpf
    def B(v):
        n,d=v.as_integer_ratio();return I(n)/I(d)
    def up(v):return float(np.nextafter(float(v.b),np.inf))
    def lo(v):return float(np.nextafter(float(v.a),-np.inf))
    T=B(float(d['t'][-1]));cc=B(float(C));gg=B(float(G));M=B(float(d['M_cm']));K=B(abs(float(d['K_cm'])))
    r0=B(model.bg.rend);rmin=I(previous['causal_radius_min_cm']);bounds=previous['frozen_polynomial_bounds']
    Nmin=I(bounds['lapse']['lower']);phi=I(bounds['phi']['absolute']);alpha=4*phi;Phi=I(bounds['Phi']['absolute'])
    amin=Nmin*iv.exp(-2*phi*phi);bmin=1-2*I(bounds['mass']['upper'])/rmin
    length=(r0-rmin)/(Nmin*iv.sqrt(bmin));ke=I(previous['mass_source_coefficient_absolute_per_cm2']);eta=I(previous['global_potential_contraction'])
    density=model.source.reshape(len(model.t),-1,8)/model.dx.reshape(-1,8)
    coefficients=np.einsum('cij,tcj->tci',model.inverse,density)
    maxima=np.max(abs(coefficients),axis=0);integral=I(0)
    for j,row in enumerate(maxima):
        if model.faces[j+1]<float(rmin.a):continue
        integral+=sum((B(v) for v in row),I(0))*2*B(model.half[j])
    source_norm=cc*T/2*integral
    assert np.all(np.diff(d['outer_cumulative_energy_erg'])>=0)
    emitted=B(float(d['outer_cumulative_energy_erg'][-1]));b0=1-2*M/r0
    N0=B(float(model.bg.metric(np.array([model.bg.rend/model.model.m.R]))[1][0]));c0=N0*iv.sqrt(b0)
    kappa=M/(r0*b0)+K*K/(2*c0*c0*r0*r0);assert up(kappa)<1
    e=gg/cc**4*emitted;ext_stress=e*K*iv.pi/(4*c0*c0*r0*iv.sqrt(1-kappa));ext_mass=cc*T/2*e*K/(c0*b0*r0*r0)
    disc=history['discard'];scale=LD(model.model.gas_scale);cx=LD(d['cx']);a0=LD(model.model.m.a0)
    floor_K=(disc[:,2].astype(LD)+a0*cx*disc[:,0])*scale
    assert floor_K.min()>=0 and np.all(np.diff(floor_K)>=0)
    floor_max=B(float(np.max(floor_K)));escape_max=B(float(np.max(abs(history['escape']))))
    q=flow.initial.Quadrature(d['edges'],8);_,_,a,BB,_=model.geo(q.r.ravel()-model.model.m.RJ)
    measure=q.r*q.r*BB.reshape(q.r.shape);mean_a=((measure*a.reshape(q.r.shape))@q.w)/(measure@q.w)
    energy=d['baryon_g'].astype(LD)*cx*LD(C)**2+d['gas_nonrest_energy_erg']+d['photon_energy_erg']
    mismatch=np.sum(energy*mean_a,axis=1,dtype=LD)+floor_K+history['escape']-d['inner_cumulative_energy_erg']+d['outer_cumulative_energy_erg']
    residual=B(float(np.max(abs(mismatch))));mass_uncertainty=floor_max+escape_max+residual
    mass_error=cc*T/2*length*ke*gg/cc**4/Nmin*mass_uncertainty
    # Conditional DIRECT source envelope for the discarded positive matter:
    # 0<=principal pressures<=E implies |trace|<=2E, |E-Pr|<=E.
    # Subsequent opacity/chemistry/transport feedback is not enclosed here.
    missing_proper=floor_max/amin
    floor_direct=gg/(2*cc**3)*T*alpha/rmin*2*missing_proper
    floor_stress=gg/(2*cc**3)*T*Phi*missing_proper
    potential=eta/(1-eta)*(source_norm+ext_stress+ext_mass+mass_error+floor_direct+floor_stress)
    unknown=(ext_stress+ext_mass+mass_error+floor_direct+floor_stress+potential)/M
    np.savez_compressed(OUT/'bound-inputs.npz',source_coefficient_maxima=maxima,source_optical_cell_widths=2*model.half,
        mean_lapse=mean_a,continuous_piecewise_linear_port_residual=mismatch,floor_Killing_energy=floor_K)
    return dict(classification='Counterexample candidate',global_potential_contraction=up(eta),absolute_source_norm_cm=up(source_norm),
        potential_normalized_bound=up(potential/M),exterior_stress_bound=up(ext_stress/M),exterior_mass_bound=up(ext_mass/M),
        floor_direct_trace_bound=up(floor_direct/M),floor_metric_stress_bound=up(floor_stress/M),mass_constraint_uncertainty_bound=up(mass_error/M),
        actual_emitted_energy_erg=float(emitted.a),floor_Killing_energy_erg=up(floor_max),remaining_source_port_error_erg=up(residual),
        remaining_source_port_relative=up(residual)/float(d['inner_cumulative_energy_erg'][-1]),total_conditional_normalized_bound=up(unknown),
        exact_stage_ports_used=True,piecewise_linear_cumulative_histories=True,
        scope='Current saved first-variation operator and reconstructed source. Same geometry envelope; new source and loss norms. Conditional on positive discarded matter with principal pressures between0 andE and no larger subsequent energy. This bounds its direct GR omission, not changes to opacity/chemistry/transport.',
        full_floor_feedback_enclosed=False,full_source_error_enclosed=False,nonlinear_GR=False)


def port_origin(d):
    old=np.load(gr.base.OUT/'source-128.npz');t=d['t'];h=t[1]-t[0]
    cumulative=d['inner_cumulative_energy_erg']-d['outer_cumulative_energy_erg']
    rates=np.r_[old['inner_luminosity'][0]-old['outer_luminosity'][0],np.diff(cumulative)/h]
    trap=flow.green.polynomial(t,rates).antiderivative()(t)
    method=cumulative-trap;exact=h/2*(rates-rates[0]);scale=float(d['inner_cumulative_energy_erg'][-1])
    identity=float(np.max(abs(method-exact))/scale);assert identity<1e-12
    old_net=old['inner_cumulative_energy_erg']-old['outer_cumulative_energy_erg']
    cadence=trap[::8]-old_net
    return dict(classification='Counterexample candidate',implicit_stage_vs_trapezoid_identity=identity,
        maximum_implicit_stage_rule_offset=float(np.max(abs(method))/scale),
        maximum_remaining_sparse_history_offset=float(np.max(abs(cadence))/scale),
        endpoint_implicit_stage_rule_offset=float(method[-1]/scale),endpoint_sparse_history_offset=float(cadence[-1]/scale),
        endpoint_total_recovered_port_difference=float((cumulative[-1]-old_net[-1])/scale),
        explanation='Accepted backward-Euler stage integral differs from trapezoids even at every original time step. The remaining difference is sparse history reconstruction, including any local-half-stage effect on sampled luminosity. Neither is a failure of conserved-state energy conversion.')


def main():
    assert not OUT.exists();OUT.mkdir();start=time.monotonic();assert json.loads((run.OUT/'result.json').read_text())['passed']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',claim='Apply recovered exact stage ports and original-step source moments to retarded scalar charge; rebuild all source-dependent GR bounds and explicitly carry discarded-matter trace/stress uncertainty.',
        decision='Does the remaining scalar charge survive removing the artificial sparse port mismatch and the newly bounded direct floor terms?',
        reuse='Same531 cells and original64/128 clocks. One outer observer instead of unnecessary full field output. Reuse only verified same-background coefficient envelopes, not old source norms or a lower bound.',
        gates=dict(time=.02,quadrature=.002,independent=1e-9,energy=1e-8),budget_seconds=90,CPU_threads=1,memory_GB=3,
        stop='No automatic resolution or horizon increase. Keep empirical time comparisons separate from the conditional uniform GR/source-loss bound.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(run.__file__),Path(gr.__file__),run.OUT/'result.json',run.OUT/'source-64.npz',run.OUT/'source-128.npz',exterior.OUT/'bound.json',flow.INPUT/'balanced-20.npz']}))
    signal.signal(signal.SIGALRM,flow.old.optical.timeout);signal.alarm(90);m=gr.Response();paths={}
    for steps,order in [(128,8),(64,8),(128,4)]:
        data=dict(np.load(run.OUT/f'source-{steps}.npz'));wave=read(m,data,order);paths[steps,order]=wave
        np.savez_compressed(OUT/f'wave-{steps}-g{order}.npz',**wave)
    d=dict(np.load(run.OUT/'source-128.npz'));fine=paths[128,8];coarse=paths[64,8];quad=paths[128,4]
    assert fine['free_scalar'][0]==0
    norm=max(np.max(abs(fine['free_scalar'])),1e-300)
    time_error=float(np.max(abs(fine['free_scalar'][::2]-coarse['free_scalar']))/norm);quadrature=float(np.max(abs(fine['free_scalar']-quad['free_scalar']))/norm)
    direct,error=independent.direct(m,d,8);agreement=abs(direct/fine['direct_scalar'][-1]-1)
    sparse={k:(v[::8] if isinstance(v,np.ndarray) and v.ndim and v.shape[0]==len(d['t']) else v) for k,v in d.items()}
    sparse_wave=read(m,sparse,8);np.savez_compressed(OUT/'wave-same-path-17-knots.npz',**sparse_wave)
    cadence=float(np.max(abs(sparse_wave['free_scalar']-fine['free_scalar'][::8]))/norm)
    bound=enclosure(m,d,np.load(run.OUT/'history-128.npz'));origin=port_origin(d);end=float(fine['free_scalar'][-1]);uncertainty=bound['total_conditional_normalized_bound']
    lower=float(np.nextafter(end-uncertainty,-np.inf));upper=float(np.nextafter(end+uncertainty,np.inf))
    old=np.load(gr.OUT/'fields-128-g8.npz');old_end=-float(old['direct_and_mass_stress_U'][-1,-1])/float(d['M_cm'])
    # A cumulative history, rather than endpoint-rate trapezoids, represents
    # the accepted implicit Euler transfer. Test the nonconstant-rate case.
    h,L0,L1=sp.symbols('h L0 L1',positive=True)
    defect=sp.expand(h*L1-h*(L0+L1)/2);assert sp.simplify(defect-h*(L1-L0)/2)==0
    passed=time_error<.02 and quadrature<.002 and agreement<1e-9 and bound['remaining_source_port_relative']<1e-8 and lower>0
    result=dict(classification='Counterexample candidate',passed=bool(passed),time_relative=time_error,quadrature_relative=quadrature,
        independent_direct_relative=agreement,retarded_coordinate_inverse_residual_cm=error,captured17_vs_all_knots_relative=cadence,
        endpoint_free_scalar=end,previous_sparse_endpoint=old_end,endpoint_relative_change=end/old_end-1,
        endpoint_direct=float(fine['direct_scalar'][-1]),conditional_scalar_interval=[lower,upper],bound=bound,port_origin=origin,
        nonconstant_stage_flux_symbolic_check=True,actual_stage_histories_applied_to_GR=True,
        interpretation='Artificial sparse-port mismatch resolved at its owner. Conditional scalar positivity includes exterior/potential and a direct discarded-matter stress envelope, not a full EOS, time-continuum, floor feedback, nonlinear or observational certificate.',
        coupled_fixed_point_verified=False,full_source_error_enclosed=False,final_charge_solved=False,full_goal_complete=False,seconds=time.monotonic()-start)
    write(OUT/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
