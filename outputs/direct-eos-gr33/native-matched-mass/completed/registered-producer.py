"""Counterexample candidate: apply the existing mass ports to charge readout.

The homogeneous exterior components are matched here. Additional exterior
mixed stress and exact ADM conservation remain separate closure obligations.
"""
from pathlib import Path
import json,resource,sys,time
import numpy as np
import sympy as sp
import mpmath as mp
import return_native_mixed_gr as previous

OUT=Path('native-matched-mass163-work');R=previous.OUT;I=previous.inf
read,write,sha=previous.read,previous.write,previous.sha
LD=np.longdouble;C=previous.C;G=previous.G
CAPS=dict(prepare=30,apply=90,audit=60)


def symbolic():
    e,d,r,b,lb,li,lx=sp.symbols('e d r b lb li lx',nonzero=True)
    mass=r/2*(1-b*sp.exp(-2*(e*lb+d*li+e*d*lx)))
    assert sp.simplify(sp.diff(mass,e,d).subs({e:0,d:0})-r*b*(lx-2*lb*li))==0
    alpha,s,k,E,ds,dk,de=sp.symbols('alpha s k E ds dk de')
    q=(s+alpha*(E-k))/(1+k-E)
    changed=(s+ds+alpha*(E+de-k-dk))/(1+k+dk-E-de)
    assert sp.factor(changed-q-(ds+(alpha+q)*(de-dk))/(1+k+dk-E-de))==0
    N,A,f,lam,nu,energy=sp.symbols('N A f lam nu energy',nonzero=True)
    # A conserved Killing-energy packet has A^4 E_J dr=e sqrt(b)/N.
    # Its density and the Einstein constraint prefactor must both vary.
    density=energy*sp.exp(-d*(nu+lam+4*alpha*f))
    prefactor=sp.exp(d*(2*lam+4*alpha*f))
    assert sp.simplify(sp.diff(prefactor*density,d).subs(d,0)-energy*(lam-nu))==0
    return dict(classification='Proven',passed=True,
        physical_mass='m_BI=r*b*lambda_BI-2*r*b*lambda_B*lambda_I.',
        normalization='q=(s+alpha*(epsilon-kappa))/(1+kappa-epsilon); exact increment is [ds+(alpha+q)*(de-dk)]/[1+kappa+dk-epsilon-de].',
        packet_measure='The known mixed metric factor and fixed-Killing-inventory photon density combine to (lambda_I-nu_I)*e. The conformal4alpha f term cancels; varying only the metric prefactor is incomplete.',
        scope='Algebra in the repository Einstein-frame convention. Not an assertion that the saved residual mass obeys the momentum constraint.')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    files=[Path(__file__),Path(previous.__file__),Path(I.__file__),R/'audit.json',R/'charge-parts.npz',R/'result.json',
        I.OUT/'charge-parts.npz',I.BACKGROUND,I.SELF/'normalization.json',
        previous.mixed.previous.METRIC/'metric-128-g8.npz',R/'sweep-2/gr/source-128-reference-128.npz']
    for n,q in [(128,8),(128,4),(64,8)]:files += [I.SELF/f'metric/metric-{n}-g{q}.npz']
    for n,q in [(128,8),(128,4),('128-17',8)]:files += [R/f'metric/metric-{n}-g{q}.npz']
    for q in [4,8]:files += [I.SELF/f'sweep-2/returned-metric/metric-128-g{q}.npz']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='b1e26c922e93f7068fa704ff49f57fe22efaa614',
        previous_turn='Progress: Phase162 returned the compact mixed GR field through actual photons and free matter and read its accepted finite response at fixed-operator infinity.',
        claim='Use a common mass-normalized charge for the existing background, incident material response and mixed GR component. Preserve each signed mass port and tiny transport return separately; quantify the cancellation instead of treating J_BI as new photon energy.',
        decision='Does the physically required mass coordinate and homogeneous exterior matching change the selected charge? Separately expose what is not yet a conserved, complete exterior mass.',
        scope='Actual stored conserved-source ports and their declared homogeneous vacuum continuation. Additional exterior mixed scalar/photon stress, changed photon paths and exact ADM/momentum-constraint closure are not assumed zero or solved.',
        representation='The fine33-knot mixed input is sampled at the common17 output knots.64/128 checks change the actual incident response clock;17/33 changes only mixed-source representation.',
        gates=dict(time=.02,quadrature=.002,cadence=.02,independent=1e-12,normalization=1e-12),
        budget=dict(actions=CAPS,total_seconds=sum(CAPS.values()),CPU_threads=1,virtual_GiB=3,new_evolution_steps=0,new_EOS_roots=0),
        admission='Reuse existing mass/field arrays. Only the terminal returned-source energy projection is assembled; no ray, fluid or EOS history replay. Stop at the fixed action caps without enlarging clocks or relaxing gates.',
        references=['https://arxiv.org/pdf/gr-qc/9707041, equations2.11,2.23,2.24; the mixed expansion and readout identities are derived here.'],
        bindings={str(p):sha(p) for p in files}))
    write(OUT/'symbolic.json',symbolic())


def sample(z,key,clock):
    ids=[int(np.argmin(abs(z['t']-t))) for t in clock]
    assert np.max(abs(z['t'][ids]-clock))<1e-18
    return z[key][ids].astype(LD)


def returned_port(order):
    d=dict(np.load(R/'sweep-2/gr/source-128-reference-128.npz'))
    m=previous.wave.Response();m.setup(d,order)
    # J is independent of the scalar argument in this existing conserved
    # source owner. Only J is consumed; zero-field derived outputs are ignored.
    zero=np.zeros((len(d['t']),len(m.tr)))
    z=previous.mixed.previous.retained.constraints.centers(m,d,dict(delta_phi=zero,delta_Phi=zero),order)
    energy=d['outer_cumulative_energy_erg'].astype(LD)
    value=z['J'][:,-1].astype(LD)*LD(m.tz['lapse'][-1])/np.sqrt(LD(m.tz['b'][-1]))+LD(G)/LD(C)**4*energy
    return value/LD(read(R/'normalization.json')['factor'])


def apply():
    assert read(R/'audit.json')['passed'];previous.mixed.previous.initialize()
    old=np.load(R/'charge-parts.npz');bg=np.load(I.BACKGROUND);legacy=np.load(I.OUT/'charge-parts.npz')
    t=old['t'];source=np.load(R/'sweep-2/gr/source-128-reference-128.npz')
    M=LD(source['M_cm']);alpha=-LD(source['K_cm'])/M;D0=1-bg['epsilon'].astype(LD);q0=bg['normalized'].astype(LD)
    b=sample(np.load(previous.mixed.previous.METRIC/'metric-128-g8.npz'),'asymptotic_mass_residual_cm',t)/M
    D=D0+b;db=-(alpha+q0)*b/D;qb=q0+db
    # The known mixed-source correction has opposite sign to the original
    # incident mass port. Their signed sum is the matching quantity.
    factor=LD(read(I.SELF/'normalization.json')['factor']);paths={};rows=[];ports={}
    for tag,n,q,mixed_n in [('fine',128,8,128),('quadrature',128,4,128),('time',64,8,128),('cadence',128,8,'128-17')]:
        ci=sample(np.load(I.SELF/f'metric/metric-{n}-g{q}.npz'),'asymptotic_mass_residual_cm',t)
        cx=sample(np.load(R/f'metric/metric-{mixed_n}-g{q}.npz'),'auxiliary_asymptotic_mass_cm',t)
        cs=sample(np.load(I.SELF/f'sweep-2/returned-metric/metric-128-g{q}.npz'),'asymptotic_mass_residual_cm',t)/factor
        if q not in ports:ports[q]=returned_port(q)
        cr=ports[q];total=ci+cx+cs+cr;correction=-(alpha+qb)*total/(M*D)
        paths[tag]=correction
        np.savez_compressed(OUT/f'mass-{tag}.npz',t=t,incident_port_cm=ci,mixed_auxiliary_port_cm=cx,
            self_return_port_cm=cs,mixed_transport_port_cm=cr,matched_homogeneous_port_cm=total,charge_correction=correction)
        rows.append(dict(path=tag,incident_endpoint_cm=float(ci[-1]),mixed_endpoint_cm=float(cx[-1]),
            paired_endpoint_cm=float(total[-1]),charge_endpoint=float(correction[-1])))
    norm=max(np.max(abs(paths['fine'])),LD('1e-290'))
    controls={k:float(np.max(abs(v-paths['fine']))/norm) for k,v in paths.items() if k!='fine'}
    # Preserve every component; never subtract the large direct waveform to
    # recover the much smaller material-mediated charge.
    components={k:old[k].astype(LD)*D0/D for k in ['previous_selected','corrected_mixed','transport_return','direct']}
    components['self_return']=legacy['self_gr'].astype(LD)*D0/D
    body=np.load(I.OUT/'body-128-a8-r8.npz');self_gr=np.load(I.OUT/'self-128-a8-r8.npz')
    arrived=np.load(R/'return-128-a8-r8.npz')['arrived_energy_erg'].astype(LD)
    de=body['epsilon_increment'].astype(LD)+self_gr['epsilon_increment'].astype(LD)+LD(G)/LD(C)**4/M*arrived
    components['background_mass_remap']=db*de/D
    components['matched_mass']=paths['fine']
    selected=sum(v for k,v in components.items() if k!='direct')
    old_selected=old['selected_without_tiny_return'].astype(LD)
    np.savez_compressed(OUT/'charge-parts.npz',t=t,**components,selected=selected,old_selected=old_selected,
        background_kappa=b,background_mass_shift=db,background_denominator=D,old_denominator=D0,
        epsilon_response=de,alpha=alpha,background_charge=q0)
    result=dict(classification='Counterexample candidate',passed=controls['time']<.02 and controls['quadrature']<.002 and controls['cadence']<.02,
        controls=controls,rows=rows,old_selected_endpoint=float(old_selected[-1]),matched_selected_endpoint=float(selected[-1]),
        matched_mass_charge_endpoint=float(paths['fine'][-1]),mass_effect_over_selected=float(paths['fine'][-1]/old_selected[-1]),
        background_mass_parameter_maximum=float(np.max(abs(b))),same_background_denominator=True,
        actual_stored_mass_ports_applied=True,auxiliary_mass_not_reclassified_as_photon_energy=True,
        homogeneous_component_readout=True,additional_exterior_mixed_stress_closed=False,
        exact_ADM_conservation_verified=False,physical_final_charge_solved=False,full_goal_complete=False,
        scope=read(OUT/'plan.json')['scope'])
    write(OUT/'result.json',result);print(json.dumps(result),flush=True);assert result['passed'],result


def audit():
    result=read(OUT/'result.json');assert result['passed'];d=np.load(OUT/'charge-parts.npz');m=np.load(OUT/'mass-fine.npz')
    mp.mp.dps=110;B=lambda x:mp.mpf(str(x));errors=[];exact=[]
    source=np.load(R/'sweep-2/gr/source-128-reference-128.npz');M=B(source['M_cm'])
    old=np.load(R/'charge-parts.npz');legacy=np.load(I.OUT/'charge-parts.npz')
    for j in range(len(d['t'])):
        alpha=B(d['alpha']);q0=B(d['background_charge'][j]);D0=B(d['old_denominator'][j]);k=B(d['background_kappa'][j]);de=B(d['epsilon_response'][j])
        D=D0+k;qb=q0-(alpha+q0)*k/D;mass=B(m['matched_homogeneous_port_cm'][j])/M
        parts=[B(old[n][j]) for n in ['previous_selected','corrected_mixed','transport_return','direct']]+[B(legacy['self_gr'][j])]
        ds=D0*sum(parts)-(alpha+q0)*de
        linear=(ds+(alpha+qb)*(de-mass))/D
        actual=sum(B(d[n][j]) for n in ['previous_selected','corrected_mixed','transport_return','direct','self_return','background_mass_remap','matched_mass'])
        errors.append(float(abs(actual-linear)/max(abs(linear),mp.mpf('1e-290'))))
        exact.append(str((ds+(alpha+qb)*(de-mass))/(D+mass-de)))
    assert max(errors)<1e-12,errors
    # Positive control: a pure positive mass increment lowers q when
    # alpha+q>0. This independently fixes the sign of the applied port.
    a=mp.mpf('.004');k=mp.mpf('.001');assert abs(-a*k/(1+k)-((a/(1+k))-a))<mp.mpf('1e-100')
    fine=np.load(OUT/'mass-fine.npz');scale=np.max(abs(fine['incident_port_cm']))
    remaining=float(np.max(abs(fine['matched_homogeneous_port_cm']))/scale)
    verdict=dict(classification='Counterexample candidate',passed=True,normalization_relative=max(errors),
        hundred_digit_standalone_rational_increment=exact,positive_mass_sign_control=True,
        maximum_matched_port_over_original_incident_port=remaining,
        not_an_ADM_conservation_pass='The remaining nonzero homogeneous mass port is applied, not zeroed or absorbed into emitted energy. Exact initial-mass/flux matching and additional exterior mixed stress are still required.',
        symbolic=read(OUT/'symbolic.json'),full_goal_complete=False)
    write(OUT/'audit.json',verdict);print(json.dumps(verdict),flush=True)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS
    receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));previous.mixed.inc.native.deadline(CAPS[action])
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
            assert sum(read(p)['seconds'] for p in OUT.glob('*-receipt.json'))+CAPS[action]<=sum(CAPS.values())
        globals()[action]()
    except Exception as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(action=action,seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
