"""Request32 independent checks; never change an active GR trajectory."""
import json, shutil, sys
import numpy as np
from scipy.linalg import eigvalsh_tridiagonal, solve_banded
import structured_enthalpy as e

ROOT_STATS=dict(calls=0,evaluations=0,maximum_score=0.,label='strict-entropy')


def strict_invert(eos,lp,entropy,eps,guess):
    """Same EOS, an energy-unit entropy residual; no analytic EOS guarantee."""
    t=float(guess);best=None;started=ROOT_STATS['evaluations']
    def sample(lt):
        a=eos(1,lp,lt,eps);temperature=np.exp(lt)
        H=a[2]+a[1]/a[0];C=a[10]-a[1]/a[0]*a[8]
        assert C>0 and np.isfinite(C)
        error=a[3]-entropy;budget=max(2.,32*np.spacing(abs(H)))
        ROOT_STATS['evaluations']+=1
        return (abs(temperature*error)/budget,a,lt,error,C/temperature,budget)
    for _ in range(15):
        trial=sample(t)
        if best is None or trial[0]<best[0]: best=trial
        if trial[0]<.25: break
        new=t-np.clip(trial[3]/trial[4],-.15,.15)
        if new==t: break
        t=new
    if best[0]>.25:
        low=high=best[2]
        for _ in range(3):
            low=np.nextafter(low,-np.inf);high=np.nextafter(high,np.inf)
            for value in [low,high]:
                trial=sample(value)
                if trial[0]<best[0]: best=trial
    ROOT_STATS['calls']+=1;ROOT_STATS['maximum_score']=max(ROOT_STATS['maximum_score'],float(best[0]))
    if best[0]>1:
        path=e.OUT/(ROOT_STATS['label']+'-failure.json')
        if not path.exists():
            e.save(path.name,dict(classification='Counterexample candidate',passed=False,
                logP=float(lp),entropy=float(entropy),X=eps.tolist(),guess=float(guess),
                best_lnT=float(best[2]),score=float(best[0]),budget_erg_g=float(best[5])))
        raise ValueError(('energy-unit entropy inverse',lp,best[0],best[5]))
    return best[1],best[2],ROOT_STATS['evaluations']-started


def strict_plan():
    assert not (e.OUT/'strict-entropy-plan.json').exists()
    e.save('strict-entropy-plan.json',dict(classification='Counterexample candidate',checkpoint='7c14403',
        reason='The direct one-step energy score passes at 9.95186e-7, but its 6.42e32 erg margin is smaller than the linearized entropy-inverse residual estimates. Keep that original numerical pass and the original running 2/4-step paths unchanged.',
        residual='At every isentropic EOS inversion require abs(T*(s-s_target)) <= max(2 erg/g,32 ulp(abs(H))). This changes numerical root termination, not the EOS, matter, boundary luminosity, or original energy gate.',
        solver='Fifteen positive-capacity Newton trials and up to three neighboring binary64 lnT values on either side. Record and stop on failure; do not raise the tolerance.',
        scope='Independent re-projection of the already accepted one-step reference state, rebuilding its EOS table and fitting GR constraints. This does not yet re-integrate a strict-root trajectory or certify continuous/physical EOS error.',
        reference_sha256=e.c.sha(e.OUT/'E-1-0-corrector-reference.npz'),
        original_state_sha256=e.c.sha(e.OUT/'E-1-0-corrector-state-17-4.npz')))


def strict_control():
    assert (e.OUT/'strict-entropy-plan.json').exists()
    b=dict(np.load(e.OUT/'initial-state.npz'));eos=e.s.v.ColdEOS();errors=[]
    ROOT_STATS.update(calls=0,evaluations=0,maximum_score=0.,label='strict-entropy-control')
    for i,(lp,t,x) in enumerate(zip(b['logP'],b['lnT'],b['X'])):
        known=eos(1,lp,t,x);_,actual,_=strict_invert(eos,lp,known[3],x,t+(.01 if i%2 else -.01))
        errors.append(abs(actual-t))
    assert max(errors)<1e-10
    e.save('strict-entropy-control.json',dict(classification='Counterexample candidate',passed=True,
        cells=len(errors),known_temperature_maximum_log_error=float(max(errors)),**ROOT_STATS,
        scope='All actual initial P,X states with fresh known-target entropy and alternating perturbed guesses. It tests root recovery, not the EOS physical error.'))
    print('STRICT ENTROPY CONTROL',ROOT_STATS,max(errors),flush=True)


def strict_projection():
    assert json.loads((e.OUT/'strict-entropy-control.json').read_text())['passed']
    plan=json.loads((e.OUT/'strict-entropy-plan.json').read_text());label='strict-E-1-0-corrector'
    path=e.OUT/'E-1-0-corrector-reference.npz';assert e.c.sha(path)==plan['reference_sha256']
    reference=dict(np.load(path));base=dict(np.load(e.OUT/'initial-state.npz'))
    prior=json.loads((e.OUT/'path-1.json').read_text());old=dict(np.load(e.OUT/'E-1-0-corrector-state-17-4.npz'))
    parameters=prior['records'][0]['corrector']['parameters'];e.configure()
    ROOT_STATS.update(calls=0,evaluations=0,maximum_score=0.,label='strict-entropy-projection')
    original=e.c.be.invert
    try:
        e.c.be.invert=strict_invert
        state,record=e.s.structure(label,reference,parameters)
    finally: e.c.be.invert=original
    assert np.array_equal(state['X'],old['X']) and np.array_equal(state['dm'],old['dm'])
    r0=base['radius_faces_m'][0]*100;r1=state['radius_faces_m'][0]*100
    work=np.exp(float(base['boundary_logP']))*4*np.pi*(r1-r0)*(r1*r1+r1*r0+r0*r0)/3
    expected=prior['expected_energy_erg']+prior['pressure_work_erg']-work
    measured=e.c.stable_mass_change(base,state)*(e.c.gr.C*100)**2
    score=float(abs(measured-expected)/prior['release_erg'])
    e.save('strict-entropy-projection.json',dict(classification='Counterexample candidate',completed=True,
        structure=record,root_statistics=ROOT_STATS,energy_score=score,energy_passed=score<1e-6,
        mass_change_energy_erg=float(measured),expected_energy_erg=float(expected),
        previous_energy_score=prior['energy_score'],previous_energy_passed=prior['energy_passed'],
        lnT_difference=float(abs(state['lnT']-old['lnT']).max()),
        mass_difference_energy_erg=float(e.c.stable_mass_change(old,state)*(e.c.gr.C*100)**2),
        full_strict_root_time_path_recomputed=False,physical_or_continuous_certificate=False))
    print('STRICT ENTROPY PROJECTION',score,ROOT_STATS,flush=True)


def transport_control():
    """Frozen interior diffusion: conservation, dissipation and time scale."""
    b=dict(np.load(e.OUT/'initial-state.npz'));basis=dict(np.load(e.OUT/'basis.npz'))
    _,p=e.c.mesa(e.OUT/'initial-profile.data.gz')
    dm=b['dm'];N=np.exp(b['nu']);q=np.exp(4*(b['lnT']+b['nu']))
    weight=dm[:-1]/(dm[:-1]+dm[1:])
    kap=(1-weight)*p['opacity'][:-1]+weight*p['opacity'][1:]
    area=4*np.pi*(b['radius_faces_m'][1:-1]*100)**2
    B=area**2*(4*5.670400e-5)/(3*kap*np.exp(b['nu_faces'][1:-1])**2)/((dm[:-1]+dm[1:])/2)
    w=dm*N*basis['capacity'];a=4*q/w
    assert np.all(B>0) and np.all(w>0) and np.all(np.isfinite(a))
    degree=np.r_[B,0]+np.r_[0,B]
    diagonal=-a*degree;off=B*np.sqrt(a[:-1]*a[1:])
    fastest=float(eigvalsh_tridiagonal(diagonal,off,select='i',select_range=(0,0))[0])
    # With frozen opacity/geometry, dq/dt=-diag(a) K q. Similar to
    # -sqrt(a) K sqrt(a), with K the positive weighted path Laplacian.
    # No opacity, geometry, composition or atmosphere Jacobian is included.
    rng=np.random.default_rng(3201);z=rng.normal(size=len(dm))
    flux=-B*np.diff(z);div=np.r_[flux,0]-np.r_[0,flux]
    identity=abs(float(z@div)-float(B@(np.diff(z)**2)))/float(B@(np.diff(z)**2))
    balance=abs(float(div.sum()))/float(abs(div).sum())
    assert identity<1e-14 and balance<1e-14,(identity,balance)
    # Independent form of the same face law; cancellation is avoided in the
    # production expm1 form, so compare with a floating reconstruction budget.
    actual=e.s.baryon_face_diffusion(b,p['opacity'])
    direct=B*(q[1:]-q[:-1])
    budget=64*np.finfo(float).eps*B*(q[1:]+q[:-1])
    score=float(np.max(abs(actual-direct)/budget));assert score<1,score
    tolman={**b,'lnT':-b['nu']}
    assert np.all(e.s.baryon_face_diffusion(tolman,p['opacity'])==0)
    paths=['star/private/opacities.f90','kap/public/kap_lib.f90',
           'kap/private/kap_eval.f90','star/defaults/controls.defaults']
    hashes={}
    for rel in paths:
        src=e.c.fresh.MESA/rel;dest=e.OUT/'sources'/rel;dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(src,dest);hashes[rel]=e.c.sha(dest)
    source=(e.OUT/'sources/kap/private/kap_eval.f90').read_text()
    defaults=(e.OUT/'sources/star/defaults/controls.defaults').read_text()
    assert 'kap = 1d0 / (1d0/kap_rad + 1d0/kap_ec)' in source
    assert 'include_electron_conduction = .true.' in defaults
    dt=float(np.min(-1/diagonal));cell=int(np.argmin(-1/diagonal))
    e.save('transport-control.json',dict(classification='Counterexample candidate',passed=True,
        initial_profile_sha256=e.c.sha(e.OUT/'initial-profile.data.gz'),source_sha256=hashes,
        frozen_interior_diffusion=True,conservation_relative_residual=balance,
        dissipation_identity_relative_residual=identity,face_reconstruction_score=score,
        constant_redshifted_temperature_zero_flux=True,
        forward_Euler_positivity_step_seconds=dt,limiting_cell=cell,
        fastest_frozen_eigenvalue_per_second=fastest,linear_Euler_stability_step_seconds=-2/fastest,
        registered_step_to_positivity_step=float(json.loads((e.OUT/'plan.json').read_text())['duration_coordinate_seconds']/dt),
        initial_maximum_log10_temperature=float(b['lnT'].max()/np.log(10)),
        opacity_interpretation='The source default includes electron conduction and combines inverse radiative and conductive opacities. Source inspection is not a new compiled-runtime proof. Do not add a second conductive flux to this total-opacity control.',
        boundary='Both boundary perturbation fluxes are fixed to zero for this interior operator. No convection, radiating atmosphere, opacity/metric derivative or physical opacity certificate.',
        full_transport_or_GR_closed=False))
    np.savez_compressed(e.OUT/'transport-control.npz',face_coefficient=B,capacity_weight=w,
        redshifted_temperature_fourth=q,symmetric_diagonal=diagonal,symmetric_off_diagonal=off)
    print('TRANSPORT CONTROL',dt,cell,fastest,score,flush=True)


def stages():
    """Recompute energy and input bindings from completed paths, not summaries."""
    plan=json.loads((e.OUT/'plan.json').read_text());base=dict(np.load(e.OUT/'initial-state.npz'))
    basis=dict(np.load(e.OUT/'basis.npz'));rest=(e.c.W/e.c.A-1)*(e.c.gr.C*100)**2
    flux=np.load(e.s.OLD/'registered-face-flux.npz')['canonical'];rows=[];previous=None
    for count in plan['step_counts']:
        path=e.OUT/f'path-{count}.json'
        if not path.exists(): continue
        report=json.loads(path.read_text());assert report['completed']
        current=base;release=0.;details=[]
        for step in range(count):
            prefix=f'E-{count}-{step}'
            for role in ['predictor','corrector']:
                label=prefix+'-'+role;ref=dict(np.load(e.OUT/f'{label}-reference.npz'))
                state=dict(np.load(e.OUT/f'{label}-state-17-4.npz'))
                binding=json.loads((e.OUT/f'{label}-projection-binding.json').read_text())
                assert binding['reference_binding']==e.fingerprint(ref)
                assert binding['state_sha256']==e.c.sha(e.OUT/f'{label}-state-17-4.npz')
                assert np.array_equal(state['dm'],base['dm']) and np.array_equal(state['X'],ref['X'])
                assert np.array_equal(state['chart_total'],ref['chart_total'])
                reference=json.loads((e.OUT/f'{label}-reference.json').read_text())
                assert reference['previous_physical_binding']==e.fingerprint(current)
                assert reference['passed'] and reference['maximum_inverse_score']<=1
                struct=json.loads((e.OUT/f'{label}-structure-17-4.json').read_text())
                assert struct['interface_max']<1e-8
                details.append(dict(label=label,inverse_score=reference['maximum_inverse_score'],
                    interface=struct['interface_max'],minimum_X=float(state['X'].min()),
                    baryon_gauge_residual=float(abs(state['X'].sum(1)-1).max()),
                    entropy_preservation_maximum=float(abs(state['s_B']-ref['s_B']).max())))
                if role=='predictor' or step+1<count:
                    source=json.loads((e.OUT/f'{label}-chart.json').read_text())
                    assert source['passed'] and source['state_binding']==e.fingerprint(ref,state)
                    assert source['native_source_sha256']==e.c.sha(e.OUT/f'{label}-corrected.npz')
            dx=state['X'].astype(np.longdouble)-current['X'].astype(np.longdouble)
            release+=float(base['dm'].astype(np.longdouble)@abs(dx@rest.astype(np.longdouble)))
            current=state
        loss=float(base['dm'].astype(np.longdouble)@(basis['capacity']*state['chart_total'][:,27]))
        r0=base['radius_faces_m'][0]*100;r1=state['radius_faces_m'][0]*100
        work=np.exp(float(base['boundary_logP']))*4*np.pi*(r1-r0)*(r1*r1+r1*r0+r0*r0)/3
        expected=-flux[0]*plan['duration_coordinate_seconds']-loss-work
        measured=e.c.stable_mass_change(base,state)*(e.c.gr.C*100)**2
        score=float(abs(measured-expected)/release)
        assert release==report['release_erg'] and score==report['energy_score']
        assert report['energy_passed']==(score<plan['controls']['global_energy_relative_to_release'])
        row=dict(steps=count,energy_score=score,energy_passed=report['energy_passed'],stages=details)
        if previous is not None:
            budget=1e-16+1e-3*abs(state['X']-base['X']);diff=state['X']-previous['X']
            ix=np.unravel_index(np.argmax(abs(diff)/budget),diff.shape)
            accumulated=state['chart_total'][:,:26]-previous['chart_total'][:,:26]
            row['refinement']=dict(composition_score=float(np.max(abs(diff)/budget)),
                accumulated_composition_score=float(np.max(abs(accumulated)/budget)),
                worst_cell=int(ix[0]),worst_species=e.c.NAMES[ix[1]],
                worst_difference=float(diff[ix]),worst_budget=float(budget[ix]),
                lnT_difference=float(abs(state['lnT']-previous['lnT']).max()))
        rows.append(row);previous=state
    e.save('stage-audit.json',dict(classification='Counterexample candidate',rows=rows,
        complete=len(rows)==len(plan['step_counts']),physical_or_continuous_certificate=False))
    print('STAGE AUDIT',[(r['steps'],r['energy_score'],r.get('refinement')) for r in rows],flush=True)


def implicit_transport_control():
    """Solve the specified frozen-coefficient interior control, not a star."""
    b=dict(np.load(e.OUT/'transport-control.npz'));q=b['redshifted_temperature_fourth']
    B=b['face_coefficient'];a=4*q/b['capacity_weight'];n=len(q)
    duration=json.loads((e.OUT/'plan.json').read_text())['duration_coordinate_seconds']
    rows=[];scale=q.max();initial=q/scale
    for dt in [duration/4,duration/2,duration]:
        ab=np.zeros((3,n));ab[1]=1+dt*a*(np.r_[B,0]+np.r_[0,B])
        ab[0,1:]=-dt*a[:-1]*B;ab[2,:-1]=-dt*a[1:]*B
        final=solve_banded((1,1),ab,initial)
        uniform=solve_banded((1,1),ab,np.ones(n))
        ld=np.longdouble;edge=B.astype(ld)*(final[:-1].astype(ld)-final[1:].astype(ld))
        residual=final.astype(ld)+dt*a.astype(ld)*(np.r_[edge,ld(0)]-np.r_[ld(0),edge])-initial.astype(ld)
        conservation=abs(np.sum((final.astype(ld)-initial.astype(ld))/a.astype(ld)))/np.sum(initial.astype(ld)/a.astype(ld))
        constant=float(max(abs(uniform-1)))
        minimum=float(final.min());maximum=float(final.max())
        assert minimum>=initial.min()*(1-1e-6) and maximum<=initial.max()*(1+1e-6)
        assert conservation<1e-6 and constant<1e-6,(conservation,constant)
        rows.append(dict(step_seconds=dt,weighted_invariant_relative_defect=float(conservation),
            constant_state_maximum_defect=constant,maximum_normalized_equation_residual=float(abs(residual).max()),
            temperature_log_change_range=[float(.25*np.log(final/initial).min()),float(.25*np.log(final/initial).max())],
            minimum_normalized_q=minimum,maximum_normalized_q=maximum))
    e.save('implicit-transport-control.json',dict(classification='Counterexample candidate',passed=True,rows=rows,
        meaning='An O(n) band solve on the actual frozen 5735-cell coefficients. It preserves the linear weighted invariant and positivity at the listed steps. The invariant is not the nonlinear physical stellar energy.',
        missing='Temperature/EOS capacity, opacity, GR metric and composition feedback, nonzero surface transport, convection and atmosphere. No full stellar evolution or physical conservation verdict.'))
    print('IMPLICIT TRANSPORT CONTROL',rows,flush=True)


def embedding_moment_control():
    """Expose an exact non-ideal boundary of the positive element embedding."""
    import sympy as sy
    lo,z,hi=sy.symbols('lo z hi',real=True)
    variance=(hi-z)/(hi-lo)*lo**2+(z-lo)/(hi-lo)*hi**2-z**2
    assert sy.factor(variance-(z-lo)*(hi-z))==0
    matrix=e.c.embedding();excess=matrix@(e.c.EZ**2)-e.c.Z**2
    missing=~np.isin(e.c.Z,e.c.EZ)
    assert np.all(excess[missing]>0) and np.all(excess[~missing]==0)
    coefficient=np.where(missing,excess,0);rows=[]
    labels=['initial-state']+[f'E-{n}-{n-1}-corrector-state-17-4' for n in [1,2,4]]
    for label in labels:
        path=e.OUT/(label+'.npz')
        if not path.exists(): continue
        b=dict(np.load(path));y=b['X']/e.c.A
        delta=y@coefficient;real=y@(e.c.Z**2);ratio=delta/real;i=int(ratio.argmax())
        rows.append(dict(label=label,state_sha256=e.c.sha(path),maximum_fractional_second_moment_excess=float(ratio[i]),
            worst_cell=i,maximum_absolute_moment_excess=float(delta.max()),
            affected_cells=int((delta>0).sum()),
            maximum_excess_fraction_of_classical_screening_sum=float(np.max(delta/(real+y@e.c.Z)))))
    hashes={}
    for rel in ['src/master_coulomb.f90','src/ionize.f90','src/mod_free_eos.f90']:
        target=e.OUT/'sources/free_eos'/rel;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(e.c.gr.SOURCE/rel,target);hashes[rel]=e.c.sha(target)
    bridge=e.ROOT/'verification/gr_eos_bridge.f90'
    assert 'call free_eos(0,3,1,-2,kif' in bridge.read_text()
    src=(e.c.gr.SOURCE/'src/mod_free_eos.f90').read_text().split('elseif(ifoption.eq.3) then',1)[1]
    branch=src.split('elseif(ifmodified.eq.1) then',1)[1].split('elseif',1)[0]
    assert 'ifcoulomb = 5' in branch
    e.save('embedding-moment-control.json',dict(classification='Proven',passed=True,
        statement='No positive representation of an absent nuclear charge Z by other charges can preserve ion number, mean charge and charge squared simultaneously. Unit weights and preserved mean imply E[(Zprime-Z)^2]=E[Zprime^2]-Z^2; a zero value requires every nonzero weight at Z itself.',
        two_neighbor_excess='(Z-Zlo)*(Zhi-Z)',
        unsupported_species_excess={name:float(excess[i]) for i,name in enumerate(e.c.NAMES) if missing[i]},
        source_sha256=hashes,rows=rows,bridge_source_sha256=e.c.sha(bridge),
        declared_EOS_options=dict(ifoption=3,ifmodified=1,ifion=-2,model='EOS1',ifcoulomb=5),
        consequence='Even in the completely ionized limit the present embedding changes the ion charge-squared input to Coulomb free energy. A rho,T-independent atomic energy offset cannot correct a rho,T-dependent Coulomb term. This disproves exact non-ideal EOS equivalence; it is not a lower or upper bound on the physical stellar pressure error.',
        scope='The theorem concerns positive composition representations and the stated charge moments. Saved-state numbers evaluate the fully stripped moments only; they do not assert actual complete ionization, a Debye-Huckel validity regime or a physical EOS error certificate.',
        primary_reference='https://freeeos.sourceforge.net/coulomb.pdf',
        physical_EOS_certified=False))
    print('EMBEDDING MOMENT CONTROL',rows,flush=True)


def predictor_comparison():
    """Quantify the old trial's substitution; never reuse its state here."""
    base=dict(np.load(e.OUT/'initial-state.npz'));rows=[]
    for count in [1,2,4]:
        path=e.OUT/f'E-{count}-0-predictor-state-17-4.npz'
        if not path.exists(): continue
        old_path=e.OLD/f'P-{count}-0-state-17-4.npz'
        new=dict(np.load(path));old=dict(np.load(old_path))
        energy=e.c.stable_mass_change(old,new)*(e.c.gr.C*100)**2
        score=float(np.max(abs(new['X']-old['X'])/(1e-16+1e-3*abs(new['X']-base['X']))))
        rows.append(dict(steps=count,new_state_sha256=e.c.sha(path),old_state_sha256=e.c.sha(old_path),
            composition_score=score,lnT_difference=float(abs(new['lnT']-old['lnT']).max()),
            lnd_difference=float(abs(new['lnd']-old['lnd']).max()),
            maximum_midpoint_radius_difference_m=float(abs(new['r_mid_m']-old['r_mid_m']).max()),
            mass_function_difference_energy_erg=float(energy)))
    e.save('predictor-comparison.json',dict(classification='Counterexample candidate',rows=rows,
        complete=len(rows)==3,
        meaning='Direct predictors versus the old source-first GR predictors at the same first-step duration. Tiny state differences can be material to an energy residual close to its gate. Every Request32 predictor was independently solved; these differences do not certify equivalence of the original substituted trial.'))
    print('PREDICTOR COMPARISON',rows,flush=True)


def inverse_energy_audit():
    """Separate actual finite inversion residuals from a passed energy gate."""
    base=dict(np.load(e.OUT/'initial-state.npz'));rows=[];eos=e.s.v.ColdEOS()
    for count in [1,2,4]:
        path=e.OUT/f'path-{count}.json'
        if not path.exists(): continue
        report=json.loads(path.read_text());release=report['release_erg'];details=[]
        for step in range(count):
            label=f'E-{count}-{step}-corrector'
            ref=dict(np.load(e.OUT/f'{label}-reference.npz'))
            state=dict(np.load(e.OUT/f'{label}-state-17-4.npz'))
            table=dict(np.load(e.OUT/f'{label}-adiabats-17.npz'))
            actual=np.array([eos(1,lp,t,x) for lp,t,x in zip(ref['logP'],ref['lnT'],ref['X'])])
            dh=actual[:,2]+actual[:,1]/actual[:,0]-ref['reference_enthalpy']
            weight=base['dm']*np.exp(state['nu']);temperature=np.exp(state['lnT'])
            terms=dict(reference_enthalpy_inverse=weight*np.exp(state['lnT']-ref['lnT'])*dh,
                reference_entropy_reconstruction=weight*temperature*(table['reference'][:,3]-actual[:,3]),
                structural_entropy_inverse=weight*temperature*(state['s_B']-table['reference'][:,3]))
            details.append(dict(label=label,terms={key:dict(signed_energy_estimate_erg=float(value.sum()),
                unsigned_energy_estimate_erg=float(abs(value).sum()),unsigned_relative_to_release=float(abs(value).sum()/release))
                for key,value in terms.items()}))
        rows.append(dict(steps=count,energy_score=report['energy_score'],energy_passed=report['energy_passed'],
            energy_gate_margin_erg=float(1e-6*release-abs(report['mass_change_energy_erg']-report['expected_energy_erg'])),
            steps_detail=details))
    e.save('inverse-energy-audit.json',dict(classification='Counterexample candidate',rows=rows,
        complete=len(rows)==3,
        scope='Actual saved-state finite residuals with linearized temperature/redshift weights. These are neither rigorous propagated mass-error bounds nor independent terms guaranteed to add to the observed mass defect. A positive gate margin does not certify robustness to EOS/structure error.',
        reference_entropy_gauge='The inherited structural inverse uses a relative entropy residual; it is not an energy-unit root error certificate.',
        continuous_or_physical_certificate=False))
    print('INVERSE ENERGY AUDIT',rows,flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
