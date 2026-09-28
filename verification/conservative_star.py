"""Request31: corrected source derivatives before conservative GR continuation.

Counterexample candidate. A new forced quasistatic evolution must be re-computed;
neither finite derivative checks nor an EOS seed certify a physical star.
"""
from pathlib import Path
import ctypes, json, shutil, sys
import numpy as np
import closure_precision as v

c=v.c
ROOT=v.ROOT
OLD=v.OUT
OUT=ROOT/'outputs/conservative-star31'
CACHE=Path('/home/lpaiu/work/conservative-star31')


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists() and not CACHE.exists();OUT.mkdir();CACHE.mkdir()
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='d7a8331',
        previous_manifest_sha256=c.sha(OLD/'manifest.json'),
        objective='Continue the full original physical EOS, GR evolution, scalar driving and observational-inference goal. This stage closes source derivative and conservative time-integration prerequisites without redefining that goal.',
        ordered_work=['All 25 independent composition directions of the same corrected 26-species source',
            'Common-EOS density and temperature directions including neutrino loss',
            'Increment accumulation with extended precision and an explicit conservation gauge',
            'New fixed-baryon GR time paths, original energy/composition/temperature refinement criteria',
            'Self-consistent transport, physical EOS/GR, scalar and observation closure remain required'],
        composition=dict(pivot='he4',steps=[1e-6,5e-7],directions=[n for n in c.NAMES if n!='he4'],
            method='One-sided simplex directions e_i-e_He4; second-order extrapolation of two step sizes. '
                'Compare each with the actual returned Jacobian and record discrepancies. Do not silently promote a finite derivative to a rigorous certificate.'),
        density_temperature_steps=[1e-4,5e-5],
        controls=dict(source_derivative_relative=1e-3,composition_absolute=1e-16,
            composition_relative_to_total_change=1e-3,temperature_log_absolute=2e-6,
            global_energy_relative_to_release=1e-6),
        time=dict(coordinate_seconds=140780.16,step_counts=[1,2,4],
            flux='Exact registered Request30 face luminosities for an explicit forced control. A later transport closure must replace this boundary prescription.'),
        policy='Use the unchanged native binary and the declared Request30 double-table correction. '
            'Use seeded FreeEOS for all fresh auxiliary inputs. Preserve all prior raw verdicts and any new failed directions. '
            'Do not equate same-time endpoint inverses or fixed-flux TOV sequences with full physical GR evolution.'))
    bindings={}
    for rel in ['docs/'+n+'.md' for n in c.cell.DOCS]+['paper/revision-manifest.json']:
        p=ROOT/rel;target=OUT/'previous-notes'/p.name;target.parent.mkdir(exist_ok=True);shutil.copy2(p,target)
        bindings[rel]=dict(snapshot=target.relative_to(ROOT).as_posix(),sha256=c.sha(p),historical_manifest=rel.startswith('paper/'))
    save('historical-note-bindings.json',bindings)
    shutil.copy2(v.OLD/'restored-state-17-4.npz',OUT/'initial-state.npz')


def auxiliary(data):
    lib=ctypes.CDLL(str(c.CACHE/'common_eos_aux.so'));call=lib.common_eos_aux
    call.argtypes=c.gr.EOS().call.argtypes;call.restype=None;eos=v.ColdEOS()
    def raw(r,t,x):
        # Same declared fixed seed, then the unchanged 22-output bridge.
        seed=np.zeros(26);seed[c.NAMES.index('c12') if x[c.NAMES.index('c12')]<.5 else c.NAMES.index('h1')]=1
        c.EOS.__call__(eos,0,-20.,np.log(1e6),seed)
        _,ym=eos.composition(x);cx=ym@eos.weights;eps=np.ascontiguousarray(ym/cx)
        result=np.full(22,np.nan);info=ctypes.c_int(-999)
        call(2,float(r+np.log(cx)),float(t),eps,result,ctypes.byref(info))
        assert info.value==0 and np.all(np.isfinite(result));return result
    rows=[];h=5e-5
    for r,t,x in zip(data['lnd'],data['lnT'],data['X']):
        a=raw(r,t,x);at=(raw(r,t+h,x)[12]-raw(r,t-h,x)[12])/(2*h)
        ar=(raw(r+h,t,x)[12]-raw(r-h,t,x)[12])/(2*h)
        rows.append([a[13]/np.exp(r),a[12],at,ar])
    return np.array(rows)


def evaluate(label,data):
    c.native.OUT=OUT;c.native.CACHE=CACHE;c.native.context()
    v.shim.OUT=OUT;v.shim.CACHE=CACHE;v.shim.save=save
    path=OUT/f'{label}-native.npz'
    if not path.exists():
        c.native.setup(label,data,species=c.NAMES,network=(c.cell.OLD/'inputs/explicit-PP/cno_extras.net').read_text())
        aux=auxiliary(data);np.save(OUT/f'{label}-replacement.npy',aux)
        v.shim.trace(label,aux)
    inp=np.load(OUT/f'{label}-input.npz')
    for key in ['lnd','lnT','X','dm']: assert np.array_equal(inp[key],data[key]),(label,key)
    native=dict(np.load(path));result=v.corrected_values(native)
    np.savez_compressed(OUT/f'{label}-corrected.npz',**result)
    return native,result


def derivatives():
    plan=json.loads((OUT/'plan.json').read_text());base=dict(np.load(OUT/'initial-state.npz'))
    native,reference=evaluate('initial',base);active=abs(reference['heat'])>1;pivot=c.NAMES.index(plan['composition']['pivot'])
    old=np.load(OLD/'unified-source-audit.npz')
    errors={key:float(np.max(abs(reference[key]-old[key])/np.maximum(1e-30,abs(old[key])))) for key in ['dxdt','heat','neutrino']}
    save('initial-control.json',dict(classification='Counterexample candidate',same_state_corrected_relative_errors=errors,
        seeded_auxiliaries_used=True,inputs_bound=True,physical_EOS_certified=False))
    rows=[]
    for species in plan['composition']['directions']:
        j=c.NAMES.index(species);sides=[]
        for h in plan['composition']['steps']:
            data={k:v.copy() for k,v in base.items()};data['X'][:,j]+=h;data['X'][:,pivot]-=h
            assert data['X'].min()>=0
            _,result=evaluate(f'X-{species}-{h}',data)
            sides.append({key:(result[key]-reference[key])/h for key in reference})
        direction={key:2*sides[1][key]-sides[0][key] for key in reference}
        reported=native['jacobian'][:,:,j]-native['jacobian'][:,:,pivot]
        scale=np.maximum(1e-30,np.maximum(abs(direction['dxdt']).max(1),abs(reported).max(1)))
        disagreement=abs(direction['dxdt']-reported).max(1)/scale
        refinement=abs(sides[1]['dxdt']-sides[0]['dxdt']).max(1)/scale
        row=dict(species=species,returned_Jacobian_score=float(max(disagreement[active])),
            finite_step_difference_score=float(max(refinement[active])),
            returned_Jacobian_passed=bool(max(disagreement[active])<1e-3),
            finite_step_passed=bool(max(refinement[active])<1e-3),
            minimum_X=float(base['X'][:,j].min()))
        rows.append(row);np.savez_compressed(OUT/f'X-{species}-derivative.npz',**direction)
        save('composition-derivatives.json',dict(classification='Counterexample candidate',rows=rows,
            complete=len(rows)==len(plan['composition']['directions']),
            full_domain_rigorous_certificate=False,GR_trajectory_recomputed=False))
        print('COMPOSITION',row,flush=True)


def precompute_tail():
    """Evaluate disjoint late inputs; the running serial audit owns all verdicts."""
    plan=json.loads((OUT/'plan.json').read_text());base=dict(np.load(OUT/'initial-state.npz'))
    species=plan['composition']['directions'][-6:];pivot=c.NAMES.index(plan['composition']['pivot'])
    assert species==['ne18','ne19','ne20','mg22','mg24','ca40']
    save('tail-input-plan.json',dict(classification='Counterexample candidate',species=species,
        original_plan_sha256=c.sha(OUT/'plan.json'),steps=plan['composition']['steps'],
        purpose='Precompute the exact later inputs under disjoint labels. The original serial audit alone writes derivative verdicts. No tolerance, source function or input changes.'))
    for name in species:
        for h in plan['composition']['steps']:
            data={k:a.copy() for k,a in base.items()};data['X'][:,c.NAMES.index(name)]+=h;data['X'][:,pivot]-=h
            evaluate(f'X-{name}-{h}',data)
    save('tail-input-complete.json',dict(classification='Counterexample candidate',species=species,complete=True))


def evaluate_replaced(label,data,values):
    """Declared seven-input control on the same unmodified native executable."""
    c.native.OUT=OUT;c.native.CACHE=CACHE;c.native.context();v.shim.OUT=OUT;v.shim.CACHE=CACHE;v.shim.save=save
    if not (OUT/f'{label}-native.npz').exists():
        c.native.setup(label,data,species=c.NAMES,network=(c.cell.OLD/'inputs/explicit-PP/cno_extras.net').read_text())
        np.save(OUT/f'{label}-replacement.npy',values);v.shim.trace(label,values)
    inp=np.load(OUT/f'{label}-input.npz')
    for key in ['X','lnd','lnT','dm']: assert np.array_equal(inp[key],data[key]),(label,key)
    result=dict(np.load(OUT/f'{label}-native.npz'))
    assert np.array_equal(result['moments_used'],values[:,:3]) and np.array_equal(result['aux_used'],values[:,3:])
    return result


def composition_partial_control():
    plan=json.loads((OUT/'composition-partial-plan.json').read_text());base=dict(np.load(OUT/'initial-state.npz'))
    _,p=c.mesa(OUT/'initial-profile.data.gz');aux=np.load(OUT/'initial-replacement.npy')
    fixed=np.c_[p['abar'],p['zbar'],p['z2bar'],aux];reference=dict(np.load(OUT/'initial-native.npz'))
    control=evaluate_replaced('partial-initial',base,fixed)
    same={key:bool(np.array_equal(control[key],reference[key])) for key in ['dxdt','heat','neutrino','jacobian']}
    save('composition-partial-baseline.json',dict(classification='Counterexample candidate',bitwise_equal=same,passed=all(same.values())))
    assert all(same.values()),same
    rows=[];pivot=c.NAMES.index('he4');active=abs(reference['heat'])>1
    for name in plan['directions']:
        j=c.NAMES.index(name);values=[]
        for h in plan['steps']:
            data={k:a.copy() for k,a in base.items()};data['X'][:,j]+=h;data['X'][:,pivot]-=h
            result=evaluate_replaced(f'partial-{name}-{h}',data,fixed);values.append((result['dxdt']-reference['dxdt'])/h)
        derivative=2*values[1]-values[0];reported=reference['jacobian'][:,:,j]-reference['jacobian'][:,:,pivot]
        scale=np.maximum(1e-30,np.maximum(abs(derivative).max(1),abs(reported).max(1)))
        score=float(max((abs(derivative-reported).max(1)/scale)[active]))
        rows.append(dict(species=name,fixed_auxiliary_partial_score=score,passed=score<plan['criterion']))
        np.savez_compressed(OUT/f'partial-{name}-derivative.npz',dxdt=derivative)
        save('composition-partial-control.json',dict(classification='Counterexample candidate',rows=rows,
            complete=len(rows)==len(plan['directions']),passed=all(r['passed'] for r in rows),
            baseline_bitwise_equal=True,full_composed_source_derivative_certified=False))
        print('FIXED AUXILIARY PARTIAL',rows[-1],flush=True)


def composition_chain_boundary():
    import sympy as sy
    result=json.loads((OUT/'composition-partial-control.json').read_text())
    assert result['complete'] and result['passed'] and result['baseline_bitwise_equal']
    paths=['net/public/net_lib.f90','net/private/net_derivs.f90','net/private/net_screen.f90']
    for rel in paths:
        target=OUT/'sources'/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(c.fresh.MESA/rel,target)
    source=(OUT/'sources/net/public/net_lib.f90').read_text()
    assert 'abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho' in source
    labels=['initial']+[f'X-{name}-{h}' for name in ['h1','n14'] for h in [1e-6,5e-7]]
    exact={}
    for label in labels:
        raw=dict(np.load(OUT/f'{label}-native.npz'));value=v.corrected_values(raw)
        exact[label]=bool(np.array_equal(raw['dxdt'],value['dxdt']))
    save('composition-zero-correction-hypothesis.json',dict(classification='Counterexample candidate',
        raw_and_corrected_vectors_bitwise_equal=exact,passed=all(exact.values()),
        meaning='A zero value correction is not required. The relevant fixed-auxiliary directional correction must be constant when only H1/N14 and He4 are varied.'))
    baseline=dict(np.load(OUT/'partial-initial-native.npz'))
    before=v.weak_vector(baseline,False)[0]-v.weak_vector(baseline,True)[0];constant={}
    for label in [f'partial-{name}-{h}' for name in ['h1','n14'] for h in [1e-6,5e-7]]:
        raw=dict(np.load(OUT/f'{label}-native.npz'))
        delta=v.weak_vector(raw,False)[0]-v.weak_vector(raw,True)[0]
        constant[label]=bool(np.array_equal(delta,before))
    assert all(constant.values()),constant
    x,z=sy.symbols('x z');a=sy.Function('a')(x);F=x*x*z+sy.exp(x*z)
    expression=sy.diff(F.subs(z,a),x)
    expected=sy.diff(F,x).subs(z,a)+sy.diff(F,z).subs(z,a)*sy.diff(a,x)
    assert sy.simplify(expression-expected)==0
    save('composition-chain-boundary.json',dict(classification='Counterexample candidate',
        baseline_seven_input_intervention_bitwise_equal=True,fixed_auxiliary_partial_checks=result['rows'],
        weak_vector_correction_constant_in_fixed_auxiliary_controls=constant,
        conclusion='For these H1 and N14 controls, the returned native fixed-auxiliary partial is supported. Its difference from the composed source derivative is an auxiliary-chain contribution, not evidence that this fixed-auxiliary API partial itself is wrong. This finding does not repair the separately observed weak temperature/density derivative omissions.',
        next_derivative_requirement='Differentiate the independent composition/EOS auxiliaries and assemble the full chain, or certify a direct composed derivative with a suitable boundary treatment. The failed all-zone H2 extrapolation is retained.',
        source_sha256={rel:c.sha(OUT/'sources'/rel) for rel in paths}))
    save('composition-chain-symbolic.json',dict(classification='Proven',chain_rule=True,nonlinear_symbolic_control=True,
        identity='D_X F(X,a(X)) = partial_X F + partial_a F D_X a',
        scope='Differentiable composed functions on their stated domain. No assertion of a globally differentiable EOS at all zero-abundance boundaries.'))
    print('PASS fixed-auxiliary composition boundary and chain identity',flush=True)


def accumulate(initial,total,delta):
    """Keep physical increments apart from the rounded native EOS input."""
    assert initial.shape==total.shape==delta.shape
    assert np.all(np.isfinite(delta)) and total.dtype==np.longdouble
    zones=np.arange(len(initial));pivot=np.argmax(initial,axis=1)
    addition=delta.astype(np.longdouble)
    correction=-addition.sum(1,dtype=np.longdouble);addition[zones,pivot]+=correction
    updated=total+addition;physical=initial.astype(np.longdouble)+updated
    minimum=float(physical.min());assert minimum>=-1e-22,('negative abundance',minimum)
    projection=np.maximum(-physical,0);updated+=projection
    updated[zones,pivot]-=projection.sum(1,dtype=np.longdouble)
    inputs=np.asarray(initial.astype(np.longdouble)+updated,dtype=float)
    assert inputs.min()>=0 and max(abs(updated.sum(1,dtype=np.longdouble)))<1e-16
    return inputs,updated,dict(minimum_unprojected=minimum,
        largest_increment_mass_correction=float(max(abs(correction))),largest_positive_projection=float(projection.max()),
        maximum_accumulated_baryon_error=float(max(abs(updated.sum(1,dtype=np.longdouble)))))


def compensation_check():
    assert np.finfo(np.longdouble).nmant>np.finfo(float).nmant
    initial=np.array([[.8,.2]]);total=np.zeros_like(initial,dtype=np.longdouble);naive=initial.copy()
    delta=np.array([[-1e-19,1e-19]])
    for _ in range(1000):
        native,total,row=accumulate(initial,total,delta);naive+=delta
    expected=1000*delta.astype(np.longdouble)
    assert np.array_equal(naive,initial)
    assert np.max(abs(total-expected))<np.longdouble('1e-32')
    assert np.max(abs(native-initial))>0 and np.max(abs(total.sum(1)))==0
    from scipy.linalg import expm
    rate=.7;dt=2.;J=np.array([[-rate,0],[rate,0]]);aug=np.zeros((3,3));aug[:2,:2]=J;aug[:2,2]=J@initial[0]
    increment=expm(dt*aug)[:2,2][None,:]
    answer,_,_=accumulate(initial,np.zeros_like(initial,dtype=np.longdouble),increment)
    exact=np.array([[.8*np.exp(-rate*dt),1-.8*np.exp(-rate*dt)]])
    assert np.max(abs(answer-exact))<1e-14
    save('compensation-control.json',dict(classification='Counterexample candidate',
        longdouble_mantissa_bits=int(np.finfo(np.longdouble).nmant),float64_mantissa_bits=int(np.finfo(float).nmant),
        tiny_increment_accumulation_error=float(np.max(abs(total-expected))),
        ordinary_addition_loses_all_test_increments=True,analytic_conversion_error=float(np.max(abs(answer-exact))),
        projection=row,passed=True,stellar_time_refinement_passed=False))
    print('PASS extended-precision increments and analytic conversion control',flush=True)


def thermodynamic_derivatives():
    plan=json.loads((OUT/'plan.json').read_text());base=dict(np.load(OUT/'initial-state.npz'))
    native=dict(np.load(OUT/'initial-native.npz'));reference=dict(np.load(OUT/'initial-corrected.npz'))
    active=abs(reference['heat'])>1;rows=[];arrays={}
    for field in ['lnT','lnd']:
        results=[]
        for h in plan['density_temperature_steps']:
            sides=[]
            for sign in [-1,1]:
                data={k:v.copy() for k,v in base.items()};data[field]+=sign*h
                _,result=evaluate(f'{field}-{h}-{sign}',data);sides.append(result)
            derivative={key:(sides[1][key]-sides[0][key])/(2*h) for key in reference}
            results.append(derivative)
            np.savez_compressed(OUT/f'{field}-{h}-derivative.npz',**derivative)
        for key in reference:
            a,b=results;scale=np.maximum(1e-30 if key=='dxdt' else 1,
                np.maximum(abs(a[key]).max(1),abs(b[key]).max(1)) if key=='dxdt' else np.maximum(abs(a[key]),abs(b[key])))
            difference=abs(b[key]-a[key]);score=(difference.max(1) if key=='dxdt' else difference)/scale
            rows.append(dict(direction=field,output=key,step_difference_score=float(max(score[active])),
                passed=bool(max(score[active])<1e-3),worst_active_cell=int(np.argmax(np.where(active,score,-1)))))
            arrays[key+'_'+field]=b[key]
        save('thermodynamic-derivatives.json',dict(classification='Counterexample candidate',rows=rows,
            complete=field=='lnd',physical_or_continuous_error_certified=False))
        np.savez_compressed(OUT/'thermodynamic-derivatives.npz',**arrays)
        print('THERMODYNAMIC',rows[-3:],flush=True)


def eos_directions(data):
    """Fixed-P energy/density derivatives on the same 25-dimensional simplex."""
    eos=v.ColdEOS();pivot=c.NAMES.index('he4');n=len(data['X'])
    e0=np.empty(n);cap=np.empty(n);rt=np.empty(n);ex=np.zeros((n,26));rx=np.zeros_like(ex)
    for i in range(n):
        lp=data['logP'][i];t=data['lnT'][i];x=data['X'][i];a=eos(1,lp,t,x)
        e0[i]=a[2]+a[1]/a[0];cap[i]=a[10]-a[1]/a[0]*a[8];rt[i]=a[8]
        assert cap[i]>0
        for k in range(26):
            if k==pivot: continue
            slopes=[]
            for h in [1e-5,5e-6]:
                z=x.copy();z[k]+=h;z[pivot]-=h;az=eos(1,lp,t,z)
                slopes.append([(az[2]+az[1]/az[0]-e0[i])/h,np.log(az[0]/a[0])/h])
            ex[i,k],rx[i,k]=2*np.array(slopes[1])-np.array(slopes[0])
        if i%1000==0: print('EOS DIRECTIONS',i,flush=True)
    return dict(energy=e0,capacity=cap,rho_lnT=rt,energy_X=ex,rho_X=rx)


def initial_eos_directions():
    result=eos_directions(dict(np.load(OUT/'initial-state.npz')))
    np.savez_compressed(OUT/'initial-eos-directions.npz',**result)
    save('initial-eos-directions.json',dict(classification='Counterexample candidate',cells=len(result['energy']),
        minimum_capacity=float(min(result['capacity'])),pivot='he4',composition_steps=[1e-5,5e-6],
        fixed_pressure=True,physical_EOS_certified=False,continuous_derivative_certified=False))
    print('PASS initial EOS energy/density directions',flush=True)


def tangent_plan():
    save('source-tangent-plan.json',dict(classification='Counterexample candidate',
        reason='The first H1 simplex direction has a 0.03484 native Jacobian discrepancy while the two finite steps agree to 5.84e-7. Use independently measured derivatives, not the known incomplete returned Jacobian.',
        inputs='All 25 composition directions at two sizes, both temperature/density central steps, and thermal-neutrino profiles from those exact runs.',
        derivative_normalization='Vector max-norm per cell with absolute floor 1e-30; scalar heat/neutrino derivatives with absolute floor 1 erg g^-1 s^-1. Keep the existing 1e-3 relative gate.',
        operator='Use the measured initial tangent as an explicit fixed exponential-integrator preconditioner. Reevaluate the full nonlinear source and EOS at each time step. Do not call the fixed preconditioner a current-state exact Jacobian. Step refinement must demonstrate the accuracy of the new trajectory.',
        energy='The coupled temperature equation uses directly differentiated total neutrino loss. It does not infer that derivative by subtracting native heat from a rest-mass source; the latter can omit chemical work.',
        conservation='Carry abundance increments in longdouble, round only native/EOS inputs, and record the increment baryon correction and any positive projection. Evaluate the discrete rest-energy change from the actual rounded endpoint compositions used by the EOS and GR structure.',
        thermal_neutrino_boundary='Check the archived neu_get signature and star caller before including the independently differentiated thermal-neutrino profile.',
        physical_EOS_certified=False,full_GR_dynamics=False))


def thermal_loss(label):
    return c.mesa(OUT/f'{label}-profile.data.gz')[1]['non_nuc_neu']


def collect_tangent():
    assert (OUT/'source-tangent-plan.json').exists()
    comp=json.loads((OUT/'composition-derivatives.json').read_text());thermo=json.loads((OUT/'thermodynamic-derivatives.json').read_text())
    assert comp['complete'] and thermo['complete']
    selected={row['species']:OUT/f"X-{row['species']}-derivative.npz" for row in comp['rows']}
    for row in comp['rows']:
        if not row['finite_step_passed']:
            refined=json.loads((OUT/f"refined-{row['species']}.json").read_text())
            assert refined['finite_step_passed'],refined
            assert refined['original_derivative_sha256']==c.sha(selected[row['species']])
            selected[row['species']]=OUT/f"refined-{row['species']}-derivative.npz"
    chosen_fields={}
    for field in ['lnT','lnd']:
        field_rows=[r for r in thermo['rows'] if r['direction']==field]
        chosen_fields[field]=5e-5
        if not all(r['passed'] for r in field_rows):
            refined=json.loads((OUT/f'refined-{field}.json').read_text())
            assert refined['passed'],refined
            chosen_fields[field]=refined['steps'][-1]
    base=dict(np.load(OUT/'initial-state.npz'));n=len(base['X']);pivot=c.NAMES.index('he4')
    fx=np.zeros((n,26,26));nx=np.zeros((n,26));thermal0=thermal_loss('initial');thermal_rows=[]
    for species in c.NAMES:
        if species=='he4': continue
        j=c.NAMES.index(species);gradient=np.load(selected[species])
        fx[:,:,j]=gradient['dxdt'];values=[]
        refined=selected[species].name.startswith('refined-')
        sizes=[1e-8,5e-9] if refined else [1e-6,5e-7];prefix='X2' if refined else 'X'
        for h in sizes: values.append((thermal_loss(f'{prefix}-{species}-{h}')-thermal0)/h)
        nx[:,j]=gradient['neutrino']+2*values[1]-values[0]
        nuclear0=np.load(OUT/'initial-corrected.npz')['neutrino']
        full=[(np.load(OUT/f'{prefix}-{species}-{h}-corrected.npz')['neutrino']-nuclear0)/h+thermal for h,thermal in zip(sizes,values)]
        active=abs(np.load(OUT/'initial-corrected.npz')['heat'])>1
        loss_score=float(max((abs(full[1]-full[0])/np.maximum(1,np.maximum(abs(full[0]),abs(full[1]))))[active]))
        thermal_rows.append(dict(species=species,maximum_thermal_step_difference=float(np.max(abs(values[1]-values[0]))),
            full_loss_step_difference_score=loss_score,passed=loss_score<1e-3))
    save('loss-composition-check.json',dict(classification='Counterexample candidate',rows=thermal_rows,passed=all(row['passed'] for row in thermal_rows)))
    assert all(row['passed'] for row in thermal_rows),thermal_rows
    h2=json.loads((OUT/'h2-extrapolation.json').read_text())
    h2_used=h2['vector_passed'] and h2['loss_passed']
    if h2_used:
        refined_h2=np.load(OUT/'h2-extrapolation-derivative.npz');ih2=c.NAMES.index('h2')
        fx[:,:,ih2]=refined_h2['dxdt'];nx[:,ih2]=refined_h2['total_loss']
    chain=json.loads((OUT/'h2-chain-tangent.json').read_text())
    assert chain['all_zone_finite_chain_passed'] and chain['burning_direct_passed'],chain
    assert json.loads((OUT/'h2-original-loss-all-zones.json').read_text())['passed']
    fx[:,:,c.NAMES.index('h2')]=np.load(OUT/'h2-chain-tangent.npz')['dxdt']
    arrays=dict(f_X=fx,loss_X=nx)
    for field in ['lnT','lnd']:
        h=chosen_fields[field];d=np.load(OUT/f'{field}-{h}-derivative.npz')
        extra=(thermal_loss(f'{field}-{h}-1')-thermal_loss(f'{field}-{h}--1'))/(2*h)
        arrays['f_'+field]=d['dxdt'];arrays['loss_'+field]=d['neutrino']+extra
    assert np.all(fx[:,:,pivot]==0) and np.all(nx[:,pivot]==0)
    np.savez_compressed(OUT/'source-tangent.npz',**arrays)
    save('source-tangent.json',dict(classification='Counterexample candidate',
        composition_directions=25,fixed_initial_preconditioner=True,thermal_neutrino_rows=thermal_rows,
        H2_separate_all_zone_extrapolation_sha256=c.sha(OUT/'h2-extrapolation.json'),
        H2_small_step_candidate_used=h2_used,
        H2_chain_candidate_used=True,H2_chain_sha256=c.sha(OUT/'h2-chain-tangent.json'),
        preconditioner_scope='Use the separately passing native partial plus auxiliary chain for H2, with the passing original total-loss direction. Preserve failed direct H2 estimates. This initial numerical operator is not a continuous-domain certificate; require independent trajectory refinement.',
        independent_source_derivatives_used=True,physical_or_continuous_certificate=False))


def all_zone_derivative_audit():
    """Expose the original burning-zone mask; never imply a full-domain gate."""
    comp=json.loads((OUT/'composition-derivatives.json').read_text())
    reference=dict(np.load(OUT/'initial-corrected.npz'));J=np.load(OUT/'initial-native.npz')['jacobian']
    pivot=c.NAMES.index('he4');rows=[]
    for row in comp['rows']:
        name=row['species'];refined=(OUT/f'refined-{name}.json').exists()
        sizes=[1e-8,5e-9] if refined else [1e-6,5e-7];prefix='X2' if refined else 'X';values=[]
        for h in sizes: values.append((np.load(OUT/f'{prefix}-{name}-{h}-corrected.npz')['dxdt']-reference['dxdt'])/h)
        derivative=2*values[1]-values[0];reported=J[:,:,c.NAMES.index(name)]-J[:,:,pivot]
        scale=np.maximum(1e-30,np.maximum(abs(derivative).max(1),abs(reported).max(1)))
        difference=abs(values[1]-values[0]).max(1);score=difference/scale;i=int(np.argmax(score))
        rows.append(dict(species=name,all_zone_step_difference_score=float(score[i]),worst_cell=i,
            worst_cell_heat_erg_g_s=float(reference['heat'][i]),worst_cell_absolute_vector_step_difference=float(difference[i]),
            worst_cell_extrapolated_derivative_norm=float(abs(derivative[i]).max()),
            worst_cell_returned_derivative_norm=float(abs(reported[i]).max()),
            passed_same_relative_gate=bool(score[i]<1e-3)))
    save('all-zone-derivative-audit.json',dict(classification='Counterexample candidate',complete=comp['complete'],rows=rows,
        passed=all(r['passed_same_relative_gate'] for r in rows),
        original_gate='The registered derivative gate used abs(initial heat)>1 erg/g/s. It is not an all-zone or continuous-domain certificate.',
        interpretation='Small physical source values do not waive a failed relative derivative gate. Keep absolute values visible; the initial tangent remains an explicitly numerical preconditioner whose trajectories need separate time checks.'))
    print('ALL ZONES',[(r['species'],r['all_zone_step_difference_score']) for r in rows],flush=True)


def refine_h2_extrapolation():
    base=dict(np.load(OUT/'initial-state.npz'));reference=dict(np.load(OUT/'initial-corrected.npz'))
    j=c.NAMES.index('h2');pivot=c.NAMES.index('he4');sizes=[1e-9,5e-10,2.5e-10]
    snapshot=OUT/'h2-original-all-zone-audit.json'
    if not snapshot.exists(): shutil.copy2(OUT/'all-zone-derivative-audit.json',snapshot)
    save('h2-extrapolation-plan.json',dict(classification='Counterexample candidate',
        original_all_zone_audit_snapshot=snapshot.name,original_all_zone_audit_sha256=c.sha(snapshot),steps=sizes,criterion=1e-3,
        observation='In cell 3184 the original H2 one-sided steps differ by about 2.05e-7 while their extrapolated derivative norm is 8.08e-15. The initial heat is 8.74e-15 erg/g/s, so the registered burning-zone mask excluded this cell.',
        method='Compare two second-order extrapolated derivatives from three new positive simplex steps. This tests the actual extrapolated estimator; it is distinct from and does not repair the original first-order step-difference verdict. Require the same 1e-3 vector-relative criterion over all 5735 cells.',
        physical_or_continuous_certificate=False))
    values=[];thermal0=thermal_loss('initial')
    for h in sizes:
        data={k:a.copy() for k,a in base.items()};data['X'][:,j]+=h;data['X'][:,pivot]-=h
        label=f'X3-h2-{h}';_,result=evaluate(label,data)
        row={key:(result[key]-reference[key])/h for key in reference}
        row['total_loss']=row['neutrino']+(thermal_loss(label)-thermal0)/h;values.append(row)
    estimates=[{key:2*values[i+1][key]-values[i][key] for key in values[0]} for i in range(2)]
    a,b=estimates;scale=np.maximum(1e-30,np.maximum(abs(a['dxdt']).max(1),abs(b['dxdt']).max(1)))
    score=abs(a['dxdt']-b['dxdt']).max(1)/scale
    loss_score=abs(a['total_loss']-b['total_loss'])/np.maximum(1,np.maximum(abs(a['total_loss']),abs(b['total_loss'])))
    result=dict(classification='Counterexample candidate',steps=sizes,
        maximum_all_zone_vector_estimator_score=float(score.max()),worst_vector_cell=int(np.argmax(score)),
        maximum_all_zone_loss_estimator_score=float(loss_score.max()),worst_loss_cell=int(np.argmax(loss_score)),
        vector_passed=bool(score.max()<1e-3),loss_passed=bool(loss_score.max()<1e-3),
        original_one_sided_failure_preserved=True,physical_or_continuous_certificate=False)
    np.savez_compressed(OUT/'h2-extrapolation-derivative.npz',**b)
    save('h2-extrapolation.json',result);print('H2 EXTRAPOLATION',result,flush=True)


def h2_chain_tangent():
    base=dict(np.load(OUT/'initial-state.npz'));native=dict(np.load(OUT/'initial-native.npz'))
    _,p=c.mesa(OUT/'initial-profile.data.gz');aux=np.load(OUT/'initial-replacement.npy')
    fixed=np.c_[p['abar'],p['zbar'],p['z2bar'],aux];j=c.NAMES.index('h2');pivot=c.NAMES.index('he4')
    sizes=[1e-6,5e-7];slopes=[(np.load(OUT/f'X-h2-{h}-replacement.npy')-aux)/h for h in sizes]
    da=2*slopes[1]-slopes[0];dS=1/c.A[j]-1/c.A[pivot];dQ=c.Z[j]/c.A[j]-c.Z[pivot]/c.A[pivot]
    dQ2=c.Z[j]**2/c.A[j]-c.Z[pivot]**2/c.A[pivot]
    direction=np.c_[-p['abar']**2*dS,p['abar']*(dQ-p['zbar']*dS),p['abar']*(dQ2-p['z2bar']*dS),da]
    save('h2-chain-plan.json',dict(classification='Counterexample candidate',steps=sizes,criterion=1e-3,
        source_partial_control_sha256=c.sha(OUT/'composition-chain-boundary.json'),
        failed_direct_extrapolation_sha256=c.sha(OUT/'h2-extrapolation.json'),
        method='Keep every nuclear abundance at its actual initial value while perturbing the seven independent auxiliary inputs along D_X a for H2-He4. Add the returned fixed-auxiliary H2-He4 partial. The direct weak-table correction is constant in this direction at fixed auxiliaries because neither H2 nor He4 is a corrected parent.',
        auxiliary_direction='Analytic three composition moments plus the already measured two-step extrapolated seeded-EOS auxiliary direction.',
        comparison='Compare total reconstructed vectors at two auxiliary directional step sizes over all 5735 cells. Separately compare against the registered direct derivative in burning cells. Preserve both failed direct all-zone controls.',
        physical_or_continuous_certificate=False))
    parents=json.loads((OLD/'linear-table-audit.json').read_text())['parent_species']
    assert 'h2' not in parents and 'he4' not in parents
    partial=native['jacobian'][:,:,j]-native['jacobian'][:,:,pivot];estimates=[]
    for h in sizes:
        sides=[]
        for sign in [-1,1]:
            label=f'H2-chain-{h}-{sign}';raw=evaluate_replaced(label,base,fixed+sign*h*direction)
            sides.append(v.corrected_values(raw)['dxdt'])
        estimates.append(partial+(sides[1]-sides[0])/(2*h))
    a,b=estimates;scale=np.maximum(1e-30,np.maximum(abs(a).max(1),abs(b).max(1)))
    score=abs(a-b).max(1)/scale
    direct=np.load(OUT/'X-h2-derivative.npz')['dxdt'];active=abs(native['heat'])>1
    direct_score=abs(b-direct).max(1)/np.maximum(1e-30,np.maximum(abs(b).max(1),abs(direct).max(1)))
    result=dict(classification='Counterexample candidate',steps=sizes,maximum_all_zone_vector_step_score=float(score.max()),
        worst_cell=int(score.argmax()),all_zone_finite_chain_passed=bool(score.max()<1e-3),
        maximum_burning_direct_comparison_score=float(max(direct_score[active])),burning_direct_passed=bool(max(direct_score[active])<1e-3),
        native_partial_plus_auxiliary_chain=True,direct_failures_preserved=True,physical_or_continuous_certificate=False)
    np.savez_compressed(OUT/'h2-chain-tangent.npz',dxdt=b,auxiliary_direction=direction,
        all_zone_step_score=score,direct_comparison_score=direct_score)
    save('h2-chain-tangent.json',result);print('H2 CHAIN',result,flush=True)


def enthalpy_inverse(eos,lp,x,target,guess,budget):
    def sample(t):
        a=eos(1,lp,t,x);return a[2]+a[1]/a[0]-target,a
    t=guess;best=(float('inf'),t,None)
    for _ in range(15):
        defect,a=sample(t)
        if abs(defect)<abs(best[0]): best=(defect,t,a)
        if abs(defect)<budget*.25: break
        trial=t-np.clip(defect/(a[10]-a[1]/a[0]*a[8]),-.1,.1)
        if trial==t: break
        t=trial
    if abs(best[0])>budget*.25:
        lo=hi=best[1]
        for _ in range(3):
            lo=np.nextafter(lo,-np.inf);hi=np.nextafter(hi,np.inf)
            for t in [lo,hi]:
                defect,a=sample(t)
                if abs(defect)<abs(best[0]): best=(defect,t,a)
    return best


def source_step(data,source,thermal,tangent,eosgrad,dt,flux,label):
    from scipy.linalg import expm
    initial=data.get('accumulation_base',data['X']);total=data.get('accumulated_X',np.zeros_like(data['X'],dtype=np.longdouble))
    eos=v.ColdEOS();N=np.exp(data['nu']);qflux=(flux[:-1]-flux[1:])/(data['dm']*N*N)
    rest=(c.W/c.A-1)*(c.gr.C*100)**2;deltaX=[];deltaT=[];loss=[]
    for i in range(len(N)):
        J=tangent['f_X'][i]+np.outer(tangent['f_lnd'][i],eosgrad['rho_X'][i])
        ft=tangent['f_lnT'][i]+tangent['f_lnd'][i]*eosgrad['rho_lnT'][i]
        nx=tangent['loss_X'][i]+tangent['loss_lnd'][i]*eosgrad['rho_X'][i]
        nt=tangent['loss_lnT'][i]+tangent['loss_lnd'][i]*eosgrad['rho_lnT'][i]
        f=source['dxdt'][i];neutrino=source['neutrino'][i]+thermal[i]
        energy=eosgrad['energy'][i];capacity=eosgrad['capacity'][i];ux=eosgrad['energy_X'][i]
        mat=np.zeros((29,29));mat[:26,:26]=J;mat[:26,26]=ft;mat[:26,-1]=f
        mat[26,:26]=(-(rest+ux)@J-nx)/capacity
        mat[26,26]=(-(rest+ux)@ft-nt)/capacity
        mat[26,-1]=(-(rest+ux)@f-neutrino-qflux[i])/capacity
        mat[27,:26]=nx;mat[27,26]=nt;mat[27,-1]=neutrino
        delta=expm(dt*N[i]*mat)[:,-1]
        deltaX.append(delta[:26]);deltaT.append(delta[26]);loss.append(delta[27])
    xx,updated,projection=accumulate(initial,total,np.array(deltaX));loss=np.array(loss)
    assert loss.min()>=-1e-12,('negative integrated neutrino loss',loss.min())
    temperatures=[];densities=[];defects=[];released=[]
    for i in range(len(N)):
        restchange=float(rest.astype(np.longdouble)@(xx[i].astype(np.longdouble)-data['X'][i].astype(np.longdouble)))
        e0=eosgrad['energy'][i];target=e0-restchange-loss[i]-dt*N[i]*qflux[i]
        budget=max(2.,32*np.spacing(abs(target)),abs(target-e0)*1e-8)
        best=enthalpy_inverse(eos,data['logP'][i],xx[i],target,data['lnT'][i]+deltaT[i],budget)
        if abs(best[0])>budget:
            save(label+'-failure.json',dict(classification='Counterexample candidate',cell=i,
                reason='energy inversion',defect_erg_g=best[0],budget_erg_g=budget))
            raise ValueError(('energy inversion',i,best[0],budget))
        temperatures.append(best[1]);densities.append(np.log(best[2][0]));defects.append(best[0]/budget);released.append(abs(restchange))
    result={**data,'X':xx,'accumulation_base':initial,'accumulated_X':updated,'lnT':np.array(temperatures),'lnd':np.array(densities)}
    np.savez_compressed(OUT/f'{label}-source.npz',**result,loss_per_baryon_gram=loss)
    save(label+'-source.json',dict(classification='Counterexample candidate',dt=dt,projection=projection,
        maximum_inverse_score=float(max(abs(np.array(defects)))),
        absolute_rest_release_erg=float(data['dm']@np.array(released)),
        neutrino_energy_infinity_erg=float(data['dm']@(N*loss)),fixed_initial_source_tangent=True))
    return result,loss


def source_control():
    import sympy as s
    q,u,C,fx,nu,flux=s.symbols('q u C fx nu flux',real=True)
    temperature_rate=(-(q+u)*fx-nu-flux)/C
    assert s.simplify(C*temperature_rate+(q+u)*fx+nu+flux)==0
    # Synthetic linear conversion exercises the actual source_step routine,
    # including its exponential, abundance accumulation and energy inversion.
    initial=np.zeros((1,26));ih=c.NAMES.index('h1');ia=c.NAMES.index('he4');initial[0,ih]=.8;initial[0,ia]=.2
    rate=.7;dt=2.;capacity=1e18;lossrate=1e10
    data=dict(X=initial,lnT=np.zeros(1),lnd=np.zeros(1),logP=np.zeros(1),nu=np.zeros(1),dm=np.ones(1))
    f=np.zeros_like(initial);f[0,ih]=-.8*rate;f[0,ia]=.8*rate
    J=np.zeros((1,26,26));J[0,ih,ih]=-rate;J[0,ia,ih]=rate
    tangent=dict(f_X=J,f_lnT=np.zeros_like(f),f_lnd=np.zeros_like(f),
        loss_X=np.zeros_like(f),loss_lnT=np.zeros(1),loss_lnd=np.zeros(1))
    grad=dict(energy=np.ones(1),capacity=np.array([capacity]),rho_lnT=np.zeros(1),energy_X=np.zeros_like(f),rho_X=np.zeros_like(f))
    class LinearControlEOS:
        def __call__(self,mode,lp,t,x):
            a=np.zeros(12);a[0]=a[1]=1.;a[2]=capacity*t;a[10]=capacity;return a
    old=v.ColdEOS
    try:
        v.ColdEOS=LinearControlEOS
        answer,loss=source_step(data,dict(dxdt=f,neutrino=np.array([lossrate])),np.zeros(1),tangent,grad,dt,np.zeros(2),'control-linear')
    finally: v.ColdEOS=old
    exactX=initial.copy();exactX[0,ih]=.8*np.exp(-rate*dt);exactX[0,ia]=1-exactX[0,ih]
    rest=(c.W/c.A-1)*(c.gr.C*100)**2
    theta=(-rest@(exactX[0]-initial[0])-lossrate*dt)/capacity
    ex=float(np.max(abs(answer['X']-exactX)));et=float(abs(answer['lnT'][0]-theta));en=float(abs(loss[0]-lossrate*dt))
    assert ex<1e-14 and et<1e-12 and en<1e-4,(ex,et,en)
    save('source-control.json',dict(classification='Proven',symbolic_energy_identity=True,
        analytic_linear_conversion_errors=dict(composition=ex,lnT=et,integrated_loss=en),passed=True,
        scope='Algebraic identity and the implemented integrator on a specified synthetic linear conversion with a constant heat capacity. No physical EOS or stellar evolution claim.'))
    print('PASS coupled source identity and analytic linear control',ex,et,en,flush=True)


def source_pilot():
    original=dict(np.load(OUT/'initial-state.npz'));snapshot=OUT/'source-pilot-tangent.npz'
    assert not snapshot.exists();shutil.copy2(OUT/'source-tangent.npz',snapshot)
    dt=json.loads((OUT/'plan.json').read_text())['time']['coordinate_seconds']
    save('source-pilot-plan.json',dict(classification='Counterexample candidate',dt=dt,
        tangent_snapshot=snapshot.name,tangent_sha256=c.sha(snapshot),
        scope='First actual 5735-cell source step under the originally registered fixed preconditioner, before structural relaxation. It does not count as a GR time path or a derivative certificate.'))
    source=dict(np.load(OUT/'initial-corrected.npz'));grad=dict(np.load(OUT/'initial-eos-directions.npz'))
    flux=np.load(OLD/'registered-face-flux.npz')['canonical']
    answer,loss=source_step(original,source,thermal_loss('initial'),dict(np.load(snapshot)),grad,dt,flux,'source-pilot')
    print('ACTUAL SOURCE PILOT',json.loads((OUT/'source-pilot-source.json').read_text()),flush=True)


def rest_mass_rounding_audit():
    old=np.load(OUT/'initial-state.npz');new=np.load(OUT/'source-pilot-source.npz')
    dm=old['dm'].astype(np.longdouble);dx=new['X'].astype(np.longdouble)-old['X'].astype(np.longdouble)
    error=dx.sum(1,dtype=np.longdouble);c2=(c.gr.C*100)**2
    direct=(new['X']/c.A)@c.W-(old['X']/c.A)@c.W
    anchored=dx@(c.W.astype(np.longdouble)/c.A-1)
    release=json.loads((OUT/'source-pilot-source.json').read_text())['absolute_rest_release_erg']
    result=dict(classification='Counterexample candidate',
        maximum_rounded_baryon_fraction_change=float(abs(error).max()),
        signed_unit_baryon_rest_energy_change_erg=float(dm@error*c2),
        sum_absolute_unit_baryon_rest_energy_changes_erg=float(dm@abs(error)*c2),
        direct_coefficient_subtraction_minus_anchored_increment_erg=float(dm@(direct-anchored)*c2),
        difference_relative_to_release=float(dm@(direct-anchored)*c2/release),
        inherited_initial_CX_metadata_max_difference=float(abs(old['CX']-(old['X']/c.A)@c.W).max()),
        treatment='Source energy and the established stable mass-change diagnostic use rest-mass excess increments, avoiding subtraction of the unit baryon rest energy. The GR structure still evaluates its declared neutral-atom mass coefficient directly. New structural state CX metadata is recomputed from its actual X; intermediate source states retain pre-relaxation geometry.',
        limitation='This finite rounding audit is not a rigorous ADM mass-difference error enclosure. The unsigned contribution cannot be discarded on the strength of signed cancellation. Prior states and their metadata are preserved.')
    save('rest-mass-rounding-audit.json',result);print('REST MASS ROUNDING',result,flush=True)


def energy_quadrature_audit():
    """Local shell quadrature anchored to the two saved GR midpoint states."""
    from numpy.polynomial.legendre import leggauss
    old=dict(np.load(OUT/'initial-state.npz'));new=dict(np.load(OUT/'P-1-0-state-17-4.npz'))
    source=dict(np.load(OUT/'P-1-0-source.npz'));run=json.loads((OUT/'evolution.json').read_text())['rows'][0]
    save('energy-quadrature-plan.json',dict(classification='Counterexample candidate',orders=[1,2,4],
        old_state_sha256=c.sha(OUT/'initial-state.npz'),new_state_sha256=c.sha(OUT/'P-1-0-state-17-4.npz'),
        original_one_step_energy_score=run['energy_score'],
        target='Separate the saved midpoint mass-functional quadrature from the physical source-step error. Anchor every local integration to the saved r,m,P midpoint and preserve each shell composition and entropy.',
        method='Use the existing local baryon-coordinate RK4 EOS-table path to Gauss nodes, then query seeded FreeEOS at the saved shell entropy. Compare one, two and four nodes without changing the original trajectory or its failed verdict.',
        checks=dict(midpoint_reconstruction_relative_to_release=1e-7,quadrature_refinement_relative_to_release=1e-7),
        scope='A finite local-quadrature diagnostic. It is not a new globally solved GR trajectory, a rigorous interpolation enclosure or a physical EOS certificate.'))
    settings=[(v.OLD,'restored',dict(np.load(v.OLD/'restored-input.npz')),old),(OUT,'P-1-0',source,new)]
    solvers=[];saved_out,saved_eos=c.OUT,c.EOS
    try:
        c.EOS=v.ColdEOS
        for folder,label,data,state in settings:
            c.OUT=folder;solver=c.Structure(label,data,17,4)
            # The diagnostic follows the actual saved shell entropy. This
            # also makes the one-node comparison a direct reconstruction.
            solver.ref=solver.ref.copy();solver.ref[:,3]=state['s_B']
            solver.mat.eps=state['X'];solver.mat.lt=state['lnT'];solvers.append(solver)
        c.OUT=saved_out
        def nodes(solver,state,z):
            m=solver.mat;rows=[]
            for i in range(len(m.dm)):
                outside=i<m.split
                start=(m.outer[i]+m.outer[i+1])/2 if outside else (m.inner[i]+m.inner[i+1])/2
                target=m.outer[i]+z*m.dm[i] if outside else m.inner[i+1]+(1-z)*m.dm[i]
                y=np.array([state['r_mid_m'][i]/m.R,state['m_mid_geom'][i]/m.B,state['logP'][i]])
                if z!=.5: y=solver.step(np.log(start),np.log(target),y,i,m.B,outside)
                a,t,_=c.be.invert(solver.eos,y[2],state['s_B'][i],state['X'][i],state['lnT'][i])
                rows.append([y[0]*m.R,y[1]*m.B,a[2],a[3]-state['s_B'][i]])
            return np.array(rows)
        dcx=((new['X']-old['X'])/c.A)@(c.W-c.A);cx=(old['X']/c.A)@c.W;c2=(c.gr.C*100)**2
        result=[];arrays={}
        for order in [1,2,4]:
            xx,ww=leggauss(order);total=0.;entropy=[]
            for j,(z,w) in enumerate(zip((xx+1)/2,ww/2)):
                a=nodes(solvers[0],old,float(z));b=nodes(solvers[1],new,float(z))
                f0=np.sqrt(1-2*a[:,1]/a[:,0]);f1=np.sqrt(1-2*b[:,1]/b[:,0])
                df=-2*((b[:,1]-a[:,1])-a[:,1]/a[:,0]*(b[:,0]-a[:,0]))/b[:,0]/(f0+f1)
                value=(dcx*c2+b[:,2]-a[:,2])*f1+(cx*c2+a[:,2])*df
                total+=float(w*(old['dm']@value));arrays[f'order{order}_node{j}']=value
                entropy.append(float(max(abs(a[:,3]).max(),abs(b[:,3]).max())))
                print('MASS QUADRATURE',order,j,'node',float(z),flush=True)
            row=dict(order=order,mass_change_energy_erg=total,
                energy_score=abs(total-run['expected_energy_erg'])/run['release_erg'],
                maximum_entropy_residual=max(entropy))
            if order==1:
                row['saved_midpoint_reconstruction_score']=abs(total-run['mass_change_energy_erg'])/run['release_erg']
                row['midpoint_reconstruction_passed']=row['saved_midpoint_reconstruction_score']<1e-7
            else:
                row['previous_order_difference_score']=abs(total-result[-1]['mass_change_energy_erg'])/run['release_erg']
            result.append(row)
            save('energy-quadrature-audit.json',dict(classification='Counterexample candidate',rows=result,
                complete=order==4,original_one_step_failure_preserved=True,
                globally_resolved_GR_path=False,rigorous_EOS_or_ADM_certificate=False))
            np.savez_compressed(OUT/'energy-quadrature-audit.npz',**arrays)
            print('MASS QUADRATURE RESULT',row,flush=True)
    finally: c.OUT,c.EOS=saved_out,saved_eos


def energy_error_budget():
    """Independent, finite attribution checks; no correction of raw verdicts."""
    from scipy.interpolate import PchipInterpolator
    old=dict(np.load(OUT/'initial-state.npz'));N=np.exp(old['nu']);dm=old['dm'];c2=(c.gr.C*100)**2
    flux=np.load(OLD/'registered-face-flux.npz')['canonical'];rows=[]
    for count in [1,2,4]:
        label=f'P-{count}-0';path=OUT/f'{label}-state-17-4.npz'
        if not path.exists(): continue
        end=dict(np.load(path));src=dict(np.load(OUT/f'{label}-source.npz'))
        source=json.loads((OUT/f'{label}-source.json').read_text());tab=np.load(OUT/f'{label}-adiabats-17.npz')
        r0=old['radius_faces_m'][0]*100;r1=end['radius_faces_m'][0]*100
        work=np.exp(float(old['boundary_logP']))*4*np.pi*(r1-r0)*(r1*r1+r1*r0+r0*r0)/3
        measured=c.stable_mass_change(old,end)*c2
        expected=-flux[0]*source['dt']-source['neutrino_energy_infinity_erg']-work
        ds=end['s_B']-tab['reference'][:,3];ds_energy=end['dm']*np.exp(end['nu']+end['lnT'])*ds
        dp=np.exp(old['logP'])*np.expm1(end['logP']-old['logP'])
        dv=np.exp(-old['lnd'])*np.expm1(old['lnd']-src['lnd'])
        dh=-src['loss_per_baryon_gram']-source['dt']*(flux[:-1]-flux[1:])/(dm*N)
        rows.append(dict(dt=source['dt'],source_steps=count,defect_erg=measured-expected,
            release_erg=source['absolute_rest_release_erg'],
            energy_score=abs(measured-expected)/source['absolute_rest_release_erg'],
            entropy_inverse_signed_energy_estimate_erg=float(ds_energy.sum()),
            entropy_inverse_unsigned_energy_estimate_erg=float(abs(ds_energy).sum()),
            leading_pressure_cross_estimate_erg=float(.5*dm@(N*dp*dv)),
            leading_lapse_cross_estimate_erg=float(.5*dm@(N*np.expm1(end['nu']-old['nu'])*dh))))
    grad=np.load(OUT/'initial-eos-directions.npz')
    dh=grad['energy']-old['u_W']-np.exp(old['logP']-old['lnd'])
    interp=[]
    for label,state,path in [('old',old,v.OLD/'restored-adiabats-17.npz'),
            ('new',dict(np.load(OUT/'P-1-0-state-17-4.npz')),OUT/'P-1-0-adiabats-17.npz')]:
        table=np.load(path);offset=state['logP']-np.log(table['reference'][:,1]);grid=table['offset']
        coef=PchipInterpolator(grid,table['values'],axis=1).c;y=[];dy=[]
        assert min(offset)>=grid[0] and max(offset)<=grid[-1]
        for i,point in enumerate(offset):
            j=int(np.clip(np.searchsorted(grid,point)-1,0,len(grid)-2));z=point-grid[j];a=coef[:,j,i]
            y.append(((a[0]*z+a[1])*z+a[2])*z+a[3]);dy.append((3*a[0]*z+2*a[1])*z+a[2])
        y=np.array(y);dy=np.array(dy);u=y[:,2]*c2;rhs=np.exp(state['logP']-y[:,0])*dy[:,0]
        weight=state['dm']*np.sqrt(1-2*state['m_mid_geom']/state['r_mid_m'])
        interp.append(dict(state=label,weighted_direct_minus_interpolated_internal_energy_erg=float(weight@(state['u_W']-u)),
            maximum_adiabatic_first_law_relative_residual=float(max(abs(dy[:,2]*c2-rhs)/np.maximum(abs(rhs),1)))))
    save('energy-error-budget.json',dict(classification='Counterexample candidate',first_step_duration_controls=rows,
        full_to_half_step_defect_ratio=rows[0]['defect_erg']/rows[1]['defect_erg'] if len(rows)>1 else None,
        cold_initial_enthalpy_signed_difference_erg=float((dm*N)@dh),
        cold_initial_enthalpy_unsigned_difference_erg=float((dm*N)@abs(dh)),
        EOS_interpolation_diagnostics=interp,
        conclusion='The finite entropy-inverse, initial-EOS, local pressure/lapse and midpoint-quadrature checks do not individually explain the leading one-step defect. The first shorter-duration control tests time-order scaling; it is not a same-final-time convergence test or a rigorous global error bound.',
        original_failed_verdicts_preserved=True,physical_or_continuous_certificate=False))
    print('ENERGY ERROR BUDGET',rows,flush=True)


def stored_mass_interval_audit():
    """Enclose the declared mass functional on exact saved binary64 inputs."""
    from mpmath import iv
    from fractions import Fraction
    iv.dps=40
    def endpoint(x):
        sign,man,exponent,_=x._mpi_[0]
        return Fraction((-1 if sign else 1)*man)*Fraction(2)**exponent
    exact=lambda x:iv.mpf(float(x))
    initial=dict(np.load(OUT/'initial-state.npz'));runs=json.loads((OUT/'evolution.json').read_text())['rows']
    save('stored-mass-interval-plan.json',dict(classification='Proven',decimal_precision=40,
        input_interpretation='Every saved binary64 value and isotope weight is an exact point input to the declared excess-rest-mass functional.',
        target='Arithmetic evaluation error only. No enclosure of the true GR state, EOS interpolation, entropy inversion, source time discretization or physical model error.',
        diagnostic_criterion_relative_to_release=1e-8,
        interval_source_binding_sha256=c.sha(OLD/'provenance.json')))
    weights=[exact(w) for w in c.W];baryons=[exact(a) for a in c.A];c2=(exact(c.gr.C)*100)**2
    def cell(i,new):
        r0,r1=exact(initial['r_mid_m'][i]),exact(new['r_mid_m'][i])
        m0,m1=exact(initial['m_mid_geom'][i]),exact(new['m_mid_geom'][i])
        f0=iv.sqrt(1-2*m0/r0);f1=iv.sqrt(1-2*m1/r1)
        df=-2*((m1-m0)-m0/r0*(r1-r0))/r1/(f0+f1)
        x0=[exact(x) for x in initial['X'][i]];x1=[exact(x) for x in new['X'][i]]
        cx=sum(x*w/a for x,w,a in zip(x0,weights,baryons))
        dcx=sum((b-a)*(w-n)/n for a,b,w,n in zip(x0,x1,weights,baryons))
        u0,u1=exact(initial['u_W'][i]),exact(new['u_W'][i])
        return exact(initial['dm'][i])*((dcx*c2+u1-u0)*f1+(cx*c2+u0)*df)
    # Exact identity through the implemented rationalized expression.
    control=cell(2501,initial);assert endpoint(control.a)<=0<=endpoint(control.b)
    rows=[]
    for run in runs:
        path=OUT/f"P-{run['steps']}-{run['steps']-1}-state-17-4.npz";state=dict(np.load(path))
        total=iv.mpf(0)
        for i in range(len(initial['dm'])): total+=cell(i,state)
        lo,hi=endpoint(total.a),endpoint(total.b);stored=Fraction(run['mass_change_energy_erg'])
        bound=max(abs(lo-stored),abs(hi-stored));score=bound/Fraction(run['release_erg'])
        rows.append(dict(steps=run['steps'],mass_energy_lower_exact=str(lo),mass_energy_upper_exact=str(hi),
            mass_energy_midpoint_erg=float((lo+hi)/2),arithmetic_error_upper_exact=str(bound),
            arithmetic_error_score_upper_exact=str(score),arithmetic_error_score_upper=float(score),
            passed=score<Fraction(1,10**8),state_sha256=c.sha(path)))
    save('stored-mass-interval-audit.json',dict(classification='Proven',rows=rows,
        initial_state_sha256=c.sha(OUT/'initial-state.npz'),zero_change_interval_contains_zero=True,
        interpretation='Rational endpoints enclose the mathematical excess-rest-mass functional of the exact stored point data. A small arithmetic bound cannot validate the input trajectories or reinterpret their energy failures.',
        full_GR_or_physical_certificate=False))
    print('STORED MASS INTERVAL',[(r['steps'],r['arithmetic_error_score_upper'],r['passed']) for r in rows],flush=True)


def structure_arithmetic_audit():
    """Same saved shooting parameters, EOS table and inputs; vary arithmetic."""
    from types import SimpleNamespace
    label='P-1-0';data=dict(np.load(OUT/f'{label}-source.npz'))
    data['CX']=(data['X']/c.A)@c.W
    native=dict(np.load(OUT/f'{label}-state-17-4.npz'));initial=dict(np.load(OUT/'initial-state.npz'))
    parameters=json.loads((OUT/f'{label}-structure-17-4.json').read_text())['parameters']
    run=json.loads((OUT/'evolution.json').read_text())['rows'][0]
    save('structure-arithmetic-plan.json',dict(classification='Counterexample candidate',
        input_sha256=c.sha(OUT/f'{label}-source.npz'),EOS_table_sha256=c.sha(OUT/f'{label}-adiabats-17.npz'),
        shooting_parameters=parameters,parameters_refitted=False,arithmetics=['binary64','longdouble'],
        scope='Replay the same declared fixed shooting point with the same EOS table. Compare material-coordinate integration arithmetic and direct EOS evaluations. Extended precision does not improve the native binary64 EOS or certify physical/model error.',
        baseline='Require the binary64 replay to reproduce saved geometry and thermodynamics bitwise before interpreting the arithmetic difference.'))
    original_out,original_eos,original_solver,original_fit,original_save=c.OUT,c.EOS,c.Structure,c.least_squares,c.save
    class ExtendedStructure(original_solver):
        def __init__(self,*args,**kwargs):
            super().__init__(*args,**kwargs)
            self.mat.R=np.longdouble(self.mat.R);self.mat.B=np.longdouble(self.mat.B)
        def step(self,a,b,y,i,B,outer):
            return super().step(a,b,np.asarray(y,dtype=np.longdouble),i,B,outer)
    rows=[];states=[]
    try:
        c.OUT=OUT;c.EOS=v.ColdEOS
        for tag,dtype,solver in [('64',float,original_solver),('extended',np.longdouble,ExtendedStructure)]:
            name='arithmetic-'+tag
            path=OUT/f'{name}-adiabats-17.npz'
            if not path.exists(): shutil.copy2(OUT/f'{label}-adiabats-17.npz',path)
            c.Structure=solver
            c.least_squares=lambda fun,x0,**kwargs:SimpleNamespace(success=True,x=np.asarray(parameters,dtype=dtype))
            c.save=lambda name,value:(OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2,
                default=lambda number:float(number))+'\n')
            state,record=c.structure(name,data,17,4,parameters);states.append(state)
            if tag=='64':
                matched={key:bool(np.array_equal(state[key],native[key])) for key in ['lnd','lnT','r_mid_m','m_mid_geom','logP','nu','u_W','s_B']}
                save('structure-arithmetic-baseline.json',dict(classification='Counterexample candidate',bitwise_equal=matched,passed=all(matched.values())))
                assert all(matched.values()),matched
            energy=c.stable_mass_change(initial,state)*(c.gr.C*100)**2
            rows.append(dict(arithmetic=tag,interface_max=float(record['interface_max']),
                mass_change_energy_erg=energy,score_against_original_expected_energy=abs(energy-run['expected_energy_erg'])/run['release_erg'],
                midpoint_radius_dtype=str(state['r_mid_m'].dtype),shooting_parameters_refitted=False))
            print('STRUCTURE ARITHMETIC',rows[-1],flush=True)
    finally: c.OUT,c.EOS,c.Structure,c.least_squares,c.save=original_out,original_eos,original_solver,original_fit,original_save
    shifts={key:float(abs(states[1][key]-states[0][key]).max()) for key in ['r_mid_m','m_mid_geom','logP','lnT','u_W']}
    save('structure-arithmetic-audit.json',dict(classification='Counterexample candidate',rows=rows,
        baseline_bitwise_equal=True,maximum_state_shifts=shifts,
        mass_energy_shift_erg=rows[1]['mass_change_energy_erg']-rows[0]['mass_change_energy_erg'],
        shift_relative_to_release=(rows[1]['mass_change_energy_erg']-rows[0]['mass_change_energy_erg'])/run['release_erg'],
        original_failed_verdict_preserved=True,new_fitted_GR_path=False,continuous_or_physical_certificate=False))


def zero_source_projection():
    """Measure the constraint solver's drift without any source or elapsed time."""
    old=dict(np.load(OUT/'initial-state.npz'))
    release=json.loads((OUT/'evolution.json').read_text())['rows'][0]['release_erg']
    save('zero-source-projection-plan.json',dict(classification='Counterexample candidate',
        initial_state_sha256=c.sha(OUT/'initial-state.npz'),elapsed_time_seconds=0,
        composition_and_baryon_inventories='Exactly unchanged; no source integration or native reaction call.',
        expected_mass_energy_change_erg=0,criterion_erg=release*1e-6,
        criterion_reference='The original full-duration one-step rest release; a zero source has no release denominator of its own.',
        target='Reinsert the saved midpoint state into the same seeded-EOS GR projection. Test numerical idempotence, including entropy inversion and rebasing of the local EOS table.',
        original_trajectory_modified=False,physical_or_continuous_certificate=False))
    parameters=json.loads((v.OLD/'restored-structure-17-4.json').read_text())['parameters']
    state,record=structure('zero-source',old,parameters)
    assert np.array_equal(state['X'],old['X']) and np.array_equal(state['dm'],old['dm'])
    c2=(c.gr.C*100)**2;r0=old['radius_faces_m'][0]*100;r1=state['radius_faces_m'][0]*100
    work=np.exp(float(old['boundary_logP']))*4*np.pi*(r1-r0)*(r1*r1+r1*r0+r0*r0)/3
    energy=c.stable_mass_change(old,state)*c2
    result=dict(classification='Counterexample candidate',structure=record,
        mass_change_energy_erg=energy,pressure_work_erg=work,defect_erg=energy+work,
        energy_score=abs(energy+work)/release,passed=bool(abs(energy+work)<release*1e-6),
        maximum_lnT_drift=float(abs(state['lnT']-old['lnT']).max()),
        maximum_lnd_drift=float(abs(state['lnd']-old['lnd']).max()),
        maximum_entropy_drift=float(abs(state['s_B']-old['s_B']).max()),
        composition_and_dm_bitwise_unchanged=True,original_failed_paths_preserved=True,
        scope='Finite zero-source projection control, not a time path or a rigorous error enclosure.')
    save('zero-source-projection.json',result);print('ZERO SOURCE PROJECTION',result,flush=True)


def time_refinement_audit():
    """Locate endpoint failures without replacing the registered float64 gate."""
    initial=dict(np.load(OUT/'initial-state.npz'))
    runs=json.loads((OUT/'evolution.json').read_text())['rows'];rows=[]
    for coarse,fine in zip(runs,runs[1:]):
        paths=[OUT/f"P-{r['steps']}-{r['steps']-1}-state-17-4.npz" for r in [coarse,fine]]
        a,b=[dict(np.load(p)) for p in paths]
        denom=1e-16+1e-3*abs(b['X']-initial['X']);score=abs(b['X']-a['X'])/denom
        assert float(score.max())==fine['refinement']['composition_score']
        indices=np.argsort(score.ravel())[-8:][::-1];worst=[]
        for flat in indices:
            i,j=np.unravel_index(flat,score.shape)
            worst.append(dict(cell=int(i),species=c.NAMES[j],score=float(score[i,j]),
                coarse_change=float(a['X'][i,j]-initial['X'][i,j]),
                fine_change=float(b['X'][i,j]-initial['X'][i,j]),
                endpoint_difference=float(b['X'][i,j]-a['X'][i,j]),budget=float(denom[i,j]),
                accumulated_endpoint_difference=float(b['accumulated_X'][i,j]-a['accumulated_X'][i,j]),
                initial_temperature_K=float(np.exp(initial['lnT'][i])),
                initial_density_g_cm3=float(np.exp(initial['lnd'][i]))))
        compensated=abs(b['accumulated_X']-a['accumulated_X'])/(1e-16+1e-3*abs(b['accumulated_X']))
        rows.append(dict(coarse_steps=coarse['steps'],fine_steps=fine['steps'],
            registered_composition_score=float(score.max()),
            cells_failing=int(np.any(score>=1,axis=1).sum()),
            species_failing=[c.NAMES[j] for j in np.where(np.any(score>=1,axis=0))[0]],
            compensated_composition_diagnostic=float(compensated.max()),worst=worst,
            state_sha256={p.name:c.sha(p) for p in paths}))
    save('time-refinement-audit.json',dict(classification='Counterexample candidate',rows=rows,
        scope='Saved same-final-time endpoint localization; compensated increments diagnose input rounding only. The original gate and every failed trajectory remain unchanged.',
        continuous_time_error_certificate=False))
    print('TIME REFINEMENT',rows,flush=True)


def structural_source_diagnostic():
    """Measure compression contributions at the worst saved Be7 time control."""
    from scipy.linalg import expm
    base=dict(np.load(OUT/'initial-state.npz'));src=dict(np.load(OUT/'P-2-0-source.npz'))
    relaxed=dict(np.load(OUT/'P-2-0-state-17-4.npz'))
    before=np.load(OUT/'initial-corrected.npz')['dxdt'];after=np.load(OUT/'P-2-1-corrected.npz')['dxdt']
    tangent=np.load(OUT/'source-tangent.npz');grad=np.load(OUT/'initial-eos-directions.npz')
    audit=json.loads((OUT/'time-refinement-audit.json').read_text());worst=audit['rows'][0]['worst'];rows=[]
    dt=json.loads((OUT/'plan.json').read_text())['time']['coordinate_seconds']/2
    for point in worst[:3]:
        i=point['cell'];j=c.NAMES.index(point['species'])
        dx=relaxed['X'][i]-base['X'][i];dr=relaxed['lnd'][i]-base['lnd'][i];dT=relaxed['lnT'][i]-base['lnT'][i]
        terms=[tangent['f_X'][i]@dx,tangent['f_lnd'][i]*dr,tangent['f_lnT'][i]*dT]
        compression=tangent['f_lnd'][i]*(relaxed['lnd'][i]-src['lnd'][i])+tangent['f_lnT'][i]*(relaxed['lnT'][i]-src['lnT'][i])
        J=tangent['f_X'][i]+np.outer(tangent['f_lnd'][i],grad['rho_X'][i])
        matrix=np.zeros((27,27));matrix[:26,:26]=J;matrix[:26,-1]=compression
        response=expm(dt*np.exp(base['nu'][i])*matrix)[:26,-1]
        rows.append(dict(cell=i,species=c.NAMES[j],actual_source_change=float(after[i,j]-before[i,j]),
            initial_source=float(before[i,j]),after_half_step_source=float(after[i,j]),
            linear_composition_contribution=float(terms[0][j]),linear_density_contribution=float(terms[1][j]),
            linear_temperature_contribution=float(terms[2][j]),
            linearization_residual=float(after[i,j]-before[i,j]-sum(terms)[j]),
            structural_lnT_change=float(relaxed['lnT'][i]-src['lnT'][i]),
            structural_lnd_change=float(relaxed['lnd'][i]-src['lnd'][i]),
            compression_source_contribution=float(compression[j]),
            frozen_temperature_compression_response=float(response[j]),
            measured_full_vs_two_step_endpoint_difference=point['endpoint_difference']))
    save('structural-source-diagnostic.json',dict(classification='Counterexample candidate',rows=rows,
        scope='Local initial-tangent decomposition at saved endpoints. The 26-species response holds temperature fixed and is an attribution diagnostic, not an improved time path or a global derivative certificate.',
        original_failed_verdicts_preserved=True))
    print('STRUCTURAL SOURCE',rows,flush=True)


def structural_split_control():
    """Classical local order boundary for a missing algebraic-structure chain."""
    import sympy as sy
    from scipy.linalg import expm
    h,w,k,y=sy.symbols('h w k y',real=True)
    # F(y,G(y))=(w+k)y; the local source operator includes w alone.
    delta=sum(h**j*w**(j-1)/sy.factorial(j) for j in range(1,4))*(w+k)*y
    exact=sy.series(y*sy.exp(h*(w+k)),h,0,4).removeO()
    defect=sy.expand(exact-y-delta)
    assert sy.simplify(defect.coeff(h,2)-k*(w+k)*y/2)==0
    remainder=k*delta
    corrected=y+delta+(h/2+h*h*w/6)*remainder
    assert sy.expand(exact-corrected).coeff(h,2)==0
    # Independent nonlinear scalar control with an intentionally partial W.
    rate=.7;beta=.2;duration=2.;initial=.8;operator=-rate
    exact_value=initial*np.exp(-rate*duration)/(1-beta*initial/rate*(1-np.exp(-rate*duration)))
    f=lambda value:-rate*value+beta*value*value
    rows=[]
    for count in [16,32,64,128]:
        step=duration/count;first=second=initial
        for _ in range(count):
            first+=expm(step*np.array([[operator,f(first)],[0.,0.]]))[0,1]
            d=expm(step*np.array([[operator,f(second)],[0.,0.]]))[0,1]
            residual=f(second+d)-f(second)-operator*d
            aug=np.zeros((3,3));aug[0,0]=operator;aug[0,1]=residual;aug[1,2]=1/step
            second+=d+expm(step*aug)[0,2]
        rows.append(dict(steps=count,partial_operator_error=abs(first-exact_value),
            residual_corrected_error=abs(second-exact_value)))
    ratios=[rows[j]['residual_corrected_error']/rows[j+1]['residual_corrected_error'] for j in range(len(rows)-1)]
    assert min(ratios)>3.8 and max(ratios)<4.2,ratios
    save('structural-split-control.json',dict(classification='Proven',
        missing_second_order_term='h^2/2 * (F_g G_y) F for a source-local W=F_y, evaluated at the starting state',
        linear_symbolic_control=True,nonlinear_control=rows,residual_corrected_error_ratios=ratios,passed=True,
        candidate_formula='d=h phi1(hW) F(y); y_next=y+d+h phi2(hW)[F(y+d)-F(y)-Wd]',
        scope='Classical local Taylor order for smooth finite-dimensional reduced functions and the stated controls. The complete F must include the algebraic structure response in a consistent independent-variable chart. No uniform stiff bound, EOS/GR accuracy, positivity or actual corrected stellar path is established.',
        stellar_candidate_implemented=False,
        related_method_source='https://doi.org/10.1137/080717717',
        source_boundary='The cited exponential Rosenbrock work uses flow linearization and its own order assumptions; its higher-order theorems are not transferred to the inexact stellar operator.'))
    print('PASS structural split boundary and nonlinear residual control',rows,ratios,flush=True)


def enthalpy_chart_control():
    """Reference-pressure enthalpy chart for a thermodynamically consistent EOS."""
    import sympy as sy
    T,Tr,mu,mur,r,f,q=sy.symbols('T Tr mu mur r f q',positive=True)
    Sdot=(-q-(r+mu)*f)/T
    Hdot=Tr*Sdot+mur*f
    expected=Tr/T*(-q-r*f)+(mur-Tr/T*mu)*f
    assert sy.simplify(Hdot-expected)==0
    gauge=sy.symbols('gauge',real=True)
    assert sy.simplify((mur-Tr*gauge)-Tr/T*(mu-T*gauge)-(mur-Tr/T*mu))==0
    x,cp,gas,P,Pr,u=sy.symbols('x cp gas P Pr u',positive=True)
    # A composition entropy singularity cancels in the reference difference.
    entropy=cp*sy.log(T)-gas*sy.log(P)+x*sy.log(x)
    enthalpy=cp*T+u*x
    chemical=sy.diff(enthalpy,x)-T*sy.diff(entropy,x)
    ratio=(Pr/P)**(gas/cp)
    difference=sy.simplify(chemical.subs(T,T*ratio)-ratio*chemical)
    assert sy.simplify(difference-(1-ratio)*u)==0
    assert sy.simplify(sy.limit(difference,x,0,dir='+')-(1-ratio)*u)==0
    save('enthalpy-chart-control.json',dict(classification='Proven',passed=True,
        independent_coordinates='X and H_ref(X,s), where the reference baryon-shell pressure is fixed throughout a path.',
        local_proper_time_equation='dH_ref/dtau=(T_ref/T)(-r dot f-epsilon_nu-q_flux)+(mu_ref-(T_ref/T)mu) dot f',
        chemical_difference='H_X_ref-(T_ref/T)H_X+T_ref(s_X-s_X_ref), with each derivative at its own fixed pressure and temperature.',
        entropy_composition_gauge_invariant=True,
        ideal_gas_composition_log_singularity_cancels=True,
        limits='Requires the same differentiable thermodynamic EOS, positive heat capacity, a common composition/entropy branch and a differentiable constraint map. Cancellation of an additive composition-only entropy reference does not prove smoothness across physical EOS phase or ionization boundaries. No stellar residual-corrected path is implemented by this identity.'))
    print('PASS reference-pressure enthalpy chart and entropy-reference cancellation',flush=True)


def enthalpy_chart_eos_audit(full=False):
    assert json.loads((OUT/'enthalpy-chart-control.json').read_text())['passed']
    prefix='P-1-0' if full else 'P-2-0';native='residual-predictor' if full else 'P-2-1'
    stem='residual-chart' if full else 'enthalpy-chart-eos'
    ref=dict(np.load(OUT/f'{prefix}-source.npz'));actual=dict(np.load(OUT/f'{prefix}-state-17-4.npz'))
    source=dict(np.load(OUT/f'{native}-corrected.npz'));thermal=thermal_loss(native)
    flux=np.load(OLD/'registered-face-flux.npz')['canonical'];eos=v.ColdEOS()
    selected=list(range(len(actual['X']))) if full else [0,2501,len(actual['X'])-1]
    sizes=[1e-5,5e-6];pivot=c.NAMES.index('he4');capacity=np.load(OUT/'initial-eos-directions.npz')['capacity']
    save(stem+'-plan.json',dict(classification='Counterexample candidate',cells=selected,
        composition_steps=sizes,pivot='he4',criterion_relative_to_chart_source=1e-3,
        input_sha256={name:c.sha(OUT/name) for name in [f'{prefix}-source.npz',f'{prefix}-state-17-4.npz',f'{native}-corrected.npz']},
        target='Check the chemical-potential difference entering the reference-pressure enthalpy chart using two finite simplex sizes at the listed saved physical/reference pairs. This is not a continuous-domain derivative certificate.'))
    rows=[];vectors=[];rest=(c.W/c.A-1)*(c.gr.C*100)**2
    for i in selected:
        x=actual['X'][i];assert np.array_equal(x,ref['X'][i])
        tr=np.exp(ref['lnT'][i]);ta=np.exp(actual['lnT'][i]);ratio=tr/ta
        ar=eos(1,ref['logP'][i],ref['lnT'][i],x);aa=eos(1,actual['logP'][i],actual['lnT'][i],x)
        hr=ar[2]+ar[1]/ar[0];ha=aa[2]+aa[1]/aa[0];directions=[]
        for h in sizes:
            derivative=np.zeros(26)
            for j in range(26):
                if j==pivot: continue
                z=x.copy();z[j]+=h;z[pivot]-=h
                br=eos(1,ref['logP'][i],ref['lnT'][i],z);ba=eos(1,actual['logP'][i],actual['lnT'][i],z)
                derivative[j]=((br[2]+br[1]/br[0]-hr)-ratio*(ba[2]+ba[1]/ba[0]-ha)
                    +tr*((ba[3]-aa[3])-(br[3]-ar[3])))/h
            directions.append(derivative)
        N=np.exp(actual['nu'][i]);qflux=(flux[i]-flux[i+1])/(actual['dm'][i]*N*N)
        f=source['dxdt'][i];loss=source['neutrino'][i]+thermal[i]
        reference_part=ratio*(-rest@f-loss-qflux)
        chemical=[float(d@f) for d in directions];chart=[N*(reference_part+term) for term in chemical]
        score=abs(chart[1]-chart[0])/max(1,abs(chart[0]),abs(chart[1]))
        vectors.append(np.r_[N*f,(2*chart[1]-chart[0])/capacity[i],N*N*loss/capacity[i]])
        rows.append(dict(cell=i,T_ref_over_T=float(ratio),entropy_pair_residual=float(aa[3]-ar[3]),
            chemical_source_corrections_erg_g_s=chemical,coordinate_chart_source_erg_g_s=chart,
            finite_step_source_score=float(score),passed=bool(score<1e-3)))
        if full and i%1000==0: print('RESIDUAL CHART EOS',i,flush=True)
    np.savez_compressed(OUT/(stem+'-source.npz'),F=np.array(vectors),cells=np.array(selected),capacity=capacity[selected])
    save(stem+'-audit.json',dict(classification='Counterexample candidate',rows=rows,complete=True,
        passed=all(row['passed'] for row in rows),all_zone_or_continuous_certificate=False,
        corrected_stellar_path_completed=False))
    print('ENTHALPY CHART EOS',len(rows),'max source score',max(row['finite_step_source_score'] for row in rows),
        'passed',all(row['passed'] for row in rows),flush=True)


def residual_chart_source(): enthalpy_chart_eos_audit(full=True)


def enthalpy_coordinate_predictor():
    """Re-express the initial exponential predictor in dimensionless chart units."""
    from scipy.linalg import expm
    data=dict(np.load(OUT/'initial-state.npz'));old=dict(np.load(OUT/'P-1-0-source.npz'))
    tangent=dict(np.load(OUT/'source-tangent.npz'));grad=dict(np.load(OUT/'initial-eos-directions.npz'))
    source=dict(np.load(OUT/'initial-corrected.npz'));thermal=thermal_loss('initial')
    flux=np.load(OLD/'registered-face-flux.npz')['canonical'];N=np.exp(data['nu'])
    dt=json.loads((OUT/'plan.json').read_text())['time']['coordinate_seconds'];capacity=grad['capacity']
    release=json.loads((OUT/'P-1-0-source.json').read_text())['absolute_rest_release_erg']
    save('enthalpy-coordinate-plan.json',dict(classification='Counterexample candidate',dt=dt,
        coordinates='26 composition increments, reference-pressure enthalpy increment/C0, neutrino energy at infinity/C0; C0 is the initial positive dH/dlnT at fixed P,X.',
        initial_tangent_sha256=c.sha(OUT/'source-tangent.npz'),
        original_source_sha256=c.sha(OUT/'P-1-0-source.npz'),
        gates=dict(composition_absolute=1e-16,composition_relative_to_total_change=1e-3,
            unsigned_thermal_energy_difference_relative_to_release=1e-6,
            unsigned_neutrino_energy_difference_relative_to_release=1e-6),
        scope='Finite covariance/replay control before reusing a saved GR predictor. A successful coordinate transformation does not turn the failed original GR path into a pass.'))
    rest=(c.W/c.A-1)*(c.gr.C*100)**2;qflux=(flux[:-1]-flux[1:])/(data['dm']*N*N)
    matrices=[];vectors=[];increments=[]
    for i in range(len(N)):
        J=tangent['f_X'][i]+np.outer(tangent['f_lnd'][i],grad['rho_X'][i])
        ft=tangent['f_lnT'][i]+tangent['f_lnd'][i]*grad['rho_lnT'][i]
        nx=tangent['loss_X'][i]+tangent['loss_lnd'][i]*grad['rho_X'][i]
        nt=tangent['loss_lnT'][i]+tangent['loss_lnd'][i]*grad['rho_lnT'][i]
        ux=grad['energy_X'][i]/capacity[i];J-=np.outer(ft,ux);nx-=nt*ux
        W=np.zeros((28,28));W[:26,:26]=N[i]*J;W[:26,26]=N[i]*ft
        W[26,:26]=-N[i]*(rest@J+nx)/capacity[i];W[26,26]=-N[i]*(rest@ft+nt)/capacity[i]
        W[27,:26]=N[i]**2*nx/capacity[i];W[27,26]=N[i]**2*nt/capacity[i]
        f=source['dxdt'][i];loss=source['neutrino'][i]+thermal[i]
        F=np.r_[N[i]*f,N[i]*(-rest@f-loss-qflux[i])/capacity[i],N[i]**2*loss/capacity[i]]
        aug=np.zeros((29,29));aug[:28,:28]=W;aug[:28,-1]=F
        increments.append(expm(dt*aug)[:28,-1]);matrices.append(W);vectors.append(F)
    delta=np.array(increments);xx,total,projection=accumulate(data['X'],np.zeros_like(data['X'],dtype=np.longdouble),delta[:,:26])
    dx=old['X'].astype(np.longdouble)-data['X'].astype(np.longdouble)
    oldtheta=np.asarray((-dx@rest.astype(np.longdouble)-old['loss_per_baryon_gram']-dt*N*qflux)/capacity,dtype=float)
    oldloss=N*old['loss_per_baryon_gram']/capacity
    xs=float(max((abs(xx-old['X'])/(1e-16+1e-3*abs(xx-data['X']))).ravel()))
    hs=float(data['dm']@(N*capacity*abs(delta[:,26]-oldtheta))/release)
    ns=float(data['dm']@(capacity*abs(delta[:,27]-oldloss))/release)
    np.savez_compressed(OUT/'enthalpy-coordinate-predictor.npz',W=np.array(matrices),F0=np.array(vectors),delta=delta,
        input_X=xx,accumulated_X=total,capacity=capacity,reference_pressure=data['logP'])
    result=dict(classification='Counterexample candidate',composition_score=xs,
        unsigned_thermal_energy_difference_score=hs,unsigned_neutrino_energy_difference_score=ns,
        maximum_thermal_coordinate_difference=float(max(abs(delta[:,26]-oldtheta))),
        projection=projection,passed=bool(xs<1 and hs<1e-6 and ns<1e-6),
        original_failed_GR_verdict_preserved=True,corrected_stellar_path_completed=False)
    save('enthalpy-coordinate-predictor.json',result);print('ENTHALPY COORDINATE PREDICTOR',result,flush=True)


def residual_predictor_source():
    path=OUT/'P-1-0-state-17-4.npz';data=dict(np.load(path))
    save('residual-predictor-plan.json',dict(classification='Counterexample candidate',
        predictor_state_sha256=c.sha(path),initial_tangent_sha256=c.sha(OUT/'source-tangent.npz'),
        target='Evaluate the complete corrected nonlinear reaction/loss source at the saved first GR predictor for a reference-pressure enthalpy residual correction.',
        proposed_chart='Composition, reference-pressure enthalpy increment divided by the initial enthalpy log-temperature derivative, and neutrino energy at infinity divided by the same scale.',
        source='Same unchanged native executable and declared Request30 four-weak correction; fresh seeded common-EOS auxiliary inputs.',
        boundary='The existing predictor is a failed one-step control. Reusing its physical endpoint to evaluate a residual does not change that verdict. A new complete corrected path requires its own source/coordinate consistency, GR projection and time/energy tests.',
        corrected_stellar_path_completed=False))
    evaluate('residual-predictor',data)
    save('residual-predictor-source.json',dict(classification='Counterexample candidate',completed=True,
        original_predictor_sha256=c.sha(path),new_source_sha256=c.sha(OUT/'residual-predictor-corrected.npz'),
        corrected_stellar_path_completed=False))
    print('PASS fresh physical predictor source for residual correction',flush=True)


def residual_corrected_trial():
    """One new ETD residual trial; a full refinement series remains separate."""
    from scipy.linalg import expm
    assert json.loads((OUT/'enthalpy-coordinate-predictor.json').read_text())['passed']
    chart=json.loads((OUT/'residual-chart-audit.json').read_text());assert chart['complete'] and chart['passed']
    data=dict(np.load(OUT/'initial-state.npz'));predictor=dict(np.load(OUT/'enthalpy-coordinate-predictor.npz'))
    grad=dict(np.load(OUT/'initial-eos-directions.npz'));F1=np.load(OUT/'residual-chart-source.npz')['F']
    dt=json.loads((OUT/'plan.json').read_text())['time']['coordinate_seconds'];eos=v.ColdEOS()
    save('residual-corrected-plan.json',dict(classification='Counterexample candidate',dt=dt,
        inputs_sha256={name:c.sha(OUT/name) for name in ['enthalpy-coordinate-predictor.npz','residual-chart-source.npz','initial-state.npz']},
        formula='d_corrected=d+h phi2(hW)[F_at_GR_predictor-F0-Wd]',
        predictor='Reuse the original failed GR endpoint only after the separate coordinate covariance gate passed. That finite reuse error is recorded and is not a rigorous stage-error bound.',
        energy_coordinates='Reference-pressure enthalpy and neutrino energy at infinity, both divided by C0. Invert the corrected reference enthalpy using the existing EOS inverse budget, then solve the GR constraints.',
        criteria=dict(composition_minimum_before_tiny_projection=-1e-22,neutrino_loss_minimum_erg_g=-1e-12,
            EOS_inverse='max(2 erg/g,32 ulp(target),1e-8 abs(target-initial))',global_energy_relative_to_release=1e-6),
        time_refinement_completed=False,full_physical_GR_evolution=False))
    delta=[]
    for W,F0,d,F in zip(predictor['W'],predictor['F0'],predictor['delta'],F1):
        remainder=F-F0-W@d;aug=np.zeros((30,30));aug[:28,:28]=W;aug[:28,28]=remainder;aug[28,29]=1/dt
        delta.append(d+expm(dt*aug)[:28,29])
    delta=np.array(delta);np.savez_compressed(OUT/'residual-corrected-increments.npz',delta=delta)
    try:
        xx,total,projection=accumulate(data['X'],np.zeros_like(data['X'],dtype=np.longdouble),delta[:,:26])
    except AssertionError as error:
        save('residual-corrected-failure.json',dict(classification='Counterexample candidate',stage='composition',reason=str(error)))
        raise
    loss=grad['capacity']*delta[:,27]
    if loss.min()<-1e-12:
        save('residual-corrected-failure.json',dict(classification='Counterexample candidate',stage='neutrino',minimum=float(loss.min())))
        raise ValueError(('negative residual-corrected neutrino energy',loss.min()))
    target=grad['energy']+grad['capacity']*delta[:,26];temperatures=[];densities=[];scores=[]
    guess=data['lnT']+delta[:,26]-np.sum(grad['energy_X']*delta[:,:26],axis=1)/grad['capacity']
    for i in range(len(xx)):
        budget=max(2.,32*np.spacing(abs(target[i])),abs(target[i]-grad['energy'][i])*1e-8)
        defect,t,a=enthalpy_inverse(eos,data['logP'][i],xx[i],target[i],guess[i],budget)
        if abs(defect)>budget:
            save('residual-corrected-failure.json',dict(classification='Counterexample candidate',stage='EOS inverse',cell=i,defect=float(defect),budget=float(budget)))
            raise ValueError(('residual-corrected EOS inverse',i,defect,budget))
        temperatures.append(t);densities.append(np.log(a[0]));scores.append(abs(defect)/budget)
    reference={**data,'X':xx,'accumulation_base':data['X'],'accumulated_X':total,
        'lnT':np.array(temperatures),'lnd':np.array(densities),'reference_enthalpy':target,
        'chart_increment':delta,'neutrino_infinity_per_baryon_gram':loss}
    np.savez_compressed(OUT/'residual-corrected-reference.npz',**reference)
    save('residual-corrected-source.json',dict(classification='Counterexample candidate',maximum_inverse_score=float(max(scores)),projection=projection,
        neutrino_infinity_erg=float(data['dm']@loss),new_reference_state_completed=True,GR_projection_completed=False))
    parameters=json.loads((OUT/'P-1-0-structure-17-4.json').read_text())['parameters']
    state,record=structure('residual-corrected',reference,parameters)
    rest=(c.W/c.A-1)*(c.gr.C*100)**2
    release=float(data['dm'].astype(np.longdouble)@abs((xx.astype(np.longdouble)-data['X'].astype(np.longdouble))@rest.astype(np.longdouble)))
    r0=data['radius_faces_m'][0]*100;r1=state['radius_faces_m'][0]*100
    work=np.exp(float(data['boundary_logP']))*4*np.pi*(r1-r0)*(r1*r1+r1*r0+r0*r0)/3
    flux=np.load(OLD/'registered-face-flux.npz')['canonical'];expected=-flux[0]*dt-float(data['dm']@loss)-work
    measured=c.stable_mass_change(data,state)*(c.gr.C*100)**2;score=abs(measured-expected)/release
    result=dict(classification='Counterexample candidate',structure=record,energy_score=float(score),energy_passed=bool(score<1e-6),
        mass_change_energy_erg=measured,expected_energy_erg=float(expected),absolute_rest_release_erg=release,
        maximum_inverse_score=float(max(scores)),new_forced_quasistatic_trial_completed=True,
        time_refinement_completed=False,original_failed_paths_preserved=True,
        physical_EOS_certified=False,full_GR_thermal_fluid_metric_evolution=False,nonlinear_observation_completed=False)
    save('residual-corrected-trial.json',result);print('RESIDUAL CORRECTED TRIAL',result,flush=True)


def preconditioner_consistency():
    import sympy as sy
    from scipy.linalg import expm
    h,lam,mu,y=sy.symbols('h lam mu y',positive=True)
    step=y*(1+lam/mu*(sy.exp(-mu*h)-1))
    assert sy.simplify(step.subs(h,0)-y)==0
    assert sy.simplify(sy.diff(step,h).subs(h,0)+lam*y)==0
    # The same augmented exponential with a deliberately inexact operator.
    errors=[];duration=2.;rate=.7;operator=.3
    for count in [16,32,64,128]:
        value=.8;dt=duration/count
        for _ in range(count):
            matrix=np.array([[-operator,-rate*value],[0.,0.]])
            value+=expm(dt*matrix)[0,1]
        errors.append(abs(value-.8*np.exp(-rate*duration)))
    ratios=np.array(errors[:-1])/errors[1:];assert min(ratios)>1.9 and max(ratios)<2.2,ratios
    sdata=dict(classification='Proven',identity_first_derivative=True,
        linear_control_errors=errors,successive_error_ratios=ratios.tolist(),passed=True,
        mathematical_statement='For a fixed finite matrix W, h*phi_1(hW)*F(y)=h*F(y)+O(h^2); choosing W distinct from the exact Jacobian does not by itself destroy consistency.',
        limits='This gives neither a uniform stiff error bound nor stability/positivity for arbitrary W or step size. It does not repair any failed derivative certificate. The nonlinear GR trajectories retain independent energy, composition and temperature gates.',
        failure_policy='The failed small-step all-zone H2 estimator is not adopted. The pilot preserves the original numerical operator in its own snapshot; the final source-tangent record identifies the separately gated H2 chain replacement. Neither is a continuous-domain derivative certificate.')
    save('preconditioner-consistency.json',sdata);print('PRECONDITIONER CONSISTENCY',errors,ratios,flush=True)


def structure(label,data,parameters):
    oldout,oldeos=c.OUT,c.EOS
    try:
        c.OUT=OUT;c.EOS=v.ColdEOS
        for name in ['atomic-binding.json','freeeos-weights.json']:
            if not (OUT/name).exists(): shutil.copy2(v.OLD/name,OUT/name)
        # The legacy CX metadata was inherited unchanged as X evolved. Store
        # exactly the coefficient this existing GR solver actually computes.
        data={**data,'CX':(data['X']/c.A)@c.W}
        return c.structure(label,data,17,4,parameters)
    finally: c.OUT,c.EOS=oldout,oldeos


def evolution():
    plan=json.loads((OUT/'plan.json').read_text());assert (OUT/'source-tangent.npz').exists()
    assert json.loads((OUT/'source-control.json').read_text())['passed']
    assert json.loads((OUT/'preconditioner-consistency.json').read_text())['passed']
    original=dict(np.load(OUT/'initial-state.npz'));tangent=dict(np.load(OUT/'source-tangent.npz'))
    original['CX']=(original['X']/c.A)@c.W
    flux=np.load(OLD/'registered-face-flux.npz')['canonical'];dtall=plan['time']['coordinate_seconds']
    previous=None;rows=[]
    for steps in plan['time']['step_counts']:
        state={k:v.copy() for k,v in original.items()};loss_total=0.;release=0.;records=[]
        parameters=json.loads((v.OLD/'restored-structure-17-4.json').read_text())['parameters']
        for step in range(steps):
            label=f'P-{steps}-{step}';input_label='initial' if step==0 else label
            _,source=evaluate(input_label,state)
            thermal=thermal_loss(input_label)
            if step==0: grad=dict(np.load(OUT/'initial-eos-directions.npz'))
            else: grad=eos_directions(state)
            next_state,loss=source_step(state,source,thermal,tangent,grad,dtall/steps,flux,label)
            loss_total+=float(state['dm']@(np.exp(state['nu'])*loss))
            release+=json.loads((OUT/f'{label}-source.json').read_text())['absolute_rest_release_erg']
            state,record=structure(label,next_state,parameters);parameters=record['parameters'];records.append(record)
        mass_energy=c.stable_mass_change(original,state)*(c.gr.C*100)**2
        r0=original['radius_faces_m'][0]*100;r1=state['radius_faces_m'][0]*100
        work=np.exp(float(original['boundary_logP']))*4*np.pi*(r1-r0)*(r1*r1+r1*r0+r0*r0)/3
        expected=-flux[0]*dtall-loss_total-work;score=abs(mass_energy-expected)/release
        row=dict(steps=steps,records=records,energy_score=score,energy_passed=bool(score<1e-6),
            mass_change_energy_erg=mass_energy,expected_energy_erg=expected,release_erg=release,
            integrated_neutrino_infinity_erg=loss_total,pressure_work_erg=work)
        if previous is not None:
            sx=float(np.max(abs(state['X']-previous['X'])/(1e-16+1e-3*abs(state['X']-original['X']))))
            st=float(np.max(abs(state['lnT']-previous['lnT'])))
            row['refinement']=dict(composition_score=sx,lnT_difference=st,passed=bool(sx<1 and st<2e-6))
        rows.append(row);previous=state
        save('evolution.json',dict(classification='Counterexample candidate',rows=rows,
            complete=len(rows)==len(plan['time']['step_counts']),fixed_initial_tangent=True,
            fixed_flux=True,physical_EOS_certified=False,full_GR_evolution=False))
    print('EVOLUTION',steps,'energy',score,'refinement',row.get('refinement'),flush=True)


def refinement_plan():
    save('composition-refinement-plan.json',dict(classification='Counterexample candidate',
        observation='He3 has a 0.04914 first-order step difference at 1e-6/5e-7 while the extrapolated derivative is close to the returned value. The original finite-step gate fails and is preserved.',
        steps=[1e-8,5e-9],criterion=1e-3,
        method='Reevaluate failed directions at a hundredfold smaller pair of positive simplex steps. Require the same two-step relative criterion; save a separate derivative and bind it to the original derivative SHA. This is an explicit refinement, not a reinterpretation of the old pair.',
        automatic_selection='Only a newly measured passing refined direction may replace its failed original direction in the numerical preconditioner.',
        continuum_certificate=False))


def refine_he3():
    refine_direction('he3')


def refine_direction(species):
    plan=json.loads((OUT/'composition-refinement-plan.json').read_text())
    base=dict(np.load(OUT/'initial-state.npz'));reference=dict(np.load(OUT/'initial-corrected.npz'))
    j=c.NAMES.index(species);pivot=c.NAMES.index('he4');active=abs(reference['heat'])>1;sides=[]
    for h in plan['steps']:
        data={k:v.copy() for k,v in base.items()};data['X'][:,j]+=h;data['X'][:,pivot]-=h
        _,result=evaluate(f'X2-{species}-{h}',data)
        sides.append({key:(result[key]-reference[key])/h for key in reference})
    direction={key:2*sides[1][key]-sides[0][key] for key in reference}
    native=np.load(OUT/'initial-native.npz');reported=native['jacobian'][:,:,j]-native['jacobian'][:,:,pivot]
    scale=np.maximum(1e-30,np.maximum(abs(direction['dxdt']).max(1),abs(reported).max(1)))
    discrepancy=float(max((abs(direction['dxdt']-reported).max(1)/scale)[active]))
    difference=float(max((abs(sides[1]['dxdt']-sides[0]['dxdt']).max(1)/scale)[active]))
    np.savez_compressed(OUT/f'refined-{species}-derivative.npz',**direction)
    result=dict(classification='Counterexample candidate',species=species,steps=plan['steps'],
        original_derivative_sha256=c.sha(OUT/f'X-{species}-derivative.npz'),
        returned_Jacobian_score=discrepancy,finite_step_difference_score=difference,
        returned_Jacobian_passed=discrepancy<1e-3,finite_step_passed=difference<1e-3,
        original_failed_verdict_preserved=True,continuum_certificate=False)
    save(f'refined-{species}.json',result);print('REFINED COMPOSITION',result,flush=True)


def neutrino_boundary():
    sources=['neu/public/neu_lib.f90','star/private/neu.f90']
    for rel in sources:
        dest=OUT/'sources'/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(c.fresh.MESA/rel,dest)
    _,base=c.mesa(v.OLD/'EOS-pass-through-profile.data.gz');_,changed=c.mesa(v.OLD/'EOS-common-profile.data.gz')
    same=np.array_equal(base['non_nuc_neu'],changed['non_nuc_neu']);assert same
    save('neutrino-boundary.json',dict(classification='Counterexample candidate',
        archived_inputs=['T','log10_T','rho','log10_rho','abar','zbar','z2bar','log10_Tlim','flags'],
        independent_eta_or_free_electron_input=False,
        same_state_nuclear_EOS_auxiliary_intervention_leaves_thermal_loss_bitwise_equal=same,
        cells=len(base['non_nuc_neu']),
        source_sha256={rel:c.sha(OUT/'sources'/rel) for rel in sources},
        source_caveat='Source signatures alone do not certify compilation; the recorded independent same-state intervention supplies the stated functional evidence.',
        physical_fit_error_certified=False))
    print('PASS thermal-neutrino input boundary',len(base['non_nuc_neu']),flush=True)


def baryon_face_diffusion(data,opacity):
    """Interior outward L_infinity faces, with baryon mass increasing outward.

    Counterexample candidate: the caller supplies opacity per baryon gram.
    No convection, surface atmosphere or flux inertia is silently included.
    """
    dm=data['dm'];theta_log=data['lnT']+data['nu'];nf=np.exp(data['nu_faces'][1:-1])
    weight=dm[:-1]/(dm[:-1]+dm[1:]);kap=weight*opacity[1:]+(1-weight)*opacity[:-1]
    assert np.all(dm>0) and np.all(kap>0)
    area=4*np.pi*(data['radius_faces_m'][1:-1]*100)**2
    difference=np.exp(4*theta_log[1:])*np.expm1(4*(theta_log[:-1]-theta_log[1:]))
    return -area**2*(4*5.670400e-5)/(3*kap*nf*nf)*difference/((dm[:-1]+dm[1:])/2)


def transport_boundary_audit():
    import sympy as sy
    data=dict(np.load(OUT/'initial-state.npz'));_,p=c.mesa(OUT/'initial-profile.data.gz')
    for rel in ['star/private/star_utils.f90','star/private/profile_getval.f90']:
        target=OUT/'sources'/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(c.fresh.MESA/rel,target)
    getter=(OUT/'sources/star/private/star_utils.f90').read_text()
    exporter=(OUT/'sources/star/private/profile_getval.f90').read_text()
    assert 'val = (s% L(k) - L_rad)/max(1d0,s% L(k))' in exporter
    assert 'del_m = s% dm_bar(j)' in getter and 'del_T4 = pow4(s% T(j-1)) - pow4(s% T(j))' in getter
    L=p['luminosity']*c.fresh.LSUN
    exported_rad=L-p['lum_conv_div_L']*np.maximum(1.,L)
    # Reconstruct the archived getter on the actual native face grid first.
    flat=dict(dm=p['dm'],lnT=p['logT']*np.log(10.),nu=np.zeros(len(L)),
        nu_faces=np.zeros(len(L)+1),radius_faces_m=np.r_[p['radius']*c.gr.RSUN,0.])
    native_interior=baryon_face_diffusion(flat,p['opacity'])
    native=np.r_[native_interior[0],native_interior] # getter uses face 2 at the outermost output
    native_score=abs(native-exported_rad)/np.maximum(1.,np.maximum(abs(native),abs(exported_rad)))
    if not (OUT/'transport-getter-control.json').exists():
        save('transport-getter-control.json',dict(classification='Counterexample candidate',
            maximum_relative_difference=float(native_score.max()),worst_cell=int(np.argmax(native_score)),
            passed=bool(native_score.max()<1e-7),original_ratio_is_signed_L_with_floor_one=True))
    # The native final thermal node is placed at the centre: its dm_bar is
    # dm[-2]/2 + dm[-1], unlike this GR solver's baryon midpoint nodes.
    assert 'dm_bar(nz) = 0.5d0*dm(nz-1) + dm(nz)' in getter
    native[-1]*=(p['dm'][-2]+p['dm'][-1])/(p['dm'][-2]+2*p['dm'][-1])
    native_score=abs(native-exported_rad)/np.maximum(1.,np.maximum(abs(native),abs(exported_rad)))
    save('transport-getter-boundary-control.json',dict(classification='Counterexample candidate',
        maximum_relative_difference=float(native_score.max()),worst_cell=int(np.argmax(native_score)),
        passed=bool(native_score.max()<1e-7),original_failure_sha256=c.sha(OUT/'transport-getter-control.json'),
        correction='Reproduce the source-defined special central dm_bar only for the native getter audit; the GR baryon-midpoint stencil is unchanged.'))
    assert native_score.max()<1e-7,('native radiative face getter',native_score.max(),np.argmax(native_score))
    gr_flux=baryon_face_diffusion(data,p['opacity']);registered=np.load(OLD/'registered-face-flux.npz')['canonical']
    meaningful=abs(L)>1e-10*c.fresh.LSUN
    proper_fraction=np.divide(L-exported_rad,abs(L),out=np.zeros_like(L),where=meaningful)
    dm=data['dm'];mass=float(dm.sum());score=abs(gr_flux-registered[1:-1])/np.maximum(abs(registered[1:-1]),1e-10*c.fresh.LSUN)
    # The chain rule removes sqrt(f) in the baryon-mass form of Fourier's law.
    area,N,T,K,rho,f,dtheta=sy.symbols('area N T K rho f dtheta',positive=True)
    dm_dr=area*rho/sy.sqrt(f)
    local=-K*sy.sqrt(f)/N*dtheta*dm_dr
    assert sy.simplify(area*N**2*local+area**2*N*K*rho*dtheta)==0
    # Constant redshifted temperature is an exact zero-flux discrete control.
    tolman={**data,'lnT':-data['nu']}
    assert np.all(baryon_face_diffusion(tolman,p['opacity'])==0)
    save('transport-boundary-audit.json',dict(classification='Counterexample candidate',
        cells=len(L),negative_luminosity_cells=int((L<0).sum()),
        negative_luminosity_mass_fraction=float(dm@(L<0)/mass),
        raw_ratio_over_one_pct_mass_fraction=float(dm@(abs(p['lum_conv_div_L'])>.01)/mass),
        meaningful_signed_residual_over_one_pct_mass_fraction=float(dm@(meaningful&(abs(proper_fraction)>.01))/mass),
        native_convective_flag_mass_fraction=float(dm@(p['mixing_type']==1)/mass),
        native_overshoot_flag_mass_fraction=float(dm@(p['mixing_type']==2)/mass),
        GR_face_diffusion_mismatch_quantiles=np.quantile(score,[0,.5,.9,.99,1]).tolist(),
        native_getter_control_passed=True,all_transport_closed=False,
        statement='The imported luminosity and temperature gradient need not obey its MLT closure. A signed-L denominator clipped at one erg/s cannot be treated as a bounded convective fraction for negative luminosity. Native mixing flags are diagnostics on the imported state, not a common-EOS GR convection solution.',
        approximation='Actual state opacity evaluated by the unchanged native module; baryon-gram opacity convention assumed for this GR face control. No new physical opacity error bound, convection or atmosphere solution.',
        source_sha256={rel:c.sha(OUT/'sources'/rel) for rel in ['star/private/star_utils.f90','star/private/profile_getval.f90']}))
    np.savez_compressed(OUT/'transport-boundary-audit.npz',native_getter=native,exported_radiative_luminosity=exported_rad,
        signed_flux_residual_fraction=proper_fraction,meaningful_luminosity=meaningful,GR_interior_face_Linf=gr_flux,
        registered_interior_face_Linf=registered[1:-1],GR_diffusion_mismatch_score=score)
    save('transport-face-symbolic.json',dict(classification='Proven',baryon_coordinate_chain_rule=True,
        discrete_constant_redshifted_temperature_zero_flux=True,
        scope='Specified spherical quasistatic Fourier law and the declared positive-opacity face stencil. No full dynamical GR or physical opacity theorem.'))
    print('TRANSPORT BOUNDARY',json.loads((OUT/'transport-boundary-audit.json').read_text()),flush=True)


def thermal_refinement_plan():
    save('thermal-refinement-plan.json',dict(classification='Counterexample candidate',
        observed_failure='The original temperature pair has a 0.0098622 nuclear-neutrino derivative difference, worst cell 2274. Vector and heat directions pass. Preserve all three original verdicts.',
        diagnostic='The larger h=1e-4 stencil crosses a weak blend boundary near cell 2274; Request30 verified that the smaller 5e-5 stencils stay on one branch at all burning cells.',
        direction='lnT',steps=[5e-5,2.5e-5],criterion=1e-3,
        method='Add a smaller central stencil and compare all three output derivatives with the existing 5e-5 stencil. Do not raise the threshold or claim a derivative through the kink.',
        continuous_domain_certified=False))


def refine_temperature():
    plan=json.loads((OUT/'thermal-refinement-plan.json').read_text());field=plan['direction'];base=dict(np.load(OUT/'initial-state.npz'))
    reference=dict(np.load(OUT/'initial-corrected.npz'));active=abs(reference['heat'])>1;results=[]
    for h in plan['steps']:
        path=OUT/f'{field}-{h}-derivative.npz'
        if path.exists(): result=dict(np.load(path))
        else:
            sides=[]
            for sign in [-1,1]:
                data={k:v.copy() for k,v in base.items()};data[field]+=sign*h
                sides.append(evaluate(f'{field}-{h}-{sign}',data)[1])
            result={key:(sides[1][key]-sides[0][key])/(2*h) for key in reference}
            np.savez_compressed(path,**result)
        results.append(result)
    rows=[]
    for key in reference:
        a,b=results;scale=np.maximum(1e-30 if key=='dxdt' else 1,
            np.maximum(abs(a[key]).max(1),abs(b[key]).max(1)) if key=='dxdt' else np.maximum(abs(a[key]),abs(b[key])))
        difference=abs(b[key]-a[key]);score=(difference.max(1) if key=='dxdt' else difference)/scale
        rows.append(dict(output=key,step_difference_score=float(max(score[active])),passed=bool(max(score[active])<1e-3)))
    save('refined-lnT.json',dict(classification='Counterexample candidate',steps=plan['steps'],rows=rows,
        passed=all(r['passed'] for r in rows),original_failure_preserved=True,continuous_domain_certified=False))
    print('REFINED TEMPERATURE',rows,flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
