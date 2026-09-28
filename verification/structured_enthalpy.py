"""Request32: explicit GR stages in a reference-pressure enthalpy chart.

Counterexample candidate. Fixed flux and an initial block preconditioner remain
declared controls; they do not constitute physical EOS/transport/GR certification.
"""
from pathlib import Path
import hashlib, json, shutil, sys
import numpy as np
from scipy.linalg import expm
import conservative_star as s

c=s.c;ROOT=s.ROOT;OLD=s.OUT
OUT=ROOT/'outputs/structured-enthalpy32'
CACHE=Path('/home/lpaiu/work/structured-enthalpy32')


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def configure():
    s.OUT=OUT;s.CACHE=CACHE


def prepare():
    assert not OUT.exists() and not CACHE.exists();OUT.mkdir();CACHE.mkdir()
    prior=json.loads((OLD/'plan.json').read_text())
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='c39152f',
        previous_manifest_sha256=c.sha(OLD/'manifest.json'),
        full_goal='Resolve physical EOS certification, continuous errors, self-consistent GR thermal/fluid/metric/atmosphere evolution, nonzero scalar drive/charge and complete nonlinear observational inference. These finite time controls do not redefine completion.',
        duration_coordinate_seconds=prior['time']['coordinate_seconds'],step_counts=[1,2,4],
        controls=prior['controls'],
        coordinates='26 mass fractions, enthalpy at a fixed initial baryon-shell reference pressure divided by C0, and accumulated neutrino energy at infinity divided by C0. C0 is the initial positive fixed-P dH/dlnT.',
        integration='For a fixed initial block matrix W, d=h phi1(hW) F(y). Reconstruct the predictor reference EOS and solve its GR constraints directly, evaluate its complete chart source, then correct by h phi2(hW)[F(y+d)-F(y)-Wd]. Reconstruct and solve the corrected state independently.',
        predictor_policy='No substitution of an older failed GR endpoint for the new predictor. Request31 covariance and scalar order controls remain supporting evidence only.',
        source='Unchanged native executable, fresh seeded common-EOS auxiliaries, Request30 four-weak double-table correction, two-size finite chemical difference at fixed reference/actual P,T.',
        chemical_steps=[1e-5,5e-6],chemical_source_criterion=1e-3,
        EOS_inverse='max(2 erg/g,32 ulp(target),1e-8 abs(target-current_reference_enthalpy))',
        fixed_flux=True,fixed_initial_preconditioner=True,
        energy_release_denominator='Sum over accepted time steps and baryon cells of the absolute neutral-atom rest-excess change; no denominator from a predictor stage.',
        policy='Preserve every earlier failure and every new failed gate. A finite energy/time pass is not a uniform error enclosure, physical EOS certificate or full GR/observational closure.'))
    history={}
    for rel in ['docs/'+name+'.md' for name in c.cell.DOCS]+['paper/revision-manifest.json']:
        path=ROOT/rel;target=OUT/'previous-notes'/path.name;target.parent.mkdir(exist_ok=True);shutil.copy2(path,target)
        history[rel]=dict(snapshot=target.relative_to(ROOT).as_posix(),sha256=c.sha(path),historical_manifest=rel.startswith('paper/'))
    save('historical-note-bindings.json',history)
    for name in ['atomic-binding.json','freeeos-weights.json']:
        shutil.copy2(s.v.OLD/name,OUT/name)
    initial=dict(np.load(OLD/'initial-state.npz'));coefficient=(initial['X']/c.A)@c.W
    save('initial-metadata.json',dict(classification='Counterexample candidate',
        inherited_CX_maximum_difference=float(max(abs(initial['CX']-coefficient))),
        change='Refresh only CX metadata from the same X and the same neutral-atom weights. Physical state and baryon inventories are unchanged.'))
    initial['CX']=coefficient;np.savez_compressed(OUT/'initial-state.npz',**initial)
    old=dict(np.load(OLD/'enthalpy-coordinate-predictor.npz'));grad=dict(np.load(OLD/'initial-eos-directions.npz'))
    assert np.array_equal(old['capacity'],grad['capacity']) and np.all(grad['capacity']>0)
    np.savez_compressed(OUT/'basis.npz',W=old['W'],F0_previous=old['F0'],capacity=grad['capacity'],
        reference_pressure=old['reference_pressure'],initial_H=grad['energy'],initial_H_X=grad['energy_X'])
    save('basis-binding.json',dict(classification='Imported from prior work',
        source_sha256={name:c.sha(OLD/name) for name in ['initial-state.npz','enthalpy-coordinate-predictor.npz','initial-eos-directions.npz','structural-split-control.json','enthalpy-chart-control.json']},
        fixed_numeric_preconditioner=True,continuous_source_Jacobian_certified=False))


def fingerprint(*states):
    digest=hashlib.sha256()
    for state in states:
        for key in ['X','lnd','lnT','logP','nu','dm']:
            a=np.ascontiguousarray(state[key]);digest.update(key.encode()+str(a.shape).encode()+a.dtype.str.encode()+a.tobytes())
    return digest.hexdigest()


def phi_increment(W,F,dt,second=False):
    n=len(F);a=np.zeros((n+1+int(second),n+1+int(second)))
    a[:n,:n]=W;a[:n,n]=F
    if second: a[n,n+1]=1/dt
    return expm(dt*a)[:n,-1]


def control():
    # Actual block exponential and correction, deliberately inexact W.
    rate=.7;beta=.2;start=.8;duration=2.;W=np.array([[-rate]])
    exact=start*np.exp(-rate*duration)/(1-beta*start/rate*(1-np.exp(-rate*duration)))
    rows=[]
    for count in [16,32,64,128]:
        y=np.array([start]);h=duration/count
        f=lambda x:-rate*x+beta*x*x
        for _ in range(count):
            F=f(y);d=phi_increment(W,F,h)
            y+=d+phi_increment(W,f(y+d)-F-W@d,h,True)
        rows.append(dict(steps=count,error=float(abs(y[0]-exact))))
    ratios=[a['error']/b['error'] for a,b in zip(rows,rows[1:])]
    assert min(ratios)>3.8 and max(ratios)<4.2,ratios
    assert np.allclose(phi_increment(np.zeros((2,2)),np.array([2.,3.]),.4,True),[.4,.6],rtol=0,atol=1e-15)
    save('control.json',dict(classification='Proven',passed=True,nonlinear_scalar_errors=rows,ratios=ratios,
        scope='Implemented matrix-function/correction control with an analytic scalar solution. No claim of a uniform stiff bound or stellar convergence.'))
    print('PASS structured enthalpy scalar control',ratios,flush=True)


def chart_source(label,reference,physical,basis):
    configure();binding=fingerprint(reference,physical);path=OUT/f'{label}-chart.npz';report=OUT/f'{label}-chart.json'
    if path.exists() and report.exists():
        result=json.loads(report.read_text());assert result['state_binding']==binding and result['passed']
        return np.load(path)['F']
    _,source=s.evaluate(label,physical);thermal=s.thermal_loss(label)
    flux=np.load(s.OLD/'registered-face-flux.npz')['canonical'];eos=s.v.ColdEOS()
    rest=(c.W/c.A-1)*(c.gr.C*100)**2;pivot=c.NAMES.index('he4');rows=[];vectors=[]
    plan=json.loads((OUT/'plan.json').read_text());sizes=plan['chemical_steps']
    for i,x in enumerate(physical['X']):
        assert np.array_equal(x,reference['X'][i])
        tr=np.exp(reference['lnT'][i]);ta=np.exp(physical['lnT'][i]);ratio=tr/ta
        f=source['dxdt'][i];loss=source['neutrino'][i]+thermal[i];N=np.exp(physical['nu'][i])
        qflux=(flux[i]-flux[i+1])/(physical['dm'][i]*N*N)
        if reference['lnT'][i]==physical['lnT'][i] and reference['logP'][i]==physical['logP'][i]:
            chemical=[0.,0.];ds=0.
        else:
            ar=eos(1,reference['logP'][i],reference['lnT'][i],x)
            aa=eos(1,physical['logP'][i],physical['lnT'][i],x)
            assert ar[10]-ar[1]/ar[0]*ar[8]>0 and aa[10]-aa[1]/aa[0]*aa[8]>0,(label,i,'nonpositive heat capacity')
            hr=ar[2]+ar[1]/ar[0];ha=aa[2]+aa[1]/aa[0];ds=aa[3]-ar[3];chemical=[]
            for h in sizes:
                derivative=np.zeros(26)
                for j in range(26):
                    if j==pivot: continue
                    z=x.copy();z[j]+=h;z[pivot]-=h;assert z.min()>=0
                    br=eos(1,reference['logP'][i],reference['lnT'][i],z)
                    ba=eos(1,physical['logP'][i],physical['lnT'][i],z)
                    derivative[j]=((br[2]+br[1]/br[0]-hr)-ratio*(ba[2]+ba[1]/ba[0]-ha)
                        +tr*((ba[3]-aa[3])-(br[3]-ar[3])))/h
                chemical.append(float(derivative@f))
        values=[N*(ratio*(-rest@f-loss-qflux)+term) for term in chemical]
        score=abs(values[1]-values[0])/max(1.,abs(values[0]),abs(values[1]))
        rows.append(dict(cell=i,finite_source_score=float(score),entropy_pair_residual=float(ds)))
        vectors.append(np.r_[N*f,(2*values[1]-values[0])/basis['capacity'][i],N*N*loss/basis['capacity'][i]])
        if i%1000==0: print('CHART SOURCE',label,i,flush=True)
    maximum=max(row['finite_source_score'] for row in rows);passed=maximum<plan['chemical_source_criterion']
    np.savez_compressed(path,F=np.array(vectors))
    save(report.name,dict(classification='Counterexample candidate',state_binding=binding,
        native_source_sha256=c.sha(OUT/f'{label}-corrected.npz'),rows=rows,maximum_source_score=maximum,
        passed=bool(passed),continuous_or_physical_certificate=False))
    assert passed,(label,'chemical chart source',maximum)
    return np.array(vectors)


def initial():
    assert json.loads((OUT/'control.json').read_text())['passed']
    base=dict(np.load(OUT/'initial-state.npz'));basis=dict(np.load(OUT/'basis.npz'))
    F=chart_source('initial',base,base,basis)
    previous=dict(np.load(OLD/'initial-corrected.npz'));current=dict(np.load(OUT/'initial-corrected.npz'))
    equal={k:bool(np.array_equal(previous[k],current[k])) for k in ['dxdt','heat','neutrino']}
    error=float(np.max(abs(F-basis['F0_previous'])/np.maximum(1e-30,abs(basis['F0_previous']))))
    passed=all(equal.values()) and error<1e-12
    save('initial-control.json',dict(classification='Counterexample candidate',native_corrected_bitwise_equal=equal,
        chart_relative_difference=error,passed=passed))
    assert passed,(equal,error);print('PASS fresh initial source and chart control',equal,error,flush=True)


def reference_state(label,base,current,previous_reference,total,increment,basis):
    target_total=total+increment.astype(np.longdouble)
    xx,accumulated,projection=s.accumulate(base['X'],total[:,:26],increment[:,:26])
    target_total[:,:26]=accumulated
    if np.min(target_total[:,27]*basis['capacity'])<-1e-12:
        raise ValueError(('negative integrated neutrino energy',label,float(np.min(target_total[:,27]*basis['capacity']))))
    path=OUT/f'{label}-reference.npz';report=OUT/f'{label}-reference.json'
    if path.exists() and report.exists():
        saved=dict(np.load(path));assert np.array_equal(saved['chart_total'],target_total) and np.array_equal(saved['X'],xx)
        assert json.loads(report.read_text())['previous_physical_binding']==fingerprint(current)
        return saved
    target=np.asarray(basis['initial_H'].astype(np.longdouble)+basis['capacity']*target_total[:,26],float)
    oldH=previous_reference.get('reference_enthalpy',basis['initial_H'])
    guess=previous_reference['lnT']+increment[:,26]-np.sum(basis['initial_H_X']*increment[:,:26],axis=1)/basis['capacity']
    temperatures=[];densities=[];energies=[];entropies=[];scores=[];eos=s.v.ColdEOS()
    for i in range(len(xx)):
        budget=max(2.,32*np.spacing(abs(target[i])),abs(target[i]-oldH[i])*1e-8)
        defect,t,a=s.enthalpy_inverse(eos,basis['reference_pressure'][i],xx[i],target[i],guess[i],budget)
        assert a is not None and a[10]-a[1]/a[0]*a[8]>0,(label,i,'nonpositive heat capacity')
        if abs(defect)>budget:
            save(f'{label}-failure.json',dict(classification='Counterexample candidate',stage='reference enthalpy inverse',cell=i,defect=float(defect),budget=float(budget)))
            raise ValueError((label,i,defect,budget))
        temperatures.append(t);densities.append(np.log(a[0]));energies.append(a[2]);entropies.append(a[3]);scores.append(abs(defect)/budget)
    result={**current,'X':xx,'lnd':np.array(densities),'lnT':np.array(temperatures),
        'logP':basis['reference_pressure'],'reference_enthalpy':target,'chart_total':target_total,
        'u_W':np.array(energies),'s_B':np.array(entropies),'CX':(xx/c.A)@c.W,
        'accumulation_base':base['X'],'accumulated_X':accumulated}
    np.savez_compressed(path,**result)
    save(report.name,dict(classification='Counterexample candidate',previous_physical_binding=fingerprint(current),
        maximum_inverse_score=float(max(scores)),projection=projection,passed=True))
    return result


def project(label,reference,parameters):
    configure();path=OUT/f'{label}-state-17-4.npz';report=OUT/f'{label}-structure-17-4.json'
    binding=OUT/f'{label}-projection-binding.json';digest=fingerprint(reference)
    if path.exists() and report.exists() and binding.exists():
        old=json.loads(binding.read_text())
        assert old['reference_binding']==digest and old['initial_parameters']==parameters
        assert old['boundary_logP']==float(reference['boundary_logP'])
        return dict(np.load(path)),json.loads(report.read_text())
    state,record=s.structure(label,reference,parameters)
    assert np.array_equal(state['X'],reference['X']) and np.array_equal(state['dm'],reference['dm'])
    save(binding.name,dict(classification='Counterexample candidate',reference_binding=digest,
        initial_parameters=list(parameters),boundary_logP=float(reference['boundary_logP']),
        state_sha256=c.sha(path),composition_and_baryon_masses_bitwise_unchanged=True))
    return state,record


def reference_control():
    base=dict(np.load(OUT/'initial-state.npz'));basis=dict(np.load(OUT/'basis.npz'))
    total=np.zeros((len(base['X']),28),dtype=np.longdouble);delta=np.zeros_like(total,dtype=float)
    state=reference_state('zero-control',base,base,base,total,delta,basis)
    same={key:bool(np.array_equal(state[key],base[key])) for key in ['X','lnT','dm']}
    assert all(same.values()) and np.array_equal(state['reference_enthalpy'],basis['initial_H']),same
    assert np.all(state['chart_total']==0)
    residual=state['u_W']+np.exp(state['logP']-state['lnd'])-state['reference_enthalpy']
    # Reconstructing density from its stored logarithm incurs its own rounding.
    assert np.all(abs(residual)<128*np.spacing(abs(state['reference_enthalpy']))+4),float(abs(residual).max())
    save('reference-control.json',dict(classification='Counterexample candidate',passed=True,bitwise_equal=same,
        maximum_lnd_change=float(abs(state['lnd']-base['lnd']).max()),
        maximum_reconstructed_enthalpy_residual_erg_g=float(abs(residual).max()),
        meaning='Zero-increment reference-coordinate reconstruction on all actual initial cells. This does not certify GR projection or continuous EOS error.'))
    print('PASS zero-increment reference reconstruction',same,flush=True)


def evolve(steps):
    configure();plan=json.loads((OUT/'plan.json').read_text());assert steps in plan['step_counts']
    assert json.loads((OUT/'initial-control.json').read_text())['passed']
    assert json.loads((OUT/'reference-control.json').read_text())['passed']
    base=dict(np.load(OUT/'initial-state.npz'));basis=dict(np.load(OUT/'basis.npz'))
    state={k:a.copy() for k,a in base.items()};reference=state.copy();total=np.zeros((len(base['X']),28),dtype=np.longdouble)
    F=np.load(OUT/'initial-chart.npz')['F'];dt=plan['duration_coordinate_seconds']/steps
    parameters=json.loads((s.v.OLD/'restored-structure-17-4.json').read_text())['parameters']
    release=0.;records=[];rest=(c.W/c.A-1)*(c.gr.C*100)**2
    for step in range(steps):
        label=f'E-{steps}-{step}'
        d=np.array([phi_increment(W,f,dt) for W,f in zip(basis['W'],F)])
        predref=reference_state(label+'-predictor',base,state,reference,total,d,basis)
        predicted,predrecord=project(label+'-predictor',predref,parameters)
        Fpred=chart_source(label+'-predictor',predref,predicted,basis)
        change=np.array([v+phi_increment(W,fp-f-W@v,dt,True) for W,v,fp,f in zip(basis['W'],d,Fpred,F)])
        corrected=reference_state(label+'-corrector',base,state,reference,total,change,basis)
        next_state,record=project(label+'-corrector',corrected,predrecord['parameters'])
        dx=next_state['X'].astype(np.longdouble)-state['X'].astype(np.longdouble)
        release+=float(base['dm'].astype(np.longdouble)@abs(dx@rest.astype(np.longdouble)))
        records.append(dict(step=step,predictor=predrecord,corrector=record))
        state=next_state;reference=corrected;total=corrected['chart_total'];parameters=record['parameters']
        save(f'path-{steps}-progress.json',dict(classification='Counterexample candidate',completed_steps=step+1,total_steps=steps,records=records))
        if step+1<steps: F=chart_source(label+'-corrector',reference,state,basis)
    loss=float(base['dm'].astype(np.longdouble)@(basis['capacity']*total[:,27]))
    r0=base['radius_faces_m'][0]*100;r1=state['radius_faces_m'][0]*100
    work=np.exp(float(base['boundary_logP']))*4*np.pi*(r1-r0)*(r1*r1+r1*r0+r0*r0)/3
    flux=np.load(s.OLD/'registered-face-flux.npz')['canonical'];expected=-flux[0]*plan['duration_coordinate_seconds']-loss-work
    measured=c.stable_mass_change(base,state)*(c.gr.C*100)**2;score=abs(measured-expected)/release
    save(f'path-{steps}.json',dict(classification='Counterexample candidate',steps=steps,completed=True,records=records,
        mass_change_energy_erg=measured,expected_energy_erg=float(expected),release_erg=release,
        neutrino_infinity_erg=loss,pressure_work_erg=float(work),energy_score=float(score),
        energy_passed=bool(score<plan['controls']['global_energy_relative_to_release']),
        fixed_flux=True,fixed_initial_preconditioner=True,physical_EOS_certified=False,full_GR_evolution=False))
    print('ENTHALPY PATH',steps,'energy',score,flush=True)


def collect():
    plan=json.loads((OUT/'plan.json').read_text());base=np.load(OUT/'initial-state.npz')['X'];rows=[];previous=None
    for steps in plan['step_counts']:
        path=OUT/f'path-{steps}.json'
        if not path.exists(): continue
        row=json.loads(path.read_text());assert row['completed']
        state=dict(np.load(OUT/f'E-{steps}-{steps-1}-corrector-state-17-4.npz'))
        if previous is not None:
            sx=float(np.max(abs(state['X']-previous['X'])/(1e-16+1e-3*abs(state['X']-base))))
            st=float(np.max(abs(state['lnT']-previous['lnT'])))
            row['refinement']=dict(composition_score=sx,lnT_difference=st,passed=bool(sx<1 and st<2e-6))
        rows.append(row);previous=state
    save('evolution.json',dict(classification='Counterexample candidate',rows=rows,
        complete=len(rows)==len(plan['step_counts']),fixed_flux=True,fixed_initial_preconditioner=True,
        physical_EOS_certified=False,full_GR_evolution=False))
    print('COLLECT',[(row['steps'],row['energy_score'],row.get('refinement')) for row in rows],flush=True)


if __name__=='__main__':
    if sys.argv[1]=='evolve': evolve(int(sys.argv[2]))
    else: globals()[sys.argv[1]]()
