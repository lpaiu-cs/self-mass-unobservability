"""Request 20: frozen-window thermal WD robustness; reuses Request 19 physics."""
from pathlib import Path
import gzip, hashlib, json, os, re, shutil, subprocess, sys, time
import numpy as np
from thermal_wd import mesa
from thermal_restart import sha

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'outputs/thermal-robustness20'
OLD = ROOT/'outputs/thermal-restart19'
CACHE = Path('/home/lpaiu/work/thermal-robustness20')
RUNTIME = Path('/home/lpaiu/work/thermal-restart19')
CASES = {'control': (-6, .8, .001, '2d-5'),
         'fast': (-5, .8, .001, '2d-5'),
         'slow': (-7, .8, .001, '2d-5'),
         'mesh': (-6, .4, .0005, '2d-5'),
         'time': (-6, .8, .001, '1d-5')}
PHASE_CASES = ['fast_phase','time_phase','fast_phase_replay','time_phase_replay']


def save(name, obj):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/name).write_text(json.dumps(obj, ensure_ascii=False, indent=2)+'\n')


def gzcopy(src, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    with src.open('rb') as inp, dest.open('wb') as out:
        with gzip.GzipFile(filename='', mode='wb', fileobj=out, mtime=0) as z:
            shutil.copyfileobj(inp, z)


def prepare():
    assert not any((CACHE/name).exists() for name in CASES)
    OUT.mkdir(parents=True, exist_ok=True); CACHE.mkdir(exist_ok=True)
    _, hist = mesa(OLD/'runs/mass18000/history.data.gz')
    age_end = float(hist['star_age'][-1])
    original = json.loads((OLD/'published-reference.json').read_text())
    age_start = original['history']['18000']['star_age']
    plan = dict(classification='Counterexample candidate',
        age_window_yr=[age_start, age_end], cases=CASES,
        mass_GM_sun=.197536385307, mass_relative_tolerance=1e-8,
        Teff_center_K=15800, Teff_halfwidth_K=300,
        logg_center=5.82, logg_halfwidth=.15,
        candidate_rule='Preserve every actual evolved row in the fixed age window satisfying all cuts; choose minimum squared standardized distance among passing rows, otherwise among mass-passing rows. No interpolated structure.',
        comparison_rule='Also compare each trajectory to baseline on 501 equally spaced common ages using linear history interpolation, and enumerate all Teff=15800 crossings. Interpolation is diagnostic only.',
        success_rule='Each of four variants must terminate normally at the frozen end age and contain at least one passing row. This is finite sensitivity testing, not continuum convergence or independent formation histories.',
        control_rule='First ten compact-output steps must reproduce the frozen baseline history and the shared profile columns exactly.',
        guards=dict(max_model_number=40000, timeout_s=7200),
        assumptions='One initial photo and reconstructed network; relax_mass resets age and is model construction, not an observed mass-loss history.',
        baseline_history_sha256=sha(OLD/'runs/mass18000/history.data.gz'),
        previous_manifest_sha256=sha(OLD/'manifest.json'))
    save('preregistered-plan.json', plan)
    shutil.copy2(OLD/'source-bindings.json', OUT/'source-bindings.json')
    baseline = RUNTIME/'mass18000'
    expected = json.loads((OLD/'mass18000-inputs.json').read_text())['sha256']
    for rel, digest in expected.items():
        # MESA advances .restart at runtime; restore the frozen input, not that cursor.
        path=OLD/'runs/mass18000/restart-input.txt' if rel=='.restart' else baseline/rel
        assert sha(path) == digest, rel
    cols = ['zone', 'logT', 'logRho', 'logP', 'logR', 'mass', 'radius', 'omega',
            'dm', 'energy', 'pressure', 'rho', 'm_grav', 'mass_correction_factor',
            'm_grav_div_m_baryonic', 'add_abundances']
    for name, (rate, mesh, dq, target) in CASES.items():
        folder = CACHE/name; folder.mkdir(); (folder/'photos1').mkdir()
        for p in baseline.iterdir():
            if p.is_file() and (p.name.startswith('inlist') or p.suffix in ('.list', '.net') or p.name=='binary'):
                shutil.copy2(p, folder/p.name)
        shutil.copy2(baseline/'photos1/18000', folder/'photos1/18000')
        (folder/'.restart').write_text('18000\n')
        p = folder/'inlist1'; s = p.read_text()
        for old, new in [('lg_max_abs_mdot = -6', f'lg_max_abs_mdot = {rate}'),
                         ('mesh_delta_coeff = 0.8', f'mesh_delta_coeff = {mesh}'),
                         ('max_dq = 0.001', f'max_dq = {dq}'),
                         ('max_model_number = 19500', f'max_model_number = {18010 if name=="control" else 40000}'),
                         ('max_age =1.4d10', f'max_age = {age_end:.17e}')]:
            assert s.count(old)==1, old
            s = s.replace(old, new)
        p.write_text(s)
        p = folder/'inlist_project'; s = p.read_text()
        s, count = re.subn(r'(varcontrol_(?:case_a|case_b|ms|post_ms)\s*=\s*)2d-5',
                          lambda m: m[1]+target, s)
        assert count==4; p.write_text(s)
        (folder/'profile_columns.list').write_text('! Output-only compact structure and mass audit.\n'+'\n'.join(cols)+'\n')
        inputs = {str(p.relative_to(folder)): sha(p) for p in sorted(folder.rglob('*')) if p.is_file()}
        dest = OUT/'runs'/name; dest.mkdir(parents=True)
        for rel in inputs:
            if rel in ('binary', 'photos1/18000'): continue
            shutil.copy2(folder/rel, dest/('restart-input.txt' if rel=='.restart' else rel))
        save(name+'-inputs.json', dict(classification='Counterexample candidate',
             run_directory=str(folder), sha256=inputs,
             plan_sha256=sha(OUT/'preregistered-plan.json'),
             original_runtime_manifest_sha256=sha(OLD/'environment.json')))
    print('Preregistered fixed age window', age_start, age_end, 'years;', CASES)


def run(name):
    assert name in CASES or name in PHASE_CASES
    folder=CACHE/name; logpath=folder/'execution.log'
    record=json.loads((OUT/(name+'-inputs.json')).read_text())
    plan_file=record.get('plan_file','preregistered-plan.json')
    assert plan_file in ('preregistered-plan.json','phase-plan.json','phase-restart-repair.json')
    assert sha(OUT/plan_file)==record['plan_sha256']
    assert not logpath.exists()
    for rel, digest in record['sha256'].items(): assert sha(folder/rel)==digest, rel
    env=os.environ.copy(); env.update(MESA_DIR=str(RUNTIME/'mesa-r7624'),
        LD_LIBRARY_PATH=str(RUNTIME/'mesasdk/lib')+':'+str(RUNTIME/'mesasdk/lib64'),
        OMP_NUM_THREADS='4')
    start=time.monotonic()
    with logpath.open('w') as log:
        try:
            result=subprocess.run(['./binary'], cwd=folder, env=env, stdout=log,
                stderr=subprocess.STDOUT, timeout=7200)
            code=result.returncode
        except subprocess.TimeoutExpired:
            code='timeout'
    shutil.copy2(logpath, OUT/(name+'-execution.log'))
    save(name+'-execution.json', dict(classification='Proven', returncode=code,
        elapsed_s=time.monotonic()-start, environment={k:env[k] for k in ['MESA_DIR','LD_LIBRARY_PATH','OMP_NUM_THREADS']}))
    print(name, code, time.monotonic()-start, flush=True)


def prepare_phase():
    import struct
    assert not (OUT/'phase-plan.json').exists()
    primary=json.loads((OUT/'preregistered-plan.json').read_text())
    end=primary['age_window_yr'][1]; parents={}
    for name in ['fast','time']:
        record=json.loads((OUT/(name+'-execution.json')).read_text())
        assert record['returncode']==0
        assert 'termination code: max_age' in (CACHE/name/'execution.log').read_text()
        _,d=mesa(CACHE/name/'LOGS1/history.data')
        selected=select(d,primary); assert selected['matching_rows']==0
        model=int(d['model_number'][-1]); photo=f'x{model%1000:03d}'
        path=CACHE/name/'photos1'/photo
        with path.open('rb') as f:
            for _ in range(3):
                marker=f.read(4); size=struct.unpack('<I',marker)[0]
                payload=f.read(size); assert f.read(4)==marker
        species,reactions,photo_model=struct.unpack_from('<iii',payload,60)
        assert photo_model==model and species==22 and reactions==72
        parents[name]=dict(model=model,photo=photo,photo_sha256=sha(path),
            history_raw_sha256=sha(CACHE/name/'LOGS1/history.data'),primary_matching_rows=0)
    save('phase-plan.json',dict(classification='Conjectural',before_task_checkpoint='8b95a41',
        scope='Post-primary follow-up: test existence at the next cooling crossing; preserve the failed fixed-window verdict.',
        registered_after='fast and time primary runs completed with zero matching rows',parents=parents,
        age_window_yr=[end,end+1e7],Teff_stop_K=15500,matching_rules=primary,
        success_rule='Actual mass/Teff/logg-passing rows after the original end age in both continuations. This tests age-shift compatibility only, not continuum accuracy.',
        guards=dict(max_model_number=40000,timeout_s=7200),
        original_verdict_may_be_changed=False))
    for parent,item in parents.items():
        name=parent+'_phase'; source=CACHE/parent; folder=CACHE/name
        folder.mkdir(exist_ok=False); (folder/'photos1').mkdir()
        for p in source.iterdir():
            if p.is_file() and (p.name.startswith('inlist') or p.suffix in ('.list','.net') or p.name=='binary'):
                shutil.copy2(p,folder/p.name)
        shutil.copy2(source/'photos1'/item['photo'],folder/'photos1'/item['photo'])
        (folder/'.restart').write_text(item['photo']+'\n')
        p=folder/'inlist1'; s=p.read_text()
        assert s.count('relax_mass = .true.')==1 and 'Teff_lower_limit' not in s
        s=s.replace('relax_mass = .true.','relax_mass = .false.')
        s,count=re.subn(r'max_age\s*=\s*\S+',f'max_age = {end+1e7:.17e}',s); assert count==1
        s=s.replace('&controls','&controls\n      Teff_lower_limit = 15500',1); p.write_text(s)
        dest=OUT/'runs'/name; dest.mkdir(parents=True)
        hashes={str(p.relative_to(folder)):sha(p) for p in sorted(folder.rglob('*')) if p.is_file()}
        for rel in hashes:
            if rel=='binary' or rel.startswith('photos1/'): continue
            shutil.copy2(folder/rel,dest/('restart-input.txt' if rel=='.restart' else rel))
        gzcopy(folder/'photos1'/item['photo'],dest/'restart-photo.gz')
        save(name+'-inputs.json',dict(classification='Conjectural',run_directory=str(folder),
            parent=parent,sha256=hashes,plan_file='phase-plan.json',plan_sha256=sha(OUT/'phase-plan.json')))
    print('Preregistered two phase continuations from exact final photos')


def prepare_phase_replay():
    import struct
    assert not (OUT/'phase-restart-repair.json').exists()
    plan=json.loads((OUT/'phase-plan.json').read_text()); repairs={}
    for parent,item in plan['parents'].items():
        failed=CACHE/(parent+'_phase')
        assert 'termination code: dt_is_zero' in (failed/'execution.log').read_text()
        model=(item['model']-1)//50*50; photo=f'x{model%1000:03d}'
        path=CACHE/parent/'photos1'/photo
        with path.open('rb') as f:
            for _ in range(3):
                marker=f.read(4); n=struct.unpack('<I',marker)[0]
                payload=f.read(n); assert f.read(4)==marker
        assert struct.unpack_from('<iii',payload,60)==(22,72,model)
        repairs[parent]=dict(replay_model=model,photo=photo,photo_sha256=sha(path))
    save('phase-restart-repair.json',dict(classification='Proven',original_plan_sha256=sha(OUT/'phase-plan.json'),
        reason='Exact max_age termination saves dt_next=0; both direct continuations stopped without taking a step.',
        correction='Use the last regular pre-termination photo, retain all physical controls, and compare shared replay history. Do not edit binary photo fields.',
        parents=repairs,age_window_yr=plan['age_window_yr'],Teff_stop_K=plan['Teff_stop_K'],
        physical_selection_or_success_rule_changed=False))
    for parent,item in repairs.items():
        source=CACHE/(parent+'_phase'); name=parent+'_phase_replay'; folder=CACHE/name
        folder.mkdir(exist_ok=False); (folder/'photos1').mkdir()
        for p in source.iterdir():
            if p.is_file() and (p.name.startswith('inlist') or p.suffix in ('.list','.net') or p.name=='binary'):
                shutil.copy2(p,folder/p.name)
        shutil.copy2(CACHE/parent/'photos1'/item['photo'],folder/'photos1'/item['photo'])
        (folder/'.restart').write_text(item['photo']+'\n')
        dest=OUT/'runs'/name; dest.mkdir(parents=True)
        hashes={str(p.relative_to(folder)):sha(p) for p in sorted(folder.rglob('*')) if p.is_file()}
        for rel in hashes:
            if rel=='binary' or rel.startswith('photos1/'): continue
            shutil.copy2(folder/rel,dest/('restart-input.txt' if rel=='.restart' else rel))
        gzcopy(folder/'photos1'/item['photo'],dest/'restart-photo.gz')
        save(name+'-inputs.json',dict(classification='Conjectural',run_directory=str(folder),
            parent=parent,sha256=hashes,plan_file='phase-restart-repair.json',plan_sha256=sha(OUT/'phase-restart-repair.json')))
    print('Prepared regular-photo replay; failed exact-endpoint attempts retained')


def control():
    _, a=mesa(OLD/'runs/mass18000/history.data.gz')
    _, b=mesa(CACHE/'control/LOGS1/history.data')
    assert np.array_equal(b['model_number'], np.arange(18001,18011))
    errors={k:float(np.max(abs(a[k][:10]-v))) for k,v in b.items()}
    # Wall-clock runtime is an execution measurement, not a stellar quantity.
    assert all(v==0 for k,v in errors.items() if k!='runtime_minutes')
    idx=np.loadtxt(CACHE/'control/LOGS1/profiles.index', skiprows=1, dtype=int, ndmin=2)
    old_idx=np.loadtxt(RUNTIME/'mass18000/LOGS1/profiles.index', skiprows=1, dtype=int, ndmin=2)
    shared={}
    for model in (18000,18010):
        i=idx[idx[:,0]==model]; j=old_idx[old_idx[:,0]==model]
        assert len(i)==len(j)==1
        _, x=mesa(CACHE/'control/LOGS1'/f'profile{i[0,2]}.data')
        _, y=mesa(RUNTIME/'mass18000/LOGS1'/f'profile{j[0,2]}.data')
        shared[str(model)]=sorted(x.keys() & y.keys())
        assert all(np.array_equal(x[k], y[k]) for k in shared[str(model)])
    assert 'termination code: max_model_number' in (CACHE/'control/execution.log').read_text()
    save('control-audit.json', dict(classification='Proven', all_physical_history_columns_identical=True,
        excluded_execution_counter='runtime_minutes', all_history_max_absolute_differences=errors,
        history_rows=10, history_columns=len(errors), profile_shared_columns=shared,
        shared_structure_identical=True))
    print('PASS: output-only control reproduces all physical history and shared structure columns')


def select(d, plan):
    units=json.loads((OUT/'source-bindings.json').read_text())['constants']
    mass=d['star_mass']*units['G_SI']*units['Msun_kg']/1.3271244e20
    teff=10**d['log_Teff']; logg=d['log_g']; age=d['star_age']
    window=(age>=plan['age_window_yr'][0]) & (age<=plan['age_window_yr'][1])
    mass_ok=window & (abs(mass/plan['mass_GM_sun']-1)<=plan['mass_relative_tolerance'])
    passes=mass_ok & (abs(teff-15800)<=300) & (abs(logg-5.82)<=.15)
    score=((teff-15800)/100)**2+((logg-5.82)/.05)**2
    eligible=np.flatnonzero(passes if passes.any() else mass_ok)
    assert len(eligible)>0
    def row(i):
        return dict(model=int(d['model_number'][i]), age_yr=float(age[i]),
            Teff_K=float(teff[i]), logg=float(logg[i]), mass_GM_sun=float(mass[i]),
            cells=int(d['num_zones'][i]), radius_source=float(10**d['log_R'][i]),
            H_source_mass=float(d['total_mass_h1'][i]), period_days=float(d['period_days'][i]),
            score=float(score[i]), candidate=bool(passes[i]))
    best=row(int(eligible[np.argmin(score[eligible])]))
    crossings=[]
    for i in np.flatnonzero((teff[:-1]-15800)*(teff[1:]-15800)<0):
        f=(15800-teff[i])/(teff[i+1]-teff[i])
        t=float(age[i]+f*(age[i+1]-age[i]))
        if not plan['age_window_yr'][0]<=t<=plan['age_window_yr'][1]: continue
        crossings.append(dict(age_yr=t, logg=float(logg[i]+f*(logg[i+1]-logg[i])),
            direction='heating' if teff[i+1]>teff[i] else 'cooling',
            bracketing_models=[int(d['model_number'][i]),int(d['model_number'][i+1])]))
    return dict(classification='Counterexample candidate', rows=len(age),
        first_age_yr=float(age[0]), last_age_yr=float(age[-1]),
        matching_rows=int(passes.sum()), best=best,
        candidates=[row(int(i)) for i in np.flatnonzero(passes)], crossings=crossings)


def analyze():
    import contextlib
    import thermal_restart
    plan=json.loads((OUT/'preregistered-plan.json').read_text())
    assert json.loads((OUT/'control-audit.json').read_text())['shared_structure_identical']
    _, baseline=mesa(OLD/'runs/mass18000/history.data.gz')
    results={'baseline':select(baseline, plan)}; histories={'baseline':baseline}
    for name in CASES:
        folder=CACHE/name; dest=OUT/'runs'/name
        record=json.loads((OUT/(name+'-execution.json')).read_text())
        assert record['returncode']==0
        gzcopy(folder/'LOGS1/history.data', dest/'history.data.gz')
        shutil.copy2(folder/'LOGS1/profiles.index', dest/'profiles.index')
        if name=='control': continue
        log=(folder/'execution.log').read_text()
        assert 'finished doing relax mass' in log
        _, hist=mesa(dest/'history.data.gz'); histories[name]=hist
        result=select(hist, plan)
        result['normal_end_age_termination']='termination code: max_age' in log
        result['end_age_reached']=bool(hist['star_age'][-1]>=plan['age_window_yr'][1]-.01)
        index=np.loadtxt(dest/'profiles.index', skiprows=1, dtype=int, ndmin=2)
        at=index[index[:,0]==result['best']['model']]; assert len(at)==1
        path=folder/'LOGS1'/f'profile{at[0,2]}.data'
        gzcopy(path, dest/'selected.data.gz')
        result['selected_profile_raw_sha256']=sha(path)
        result['selected_profile_gzip_sha256']=sha(dest/'selected.data.gz')
        result['elapsed_s']=record['elapsed_s']; results[name]=result
        old_out=thermal_restart.OUT
        try:
            thermal_restart.OUT=OUT
            with (OUT/(name+'-response.log')).open('w') as log, contextlib.redirect_stdout(log):
                thermal_restart.response_profile(path,result['best']['logg'],name,result['best']['candidate'])
        finally: thermal_restart.OUT=old_out
        save(name+'-mass-diagnostic.json', mass_diagnostic(dest/'selected.data.gz'))
    low=max(d['star_age'][0] for d in histories.values())
    high=min(plan['age_window_yr'][1], *(d['star_age'][-1] for d in histories.values()))
    assert high>low
    grid=np.linspace(low,high,501)
    comparison={}
    for name,d in histories.items():
        assert np.all(np.diff(d['star_age'])>0), name
        if name=='baseline': continue
        fields={}
        for key in ['log_Teff','log_R','log_g','total_mass_h1','period_days']:
            delta=np.interp(grid,d['star_age'],d[key])-np.interp(grid,baseline['star_age'],baseline[key])
            fields[key]=dict(max_absolute_difference=float(np.max(abs(delta))),
                rms_difference=float(np.sqrt(np.mean(delta**2))))
        comparison[name]=fields
    save('selection.json',results)
    save('common-age-comparison.json',dict(classification='Proven',diagnostic_interpolation_only=True,
        grid_rows=501, common_age_window_yr=[low,high], differences=comparison))
    print(json.dumps({k:{x:v[x] for x in ['rows','matching_rows','best','crossings']} for k,v in results.items()},indent=2))


def analyze_phase():
    import contextlib,thermal_restart
    plan=json.loads((OUT/'phase-plan.json').read_text())
    repairs=json.loads((OUT/'phase-restart-repair.json').read_text())['parents']
    rules=dict(plan['matching_rules'],age_window_yr=plan['age_window_yr'])
    results={}
    for parent,repair in repairs.items():
        name=parent+'_phase_replay'; folder=CACHE/name; dest=OUT/'runs'/name
        assert json.loads((OUT/(name+'-execution.json')).read_text())['returncode']==0
        log=(folder/'execution.log').read_text()
        gzcopy(folder/'LOGS1/history.data',dest/'history.data.gz')
        shutil.copy2(folder/'LOGS1/profiles.index',dest/'profiles.index')
        _,d=mesa(dest/'history.data.gz'); _,original=mesa(CACHE/parent/'LOGS1/history.data')
        r=select(d,rules)
        r['normal_temperature_termination']='termination code: Teff_lower_limit' in log
        r['original_fixed_window_verdict_changed']=False
        # The final primary step was clipped at max_age and is not a replay target.
        models=d['model_number']; take=models<original['model_number'][-1]
        source=np.searchsorted(original['model_number'],models[take])
        keys=['star_age','star_mass','log_R','log_Teff','log_g','log_center_T',
              'log_center_Rho','log_center_P','total_mass_h1','period_days','log_dt']
        errors={k:float(np.max(abs(d[k][take]-original[k][source]))) for k in keys}
        r['replay_comparison']=dict(rows=int(take.sum()),max_absolute_differences=errors,
            physical_history_identical=all(v==0 for v in errors.values()))
        assert r['replay_comparison']['physical_history_identical']
        index=np.loadtxt(dest/'profiles.index',skiprows=1,dtype=int,ndmin=2)
        old_index=np.loadtxt(CACHE/parent/'LOGS1/profiles.index',skiprows=1,dtype=int,ndmin=2)
        at=index[index[:,0]==repair['replay_model']]; before=old_index[old_index[:,0]==repair['replay_model']]
        assert len(at)==len(before)==1
        start=folder/'LOGS1'/f'profile{at[0,2]}.data'
        old_start=CACHE/parent/'LOGS1'/f'profile{before[0,2]}.data'
        _,x=mesa(start); _,y=mesa(old_start)
        assert x.keys()==y.keys() and all(np.array_equal(x[k],y[k]) for k in x)
        r['restored_structure_columns_identical']=len(x)
        gzcopy(start,dest/'restored.data.gz'); gzcopy(old_start,dest/'reference.data.gz')
        at=index[index[:,0]==r['best']['model']]; assert len(at)==1
        path=folder/'LOGS1'/f'profile{at[0,2]}.data'
        gzcopy(path,dest/'selected.data.gz'); r['selected_profile_raw_sha256']=sha(path)
        old_out=thermal_restart.OUT
        try:
            thermal_restart.OUT=OUT
            with (OUT/(name+'-response.log')).open('w') as stream,contextlib.redirect_stdout(stream):
                thermal_restart.response_profile(path,r['best']['logg'],name,r['best']['candidate'])
        finally: thermal_restart.OUT=old_out
        save(name+'-mass-diagnostic.json',mass_diagnostic(dest/'selected.data.gz'))
        results[parent]=r
    save('phase-selection.json',results)
    print(json.dumps({k:{x:v[x] for x in ['rows','matching_rows','best','crossings','replay_comparison']} for k,v in results.items()},indent=2))


def shell_energies(r, m, omega):
    """Exact Newtonian W and rotational energy of constant-density spherical shells."""
    G=6.67428e-11
    a,b=r[:-1],r[1:]; dm=np.diff(m)
    # Factored differences avoid cancellation in thin shells.
    d2=(b-a)*(b+a); d3=(b-a)*(b*b+a*b+a*a)
    d5=(b-a)*(b**4+b**3*a+b*b*a*a+b*a**3+a**4)
    B=dm/d3
    width=b-a
    self_term=width**2*(1.5*a**3+2*a*a*width+a*width**2+width**3/5)
    w=-3*G*np.sum(B*(m[:-1]*d2/2+B*self_term))
    rotation=np.sum(omega**2*B*d5/5)
    return float(w),float(rotation)


def mass_diagnostic(path):
    h,d=mesa(path)
    units=json.loads((OUT/'source-bindings.json').read_text())['constants']
    r=np.r_[0,d['radius'][::-1]*units['Rsun_m']]
    m=np.r_[0,d['mass'][::-1]*units['Msun_kg']]; dm=np.diff(m)
    assert np.all(np.diff(r)>0) and np.all(dm>0)
    native_dm=d['dm'][::-1]*.001
    assert np.max(abs(dm-native_dm))/m[-1]<1e-14
    assert abs(np.sum(native_dm)/m[-1]-1)<1e-13
    assert np.array_equal(d['m_grav'],d['mass'])
    assert np.all(d['m_grav_div_m_baryonic']==1)
    C=299792458.; scale=m[-1]*C*C
    u=float(np.dot(native_dm,d['energy'][::-1]*1e-4))
    w,t=shell_energies(r,m,d['omega'][::-1])
    correction=float(np.dot(native_dm,d['mass_correction_factor'][::-1])/m[-1]-1)
    lines=(OUT/'sources/data/chem_data/isotopes.data').read_text().splitlines()[1:]
    weights={}
    for line in lines[::4]:
        words=line.split()
        if words: weights[words[0]]=float(words[1])/(int(words[2])+int(words[3]))
    isos=json.loads((OLD/'network-diagnosis.json').read_text())['profile_isotopes']
    reconstructed=sum(d[k]*weights[k] for k in isos)
    weight_error=float(np.max(abs(reconstructed-d['mass_correction_factor'])))
    assert weight_error<3e-15
    return dict(classification='Proven',scope='fixed exported Newtonian shells; diagnostic energy scales only',
        model=int(h['model_number']), nuclear_weight_mass_fraction_shift=correction,
        independently_reconstructed_weight_max_error=weight_error,
        internal_energy_over_Mc2=u/scale, Newtonian_binding_over_Mc2=w/scale,
        spherical_rotation_over_Mc2=t/scale,
        surface_GM_Rc2=float(units['G_SI']*m[-1]/r[-1]/C**2),
        maximum_shell_GM_rc2=float(np.max(units['G_SI']*m[1:]/r[1:]/C**2)),
        max_pressure_over_baryon_rest_energy=float(np.max(d['pressure']/d['rho']/(C*100)**2)),
        internal_plus_binding_plus_rotation_over_Mc2=(u+w+t)/scale,
        physical_ADM_mass_correction_computed=False,
        warning='Do not add these diagnostics and call the result ADM mass: absolute EOS/rest-energy convention, proper volume, structural response, and rotating metric remain to be completed.')


def symbolic():
    import sympy as s
    r,a,b,B,m0,G,M,R,O,rho,C,delta=s.symbols('r a b B m0 G M R O rho C delta',positive=True)
    mass=m0+B*(r**3-a**3)
    energy=-3*G*B*((m0-B*a**3)*(b**2-a**2)/2+B*(b**5-a**5)/5)
    assert s.simplify(s.integrate(-G*mass/r*s.diff(mass,r),(r,a,b))-energy)==0
    width=b-a
    stable=-3*G*B*(m0*(b*b-a*a)/2+B*width**2*(s.Rational(3,2)*a**3+2*a*a*width+a*width**2+width**3/5))
    assert s.simplify(stable-energy)==0
    assert s.simplify(energy.subs({a:0,m0:0,B:M/R**3,b:R})+3*G*M*M/(5*R))==0
    rotation=O**2*B*(b**5-a**5)/5
    assert s.simplify(rotation.subs({a:0,B:M/R**3,b:R})-M*R**2*O**2/5)==0
    e=s.Function('e')(rho)
    assert s.simplify(rho**2*s.diff(e+delta,rho)-rho**2*s.diff(e,rho))==0
    assert s.simplify(rho*(C*C+e+delta)-rho*(C*C+e)-rho*delta)==0
    pressure,cx,ee=s.symbols('pressure cx ee',positive=True)
    epsilon=rho*(cx*C*C+ee); metric=1-2*G*M/(r*C*C)
    pprime=-G*(epsilon+pressure)*(M+4*s.pi*r**3*pressure/C**2)/(C*C*r*r*metric)
    baryon_prime=4*s.pi*r*r*rho/s.sqrt(metric)
    expected=-G*(cx+(ee+pressure/rho)/C**2)*(M+4*s.pi*r**3*pressure/C**2)/(4*s.pi*r**4*s.sqrt(metric))
    assert s.simplify(pprime/baryon_prime-expected)==0
    # Independent uniform-sphere numerical control of the integration code.
    rr=np.linspace(0,2e7,127); mm=3e29*(rr/rr[-1])**3
    w,t=shell_energies(rr,mm,np.full(126,1e-4))
    assert abs(w/(-3*6.67428e-11*mm[-1]**2/(5*rr[-1]))-1)<1e-12
    assert abs(t/(mm[-1]*rr[-1]**2*1e-8/5)-1)<1e-12
    save('symbolic-audit.json',dict(classification='Proven',shell_binding_antiderivative=True,
        uniform_sphere_binding_and_rotation_controls=True,
        Newtonian_pressure_invariant_under_specific_energy_constant=True,
        total_energy_source_changes_when_rest_convention_fixed=True,
        TOV_pressure_in_proper_baryon_mass_coordinate=True,
        implication='Newtonian pressure and its thermodynamic derivatives alone do not specify the relativistic energy source. This is not a proof that source-normalized EOS completion is impossible.'))
    print('PASS: shell antiderivatives, uniform sphere and absolute-energy completion boundary')


def profile_bound():
    """Outward perturbation bound for two fixed, positive, spherical potentials."""
    import mpmath as mp
    mp.iv.dps=65
    def potential(r,m):
        rr=[mp.iv.mpf(float(v)) for v in r]
        mm=[mp.iv.mpf(float(v)) for v in m]
        return [12*(y-x)/(b**3-a**3) for a,b,x,y in zip(rr,rr[1:],mm,mm[1:])]
    base=np.load(OLD/'mass-shells.npz')
    r1,m1=base['r_m'],base['m_geom_m']; V1=potential(r1,m1)
    reference=json.loads((OLD/'mass-response.json').read_text())
    # Recompute eta intervals; recorded eta was a point value, not an enclosure.
    from thermal_wd import shell_response
    chi1,eta1=shell_response(r1,m1,ctx=mp.iv)
    Q1=4*mp.iv.mpf(float(m1[-1]))
    results={}
    for name in [k for k in CASES if k!='control']+['fast_phase_replay','time_phase_replay']:
        data=np.load(OUT/(name+'-shells.npz')); r2,m2=data['r_m'],data['m_geom_m']
        chi2,eta2=shell_response(r2,m2,ctx=mp.iv)
        assert eta1.b<1 and eta2.b<1
        V2=potential(r2,m2); Q2=4*mp.iv.mpf(float(m2[-1]))
        grid=np.unique(np.r_[r1,r2]); delta_eta=mp.iv.mpf(0); delta_Q=mp.iv.mpf(0)
        for left,right in zip(grid,grid[1:]):
            i=np.searchsorted(r1,left,side='right')-1
            j=np.searchsorted(r2,left,side='right')-1
            diff=abs((V1[i] if i<len(V1) else 0)-(V2[j] if j<len(V2) else 0))
            a,b=mp.iv.mpf(float(left)),mp.iv.mpf(float(right))
            delta_eta+=diff*(b*b-a*a)/2
            delta_Q+=diff*(b**3-a**3)/3
        bound1=abs(Q1-Q2)+delta_Q*eta1/(1-eta1)+Q2*delta_eta/((1-eta1)*(1-eta2))
        bound2=abs(Q1-Q2)+delta_Q*eta2/(1-eta2)+Q1*delta_eta/((1-eta1)*(1-eta2))
        bound=bound1 if bound1.b<=bound2.b else bound2
        observed=abs(chi1-chi2)
        assert observed.b<=bound.a
        results[name]=dict(classification='Proven',
            delta_eta_interval=str(delta_eta), delta_Q_interval_m=str(delta_Q),
            absolute_response_bound_interval_m=str(bound),
            absolute_response_bound_upper_m=float(np.nextafter(float(bound.b),np.inf)),
            observed_response_difference_interval_m=str(observed),
            reference_profile_sha256=reference['raw_profile_sha256'],
            candidate_shell_sha256=sha(OUT/(name+'-shells.npz')),
            certifies_MESA_discretization_error=False)
    save('profile-perturbation-bounds.json',results)
    print('PASS:',len(results),'outward frozen-profile perturbation bounds')


def sources():
    import zipfile
    wanted=['chem/public/chem_lib.f90','chem/private/chem_isos_io.f90',
        'data/chem_data/isotopes.data','eos/public/eos_def.f',
        'star/private/star_utils.f90','star/private/hydro_vars.f90',
        'star/private/profile_getval.f90','star/defaults/controls.defaults',
        'star/defaults/profile_columns.list','star/private/relax.f90',
        'star/private/evolve.f90','star/job/run_star_support.f90']
    bindings={}
    with zipfile.ZipFile(RUNTIME/'mesa-r7624.zip') as archive:
        for rel in wanted:
            src=RUNTIME/'mesa-r7624'/rel
            member='chem/data/isotopes.data' if rel=='data/chem_data/isotopes.data' else rel
            raw=src.read_bytes(); assert raw==archive.read('mesa-r7624/'+member)
            dest=OUT/'sources'/rel; dest.parent.mkdir(parents=True,exist_ok=True)
            dest.write_bytes(raw); bindings[rel]=sha(dest)
    save('mass-definition-sources.json',dict(classification='Imported from prior work',
        sources_equal_verified_release_archive=True,sha256=bindings,
        MESA_release='https://zenodo.org/records/2630796',
        TOV_reference='https://doi.org/10.1093/mnras/staa1493',
        TOV_equations='Pretel and da Silva 2020, equations 4-9; inspected 2026-09-10',
        rotating_mass_reference='https://arxiv.org/abs/1905.03784',
        rotating_equations='Gao et al. 2019, equations 1-3; not importing empirical NS fits into the WD calculation'))
    print('PASS: exact release source for nuclear-weight mass, EOS energy and GR-factor audit')


def runtime_audit():
    record=json.loads((OLD/'environment.json').read_text())
    for item in record['archives']:
        assert sha(RUNTIME/item['name'])==item['sha256'],item['name']
    for item in record['dynamic_libraries']:
        assert sha(item['path'])==item['sha256'],item['path']
    for path,expected in record['data_sha256'].items():
        assert sha(RUNTIME/'mesa-r7624'/path)==expected,path
    save('runtime-audit.json',dict(classification='Proven',
        verified_archives=len(record['archives']),verified_libraries=len(record['dynamic_libraries']),
        verified_static_data_files=len(record['data_sha256']),
        matches_frozen_request19_environment=True,excluded_native_generated_caches=True))
    print('PASS: active runtime archives, libraries and',len(record['data_sha256']),'static data files')


def git_bytes():
    """Check actual committed bytes, which a normalization-aware git diff can hide."""
    entries=subprocess.check_output(['git','ls-tree','-r','-z','HEAD'],cwd=ROOT).decode().split('\0')
    records=[]; count=0; skipped=[]
    with subprocess.Popen(['git','cat-file','--batch'],cwd=ROOT,stdin=subprocess.PIPE,stdout=subprocess.PIPE) as proc:
        for entry in entries:
            if not entry: continue
            info,path=entry.split('\t',1); mode,kind,oid=info.split()
            if mode not in ('100644','100755'):
                skipped.append(dict(path=path,mode=mode)); continue
            proc.stdin.write((oid+'\n').encode()); proc.stdin.flush()
            header=proc.stdout.readline().split()
            assert len(header)==3 and header[1]==b'blob',header
            raw=proc.stdout.read(int(header[2])); assert proc.stdout.read(1)==b'\n'
            actual=(ROOT/path).read_bytes(); count+=1
            if raw==actual: continue
            newline_only=raw.replace(b'\r\n',b'\n')==actual.replace(b'\r\n',b'\n')
            assert newline_only,path
            records.append(dict(path=path,git_blob_sha256=hashlib.sha256(raw).hexdigest(),
                audited_worktree_sha256=sha(ROOT/path),newline_only=True))
        proc.stdin.close(); assert proc.wait()==0
    save('git-byte-audit-before.json',dict(classification='Proven',
        before_repair_checkpoint='be010c7',tracked_regular_files=count,skipped_nonregular_entries=skipped,mismatches=records,
        conclusion='Local hash checks alone did not certify Git-exported bytes. All detected differences are newline conversion.'))
    print(count,'tracked files;',len(records),'newline-only Git/worktree differences')


def retain_control():
    for label,folder in [('control',CACHE/'control'),('reference',RUNTIME/'mass18000')]:
        index=np.loadtxt(folder/'LOGS1/profiles.index',skiprows=1,dtype=int,ndmin=2)
        for model in (18000,18010):
            at=index[index[:,0]==model]; assert len(at)==1
            gzcopy(folder/'LOGS1'/f'profile{at[0,2]}.data',OUT/f'{label}-{model}.data.gz')


def plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    fig,axes=plt.subplots(1,2,figsize=(11,4.4),constrained_layout=True)
    selections=json.loads((OUT/'selection.json').read_text())
    for name,color in zip(['baseline','fast','slow','mesh','time'],['black','#0072B2','#D55E00','#009E73','#CC79A7']):
        path=OLD/'runs/mass18000/history.data.gz' if name=='baseline' else OUT/'runs'/name/'history.data.gz'
        _,d=mesa(path)
        axes[0].plot(10**d['log_Teff'],d['log_g'],color=color,lw=1,label=name)
        best=selections[name]['best']
        axes[0].scatter(best['Teff_K'],best['logg'],color=color,s=27,marker='o' if best['candidate'] else 'x',zorder=3)
        path=OLD/'mass-shells.npz' if name=='baseline' else OUT/(name+'-shells.npz')
        d=np.load(path)
        axes[1].plot(d['r_m'][1:]/1000,d['m_geom_m'][1:]/d['m_geom_m'][-1],color=color,lw=1,label=name)
    axes[0].add_patch(Rectangle((15500,5.67),600,.3,facecolor='#FFD966',alpha=.3,edgecolor='#8A6D00'))
    axes[0].set(xlim=(14000,21000),ylim=(5.3,6.4),xlabel='Effective temperature (K)',ylabel='log10 g (cgs)',title='Actual evolution near the optical target')
    axes[0].legend(frameon=False,fontsize=8)
    axes[1].set(xscale='log',xlim=(100,120000),ylim=(0,1.02),xlabel='Radius (km)',ylabel='Enclosed model mass / total model mass',title='Selected fixed Newtonian profiles')
    for ax in axes: ax.grid(alpha=.18)
    fig.savefig(OUT/'thermal-comparison.png',dpi=170)
    fig.savefig(OUT/'thermal-comparison.pdf')
    plt.close(fig)
    save('plot-provenance.json',dict(classification='Proven',rendering_only=True,
        versions=dict(numpy=np.__version__,matplotlib=matplotlib.__version__),
        module_file_sha256={str(m.__file__):sha(m.__file__) for m in [np,matplotlib]},
        command='PYTHONPATH= python3 -s verification/thermal_robustness.py plot',
        reason='System matplotlib uses NumPy 1.x ABI; scientific calculations retain their separately bound NumPy 2 environment.'))


def check():
    import mpmath as mp
    from thermal_wd import shell_response
    mp.mp.dps=80
    _,a=mesa(OLD/'runs/mass18000/history.data.gz')
    _,b=mesa(OUT/'runs/control/history.data.gz')
    assert np.array_equal(b['model_number'],np.arange(18001,18011))
    assert all(np.array_equal(a[k][:10],v) for k,v in b.items() if k!='runtime_minutes')
    for model in (18000,18010):
        _,x=mesa(OUT/f'control-{model}.data.gz'); _,y=mesa(OUT/f'reference-{model}.data.gz')
        assert all(np.array_equal(x[k],y[k]) for k in x.keys() & y.keys())
    selections=json.loads((OUT/'selection.json').read_text())
    plan=json.loads((OUT/'preregistered-plan.json').read_text())
    bound_records=json.loads((OUT/'profile-perturbation-bounds.json').read_text())
    summary={}
    for name in CASES:
        if name=='control': continue
        _,d=mesa(OUT/'runs'/name/'history.data.gz')
        r=selections[name]
        assert r['rows']==len(d['star_age'])
        c=json.loads((OUT/'source-bindings.json').read_text())['constants']
        eligible=[]
        for i in range(len(d['star_age'])):
            mass=d['star_mass'][i]*c['G_SI']*c['Msun_kg']/1.3271244e20
            if (plan['age_window_yr'][0]<=d['star_age'][i]<=plan['age_window_yr'][1]
                and abs(mass/.197536385307-1)<=1e-8
                and 15500<=10**d['log_Teff'][i]<=16100 and 5.67<=d['log_g'][i]<=5.97):
                eligible.append(int(d['model_number'][i]))
        assert eligible==[row['model'] for row in r['candidates']]
        assert len(eligible)==r['matching_rows']
        p=OUT/'runs'/name/'selected.data.gz'; h,structure=mesa(p)
        assert hashlib.sha256(gzip.decompress(p.read_bytes())).hexdigest()==r['selected_profile_raw_sha256']
        assert int(h['model_number'])==r['best']['model']
        assert abs(float(h['Teff'])/r['best']['Teff_K']-1)<1e-14
        assert abs(float(h['photosphere_r'])/r['best']['radius_source']-1)<1e-14
        assert int(h['num_zones'])==r['best']['cells']
        shell=np.load(OUT/(name+'-shells.npz'))
        actual,eta=shell_response(shell['r_m'],shell['m_geom_m'])
        response=json.loads((OUT/(name+'-response.json')).read_text())
        lo,hi=map(mp.mpf,response['susceptibility_interval_m'].strip('[]').split(','))
        assert lo<=actual<=hi and eta<1
        lo,hi=map(mp.mpf,response['relative_bound_interval'].strip('[]').split(','))
        assert mp.mpf(response['relative_bound_upper'])>=hi
        rb=bound_records[name]
        lo,hi=map(mp.mpf,rb['absolute_response_bound_interval_m'].strip('[]').split(','))
        assert mp.mpf(rb['absolute_response_bound_upper_m'])>=hi
        baseline_response=json.loads((OLD/'mass-response.json').read_text())
        assert abs(actual-mp.mpf(baseline_response['susceptibility_m']))<=lo
        diagnostic=mass_diagnostic(p)
        assert diagnostic==json.loads((OUT/(name+'-mass-diagnostic.json')).read_text())
        summary[name]=dict(rows=r['rows'],matching_rows=len(eligible),
            normal_end_age_termination=r['normal_end_age_termination'],end_age_reached=r['end_age_reached'])
    z=mp.sqrt(12*mp.mpf(.001))
    chi,_=shell_response(np.array([0.,1.]),np.array([0.,.001]))
    assert abs(chi-(mp.tan(z)/z-1))<mp.mpf('1e-70')
    # Opposite-radius uniform spheres independently exercise equal-mass response changes.
    chi2,_=shell_response(np.array([0.,2.]),np.array([0.,.001]))
    z2=mp.sqrt(12*mp.mpf(.001)/2)
    assert abs(chi2-2*(mp.tan(z2)/z2-1))<mp.mpf('1e-70')
    symbolic()
    save('audit.json',dict(classification='Proven',physical_control_history_and_structure_identical=True,
        candidate_rows_independently_selected=True,profile_history_and_raw_sha_checked=True,
        scalar_enclosures_recomputed_at_80_digits=True,independent_uniform_sphere_controls=True,
        outward_scalar_frequency_and_profile_bounds_checked=True,
        nuclear_weights_independently_reconstructed=True,results=summary))
    print('PASS: retained controls, independent row selection, structures, enclosures and nuclear weights')


def check_phase():
    import mpmath as mp
    from thermal_wd import shell_response
    mp.mp.dps=80
    results=json.loads((OUT/'phase-selection.json').read_text())
    plan=json.loads((OUT/'phase-plan.json').read_text())
    bounds=json.loads((OUT/'profile-perturbation-bounds.json').read_text())
    verified={}
    for parent,r in results.items():
        name=parent+'_phase_replay'; dest=OUT/'runs'/name
        assert r['normal_temperature_termination'] and r['replay_comparison']['physical_history_identical']
        _,a=mesa(dest/'restored.data.gz'); _,b=mesa(dest/'reference.data.gz')
        assert a.keys()==b.keys() and all(np.array_equal(a[k],b[k]) for k in a)
        _,d=mesa(dest/'history.data.gz'); _,original=mesa(OUT/'runs'/parent/'history.data.gz')
        take=d['model_number']<original['model_number'][-1]
        at=np.searchsorted(original['model_number'],d['model_number'][take])
        for k in r['replay_comparison']['max_absolute_differences']:
            assert np.array_equal(d[k][take],original[k][at])
        candidates=[]
        c=json.loads((OUT/'source-bindings.json').read_text())['constants']
        for i in range(len(d['star_age'])):
            mass=d['star_mass'][i]*c['G_SI']*c['Msun_kg']/1.3271244e20
            if (plan['age_window_yr'][0]<=d['star_age'][i]<=plan['age_window_yr'][1]
                and abs(mass/.197536385307-1)<=1e-8
                and 15500<=10**d['log_Teff'][i]<=16100 and 5.67<=d['log_g'][i]<=5.97):
                candidates.append(int(d['model_number'][i]))
        assert candidates==[v['model'] for v in r['candidates']]
        h,_=mesa(dest/'selected.data.gz')
        assert int(h['model_number'])==r['best']['model']
        assert abs(float(h['Teff'])/r['best']['Teff_K']-1)<1e-14
        assert hashlib.sha256(gzip.decompress((dest/'selected.data.gz').read_bytes())).hexdigest()==r['selected_profile_raw_sha256']
        shells=np.load(OUT/(name+'-shells.npz')); chi,_=shell_response(shells['r_m'],shells['m_geom_m'])
        response=json.loads((OUT/(name+'-response.json')).read_text())
        lo,hi=map(mp.mpf,response['susceptibility_interval_m'].strip('[]').split(','))
        assert lo<=chi<=hi
        lo,hi=map(mp.mpf,response['relative_bound_interval'].strip('[]').split(','))
        assert mp.mpf(response['relative_bound_upper'])>=hi
        lo,hi=map(mp.mpf,bounds[name]['absolute_response_bound_interval_m'].strip('[]').split(','))
        assert mp.mpf(bounds[name]['absolute_response_bound_upper_m'])>=hi
        base=mp.mpf(json.loads((OLD/'mass-response.json').read_text())['susceptibility_m'])
        assert abs(chi-base)<=lo
        assert mass_diagnostic(dest/'selected.data.gz')==json.loads((OUT/(name+'-mass-diagnostic.json')).read_text())
        verified[parent]=dict(replayed_rows=int(take.sum()),identical_structure_columns=len(a),
            phase_candidate_rows=len(candidates))
    save('phase-audit.json',dict(classification='Proven',results=verified,
        fixed_window_failures_preserved=True,phase_candidate_selection_independently_checked=True,
        scalar_enclosures_and_response_difference_bounds_checked=True))
    print('PASS: independent phase candidates, exact replay, structures and scalar enclosures')


def report():
    results=json.loads((OUT/'selection.json').read_text())
    bounds=json.loads((OUT/'profile-perturbation-bounds.json').read_text())
    success=all(r['matching_rows']>0 and r['normal_end_age_termination'] and r['end_age_reached']
        for name,r in results.items() if name!='baseline')
    text='\n## 실행 결과\n\n분류: Counterexample candidate. 고정 나이 구간의 유한 변형 검사 '+('통과' if success else '미통과')+'. '
    text+='아래 값은 각 실행의 실제 선택 구조다. 통과 행이 없는 실행은 질량 조건 안의 최선 행을 표시하며 후보로 승격하지 않는다.\n\n'
    text+='| 실행 | 진화 행 | 조건 통과 행 | 모델 | 셀 수 | Teff (K) | logg | 나이 (년) |\n|---|---:|---:|---:|---:|---:|---:|---:|\n'
    for name,r in results.items():
        b=r['best']; text+=f"| {name} | {r['rows']} | {r['matching_rows']} | {b['model']} | {b['cells']} | {b['Teff_K']:.4f} | {b['logg']:.7f} | {b['age_yr']:.6f} |\n"
    text+='\n분류: Proven. 각 실행의 Teff=15800 K 교차는 `selection.json`에 가열·냉각 방향과 두 실제 모델 번호를 함께 보존했다. 공통 나이 501개 표본의 보간 차이는 `common-age-comparison.json`에 있다. 이 표본의 최대 차이를 전구간 엄밀 최대 오차로 해석하지 않는다.\n'
    text+='\n분류: Counterexample candidate. 다음 그림은 원래 고정 구간 실행만 표시한다. 음영은 사전 탐색 사각형이며 확률 등고선이 아니다. 원은 통과한 선택 구조, X는 통과하지 못한 최선 구조다.\n\n![열 진화와 선택 구조](../outputs/thermal-robustness20/thermal-comparison.png)\n'
    text+='\n분류: Proven. 각 선택 구조의 영 배경·고정 밀도 scalar 결과는 다음과 같다. 구조 차이 상계는 기존 기준 구조와의 차이이며, 시간 분해능 오차 인증이 아니다. 표는 표시용 반올림 값이며 엄밀한 바깥 방향 상계는 JSON 원문에 보존했다.\n\n'
    text+='| 실행 | chi (m) | 궤도 주파수 상대 차이 상계 | 기준 구조와 chi 차이 상계 (m) |\n|---|---:|---:|---:|\n'
    for name in CASES:
        if name=='control': continue
        r=json.loads((OUT/(name+'-response.json')).read_text())
        text+=f"| {name} | {float(r['susceptibility_m']):.9f} | {r['relative_bound_upper']:.9e} | {bounds[name]['absolute_response_bound_upper_m']:.9e} |\n"
    text+='\n분류: Proven. 선택 구조의 질량·에너지 진단은 다음과 같다. 모든 수는 무차원이며 마지막 두 에너지는 Mc²로 나눈 값이다. C_X는 배포 원자량과 22개 조성에서 독립 재계산했다.\n\n'
    text+='| 실행 | 가중 평균 C_X−1 | U/(Mc²) | W/(Mc²) |\n|---|---:|---:|---:|\n'
    for name in CASES:
        if name=='control': continue
        r=json.loads((OUT/(name+'-mass-diagnostic.json')).read_text())
        text+=f"| {name} | {r['nuclear_weight_mass_fraction_shift']:.9e} | {r['internal_energy_over_Mc2']:.9e} | {r['Newtonian_binding_over_Mc2']:.9e} |\n"
    text+='\n분류: Conjectural. 다음 완료 경계는 절대 에너지 기준을 결속한 열 EOS와 상대론적 구조·질량 재계산이다. 한 초기 모형의 유한 변형 시험은 전체 형성 이력의 강건성이나 연속 격자 극한을 인증하지 않는다. 실제 궤도에 맞는 형성 모형, 유체·metric·scalar 결합, 전구간 미분 오차, 전체 관측 비선형 추론도 미완료다. 이번 결과를 논문 제출 준비 완료로 판정하지 않는다.\n'
    text+='\n재검증은 저장된 증거만 사용하여 `thermal_robustness.py check`와 `thermal_robustness.py verify`를 실행한다. `check`는 감사 파일을 같은 내용으로 재생성하며 `verify`는 현재 및 역사적 해시를 읽기만 한다. 기존 원고 PDF·ZIP은 Request 12 동결본이며 이번 추가 연구를 아직 포함하지 않는다.\n'
    path=ROOT/'notes/REQUEST20_THERMAL_ROBUSTNESS_KO.md'
    before=path.read_text().split('\n## 실행 결과\n')[0]
    path.write_text(before+text)
    save('gates.json',dict(classification='Proven',theorem_progress=True,
        frozen_profile_perturbation_bound_proven=True,
        four_preregistered_variants_completed=all(r['normal_end_age_termination'] and r['end_age_reached'] for k,r in results.items() if k!='baseline'),
        all_four_variants_have_mass_and_optical_candidate=success,
        mass_definition_source_audit_completed=True,
        git_export_byte_preservation_repaired=True,
        exact_J0337_GR_mass_and_structure_matched=False,
        full_stellar_interval_certificate=False,continuum_mesh_time_convergence_certified=False,
        independent_formation_histories_validated=False,genuine_orbital_timescale_state_established=False,
        full_physical_dynamic_force_and_readout=False,full_span_variational_certificate=False,
        full_28_parameter_initialization_certificate=False,complete_nonlinear_observational_inference=False))


def report_phase():
    results=json.loads((OUT/'phase-selection.json').read_text())
    baseline=json.loads((OLD/'mass-selection.json').read_text())['best']['star_age_yr']
    text='\n## 고정 구간 실패 뒤 사전 등록한 후속 나이 검사\n\n'
    text+='분류: Proven. 빠른 제거율과 세밀한 시간 설정의 고정 구간 실패를 확인한 뒤 별도 `phase-plan.json`을 등록했다. 기준 질량·온도·표면중력은 그대로 두고, 원래 종료 나이 이후 다음 냉각에서 Teff<15500 K가 될 때까지 검사한다. 추가 1천만 년과 최대 모델 40000의 중단 한도를 두었다. 이 후속 결과는 원래 유한 변형 검사의 실패 판정을 바꿀 수 없다. 후속 작업 전 체크포인트는 `8b95a41`이다.\n\n'
    text+='분류: Proven. 정확한 종료 사진은 dt_next=0을 저장하여 두 직접 재시작이 진화 없이 종료했다. 실패 로그를 보존하고, 첫 후속 진화 결과를 얻기 전에 `phase-restart-repair.json`에 직전 정상 사진 사용을 명시했다. 바이너리 사진의 숫자를 직접 고치지 않았다. 기존 종료 나이에 맞춰 잘린 마지막 단계는 재생 대조에서 제외했다.\n\n'
    for parent,r in results.items():
        text+=f"분류: Proven. {parent}의 직전 사진 재생에서 구조 {r['restored_structure_columns_identical']}개 열과 기존 이력 {r['replay_comparison']['rows']}개 행의 질량·열 구조·나이·dt·주기를 동일하게 재현했다.\n\n"
    text+='분류: Counterexample candidate. 아래는 후속 나이 구간의 실제 후보다. 행 수에는 기준 구간을 재생한 단계도 포함되지만, 통과 행은 원래 종료 나이 뒤의 시점만 센다.\n\n'
    text+='| 후속 실행 | 진화 행 | 조건 통과 행 | Teff (K) | logg | 기준 후보보다 늦은 나이 (년) | 궤도 주파수 상대 scalar 상계 |\n|---|---:|---:|---:|---:|---:|---:|\n'
    for parent,r in results.items():
        b=r['best']; response=json.loads((OUT/(parent+'_phase_replay-response.json')).read_text())
        text+=f"| {parent} | {r['rows']} | {r['matching_rows']} | {b['Teff_K']:.4f} | {b['logg']:.7f} | {b['age_yr']-baseline:.3f} | {response['relative_bound_upper']:.9e} |\n"
    spans={name:r['candidates'][-1]['age_yr']-r['candidates'][0]['age_yr'] for name,r in results.items() if r['candidates']}
    original=json.loads((OLD/'mass-selection.json').read_text())['candidates']
    original_span=original[-1]['star_age_yr']-original[0]['star_age_yr']
    text+=f'\n분류: Proven. 기존 후보의 첫·마지막 통과 표본 사이 나이는 {original_span:.3f}년이다. '
    text+=', '.join(f'{name} 후속은 {span:.3f}년' for name,span in spans.items())+'이다. 이 표본 간격보다 진화 시점 이동이 훨씬 크다. 연속 통과 구간 길이의 엄밀한 측정이나 관측 발생 확률로 해석하지 않는다.\n'
    text+='\n분류: Counterexample candidate. 후속 냉각에서 후보가 다시 나타나는 경우는 진화 시점 이동과 양립한다. 이는 나이를 고정하지 않은 후보 존재에 대한 제한된 근거다. 같은 나이의 수치 수렴, 동일한 형성 이력, GR 질량 일치, 관측 posterior의 강건성을 증명하지 않는다. 후속 구조도 구각 응답 차이 상계와 질량 정의 감사에 포함했다.\n'
    bounds=json.loads((OUT/'profile-perturbation-bounds.json').read_text())
    largest=max(bounds[k]['absolute_response_bound_upper_m'] for k in ['slow','mesh','fast_phase_replay','time_phase_replay'])
    ceiling=float(np.ceil(largest*1e6)/1000)
    text+=f'\n분류: Proven. 질량·광학 기준을 통과한 네 변형 구조와 기존 기준 구조 사이 정적 scalar 응답 차이의 계산된 상계는 각각 {ceiling:.3f} mm 이하이다. 이 수치는 고정 구각 모형 사이의 차이에 대한 보장이며 실제 별의 수치·물리 오차 상한이 아니다.\n'
    path=ROOT/'notes/REQUEST20_THERMAL_ROBUSTNESS_KO.md'
    path.write_text(path.read_text().split('\n## 고정 구간 실패 뒤 사전 등록한 후속 나이 검사\n')[0]+text)
    gates=json.loads((OUT/'gates.json').read_text())
    gates['post_primary_age_followup_completed']=all(r['normal_temperature_termination'] for r in results.values())
    gates['both_failed_cases_recover_candidates_at_later_ages']=all(r['matching_rows']>0 for r in results.values())
    gates['original_fixed_window_verdict_preserved']=True
    save('gates.json',gates)


def maintain():
    results=json.loads((OUT/'selection.json').read_text())
    counts=', '.join(f"{name} {r['matching_rows']}개" for name,r in results.items() if name!='baseline')
    gates=json.loads((OUT/'gates.json').read_text())
    status='통과' if gates['all_four_variants_have_mass_and_optical_candidate'] else '미통과'
    additions={
        'model-definition':f'분류: Counterexample candidate. 같은 초기 사진의 제거율 두 가지 및 공간·시간 분해능 변형 네 가지를 고정 나이 구간에서 검사했다. 실제 질량·광학 조건 통과 행은 {counts}이며 유한 변형 판정은 {status}다. 독립 형성 모형이나 GR 질량 일치로 승격하지 않는다.',
        'observable-targets':f'분류: Counterexample candidate. Request 20의 실제 진화 행에 동일한 질량·Teff·logg 선택 규칙을 적용했다: {counts}. 모든 교차와 통과 행을 보존했다. 같은 물리적 나이의 표본 차이와 온도 교차를 구분했으며, 보간값을 실제 내부 구조나 관측 우도로 쓰지 않았다.',
        'adiabatic-limit':'분류: Proven. 두 고정 양의 구대칭 퍼텐셜의 eta_i<1 조건에서 정적 응답 차이의 해석 상계를 도출했다. 질량이 같으면 선도 Q 차이가 소거된다. 실제 변형 구조의 구각 합집합에서 65자리 구간 연산으로 상계를 계산했다. 이는 물리적 항성 오차 또는 전구간 미분 오차 보장이 아니다.',
        'nonadiabatic-regime':'분류: Proven. 네 변형 실행의 선택 구조에 기존 영 배경·고정 밀도·평탄 시공간 scalar 응답과 궤도 주파수 상계를 다시 적용했다. 별의 열 진화 시점 민감도와 독립 scalar 모형의 정적 근사는 다른 결과다. 유체·회전·metric 결합의 궤도 시간척도 상태를 배제하지 않는다.',
        'failure-ledger-dynamic-chi':f'분류: Counterexample candidate. 고정 구간의 네 변형 검사 판정은 {status}이며 통과 행은 {counts}다. 검사 밖 나이 구간이나 다른 형성 이력에 대한 결론으로 확장하지 않는다.\n\n분류: Proven. 원래 MESA mass_correction은 배포 원자량과 바리온 조성의 비율이며 내부에너지·결합에너지까지 포함한 ADM 질량이 아니다. 정확한 소스와 조성으로 이를 대조했다. 절대 에너지 기준·고유 부피·구조 반작용이 누락된 단계가 GR 질량 완성의 경계다.\n\n분류: Conjectural. 열 EOS의 절대 기준과 GR 구조를 먼저 완성한 후, 전체 force/readout 및 엄밀한 미분·관측 추론을 연결해야 한다.'}
    phase=json.loads((OUT/'phase-selection.json').read_text())
    later=', '.join(f"{name} {r['matching_rows']}개" for name,r in phase.items())
    for key in ['model-definition','observable-targets','failure-ledger-dynamic-chi']:
        additions[key]+=f'\n\n분류: Counterexample candidate. 원래 구간에서 실패한 두 실행을 대상으로 별도 등록한 다음 냉각 구간 검사에서는 {later}의 실제 후보를 확보했다. 원래 실패 판정은 유지한다. 나이 이동과 양립하는 후보 회복이며, 고정 나이의 수렴 인증은 아니다.'
    manifest=json.loads((OLD/'manifest.json').read_text())['sha256']
    dest=OUT/'request19-notes'; dest.mkdir(exist_ok=False); bindings={}
    rels=['docs/'+k+'.md' for k in additions]+['paper/revision-manifest.json']
    for rel in rels:
        assert sha(ROOT/rel)==manifest[rel]
        snap=dest/Path(rel).name; shutil.copy2(ROOT/rel,snap)
        bindings[rel]=dict(snapshot=snap.relative_to(ROOT).as_posix(),sha256=manifest[rel],historical_manifest=rel.startswith('paper/'))
    save('historical-note-bindings.json',bindings)
    revision=json.loads((ROOT/'paper/revision-manifest.json').read_text())
    for stem,body in additions.items():
        rel='docs/'+stem+'.md'
        with (ROOT/rel).open('ab') as f:
            f.write(('\n\n## Request 20 열 구조 민감도와 질량 정의\n\n'+body+'\n\n세부 근거: [한글 보고서](../notes/REQUEST20_THERMAL_ROBUSTNESS_KO.md).\n').encode())
        revision['sha256'][rel]=sha(ROOT/rel)
    revision['request20_supporting_note_update']=dict(evidence_manifest='outputs/thermal-robustness20/manifest.json',
        historical_notes='outputs/thermal-robustness20/historical-note-bindings.json',
        status='유한 열 구조 민감도 검사·질량 정의 감사·조건부 응답 차이 정리; GR matching 및 전체 인증·추론 미완료',
        artifact_status='원고 PDF와 ZIP은 Request 12 동결본이다.')
    (ROOT/'paper/revision-manifest.json').write_text(json.dumps(revision,ensure_ascii=False,indent=2)+'\n')


def seal():
    import mpmath,sympy,thermal_wd,thermal_restart,stellar_matching,nonzero_drive
    modules=[np,mpmath,sympy,thermal_wd,thermal_restart,stellar_matching,nonzero_drive]
    save('provenance.json',dict(before_task_checkpoint='389eaef',before_perturbation_theorem_checkpoint='9b006eb',
        before_post_primary_phase_followup_checkpoint='8b95a41',
        byte_preservation_repair_checkpoint='b3a8944',
        interpreter=sys.executable,versions={m.__name__:m.__version__ for m in [np,mpmath,sympy]},
        module_file_sha256={str(m.__file__):sha(m.__file__) for m in modules},
        prior_runtime_environment_sha256=sha(OLD/'environment.json'),
        input_sha256={'outputs/validated-variational/ivp.hex':sha(ROOT/'outputs/validated-variational/ivp.hex')},
        producer='verification/thermal_robustness.py'))
    paths=[p for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json']
    paths += [ROOT/'verification/thermal_robustness.py',ROOT/'notes/REQUEST20_THERMAL_ROBUSTNESS_KO.md',ROOT/'.gitattributes']
    paths += [ROOT/rel for rel in json.loads((OUT/'historical-note-bindings.json').read_text())]
    save('manifest.json',dict(classification='Proven',scope='Request 20 유한 열 구조 변형·질량 정의·조건부 scalar 응답 상계',
        sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(paths)}))


def verify():
    histories={
        'outputs/validated-variational/manifest.json':'outputs/remaining-levers15/historical-note-bindings.json',
        'outputs/remaining-levers15/manifest.json':'outputs/nbody-readout16/historical-note-bindings.json',
        'outputs/nbody-readout16/manifest.json':'outputs/nonzero-drive17/historical-note-bindings.json',
        'outputs/nonzero-drive17/manifest.json':'outputs/thermal-wd18/historical-note-bindings.json',
        'outputs/thermal-wd18/manifest.json':'outputs/thermal-restart19/historical-note-bindings.json',
        'outputs/thermal-restart19/manifest.json':'outputs/thermal-robustness20/historical-note-bindings.json'}
    count=0
    for label in ['outputs/research-remediation/manifest.json',*histories,'paper/revision-manifest.json','outputs/thermal-robustness20/manifest.json']:
        old=json.loads((ROOT/histories[label]).read_text()) if label in histories else {}
        for name,expected in json.loads((ROOT/label).read_text())['sha256'].items():
            path=ROOT/name
            if name in old:
                bind=old[name]; path=ROOT/bind['snapshot']; assert expected==bind['sha256']
                if bind.get('historical_manifest'):
                    before=json.loads(path.read_text()); after=json.loads((ROOT/name).read_text())
                    for k,value in before.items():
                        if k!='sha256': assert after[k]==value,k
                    for k,value in before['sha256'].items():
                        if k not in old: assert after['sha256'][k]==value,k
                else: assert (ROOT/name).read_bytes().startswith(path.read_bytes())
            assert sha(path)==expected,(label,name)
            count+=1
    provenance=json.loads((OUT/'provenance.json').read_text())
    for path,expected in provenance['module_file_sha256'].items(): assert sha(path)==expected,path
    for path,expected in provenance['input_sha256'].items(): assert sha(ROOT/path)==expected,path
    assert not json.loads((OUT/'gates.json').read_text())['complete_nonlinear_observational_inference']
    print('PASS:',count,'현재·역사적 SHA 및 동결 원고·미완료 경계 보존')


def progress():
    for name in [*CASES,*PHASE_CASES]:
        path=CACHE/name/'LOGS1/history.data'
        if not path.exists(): print(name,'no history'); continue
        with path.open() as f:
            head=[next(f) for _ in range(6)]; last=''
            for line in f:
                if len(line.split())==len(head[5].split()): last=line
        d=dict(zip(head[5].split(), last.split()))
        print(name,{k:d.get(k) for k in ['model_number','star_age','log_Teff','num_zones']},
            'finished' if (OUT/(name+'-execution.json')).exists() else 'running')


if __name__=='__main__':
    command=sys.argv[1]
    if command=='run': run(sys.argv[2])
    else: globals()[command]()
