"""Keep the failed small-box spectrum and resolve it with separate mesh/range tests."""
from concurrent.futures import ProcessPoolExecutor, as_completed
from decimal import Decimal
from types import FunctionType
import json, shutil, sys
import mpmath as mp
import numpy as np
import h2_spectre_data as original
import molecular_partition_data as partition

g=original.g;OUT=g.OUT/'gr-h2-spectre-refinement'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();original.bindings()
    baseline=json.loads((original.OUT/'baseline-result.json').read_text())
    keys=sorted({(r['v'],r['J']) for r in baseline['failed_levels']}|{(0,0),(0,31),(7,20)})
    shutil.copy2(original.OUT/'build.json',OUT/'build.json')
    save('baseline-manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p)
        for p in original.OUT.iterdir() if p.is_file()}))
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='3d70a3b',
        discovery='The baseline has 294/302 finite passes; eight failed small-box levels and the upstream weak-state warnings motivate this separately declared refinement. The baseline failures and its original gate remain frozen.',
        keys=keys,grids=[[800,40],[1600,40],[1200,60],[2400,60]],
        finite_energy_difference_gate_cm_inverse=1e-4,
        selection='Retain baseline 400/20 for levels passing both baseline controls. For all failed keys and three positive controls, require both spacing comparisons (800/40 vs 1600/40;1200/60 vs 2400/60) and both range comparisons (800/40 vs1200/60;1600/40 vs2400/60) below the original 1e-4 cm^-1 gate, and positive dissociation energy on all four grids. Use 2400/60 only after all four tests pass.',
        uniform_T_K=[999999,1000001],interval_digits=70,maximum_logT_derivative=10,
        printed_level_halfwidth_cm_inverse='0.000000000001',
        interval_scope='Conditional enclosure of the declared computed 302-level fixed spectrum. Excitation E=D00-DvJ; +/-1e-12 cm^-1 encloses two half-final printed digits. This does not bound DVR error, potential error, QED truncation, constants, excited electronic levels, continuum or plasma occupation. Upstream estimated uncertainties are reported separately and are not hard bounds.',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in [g.ROOT/'verification/h2_spectre_refinement.py',
            OUT/'baseline-manifest.json',OUT/'build.json',partition.OUT/'manifest.json']}))


def job(args):
    fn=original.job
    return FunctionType(fn.__code__,dict(fn.__globals__,OUT=OUT,save=save))(*[args])


def run():
    plan=json.loads((OUT/'plan.json').read_text());jobs=[]
    for N,R in plan['grids']:
        for v,J in plan['keys']:jobs.append((N,R,[(v,J)],f'fine-{N}-{R}-v{v:02}-j{J:02}'))
    with ProcessPoolExecutor(max_workers=4) as pool:
        for done in as_completed([pool.submit(job,args) for args in jobs]):done.result()
    collect()


def collect():
    plan=json.loads((OUT/'plan.json').read_text());rows=[];selected={}
    for v in range(15):
        for row in json.loads((original.OUT/f'grid-400-20-v{v:02}.json').read_text())['rows']:
            selected[row['v'],row['J']]=dict(row,grid=[400,20])
    for v,J in plan['keys']:
        series=[json.loads((OUT/f'fine-{N}-{R}-v{v:02}-j{J:02}.json').read_text())['rows'][0] for N,R in plan['grids']]
        D=[Decimal(r['dissociation_cm_inverse']) for r in series]
        differences=[float(abs(D[a]-D[b])) for a,b in [(0,1),(2,3),(0,2),(1,3)]]
        passed=all(r['positive'] for r in series) and max(differences)<plan['finite_energy_difference_gate_cm_inverse']
        rows.append(dict(v=v,J=J,spacing_40= differences[0],spacing_60=differences[1],
            range_coarse=differences[2],range_fine=differences[3],all_four_passed=passed))
        if passed:selected[v,J]=dict(series[-1],grid=plan['grids'][-1])
    result=dict(classification='Counterexample candidate',finite_refined_controls=rows,
        all_finite_controls_passed=all(r['all_four_passed'] for r in rows),physical_EOS_certified=False)
    save('result.json',result)
    if not result['all_finite_controls_passed']:return
    D0=Decimal(selected[0,0]['dissociation_cm_inverse']);levels=[]
    for key,r in sorted(selected.items()):
        levels.append(dict(r,excitation_cm_inverse=str(D0-Decimal(r['dissociation_cm_inverse']))))
    assert len(levels)==302 and all(Decimal(r['excitation_cm_inverse'])>=0 for r in levels)
    save('levels.json',dict(classification='Counterexample candidate',rows=levels,
        boundary='Energies calculated here with the unmodified upstream model; support imported from H2SPECTRE 7.4. Finite mesh/range agreement is not a rigorous physical error bound.'))
    spectrum=[(r['v'],r['J'],r['excitation_cm_inverse']) for r in levels]
    boxes=partition.interval_audit(spectrum,partition.coefficients(),plan)
    subset=boxes.pop('fixed_348_level_model');boxes['fixed_302_level_model']=subset;save('interval.json',boxes)
    mp.mp.dps=90;c2=mp.mpf('6.62607015e-34')*299792458*100/mp.mpf('1.380649e-23')
    args=[(mp.mpf(e),mp.mpf((2*J+1)*(1 if J%2==0 else 3))/4) for v,J,e in spectrum]
    def logQ(t):return mp.log(mp.fsum(w*mp.exp(-c2*e/mp.exp(t)) for e,w in args))
    controls=[]
    for T in [999999,1000000,1000001]:
        for n,box in enumerate(subset['logQ_derivatives']):
            value=mp.diff(logQ,mp.log(T),n);lo,hi=[mp.mpf(tuple(x)) for x in box['binary_endpoints']]
            assert lo<=value<=hi,(T,n);controls.append(dict(T_K=T,order=n,value=mp.nstr(value,90),contained=True))
    save('independent-controls.json',dict(classification='Counterexample candidate',rows=controls,all_contained=True))
    E=np.array([float(e) for v,J,e in spectrum]);weights=np.array([float(w) for e,w in args]);comparisons=[]
    for T in [1000,3000,9000,18000,20000,100000,1000000,32000000]:
        x=float(c2)*E/T;w=weights*np.exp(-x);p=w/w.sum();mean=p@x
        comparisons.append(dict(T_K=T,Q=float(w.sum()),DlogQ=float(mean),internal_Cv_over_kB=float(p@((x-mean)**2))))
    result.update(completed=True,computed_level_keys=302,ground_dissociation_cm_inverse=str(D0),
        weak_14_4_dissociation_cm_inverse=selected[14,4]['dissociation_cm_inverse'],
        maximum_upstream_estimated_uncertainty_cm_inverse=max(float(r['upstream_estimated_uncertainty_cm_inverse']) for r in levels),
        upstream_full_four_body_E2_keys=sum(r['E2_method']=='FULL' for r in levels),comparisons=comparisons,
        independent_controls=33,all_independent_controls_contained=True)
    save('result.json',result)
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify();print('H2 REFINED SPECTRUM',result['computed_level_keys'],result['weak_14_4_dissociation_cm_inverse'],comparisons[-2],flush=True)


def verify():
    original.bindings();plan=json.loads((OUT/'plan.json').read_text())
    for name,key in [('plan.json','bindings'),('baseline-manifest.json','sha256'),('manifest.json','sha256')]:
        for rel,digest in json.loads((OUT/name).read_text())[key].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    build=json.loads((OUT/'build.json').read_text());assert g.c.sha(build['executable'])==build['sha256']
    for folder in [original.OUT,OUT]:
        for path in folder.glob('*-execution.json'):
            r=json.loads(path.read_text());stem=path.name.removesuffix('-execution.json')
            assert r['returncode']==0 and g.c.sha(folder/(stem+'.log'))==r['raw_sha256']
            assert g.c.sha(folder/(stem+'.input'))==r['input_sha256']
            stdin=(folder/(stem+'.input')).read_text();keys=[tuple(map(int,line.split())) for line in stdin.splitlines()[1:]]
            parsed=original.parse((folder/(stem+'.log')).read_text(),keys)
            assert parsed==json.loads((folder/(stem+'.json')).read_text())['rows']
    result=json.loads((OUT/'result.json').read_text())
    assert result['completed'] and result['computed_level_keys']==302 and result['all_finite_controls_passed']
    assert result['all_independent_controls_contained']
    print('PASS 302 computed H2 levels, preserved baseline failures, independent finite grids and interval controls',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
