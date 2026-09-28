"""Extend the frozen product outer evaluator to every neutral state."""
from concurrent.futures import ProcessPoolExecutor,as_completed
from fractions import Fraction as F
from pathlib import Path
from types import ModuleType,SimpleNamespace
import gzip,hashlib,json,math,shutil,subprocess,sys
from mpmath import iv
import gr_outer_product_pilot as pilot
import gr_response_H_table_runner as H_table

ROOT=pilot.ROOT;OUT=pilot.OUT.parent/'gr-outer-product-table';CACHE=pilot.CACHE.parent/'outer-product-table'
sha=pilot.sha;I=pilot.I;high=pilot.high


def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2)+'\n')


def schedule(position,row,state,data):
    """The original exact dyadic subdivision and 128-bit outward remainder."""
    precision=iv.prec;iv.prec=128
    try:
        qD=F(row['screening_disk_radius']);s=F(row['strip_radius_lower'])
        a=F.from_float(float(state['scale'][position]))/2**24;b=4*F.from_float(float(data['pcut'][position]));length=b-a
        norm=I(1)/((2*16+1)*math.comb(32,16)**2);field=I(max(map(F,row['physical_integrand_majorants'])))
        stack=[(a,b)];panels=[];error=I(0)
        while stack:
            lo,hi=stack.pop();center=(lo+hi)/2;h=(hi-lo)/2;radius=max(qD if center<=2*qD else min(s,center/2),center/256)
            if radius>h:
                remainder=I(hi-lo)*norm*(I(hi-lo)/I(radius-h))**32*field
                if high(remainder)<=F('2e-8')*(hi-lo)/length:
                    panels.append([str(lo),str(hi)]);error+=remainder;continue
            assert h>0 and len(panels)<100000;stack.extend([(center,hi),(lo,center)])
        panels.sort(key=lambda x:F(x[0]));assert F(panels[0][0])==a and F(panels[-1][1])==b
        assert all(F(x[1])==F(y[0]) for x,y in zip(panels,panels[1:]));assert high(error)<F('2e-8')
        return dict(position=position,cell=row['cell'],panels=panels,nodes=16*len(panels),field_remainder_upper=str(high(error)),passed=True)
    finally:iv.prec=precision


def source():
    before=(ROOT/'verification/gr_outer_product_pilot.py').read_text()
    changes={
        "meta=json.loads((outer.domain.OUT/'result.json').read_text());schedule=next(x for x in meta['schedules'] if x['position']==position);row=meta['rows'][position]":
        "meta=json.loads((outer.domain.OUT/'result.json').read_text());row=meta['rows'][position];schedule=full_schedule(position,row,state,data)",
        "cmd=[binary,'response',str(product.OUT/'inputs.tsv'),str(product.OUT/'rule.tsv'),str(target),str(position),str(product.OUT/f'h-{position:04d}.tsv'),str(qpath)]":
        "cmd=[binary,'response',str(FULL_INPUTS),str(product.OUT/'rule.tsv'),str(target),str(position),str(CURRENT_H),str(qpath)]"}
    after=before
    for old,new in changes.items():assert after.count(old)==1;after=after.replace(old,new)
    return after,changes


def prepare():
    assert not OUT.exists() and not CACHE.exists();pilot.verify_preflight();H_table.verify_preflight()
    OUT.mkdir();CACHE.mkdir();text,changes=source();(OUT/'candidate.py').write_text(text)
    files=[ROOT/'verification/gr_outer_product_table.py',ROOT/'verification/gr_outer_product_pilot.py',OUT/'candidate.py',
        pilot.OUT/'preflight-manifest.json',pilot.OUT/'plan.json',pilot.OUT/'runtime.json',
        pilot.outer.domain.OUT/'manifest.json',pilot.outer.domain.OUT/'result.json',H_table.OUT/'plan.json',
        H_table.OUT/'inputs.tsv',H_table.OUT/'preflight-manifest.json',ROOT/'verification/gr_response_H_table_runner.py']
    p=pilot.bindings();p.update(checkpoint=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        cells=3206,positions=list(range(3206)),processes=4,preflight_position=1,pilot_positions=[0,1603,3205],
        substitutions=changes,bindings={x.relative_to(ROOT).as_posix():sha(x) for x in files},
        target='Every existing state, with the same native H/product response, exact Q nodes, joint Cauchy domain, true-neutral displacement and 2e-7 final field gate.',
        schedule='Reproduce the original 128-bit outward subdivision at all states. Require exact replay of all three original schedules before a new state is evaluated.',
        H_binding='Each worker requires its H block result to exist and pass. Bind the block bytes, compressed H hash and decompressed H hash before use. No missing, failed or incomplete H table is skipped or replaced. Full finalization additionally requires the complete H table manifest and independent H controls.',
        reuse='Only the three completed, verified original product pilots and the single new preflight state are reused. All other outputs must be new. No generic resume, adaptive gate change or retry after a scientific failure.',
        preservation='Compress new node schedules, Q inputs and native responses losslessly, replay their uncompressed hashes, then remove only the newly created plain copies. Retain point intervals, result, native log, H dependencies and all hashes. Cache H copies stay inside this new worker directory.',
        boundary='A complete pass certifies this ideal-electron RPA self term on the specified finite 3206-state inventory only. It does not supply missing correlations, a partially ionized free energy, transport calibration, self-consistent GR dynamics or observations.')
    save('plan.json',p)


def bindings():
    p=json.loads((OUT/'plan.json').read_text());text,_=source();assert text==(OUT/'candidate.py').read_text()
    for rel,digest in p['bindings'].items():assert sha(ROOT/rel)==digest,rel
    return p


def schedule_controls():
    bindings();state,_,_,data,_,_=pilot.outer.highq.context();meta=json.loads((pilot.outer.domain.OUT/'result.json').read_text());checks=[]
    for original in meta['schedules']:
        i=original['position'];current=schedule(i,meta['rows'][i],state,data);assert current==original
        checks.append(dict(position=i,nodes=current['nodes'],exact_schedule_replay=True))
    save('schedule-controls.json',dict(classification='Proven',passed=True,checks=checks));print('PASS three exact full schedule replays',flush=True)


def evaluate(position):
    p=bindings();folder=OUT/f'cell-{position:04d}';scratch=CACHE/f'cell-{position:04d}';assert not folder.exists() and not scratch.exists()
    hp=json.loads((H_table.OUT/'plan.json').read_text());block=H_table.OUT/f'block-{position//hp["block_size"]*hp["block_size"]:04d}'/'result.json'
    block_sha=sha(block);record=json.loads(block.read_text());assert record['passed'];hr=next(x for x in record['rows'] if x['position']==position);assert hr['passed']
    compressed=block.parent/f'h-{position:04d}.tsv.gz';assert sha(compressed)==hr['gzip_sha256'];folder.mkdir();scratch.mkdir();hpath=scratch/'H.tsv'
    with gzip.open(compressed,'rb') as src,hpath.open('xb') as dst:shutil.copyfileobj(src,dst)
    assert sha(hpath)==hr['uncompressed_sha256']
    obj=ModuleType(f'gr_outer_product_full_{position}');exec(compile((OUT/'candidate.py').read_text(),str(OUT/'candidate.py'),'exec'),obj.__dict__)
    obj.OUT=folder;obj.bindings=lambda:p;obj.engine=lambda:SimpleNamespace(runtime=pilot.engine().runtime)
    obj.full_schedule=schedule;obj.FULL_INPUTS=H_table.OUT/'inputs.tsv';obj.CURRENT_H=hpath
    result=obj.evaluate(position);assert result['cell']==hr['cell']
    replays={}
    for suffix in ['-nodes.json','-q.tsv','-native.jsonl']:
        path=folder/f'cell-{position:04d}{suffix}';digest=sha(path);target=path.with_name(path.name+'.gz')
        with path.open('rb') as src,gzip.open(target,'xb',compresslevel=1) as dst:shutil.copyfileobj(src,dst)
        with gzip.open(target,'rb') as src:assert hashlib.sha256(src.read()).hexdigest()==digest
        replays[path.name]=dict(uncompressed_sha256=digest,gzip_sha256=sha(target));path.unlink()
    # Both resolved paths are worker-owned descendants of the declared cache.
    assert hpath.resolve().parent==scratch.resolve() and scratch.resolve().parent==CACHE.resolve();hpath.unlink();scratch.rmdir()
    dependency=dict(classification='Proven',block=block.relative_to(ROOT).as_posix(),block_sha256=block_sha,H_table=compressed.relative_to(ROOT).as_posix(),H_record=hr,compressed_replays=replays)
    (folder/'dependencies.json').write_text(json.dumps(dependency,indent=2)+'\n')
    manifest={x.relative_to(ROOT).as_posix():sha(x) for x in folder.iterdir() if x.is_file()}
    (folder/'manifest.json').write_text(json.dumps(dict(sha256=manifest),indent=2)+'\n')
    verify_cell(position);return result


def verify_cell(position):
    folder=OUT/f'cell-{position:04d}';d=json.loads((folder/'dependencies.json').read_text())
    for rel,digest in json.loads((folder/'manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    assert sha(ROOT/d['block'])==d['block_sha256'] and sha(ROOT/d['H_table'])==d['H_record']['gzip_sha256']
    r=json.loads((folder/f'cell-{position:04d}-result.json').read_text());assert r['passed'] and r['point_budget_passed'] and r['complete_finite_interval']
    assert r['finite_reference_control']['passed'] and not r['physical_EOS_certified'] and F(r['response_midpoint_error_upper'])<F('1e-15')
    for text,point,error in zip(r['complete_field_enclosures'],r['field_approximation'],r['field_error_upper'],strict=True):
        a,b=pilot.cusp.endpoints(text);assert max(abs(a-F.from_float(point)),abs(b-F.from_float(point)))==F(error)<F('2e-7')
    return r


def preflight():
    schedule_controls();p=bindings();evaluate(p['preflight_position'])
    files=[OUT/'plan.json',OUT/'candidate.py',OUT/'schedule-controls.json']+list((OUT/f'cell-{p["preflight_position"]:04d}').iterdir())
    save('preflight-manifest.json',dict(sha256={x.relative_to(ROOT).as_posix():sha(x) for x in files}));verify_preflight()


def verify_preflight():
    p=bindings()
    for rel,digest in json.loads((OUT/'preflight-manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    assert json.loads((OUT/'schedule-controls.json').read_text())['passed'];r=verify_cell(p['preflight_position'])
    print('PASS new full-table preflight state',r['position'],r['nodes'],flush=True)


def run():
    verify_preflight();pilot.verify();p=bindings();excluded=set(p['pilot_positions']+[p['preflight_position']]);completed=sorted(excluded);failures=[]
    with ProcessPoolExecutor(max_workers=p['processes']) as pool:
        jobs={pool.submit(evaluate,i):i for i in p['positions'] if i not in excluded}
        for future in as_completed(jobs):
            i=jobs[future]
            try:r=future.result();completed.append(i);print('FULL PRODUCT',i,len(completed),r['passed'],flush=True)
            except Exception as error:failures.append(dict(position=i,type=type(error).__name__,message=str(error)));print('FULL PRODUCT FAILURE',i,repr(error),flush=True)
            save('progress.json',dict(classification='Counterexample candidate',completed_positions=sorted(completed),failures=failures))
    assert not failures


def finalize():
    p=bindings();verify_preflight();pilot.verify();H_table.verify();rows=[];files=[OUT/'plan.json',OUT/'candidate.py',OUT/'schedule-controls.json',OUT/'preflight-manifest.json',pilot.OUT/'manifest.json',H_table.OUT/'manifest.json']
    for i in p['positions']:
        if i in p['pilot_positions']:
            path=pilot.OUT/f'cell-{i:04d}-result.json';r=json.loads(path.read_text());files.append(path)
        else:r=verify_cell(i);files.append(OUT/f'cell-{i:04d}'/'manifest.json')
        assert r['position']==i and r['passed'];rows.append(r)
    assert len({r['cell'] for r in rows})==p['cells']
    save('result.json',dict(classification='Proven',passed=True,cells=len(rows),rows=rows,
        maximum_field_error_upper=str(max(F(e) for r in rows for e in r['field_error_upper'])),physical_EOS_certified=False))
    files.append(OUT/'result.json');save('manifest.json',dict(sha256={x.relative_to(ROOT).as_posix():sha(x) for x in files}));verify()


def verify():
    p=bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());assert r['passed'] and r['cells']==p['cells'] and F(r['maximum_field_error_upper'])<F('2e-7')
    for i in p['positions']:
        if i not in p['pilot_positions']:verify_cell(i)
    print('PASS all 3206 full-Q neutral self integrals; full physical EOS, GR dynamics and observational closure remain separate',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
