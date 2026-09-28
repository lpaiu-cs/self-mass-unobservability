"""Build all reusable reference-center H tables with certified uniform errors."""
from concurrent.futures import ProcessPoolExecutor,as_completed
from fractions import Fraction as F
from types import FunctionType
import gzip,hashlib,json,shutil,subprocess,sys
import numpy as np
import gr_response_product_128 as product

g=product.g;ROOT=g.ROOT;OUT=product.OUT.parent/'gr-response-H-table-defined'


def save(name,value):(OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    product.verify();assert not OUT.exists();OUT.mkdir();data=dict(np.load(g.rule.highq.moments.OUT/'inputs.npz'));state=dict(np.load(g.cusp.THERMO/'states.npz'))
    assert len(state['cells'])==3206 and len(set(map(int,state['cells'])))==3206 and np.array_equal(data['cells'],state['cells'])
    (OUT/'inputs.tsv').write_text(''.join(' '.join(map(str,[i,int(data['cells'][i]),float(data['beta'][i]).hex(),float(data['Sref'][i]).hex(),float(state['eta'][i]).hex(),float(data['pcut'][i]).hex()]))+'\n' for i in range(3206)))
    files=[ROOT/'verification/gr_response_H_table.py',product.OUT/'manifest.json',product.OUT/'runtime.json',product.OUT/'effective-source.cpp',
        product.OUT/'rule.tsv',OUT/'inputs.tsv',g.rule.highq.moments.OUT/'inputs.npz',g.cusp.THERMO/'states.npz']
    files+=list(product.OUT.glob('h-*.tsv'))
    save('plan.json',dict(classification='Proven',checkpoint=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        cells=3206,block_size=32,processes=4,coefficient_bits=128,analytic_response_budget='1e-18',total_response_budget='1e-17',
        H_control_digits=80,H_control_x=['-1/2','0','1/2'],control_positions=[0,256,1024,1603,2048,3072,3205],
        bindings={x.relative_to(ROOT).as_posix():g.sha(x) for x in files},
        target='Build actual piecewise H/Sref and five eta/tau jets at every stored reference center. Reuse the unchanged certified native build mode. Audit exact coverage, six analytic errors, coefficient endpoint ordering and uniform error of the exact coefficient-midpoint polynomial.',
        interval_error='For every panel and component, |P_n(x)|<=1 on [-1,1] bounds coefficient-midpoint rounding by half the sum of coefficient widths. Add the certified interpolation/omission error. The positive kernel gives response error <=P*maximum_panel_error. Require analytic <=1e-18 and combined <=1e-17; this does not relax the separate final 1e-15 response gate.',
        preservation='Compress each new complete native table losslessly. Record the SHA of its uncompressed bytes and compressed artifact. The original three reference tables must replay bitwise. Only the runner-created temporary uncompressed copy is removed after compression and hash replay.',
        boundary='All stored eta centers, not a new true-root response or outer-integral certificate. No interacting/partially-ionized EOS, physical transport, GR or observation claim.'))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert g.sha(ROOT/rel)==digest,rel
    return p


def native():
    text,_,_=product.sources();return product.module(text).runtime()


def audit_table(path,pcut,plan):
    maximum=[F(0)]*6;analytic=[F(0)]*6;end=F(0);panels=0;omitted=0
    for line in path.read_text().splitlines():
        words=line.split();a,b=map(F,words[:2]);skip=int(words[2]);assert a==end<b<=pcut and skip in [0,1] and len(words)==(9 if skip else 393)
        errors=list(map(F,words[3:9]));assert min(errors)>=0;radii=[F(0)]*6
        if not skip:
            for n in range(32):
                for k in range(6):
                    lo,hi=map(F,words[9+12*n+2*k:11+12*n+2*k]);assert lo<=hi;radii[k]+=(hi-lo)/2
        for k in range(6):analytic[k]=max(analytic[k],errors[k]);maximum[k]=max(maximum[k],errors[k]+radii[k])
        end=b;panels+=1;omitted+=skip
    assert end==pcut;analytic=[pcut*x for x in analytic];maximum=[pcut*x for x in maximum]
    return dict(panels=panels,bounded_omissions=omitted,analytic_response_errors=list(map(str,analytic)),combined_response_errors=list(map(str,maximum)),
        passed=max(analytic)<F(plan['analytic_response_budget']) and max(maximum)<F(plan['total_response_budget']))


def block(first):
    p=bindings();binary=native();data=dict(np.load(g.rule.highq.moments.OUT/'inputs.npz'));folder=OUT/f'block-{first:04d}';assert not folder.exists();folder.mkdir();rows=[]
    for i in range(first,min(first+p['block_size'],p['cells'])):
        target=folder/f'h-{i:04d}.tsv';done=subprocess.run([binary,'build',str(OUT/'inputs.tsv'),str(product.OUT/'rule.tsv'),str(target),str(i)],capture_output=True,text=True)
        (folder/f'h-{i:04d}-native.log').write_text(done.stdout+done.stderr);assert done.returncode==0,done.stderr
        summary=json.loads(done.stdout);audit=audit_table(target,F.from_float(float(data['pcut'][i])),p);assert summary['panels']==audit['panels'] and summary['bounded_omissions']==audit['bounded_omissions']
        original=product.OUT/f'h-{i:04d}.tsv';replay=None
        if original.exists():replay=target.read_bytes()==original.read_bytes();assert replay
        digest=g.sha(target);compressed=target.with_suffix('.tsv.gz')
        with target.open('rb') as source,gzip.open(compressed,'xb',compresslevel=1) as sink:shutil.copyfileobj(source,sink)
        with gzip.open(compressed,'rb') as source:assert hashlib.sha256(source.read()).hexdigest()==digest
        size=target.stat().st_size;target.unlink()
        rows.append(dict(position=i,cell=int(data['cells'][i]),**audit,coefficient_seconds=summary['seconds'],uncompressed_sha256=digest,
            gzip_sha256=g.sha(compressed),uncompressed_bytes=size,compressed_bytes=compressed.stat().st_size,prior_bitwise_replay=replay))
    result=dict(classification='Proven',first=first,rows=rows,passed=all(x['passed'] for x in rows));(folder/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def preflight():
    r=block(0);assert r['passed']
    files=[OUT/'plan.json',OUT/'inputs.tsv']+list((OUT/'block-0000').iterdir())
    save('preflight-manifest.json',dict(sha256={x.relative_to(ROOT).as_posix():g.sha(x) for x in files}));verify_preflight()


def verify_preflight():
    bindings();native()
    for rel,digest in json.loads((OUT/'preflight-manifest.json').read_text())['sha256'].items():assert g.sha(ROOT/rel)==digest,rel
    r=json.loads((OUT/'block-0000/result.json').read_text());assert r['passed'] and [x['position'] for x in r['rows']]==list(range(32))
    assert r['rows'][0]['prior_bitwise_replay']
    print('PASS first 32 full H tables and original state-0 byte replay',flush=True)


def run():
    verify_preflight();p=bindings();completed=[0];failures=[]
    with ProcessPoolExecutor(max_workers=p['processes']) as pool:
        jobs=[pool.submit(block,first) for first in range(p['block_size'],p['cells'],p['block_size'])]
        for future in as_completed(jobs):
            r=future.result();completed.append(r['first']);failures.extend(x['position'] for x in r['rows'] if not x['passed'])
            save('progress.json',dict(classification='Counterexample candidate',completed_blocks=sorted(completed),failed_positions=failures));print('H TABLE BLOCK',r['first'],len(completed),r['passed'],flush=True)


def controls():
    p=bindings();folder=OUT/'controls';assert not folder.exists();folder.mkdir()
    for i in p['control_positions']:
        source=OUT/f'block-{i//p["block_size"]*p["block_size"]:04d}'/f'h-{i:04d}.tsv.gz'
        with gzip.open(source,'rb') as stream:(folder/f'h-{i:04d}.tsv').write_bytes(stream.read())
    ns=dict(g.__dict__,OUT=folder,bindings=lambda:dict(p,positions=p['control_positions']))
    ns['panels']=FunctionType(g.panels.__code__,ns)
    ns['save']=lambda name,value:(folder/name).write_text(json.dumps(value,indent=2)+'\n')
    FunctionType(g.H_controls.__code__,ns)()
    for i in p['control_positions']:(folder/f'h-{i:04d}.tsv').unlink()


def finalize():
    p=bindings();rows=[]
    for first in range(0,p['cells'],p['block_size']):rows.extend(json.loads((OUT/f'block-{first:04d}/result.json').read_text())['rows'])
    assert [r['position'] for r in rows]==list(range(p['cells']))
    controls=json.loads((OUT/'controls/H-controls.json').read_text());assert controls['passed']
    maxima=[max(F(r['combined_response_errors'][k]) for r in rows) for k in range(6)]
    save('result.json',dict(classification='Proven',cells=len(rows),passed=all(r['passed'] for r in rows),rows=rows,maximum_combined_response_errors=list(map(str,maxima)),
        independent_H_components=controls['components'],total_compressed_bytes=sum(r['compressed_bytes'] for r in rows),
        total_panels=sum(r['panels'] for r in rows),true_neutral_outer_integrals_certified=False,physical_EOS_certified=False))
    files=[x for x in OUT.rglob('*') if x.is_file()];save('manifest.json',dict(sha256={x.relative_to(ROOT).as_posix():g.sha(x) for x in files}));verify()


def verify():
    verify_preflight();p=bindings();native()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.sha(ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());assert r['passed'] and r['cells']==p['cells'] and max(map(F,r['maximum_combined_response_errors']))<F(p['total_response_budget'])
    assert r['independent_H_components']==6*9*len(p['control_positions'])
    print('PASS all 3206 reusable reference-center H tables and independent controls; outer/EOS/GR closure remains separate',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
