"""Apply the certified degree remainder to the unchanged three full-Q pilots."""
from fractions import Fraction as F
from types import FunctionType,SimpleNamespace
import json,sys
import gr_outer_product_pilot as previous
import gr_response_product_adaptive as adaptive

ROOT=previous.ROOT;OUT=previous.OUT.parent/'gr-outer-product-adaptive-pilot';sha=previous.sha


def save(name,value):(OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    adaptive.verify();previous.verify();assert not OUT.exists()
    timing=json.loads((adaptive.OUT/'benchmark.json').read_text())
    assert all(x['median_ratio']>1 for x in timing['rows'])
    OUT.mkdir();p=previous.bindings();p.update(checkpoint='ad93dce4',classification='Counterexample candidate',
        scope='Repeat all three exact full-Q integral pilots with the separately certified adaptive Legendre degree binary. All prior inputs,80-bit Q nodes,128-bit native/256-bit outer arithmetic, neutral-root corrections, IR/high-Q endpoints and1e-15/2e-7 gates are unchanged. The old completed pilot and whole-table producer remain untouched.',
        runtime_change='Only engine() now uses gr_response_product_adaptive.engine().runtime(), whose binary and shared libraries are hash-bound. The original evaluate/preflight/run/finalize/verify Python code objects are reused with this OUT,bindings and engine; no equation or gate substitutions.',
        extra_gate='Every complete physical-field interval must overlap its prior full-degree counterpart. Three paired finite benchmark medians must improve before this pilot begins; no all-state speed claim follows.',
        boundary='Three repeated states only. No whole3206-state completion, full physical EOS, nonlinear GR evolution or actual observation inference.')
    files=[ROOT/'verification/gr_outer_product_adaptive_pilot.py',ROOT/'verification/gr_outer_product_pilot.py',
        ROOT/'verification/gr_response_product_adaptive.py',adaptive.OUT/'manifest.json',adaptive.OUT/'runtime.json',adaptive.OUT/'effective-source.cpp',
        previous.OUT/'manifest.json']
    p['bindings'].update({f.relative_to(ROOT).as_posix():sha(f) for f in files});save('plan.json',p)


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert sha(ROOT/rel)==digest,rel
    return p


def engine():
    ns=dict(previous.__dict__,OUT=OUT,bindings=bindings,engine=adaptive.engine)
    for name in ['save','evaluate','preflight','verify_preflight','run','finalize','verify']:
        fn=getattr(previous,name);ns[name]=FunctionType(fn.__code__,ns,name,fn.__defaults__)
    return SimpleNamespace(**{name:ns[name] for name in ['preflight','verify_preflight','run','finalize','verify']})


def preflight():engine().preflight()


def run():engine().run()


def comparisons():
    rows=[]
    for i in bindings()['positions']:
        old=json.loads((previous.OUT/f'cell-{i:04d}-result.json').read_text());new=json.loads((OUT/f'cell-{i:04d}-result.json').read_text())
        assert old['passed'] and new['passed'] and old['nodes']==new['nodes'] and old['cell']==new['cell']
        assert (previous.OUT/f'cell-{i:04d}-q.tsv').read_bytes()==(OUT/f'cell-{i:04d}-q.tsv').read_bytes()
        widths=[]
        for x,y in zip(old['complete_field_enclosures'],new['complete_field_enclosures'],strict=True):
            a,b=previous.cusp.endpoints(x);c,d=previous.cusp.endpoints(y);width=min(b,d)-max(a,c);assert width>=0;widths.append(str(width))
        difference=max(abs(F.from_float(x)-F.from_float(y)) for x,y in zip(old['field_approximation'],new['field_approximation'],strict=True))
        rows.append(dict(position=i,cell=old['cell'],nodes=old['nodes'],exact_Q_byte_replay=True,intersection_widths=widths,maximum_midpoint_difference=str(difference),
            prior_native_seconds=old['native_seconds'],adaptive_native_seconds=new['native_seconds'],maximum_new_field_error=str(max(map(F,new['field_error_upper'])))))
    save('comparisons.json',dict(classification='Counterexample candidate',passed=True,components=21,rows=rows,
        scope='Same three physical model states, shared domains and endpoint certificates. Full-degree versus bounded-degree algorithm compatibility. Historical full-job times are descriptive, not the paired benchmark.'))


def finalize():comparisons();engine().finalize();verify()


def verify_preflight():engine().verify_preflight()


def verify():
    engine().verify();r=json.loads((OUT/'comparisons.json').read_text());assert r['passed'] and r['components']==21
    for row in r['rows']:assert row['exact_Q_byte_replay'] and min(map(F,row['intersection_widths']))>=0 and F(row['maximum_new_field_error'])<F('2e-7')
    print('PASS three full-Q adaptive-degree pilots at original response and physical-field gates',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
