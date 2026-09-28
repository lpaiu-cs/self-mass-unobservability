"""Compare completed true-neutral field intervals without conflating eta centers."""
from fractions import Fraction as F
import json,sys
import gr_outer_product_pilot as product
import gr_outer_self_refined as direct

ROOT=product.ROOT;OUT=product.OUT.parent/'gr-outer-path-comparison';sha=product.sha


def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2)+'\n')


def compare(a,b):
    assert a['passed'] and b['passed'] and a['complete_finite_interval'] and b['complete_finite_interval']
    assert (a['position'],a['cell'],a['nodes'])==(b['position'],b['cell'],b['nodes'])
    fields=product.outer.highq.FIELDS;assert len(a['complete_field_enclosures'])==len(b['complete_field_enclosures'])==len(fields)==7
    rows=[]
    for name,x,y in zip(fields,a['complete_field_enclosures'],b['complete_field_enclosures'],strict=True):
        al,ah=product.cusp.endpoints(x);bl,bh=product.cusp.endpoints(y);lo=max(al,bl);hi=min(ah,bh)
        rows.append(dict(field=name,intersection=[str(lo),str(hi)],overlap=lo<=hi,direct_subset_product=bl<=al<=ah<=bh))
    return dict(position=a['position'],cell=a['cell'],nodes=a['nodes'],rows=rows,passed=all(x['overlap'] for x in rows),
        maximum_midpoint_difference=str(max(abs(F.from_float(x)-F.from_float(y)) for x,y in zip(a['field_approximation'],b['field_approximation'],strict=True))))


def run():
    assert not OUT.exists();product.verify_preflight()
    text,_=direct.source();assert text==(direct.OUT/'candidate.py').read_text();obj=direct.module(text);obj.runtime()
    for rel,digest in json.loads((direct.OUT/'preflight-manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    OUT.mkdir()
    files=[ROOT/'verification/gr_outer_path_comparison.py',ROOT/'verification/gr_outer_product_pilot.py',
        ROOT/'verification/gr_outer_self_refined.py',product.OUT/'preflight-manifest.json',direct.OUT/'preflight-manifest.json']
    # Only completed state 0 is frozen here; the other direct workers keep their files.
    files+=list(product.OUT.glob('cell-0000-*'))+list(direct.OUT.glob('cell-0000-*'))
    assert len(list(product.OUT.glob('cell-0000-*')))==6 and len(list(direct.OUT.glob('cell-0000-*')))==6
    save('plan.json',dict(classification='Proven',positions=[0],bindings={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        gate='All seven true-neutral full-Q field intervals must intersect. Each path independently retains its original 2e-7 field gate.',
        scope='Two different momentum representations, not an independent end-to-end proof: they share the outer Gauss rule, analytic domains, IR/high-Q certificates, physical field definitions and some interval infrastructure. Raw middle jets use different eta centers and must not be required to intersect.'))
    a=json.loads((direct.OUT/'cell-0000-result.json').read_text());b=json.loads((product.OUT/'cell-0000-result.json').read_text())
    r=compare(a,b);assert r['passed']
    for item in [a,b]:
        for text,point,error in zip(item['complete_field_enclosures'],item['field_approximation'],item['field_error_upper'],strict=True):
            lo,hi=product.cusp.endpoints(text);assert max(abs(lo-F.from_float(point)),abs(hi-F.from_float(point)))==F(error)<F('2e-7')
    bad=json.loads(json.dumps(b));bad['complete_field_enclosures'][0]='[1, 2]';assert not compare(a,bad)['passed']
    save('result.json',dict(classification='Proven',passed=r['passed'],rows=[r],disjoint_negative_control_passed=True,physical_EOS_certified=False))
    save('manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    p=json.loads((OUT/'plan.json').read_text());m=json.loads((OUT/'manifest.json').read_text())
    for rel,digest in {**p['bindings'],**m['sha256']}.items():assert sha(ROOT/rel)==digest,rel
    a=json.loads((direct.OUT/'cell-0000-result.json').read_text());b=json.loads((product.OUT/'cell-0000-result.json').read_text())
    r=json.loads((OUT/'result.json').read_text());assert r['passed'] and r['rows']==[compare(a,b)] and not r['physical_EOS_certified']
    print('PASS two full-Q momentum paths: all seven true-neutral field intervals overlap; shared endpoint/domain dependencies remain',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
