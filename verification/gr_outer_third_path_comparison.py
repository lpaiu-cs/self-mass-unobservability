"""Reuse the frozen center-adjusted comparison for the completed central state."""
import json,sys
from fractions import Fraction as F
import gr_outer_second_path_comparison as checks

ROOT=checks.ROOT;sha=checks.sha;OUT=checks.OUT.parent/'gr-outer-third-path-comparison'


def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2)+'\n')


def prepare():
    assert not OUT.exists();checks.verify();checks.paths.direct.verify()
    stem='cell-3205';files=[ROOT/'verification/gr_outer_third_path_comparison.py',
        ROOT/'verification/gr_outer_second_path_comparison.py',checks.OUT/'manifest.json',
        checks.paths.direct.OUT/'manifest.json',checks.paths.product.OUT/'manifest.json']
    for folder in [checks.paths.direct.OUT,checks.paths.product.OUT]:
        current=list(folder.glob(stem+'-*'));assert len(current)==6;files+=current
        row=json.loads((folder/(stem+'-result.json')).read_text())
        assert row['passed'] and row['position']==3205 and row['nodes']==19248
    prior=json.loads((checks.OUT/'plan.json').read_text())
    OUT.mkdir();save('plan.json',dict(prior,classification='Proven',position=3205,expected_nodes=19248,
        checkpoint='b173df96',bindings={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        node_gate='All 19248 Q values agree exactly; all three eta responses intersect after the same proved neutral-center shift. No fitted tolerance.',
        comparison='Reuse the unchanged second-path calculate function with only the output directory, position and declared node count changed.'))


def engine():
    checks.OUT=OUT;checks.POSITION=3205;return checks


def run():
    assert not (OUT/'result.json').exists();result=engine().calculate();save('result.json',result)
    assert result['passed'];save('manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    plan=engine().bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());assert r['passed'] and r['nodes']==plan['expected_nodes']==19248
    assert r['components']==57744 and r['disjoint_counts']==[0,0,0] and min(map(F,r['minimum_intersection_widths']))>=0
    assert r['field_comparison']['passed'] and r['disjoint_field_negative_control_passed']
    assert not r['physical_EOS_certified'] and not r['nonlinear_GR_evolution']
    print('PASS state3205: seven true-neutral fields and 57744 center-adjusted exact-Q components',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
