"""Correct symbolic equality and bind the actual complete-response schema."""
from types import ModuleType
import json,sys
import gr_response_regular_identity as original

OUT=original.OUT.parent/'gr-response-regular-defined'


def source():
    before=(original.ROOT/'verification/gr_response_regular_identity.py').read_text()
    changes={"OUT=sector.OUT.parent/'gr-response-regular-identity'":"OUT=sector.OUT.parent/'gr-response-regular-defined'",
        "assert sp.cancel(2*x-(2*x+x*x)/(1+x)**2)==x*x*(2*x+3)/(x+1)**2":"assert sp.cancel(2*x-(2*x+x*x)/(1+x)**2-x*x*(2*x+3)/(x+1)**2)==0",
        "if row['case'] in wanted:source[row['case']]=row":"case=12*row['position']+row['z_index']\n            if case in wanted:source[case]=row",
        "texts=row['enclosures'];inside=all(M(a)<=value<=M(b) for (a,b),value in zip(texts,values))":"texts=row['enclosures'];inside=all(M(a)<=value<=M(b) for (a,b),value in zip(map(cusp.endpoints,texts),values))"}
    after=before
    for a,b in changes.items():assert after.count(a)==1;after=after.replace(a,b)
    reverse=after
    for a,b in reversed(list(changes.items())):reverse=reverse.replace(b,a)
    assert reverse==before;return after,changes


def module(text):
    obj=ModuleType('gr_response_regular_defined');exec(compile(text,str(OUT/'candidate.py'),'exec'),obj.__dict__);return obj


def prepare():
    assert not OUT.exists();text,changes=source();obj=module(text);obj.prepare();(OUT/'candidate.py').write_text(text)
    p=json.loads((OUT/'plan.json').read_text());p.update(checkpoint='0ca89c3',execution_corrections=changes,
        correction_scope='The first symbolic run stopped at a structurally different but algebraically equal SymPy expression. Compare the simplified difference to zero. The complete-response records have position/z_index and interval-text fields; bind those actual keys before running controls. Equations, mathematical bounds, state selection and gates are unchanged.')
    files=[original.ROOT/'verification/gr_response_regular_runner.py',OUT/'candidate.py',original.OUT/'plan.json',original.OUT/'prove-failure.log']
    p['bindings'].update({f.relative_to(original.ROOT).as_posix():original.g.sha(f) for f in files});obj.save('plan.json',p)


def run(name):
    text,_=source();assert text==(OUT/'candidate.py').read_text();obj=module(text);obj.bindings();getattr(obj,name)()


def verify():run('verify')


if __name__=='__main__':
    if sys.argv[1]=='prepare':prepare()
    else:run(sys.argv[1])
