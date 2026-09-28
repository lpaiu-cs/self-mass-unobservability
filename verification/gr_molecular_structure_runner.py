"""Serialize the measured NumPy predicate without changing any scientific gate."""
from types import ModuleType
import json, shutil, sys
import gr_molecular_structure as original

g=original.g;OUT=g.OUT/'gr-molecular-structure-defined'


def source():
    before=(g.ROOT/'verification/gr_molecular_structure.py').read_text()
    changes={"OUT=g.OUT/'gr-molecular-structure'":"OUT=g.OUT/'gr-molecular-structure-defined'",
        "passed=abs(actual-t)<plan['known_root_log_temperature_tolerance']":
        "passed=bool(abs(actual-t)<plan['known_root_log_temperature_tolerance'])"}
    after=before
    for a,b in changes.items():assert after.count(a)==1;after=after.replace(a,b)
    reverse=after
    for a,b in reversed(list(changes.items())):assert reverse.count(b)==1;reverse=reverse.replace(b,a)
    assert reverse==before;return after,changes


def module(text):
    # Register the reused worker functions so the existing process pool can pickle them.
    name='gr_molecular_structure_defined';obj=ModuleType(name);sys.modules[name]=obj
    exec(compile(text,str(OUT/'candidate.py'),'exec'),obj.__dict__);return obj


def prepare():
    assert not OUT.exists();text,changes=source();obj=module(text);obj.prepare()
    (OUT/'candidate.py').write_text(text)
    shutil.copy2(g.ROOT/'outputs/gr-molecular-structure33-table.log',OUT/'original-serialization-failure.log')
    plan=json.loads((OUT/'plan.json').read_text())
    plan.update(checkpoint='37389db',serialization_correction='The comparison produced a NumPy bool that stdlib JSON rejected before any table block started. Explicitly convert that predicate to a Python bool; preserve the original code/plan/startup/error log. All inputs, comparisons, seeds, algorithms and gates are unchanged.',substitutions=changes)
    plan['bindings'].update({p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in [
        g.ROOT/'verification/gr_molecular_structure_runner.py',OUT/'candidate.py',OUT/'original-serialization-failure.log',
        original.OUT/'plan.json',original.OUT/'startup-binding.json']})
    obj.save('plan.json',plan)
    # The specific boundary check that failed, now independently JSON round-tripped.
    assert json.loads(json.dumps(bool(original.np.float64(0)<1))) is True


def run(name):
    text,_=source();assert text==(OUT/'candidate.py').read_text();obj=module(text)
    obj.bindings();getattr(obj,name)()


if __name__=='__main__':
    if sys.argv[1]=='prepare':prepare()
    else:run(sys.argv[1])
