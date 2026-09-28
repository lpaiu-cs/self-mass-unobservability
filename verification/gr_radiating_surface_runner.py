"""Expand the Lorentz polynomial before the exact W-squared substitution."""
from types import ModuleType
import sys
import gr_radiating_surface_connection as original

OUT=original.OUT.parent/'gr-radiating-surface-defined'


def source():
    before=(original.g.ROOT/'verification/gr_radiating_surface_connection.py').read_text()
    changes={"OUT=g.OUT/'gr-radiating-surface-connection'":"OUT=g.OUT/'gr-radiating-surface-defined'",
        'ratio=((gamma+ur)**2).subs(W**2,1/(1-v*v))':'ratio=sp.expand((gamma+ur)**2).subs(W**2,1/(1-v*v))',
        "with urllib.request.urlopen(URL,timeout=60) as response:(OUT/'maharaj-govender-govender2013.pdf').write_bytes(response.read())":
        "(OUT/'maharaj-govender-govender2013.pdf').write_bytes((g.OUT/'gr-radiating-surface-connection/maharaj-govender-govender2013.pdf').read_bytes())"}
    after=before
    for old,new in changes.items():assert after.count(old)==1;after=after.replace(old,new)
    return after,changes


def module(text):
    obj=ModuleType('gr_radiating_surface_defined');exec(compile(text,str(OUT/'candidate.py'),'exec'),obj.__dict__);return obj


def run():
    assert not OUT.exists();text,changes=source();obj=module(text);old_save=obj.save
    def save(name,value):
        if name=='plan.json':
            (OUT/'candidate.py').write_text(text)
            error=original.g.ROOT/'outputs/gr-radiating-surface-connection33-run.log';assert 'AssertionError' in error.read_text()
            files=[original.g.ROOT/'verification/gr_radiating_surface_runner.py',OUT/'candidate.py',original.OUT/'plan.json',error]
            value.update(checkpoint='0f63f349',substitutions=changes,
                correction='The first check tried to replace W^2 inside the unexpanded square (W/a+W*v/a)^2; no W^2 node existed in that expression tree. Expand the polynomial before the identical substitution W^2=1/(1-v^2). All metric, junction and luminosity equations, saved inputs and gates are unchanged. Replay the already-read source PDF bytes without a new download.')
            value['bindings'].update({p.relative_to(original.g.ROOT).as_posix():original.g.c.sha(p) for p in files})
        old_save(name,value)
    obj.save=save;obj.run()


def verify():
    text,_=source();assert text==(OUT/'candidate.py').read_text();module(text).verify()


if __name__=='__main__':globals()[sys.argv[1]]()
