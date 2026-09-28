"""Change only TOPS display names after isolating the server's hyphen failure."""
from types import FunctionType
import json, shutil, sys
import lanl_tops_stellar_probe as original

g=original.g;OUT=g.OUT/'lanl-tops-stellar-names'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();original.verify()
    for name in ['lanl-tops-form.html','selected-inputs.json','stellar-inputs.npz']:
        shutil.copy2(original.OUT/name,OUT/name)
    for label in ['RealAlnum','HyphenOnly']:
        for suffix in ['.html','.request.json']:
            name='lanl-tops-interface-'+label+suffix;shutil.copy2(g.ROOT/'outputs'/name,OUT/name)
    rows=json.loads((original.OUT/'requests.json').read_text());before=json.loads(json.dumps(rows))
    for a,b in zip(rows,before):
        a['fields']['mixname']=b['fields']['mixname'].replace('-','')
        assert a['fields']['mixname'].isalnum() and len(a['fields']['mixname'])<=15
        assert dict(a['fields'],mixname=b['fields']['mixname'])==b['fields']
    save('requests.json',rows)
    plan=json.loads((original.OUT/'plan.json').read_text())
    plan['display_name_repair']='All seven original queries returned HTTP 500 at /submit. A one-field intervention made the actual outer-cell input return HTTP 200 when only its mixture display name lost hyphens. Adding a hyphen to the otherwise working aluminum control returned HTTP 500. Change only mixname to alphanumeric. All compositions, isotope masses, density, tabulated temperatures, library, spectral and cutoff settings are byte-for-byte unchanged. Preserve the original failed retrieval manifest.'
    paths=[g.ROOT/'verification/lanl_tops_name_runner.py',original.OUT/'manifest.json']+list(OUT.iterdir())
    plan['bindings'].update({p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths if p.is_file()})
    save('plan.json',plan)


def namespace():
    env=dict(vars(original),OUT=OUT)
    for name in ['save','bindings','fetch','run','verify']:
        fn=getattr(original,name);env[name]=FunctionType(fn.__code__,env,argdefs=fn.__defaults__)
    return env


def run():namespace()['run']()


def verify():namespace()['verify']()


if __name__=='__main__':globals()[sys.argv[1]]()
