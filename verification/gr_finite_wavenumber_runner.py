"""Preserve the structural-equality failure; compare the symbolic difference."""
import json, shutil, sys
import gr_finite_wavenumber_response as original

g=original.g;OUT=g.OUT/'gr-finite-wavenumber-defined'
OLD="assert sp.limit((x*x/gamma+x*(1+x*x-Q*Q)/(2*Q*gamma)*sp.log((x+Q)/(x-Q))),Q,0)==(2*x*x+1)/gamma"
NEW="assert sp.simplify(sp.limit((x*x/gamma+x*(1+x*x-Q*Q)/(2*Q*gamma)*sp.log((x+Q)/(x-Q))),Q,0)-(2*x*x+1)/gamma)==0"


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2,default=original.previous.scalar)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();original.bindings()
    source=(g.ROOT/'verification/gr_finite_wavenumber_response.py').read_text()
    assert source.count(OLD)==1
    (OUT/'effective-source.py').write_text(source.replace(OLD,NEW))
    error=g.ROOT/'outputs/gr-finite-wavenumber-response33-run.log'
    assert 'AssertionError' in error.read_text() and not (original.OUT/'states.npz').exists()
    shutil.copy2(error,OUT/'original-error.log')
    paths=[g.ROOT/'verification/gr_finite_wavenumber_runner.py',original.OUT/'plan.json',OUT/'effective-source.py',OUT/'original-error.log']
    save('plan.json',dict(classification='Counterexample candidate',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        original_failure='SymPy returned 2*x^2/sqrt(1+x^2)+1/sqrt(1+x^2), which is not structurally equal to (2*x^2+1)/sqrt(1+x^2). Their simplified difference is exactly zero. Failure occurred before numerical controls or stellar quadrature.',
        change='Change only structural == to simplify(lhs-rhs)==0 in the symbolic k=0 limit assertion. All physics, formulae, controls and numerical gates are identical.'))


def bindings():
    for rel,digest in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    source=(g.ROOT/'verification/gr_finite_wavenumber_response.py').read_text()
    assert source.count(OLD)==1 and (OUT/'effective-source.py').read_text()==source.replace(OLD,NEW)
    return original.bindings()


def module():
    ns={'__name__':'finite_response_defined'}
    exec(compile((OUT/'effective-source.py').read_text(),str(OUT/'effective-source.py'),'exec'),ns)
    ns.update(OUT=OUT,save=save,bindings=bindings)
    return ns


def run():module()['run']()
def verify():module()['verify']()


if __name__=='__main__':globals()[sys.argv[1]]()
