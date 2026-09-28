"""Remove the integrable logarithm analytically in cold reference controls."""
import json, shutil, sys
import gr_finite_wavenumber_runner as previous

g=previous.g;OUT=g.OUT/'gr-finite-wavenumber-regular'
OLD='''            def f(p):
                if p==0:return mp.mpf(0)
                return p*p/mp.sqrt(1+p*p)+p*(1+p*p-Q*Q)/(2*Q*mp.sqrt(1+p*p))*mp.log(abs((p+Q)/(p-Q)))
            integral=mp.quad(f,points);source=cold(x,Q)'''
NEW='''            bq=1/(2*mp.sqrt(1+Q*Q))
            def f(p):
                gamma=mp.sqrt(1+p*p);b=p*(1+p*p-Q*Q)/(2*Q*gamma)
                regular=p*p/gamma+b*mp.log(p+Q)
                return regular if p==Q else regular-(b-bq)*mp.log(abs(p-Q))
            def primitive(t):return t*(mp.log(abs(t))-1) if t else mp.mpf(0)
            integral=mp.quad(f,points)-bq*(primitive(x-Q)-primitive(-Q));source=cold(x,Q)'''


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2,default=previous.original.previous.scalar)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();previous.bindings()
    source=(previous.OUT/'effective-source.py').read_text();assert source.count(OLD)==1
    (OUT/'effective-source.py').write_text(source.replace(OLD,NEW))
    error=g.ROOT/'outputs/gr-finite-wavenumber-defined33-run.log'
    assert 'ZeroDivisionError' in error.read_text() and not (previous.OUT/'states.npz').exists()
    shutil.copy2(error,OUT/'original-error.log')
    paths=[g.ROOT/'verification/gr_finite_wavenumber_regular.py',previous.OUT/'plan.json',
        previous.OUT/'symbolic.json',OUT/'effective-source.py',OUT/'original-error.log']
    save('plan.json',dict(classification='Counterexample candidate',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        original_failure='High-precision tanh-sinh control nodes rounded onto p=Q, producing a ZeroDivisionError at an integrable logarithmic cusp, before stellar calculation.',
        change='For the independent cold reference only, subtract -b(Q)*log|p-Q| analytically, b(Q)=1/(2*sqrt(1+Q^2)). Integrate its primitive (p-Q)*(log|p-Q|-1) exactly. The remaining (b(p)-b(Q))*log|p-Q| has the continuous value zero at p=Q. No changed kernel or tolerance.'))


def bindings():
    for rel,digest in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    source=(previous.OUT/'effective-source.py').read_text()
    assert source.count(OLD)==1 and (OUT/'effective-source.py').read_text()==source.replace(OLD,NEW)
    return previous.bindings()


def module():
    ns={'__name__':'finite_response_regular'}
    exec(compile((OUT/'effective-source.py').read_text(),str(OUT/'effective-source.py'),'exec'),ns)
    ns.update(OUT=OUT,save=save,bindings=bindings)
    return ns


def run():module()['run']()
def verify():module()['verify']()


if __name__=='__main__':globals()[sys.argv[1]]()
