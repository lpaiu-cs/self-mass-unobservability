"""Supply the missing dielectric constant from its bound source, outwardly."""
import inspect,json,shutil,sys
from fractions import Fraction as F
import numpy as np
from mpmath import iv
import gr_response_complex_domain as g
import gr_polarization_thermodynamics as thermo

ORIGINAL=g.OUT;OUT=ORIGINAL.parent/'gr-response-complex-defined'
OLD="stored=dict(np.load(cusp.THERMO/'states.npz'));B=I(float(stored['B']))"
NEW="B=cusp.interval(json.loads((OUT/'B-binding.json').read_text())['interval'])"


def prepare():
    p=g.bindings();assert not OUT.exists();OUT.mkdir()
    error=g.ROOT/'outputs/gr-response-complex-domain33-run.log';assert "KeyError: 'B'" in error.read_text() and not (ORIGINAL/'result.json').exists()
    shutil.copy2(error,OUT/'original-error.log');iv.prec=p['bits']
    constants=thermo.density.ionic.OUT/'constants.json';alpha=json.loads(constants.read_text())['native_binary64']['alpha']
    exact=g.I(alpha)/iv.pi;candidate=float(alpha/np.pi);lo=min(g.low(exact),F.from_float(candidate));hi=max(g.high(exact),F.from_float(candidate))
    (OUT/'B-binding.json').write_text(json.dumps(dict(classification='Proven',interval=f'[{lo},{hi}]',
        alpha_binary64=alpha,candidate_B_binary64=candidate,exact_alpha_over_pi=g.cusp.interval_text(exact),
        scope='The interval contains both exact declared alpha divided by mathematical pi and the binary64 constant used by the existing numerical candidate. All zero-free statements hold uniformly throughout this interval. No equivalence of the two constants is assumed.'),indent=2)+'\n')
    source=(g.ROOT/'verification/gr_response_complex_domain.py').read_text();assert source.count(OLD)==1
    (OUT/'effective-source.py').write_text(source.replace(OLD,NEW))
    p['bindings'].update({path.relative_to(g.ROOT).as_posix():g.sha(path) for path in [constants,g.ROOT/'verification/gr_response_complex_runner.py',ORIGINAL/'plan.json',OUT/'original-error.log',OUT/'B-binding.json',OUT/'effective-source.py']})
    p['execution_overlay']='The saved state NPZ omits B; the original run stopped before computing any state row. Read the declared alpha source and enclose both exact alpha/pi and the existing numerical binary B. Preserve all analytic formulas, state inputs, radius rules and thresholds.'
    (OUT/'plan.json').write_text(json.dumps(p,ensure_ascii=False,indent=2)+'\n')


def activate():
    if g.OUT==OUT:return
    g.OUT=OUT;source=inspect.getsource(g.run);assert source.count(OLD)==1
    exec(compile(source.replace(OLD,NEW),str(OUT/'effective-source.py'),'exec'),g.__dict__)


def run():activate();g.run()
def verify():activate();g.verify()


if __name__=='__main__':globals()[sys.argv[1]]()
