"""Recompute the exported coefficient and kernel budgets using exact rationals."""
from fractions import Fraction as F
import json, sys
import gr_plasma_moment_certificate as certificate

g=certificate.g;OUT=g.OUT/'gr-plasma-moment-audit'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')
def ends(s):
    a,b=map(F,s[1:-1].split(','));assert a<=b;return a,b


def prepare():
    assert not OUT.exists();OUT.mkdir();certificate.verify()
    paths=[g.ROOT/'verification/verify_plasma_moment_certificate.py',certificate.OUT/'manifest.json']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='d5246ee',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        rule='Use outward saved decimal moment-error upper bounds and denominator lower bound as exact Fractions. Recompute coefficient/kernel budgets, including the saved certified tail upper bounds. Preserve the original arithmetic outputs and record any additional outward export padding.',
        threshold='1e-10',physical_EOS_certified=False))


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    cert=json.loads((certificate.OUT/'uniform-certificate.json').read_text())
    result=json.loads((certificate.OUT/'result.json').read_text())
    D=ends(cert['lower_normalized_Ip'])[0]
    errors={r['component']:ends(r['common_rule_error'])[1] for r in cert['rows'] if r['eta_order']==r['beta_order']==0}
    ED=errors[0]+errors[1]/3;assert D>ED
    coeff=[F(0)]+[(errors[j]/(2*j+1)+errors[j+1]/(2*j+3)+F(3,2)*ED/(2*j+1))/(D-ED) for j in range(1,25)]
    maximum_padding=F(0)
    for value,stored in zip(coeff,result['coefficient_errors']):
        maximum_padding=max(maximum_padding,value-ends(stored)[1])
    T=sum(coeff)+ends(result['transverse_uniform_tail'])[1]
    L=sum((2*j+1)*v for j,v in enumerate(coeff))+ends(result['longitudinal_uniform_tail'])[1]
    maximum_padding=max(maximum_padding,T-ends(result['transverse_uniform_kernel_error'])[1],L-ends(result['longitudinal_uniform_kernel_error'])[1])
    assert max(T,L)<F(plan['threshold']) and maximum_padding<F('1e-55')
    assert len(cert['rows'])==157
    assert all(ends(r['common_rule_error'])[1]<=F('1e-14') for r in cert['rows'])
    save('result.json',dict(classification='Proven',passed=True,components=157,
        denominator_lower_rational=str(D),denominator_error_upper_rational=str(ED),
        transverse_kernel_upper_rational=str(T),longitudinal_kernel_upper_rational=str(L),
        maximum_additional_export_padding_rational=str(maximum_padding),
        display_only=dict(transverse=float(T),longitudinal=float(L),export_padding=float(maximum_padding)),
        scope='Independent exact rational propagation of the saved moment/tail certificates. Original directed-interval calculations and their analytic assumptions remain binding. This does not certify the previous adaptive native binary or physical EOS.'))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    certificate.verify()
    for name,key in [('plan.json','bindings'),('manifest.json','sha256')]:
        for rel,digest in json.loads((OUT/name).read_text())[key].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'result.json').read_text())['passed']
    print('PASS exact-rational plasma coefficient/kernel error propagation',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
