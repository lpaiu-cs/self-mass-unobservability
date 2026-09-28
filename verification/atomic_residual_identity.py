"""Symbolic reaction/Fisher identities and a gradient-shift positive control."""
from fractions import Fraction as F
import json,sys
import sympy as sp
import direct_eos_gr as g

OUT=g.OUT/'atomic-residual-identities'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def run():
    OUT.mkdir(exist_ok=True);assert not (OUT/'result.json').exists()
    x,y,N=sp.symbols('x y N',positive=True);c0,c1,delta=sp.symbols('c0 c1 delta',real=True)
    phi=sp.Function('phi')(x,y)
    energy=x*(sp.log(x)-1)+y*(sp.log(y)-1)+c0*x+c1*y+phi
    reaction=sp.log(y/x)+c1-c0+sp.diff(phi,y)-sp.diff(phi,x)
    assert sp.simplify(sp.diff(energy,y)-sp.diff(energy,x)-reaction)==0
    perturb=energy+delta*y
    assert sp.simplify(sp.diff(perturb,y)-sp.diff(perturb,x)-reaction)==delta
    assert sp.hessian(perturb,[x,y])==sp.hessian(energy,[x,y])
    n=sp.symbols('n0:3',positive=True);r=sp.symbols('r0:3',real=True);lam=sp.symbols('lam',real=True)
    norm=sum(v*(w-lam)**2 for v,w in zip(n,r));mean=sum(v*w for v,w in zip(n,r))/sum(n)
    assert sp.simplify(sp.diff(norm,lam).subs(lam,mean))==0
    assert sp.expand(sp.diff(norm,lam,2)-2*sum(n))==0
    assert sp.factor(norm.subs(lam,0)-norm.subs(lam,mean))==sp.factor(sum(v*w for v,w in zip(n,r))**2/sum(n))
    path=g.OUT/'atomic-stationarity-defined/frozen-residual-bounds.json'
    bound=json.loads(path.read_text());rmax=F(bound['maximum_absolute_atomic_residual_upper'])
    injected=F(1,10000);threshold=F(2,10**8)
    assert injected-rmax>threshold
    save('result.json',dict(classification='Proven',passed=True,
        reaction_identity='The derivative along atomic transfer x->y equals ln(y/x)+(c1-c0)+partial_y(phi)-partial_x(phi).',
        Fisher_identity='Subtracting the weighted mean minimizes sum_i nu_i*(r_i-lambda)^2; fixing the reference residual to zero gives a valid upper bound.',
        positive_control=dict(gradient_shift=str(injected),certified_shifted_residual_lower=str(injected-rmax),gate=str(threshold),detected=True,
            method='Adding delta times the receiving species population shifts its reaction residual by delta and leaves the Hessian unchanged. Even the worst frozen observed residual cannot cancel this registered-scale perturbation.'),
        scope='Algebraic identities and a manufactured additive-gradient control. This is not an injected native-library run or a bound on native/physical EOS discrepancy.'))
    paths=[p for p in OUT.iterdir() if p.is_file()]+[g.ROOT/'verification/atomic_residual_identity.py',path]
    save('manifest.json',dict(classification='Proven',sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths}))
    print('PASS ATOMIC REACTION/FISHER IDENTITIES AND GRADIENT POSITIVE CONTROL',flush=True)


def verify():
    manifest=json.loads((OUT/'manifest.json').read_text())
    for rel,digest in manifest['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    print('PASS ATOMIC RESIDUAL IDENTITIES',len(manifest['sha256']),'artifact SHA',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
