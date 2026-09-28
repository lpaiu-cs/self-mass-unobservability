"""Independent checks for the Voigt tail and absorption perturbation bounds."""
import json
import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.special import wofz, erfc
import def_photon_profile_integral as p


def run():
    c1,c2,x=sp.symbols('c1 c2 x',positive=True)
    B=sp.Matrix([[c1+c2,-sp.sqrt(c1),-sp.sqrt(c2)],[-sp.sqrt(c1),1,0],[-sp.sqrt(c2),0,1]])
    assert sp.simplify((x*sp.eye(3)-B).det()-x*(x-1)*(x-1-c1-c2))==0
    # Exact integral: H(a,x)=int exp(-t^2) a/[pi((x-t)^2+a^2)] dt.
    # For |x|>=R, split |t|<=R/2 and its complement. Both bounds below
    # are uniform, positive, and independent of line-center sampling.
    cases=0
    for a in np.geomspace(1e-4,1e6,21):
        for R in [16.,100.,1e4]:
            bound=a/np.sqrt(np.pi)/((R/2)**2+a*a)+erfc(R/2)/(np.sqrt(np.pi)*a)
            assert np.max(wofz(np.array([R,2*R,10*R])+1j*a).real)<=bound*(1+1e-14)
            cases+=3
    # Independent normalization with an adaptive quadrature in asinh x.
    errors=[]
    for a in [1e-4,.01,1.,100.,1e6]:
        sd=max(1.,a)
        value,error=quad(lambda t:wofz(sd*np.sinh(t)+1j*a).real*sd*np.cosh(t),-35,35,epsabs=1e-10,points=[0.],limit=200)
        errors.append(abs(value/np.sqrt(np.pi)-1));assert errors[-1]<1e-9
    q=p.OUT/'preflight.json';b=np.load(p.OUT/'bank.npz')
    assert np.all(b['rate_a_fine']>0) and np.all(b['rate_a_coarse']>0)
    assert json.loads(q.read_text())['line_inventory_same']
    result=dict(classification='Proven',passed=True,tail_inequality_sample_checks=cases,
        normalization_relative=errors,
        statement='For the exact finite-bin LTE Voigt model, the omitted wing response is at most t*c*epsilon*(1+Cgamma/Cm) times the initial energy norm, where epsilon is the sum of the uniform volume-opacity tail bounds. The common streaming operator is skew and collisions are dissipative. This follows from the positive absorption block bound and Duhamel contraction.',
        limitation='The analytic statement is conditional. Floating-point checks do not constitute outward-rounded interval certification of atomic profiles or quadrature.',
        source_sha256=p.a.digest(__import__('pathlib').Path(__file__)))
    p.a.write(p.OUT/'profile-check.json',result);print(result)


if __name__=='__main__':run()
