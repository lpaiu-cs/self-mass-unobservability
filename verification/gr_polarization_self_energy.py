"""All-wavenumber ideal-medium self term, with an analytic UV tail bound."""
from fractions import Fraction
import json, sys
import numpy as np
import mpmath as mp
from mpmath import iv
import sympy as sp
from scipy.integrate import quad
from interval_records import exact_endpoint, interval_text
import gr_finite_wavenumber_regular as finite

g=finite.g;OUT=g.OUT/'gr-polarization-self-energy'


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2,default=finite.previous.original.previous.scalar)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();finite.verify()
    paths=[g.ROOT/'verification/gr_polarization_self_energy.py',g.ROOT/'verification/interval_records.py',
        finite.OUT/'manifest.json',finite.OUT/'states.npz']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='5639306',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        arithmetic_sources={str(p):g.c.sha(p) for p in g.Path(mp.__file__).parent.rglob('*.py')},
        target='Compute the polarization self/volume term of the already defined finite-T finite-k RPA medium, with a rigorous all-state omitted-wavenumber tail bound and finite numerical interior comparisons.',
        z_cut=2**20,gauss_orders=[64,128],inner_tolerance=2e-11,
        outer_absolute_gate=2e-7,tail_absolute_gate='1e-9',interval_bits=128,
        control_positions=[0,801,1603,2404,3205],control_outer_tolerance=1e-9,control_gate=2e-7,
        quadrature='Set z=k/kappa_0=tan(v), integrate v from 0 to atan(z_cut). Independent outer control uses adaptive Gauss-Kronrod on the same previously verified inner response. This is a finite convergence check, not an interval enclosure of the interior integral.',
        tail='Positive S and its exact integral J=(pi^2/4)*beta*integral_0^infinity (1+beta*t)^2*f(t)dt imply tail <= J/(a*S0*z_cut^2), a=Q at k/kappa_0=1. Bound J by f<=1 below max(eta,0), f<=exp(eta-t) above it, using outward interval arithmetic.',
        global_bound='If A=J/(a*S0), positivity gives I=integral_0^infinity R(z)/(z^2+R(z))dz <= (6*A)^(1/3)/2 by a pointwise concave tangent bound. No empirical tail fit.',
        input_scope='Stored binary eta,beta,a,S0 treated as exact model parameters. Their upstream finite normalization checks are not an interval certificate for physical constants, density inversion, or the total EOS.',
        energy='C_pol/(N_i*<Z^2>*e^2*kappa_0)=-I/pi. This volume term is already part of a consistent electron-screening Hamiltonian; do not add it again to native/PC free energy.',
        physical_EOS_certified=False,native_EOS_replaced=False))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for path,digest in p['arithmetic_sources'].items():assert g.c.sha(g.Path(path))==digest,path
    finite.bindings();return p


def symbolic():
    y,z,L,R,A=sp.symbols('y z L R A',positive=True);n=sp.symbols('n',integer=True,nonnegative=True)
    ell=sp.log((y+1)/(y-1));primitive=y/2-(y*y-1)*ell/4
    assert sp.simplify(sp.diff(primitive,y)-(1-y*ell/2))==0
    assert sp.limit(primitive,y,sp.oo)==0
    assert sp.summation(1/(2*n+1)**2,(n,0,sp.oo))==sp.pi**2/8
    # Exact gap proves the tangent bound for z<L; z>=L uses R/(z^2+R)<=R/L^2.
    gap=R/L**2+(1-z/L)**2-R/(z*z+R)
    assert sp.simplify(gap-(L*z-R-z*z)**2/(L*L*(R+z*z)))==0
    objective=A/L**2+L/3
    assert sp.simplify(objective.subs(L,(6*A)**sp.Rational(1,3))-(6*A)**sp.Rational(1,3)/2)==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        kernel_integral='For fixed p>0, integral_0^infinity dQ dK(p,Q)/dp = (pi^2/4)*p*sqrt(1+p^2). Use y=Q/p, integral log|(1+y)/(1-y)|/(2y)dy=pi^2/4, and integral [1-y*log|(1+y)/(1-y)|/2]dy=0. The logarithmic cusp is integrable.',
        susceptibility_moment='Tonelli applies to the positive kernel. J=integral S(Q)dQ=(pi^2/4)*integral f(p)*p*gamma(p)dp=(pi^2*beta/4)*integral (1+beta*t)^2*f(t)dt. It is finite for every finite eta and beta>0.',
        upper_moment='Let b=max(eta,0). J <= (pi^2*beta/4)*[b+beta*b^2+beta^2*b^3/3+exp(eta-b)*((1+beta*b)^2+2*beta*(1+beta*b)+2*beta^2)]. The bound uses f<=1 on [0,b] and f<=exp(eta-t) on [b,infinity).',
        ultraviolet='With a=Q/z, R(z)=S(a*z)/S0, A=J/(a*S0), I_tail(Z)=integral_Z^infinity R/(z^2+R)dz <= A/Z^2. The inequality controls the actual infinite tail of this declared medium model, without assuming a fitted q^-4 coefficient.',
        global_bound='For any L>0: R/(z^2+R)<=R/L^2+(1-z/L)_+^2. Integrate to get I<=A/L^2+L/3, optimize L=(6*A)^(1/3): 0<I<=(6*A)^(1/3)/2. Consequently -(6*A)^(1/3)/(2*pi)<=C_pol/(N_i*<Z^2>*e^2*kappa_0)<0.',
        scope='A real all-wavenumber bound for the declared ideal-medium self term. Its numerical interior integral, thermodynamic derivatives, interacting EOS and actual GR evolution are separate obligations.'))


def tail_bounds(eta,beta,a,S0):
    plan=bindings();iv.prec=plan['interval_bits'];records=[];tails=[];globals_=[]
    for i,(ev,bv,av,sv) in enumerate(zip(eta,beta,a,S0)):
        e,b,aa,ss=[iv.mpf(float(v)) for v in [ev,bv,av,sv]];cut=iv.mpf(float(max(ev,0)))
        j=iv.pi**2*b/4*(cut+b*cut**2+b*b*cut**3/3+iv.exp(e-cut)*((1+b*cut)**2+2*b*(1+b*cut)+2*b*b))
        total=j/(aa*ss);tail=total/plan['z_cut']**2;full=iv.exp(iv.ln(6*total)/3)/2
        tails.append(exact_endpoint(tail._mpi_[1]));globals_.append(float(exact_endpoint(full._mpi_[1])))
        records.append(dict(position=i,J_upper_expression=interval_text(j),A_upper_expression=interval_text(total),tail_upper_expression=interval_text(tail),global_I_upper_expression=interval_text(full)))
    passed=max(tails)<Fraction(plan['tail_absolute_gate'])
    save('tail-bounds.json',dict(classification='Proven',passed=passed,records=records,
        maximum_tail_upper=str(max(tails)),maximum_tail_upper_float=float(max(tails)),
        scope='All z>=z_cut, exact saved input parameters, outward interval arithmetic. The lower endpoint of an upper-bound expression is NOT a lower bound on the unknown tail or full I.'))
    assert passed
    return np.array([float(x) for x in tails]),np.array(globals_)


def run():
    plan=bindings();symbolic();model=finite.module();a0=dict(np.load(finite.OUT/'states.npz'))
    eta,beta,S0=[a0[k] for k in ['eta','beta','S0']]
    j=int(np.flatnonzero(a0['k_over_kappa']==1)[0]);a=a0['Q'][:,j]
    tail,upper=tail_bounds(eta,beta,a,S0);end=np.arctan(float(plan['z_cut']))
    orders=[];positive=[]
    for order in plan['gauss_orders']:
        nodes,weights=np.polynomial.legendre.leggauss(order);v=(nodes+1)*end/2;weights=weights*end/2
        samples=[]
        for k,angle in enumerate(v):
            z=np.tan(angle);r=model['response'](eta,beta,a*z,S0,plan['inner_tolerance'])
            assert np.all(r>0)
            samples.append(r*(1+z*z)/(z*z+r))
            if (k+1)%16==0:print('SELF ENERGY',order,k+1,flush=True)
        value=weights@np.array(samples);orders.append(value)
        # Positive control: the same mapping/quadrature integrates the TF constant exactly.
        tf=float(weights@np.ones(order));positive.append(dict(order=order,TF_integral=tf,exact=end,difference=abs(tf-end),passed=abs(tf-end)<1e-13))
        np.savez_compressed(OUT/f'order-{order}.npz',v=v,weights=weights,integrand=samples,value=value)
    coarse,fine_value=orders;controls=[]
    for i in plan['control_positions']:
        def f(angle):
            z=np.tan(angle)
            r=float(model['response'](eta[i:i+1],beta[i:i+1],a[i:i+1]*z,S0[i:i+1],plan['inner_tolerance'])[0])
            return r*(1+z*z)/(z*z+r)
        value,error=quad(f,0,end,epsabs=plan['control_outer_tolerance'],epsrel=plan['control_outer_tolerance'],limit=250)
        delta=abs(value-fine_value[i]);controls.append(dict(cell=int(a0['cells'][i]),reference=value,quad_estimate=error,gauss=float(fine_value[i]),difference=float(delta),passed=delta<plan['control_gate']))
        print('SELF CONTROL',int(a0['cells'][i]),float(delta),flush=True)
    score=float(np.max(abs(coarse-fine_value)))
    np.savez_compressed(OUT/'states.npz',cells=a0['cells'],I_coarse=coarse,I_fine=fine_value,tail_upper=tail,global_I_upper=upper,
        self_coefficient=-fine_value/np.pi,magnitude_over_TF=2*fine_value/np.pi)
    passed=score<plan['outer_absolute_gate'] and all(x['passed'] for x in controls+positive) and np.all(fine_value>0) and np.all(fine_value<upper)
    save('result.json',dict(classification='Counterexample candidate',passed=passed,cells=len(eta),finite_outer_difference=score,
        independent_outer_controls=controls,TF_positive_controls=positive,
        maximum_tail_upper=float(tail.max()),I_range=[float(fine_value.min()),float(fine_value.max())],
        self_coefficient_range=[float((-fine_value/np.pi).min()),float((-fine_value/np.pi).max())],
        magnitude_over_TF_range=[float((2*fine_value/np.pi).min()),float((2*fine_value/np.pi).max())],
        physical_EOS_certified=False,native_EOS_replaced=False,
        scope='Finite numerical interior with separately rigorous omitted-tail and global bounds for exact declared input parameters. No rigorous numerical interior enclosure, thermodynamic derivative, full EOS correction or observation inference.'))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify()


def verify():
    p=bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());a=dict(np.load(OUT/'states.npz'))
    assert r['passed'] and r['cells']==len(a['cells'])==3206
    assert np.all(a['I_fine']>0) and np.all(a['I_fine']<a['global_I_upper'])
    for name in ['symbolic.json','tail-bounds.json']:assert json.loads((OUT/name).read_text())['passed']
    print('PASS all-state ideal-medium self integral finite checks and rigorous UV-tail bounds; full EOS remains open',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
