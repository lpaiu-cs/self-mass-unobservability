"""Composition derivatives of the same finite-temperature polarization self term."""
import json,sys
import numpy as np
import sympy as sp
import gr_polarization_thermodynamics as thermo

g=thermo.g;OUT=g.OUT/'gr-polarization-composition'


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();thermo.bindings()
    assert json.loads((thermo.OUT/'numerical.json').read_text())['passed']
    paths=[g.ROOT/'verification/gr_polarization_composition.py',
        g.ROOT/'verification/gr_polarization_thermodynamics.py',thermo.OUT/'plan.json',
        thermo.OUT/'states.npz',thermo.OUT/'numerical.json',thermo.density.ionic.OUT/'inputs.npz']
    save('plan.json',dict(classification='Proven',checkpoint='533da79',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        scope='Point ions at full ionization in the same ideal-electron RPA self term. Free-energy density is C_Z*g(n_e,T), C_Z=sum(n_j*Z_j^2), n_e=sum(n_j*Z_j). Ion densities are independent thermodynamic coordinates before imposing reaction or inventory constraints.',
        identities_gate=1e-12,directional_polynomial_gate=1e-11,
        controls='Exact symbolic differentiation in three independent species; Euler, mixed-temperature and rank identities. Numerical linear and quadratic directional coefficients use a separately multiplied cubic free-energy Taylor polynomial at all saved states.',
        nominal_reaction='He4+C12 -> O16: use nominal integer A and Z only to check the stoichiometric baryon/charge identities. This does not assert a physical reaction rate or its full equilibrium free energy.',
        numerical_dependency='Numerical thermodynamic primitives have passed the fixed 512/1024 quadrature comparison. Separate finite-difference and rigorous-interior acceptance is not inferred from the composition identities.',
        native_EOS_replaced=False,physical_EOS_certified=False))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    return p


def symbolic():
    n=sp.symbols('n0:3',positive=True);z=sp.symbols('z0:3',positive=True);T=sp.symbols('T',positive=True)
    e=sum(x*y for x,y in zip(n,z));C=sum(x*y*y for x,y in zip(n,z));q=sp.symbols('q',positive=True)
    fun=sp.Function('g');D=lambda k,l=0:sp.diff(fun(q,T),q,k,T,l).subs(q,e)
    free=C*fun(e,T);mu=[sp.diff(free,x) for x in n]
    for j in range(3):
        assert sp.simplify(mu[j]-z[j]**2*D(0)-C*z[j]*D(1))==0
        assert sp.simplify(sp.diff(mu[j],T)-z[j]**2*D(0,1)-C*z[j]*D(1,1))==0
        for k in range(3):
            expected=(z[j]**2*z[k]+z[j]*z[k]**2)*D(1)+C*z[j]*z[k]*D(2)
            assert sp.simplify(sp.diff(mu[j],n[k])-expected)==0
    assert sp.simplify(sum(x*m for x,m in zip(n,mu))-free-C*e*D(1))==0
    a,b,Q,d1,d2=sp.symbols('a b Q d1 d2',real=True)
    h=sp.Matrix([[2*a**3*d1+Q*a*a*d2,(a*a*b+a*b*b)*d1+Q*a*b*d2],
        [(a*a*b+a*b*b)*d1+Q*a*b*d2,2*b**3*d1+Q*b*b*d2]])
    assert sp.expand(h.det()+a*a*b*b*(a-b)**2*d1*d1)==0
    eps,dQ=sp.symbols('eps dQ',real=True)
    assert sp.diff((Q+eps*dQ)*fun(q,T),eps,2)==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        chemical='mu_j=Z_j^2*g+C_Z*Z_j*g_n. d(mu_j)/dT=Z_j^2*g_T+C_Z*Z_j*g_nT.',
        hessian='H_jk=(Z_j^2*Z_k+Z_j*Z_k^2)*g_n+C_Z*Z_j*Z_k*g_nn; it is symmetric with rank at most two.',
        pressure='sum(n_j*mu_j)-F/V=C_Z*n_e*g_n, equal to the fixed-composition pressure. The chemical caloric coefficient is mu_j-T*d(mu_j)/dT.',
        constrained='At fixed n_e,T the self free-energy density is exactly linear in C_Z. Any finite composition change conserving both sum Z_j*n_j and sum Z_j^2*n_j leaves it unchanged. Charge-conserving reaction curvature is zero and its free-energy increment is g*Delta(C_Z).',
        indefiniteness='For two distinct positive charges a,b, det(H_2x2)=-a^2*b^2*(a-b)^2*g_n^2<0 when g_n !=0. The isolated self term has an indefinite unconstrained composition Hessian. This is not instability of the complete EOS or of a fixed-charge constrained mixture.',
        sign='The positive response kernel gives S_eta>0, and the neutral eta_lnn is positive. Therefore the negative self prefactor implies g_n<0 at finite positive density/temperature in this declared model.',
        scope='An exact composition closure for one already identified Hamiltonian term, not a new correction to add on top of a screening-inclusive EOS.'))


def run():
    p=bindings();symbolic();states=dict(np.load(thermo.OUT/'states.npz'));ions=dict(np.load(thermo.density.ionic.OUT/'inputs.npz'))
    assert np.array_equal(states['cells'],ions['cells']);Z=ions['Z'];A=ions['A'];weights=ions['X']/A;weights/=weights.sum(1)[:,None]
    meanZ=weights@Z;Z2=weights@(Z*Z);assert np.allclose(Z2,states['Z2'],rtol=2e-15,atol=0)
    free,U,P,S,CV,PT,PR=states['fine_fields'];u=Z*Z;v=Z
    mu=(free/Z2)[:,None]*u+(P/meanZ)[:,None]*v
    muT=-(S/Z2)[:,None]*u+(PT/meanZ)[:,None]*v
    caloric=(U/Z2)[:,None]*u+((P-PT)/meanZ)[:,None]*v
    h_cross=P/Z2;h_charge=(PR-2*P)/meanZ
    euler=float(np.max(abs(np.sum(weights*mu,axis=1)-free-P)))
    thermal=float(np.max(abs(mu-muT-caloric)))
    # A bounded fractional composition direction stays inside the positive cone.
    direction=weights*np.linspace(-1,1,len(Z));de=direction@Z;dc=direction@u
    g0=free/Z2;g1=P/(Z2*meanZ);g2=(PR-2*P)/(Z2*meanZ**2)
    linear=Z2*g1*de+dc*g0
    quadratic=Z2*g2*de**2+2*dc*g1*de
    matrix_linear=np.sum(direction*mu,axis=1)
    matrix_quadratic=(2*h_cross*dc*de+h_charge*de**2)/meanZ
    polynomial=max(float(np.max(abs(linear-matrix_linear))),float(np.max(abs(quadratic-matrix_quadratic))))
    selected=[]
    for mass,charge in [(4,2),(12,6),(16,8)]:
        ix=np.flatnonzero((A==mass)&(Z==charge));assert len(ix)==1;selected.append(int(ix[0]))
    he,carbon,oxygen=selected;stoich=np.zeros(len(Z));stoich[[he,carbon,oxygen]]=[-1,-1,1]
    assert Z@stoich==0 and A@stoich==0 and u@stoich==24
    affinity=mu@stoich;reaction=float(np.max(abs(affinity-24*g0)))
    np.savez_compressed(OUT/'states.npz',cells=states['cells'],Z=Z,A=A,meanZ=meanZ,Z2=Z2,
        chemical_potential_kBT=mu,chemical_temperature_derivative_kB=muT,chemical_caloric_kBT=caloric,
        ne_hessian_cross=h_cross,ne_hessian_charge=h_charge,nominal_He_C_O_stoichiometry=stoich,
        nominal_He_C_O_self_affinity_kBT=affinity)
    save('result.json',dict(classification='Counterexample candidate',cells=len(meanZ),species=len(Z),
        passed=max(euler,thermal,reaction)<p['identities_gate'] and polynomial<p['directional_polynomial_gate'],
        Euler_maximum_absolute_error=euler,chemical_caloric_maximum_absolute_error=thermal,
        directional_polynomial_maximum_absolute_error=polynomial,nominal_reaction_maximum_absolute_error=reaction,
        chemical_potential_range=[float(mu.min()),float(mu.max())],
        nominal_He_C_O_self_affinity_range=[float(affinity.min()),float(affinity.max())],
        reconstructed_hessian='n_e*H = cross*(Z^2 outer Z + Z outer Z^2) + charge*(Z outer Z), in reference kBT units. No full species tensor is necessary.',
        rigorous_numerical_error_certified=False,full_physical_EOS_certified=False))
    save('manifest.json',dict(sha256={f.relative_to(g.ROOT).as_posix():g.c.sha(f) for f in OUT.iterdir() if f.is_file()}))
    verify()


def verify():
    bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());assert r['passed'] and r['cells']==3206 and r['species']==26
    assert json.loads((OUT/'symbolic.json').read_text())['passed']
    print('PASS same-self-term composition identities and 3206 x26 numerical mapping; complete EOS remains open',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
