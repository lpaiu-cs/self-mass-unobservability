"""Conditional thermodynamic boundary for replacing a resolved cell by a point."""
import json,sys
import sympy as sp
import direct_eos_gr as g

OUT=g.OUT/'gr-cell-average-identity'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def run():
    assert not OUT.exists();OUT.mkdir()
    paths=[g.ROOT/'verification/gr_cell_average_identity.py',g.OUT/'gr-spatial-preflight/manifest.json']
    save('plan.json',dict(classification='Proven',checkpoint='c8b8dad',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        scope='Fixed positive volume measure; a cell has uniform specific entropy and nuclear composition but nonuniform baryon density. Differentiable first law and strictly positive isentropic pressure derivative over its density interval are assumptions.'))
    rho,C=sp.symbols('rho C',positive=True);u=sp.Function('u')(rho);P=sp.Function('P')(rho)
    energy=rho*(C+u)
    transformed=sp.diff(energy,rho,2).subs(sp.diff(u,rho,2),sp.diff(P,rho)/rho**2-2*P/rho**3).subs(sp.diff(u,rho),P/rho**2)
    assert sp.simplify(transformed-sp.diff(P,rho)/rho)==0
    K,gamma=sp.symbols('K gamma',positive=True)
    density1,density2=sp.Rational(1),sp.Rational(3);mean=(density1+density2)/2
    # gamma=2, K=1, arbitrary affine rest-energy coefficient C.
    e=lambda r:C*r+r*r
    gap=sp.simplify((e(density1)+e(density2))/2-e(mean));assert gap==1
    specific_average=((density1**2+density2**2)/2)/mean
    assert specific_average==sp.Rational(5,2) and specific_average-mean==sp.Rational(1,2)
    save('symbolic.json',dict(classification='Proven',passed=True,
        identity='At fixed specific entropy s and composition X, du/drho=P/rho^2 gives d^2[rho*(C_X*c^2+u)]/drho^2=(1/rho)*(dP/drho)_s=Gamma1*P/rho^2.',
        Jensen='If this second derivative is positive over the cell density interval, the energy density is strictly convex. With a fixed positive normalized proper-volume measure, <epsilon(rho,s,X)> >= epsilon(<rho>,s,X), with equality only for constant rho (up to measure zero). Affine nuclear rest energy cancels from the difference.',
        uniform_state_limit='A nonuniform isentropic cell cannot in general be replaced by a single uniform state preserving proper volume, baryon mass, internal energy, composition and the original specific entropy simultaneously. At the same mean density, matching the higher averaged energy instead changes the entropy when T>0. This is a coarse-graining effect, not a physical heat source.',
        conditional_variance_bound='If 0<m<=epsilon_rhorho<=M over the entire density interval, then (m/2)*Var(rho)<= <epsilon>-epsilon(<rho>) <=(M/2)*Var(rho). This follows from Taylor remainder bounds around <rho>; the linear term averages to zero. No native continuum m,M are supplied here.',
        exact_control='For two equal-volume ideal polytropic gamma=2 states with K=1 and densities 1 and 3, mean density is 2, the energy-density gap is 1, averaged specific internal energy is 5/2, and the same-entropy uniform value is 2. Any affine rest-energy offset leaves these gaps unchanged.',
        required_reconstruction='Preserve the subcell thermodynamic profile/moments, or explicitly account for the entropy/closure error introduced by finite-volume primitive recovery. Identifying midpoint state values with cell averages without a reconstruction/error analysis is not justified.',
        limitations='No claim that all stored EOS branches satisfy the continuum convexity assumptions, no complete GR finite-volume evolution, no physical EOS certification, and no reinterpretation of numerical averaging entropy as stellar entropy production.'))
    save('manifest.json',dict(classification='Proven',sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    print('PASS isentropic cell-averaging convexity and exact control',flush=True);verify()


def verify():
    for rel,digest in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'symbolic.json').read_text())['passed']
    print('PASS cell-average identity SHA',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
