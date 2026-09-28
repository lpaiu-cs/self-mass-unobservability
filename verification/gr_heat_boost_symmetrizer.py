"""Carry the certified frozen heat symmetrizer through every radial Lorentz boost."""
from fractions import Fraction as F
import gzip,json,sys
import sympy as sp
import gr_heat_rank_one_damping as prior

g=prior.g;OUT=g.OUT/'gr-heat-boost-symmetrizer';sha=g.c.sha


def save(name,value):(OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def run():
    prior.verify();assert not OUT.exists();OUT.mkdir()
    files=[g.ROOT/'verification/gr_heat_boost_symmetrizer.py',prior.OUT/'manifest.json',prior.OUT/'symbolic.json',prior.OUT/'result.json',prior.OUT/'certificates.jsonl.gz']
    save('plan.json',dict(classification='Proven',checkpoint='404f5cc0',bindings={p.relative_to(g.ROOT).as_posix():sha(p) for p in files},
        target='For the same declared constant-coefficient four-variable heat model, prove strict damping at every finite nonzero real radial wavenumber in every constant radial Lorentz-boosted coordinate system. Audit every existing simple causal root and positive residue exactly.',
        assumptions='Light speed is one; time is scaled by the same positive proper tau. The four distinct characteristic speeds satisfy |c_i|<1, all rank-one residues u_i*v_i are positive, and the boost speed satisfies |V|<1. Perturbation amplitudes are expressed in the original constant characteristic basis. Any additional constant invertible change of field components is equivalent.',
        boundary='This is a coordinate boost of the frozen forced diagnostic, not a fresh physical linearization about an evolving star. No transverse directions, variable EOS/transport coefficients, geometry, reactions, metric perturbations, surface matching or nonlinear finite trajectory are supplied.'))
    save('sources.json',dict(classification='Imported from prior work',consulted='2026-09-11',scope='Primary abstract-level context only; the algebraic proof here is self-contained.',sources=[
        dict(url='https://arxiv.org/abs/0907.3906',title='Does stability of relativistic dissipative fluid dynamics imply causality?',authors='Shi Pu, Tomoi Koide, Dirk H. Rischke',
             statement='Their specified dissipative-fluid model is analysed in the rest and Lorentz-boosted frames, with a causality condition relevant to boosted stability. Their model is not substituted for ours.'),
        dict(url='https://arxiv.org/abs/2210.05067',title='Is Relativistic Hydrodynamics always Symmetric-Hyperbolic in the Linear Regime?',authors='Lorenzo Gavassino',
             statement='Connects Onsager-Casimir symmetry with symmetric hyperbolicity near equilibrium. No equilibrium premise from that abstract is imposed on our declared nonzero-heat frozen model.')]))
    V,c,d,gamma=sp.symbols('V c d gamma',real=True)
    velocity=lambda c:(c-V)/(1-V*c)
    assert sp.factor(1-velocity(c)**2-(1-V*V)*(1-c*c)/(1-V*c)**2)==0
    assert sp.factor(velocity(d)-velocity(c)-(1-V*V)*(d-c)/((1-V*c)*(1-V*d)))==0
    assert sp.factor(sp.diff(velocity(c),c)-(1-V*V)/(1-V*c)**2)==0
    speeds=sp.symbols('c:4',real=True);u=sp.Matrix(sp.symbols('u:4',nonzero=True,real=True));v=sp.Matrix(sp.symbols('v:4',nonzero=True,real=True))
    A=sp.diag(*speeds);B=u*v.T;H=sp.diag(*[v[i]/u[i] for i in range(4)])
    K0=gamma*H*(sp.eye(4)-V*A);K1=gamma*H*(A-V*sp.eye(4));D=H*B
    assert K0==K0.T and K1==K1.T and D==v*v.T
    # Lorentz coordinate differentiation and flux transform in the fixed amplitude basis.
    yt=sp.Matrix(sp.symbols('yt:4'));yx=sp.Matrix(sp.symbols('yx:4'))
    transformed=H*(gamma*(yt-V*yx)+A*gamma*(yx-V*yt))
    assert (transformed-K0*yt-K1*yx).applyfunc(sp.expand)==sp.zeros(4,1)
    # A timelike quadratic current bounds the boundary flux in every boost.
    y=sp.Matrix(sp.symbols('y:4',real=True));energy=(y.T*H*y)[0]/2;flux=(y.T*H*A*y)[0]/2
    assert sp.expand((y.T*K0*y)[0]/2-gamma*(energy-V*flux))==0
    assert sp.expand((y.T*K1*y)[0]/2-gamma*(flux-V*energy))==0
    # Derive the extra terms for varying symmetric principal coefficients.
    t,x=sp.symbols('t x');state=sp.Matrix([sp.Function(f'y{i}')(t,x) for i in range(4)])
    def symmetric_functions(prefix):
        return sp.Matrix(4,4,lambda i,j:sp.Function(f'{prefix}{min(i,j)}{max(i,j)}')(t,x))
    time=symmetric_functions('T');space=symmetric_functions('X');loss=symmetric_functions('D')
    divergence=sp.diff((state.T*time*state)[0]/2,t)+sp.diff((state.T*space*state)[0]/2,x)
    residual=(state.T*(time*sp.diff(state,t)+space*sp.diff(state,x)+loss*state))[0]
    expected=-(state.T*loss*state)[0]+(state.T*(sp.diff(time,t)+sp.diff(space,x))*state)[0]/2
    assert sp.expand(divergence-residual-expected)==0
    # Exact negative control: a superluminal speed can make boosted time energy negative.
    assert (1-V*c).subs({V:sp.Rational(3,4),c:2})==sp.Rational(-1,2)
    count=0;root_count=0;minimum_causal=F(1);minimum_separation=None;minimum_residue=None
    seen=set()
    with gzip.open(prior.OUT/'certificates.jsonl.gz','rt') as stream:
        for line in stream:
            row=json.loads(line);key=(row['cell'],row['model']);assert key not in seen;seen.add(key)
            assert row['verdict']=='stable' and len(row['roots'])==4
            end=None
            for root in row['roots']:
                a,b=map(F,root['bounds']);lo,hi=map(F,root['residue_interval'])
                assert -1<a<=b<1 and 0<lo<=hi and root['sign']==1
                minimum_causal=min(minimum_causal,1-max(abs(a),abs(b)))
                minimum_residue=lo if minimum_residue is None else min(minimum_residue,lo)
                if end is not None:
                    gap=a-end;assert gap>0
                    minimum_separation=gap if minimum_separation is None else min(minimum_separation,gap)
                end=b;root_count+=1
            count+=1
    p=prior.bindings();assert seen=={(i,j) for i in range(p['cells']) for j in range(p['models'])}
    save('result.json',dict(classification='Proven',passed=True,frozen_models=count,certified_characteristic_roots=root_count,
        minimum_causal_margin=str(minimum_causal),minimum_adjacent_root_separation=str(minimum_separation),minimum_residue_lower=str(minimum_residue),
        coordinate_equation='gamma*[(I-V*A)*y_tprime+(A-V*I)*y_xprime]+B*y=0. H=diag(v_i/u_i)>0, H*B=v*v^T>=0.',
        energy='K0=gamma*H*(I-V*A)>0 and K1=gamma*H*(A-V*I) are symmetric. The quadratic current satisfies Eprime=gamma*(E-V*F), Fprime=gamma*(F-V*E). Since |F|<=max_i|c_i|*E<E for nonzero amplitudes, every subluminal boost has positive time energy.',
        strict_damping='For a Fourier eigenpair, Re(s)*(z^*K0*z)=-|v^T*z|^2. Equality at nonzero real k would imply support only on one of the distinct transformed velocities (c_i-V)/(1-V*c_i). All v_i are nonzero, so this is impossible. Hence Re(s)<0 for every finite k!=0 and every |V|<1. At k=0 there are three conserved zero modes and one damped mode. No positive decay gap uniform in k or boosts is asserted.',
        causal='1-cprime_i^2=(1-V^2)*(1-c_i^2)/(1-V*c_i)^2>0, and the transformed root difference is (1-V^2)*(c_j-c_i)/[(1-V*c_i)*(1-V*c_j)]>0. The exact saved margins certify the premises for all11470 frozen models.',
        varying_coefficients='For K0*y_t+K1*y_x+D*y=0 with symmetric K0,K1, d_t(y^T*K0*y/2)+d_x(y^T*K1*y/2)=-y^T*sym(D)*y+y^T*(d_t K0+d_x K1)*y/2. The last term and boundary flux have no sign from the frozen residues. They cannot be dropped when judging the inhomogeneous GR star.',
        negative_control='c=2,V=3/4 gives1-V*c=-1/2. Positive rest-frame time energy alone therefore cannot justify the all-boost conclusion without the causal premise.',
        full_GR_stability_certified=False,physical_transport_calibrated=False,finite_nonlinear_trajectory_certified=False,observational_inference_complete=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    for name,key in [('plan.json','bindings'),('manifest.json','sha256')]:
        for rel,digest in json.loads((OUT/name).read_text())[key].items():assert sha(g.ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());assert r['passed'] and r['frozen_models']==11470 and r['certified_characteristic_roots']==45880
    assert min(F(r[k]) for k in ['minimum_causal_margin','minimum_adjacent_root_separation','minimum_residue_lower'])>0
    assert not r['full_GR_stability_certified'] and not r['finite_nonlinear_trajectory_certified']
    print('PASS all radial Lorentz boosts of11470 frozen heat models;45880 exact causal roots and positive residues; variable-coefficient GR boundary retained',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
