"""Proven conditional GR heat-current constraints and actual initial-data candidate.

The existing fixed-geometry transport operators remain controls. This supplies
the missing initial extrinsic curvature, not a solved full fluid/metric path.
"""
import json,sys,urllib.request
import numpy as np
import sympy as sp
import gr_microphysics as micro

g=micro.g;OUT=g.OUT/'gr-heat-initial-constraints'
URL='https://arxiv.org/pdf/gr-qc/0703035'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    with urllib.request.urlopen(URL,timeout=60) as response:(OUT/'gourgoulhon-0703035.pdf').write_bytes(response.read())
    assert (OUT/'gourgoulhon-0703035.pdf').read_bytes().startswith(b'%PDF')
    paths=[g.ROOT/'verification/gr_heat_initial_constraints.py',g.OUT/'initial-state-17-4.npz',
        g.OUT/'gr-microphysics/auxiliaries.npz',g.OUT/'gr-transport/diagnostics.npz',OUT/'gourgoulhon-0703035.pdf']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='a4fb34c',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        source=dict(classification='Imported from prior work',url=URL,title='3+1 Formalism and Bases of Numerical Relativity',
            author='Eric Gourgoulhon',date='2007-03-06',equations=['4.93','4.94','5.12','5.50'],
            inspected_pdf_pages_zero_based=[64,72,78]),
        convention='Signature -+++, K_ij=-(1/(2N))*partial_(ct) gamma_ij, shift zero. Polar areal 3-metric diag(a^2,r^2,r^2 sin^2 theta), K^i_j=diag(A,0,0). Initial material 4-velocity equals slice normal, but a radial heat current is nonzero.',
        material='Retain the new EOS rho_B,T,X and the supplied isotropic pressure; total energy includes nuclear rest energy. The initial interior heat current is the same actual EOS/opacity diffusion diagnostic, with both boundary heat currents zero.',
        scope='Conditional continuum identities plus a finite initial-data candidate. Hamiltonian discretization, atmosphere, diffusion validity, causal transport closure and subsequent fluid/metric evolution are not certified.',
        algebraic_relative_tolerance=1e-12,physical_EOS_certified=False,full_GR_evolution=False))
    symbolic()


def symbolic():
    r,theta,phi=sp.symbols('r theta phi',positive=True);coords=[r,theta,phi]
    a=sp.Function('a')(r);A=sp.Function('A')(r);B=sp.Function('B')(r)
    metric=sp.diag(a*a,r*r,r*r*sp.sin(theta)**2);inverse=metric.inv()
    Gamma=[[[sp.simplify(sum(inverse[k,l]*(sp.diff(metric[l,j],coords[i])+sp.diff(metric[l,i],coords[j])-sp.diff(metric[i,j],coords[l])) for l in range(3))/2)
        for j in range(3)] for i in range(3)] for k in range(3)]
    Ricci=sp.zeros(3)
    for i in range(3):
        for j in range(3):
            Ricci[i,j]=sp.simplify(sum(sp.diff(Gamma[k][i][j],coords[k])-sp.diff(Gamma[k][i][k],coords[j])+
                sum(Gamma[k][i][j]*Gamma[l][k][l]-Gamma[l][i][k]*Gamma[k][j][l] for l in range(3)) for k in range(3)))
    curvature=sp.simplify(sp.trace(inverse*Ricci));mass=r*(1-a**-2)/2
    assert sp.simplify(curvature-4*sp.diff(mass,r)/r**2)==0
    K=sp.diag(A,B,B);S=K-sp.eye(3)*sp.trace(K);momentum=[]
    for i in range(3):
        value=sum(sp.diff(S[j,i],coords[j])+sum(Gamma[j][j][k]*S[k,i]-Gamma[k][j][i]*S[j,k] for k in range(3)) for j in range(3))
        momentum.append(sp.simplify(value))
    assert sp.simplify(momentum[0]-2*(A-B)/r+2*sp.diff(B,r))==0
    assert momentum[1:]==[0,0]
    assert sp.expand(sp.trace(K)**2-sp.trace(K*K)-(4*A*B+2*B*B))==0
    q=sp.Function('q')(r);N=sp.Function('N')(r);E,P=sp.symbols('E P',real=True)
    polar_A=4*sp.pi*r*a*q
    assert sp.simplify(momentum[0].subs(B,0).doit().subs(A,polar_A)-8*sp.pi*a*q)==0
    # Differentiate the mass evolution law and compare the independently
    # projected matter energy equation, then use Hamiltonian + polar lapse.
    mdot=-4*sp.pi*r*r*N*q/a
    Edot=N*polar_A*(E+P)-N*sp.diff(r*r*q,r)/(a*r*r)-2*q*sp.diff(N,r)/a
    lapse_derivative=N*(4*sp.pi*r*a*a*(E+P)-sp.diff(a,r)/a)
    assert sp.simplify((sp.diff(mdot,r)-4*sp.pi*r*r*Edot).subs(sp.diff(N,r),lapse_derivative))==0
    assert sp.simplify(-N*a*polar_A/a+N*polar_A)==0
    rho,C,u,p,heat=sp.symbols('rho C u p heat',real=True)
    specific=N*polar_A*p/rho+heat
    assert sp.simplify(rho*specific+(C+u)*rho*N*polar_A-N*polar_A*(rho*(C+u)+p)-rho*heat)==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        spatial_curvature='R3=4*m_prime/r^2, m=r*(1-a^-2)/2.',
        general_spherical_momentum='For K^i_j=diag(A,B,B), D_j(K^j_r-delta^j_r K)=2*(A-B)/r-2*B_prime.',
        Hamiltonian_extrinsic_term='K^2-K^i_j*K^j_i=4*A*B+2*B^2; identically zero at B=0.',
        polar_initial_solution='For the proper-frame geometrized heat current q and initial material velocity zero, j_r=a*q. Polar B=0 gives A=4*pi*r*a*q. K=0 is inconsistent with any nonzero j_r. Adding this A leaves the original Hamiltonian equation unchanged.',
        metric_and_mass_rates='partial_(ct) ln a=-N*A; partial_(ct) m=-4*pi*r^2*N*q/a. The outer vacuum-normalized N*a=1 converts this to the usual luminosity-at-infinity mass-loss relation.',
        baryon_and_specific_energy='At initial normal material velocity zero, partial_t ln rho_B=c*N*A; partial_t u=(P/rho_B)*partial_t ln rho_B -(1/N)*partial_mu L_infinity, where mu is outward-increasing baryon mass and physical time is in seconds.',
        initial_Hamiltonian_propagation='partial_r(partial_(ct) m)=4*pi*r^2*partial_(ct) E follows from the matter energy equation and N_prime/N+a_prime/a=4*pi*r*a^2*(E+P). The latter uses Hamiltonian plus polar lapse. This is an exact initial continuum identity, not a discrete evolution certificate.',
        missing='A transport evolution law, material acceleration, metric/lapse evolution after the initial slice, boundary atmosphere and the existing physical EOS uncertainties remain. These identities do not keep the material velocity zero for a finite interval.'))
    print('PASS spherical heat-current constraints and initial propagation identities',flush=True)


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    state,aux=micro.inputs();aEOS=aux['eos'];c=g.c.gr.C*100;G=g.c.gr.G*1000
    rho=np.exp(state['lnd']);P=aEOS[:,1];CX=(state['X']/g.c.A)@g.c.W
    energy=rho*(CX*c*c+aEOS[:,2]);enthalpy=energy+P
    diagnostics=dict(np.load(g.OUT/'gr-transport/diagnostics.npz'))
    faces=np.r_[0.,diagnostics['interior_Linf'],0.]
    Lmid=(faces[:-1]+faces[1:])/2;r=state['r_mid_m']*100;m=state['m_mid_geom']*100;N=np.exp(state['nu'])
    radial=(1-2*m/r)**-.5;proper_heat=Lmid/(4*np.pi*r*r*N*N);q=G*proper_heat/c**5
    j_r=radial*q;A=4*np.pi*r*j_r
    momentum=2*A/r-8*np.pi*j_r;nonzero=j_r!=0
    relative=float((abs(momentum[nonzero])/(8*np.pi*abs(j_r[nonzero]))).max())
    assert relative<plan['algebraic_relative_tolerance']
    dlnrho=c*N*A;dlnradial=-dlnrho
    mdot=-4*np.pi*r*r*N*q/radial*c
    heating=(faces[1:]-faces[:-1])/(state['dm']*N)
    compression=P/rho*dlnrho
    dlnT=(heating+(P/rho-aEOS[:,9])*dlnrho)/aEOS[:,10]
    old_dlnT=heating/aEOS[:,10]
    # Algebraic dominant-energy check for the specified isotropic-stress,
    # radial-heat tensor. This is not a causal transport/evolution theorem.
    Q=proper_heat/c;discriminant=enthalpy**2-4*Q**2;assert np.all(discriminant>0)
    root=np.sqrt(discriminant);landau_energy=(energy-P+root)/2;landau_pressure=(-energy+P+root)/2
    dominant=np.minimum(landau_energy-abs(landau_pressure),landau_energy-abs(P))
    assert np.all(dominant>=0)
    beta=2*Q/(enthalpy+root);assert np.all(abs(beta)<1)
    np.savez_compressed(OUT/'initial-data.npz',r_cm=r,m_geom_cm=m,lapse=N,radial_metric_factor=radial,
        proper_heat_flux_cgs=proper_heat,normal_momentum_covariant_geom=j_r,
        radial_extrinsic_curvature_per_cm=A,dlnrho_dt=dlnrho,dlnradial_metric_dt=dlnradial,
        dmass_geom_cm_dt=mdot,compression_specific_energy_rate=compression,
        heat_specific_energy_rate=heating,dlnT_dt=dlnT,fixed_metric_dlnT_dt=old_dlnT,
        energy_frame_velocity_over_c=beta,momentum_constraint_residual=momentum)
    result=dict(classification='Counterexample candidate',completed=True,cells=len(r),
        original_state_sha256=g.c.sha(g.OUT/'initial-state-17-4.npz'),
        nonzero_normal_heat_current_cells=int(nonzero.sum()),
        original_zero_extrinsic_curvature_momentum_relative_defect=1.,
        new_pointwise_momentum_relative_residual=relative,
        added_Hamiltonian_extrinsic_term_identically_zero=True,
        maximum_abs_radial_extrinsic_curvature_per_cm=float(abs(A).max()),
        maximum_abs_dlnrho_dt_per_second=float(abs(dlnrho).max()),
        maximum_abs_mass_geom_rate_cm_per_second=float(abs(mdot).max()),
        maximum_abs_compression_energy_rate_erg_g_s=float(abs(compression).max()),
        maximum_abs_added_dlnT_dt=float(abs(dlnT-old_dlnT).max()),
        maximum_energy_frame_velocity_over_c=float(abs(beta).max()),
        algebraic_dominant_energy_condition_at_saved_points=True,
        old_spatial_Hamiltonian_error_recertified=False,continuous_constraint_error_certified=False,
        atmosphere_solved=False,transport_causality_certified=False,physical_EOS_certified=False,full_GR_evolution=False)
    save('result.json',result);print('HEAT CURRENT INITIAL GR DATA',result,flush=True);verify()


def verify():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'result.json').read_text())['completed']
    path=OUT/'manifest.json'
    if not path.exists():save('manifest.json',dict(classification='Counterexample candidate',
        sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    for rel,digest in json.loads(path.read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    print('PASS HEAT CURRENT INITIAL CONSTRAINT SHA',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
