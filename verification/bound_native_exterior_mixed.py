"""Counterexample candidate: bound the missing direct-photon mixed source.

Reuse positive saved emission and primary plus finite-Born norms. The generated
background scalar, signed later emission returns and full physical errors remain
outside this sector. No fluid, EOS-root or ray history is replayed.
"""
from pathlib import Path
import json, resource, sys, time
import numpy as np
import mpmath as mp
import sympy as sp
import bound_native_exterior_incident as prior
import read_native_matched_mass as matched

run=prior.run
OUT=Path('native-exterior-mixed164-work')
read,write,sha=run.read,run.write,run.sha
CAPS=dict(prepare=30,apply=60,audit=45)
TOTAL=sum(CAPS.values())


def symbolic():
    r,R,p,e,g,rho=sp.symbols('r R p e g rho',positive=True)
    b,a,Phi,li,zi,zb,zbr=sp.symbols('b a Phi li zi zb zbr',nonzero=True)
    # Test the distributional primitive against a decaying, nonconstant test
    # function. This checks both the moving delta sign and its smooth weight.
    z=r-R;v=sp.exp(-r)
    source=e*(sp.integrate(v,(r,R,p))+(g-(p-R))*v.subs(r,p)+rho*sp.diff(v,r).subs(r,p))
    primitive=-e*(sp.integrate(z*sp.diff(v,r),(r,R,p))+sp.integrate(g*sp.diff(v,r),(r,p,sp.oo)))
    primitive+=e*rho*sp.diff(v,r).subs(r,p)
    assert sp.simplify(source-primitive)==0
    radial=4*(-e/(a*r))*li/(r*b)
    assert sp.simplify(r*a*Phi*radial+4*e*Phi*li/(r*b))==0
    gradient=2*a*r*Phi*(zb*2*li/(r*b)+zi*zbr)
    assert sp.simplify(gradient.subs(zbr,-2*e/(a*b*r*r))-(4*a*Phi*li*zb/b-4*e*Phi*zi/(b*r)))==0
    eps,mu,dm,zet=sp.symbols('eps mu dm zet')
    Q=lambda x:(1-x*x)/x
    exact=sp.exp(eps*(g-2*zet))*Q(mu+eps*dm)
    first=Q(mu)*(g-2*zet)-(1+mu*mu)*dm/(mu*mu)
    assert sp.simplify(sp.diff(exact,eps).subs(eps,0)-first)==0
    # The static-gradient shell contributes 2*zeta times the old shell:
    # together its amplitude and the two measure factors cancel the zeta.
    assert sp.simplify(first+2*zet*Q(mu)-(g*Q(mu)-(1+mu*mu)*dm/(mu*mu)))==0
    z=sp.symbols('z',positive=True)
    F=-sp.sqrt(1-z)/sp.sqrt(1+z)
    assert sp.simplify(sp.diff(F,z)-1/(sp.sqrt(1-z)*(1+z)**sp.Rational(3,2)))==0
    tails={n:sp.integrate((r-R)/r**n,(r,R,sp.oo)) for n in [3,4]}
    assert tails[3]==1/(2*R) and tails[4]==1/(6*R*R)
    return dict(classification='Proven',passed=True,
        constraints='C_prime=a*r*q_lambda=e*[zeta_I_prime*H(rp-r)+(g-zeta_I(rp))*delta(r-rp)-rho*delta_prime(r-rp)]. Thus C=e*[(zeta_I(r)-zeta_I(R))*H(rp-r)+(g-zeta_I(R))*H(r-rp)-rho*delta(r-rp)], with zero additional inner port.',
        scalar_measure='The mass-source measure is 2*Phi*C/(r*b) dr. The smooth known constraint stress is -4*e*Phi*lambda_I/(r*b)*H(rp-r) dr.',
        static_gradient='Also include 2*a*r*Phi*(zeta_B*zeta_I_prime+zeta_I*zeta_B_prime) dr. Its smooth part is 4*a*Phi*lambda_I*zeta_B/b-4*e*Phi*zeta_I/(b*r); its shell is -2*e*Phi*zeta_I*(1-mu^2)*delta dr.',
        shell_pullback='At fixed radius, the ray time Jacobian adds a second -zeta_I. Together with the static-gradient shell it cancels: delta shell integrand is -e*Phi/(2*a)*[g*(1-mu^2)/mu-(1+mu^2)*delta_mu/mu^2] H dr, plus the arrival-front shift. delta_mu=(1-mu^2)*h/mu and |h|<=Hmax*(1-R/r).',
        uniform_domain='On t>=0,t-x/c<=T each source radius has time width at most T+x(r)/c. This also covers the noncausal-in-radius constraint tail. Reuse the old c*T potential norm only after adding its positive x-weighted vacuum integral.',
        scope='Distribution and first-variation identities in the declared sector, plus elementary tail integrals. No assertion of exact ADM conservation or physical EOS closure.')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    files=[Path(__file__),Path(prior.__file__),Path(run.__file__),Path(matched.__file__),
        Path(matched.previous.mixed.__file__),run.OUT/'bound-plan.json',run.OUT/'bound-result.json',
        run.OUT/'audit.json',run.OUT/'result.json',matched.OUT/'final-result.json',matched.OUT/'compensated-charge-parts.npz',
        matched.OUT/'charge-parts.npz',run.previous.BACKGROUND,
        run.previous.ROOT/'def-native-global-scalar-closure/bound.json']
    # Binding old plans is insufficient if their physical input files drift.
    for name in ['plan.json','bound-plan.json']:
        for path,h in read(run.OUT/name)['bindings'].items():
            assert sha(path)==h,path
            files.append(Path(path))
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='53a1718d7ec4b7c7bf4ff01152e83dab882d8ec9',
        claim='Enclose the direct emitted-photon mixed constraints, static-gradient cross term, changed ray scalar readout, and their additional homogeneous mass normalization. Apply the envelope to the Phase163 selected charge.',
        decision='Can these omitted terms erase the selected endpoint? A selected exterior envelope below1percent avoids a new ray/fluid history for this sector.',
        family='The original nonnegative64/128 accepted-port histories and every nonnegative emission satisfying their energy/luminosity caps, with the prescribed primary plus continuous finite-source first Born on the exact initial scalar vacuum. No initially exterior photons, fixed emission inventory and launch angles.',
        exclusions=['Later signed background emission corrections and generated background scalar, including its reciprocal constraints and metric action.','Returned material fields in the external input metric, additional inner/interface ports, exact initial ADM/momentum/flux conservation.','Full direct-field, EOS, uniform derivative, source interpolation, spatial, boundary, nonlinear and same-inventory static/observational errors.'],
        method='50-digit outward intervals on existing binary upper bounds. Integrate positive source envelopes over the entire observer strip. Separate packet measure, radial displacement, angular and arrival fronts; include the static-gradient cross term. Bound the infinity mass port instead of assuming it zero.',
        gates=dict(selected_envelope_fraction=.01,potential_contraction=1.,independent_relative=1e-10),
        budget=dict(actions=CAPS,total_seconds=TOTAL,CPU_threads=1,virtual_GiB=3,new_fluid_steps=0,new_EOS_roots=0,new_rays=0),
        measured_basis='Phase159 interval body4.96s, Phase163 all four action bodies6.63s. One existing-background load and scalar arithmetic;60s apply and135s total include a conservative margin. This is not a prediction for a new physical trajectory.',
        stop='No automatic history, horizon, grid or iteration enlargement. Preserve the Phase159 point failure; an inconclusive bound remains inconclusive.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)}))
    write(OUT/'symbolic.json',symbolic())


def terms(p,I,sqrt,pi):
    """The same explicit inequalities, usable by interval and scalar auditors."""
    R,M,K,c,a,b,N,T=[I(p[k]) for k in ['R','M','K','c','amin','bmin','N0','T']]
    E,Lum,U,Ux,Ut,phi=[I(p[k]) for k in ['energy_geometric','luminosity_geometric','U','Ux','Ut','phi']]
    k=I(p['ray_kappa']);lam=K*U/(a*R**2);Z=lam/b
    g=4*phi*U/R+(1+1/b)*lam+(2+1/b)*K*Ut*pi/(2*c*a*a*R*sqrt(1-k))
    H=K/(a*R)*((Ux+Ut/c)/a+U/R*(1/b-1))
    tau=R/(c*a)*(Z*pi/(2*sqrt(1-k))+H/(N*N*(1-k)*sqrt(1-k)))
    W=K/(2*a*b)*(T/R**2+1/(c*a*R))
    result=dict(
        constraint_mass=c*E/M*(g+2*Z)*W,
        constraint_smooth_stress=2*c*E/M*lam*W,
        static_gradient_smooth=c*E/M*(4*lam/a+2*Z)*W,
        moving_constraint_shell=c*E/M*tau*K/(2*a*b*R**2),
        photon_shell_amplitude_and_front=K*pi/(4*M*a*a*R*sqrt(1-k))*(E*g+Lum*tau),
        photon_shell_angle=E*K*H/(M*a*a*R*(1-k)*sqrt(1-k)))
    extra=(M/R+K*K/(3*a*a*R**2))/(2*a*sqrt(b))
    eta=I(p['old_potential_contraction'])+extra
    free=sum(result.values());potential=eta/(1-eta)*free
    result.update(additional_potential=potential,
        homogeneous_mass_cm=E*(g+Z),lambda_bound=lam,zeta_bound=Z,
        photon_energy_fraction=g,photon_delay_seconds=tau,H_bound=H,
        weighted_potential_addition=extra,potential_contraction=eta,free_scalar_norm=free)
    return result


def apply():
    assert read(run.OUT/'audit.json')['passed'] and not read(run.OUT/'result.json')['passed']
    prior_result=read(run.OUT/'bound-result.json');final=read(matched.OUT/'final-result.json')
    assert prior_result['analytic_envelope_verified'] and final['passed'] and final['audit_passed']
    run.previous.prior.initialize();d=run.previous.incident.Driver(8);m=run.previous.exterior.Exterior()
    iv=mp.iv;iv.dps=50;I=iv.mpf
    def B(x):
        n,den=float(x).as_integer_ratio();return I(n)/I(den)
    def up(x):return float(np.nextafter(float(x.b),np.inf))
    def lo(x):return float(np.nextafter(float(x.a),-np.inf))
    R=B(m.r0);M=B(m.M);K=B(abs(m.K));c=B(run.C);G=B(run.G);N=B(m.N0);T=B(d.T)
    b=1-2*M/R;a=N*iv.sqrt(b);k=M/(R*b)+K*K/(2*a*a*R**2)
    Bn=prior_result['Born_source'];caps=prior_result['emission_caps'];U0=B(run.previous.incident.ETA)*R
    Ux0=U0/(c*B(d.D))*I(2048)/70
    # Save exact enclosing endpoints as decimal strings, not float-rounded
    # derived parameters. The auditor reconstructs each interval independently.
    values=dict(R=R,M=M,K=K,c=c,N0=N,T=T,bmin=b,amin=a,ray_kappa=k,
        energy_geometric=G/c**4*B(caps['through_T_erg']),luminosity_geometric=G/c**4*B(caps['binwise_luminosity_erg_per_second']),
        U=U0+B(Bn['U']),Ux=Ux0+B(Bn['Ux']),Ut=c*Ux0+B(Bn['Ut']),
        phi=B(abs(float(d.z0['phi'][0]))),old_potential_contraction=B(prior_result['global_potential_contraction']))
    # Serialize the exact dyadic endpoints and verify the interval round trip.
    mp.mp.dps=250
    p={key:[mp.nstr(mp.mpf(x._mpi_[j]),250) for j in [0,1]] for key,x in values.items()}
    assert all(I(p[key])._mpi_==x._mpi_ for key,x in values.items())
    write(OUT/'interval-inputs.json',p)
    out=terms({k:I(v) for k,v in p.items()},I,iv.sqrt,iv.pi)
    assert up(out['potential_contraction'])<1
    bg=np.load(run.previous.BACKGROUND)
    emax=B(np.nextafter(float(np.max(bg['epsilon'])),np.inf))
    qmax=B(np.nextafter(float(np.max(abs(bg['normalized']))),np.inf))
    kb=B(np.nextafter(final['background_mass_parameter_maximum'],np.inf))
    D0=1-emax;denominator=D0-kb;assert lo(denominator)>0
    scalar=(out['free_scalar_norm']+out['additional_potential'])/denominator
    mass_k=out['homogeneous_mass_cm']/M
    # Enclose the new port's rational denominator as well; this does not
    # promote the other frozen components to nonlinear physical accuracy.
    coefficient=(K/M+qmax)*D0/denominator
    mass=coefficient*mass_k/(denominator-mass_k)
    old=B(prior_result['total_selected_envelope'])/(1-kb/D0)**2
    new=scalar+mass;total=old+new;center=I(final['selected_decimal'][-1])
    fraction=up(total/abs(center));interval=center+I([-up(total),up(total)])
    result=dict(classification='Counterexample candidate',passed=fraction<.01,
        bounds={key:up(v) for key,v in out.items()},new_scalar_charge_envelope=up(scalar),new_homogeneous_mass_charge_envelope=up(mass),
        new_sector_charge_envelope=up(new),previous_selected_exterior_envelope_remapped=up(old),
        combined_selected_exterior_envelope=up(total),combined_envelope_over_selected=fraction,
        selected_endpoint=float(center.mid),conditional_selected_interval=[lo(interval),up(interval)],
        selected_sign_survives=up(interval)<0,old_arrived_photon_mass_counted_once=True,
        original_point_verdict='FAILED and unchanged',original_point_controls=prior_result['original_point_controls'],
        emission_family=read(OUT/'plan.json')['family'],exclusions=read(OUT/'plan.json')['exclusions'],
        direct_photon_mixed_sector_enclosed=True,additional_exterior_mixed_stress_closed=False,
        exact_ADM_conservation_verified=False,physical_final_charge_solved=False,full_goal_complete=False,
        scope='A contribution envelope around the stored Phase163 selected endpoint. It is not a confidence interval or enclosure of the body response, complete exterior source or physical/observational error.')
    write(OUT/'result.json',result);print(json.dumps(result),flush=True);assert result['passed'],result


def audit():
    result=read(OUT/'result.json');assert result['passed'];mp.mp.dps=100
    p=read(OUT/'interval-inputs.json');mid={k:(mp.mpf(v[0])+mp.mpf(v[1]))/2 for k,v in p.items()}
    actual=terms(mid,mp.mpf,mp.sqrt,mp.pi);errors={}
    for name,value in actual.items():
        upper=mp.mpf(float(result['bounds'][name]));assert upper>=value,(name,upper,value)
        errors[name]=float(abs(upper-value)/max(abs(value),mp.mpf('1e-290')))
    assert max(errors.values())<1e-10
    # Independent positive-control radial integrals, including the grazing
    # endpoint, and a nonzero moving shell paired with a test function.
    integrals=[mp.quad(lambda y:mp.sin(y),[0,mp.pi/2]),
        mp.quad(lambda y:2/(1+y*y),[0,1]),
        mp.quad(lambda y:y*(1-y),[0,1])]
    assert max(abs(a-b) for a,b in zip(integrals,[1,mp.pi/2,mp.mpf(1)/6]))<mp.mpf('1e-90')
    zero=dict(mid,energy_geometric=mp.mpf(0),luminosity_geometric=mp.mpf(0))
    z=terms(zero,mp.mpf,mp.sqrt,mp.pi)
    names=['constraint_mass','constraint_smooth_stress','static_gradient_smooth','moving_constraint_shell',
           'photon_shell_amplitude_and_front','photon_shell_angle','additional_potential','homogeneous_mass_cm']
    assert all(z[k]==0 for k in names) and all(actual[k]>0 for k in names)
    double=dict(mid,energy_geometric=2*mid['energy_geometric'],luminosity_geometric=2*mid['luminosity_geometric'])
    doubled=terms(double,mp.mpf,mp.sqrt,mp.pi);assert all(doubled[k]==2*actual[k] for k in names)
    verdict=dict(classification='Counterexample candidate',passed=True,symbolic=symbolic(),
        scalar_arithmetic_relative=errors,zero_source_and_nonzero_controls=True,linear_inventory_scaling=True,
        interval_precision=50,independent_scalar_precision=100,
        limitation='Arithmetic, distribution and elementary integral checks for the conditional inequalities. They do not independently certify the saved vacuum/EOS, emission caps against a continuum evolution, or omitted scalar sectors.',
        physical_final_charge_solved=False,full_goal_complete=False)
    write(OUT/'audit.json',verdict);print(json.dumps(verdict),flush=True)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));run.previous.incident.native.deadline(CAPS[action])
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
            assert sum(read(p)['seconds'] for p in OUT.glob('*-receipt.json'))+CAPS[action]<=TOTAL
        globals()[action]()
    except Exception as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(action=action,seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
