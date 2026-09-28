"""Bound the actual exterior terms without rescuing failed point controls.

Counterexample candidate. Outward interval arithmetic on specified saved
binary inputs and an exact-vacuum idealization, not a physical EOS certificate.
"""
from pathlib import Path
import json, resource, time
import numpy as np
import mpmath as mp
import sympy as sp
from numpy.polynomial import legendre as leg
import couple_native_exterior_incident as run


def main():
    out=run.OUT;assert not (out/'bound-plan.json').exists()
    old=run.read(out/'result.json');assert not old['passed']
    boundfile=run.previous.ROOT/'def-native-global-scalar-closure/bound.json'
    bp=boundfile.with_name('bound-plan.json')
    files=[Path(__file__),out/'plan.json',out/'result.json',out/'production-receipt.json',boundfile,bp,
           run.previous.OUT/'direct-g8.npz',run.previous.OUT/'charge-parts.npz',run.previous.BACKGROUND]
    run.write(out/'bound-plan.json',dict(classification='Counterexample candidate',
        decision='Can the actual omitted exterior radiation-metric interaction erase the stored matter-mediated charge? Bound the whole selected contribution, independently of the failed point quadrature/history criterion.',
        original_point_verdict='FAILED and unchanged:2.3616percent background-history comparison and0.26206percent shell temporal-quadrature comparison.',
        scope='Primary incoming wave on the direct emitted-radiation metric; the continuous finite-source first-Born field on that metric; static potential repetitions of these sources; photon arrived-mass correction from the primary plus that finite-source Born metric. Reciprocal mixed scalar constraints, scalar readout changes from perturbed packet paths and interior evolving operators remain open.',
        family='Any nonnegative angular emission with the recorded energy caps throughD andT and recorded bin-wise luminosity caps. This contains both saved histories; it is not an error bound on an unknown continuum emission law.',
        method='Exact total variation of the C3 pulse; positive emission bounds; m<=ADM, N>=N0, b>=1-2M/r0 and a>=N0*sqrt(bmin) in the exact static scalar vacuum. Bound the finite Born source polynomial by absolute inverse-Legendre column sums using exact binary entries and50-digit outward intervals. No new rays, fluid steps, clocks or quadrature orders.',
        gate='Selected total envelope below1percent of the Phase158 stored body endpoint; this is a separate absolute-contribution conclusion and does not turn the failed point result into a pass.',
        budget_seconds=90,from_original_reserved_bounds_budget=True,CPU_threads=1,virtual_GiB=3,
        stop='No automatic tighter source/time history or new physical trajectory if the envelope is inconclusive.',
        bindings={str(p):run.sha(p) for p in files}))
    start=time.monotonic();cpu=time.process_time();run.previous.incident.native.deadline(90)
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
    run.previous.prior.initialize();d=run.previous.incident.Driver(8);m=run.previous.exterior.Exterior()
    frozen=run.read(boundfile);old_bind=run.read(bp)['bindings'];bgfile=run.previous.wave.base.flow.INPUT/'balanced-20.npz'
    assert run.sha(bgfile) in [h for p,h in old_bind.items() if Path(p).name=='balanced-20.npz']
    iv=mp.iv;iv.dps=50;I=iv.mpf
    def B(x):
        n,den=float(x).as_integer_ratio();return I(n)/I(den)
    def up(x):return float(np.nextafter(float(x.b),np.inf))
    def low(x):return float(np.nextafter(float(x.a),-np.inf))
    R=B(m.r0);M=B(m.M);K=B(abs(m.K));cc=B(run.C);gg=B(run.G);D=B(d.D);T=B(d.T);eta=B(run.previous.incident.ETA)
    N0=B(m.N0);bmin=1-2*M/R;amin=N0*iv.sqrt(bmin)
    kappa=M/(R*bmin)+K*K/(2*amin*amin*R*R);assert up(kappa)<1
    Zmax=(1+cc*D/(bmin*R))/(amin*R)
    U0=eta*R;Ux0=U0/(cc*D)*I(2048)/70;Uxx0=U0/(cc*D)**2*I(2048)/5
    # Bernstein coefficients give |p'|<=2048/70 and |p''|<=2048/5.
    # The exact integral of |p'| is2 because p rises from0 to1 then falls to0.
    energies=[];lum_caps=[];port=[]
    for n in [64,128]:
        L,h,e=run.emission(n);aw=[I(2*j+1)/32 for j in range(4)]
        energy=[]
        for count in [n//2,n]:energy.append(sum((B(h)*aw[j]*sum((B(v) for v in L[:count,j]),I(0)) for j in range(4)),I(0)))
        energies.append(energy);lum_caps.append(sum((aw[j]*B(max(L[:,j])) for j in range(4)),I(0)));port.append(e)
    ED=I(max(up(v[0]) for v in energies));ET=I(max(up(v[1]) for v in energies));Lmax=I(max(up(v) for v in lum_caps))
    eD=gg/cc**4*ED;eT=gg/cc**4*ET
    boundary=eD/M*U0*Zmax
    volume=cc/2*eD/M*D*(cc*D)*U0*(2+4*Zmax*M)/(bmin*R**3)
    shell=cc/2*eD/M*D*U0/R**2
    primary=boundary+volume+shell
    # The actual first-Born source, including its initial exterior support.
    z=np.load(run.previous.OUT/'direct-g8.npz');q=8
    faces=np.r_[d.optical(d.edges),d.ex[1:]];mid=(faces[:-1]+faces[1:])/2;half=np.diff(faces)/2
    xx=(z['source_x'].reshape(-1,q)-mid[:,None])/half[:,None]
    inverse=np.linalg.inv(leg.legvander(xx,q-1));max_source=np.max(abs(z['source']),axis=0).reshape(-1,q)
    dx=z['source_dx'].reshape(-1,q);integral=I(0);peak=I(0);cell_bounds=[]
    for j in range(len(half)):
        envelope=sum((sum((B(abs(v)) for v in inverse[j,:,k]),I(0))*B(max_source[j,k])/B(dx[j,k]) for k in range(q)),I(0))
        integral+=(B(faces[j+1])-B(faces[j]))*envelope;peak=I(max(up(peak),up(envelope)));cell_bounds.append(up(envelope))
    dtmin=min(float(b-a) for a,b in zip(z['source_t'][:-1],z['source_t'][1:]))
    # PPoly linear time source, zero before0 and afterT. Its two endpoint
    # jumps and the local source are bounded in the second derivative.
    S=integral;St=2*S/B(dtmin);UB=cc*T/2*S;UBt=cc/2*S;UBx=S/2;UBxx=St/(2*cc)+3*peak
    Vgeom=2*M/(bmin*R**3)
    # On the outgoing strip0<=u<=T, a photon outside r obeys xp-x<=cT.
    # zeta_unit=2*I(r,rp)+(1-mu_p^2)/(rp*a_p), giving an integrable tail.
    Zint=2*cc*T/(amin**2*bmin*R)+1/(2*N0**2*amin**2)
    Xint=2/(amin*bmin*R);T_int=(2/bmin+2*(1+kappa))/(amin*R);H_int=1/(bmin*R**2)
    born_smooth=cc*T/2*eT/M*(2*Zint*(UBxx+Vgeom*UB)+Xint*UBx+T_int*UBt/cc+H_int*UB)
    born_shell=eT/(2*M*N0**2*amin**2*iv.sqrt(1-kappa))*(UBx+UBt/cc+UB/R)
    born=born_smooth+born_shell
    # Uniform free-field norm for the primary source, before reducing its
    # derivatives only at the outer observer. Rebuild, never reuse old norms.
    smooth=2*Zmax*(Uxx0+Vgeom*U0)+(4/bmin+2+2*kappa)/R**2*Ux0+2/(bmin*R**3)*U0
    front=(Ux0+U0/R)/(R*amin)
    free_norm=cc/2*eD/M*D*(cc*D*smooth+front)
    contraction=I(frozen['global_potential_contraction']);assert up(contraction)<1
    potential=contraction/(1-contraction)*(free_norm+born)
    # Global photon energy and arrival response for the primary+finite Born
    # input. The 1/mu launch singularity is bounded by integrable vacuum rays.
    U=U0+UB;Ux=Ux0+UBx;Ut=cc*Ux0+UBt;phi=B(abs(float(d.z0['phi'][0])))
    lam=K*U/(amin*R**2);zeta=lam/bmin;nu=(1+1/bmin)*lam
    energy_fraction=4*phi*U/R+nu+(2+1/bmin)*K*Ut*iv.pi/(2*cc*amin**2*R*iv.sqrt(1-kappa))
    Hmax=K/(amin*R)*((Ux+Ut/cc)/amin+U/R*(1/bmin-1))
    delay=R/(cc*amin)*(zeta*iv.pi/(2*iv.sqrt(1-kappa))+Hmax/(N0**2*(1-kappa)**I('1.5')))
    arrived=ET*energy_fraction+Lmax*delay
    background=np.load(run.previous.BACKGROUND);e0=B(max(background['epsilon']));q0=B(max(abs(background['normalized'])))
    mass=(K/M+q0)*gg/cc**4/M*arrived/(1-e0-gg/cc**4/M*arrived)
    total=(primary+born+potential)/(1-e0)+mass
    target=B(abs(run.read(run.previous.OUT/'final-result.json')['endpoint']['body']))
    fraction=up(total/target);bound=up(total);old_endpoint=float(run.read(run.previous.OUT/'final-result.json')['endpoint']['body'])
    assert abs(old['endpoint_exterior'])<up(primary)
    # Symbolic identities behind the tail and near-grazing bounds.
    v=sp.symbols('v',positive=True)
    assert sp.integrate(v**4*(1-v)**4,(v,0,1))==sp.Rational(1,630)
    p=256*v**4*(1-v)**4
    assert sp.integrate(sp.diff(p,v),(v,0,sp.Rational(1,2)))-sp.integrate(sp.diff(p,v),(v,sp.Rational(1,2),1))==2
    theta=sp.symbols('theta',real=True)
    assert sp.integrate(sp.sin(theta),(theta,0,sp.pi/2))==1
    run.write(out/'bound-symbolic.json',dict(classification='Proven',passed=True,
        pulse='The input pulse has total variation2; Bernstein coefficient bounds give |pprime|<=2048/70 and |psecond|<=2048/5.',
        vacuum='a_prime/a=2m/(r^2*b), hence I(r,rp)=1/(r*a)-1/(rp*a_p) and zeta_unit=2*I+(1-mu_p^2)/(rp*a_p). Positive static scalar energy gives m<=M,N>=N0.',
        rays='For kappa>=r*nu_prime, mu^2>=(1-kappa)*(1-(r0/r)^2). With h=delta_lnE-delta_nu-delta_lnL, h(r0)=0 and h_prime=-delta_nu_prime-mu*delta_lambda_t/(c*a). This gives |h|<=Hmax*(1-r0/r).',
        delay='delta_tprime=[-zeta-(1-mu^2)*h/mu^2]/(c*a*mu). The two dimensionless integrals are pi/2 and int_0^1(1-z)/(1-z^2)^(3/2)dz=1; the grazing limit stays finite.',
        mass='Instantaneous Hamiltonian/reference-energy conversion at launch is delta_lnE=u+delta_nu; propagation adds integral(delta_nu_t-mu^2*delta_lambda_t)dt. Positive cumulative arrival changes by at most E_total*energy_fraction+Lmax*delay.',
        scope='Algebraic inequalities under the named vacuum, continuous finite-source and nonnegative emission hypotheses. They do not certify microscopic or evolving-background closure.'))
    np.savez_compressed(out/'bound-inputs.npz',source_cell_envelope=cell_bounds,source_faces=faces,
        source_inverse=inverse,source_node_absolute_maximum=max_source,source_dx=dx)
    result=dict(classification='Counterexample candidate',analytic_envelope_verified=fraction<.01,
        original_point_readout_passed=False,original_point_controls=old['controls'],original_component_controls=old['component_controls'],
        primary_boundary_bound=up(boundary),primary_volume_bound=up(volume),primary_shell_bound=up(shell),primary_total_bound=up(primary),
        finite_Born_input_interaction_bound=up(born),additional_static_potential_bound=up(potential),
        photon_arrival_mass_bound=up(mass),total_selected_envelope=bound,envelope_over_stored_body=fraction,
        conditional_interval_around_stored_body=[float(np.nextafter(old_endpoint-bound,-np.inf)),float(np.nextafter(old_endpoint+bound,np.inf))],
        same_input_radiation_metric_cannot_erase_stored_body=fraction<1,
        emission_caps=dict(through_D_erg=up(ED),through_T_erg=up(ET),binwise_luminosity_erg_per_second=up(Lmax),port_relative=port),
        photon_energy_fraction_bound=up(energy_fraction),photon_arrival_shift_seconds_bound=up(delay),
        Born_source=dict(integral=up(S),absolute=up(peak),U=up(UB),Ut=up(UBt),Ux=up(UBx),Uxx=up(UBxx)),
        global_potential_contraction=up(contraction),new_free_field_norm=up(free_norm),
        scope='Only the listed exterior contributions, for the saved emission caps and continuous finite-source input on the exact static-vacuum idealization. This interval does not enclose the previous body response error or all physical exterior effects.',
        remaining=['Reciprocal mixed scalar-metric constraints and already-generated scalar background response.','Scalar readout change from perturbed packet paths, beyond the arrived-mass bound.','Interior evolving-background operator and input interpolation errors.','Full EOS, uniform derivatives, spatial/boundary/nonlinear errors and static/observational comparison.'],
        physical_exterior_closed=False,full_goal_complete=False,seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
        peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    run.write(out/'bound-result.json',result);print(json.dumps(result),flush=True)
    assert result['analytic_envelope_verified'],result


if __name__=='__main__':main()
