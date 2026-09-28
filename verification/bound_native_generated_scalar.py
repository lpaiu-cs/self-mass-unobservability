"""Counterexample candidate: generated scalar mixed exterior charge/flux.

Use the saved current compact source and actual photon inventory. Enclose the
declared finite-source continuation; do not identify it with a complete future
stellar history or use a flux identity as evidence of exact ADM conservation.
"""
from pathlib import Path
import json,resource,sys,time
import numpy as np
import mpmath as mp
import sympy as sp
import bound_native_exterior_mixed as previous

OUT=Path('native-generated-scalar165-work')
background=previous.matched.previous.mixed.previous
inf=previous.run.previous
read,write,sha=previous.read,previous.write,previous.sha
CAPS=dict(prepare=30,norms=60,apply=30,audit=45)
TOTAL=sum(CAPS.values())
iv=mp.iv;iv.dps=50;I=iv.mpf


def B(x):
    n,d=x.as_integer_ratio();return I(n)/I(d)


def up(x):return float(np.nextafter(float(x.b),np.inf))
def lo(x):return float(np.nextafter(float(x.a),-np.inf))


def save_intervals(path,values):
    mp.mp.dps=250
    p={k:[mp.nstr(mp.mpf(v._mpi_[j]),250) for j in [0,1]] for k,v in values.items()}
    assert all(I(p[k])._mpi_==v._mpi_ for k,v in values.items())
    write(path,p)


def symbolic():
    a,b,r,Phi,ax,V,c=sp.symbols('a b r Phi ax V c',nonzero=True)
    U,W,Ux,Wx,Ut,Wt,Uxx,Wxx,S,T=sp.symbols('U W Ux Wx Ut Wt Uxx Wxx S T')
    constraint=Ux*Wx+Ut*Wt/c**2-a/r*(Ux*W+U*Wx)+a*a*U*W/r**2-2*a*a*Phi*Phi*U*W/b
    energy=Ux*Wx+Ut*Wt/c**2+V*U*W
    boundary=(ax/r-a*a/r**2)*U*W+a/r*(Ux*W+U*Wx)
    assert sp.simplify((constraint-energy+boundary).subs(V,ax/r-2*a*a*Phi*Phi/b))==0
    x,t=sp.symbols('x t',real=True);u=sp.Function('u')(t,x);w=sp.Function('w')(t,x);v=sp.Function('v')(x)
    en=sp.diff(u,x)*sp.diff(w,x)+sp.diff(u,t)*sp.diff(w,t)/c**2+v*u*w
    flux=sp.diff(u,t)*sp.diff(w,x)+sp.diff(w,t)*sp.diff(u,x)
    eq=sp.diff(en,t)-sp.diff(flux,x)
    eq=eq.subs(sp.diff(u,t,2),c*c*(sp.diff(u,x,2)-v*u+S)).subs(sp.diff(w,t,2),c*c*(sp.diff(w,x,2)-v*w+T))
    assert sp.simplify(eq-S*sp.diff(w,t)-T*sp.diff(u,t))==0
    # Oppositely directed free waves cancel their leading mixed energy.
    assert (Ux*Wx+Ut*Wt/c**2).subs({Ux:-Ut/c,Wx:Wt/c})==0
    z=sp.Function('z')(t,x);q=sp.Function('q')(t,x);vg=sp.symbols('vg')
    op=2*z*(sp.diff(q,x,2)-vg*q)+sp.diff(z,x)*(sp.diff(q,x)-a*q/r)+sp.diff(z,t)*sp.diff(q,t)/c**2
    div=sp.diff(z*sp.diff(q,x),x)+sp.diff(z*sp.diff(q,t),t)/c**2
    rhs=div-2*z*vg*q-z*(sp.diff(q,t,2)/c**2-sp.diff(q,x,2))-a*sp.diff(z,x)*q/r
    assert sp.simplify(op-rhs)==0
    mu=sp.symbols('mu');assert sp.cancel((1-mu*mu)/(1-mu))==1+mu
    return dict(classification='Proven',passed=True,
        mass='For the two scalar legs with J_B=J_I=0, C_x=U_Bx*U_Ix+U_Bt*U_It/c^2+V*U_B*U_I-d_x(a*U_B*U_I/r). Use C_part=-integral_x^infinity E_cross dy-a*U_B*U_I/r, with zero extra asymptotic mass; its generally nonzero inner port must be matched.',
        work_flux='E_cross_t=F_cross_x+U_Bt*S_I+U_It*S_B, F_cross=U_Bt*U_Ix+U_It*U_Bx. Thus source work and the inner flux are required; an algebraic identity is not a measured conservation pass.',
        infinity='For outgoing waves, F_cross=-2*U_Bt*U_It/c. The signed emitted scalar mass length is 2/c integral U_Bt*U_It du; the incoming primary vanishes at fixed outgoing u and infinity.',
        operator='O_zeta[U]=d_x(zeta*U_x)+d_t(zeta*U_t/c^2)-2*zeta*Vgeom*U-zeta*(U_tt/c^2-U_xx)-a*zeta_x*U/r. This avoids numerical second differentiation.',
        packet_derivative='At a past-cone intersection the moving-shell factor (1-mu^2)/(1-mu)=1+mu<=2. A timelike outward ray crosses that past cone at most once.',
        scope='Conditional first-variation identities on the initial static scalar vacuum. No physical continuation, EOS or exact ADM verdict follows from the algebra alone.')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    files=[Path(__file__),Path(background.__file__),Path(inf.wave.base.__file__),
        previous.OUT/'result.json',previous.OUT/'audit.json',previous.OUT/'interval-inputs.json',
        previous.run.OUT/'bound-result.json',previous.matched.OUT/'final-result.json',inf.BACKGROUND,
        background.FIELDS/'source-128.npz',background.FIELDS/'fields-128-g8.npz',
        inf.OUT/'direct-g8.npz',inf.ROOT/'def-native-global-scalar-closure/bound.json',
        inf.wave.base.flow.INPUT/'balanced-20.npz']+background.PACKETS
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='a57c0d5da',
        previous_turn='Progress: Phase164 bounded the declared direct photon mixed source and applied it to the selected charge without accepting the failed point estimate.',
        claim='Connect the actually generated scalar leg to exterior mixed Einstein constraints, both scalar metric actions, the static-gradient cross term and outgoing scalar energy normalization.',
        decision='Can this omitted sector erase the selected response? Include its signed scalar energy channel in the mass-normalized readout, and preserve unresolved ADM/source-work conditions.',
        source='The saved Phase160 current33-knot compact forcing, its declared spatial polynomial, original positive emission caps, and all three actual signed background SDIRK photon histories. Input I is the same primary plus continuous finite-source first Born.',
        continuation='Only the outgoing observer strip t>=0,t-x/c<=T is bounded. Exterior scalar fields there use the saved compact source times and paired photons. Integrate the mixed constraint from larger radii toward the inner interface, so no unmeasured later inner history is substituted. The resulting particular solution has zero extra asymptotic mass and a required inner port, whose matching is still open.',
        sector='B is the generated scalar from compact forcing and compensated radiation, including all static canonical potential repetitions. Its scalar-only metric has lambda_B=Phi*U_B. The separate homogeneous background mass-port forcing is not silently added to this leg.',
        exclusions=['Homogeneous background mass-port scalar forcing and a conserved physical future continuation.','Direct metric actions of later signed photons, exterior input from actual material returns and all mutual material feedback of the newly bounded waves.','Initial ADM/momentum/flux matching, full direct-field/EOS/derivative/spatial/boundary/nonlinear errors and static/observational comparison.'],
        method='Outward50-digit input/source envelopes. Use the exact mixed energy-plus-boundary identity and integration by parts of both metric actions. Bound the zero-asymptotic-mass particular field and outgoing scalar flux separately; compare the finite-T required port with the existing homogeneous mass. No finite difference of tiny waveforms, fluid replay or new rays.',
        gates=dict(combined_selected_fraction=.01,potential_contraction=1.,arithmetic_relative=1e-10),
        budget=dict(actions=CAPS,total_seconds=TOTAL,CPU_threads=1,virtual_GiB=3,new_fluid_steps=0,new_EOS_roots=0,new_rays=0),
        measured_basis='Previous input load and interval calculation4.60s. One current-source setup plus the same polynomial bound pattern is capped at60s; all actions165s. This is not a forecast for any longer physical history.',
        stop='No automatic additional source knots, spatial order, evolution interval or relaxed gate. An insufficient envelope remains insufficient.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)}))
    write(OUT/'symbolic.json',symbolic())


def norms():
    assert read(previous.OUT/'audit.json')['passed']
    background.initialize();d=dict(np.load(background.FIELDS/'source-128.npz'));m=background.wave.Response();m.setup(d,8)
    f=np.load(background.FIELDS/'fields-128-g8.npz');assert not np.any(f['U'][0])
    # Absolute inverse-Legendre column sums enclose every coefficient, without
    # treating the rounded product of a matrix multiply as an exact proof.
    nodes=np.max(abs(m.source),axis=0).reshape(-1,8);dx=m.dx.reshape(-1,8)
    source=I(0);cells=[]
    for j in range(len(m.half)):
        v=sum((sum((B(abs(x)) for x in m.inverse[j,:,k]),I(0))*B(nodes[j,k])/B(dx[j,k]) for k in range(8)),I(0))
        width=B(m.xfaces[j+1])-B(m.xfaces[j]);source+=width*v;cells.append(up(v))
    np.savez_compressed(OUT/'source-enclosure-inputs.npz',node_absolute_maximum=nodes,dx=dx,inverse=m.inverse,
        xfaces=m.xfaces-m.xfaces[-1],cell_absolute_envelope=cells)
    p={k:I(v) for k,v in read(previous.OUT/'interval-inputs.json').items()}
    old=read(previous.run.OUT/'bound-result.json');T=p['T'];c=p['c'];a=p['amin'];b=p['bmin'];R=p['R'];K=p['K'];M=p['M']
    signed=[];gamma=1-1/np.sqrt(2);aw=[I(2*j+1)/32 for j in range(4)]
    for path in background.PACKETS:
        z=np.load(path);L=z['accepted_angular_luminosity'];n=int(path.stem.split('-')[1]);assert L.shape==(2*n,4)
        h=T/n;value=sum((h*B(float([1-gamma,gamma][i%2]))*aw[j]*B(abs(L[i,j])) for i in range(2*n) for j in range(4)),I(0))
        signed.append(dict(path=str(path),total_variation_energy_erg=up(value)))
    E=p['energy_geometric']+B(previous.run.G)/c**4*sum((B(v['total_variation_energy_erg']) for v in signed),I(0))
    eta=p['old_potential_contraction'];Vnorm=2*eta/(c*T)
    radiation=c*T*E*K/(2*a*b*R**2)+E*K*iv.pi/(4*a*a*R*iv.sqrt(1-p['ray_kappa']))
    free=c*T/2*source+radiation;U=free/(1-eta)
    derivative=source/2+E*K/(2*a*b*R**2)+E*K/(a*R**2)+Vnorm*U/2
    # Use source CELL faces, including end-cell margins, not only Gauss nodes.
    driver=inf.incident.Driver(8);faces=np.r_[driver.optical(driver.edges),driver.ex[1:]]
    span=B(faces[-1])-B(faces[0])
    p.update(B_source_L1=source,B_U=U,B_derivative=derivative,B_radiation_free=radiation,
        B_total_photon_energy_geometric=E,I_Born_U=B(old['Born_source']['U']),I_Born_Ux=B(old['Born_source']['Ux']),
        I_Born_Ut=B(old['Born_source']['Ut']),I_Born_source_L1=B(old['Born_source']['integral']),
        I_primary_U=B(inf.incident.ETA)*R,I_primary_Ux=B(inf.incident.ETA)*R/(c*B(driver.D))*I(2048)/70,
        I_primary_duration=B(driver.D),I_Born_source_span=span,
        vacuum_V_L1=(M/R**2+2*K*K/(3*a*a*R**3))/iv.sqrt(b))
    save_intervals(OUT/'inputs.json',p)
    result=dict(classification='Counterexample candidate',current_compact_source_enclosed=True,
        scalar_U_bound=up(U),scalar_Ux_and_Ut_over_c_bound=up(derivative),compact_source_L1_bound=up(source),
        total_photon_variation_geometric_cm=up(E),signed_emission=signed,source_cells=len(cells),source_knots=len(d['t']),
        declared_source_continuation=read(OUT/'plan.json')['continuation'],full_goal_complete=False)
    write(OUT/'norms.json',result);print(json.dumps(result),flush=True)


def terms(p,Q,sqrt,pi):
    p={k:Q(v) for k,v in p.items()};R,M,K,c,a,b,T=[p[k] for k in ['R','M','K','c','amin','bmin','T']]
    U,D=p['B_U'],p['B_derivative'];F=p['I_primary_U']+p['I_Born_U'];P=p['I_primary_Ux']
    X,Y=p['I_Born_Ux'],p['I_Born_Ut']/c;AB=K*U/(a*b);AI=K*F/(a*b)
    length=2*(p['I_Born_source_span']+c*T)
    C=2*D*P*c*p['I_primary_duration']+D*(X+Y)*length+p['vacuum_V_L1']*U*F+U*F/R
    J=lambda n:T/((n-1)*R**(n-1))+1/(c*a*(n-1)*(n-2)*R**(n-2))
    w=K/(a*b)*J(3)
    value=dict(
        mixed_mass_source=c*C/M*w,
        mixed_constraint_and_static_gradient=c/(2*M)*4*K**3*U*F/(a*a*b)*(1+2/b)*J(6),
        metric_inner_boundary=c*T/(2*M*R**2)*(AB*(P+X)+AI*D),
        metric_cone_edge=(AB*(2*P+X+Y)+2*AI*D)/(2*M*a*R),
        metric_smooth=c/(2*M)*(4*M*AB*F/a*J(5)+AI*U/a*(2*M*J(5)+2*K*K/(a*a)*J(6))+4*K*U*F/b*J(4)),
        finite_Born_source=c*T/(2*M)*AB/R**2*p['I_Born_source_L1'],
        emitted_radiation_source=AI/R**2*p['B_radiation_free']/M)
    free=sum(value.values());eta=p['old_potential_contraction']+(M/R+K*K/(3*a*a*R**2))/(2*a*sqrt(b))
    value.update(additional_potential=free*eta/(1-eta),free_scalar_norm=free,potential_contraction=eta,
        intermediate_mass_constraint_cm=C,
        required_inner_port_at_T_cm=D*(X+Y)*length+p['vacuum_V_L1']*U*p['I_Born_U']+U*p['I_Born_U']/R,
        emitted_scalar_mass_cm=2*T*D*p['I_Born_Ut'])
    return value


def apply():
    p={k:I(v) for k,v in read(OUT/'inputs.json').items()};q=terms(p,I,iv.sqrt,iv.pi)
    assert up(q['potential_contraction'])<1
    prior=read(previous.OUT/'result.json');mass=read(previous.matched.OUT/'final-result.json');bg=np.load(inf.BACKGROUND)
    e=B(np.nextafter(float(max(bg['epsilon'])),np.inf));q0=B(np.nextafter(float(max(abs(bg['normalized']))),np.inf))
    kb=B(np.nextafter(mass['background_mass_parameter_maximum'],np.inf));D0=1-e;den=D0-kb
    # The particular field has zero extra asymptotic mass. Its inner port is
    # an unclosed matching condition, not a second scalar radiation inventory.
    delta_mass=q['emitted_scalar_mass_cm']/p['M']
    scalar=(q['free_scalar_norm']+q['additional_potential'])/den
    mass_bound=(p['K']/p['M']+q0)*D0/den*delta_mass/(den-delta_mass)
    added=scalar+mass_bound;total=B(prior['combined_selected_exterior_envelope'])+added
    center=I(mass['selected_decimal'][-1]);interval=center+I([-up(total),up(total)])
    fraction=up(total/abs(center))
    port=mass['rows'][0]['paired_endpoint_cm']
    port_center=I([np.nextafter(port,-np.inf),np.nextafter(port,np.inf)])
    port_bound=B(prior['bounds']['homogeneous_mass_cm'])+q['required_inner_port_at_T_cm']
    port_interval=port_center+I([-up(port_bound),up(port_bound)])
    result=dict(classification='Counterexample candidate',passed=fraction<.01,bounds={k:up(v) for k,v in q.items()},
        new_scalar_charge_envelope=up(scalar),new_mass_and_scalar_flux_charge_envelope=up(mass_bound),
        new_sector_charge_envelope=up(added),combined_selected_exterior_envelope=up(total),combined_envelope_over_selected=fraction,
        conditional_selected_interval=[lo(interval),up(interval)],selected_sign_survives=up(interval)<0,
        finite_T_selected_mass_port_interval_cm=[lo(port_interval),up(port_interval)],
        finite_T_selected_mass_port_excludes_zero=up(port_interval)<0 or lo(port_interval)>0,
        mass_port_scope='Existing stored paired port plus the Phase164 photon sector and this generated-scalar sector at the actual finite coordinate endpoint T. Other signed/input/homogeneous/initial-matching sectors remain open; this is not a full GR no-go theorem.',
        outgoing_scalar_energy_channel_included=True,generated_scalar_sector_enclosed=True,
        prescribed_continuation=read(OUT/'plan.json')['continuation'],exclusions=read(OUT/'plan.json')['exclusions'],
        original_point_verdict='FAILED and unchanged',additional_exterior_mixed_stress_closed=False,
        exact_ADM_conservation_verified=False,physical_final_charge_solved=False,full_goal_complete=False,
        scope='Contribution envelope for the zero-extra-asymptotic-mass exterior scalar particular solution and outgoing scalar flux. Its required inner port is exposed, not assumed matched. The interval excludes that matching correction, body-response and full physical errors.')
    write(OUT/'result.json',result);print(json.dumps(result),flush=True);assert result['passed'],result


def audit():
    result=read(OUT/'result.json');mp.mp.dps=100
    p={k:(mp.mpf(v[0])+mp.mpf(v[1]))/2 for k,v in read(OUT/'inputs.json').items()}
    exact=terms(p,mp.mpf,mp.sqrt,mp.pi);errors={}
    for k,v in exact.items():
        upper=mp.mpf(result['bounds'][k]);assert upper>=v,k
        errors[k]=float((upper-v)/max(abs(v),mp.mpf('1e-290')))
    assert max(errors.values())<1e-10
    # Positive controls distinguish radiation energy from opposite-direction
    # interference; a zero scalar leg must remove every new source and flux.
    zero=dict(p,B_U=mp.mpf(0),B_derivative=mp.mpf(0),B_radiation_free=mp.mpf(0))
    z=terms(zero,mp.mpf,mp.sqrt,mp.pi)
    assert all(v==0 for k,v in z.items() if k!='potential_contraction')
    assert exact['emitted_scalar_mass_cm']>0 and exact['intermediate_mass_constraint_cm']>0
    x=sp.symbols('x');pulse=x*x*(1-x)**2
    work=sp.integrate(2*sp.diff(pulse,x)**2,(x,0,1));assert work==sp.Rational(4,105)
    receipt=dict(classification='Counterexample candidate',passed=True,symbolic=symbolic(),scalar_arithmetic_relative=errors,
        zero_leg_control=True,nonzero_outgoing_flux_control=str(work),
        actual_mass_flux_balance_verified=False,
        limit='The same explicit inequalities are recalculated at100 digits; this is not an independent physical solution. Full ADM conservation requires the measured inner port and source work, not just the proven identity.',full_goal_complete=False)
    write(OUT/'audit.json',receipt);print(json.dumps(receipt),flush=True)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));inf.incident.native.deadline(CAPS[action])
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
