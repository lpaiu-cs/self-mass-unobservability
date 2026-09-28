"""Counterexample candidate: actual incident/SDIRK histories at null infinity.

Reuse the frozen finite-source and fixed-ray operators. Keep the incoming,
matter-mediated and normalized self-GR increments separate at every step.
This is not a complete evolving-exterior Einstein solution or an EFT test.
"""
from pathlib import Path
import inspect, json, resource, sys, textwrap, time
import numpy as np
import sympy as sp
import mpmath as mp
from numpy.polynomial import legendre as leg
import def_native_incident_drive as incident
import def_retained_native_return as prior
import def_native_global_scalar_closure as exterior
import def_native_characteristic_gr as wave

OUT=Path('native-incident-infinity158-work')
BEFORE=Path('native-incident-reciprocal156-work/sweep-2')
SELF=Path('native-incident-self-gr157-work')
ROOT=Path('outputs/direct-eos-gr33')
BACKGROUND=prior.CURRENT/'infinity/completed/retained-128-a8-r8.npz'
C=wave.C;G=exterior.G;LD=np.longdouble
read,write,sha=prior.read,prior.write,prior.sha
CAPS=dict(prepare=30,direct=120,packets=240,audit=60)


def inputs():
    files=[Path(__file__),Path(incident.__file__),Path(prior.__file__),Path(wave.__file__),
           Path(wave.base.__file__),Path(exterior.__file__),BACKGROUND,
           wave.base.flow.INPUT/'balanced-20.npz',prior.EV/'coupled-128.npz',
           prior.EV/'source-128.npz',SELF/'normalization.json',SELF/'block-result.json']
    for q in [4,8]:files.append(incident.FIELDS/f'born-g{q}.npz')
    for name,folder,photons in [('body',BEFORE,'photons-precise'),('self',SELF/'sweep-2','photons-common-gr')]:
        for n in [64,128]:
            files.extend([folder/photons/f'steps-{n}-reference-128.npz',folder/'gr'/f'wave-{n}-g8.npz'])
        files.append(folder/'gr'/'wave-128-g4.npz')
    return files


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='3d1976bd7',
        claim='Read the actual Phase156 incident-driven response and Phase157 self-GR correction at outgoing null infinity, with their signed SDIRK angular emission and the same retained mass normalization.',
        decision='Does the matter-mediated compact response survive the exterior photon and arrived-mass terms in the fixed initial GR readout? Separate it from the much larger directly scattered incident field and identify the unresolved error floor.',
        input='Existing eta=1e-30, C3 compact incoming pulse, 3.4344311179287023ms,531cells,17 output times,64/128 response clocks;33-time first-Born input representation at4/8 radial orders.',
        method='Reuse the existing characteristic-cut scalar propagator with independently specified observer times. A finite observer outside ALL compact source support at t=u+x_observer/c equals the outgoing characteristic limit of this declared free Green operator. The free incoming primary is zero there. Reuse the existing causal-support photon packet kernel without changes; integrate actual accepted signed SDIRK stages, never endpoint trapezoids.',
        normalization='Same saved fine retained background for both response clocks and all readout quadratures. Store separate direct, body and physical self-GR contributions; do not add a sub-ulp increment and then subtract the background. Rational normalization is an exact readout identity, not nonlinear evolution.',
        gates=dict(time=.02,quadrature=.002,port=1e-12,independent_observer=1e-10,normalization=1e-12),
        budgets=CAPS,total_seconds=sum(CAPS.values()),CPU_threads=1,virtual_GiB=3,
        new_EOS_calls=0,new_fluid_steps=0,new_background_steps=0,
        forecast='Measure three causal packet ages at the unchanged angular/radial orders. Twice the slowest cost times319 distinct ages plus5s must fit the remaining240s. Reuse each kernel for both actual histories and all output times.',
        stop='Stop on failed gates or measured budget. No longer horizon, new clocks, higher quadrature orders, additional Born feedback or fluid replay. Preserve original upstream-GR and whole-state iteration failures.',
        limits=['Initial frozen scalar coefficients and null ray geometry; response of background emitted photons to incident metric not yet included.',
                'Higher potential returns, evolving-background scalar scattering and full physical exterior coupling not enclosed.',
                'Time/4vs8 comparisons are diagnostics, not continuum/source/EOS/derivative error bounds.',
                'Same-input direct-only subtraction is bookkeeping, not a same-inventory static EFT or observational comparator.'],
        bindings={str(p):sha(p) for p in inputs()}))
    (OUT/'registered-producer.py').write_bytes(Path(__file__).read_bytes())


def observer_propagator():
    # ponytail: reuse the exact finite-source owner; only separate output times
    # from source knots. No new scalar PDE solver or input-field refinement.
    source=inspect.getsource(wave.Response.propagate)
    source=textwrap.dedent(source)
    assert source.count('for t in self.t:')==1
    source=source.replace('for t in self.t:','for t in self.output_times:')
    namespace=dict(wave.Response.propagate.__globals__)
    exec(compile(source,__file__,'exec'),namespace)
    return namespace['propagate'],source


def observe(m,source,times,observer):
    m.output_times=np.asarray(times)+observer/C;m.tx=np.array([observer])
    m.distance=abs(m.tx[:,None]-m.x[None,:])/C;m.sign=np.sign(m.tx[:,None]-m.x[None,:])
    assert observer>m.xfaces[-1]
    return observer_propagator()[0](m,source)[0][:,0]


def symbolic():
    a,s0,ds,de,e0=sp.symbols('alpha s0 ds de e0')
    q0=(s0+a*e0)/(1-e0)
    assert sp.factor((s0+ds+a*(e0+de))/(1-e0-de)-q0-(ds+(a+q0)*de)/(1-e0-de))==0
    x,y,u,c=sp.symbols('x y u c',real=True)
    # Under x>y and c>0, t-|x-y|/c at t=u+x/c is u+y/c.
    assert sp.expand(u+x/c-(x-y)/c-(u+y/c))==0
    assert sp.integrate(sp.Symbol('mu'),(sp.Symbol('mu'),sp.Symbol('lo'),sp.Symbol('hi')))==(sp.Symbol('hi')**2-sp.Symbol('lo')**2)/2
    return dict(classification='Proven',passed=True,
        scope='Exact rational charge-increment identity, exterior-observer retardation and angular arrival weight. No numerical or physical source error is certified.')


def direct():
    prior.initialize();start=time.monotonic();fn,expanded=observer_propagator()
    (OUT/'expanded-observer.py').write_text(expanded+'\n')
    box=[]
    for order in [4,8]:
        m=wave.Response.__new__(wave.Response);m.order=order;m.t=np.array([0.,.3,.7,1.])/C
        m.xfaces=np.array([-1.,0.,1.]);m.mid=np.array([-.5,.5]);m.half=np.array([.5,.5])
        gx,gw=leg.leggauss(order);m.x=(m.mid[:,None]+m.half[:,None]*gx).ravel()
        m.ids=np.repeat(np.arange(2),order);m.dx=(m.half[:,None]*np.broadcast_to(gw,(2,order))).ravel()
        m.inverse=np.linalg.inv(leg.legvander(np.broadcast_to(gx,(2,order)),order-1));m.gx,m.gw=leg.leggauss((order+3)//2)
        uu=np.array([-2.,-.5,0.,.4,1.,2.]);got=observe(m,(m.t*C)[:,None]*m.dx,uu/C,4.)
        primitive=lambda s:np.where(s<=0,0.,np.where(s<1,s**3/3,s-2/3))
        exact=(primitive(uu+1)-primitive(uu-1))/4
        box.append(float(np.max(abs(got-exact))))
    assert max(box)<3e-15,box
    rows=[];paths={}
    for order in [4,8]:
        d=incident.Driver(order);saved=np.load(incident.FIELDS/f'born-g{order}.npz')
        m=wave.Response.__new__(wave.Response);m.order=order;m.t=saved['t']
        m.xfaces=np.r_[d.optical(d.edges),d.ex[1:]];m.mid=(m.xfaces[:-1]+m.xfaces[1:])/2;m.half=np.diff(m.xfaces)/2
        m.x=saved['source_x'];m.dx=saved['dx'];m.ids=np.repeat(np.arange(len(m.mid)),order)
        xx=(m.x.reshape(-1,order)-m.mid[:,None])/m.half[:,None]
        m.inverse=np.linalg.inv(leg.legvander(xx,order-1));m.gx,m.gw=leg.leggauss((order+3)//2)
        source=-m.dx*saved['V']*incident.ETA*d.r0*incident.pulse((m.t[:,None]+m.x/C)/d.D)
        clock=np.linspace(0,d.T,17);U=observe(m,source,clock,3*C*d.T);check=observe(m,source,clock,4*C*d.T)
        error=float(np.max(abs(U-check))/max(np.max(abs(U)),1e-290));assert error<1e-10,error
        assert np.max(abs(incident.pulse((clock+6*d.T)/d.D)))==0
        paths[order]=-U/float(d.bg.z['ADM_mass'])
        np.savez_compressed(OUT/f'direct-g{order}.npz',t=clock,outgoing_U=U,normalized=paths[order],
            other_observer_U=check,source=source,source_x=m.x,source_dx=m.dx,source_t=m.t)
        rows.append(dict(order=order,independent_observer_relative=error,
            endpoint=float(paths[order][-1]),potential_norm_estimate=float(saved['eta'])))
    norm=max(float(np.max(abs(paths[8]))),1e-290)
    error=float(np.max(abs(paths[4]-paths[8]))/norm)
    eta=float(np.load(incident.FIELDS/'born-g8.npz')['eta'])
    remainder=eta**2/(1-eta)*incident.ETA*d.r0/float(d.bg.z['ADM_mass'])
    write(OUT/'symbolic.json',symbolic())
    result=dict(classification='Counterexample candidate',passed=error<.002,rows=rows,
        quadrature_relative=error,box_source_absolute_errors=box,
        incoming_primary_at_outgoing_null_infinity_zero=True,initial_exterior_data_preserved=True,
        direct_is_existing_first_potential_return=True,first_return_refed_to_fluid=False,
        higher_return_finite_estimate_normalized=remainder,
        estimate_scope='Historical finite coefficient operator estimate; not a null-infinity physical error enclosure.',
        seconds=time.monotonic()-start,full_goal_complete=False)
    write(OUT/'direct-result.json',result);assert result['passed'];print(json.dumps(result),flush=True)


def packet_factory():
    # Reuse the already verified causal-support kernel verbatim, including
    # its log(1-mu) angular coordinate and per-packet arrival front.
    source=inspect.getsource(prior.infinity)
    marker='    # Measure the new causal-support method before all saved packets.'
    assert source.count(marker)==1
    source=source.split(marker)[0]+'    return packet,m,kernels,lambda:inverse_error\n'
    namespace=dict(vars(prior),OUT=OUT)
    exec(compile(source,__file__,'exec'),namespace)
    (OUT/'expanded-packet.py').write_text(source)
    return namespace['infinity']()


def emitted(folder,n,m):
    p=np.load(folder/f'steps-{n}-reference-128.npz');h=m.T/n;gamma=1-1/np.sqrt(2)
    t=p['accepted_angular_times'];L=p['accepted_angular_luminosity']
    assert L.shape==(2*n,4)
    assert np.max(abs(t-h*(np.arange(n)[:,None]+[gamma,1.]).ravel()))<1e-18
    packets=LD(h)*np.tile([1-gamma,gamma],n)[:,None]*L.astype(LD)
    aw=np.arange(1,8,2,dtype=LD)/32
    cumulative=np.r_[LD(0),np.cumsum(packets@aw)];ids=2*np.rint(p['t']/h).astype(int)
    error=float(np.max(abs(cumulative[ids]-p['radial_ports'][:,1,1]))/max(np.sum(abs(packets)@aw),LD('1e-290')))
    assert error<1e-12,error
    return packets,error


def packets():
    prior.initialize();start=time.monotonic();packet,m,kernels,inverse_error=packet_factory()
    costs=[];gamma=1-1/np.sqrt(2)
    for age in [(1-gamma)*m.T/128,m.T/2,m.T*(128-gamma)/128]:
        mark=time.monotonic()
        for a,r in [(8,8),(4,8),(8,4)]:packet(age,a,r)
        costs.append(time.monotonic()-mark)
    upper=2*319*max(costs)+5;eligible=upper<CAPS['packets']-(time.monotonic()-start)
    write(OUT/'packet-pilot.json',dict(classification='Counterexample candidate',point_seconds=costs,
        upper_remaining_seconds=upper,eligible=eligible,distinct_age_forecast=319));assert eligible
    bg=np.load(BACKGROUND);assert np.max(abs(bg['t']-m.t))<1e-18
    alpha=LD(-m.K/m.M);e0=bg['epsilon'].astype(LD);q0=bg['normalized'].astype(LD)
    factor=LD(read(SELF/'normalization.json')['factor']);rows=[];paths={}
    for name,folder,photon,scale in [('body',BEFORE,'photons-precise',LD(1)),('self',SELF/'sweep-2','photons-common-gr',factor)]:
        emission={n:emitted(folder/photon,n,m) for n in [64,128]}
        for n,a,r in [(128,8,8),(128,4,8),(128,8,4),(64,8,8)]:
            mark=time.monotonic();pieces=[];p,error=emission[n];h=m.T/n
            for count in np.arange(17)*(n//16):
                value=np.zeros((2,4),LD)
                for j in range(count):
                    for stage,offset in enumerate([gamma,1.]):value+=packet((count-j-offset)*h,a,r)*p[2*j+stage]
                pieces.append(np.sum(value,axis=1,dtype=LD))
            pieces=np.asarray(pieces,LD)/scale;scalar=pieces[:,0];arrived=pieces[:,1]
            compact=np.load(folder/'gr'/f'wave-{n}-g8.npz')['free_scalar'].astype(LD)/scale
            deps=LD(G)/LD(C)**4/LD(m.M)*arrived
            # First derivative of the declared rational charge around the
            # SAME retained background. Cross-response terms stored below.
            mass=(alpha+q0)*deps
            delta=(compact+scalar+mass)/(1-e0)
            d=dict(t=m.t,compact=compact,exterior=scalar,arrived_energy_erg=arrived,
                epsilon_increment=deps,scalar_increment=compact+scalar,mass_numerator=mass,
                linear_charge_increment=delta,background_epsilon=e0,background_normalized=q0,
                standalone_rational_increment=(compact+scalar+mass)/(1-e0-deps))
            paths[name,n,a,r]=d;np.savez_compressed(OUT/f'{name}-{n}-a{a}-r{r}.npz',**d)
            rows.append(dict(component=name,steps=n,angular=a,radial=r,port_relative=error,seconds=time.monotonic()-mark))
    direct=np.load(OUT/'direct-g8.npz')['normalized'].astype(LD)
    body=paths['body',128,8,8];self_gr=paths['self',128,8,8]
    de=body['epsilon_increment']+self_gr['epsilon_increment']
    # The following identities preserve the tiny self correction as a
    # separate column. Summing it into body would round it away.
    denom=1-e0;full_denom=denom-de
    parts=np.column_stack([direct/denom,body['linear_charge_increment'],self_gr['linear_charge_increment']])
    rational_parts=parts*denom[:,None]/full_denom[:,None]
    np.savez_compressed(OUT/'charge-parts.npz',t=m.t,direct=parts[:,0],body=parts[:,1],self_gr=parts[:,2],
        rational_parts=rational_parts,normalization_cross_terms=parts*(de/full_denom)[:,None],
        epsilon_increment=de,background_epsilon=e0,background_normalized=q0,alpha0=alpha)
    controls={}
    for name in ['body','self']:
        controls[name]={}
        for key in ['compact','exterior','arrived_energy_erg','linear_charge_increment']:
            fine=paths[name,128,8,8][key];norm=max(float(np.max(abs(fine))),1e-290)
            controls[name][key]={tag:float(np.max(abs(paths[(name,)+p][key]-fine))/norm)
                for tag,p in [('angular',(128,4,8)),('radial',(128,8,4)),('time',(64,8,8))]}
    passed=all(v['angular']<.002 and v['radial']<.002 and v['time']<.02 for group in controls.values() for v in group.values())
    np.savez_compressed(OUT/'packet-kernels.npz',keys=np.asarray(list(kernels)),values=np.asarray(list(kernels.values())))
    result=dict(classification='Counterexample candidate',passed=passed,controls=controls,rows=rows,
        endpoint=dict(direct=float(parts[-1,0]),body=float(parts[-1,1]),self_GR=float(parts[-1,2]),
            body_compact=float(body['compact'][-1]),body_exterior=float(body['exterior'][-1]),body_mass_numerator=float(body['mass_numerator'][-1]),
            body_arrived_energy_erg=float(body['arrived_energy_erg'][-1]),self_arrived_energy_erg=float(self_gr['arrived_energy_erg'][-1]),
            body_over_compact=float(parts[-1,1]/body['compact'][-1]),self_over_body=float(parts[-1,2]/parts[-1,1])),
        computed_packet_kernels=len(kernels),delay_inverse_seconds=inverse_error(),seconds=time.monotonic()-start,
        same_background_mass_normalization=True,actual_signed_angular_packets_used=True,
        actual_response_at_fixed_operator_null_infinity=True,background_photon_metric_response_included=False,
        scalar_evolving_exterior_interaction_included=False,upstream_GR_uncertainty_closed=False,
        original_waveform_change_test_passed=False,static_EFT_compared=False,physical_final_charge_solved=False,full_goal_complete=False)
    write(OUT/'result.json',result);print(json.dumps(result),flush=True);assert passed


def audit():
    result=read(OUT/'result.json');assert result['passed'] and read(OUT/'direct-result.json')['passed']
    d=np.load(OUT/'charge-parts.npz');body=np.load(OUT/'body-128-a8-r8.npz');self_gr=np.load(OUT/'self-128-a8-r8.npz')
    mp.mp.dps=100;errors=[];mp_exact=[];cross=[]
    for j in range(len(d['t'])):
        B=lambda v:mp.mpf(str(v))
        a=B(d['alpha0']);e0=B(d['background_epsilon'][j]);q0=B(d['background_normalized'][j])
        direct=B(np.load(OUT/'direct-g8.npz')['normalized'][j]);de=B(body['epsilon_increment'][j])+B(self_gr['epsilon_increment'][j])
        ds=direct+B(body['scalar_increment'][j])+B(self_gr['scalar_increment'][j]);s0=q0*(1-e0)-a*e0
        exact=(s0+ds+a*(e0+de))/(1-e0-de)-q0
        stable=(ds+(a+q0)*de)/(1-e0-de)
        assert abs(exact-stable)<mp.mpf('1e-110'),(j,exact,stable)
        stored=sum((B(x) for x in d['rational_parts'][j]),mp.mpf(0))
        errors.append(float(abs(stored-stable)/max(abs(stable),mp.mpf('1e-290'))))
        mp_exact.append(str(exact));cross.append(float(de/(1-e0-de)))
    assert max(errors)<1e-12,errors
    # Independent angular arrival identity: exact bin integral, including a
    # front at each physical edge and strictly inside the radial bin.
    for mu0 in [0.,.25,.5,.75,.9,1.]:
        edges=np.linspace(0,1,5);arr=np.maximum(edges[1:]**2-np.maximum(edges[:-1],mu0)**2,0)/2
        assert abs(arr.sum()-(1-mu0**2)/2)<2e-16
    direct=read(OUT/'direct-result.json');endpoint=result['endpoint']
    result.update(audit=dict(classification='Counterexample candidate',passed=True,
        hundred_digit_normalization_max_relative=max(errors),hundred_digit_exact_difference=mp_exact,
        maximum_rational_cross_term_fraction=max(abs(x) for x in cross),
        conditional_input_remainder_scale_over_body=direct['higher_return_finite_estimate_normalized']/abs(endpoint['body']),
        direct_readout_quadrature_difference_over_body=direct['quadrature_relative']*max(abs(np.load(OUT/'direct-g8.npz')['normalized']))/abs(endpoint['body']),
        scope='Independent high-precision rational readout and exact arrival integral. Neither the finite Born estimate nor4vs8 difference is a complete error bound.'))
    result['contribution']='Conditional loophole progress: the actual incident-driven body response is connected to fixed-operator outgoing null infinity, with direct input, signed exterior photons, arrived-mass normalization and separated self-GR return. The physical/observational conclusion remains open.'
    write(OUT/'final-result.json',result);print(json.dumps(dict(endpoint=endpoint,audit=result['audit'])),flush=True)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS
    receipt=OUT/f'{action}-receipt.json';assert not receipt.exists(),'Preserve receipts and original verdicts'
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
    incident.native.deadline(CAPS[action]);start=time.monotonic();cpu=time.process_time();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        globals()[action]()
    except Exception as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(action=action,seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
