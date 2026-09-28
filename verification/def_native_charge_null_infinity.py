"""Apply the corrected saved emission to exterior charge and Bondi mass.

Counterexample candidate: the declared first-variation source and discrete
emission, not a nonlinear or microscopic-error certificate.
"""
from pathlib import Path
import hashlib,json,signal,sys,time

ROOT=Path('outputs/direct-eos-gr33')
THERMAL=ROOT/'def-native-conservative-rates/thermal-refined'
UPDATED=THERMAL/'updated-gr-return'
HISTORY=UPDATED/'conserved-history-charge'
APPLIED=HISTORY/'collision-forcing/response/material-charge/finite/zero-exact/compensated/motion-feedback/finite-collision/extended/applied'
OUT=ROOT/'native-charge-null-infinity'

def read(p):return json.loads(p.read_text())
def write(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024**2),b''):h.update(b)
    return h.hexdigest()


def prepare():
    assert not OUT.exists();OUT.mkdir()
    files=[Path(__file__),Path('verification/def_native_global_scalar_closure.py'),Path('verification/verify_native_global_scalar_closure.py'),
           Path('verification/verify_native_updated_gr_return.py'),Path('verification/def_native_characteristic_gr.py'),
           APPLIED/'completed-return/audit.json',APPLIED/'return/gr/combined-source.npz',APPLIED/'return/gr/combined-charge.npz',
           ROOT/'def-native-global-scalar-closure/bound.json',ROOT/'def-native-global-scalar-closure/bound-plan.json']
    for n in [64,128]:files += [UPDATED/f'background/accepted-ports-{n}.npz',THERMAL/f'evolution/source-{n}.npz',APPLIED/f'total-photons/steps-{n}-reference-128.npz']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='d08b523a3',
        previous_turn='Progress: finite collision correction actually reached photons, finite material and compact GR.',
        claim='Propagate the latest actual signed native photon increment and corrected background emission to null infinity, pair exterior stress with body energy debit, and include arrived photon energy in the mass-normalized scalar readout.',
        decision='Does the corrected compact positive endpoint survive the same-solution exterior and mass normalization? Retain full EOS, discarded-matter feedback and nonlinear limitations.',
        reuse='Saved64/128 accepted ports, latest native photon increment, combined compact charge, original4/8 exterior kernel and initial vacuum geometry. No fluid/EOS/GR trajectory replay.',
        emission='Background accepted implicit-stage luminosity is constant over its step, as its owner prescribes. The response uses its actual two SDIRK abscissae and weights as signed packets. Integrate their first/second primitives exactly; never substitute sparse endpoint trapezoids.',
        quadrature='Reuse the existing17 causal cuts and4/8 angular/radial controls. No physical bin/horizon expansion. Any unresolved new emission cuts must appear in controls; do not call these comparisons rigorous continuum bounds.',
        budgets=dict(pilot_seconds=45,production_seconds=120,audit_seconds=60,CPU_threads=1,virtual_GiB=3),
        forecast='Measure one4/4 exterior path including both emission parts. Bound four remaining settings by twice7x that path plus20s setup. Reuse the successful pilot; no dispatch if forecast exceeds120s.',
        gates=dict(port=1e-12,quadrature=.002,time=.02,normalization_identity=1e-12,potential_contraction=1.),
        stop='One pilot and one production. Preserve failures; no automatic resolution, horizon, tolerance or repeated physical path.',
        limits='Fixed initial geometry and finite source/emission representation. Potential enclosure, if used, is conditional on the same frozen coefficients. Full source error, floor feedback, uniform EOS derivatives, nonlinear GR, final charge and observability remain unproved.',
        bindings={(p.relative_to(Path.cwd()) if p.is_absolute() else p).as_posix():sha(p) for p in files}))


def execute(action):
    import resource
    assert action in ['pilot','production']
    dest=OUT/f'{action}.json';failure=OUT/f'{action}-failure.json';assert not dest.exists() and not failure.exists()
    plan=read(OUT/'plan.json')
    for p,h in plan['bindings'].items():assert sha(Path(p))==h,p
    cap=plan['budgets'][action+'_seconds'];started=time.monotonic()
    def timeout(*_):raise TimeoutError('Phase145 '+action+' wall cap')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(cap);resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,)*2)
    try:
        import numpy as np
        import sympy as sp
        import mpmath as mp
        import def_native_global_scalar_closure as old
        import verify_native_updated_gr_return as current
        from types import SimpleNamespace
        C,G,LD=old.C,old.G,np.longdouble
        previous=old.Exterior()
        original_gamma=previous.response.gamma.copy();original_ratio=previous.response.ratio4.copy()
        current.configure();m=old.Exterior()
        assert np.array_equal(original_gamma,m.response.gamma) and np.array_equal(original_ratio,m.response.ratio4)
        binding=read(old.OUT/'bound-plan.json')['bindings'];bgfile=old.flow.INPUT/'balanced-20.npz'
        assert sha(bgfile)==binding[str(bgfile)]
        del previous
        emission={}
        def inputs(n):
            if n in emission:return emission[n]
            with np.load(UPDATED/f'background/accepted-ports-{n}.npz') as z:lum=z['angular_luminosity'].astype(LD);h=float(z['h'])
            with np.load(APPLIED/f'total-photons/steps-{n}-reference-128.npz') as z:
                t=z['accepted_angular_times'];delta=z['accepted_angular_luminosity'].astype(LD);ports=z['radial_ports'];canonical=z['t']
            gamma=1-1/np.sqrt(2);assert lum.shape==(n,4) and delta.shape==(2*n,4)
            assert np.max(abs(t-(h*(np.arange(n)[:,None]+[gamma,1.])).ravel()))<1e-18
            assert abs(h*n-m.T)<1e-18 and np.min(lum)>=0
            w=np.tile([1-gamma,gamma],n).astype(LD)*LD(h);packets=w[:,None]*delta
            q=np.vstack([np.zeros(4,LD),np.cumsum(LD(h)*lum,axis=0)])
            r=np.vstack([np.zeros(4,LD),np.cumsum(LD(h)*q[:-1]+LD(h)**2/2*lum,axis=0)])
            dq=np.vstack([np.zeros(4,LD),np.cumsum(packets,axis=0)])
            dtq=np.vstack([np.zeros(4,LD),np.cumsum(packets*t[:,None],axis=0)])
            aw=np.arange(1,8,2,dtype=LD)/32
            with np.load(THERMAL/f'evolution/source-{n}.npz') as z:baseport=z['outer_cumulative_energy_erg'];bt=z['t']
            err0=float(np.max(abs(q@aw-baseport))/max(abs(baseport[-1]),1.))
            ids=np.array([int(np.argmin(abs(bt-v))) for v in canonical]);assert np.max(abs(bt[ids]-canonical))<1e-18
            err1=float(np.max(abs((dq[2*ids]@aw)-ports[:,1,1]))/max(np.sum(abs(packets)@aw),1.))
            assert max(err0,err1)<1e-12,(err0,err1)
            def primitives(at):
                x=np.clip(np.asarray(at),0,m.T);j=np.clip((x/h).astype(int),0,n-1);s=(x-j*h).astype(LD)[...,None]
                hb=q[j]+s*lum[j];hhb=r[j]+s*q[j]+s*s/2*lum[j]
                k=np.searchsorted(t,x,side='right');hd=dq[k];hhd=x[...,None]*dq[k]-dtq[k]
                return hb,hhb,hd,hhd
            full=primitives(m.T);assert np.max(abs(full[0]-q[-1]))/np.max(q[-1])<1e-14
            emitted=float((q[-1]+dq[-1])@aw);absolute=float(q[-1]@aw+np.sum(abs(packets)@aw))
            row=dict(steps=n,background_port_relative=err0,response_port_relative=err1,emitted_energy_erg=emitted,
                     absolute_emission_energy_erg=absolute,signed_increment_energy_erg=float(dq[-1]@aw))
            emission[n]=SimpleNamespace(primitives=primitives,row=row)
            return emission[n]

        def evaluate(n,a,r):
            begin=time.monotonic();source=inputs(n);k=m.kernel(a,r);mass=[];stress=[];arrival=[]
            for t in m.t:
                hb,hhb,hd,hhd=source.primitives(np.maximum(t-k['delay'],0.));ids=np.arange(len(hb));bins=k['node_bins']
                mass.append([-G/C**3*float(np.sum(k['weights']*k['mass']*v[ids,bins],dtype=LD)) for v in [hhb,hhd]])
                stress.append([-G/(2*C**4)*float(np.sum(k['weights']*k['stress']*v[ids,bins],dtype=LD)) for v in [hb,hd]])
                hb,_,hd,_=source.primitives(np.maximum(t-k['infinity'],0.));ids=np.arange(len(hb))
                arrival.append([float(np.sum(k['mw']*k['mu']*v[ids,k['bins']],dtype=LD)) for v in [hb,hd]])
            mass=np.asarray(mass);stress=np.asarray(stress);arrival=np.asarray(arrival)
            d=dict(t=m.t,mass_U=mass,stress_U=stress,normalized_exterior=-(mass+stress).sum(1)/m.M,
                   arrived_energy_erg=arrival.sum(1),arrived_parts_erg=arrival,normalized_exterior_parts=-(mass+stress)/m.M)
            row=dict(source.row,angular=a,radial=r,seconds=time.monotonic()-begin,kernel_seconds=k['seconds'],
                delay_inverse_seconds=k['inversion_error'],endpoint_exterior=float(d['normalized_exterior'][-1]),
                endpoint_arrived_energy_erg=float(d['arrived_energy_erg'][-1]),endpoint_increment_exterior=float(d['normalized_exterior_parts'][-1,1]))
            label=f'exterior-{n}-a{a}-r{r}';np.savez_compressed(OUT/f'{label}.npz',**d);write(OUT/f'{label}.json',row)
            return d,row

        if action=='pilot':
            # Symbolic numerator/denominator pairing and primitive identities.
            M,K,u,e=sp.symbols('M K u e',nonzero=True)
            assert sp.simplify(((-K-u)/(M-e)+K/M)-((-u/M-K/M*e/M)/(1-e/M)))==0
            x,L,h=sp.symbols('x L h',positive=True)
            assert sp.diff(L*x*x/2,x)==L*x and sp.diff(L*x,x)==L
            write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='Exact normalization difference and constant-stage primitives. Signed packet first/second primitives are cumulative weight and elapsed-time weighted sum. Not a nonlinear evolution theorem.'))
            d,row=evaluate(128,4,4);upper=2*(7*row['seconds']+20)
            result=dict(classification='Counterexample candidate',row=row,upper_seconds=upper,eligible=upper<120,
                same_initial_gamma_and_photon_fourth_moment=True,seconds=time.monotonic()-started)
            write(dest,result);print(json.dumps(result),flush=True);return

        pilot=read(OUT/'pilot.json');assert pilot['eligible']
        paths=[evaluate(n,a,r) for n,a,r in [(128,8,8),(128,4,8),(128,8,4),(64,8,8)]]
        fine=paths[0][0];norm=max(np.max(abs(fine['normalized_exterior'])),1e-300);arrival=max(fine['arrived_energy_erg'][-1],1.)
        controls={key:float(np.max(abs(d['normalized_exterior']-fine['normalized_exterior']))/norm) for key,(d,_) in zip(['angular','radial','time'],paths[1:])}
        arrivals={key:float(np.max(abs(d['arrived_energy_erg']-fine['arrived_energy_erg']))/arrival) for key,(d,_) in zip(['angular','radial','time'],paths[1:])}
        with np.load(APPLIED/'return/gr/combined-charge.npz') as z:
            ids=np.array([int(np.argmin(abs(z['t']-t))) for t in m.t]);assert np.max(abs(z['t'][ids]-m.t))<1e-18
            compact=z['free_scalar'][ids]
        with np.load(APPLIED/'return/gr/combined-source.npz') as z:d=dict(z)
        assert float(d['M_cm'])==m.M and abs(float(d['K_cm'])-m.K)<1e-12*abs(m.K) and m.K<0
        epsilon=G/C**4*fine['arrived_energy_erg']/m.M;alpha0=-m.K/m.M
        scalar=compact+fine['normalized_exterior'];mass_term=alpha0*epsilon
        normalized=(scalar+mass_term)/(1-epsilon)
        assert np.all(epsilon>=0) and max(epsilon)<1
        # Same frozen coefficient proof; rebuild the entire new source norm.
        prior=read(old.OUT/'bound.json');rsp=m.response;rsp.setup(d,8)
        density=rsp.source.reshape(len(rsp.t),-1,8)/rsp.dx.reshape(-1,8)
        co=np.einsum('cij,tcj->tci',rsp.inverse,density);maxima=np.max(abs(co),axis=0)
        mp.iv.dps=40;I=mp.iv.mpf
        def B(v):
            a,b=float(v).as_integer_ratio();return I(a)/I(b)
        def up(v):return float(np.nextafter(float(v.b),np.inf))
        integral=I(0)
        for j,row in enumerate(maxima):
            if rsp.faces[j+1]>=prior['causal_radius_min_cm']:integral+=sum((B(v) for v in row),I(0))*2*B(rsp.half[j])
        source_norm=B(C)*B(m.T)/2*integral;eta=I(prior['global_potential_contraction']);assert up(eta)<1
        r0=B(m.r0);M=B(m.M);K=B(abs(m.K));b0=1-2*M/r0;c0=B(m.N0)*mp.iv.sqrt(b0)
        kappa=M/(r0*b0)+K*K/(2*c0*c0*r0*r0)
        energy=B(inputs(128).row['absolute_emission_energy_erg'])*B(G)/B(C)**4
        es=energy*K*mp.iv.pi/(4*c0*c0*r0*mp.iv.sqrt(1-kappa));em=B(C)*B(m.T)/2*energy*K/(c0*b0*r0*r0)
        potential=eta/(1-eta)*(source_norm+es+em)/M
        uncertainty=up(potential/(1-B(epsilon[-1])))
        interval=[float(np.nextafter(normalized[-1]-uncertainty,-np.inf)),float(np.nextafter(normalized[-1]+uncertainty,np.inf))]
        np.savez_compressed(OUT/'charge.npz',t=m.t,compact=compact,exterior=fine['normalized_exterior'],scalar=scalar,
            arrived_energy_erg=fine['arrived_energy_erg'],epsilon=epsilon,mass_normalization_term=mass_term,normalized=normalized,
            new_source_coefficient_maxima=maxima,source_optical_half=rsp.half)
        passed=max(controls['angular'],controls['radial'],arrivals['angular'],arrivals['radial'])<.002 and max(controls['time'],arrivals['time'])<.02 and interval[0]>0
        result=dict(classification='Counterexample candidate',passed=bool(passed),controls=controls,arrival_controls=arrivals,paths=[r for _,r in paths],
            endpoint_compact=float(compact[-1]),endpoint_exterior=float(fine['normalized_exterior'][-1]),endpoint_scalar=float(scalar[-1]),
            endpoint_mass_normalization_term=float(mass_term[-1]),endpoint_normalized=float(normalized[-1]),epsilon=float(epsilon[-1]),alpha0=alpha0,
            frozen_source_potential_bound=uncertainty,conditional_frozen_source_GR_interval=interval,
            same_initial_gamma_and_photon_fourth_moment=True,actual_current_emission_at_null_infinity=True,
            actual_emission_mass_normalization=True,signed_increment_kept_separate=True,
            potential_scope='Same frozen initial coefficients; recomputed current source polynomial norm and absolute signed-emission energy. No old scalar interval or old source norm reused.',
            limitations='The interval encloses potential repetitions of the declared source only. Quadrature/time comparison, source construction, microscopic EOS, continuous coupling, discarded matter transport/trace and nonlinear corrections are not enclosed.',
            full_floor_feedback_enclosed=False,uniform_EOS_derivative_bound=False,coupled_fixed_point_verified=False,nonlinear_GR=False,final_charge_solved=False,full_goal_complete=False,
            seconds=time.monotonic()-started,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
        write(dest,result);print(json.dumps(result),flush=True);assert passed
    except Exception as exc:
        write(failure,dict(error=repr(exc),seconds=time.monotonic()-started));raise
    finally:signal.alarm(0)


if __name__=='__main__':
    prepare() if sys.argv[1]=='prepare' else execute(sys.argv[1])
