"""Bound the radiative energy acquired by the discarded baryon inventory.

Theorem: a passive tail and its tagged photons cannot gain more than the
incident positive background energy. Numerical application: the saved linear
photon reconstruction, not the unknown fully coupled radiation/fluid solution.
"""
from pathlib import Path
import hashlib,json,shutil,signal,sys,time
import def_native_charge_null_infinity as prior

OUT=prior.ROOT/'native-discarded-radiation'
ATOMIC=Path('/home/lpaiu/work/direct-eos-gr33/photon-hhe-coupled')
read,write,sha=prior.read,prior.write,prior.sha

# Coefficients of the existing GAUNT Laurent polynomials, powers -3,...,3.
# Fortran unsuffixed literals are binary32, promoted to the binary64 owner.
COEFFICIENTS=[
 [0,12.803223,-5.5759888,1.2302628,-2.9094219e-3,7.3993579e-6,-8.7356966e-9],
 [-2.0244141,2.1325684,-1.2709045,1.1595421,-2.0735860e-3,2.7033384e-6,0],
 [-.23387146,.52471924,-.55936432,1.1450949,-1.9366592e-3,2.3572356e-6,0],
 [-5.4418565e-2,.19683564,-.31190730,1.1306695,-1.3482273e-3,-4.6949424e-6,2.3548636e-8],
 [-8.9182854e-3,5.5545091e-2,-.16051018,1.1190904,-1.0401085e-3,-6.9943488e-6,2.8496742e-8],
 [-5.5303574e-3,4.1921183e-2,-.13075417,1.1168376,-8.9466573e-4,-8.8393133e-6,3.4696768e-8],
 [-2.2752881e-3,2.3350812e-2,-9.5441161e-2,1.1128632,-7.4833260e-4,-1.0244504e-5,3.8595771e-8],
 [-9.7200274e-4,1.3298411e-2,-7.1010560e-2,1.1093137,-6.2619148e-4,-1.1342068e-5,4.1477731e-8],
 [-4.9576163e-4,8.5139736e-3,-5.6046560e-2,1.1078717,-5.4837392e-4,-1.2157943e-5,4.3796716e-8],
 [-2.9467046e-4,6.1516856e-3,-4.7326370e-2,1.1052734,-4.4341570e-4,-1.3235905e-5,4.7003140e-8]]


def prepare():
    assert not OUT.exists();OUT.mkdir()
    for name in ['atomic.f','bridge.f90','PARAMS.FOR']:shutil.copyfile(ATOMIC/name,OUT/name)
    paths=[Path(__file__),Path(prior.__file__),prior.OUT/'production.json',prior.OUT/'audit.json',
           prior.ROOT/'def-native-global-scalar-closure/bound.json',
           Path('verification/def_native_material_join.py'),Path('verification/def_native_atmosphere_spectrum.py'),
           Path('verification/verify_native_updated_gr_return.py')]
    for n in [64,128]:paths += [prior.THERMAL/f'evolution/history-{n}.npz',
        prior.UPDATED/f'background/coupled-{n}.npz',prior.APPLIED/f'total-photons/steps-{n}-reference-128.npz']
    paths += list(OUT.glob('*'))
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='971b1a114',
        claim='Replace the old no-subsequent-radiative-energy-growth assumption by a conserved-baryon incident-energy bound, and propagate it into direct scalar, metric stress and photon mass normalization.',
        decision='Does the passive radiative discarded-material uncertainty fit the original2percent comparison scale of the fixed-mass scalar signal? If not, do not enlarge a physical trajectory automatically.',
        representation='Use all original-step discard ledgers and17-knot piecewise-linear saved photon fields, plus the absolute latest native photon increment. Bound every retained atomic level and arbitrary neutral fractions; do not extrapolate the cold EOS below its support.',
        source_scope='Passive reservoir remains within the original compact atmosphere during the saved horizon; its photons are freely outgoing once they cross the original boundary. This is a defined sensitivity subproblem, not a verified trajectory of the deleted fluid. The sampled photon reconstruction is an input envelope, not a uniform bound on the live coupled solution. Pressure-work, changed-photon interactions with retained gas, continuum EOS, and nonlinear GR remain outside.',
        budget=dict(seconds=90,CPU_threads=1,virtual_GiB=3,new_native_states=0,new_trajectories=0),
        forecast='Only two saved photon histories and one existing model construction. No physical integration;90s hard cap includes scientific imports, audit and saved-array reads, excludes initial SHA checks.',
        gates=dict(native_polynomial=1e-12,ledger_monotone=0.,passive_scalar_fraction=.02),
        stop='One bounded readout; preserve failed gate. No additional mesh, clock, horizon, floor or acceptance changes.',
        bindings={(p.relative_to(Path.cwd()) if p.is_absolute() else p).as_posix():sha(p) for p in paths},
        external_inputs={str(ATOMIC/'atomic.so'):sha(ATOMIC/'atomic.so')}))


def run():
    assert not (OUT/'result.json').exists() and not (OUT/'failure.json').exists()
    plan=read(OUT/'plan.json')
    for file,h in dict(plan['bindings'],**plan['external_inputs']).items():assert sha(Path(file))==h,file
    started=time.monotonic()
    def timeout(*_):raise TimeoutError('Phase14690s budget')
    def deadline():
        remaining=plan['budget']['seconds']-(time.monotonic()-started)
        if remaining<=0:timeout()
        # Imported native owners may replace/cancel the process-wide alarm.
        signal.signal(signal.SIGALRM,timeout);signal.setitimer(signal.ITIMER_REAL,remaining)
    deadline()
    import resource
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,)*2)
    try:
        import numpy as np
        import sympy as sp
        import mpmath as mp
        import def_native_hydrogen_exchange as atomic
        import verify_native_updated_gr_return as current
        current.configure();model=current.Model();b=model.bulk;f=model.flow;m=model.m
        deadline()
        C=current.C;G=current.G;LD=np.longdouble;iv=mp.iv;iv.dps=40
        def B(x):
            n,d=float(x).as_integer_ratio();return iv.mpf(n)/iv.mpf(d)
        def up(x):return float(np.nextafter(float(x.b),np.inf))
        def down(x):return float(np.nextafter(float(x.a),-np.inf))
        native=model.spectrum.native;F=float(np.float32(2.99793e14));H=float(np.float32(6.6256e-27))
        sigma=[];owner_errors=[];powers=np.arange(-3,4)
        for n,cc in enumerate(COEFFICIENTS,1):
            cc=np.array(cc,dtype=np.float32).astype(float)
            nu0=native.oldground/n**2/H;x0=nu0/F
            # Each positive term is proportional to nu**(p-3), p<=3.
            # Dropping negative terms therefore bounds all nu>=nu0.
            envelope=B(2.815e29)/n**5/B(nu0)**3*sum(B(v)*B(x0)**int(p) for v,p in zip(cc,powers) if v>0)
            sigma.append(up(envelope))
            nu=np.ascontiguousarray(nu0*np.geomspace(1,1e4,129));bf=np.zeros_like(nu);ff=np.zeros_like(nu)
            native.cross(len(nu),1,n,1.e4,nu,bf,ff)
            polynomial=2.815e29/n**5/nu**3*np.sum(cc[:,None]*(nu[None,:]/F)**powers[:,None],axis=0)
            error=float(np.max(abs(polynomial-bf))/sigma[-1]);owner_errors.append(error)
            assert error<plan['gates']['native_polynomial'] and np.max(bf)<=sigma[-1]*(1+1e-12)
        z=sp.symbols('z',positive=True)
        for p in powers:assert sp.diff(z**int(p-3),z).subs(z,1)<=0
        beta,mu,rho,W,kappa,u=sp.symbols('beta mu rho W kappa u',real=True)
        assert sp.expand(kappa*rho*W*(1-beta*mu)*u-kappa*(rho*W)*(1-beta*mu)*u)==0
        # The restricted normal-nucleus model has <=1 electron per baryon;
        # pairs, lines and omitted opacity channels are not imported here.
        ne_per_mass=B(f.eos.cx)*B(atomic.NA)
        kap_bf=B(f.eos.nH)*B(max(sigma));kap_th=ne_per_mass*B(6.6524587321e-25)
        opacity=kap_bf+kap_th
        geo=read(prior.ROOT/'def-native-global-scalar-closure/bound.json')
        bounds=geo['frozen_polynomial_bounds'];rmin=B(geo['causal_radius_min_cm'])
        phi=B(bounds['phi']['absolute']);amin=B(bounds['lapse']['lower'])*iv.exp(-2*phi*phi)
        r=np.r_[b.d['r'],m.r];edges=np.r_[b.d['edges'],m.rf[1:]]
        assert len(edges)==len(r)+1 and np.all(np.diff(edges)>0)
        causal=edges[1:]>=float(rmin.a)
        lapse=np.r_[b.d['a'],m.a];weights=b.w[:,None]*b.d['num']*b.d['Einf']
        # a<=1 in this saved geometry. With d tau=a dt and Killing energy
        # a E_local, a**2<=1 bounds absorbed Killing power in coordinate time.
        assert max(lapse)<=1
        source=np.load(prior.APPLIED/'return/gr/combined-source.npz')
        M=B(float(source['M_cm']));K=B(abs(float(source['K_cm'])));T=B(float(source['t'][-1]))
        nominal=read(prior.OUT/'production.json');target=abs(nominal['endpoint_scalar'])
        alpha0=B(nominal['alpha0']);Phi=B(bounds['Phi']['absolute']);eta=B(geo['global_potential_contraction'])
        cc=B(C);gg=B(G);r0=B(float(edges[-1]));N0=B(float(lapse[-1]));b0=1-2*M/r0;c0=N0*iv.sqrt(b0)
        ray_kappa=M/(r0*b0)+K*K/(2*c0*c0*r0*r0)
        baseline_epsilon=B(nominal['epsilon']);rows=[]
        for n in [64,128]:
            deadline()
            with np.load(prior.UPDATED/f'background/coupled-{n}.npz') as bg:
                t=bg['snapshot_t'];I=bg['snapshot_I'].sum(1);bulk=bg['snapshot_bulk_I']
            cap=np.zeros(len(t));basecap=np.zeros(len(t))
            with np.load(prior.APPLIED/f'total-photons/steps-{n}-reference-128.npz') as delta:
                assert np.max(abs(delta['t']-t))<1e-18
                dI=delta['photon_history_scaled_occupation']
            for j in range(len(t)):
                background=np.concatenate([bulk[j],I[j]])
                assert background.min()>=0
                density=np.sum(background*weights,axis=(1,2))/lapse**4
                added=np.sum(abs(dI[j])*weights,axis=(1,2))/lapse**4
                basecap[j]=np.max(density[causal]);cap[j]=np.max((density+added)[causal])
            with np.load(prior.THERMAL/f'evolution/history-{n}.npz') as hist:
                tt=hist['t'];discard=hist['discard'].astype(LD)
            mass=discard[:,0]*LD(model.gas_scale)/LD(C)**2
            energy=(discard[:,2]+LD(m.a0)*LD(model.cx)*discard[:,0])*LD(model.gas_scale)
            assert mass.min()>=0 and energy.min()>=0 and np.diff(mass).min()>=0 and np.diff(energy).min()>=0
            # Deletion times within each original step are not stored. Put
            # the entire end-of-step inventory at its start for an upper bound.
            cap_at=np.interp(tt,t,cap);interval_cap=np.maximum(cap_at[:-1],cap_at[1:])
            incident=iv.mpf(0)
            for j,ucap in enumerate(interval_cap):
                incident+=2*cc*opacity*B(float(mass[j+1]))*B(ucap)*B(float(tt[j+1]-tt[j]))
            A=incident;F0=B(float(energy[-1]));available=F0+A
            # Using the whole initial Killing energy as releasable is loose
            # but does not assume an EOS energy zero or remove binding energy.
            photon_L1=F0+2*A
            direct=gg/(2*cc**3)*T*(4*phi)/rmin*2*available/amin/M
            stress=gg/(2*cc**3)*T*Phi*(available+photon_L1)/amin/M
            e=gg/cc**4*photon_L1
            external=(e*K*iv.pi/(4*c0*c0*r0*iv.sqrt(1-ray_kappa))+cc*T/2*e*K/(c0*b0*r0*r0))/M
            mass_error=cc*T/2*B(geo['mass_source_coefficient_absolute_per_cm2'])*(r0-rmin)/(B(bounds['lapse']['lower'])*iv.sqrt(1-2*B(bounds['mass']['upper'])/rmin))*gg/cc**4/amin*photon_L1/M
            fixed=(direct+stress+external+mass_error)/(1-eta)
            eps=gg/cc**4*photon_L1/M
            normalized=fixed/(1-baseline_epsilon-eps)+(alpha0+B(target))*eps/((1-baseline_epsilon-eps)*(1-baseline_epsilon))
            np.savez_compressed(OUT/f'inputs-{n}.npz',t=t,base_radiation_cap=basecap,radiation_cap=cap,original_t=tt,
                baryon_discard_g=mass,Killing_discard_erg=energy,interval_radiation_cap=interval_cap,causal_cells=causal)
            rows.append(dict(steps=n,discarded_baryon_g=float(mass[-1]),discarded_Killing_erg=float(energy[-1]),
                radiation_cap_erg_cm3=float(max(cap)),background_cap_erg_cm3=float(max(basecap)),
                incident_background_energy_bound_erg=up(A),available_tail_energy_bound_erg=up(available),
                photon_difference_L1_bound_erg=up(photon_L1),direct_scalar_bound=up(direct),
                metric_stress_bound=up(stress),exterior_bound=up(external),mass_constraint_bound=up(mass_error),
                fixed_mass_scalar_bound=up(fixed),normalized_charge_bound=up(normalized),
                fixed_scalar_fraction=up(fixed)/target,normalized_over_scalar=up(normalized)/target,
                passive_normalized_interval=[float(np.nextafter(nominal['endpoint_normalized']-up(normalized),-np.inf)),
                                             float(np.nextafter(nominal['endpoint_normalized']+up(normalized),np.inf))]))
            del I,bulk,dI
        deadline()
        result=dict(classification='Counterexample candidate',passed=all(r['normalized_over_scalar']<.02 for r in rows),
            paths=rows,level_cross_section_bounds_cm2=sigma,native_owner_polynomial_relative=max(owner_errors),
            opacity_bound_cm2_g=up(opacity),bound_free_opacity_bound_cm2_g=up(kap_bf),Thomson_opacity_bound_cm2_g=up(kap_th),
            symbolic_checks=True,all_original_discard_steps_used=True,conserved_baryon_opacity_bound=True,
            no_subsequent_radiative_energy_growth_assumption_removed=True,
            scope='Passive compact-atmosphere reservoir plus tagged photons under the declared17-knot positive background radiation envelope; freely outgoing exterior photons; same restricted10-level H and Thomson model and frozen first-variation geometry. Reservoir confinement is assumed, not verified for the deleted fluid. No density/temperature/y extrapolation. Photon change is counted in stress, exterior and mass normalization. Budget inequalities are analytic; applying the saved caps to an unsaved coupled history is not proved.',
            limitations='Arbitrary mechanical work from retained gas and amplification of changed photons by that gas are not bounded. The saved cap is not a uniform source/EOS derivative or continuum error certificate. Other opacity channels, nonlinear gravity, static nonabsorbability and observation remain open.',
            full_floor_feedback_enclosed=False,full_source_error_enclosed=False,coupled_fixed_point_verified=False,
            nonlinear_GR=False,final_charge_solved=False,full_goal_complete=False,
            seconds=time.monotonic()-started,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
        write(OUT/'result.json',result);print(json.dumps(result),flush=True)
    except Exception as exc:
        write(OUT/'failure.json',dict(error=repr(exc),seconds=time.monotonic()-started));raise
    finally:signal.alarm(0)


if __name__=='__main__':globals()[sys.argv[1]]()
