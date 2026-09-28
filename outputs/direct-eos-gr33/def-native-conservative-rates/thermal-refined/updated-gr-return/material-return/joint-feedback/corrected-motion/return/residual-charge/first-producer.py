"""Conserved E/H mismatch -> actual native stress -> retarded charge enclosure.

The certified box is on stored linear source coordinates, not an EOS or
waveform fixed-point error certificate. No new physical trajectory is run.
"""
from pathlib import Path
import json,resource,signal,sys,time
import numpy as np
import sympy as sp
import verify_native_corrected_motion_feedback as run
import verify_native_material_response as native
import verify_native_stage_energy_charge as readout

OUT=run.OUT/'residual-charge';write=run.write;sha=run.sha
C=run.prior.C;G=readout.G;AMP=run.physical.AMP


def prepare():
    assert not OUT.exists();OUT.mkdir()
    assert json.loads((run.OUT/'audit.json').read_text())['passed']
    files=[Path(__file__),Path(run.__file__),Path(native.__file__),Path(readout.__file__),
           Path(readout.gr.__file__),Path(readout.gr.base.__file__),run.OUT/'audit.json']
    for n,r in run.PATHS:
        files.extend([run.OUT/f'steps-{n}-reference-{r}.npz',run.photon_path(n,r),run.GR/f'source-{n}-reference-{r}.npz'])
    write(OUT/'plan.json',dict(classification='Counterexample candidate',before_checkpoint='310e6e67f',
        claim='Propagate the actual saved photon-minus-free-material E/H mismatch through the same conservative native pressure map and retarded GR. Bound every signed source in a cellwise E/H box, including mass/stress and all compact potential repeats of the declared linear reconstruction.',
        decision='Determine the measured charge shift, a sign-independent conditional source bound, and the amplification required for this residual scale to threaten the existing primary charge. Do not confuse this input-box bound with an actual fixed-point error or an EOS derivative certificate.',
        source_box='Baryon,momentum,metric,photon state and ports held fixed; each material reference-E/H deviation is bounded by its per-cell maximum stored mismatch over17 common times. Native conservative pressure coefficients are held at each saved background knot; derived GR sources are interpolated linearly in time and with the existing spatial polynomial. This explicitly declares a stored-coefficient uncertainty model.',
        reuse='Completed Phase137 actual paths; same531 cells,17 knots,three64/128 cases and4/8 GR quadrature. Reuse the shared native pressure owner and one-outer-observer characteristic evaluator. No new EOS bank,fluid/photon steps,grid,horizon or waveform iteration.',
        budgets=dict(pilot_seconds=45,sources_seconds=120,readout_seconds=45,enclosure_seconds=90,CPU_threads=1,virtual_GiB=3),
        forecast='Phase137 three independent full source readouts took17.57s; this uses two shared background models and two pressure basis calls per knot, plus one direct combined check. Pilot one midpoint to require2x extrapolated remaining source work within120s. Readout uses one observer, not532 field locations.',
        gates=dict(pressure_probe=.002,linear_identity=1e-10,independent=1e-9,quadrature=.002,box_contains_actual=True,potential=1.,inactive_fraction=1e-12),
        limits='No coupling gain,uniform physical EOS/derivative error,continuous-time residual,nonlinear GR,exterior/floor feedback or final charge is certified. An arbitrary-sign stored source box need not itself solve the material equations; it is an overinclusive readout uncertainty set.',
        stop='Stop on a failed mapping,forecast,gate or cap. Preserve completed maps. No automatic new trajectory,larger budget or finer resolution.',
        bindings={str(p):sha(p) for p in files}))
    a,b,E,H=sp.symbols('a b E H',real=True);pr,pt=sp.symbols('pr pt',real=True)
    assert sp.expand((E-pr-2*pt)+pr+2*pt-E)==0
    assert sp.expand(a*(E+H)+b*(E-H)-(a+b)*E-(a-b)*H)==0
    from fractions import Fraction
    values=[0.,float(np.nextafter(0.,1.)),1e-200,1e-90,1e-30,1.,1e30,1e90]
    for x in values:
        for y in values:
            assert Fraction(float(plus(x,y)))>=Fraction(x)+Fraction(y)
            assert Fraction(float(times(x,y)))>=Fraction(x)*Fraction(y)
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        scope='Linear conserved-coordinate pressure superposition and trace=energy-radial_pressure-2*tangential_pressure. Triangle inequality then bounds any input box and eta<1 gives B/(1-eta). Positive outward arithmetic also checked against exact rational binary values, including underflow. Physical EOS derivatives and dynamical amplification are separate premises.'))


def pressure_map(m,k):
    z=np.zeros((4,m.n));field=np.zeros((3,m.n));maps=[];errors=[]
    for axis,scale in [(2,1e30),(3,1e38)]:
        z[:]=0.;z[axis]=scale
        _,p,e=native.pressure(m,k,z,field);maps.append(p/scale)
        errors.append(float(np.sum(abs(e))/max(np.sum(abs(p)),1e-300)))
    return np.stack(maps,axis=-1),max(errors)


def pilot():
    assert not (OUT/'pilot.json').exists();run.physical.configure()
    signal.signal(signal.SIGALRM,run.physical.branch.base.flow.old.optical.timeout);signal.alarm(45)
    start=time.monotonic();m=run.Material(128,128);setup=time.monotonic()-start
    start_map=time.monotonic();coeff,error=pressure_map(m,8);seconds_map=time.monotonic()-start_map
    np.savez_compressed(OUT/'map-pilot.npz',coeff=coeff)
    forecast=33*seconds_map+2*setup+10
    row=dict(classification='Counterexample candidate',seconds=time.monotonic()-start,setup_seconds=setup,map_seconds=seconds_map,
        pressure_probe=error,forecast_seconds=forecast,upper_seconds=2*forecast,eligible=error<.002 and 2*forecast<120)
    write(OUT/'pilot.json',row);signal.alarm(0);print(json.dumps(row),flush=True)


def sources():
    p=json.loads((OUT/'pilot.json').read_text());assert p['eligible'] and not (OUT/'sources.json').exists()
    for f,h in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert sha(f)==h,f
    run.physical.configure();signal.signal(signal.SIGALRM,run.physical.branch.base.flow.old.optical.timeout);signal.alarm(120)
    start=time.monotonic();maps={};rows=[]
    for reference in [128,64]:
        m=run.Material(reference,128);coeff=[];probe=0.
        for k in range(17):
            if reference==128 and k==8:
                c=np.load(OUT/'map-pilot.npz')['coeff'];e=p['pressure_probe']
            else:c,e=pressure_map(m,k)
            coeff.append(c);probe=max(probe,e)
        maps[reference]=(m,np.array(coeff),probe)
        np.savez_compressed(OUT/f'pressure-map-{reference}.npz',coeff=coeff,t=m.t,reference_lapse=m.a)
    for n,r in run.PATHS:
        m,coeff,probe=maps[r];d=np.load(run.OUT/f'steps-{n}-reference-{r}.npz');p=np.load(run.photon_path(n,r))
        ii=[int(np.argmin(abs(d['t']-t))) for t in m.t];jj=[int(np.argmin(abs(p['t']-t))) for t in m.t]
        assert np.max(abs(d['t'][ii]-m.t))<1e-18 and np.max(abs(p['t'][jj]-m.t))<1e-18
        actual=d['history_scaled'][ii][:,[2,3]]*AMP;residual=p['moments'][jj][:,[1,2]]-actual
        box=np.nextafter(np.max(abs(residual),axis=0),np.inf)
        pressure=np.einsum('tpnc,tcn->tpn',coeff,residual);energy=residual[:,0]/m.a
        # A combined native readout checks units and order independently of the
        # two basis evaluations. Fixed B/P makes this pressure map cell-local.
        z=np.zeros((4,m.n));z[[2,3]]=residual[8]/AMP
        _,combined,err=native.pressure(m,8,z,np.zeros((3,m.n)))
        identity=float(np.sum(abs(combined*AMP-pressure[8]))/max(np.sum(abs(pressure[8])),1e-300))
        inactive=max(float(np.sum(abs(residual[k,:,~m.point(k)['active']]))/max(np.sum(abs(residual[k])),1e-300)) for k in range(17))
        base=dict(np.load(run.GR/f'source-{n}-reference-{r}.npz'))
        for key in ['baryon_g','photon_energy_erg','photon_radial_pressure_erg','inner_cumulative_energy_erg','outer_cumulative_energy_erg','inner_luminosity','outer_luminosity']:
            if key in base:base[key]=np.zeros_like(base[key])
        base.update(gas_nonrest_energy_erg=energy,nonrest_trace_erg=energy-pressure[:,1]-2*pressure[:,0],
                    nonrest_stress_erg=energy-pressure[:,1],pressure_volume_erg=pressure[:,0],metric_stress_erg=energy-pressure[:,1])
        label=f'{n}-{r}';np.savez_compressed(OUT/f'source-{label}.npz',**base)
        np.savez_compressed(OUT/f'residual-{label}.npz',t=m.t,residual=residual,box=box,pressure=pressure,energy=energy)
        row=dict(steps=n,reference=r,pressure_probe=probe,superposition_relative=identity,inactive_residual_fraction=inactive,
            absolute_E_H_box_sum=np.sum(box,axis=1).tolist(),passed=probe<.002 and identity<1e-10 and inactive<1e-12)
        rows.append(row);assert row['passed'],row
    result=dict(classification='Counterexample candidate',passed=True,paths=rows,seconds=time.monotonic()-start)
    write(OUT/'sources.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


def fields():
    assert json.loads((OUT/'sources.json').read_text())['passed'];assert not (OUT/'readout.json').exists()
    run.physical.configure();signal.signal(signal.SIGALRM,run.physical.branch.base.flow.old.optical.timeout);signal.alarm(45)
    start=time.monotonic();m=run.GRResponse();rows=[]
    for n,r,order in [(128,128,8),(128,128,4),(64,128,8),(128,64,8)]:
        d=dict(np.load(OUT/f'source-{n}-{r}.npz'));wave=readout.read(m,d,order)
        np.savez_compressed(OUT/f'wave-{n}-{r}-g{order}.npz',**wave)
        direct,error=readout.independent.direct(m,d,order);value=float(wave['direct_scalar'][-1])
        norm=max(float(np.max(abs(wave['direct_scalar']))),1e-300)
        agreement=abs(direct-value)/norm
        rows.append(dict(steps=n,reference=r,order=order,endpoint_free=float(wave['free_scalar'][-1]),endpoint_direct=value,independent_relative=agreement,inverse_radius_error=error))
    a=np.load(OUT/'wave-128-128-g8.npz')['free_scalar'];b=np.load(OUT/'wave-128-128-g4.npz')['free_scalar']
    quad=float(np.max(abs(a-b))/max(np.max(abs(a)),1e-300));assert quad<.002 and max(r['independent_relative'] for r in rows)<1e-9
    result=dict(classification='Counterexample candidate',passed=True,paths=rows,quadrature_relative=quad,seconds=time.monotonic()-start,
        interpretation='Actual signed residual source applied to retarded compact free response. Other paths measure changed residuals; they are not a time-convergence certificate of a fixed residual forcing.')
    write(OUT/'readout.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


def up(x):return np.nextafter(np.asarray(x,float),np.inf)
def plus(a,b):return up(np.asarray(a)+b)
def times(a,b):return up(np.asarray(a)*b)
def total(a,axis=0):
    x=np.moveaxis(np.asarray(a),axis,0);s=np.zeros(x.shape[1:])
    for v in x:s=plus(s,v)
    return s


def box_bound(model,d,coeff,box,lapse):
    """Outward nonnegative arithmetic for frozen binary coefficient maps."""
    model.setup(d,8);ids=model.ids;z=model.z
    energy=up(box[0]/lapse)
    pp=total(times(abs(coeff),box.T[None,None]),axis=3)
    trace=np.max(plus(plus(energy,pp[:,1]),times(2.,pp[:,0])),axis=0)
    stress=np.max(plus(energy,pp[:,1]),axis=0)
    q=readout.flow.initial.Quadrature(d['edges'],8)
    w,_,a,B,_=model.geo(q.r.ravel()-model.model.m.RJ)
    measure=q.r*q.r*B.reshape(q.r.shape);norm=measure@q.w
    mean=((measure*a.reshape(q.r.shape))@q.w)/norm
    partial=((measure*a.reshape(q.r.shape))@q.Q.T)/norm[:,None]
    killing=times(energy,abs(mean));prefix=np.zeros(len(energy));s=0.
    for j,v in enumerate(killing):prefix[j]=s;s=plus(s,v)
    mass_co=G/C**4*np.sqrt(z['b'])/z['lapse']
    jb=times(abs(mass_co),plus(prefix[ids],times(energy[ids],abs(partial.ravel()))))
    direct_co=-G/C**4*model.weights*w
    stress_co=-G/C**4*model.weights*a*z['Phi']
    mass_source_co=model.dx*z['K']
    nodes=plus(plus(times(abs(direct_co),trace[ids]),times(abs(stress_co),stress[ids])),times(abs(mass_source_co),jb))
    density=up(nodes.reshape(-1,8)/model.dx.reshape(-1,8))
    spatial=total(times(abs(model.inverse),density[:,None,:]),axis=2)
    integrals=times(total(spatial,axis=1),times(2.,model.half))
    free=times(times(float(C),float(model.t[-1]))/2,total(integrals))
    # The potential's cell field and time interpolation have norm <=1.
    # Enclose its spatial interpolation too, rather than using nodal quadrature
    # as a continuum bound for a sign-changing polynomial.
    vp=total(times(abs(model.inverse),abs(z['V']).reshape(-1,1,8)),axis=2)
    vint=total(times(total(vp,axis=1),times(2.,model.half)))
    eta=times(times(float(C),float(model.t[-1]))/2,vint)
    denominator=np.nextafter(1.-eta,-np.inf);assert denominator>0
    normalized=up(up(free/denominator)/float(d['M_cm']))
    potential_only=up(times(eta,normalized))
    actual=model.source.reshape(len(model.t),-1,8)/model.dx.reshape(-1,8)
    actual_co=np.einsum('cij,tcj->tci',model.inverse,actual)
    ratio=float(np.max(np.divide(abs(actual_co),spatial[None],out=np.zeros_like(actual_co),where=spatial[None]>0)))
    assert ratio<=1 and np.isfinite(normalized)
    return dict(free_norm_cm=float(free),compact_potential_contraction=float(eta),
        normalized_box_bound=float(normalized),normalized_potential_remainder=float(potential_only),
        actual_source_coefficient_fraction=ratio),dict(pressure_coeff=coeff,box=box,reference_lapse=lapse,
        source_coefficient_bounds=spatial,potential_coefficient_bounds=vp,cell_optical_half=model.half,
        frozen_inverse=model.inverse,frozen_direct_coefficient=direct_co,frozen_stress_coefficient=stress_co,
        frozen_mass_source_coefficient=mass_source_co,frozen_mass_coefficient=mass_co,frozen_mean_lapse=mean,
        frozen_partial_lapse=partial,frozen_potential=z['V'],frozen_dx=model.dx,actual_source_coefficients=actual_co)


def enclose():
    assert json.loads((OUT/'readout.json').read_text())['passed'];assert not (OUT/'enclosure.json').exists()
    run.physical.configure();signal.signal(signal.SIGALRM,run.physical.branch.base.flow.old.optical.timeout);signal.alarm(90)
    start=time.monotonic();m=run.GRResponse();rows=[]
    for n,r in run.PATHS:
        d=dict(np.load(OUT/f'source-{n}-{r}.npz'));p=np.load(OUT/f'pressure-map-{r}.npz');e=np.load(OUT/f'residual-{n}-{r}.npz')
        assert np.all(abs(e['residual'])<=e['box'][None])
        row,inputs=box_bound(m,d,p['coeff'],e['box'],p['reference_lapse'])
        wave=np.load(OUT/f'wave-{n}-{r}-g8.npz');nominal=float(wave['free_scalar'][-1])
        assert np.max(abs(wave['free_scalar']))<row['normalized_box_bound']
        row.update(steps=n,reference=r,actual_signed_free_endpoint=nominal,
            signed_with_compact_potential_interval=[float(np.nextafter(nominal-row['normalized_potential_remainder'],-np.inf)),float(np.nextafter(nominal+row['normalized_potential_remainder'],np.inf))])
        rows.append(row);np.savez_compressed(OUT/f'enclosure-inputs-{n}-{r}.npz',**inputs)
    # A fixed linear source is enclosed exactly as declared; interpolation and
    # native coefficient generation errors are explicitly outside that theorem.
    primary=json.loads((run.prior.run.physical.GR/'result.json').read_text())
    fine=rows[1];lower=primary['conditional_scalar_interval'][0]
    threshold=float(np.nextafter(lower/fine['normalized_box_bound'],-np.inf))
    for path,h in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert sha(path)==h,path
    result=dict(classification='Counterexample candidate',passed=True,paths=rows,seconds=time.monotonic()-start,
        primary_conditional_lower=lower,residual_box_over_primary_lower=fine['normalized_box_bound']/lower,
        maximum_E_H_box_amplification_before_primary_lower_can_be_erased=threshold,
        actual_residual_applied_to_native_stress_and_GR=True,declared_source_box_and_compact_potential_enclosed=True,
        theorem_scope='Frozen binary conservative-pressure and GR coefficient maps. Arbitrary signed E/H inputs within each saved cell box at each knot, and the declared linear time/spatial source reconstruction. Outward binary64 arithmetic encloses nonnegative sums/products/divisions; the compact Volterra norm is below1, giving B/(1-eta). This excludes coefficient-generation and physical interpolation errors.',
        amplification_scope='If the eventual E/H-induced source error is contained in A times this exact box, its charge effect is at most A times this bound. The threshold compares to the old primary conditional lower only. No measured iteration ratio is assigned to A; no assertion that all other coupled channels share this box.',
        EOS_derivative_error_enclosed=False,coupled_fixed_point_verified=False,continuous_residual_enclosed=False,
        exterior_floor_feedback_enclosed=False,nonlinear_GR=False,final_charge_solved=False,full_goal_complete=False)
    write(OUT/'enclosure.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


if __name__=='__main__':
    cap=3*1024**3;resource.setrlimit(resource.RLIMIT_AS,(cap,cap));globals()[sys.argv[1]]()
