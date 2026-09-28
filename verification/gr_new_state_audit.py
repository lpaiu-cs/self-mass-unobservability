"""Counterexample candidate: independent replay of the completed GR inputs.

This audit binds the new state, reconstructs transport without calling its
implementation, and distinguishes a structure solve from evolution.
"""
from concurrent.futures import ProcessPoolExecutor
import json, sys
import numpy as np
import gr_microphysics as micro
import gr_transport_state as transport

g=micro.g;OUT=g.OUT/'gr-state-audit'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    paths=[g.ROOT/'verification'/n for n in ['gr_new_state_audit.py','direct_eos_gr.py','gr_microphysics.py','gr_transport_state.py']]
    paths += [g.OUT/n for n in ['initial-GR.json','initial-structure-17-4.json','initial-state-17-4.npz','reference-state.npz',
        'gr-microphysics/auxiliaries.json','gr-microphysics/auxiliaries.npz','gr-microphysics/source.json',
        'gr-microphysics/common-EOS-corrected.npz','gr-microphysics/common-EOS-native.npz',
        'gr-microphysics/opacity.json','gr-opacity/result.json','gr-opacity/evaluation.npz',
        'gr-opacity/new-GR-captured.npz','gr-opacity/common-EOS-captured.npz',
        'gr-transport/result.json','gr-transport/diagnostics.npz']]
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='8ab926a',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},processes=4,
        full_EOS_replay='Every new GR midpoint, independent worker order, all 21 outputs bitwise.',
        transport_relative_tolerance=1e-12,conservation_relative_tolerance=1e-12,
        entropy_score_tolerance=1,pressure_log_tolerance=1e-10,
        physical_stability_certified=False,continuous_GR_error_certified=False,full_GR_evolution=False))


def replay(start):
    state,aux=micro.inputs();stop=min(start+256,len(state['X']));eos=g.EOS();score=[]
    ref=dict(np.load(g.OUT/'reference-state.npz'));local=[]
    for i in range(start,stop):
        a=eos(2,state['lnd'][i],state['lnT'][i],state['X'][i])
        assert np.array_equal(a,aux['eos'][i]),i
        pressure=abs(np.log(a[1])-state['logP'][i]);assert pressure<1e-10,(i,pressure)
        H=a[2]+a[1]/a[0];scale=max(2.,32*np.spacing(abs(H)))
        value=abs(np.exp(state['lnT'][i])*(a[3]-ref['s_B'][i]))/scale
        # The actual GR root compares against the table reference, whose s
        # is the same EOS evaluated at the original binary64 rho,T,X.
        score.append(value);local.append(pressure)
    return start,np.array(score),np.array(local)


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    state,aux=micro.inputs();reference=dict(np.load(g.OUT/'reference-state.npz'))
    gr=json.loads((g.OUT/'initial-GR.json').read_text());assert gr['completed']
    for key in ['dm','X']:assert np.array_equal(state[key],reference[key]),key
    assert g.c.sha(g.OUT/'initial-state-17-4.npz')==gr['state_sha256']
    rows=[]
    with ProcessPoolExecutor(max_workers=plan['processes']) as pool:
        for row in pool.map(replay,range(0,len(state['X']),256)):
            rows.append(row);print('GR INDEPENDENT EOS',len(rows),flush=True)
    entropy=np.concatenate([r[1] for r in rows]);pressure=np.concatenate([r[2] for r in rows])
    # Recompute the actual solve's reference values, rather than silently
    # substituting separately saved rounded entropy if they differ.
    tab=dict(np.load(g.OUT/'initial-adiabats-17.npz'))
    assert np.array_equal(tab['reference'][:,3],reference['s_B'])
    assert entropy.max()<=plan['entropy_score_tolerance'],entropy.max()
    diagnostics=dict(np.load(g.OUT/'gr-transport/diagnostics.npz'))
    op=dict(np.load(g.OUT/'gr-opacity/evaluation.npz'))['values'][:,0]
    a=aux['eos'];rho=np.exp(state['lnd']);T=np.exp(state['lnT']);dm=state['dm'];N=np.exp(state['nu'])
    assert np.array_equal(state['CX'],(state['X']/g.c.A)@g.c.W)
    assert np.all(np.diff(state['radius_faces_m'])<0) and np.all(np.diff(state['mass_faces_geom'])<0)
    assert state['radius_faces_m'][-1]==state['mass_faces_geom'][-1]==0
    surface=.5*np.log1p(-2*state['mass_faces_geom'][0]/state['radius_faces_m'][0])
    assert surface==state['nu_faces'][0]
    # Use long double exp and sinh: this is independent of expm1 in the
    # original face routine. It is a finite arithmetic comparison.
    ld=np.longdouble;z=(state['lnT']+state['nu']).astype(ld)
    facekap=(dm[:-1].astype(ld)*op[1:]+dm[1:].astype(ld)*op[:-1])/(dm[:-1]+dm[1:])
    area=4*ld(np.pi)*(state['radius_faces_m'][1:-1].astype(ld)*100)**2
    fourth=2*np.exp(2*(z[1:]+z[:-1]))*np.sinh(2*(z[1:]-z[:-1]))
    flux=area**2*(4*ld(5.670400e-5))*fourth/(3*facekap*np.exp(2*state['nu_faces'][1:-1].astype(ld))*((dm[:-1]+dm[1:])/2))
    original=diagnostics['interior_Linf'];relative=abs(flux-original)/np.maximum(abs(flux),ld(1e-200))
    assert relative.max()<plan['transport_relative_tolerance'],relative.max()
    divergence=np.r_[flux,ld(0)]-np.r_[ld(0),flux]
    conservation=abs(divergence.sum())/abs(divergence).sum()
    assert conservation<plan['conservation_relative_tolerance']
    theta=np.exp(z);entropy_rate=flux*(1/theta[:-1]-1/theta[1:]);assert np.all(entropy_rate>=0)
    # cv*T, a[10], is the derivative of specific energy at fixed density.
    # These are initial frozen-geometry rates, not a time integration.
    heat_rate=divergence/dm.astype(ld)/N.astype(ld)
    logarithmic_T_rate=heat_rate/a[:,10]
    N2=diagnostics['proper_buoyancy_squared'][1];wide=diagnostics['wide_stencil_proper_buoyancy_squared']
    agreed=(N2<0)&(wide<0);proper=np.full(len(dm),np.inf);proper[agreed]=1/np.sqrt(-N2[agreed])
    native=dict(np.load(g.OUT/'gr-opacity/new-GR-captured.npz'));common=dict(np.load(g.OUT/'gr-opacity/common-EOS-captured.npz'))
    assert np.array_equal(native['outputs'],common['outputs'])
    assert np.array_equal(common['used'],aux['electron'])
    result=dict(classification='Counterexample candidate',passed=True,cells=len(dm),
        same_baryon_nuclear_inventories=True,all_21_EOS_outputs_bitwise=True,
        maximum_entropy_inverse_score=float(entropy.max()),maximum_log_pressure_error=float(pressure.max()),
        independent_transport_relative_error=float(relative.max()),closed_energy_relative_residual=float(conservation),
        minimum_pair_entropy_production=float(entropy_rate.min()),
        initial_closed_diffusion_maximum_abs_dlnT_dt=float(abs(logarithmic_T_rate).max()),
        agreed_local_buoyancy_negative_mass_fraction=float(dm@agreed/dm.sum()),
        shortest_agreed_local_buoyancy_proper_seconds=float(proper.min()),
        EOS_electron_replacement_opacity_outputs_bitwise_unchanged=True,
        model_replacement_not_time_evolution=True,physical_EOS_certified=False,
        physical_stability_certified=False,continuous_GR_error_certified=False,full_GR_evolution=False)
    np.savez_compressed(OUT/'audit.npz',entropy_inverse_scores=entropy,pressure_errors=pressure,
        independent_Linf=flux,initial_closed_dlnT_dt=logarithmic_T_rate,
        agreed_local_buoyancy_proper_seconds=proper)
    save('result.json',result);print('GR STATE AUDIT COMPLETE',result,flush=True)
    verify()


def verify():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'result.json').read_text())['passed']
    path=OUT/'manifest.json'
    if not path.exists():
        inputs=[p for folder in [g.OUT/'gr-opacity',micro.OUT,transport.OUT,OUT] for p in folder.rglob('*') if p.is_file()]
        inputs += [g.OUT/n for n in ['initial-GR.json','initial-state-17-4.npz','initial-structure-17-4.json']]
        save('manifest.json',dict(classification='Counterexample candidate',sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in inputs}))
    entries=json.loads(path.read_text())['sha256']
    for rel,digest in entries.items():assert g.c.sha(g.ROOT/rel)==digest,rel
    print('PASS GR STATE AUDIT',len(entries),'artifact SHA',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
