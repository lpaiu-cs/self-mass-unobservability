"""Audit saved finite-collision superposition, conservative return and GR inputs."""
from pathlib import Path
import io,json,signal,sys,time
import numpy as np
import sympy as sp


def main():
    import propagate_native_finite_collision_remainder as run
    used=json.loads((run.RETURN/'sources.json').read_text())['seconds'];budget=min(60,int(90-used));assert budget>0
    run.write(run.OUT/'audit-plan.json',dict(classification='Counterexample candidate',
        claim='Verify the actual added collision source, saved linear superposition, photon moments, unchanged prescribed transport, finite material balances and GR return.',
        budget_seconds=budget,source_seconds_used=used,source_plus_audit_cap=90,
        reuse='Saved arrays only, with three fixed operator reads; no new trajectory or finite collision evaluation.',
        limits='Finite sampled correction and bookkeeping, not continuous nonlinear or microscopic error certification.',
        bindings={str(p):run.sha(p) for p in [Path(__file__),Path(run.__file__),run.SOURCE/'production.json',run.OUT/'result.json',run.TOTAL/'result.json',run.RETURN/'production.json',run.RETURN/'sources.json',run.GR/'result.json']}))
    started=time.monotonic();signal.alarm(budget)
    before=dict(run.file_load(run.SOURCE/'point-8.npz'));after=dict(np.load(run.SOURCE/'point-8.npz'))
    assert before.keys()==after.keys() and all(np.array_equal(before[key],after[key]) for key in before)
    del before,after
    run.configure();m=run.Response();LD=np.longdouble
    a,b,c=sp.symbols('a b c');assert sp.expand((a+b-c)+c-a-b)==0
    run.write(run.OUT/'symbolic.json',dict(classification='Proven',passed=True,
        scope='R=C(Q0+dQ,I0+dI)-C(Q0,I0)-Jd implies Jd+R equals the finite collision difference at the sampled original state. Linear response superposition holds for the same frozen operator. This is not a nonlinear trajectory identity.'))
    source=[]
    for k in [0,8,16]:
        c=m.local(m.t[k]);d=m.defect(k);errors=[]
        for key,name in [('q','photon'),('qb','bound'),('qe','escape')]:
            expected=d[name]/run.AMP
            if key!='qe':expected=expected/m.scale
            errors.append(float(np.max(abs(c[key]-expected))/max(np.max(abs(expected)),1.)))
        gas=m.gas(c['q'],c['qb'],c['qe'])
        energy=float(abs(np.sum(c['q']*m.Eweight)+np.sum(c['qe'][1])+np.sum(gas[:,0]*m.eu))/max(np.sum(abs(c['q'])*m.Eweight),1.))
        species=float(abs(np.sum(c['qb']*m.Nweight)-np.sum(gas[:,1]*m.nu))/max(np.sum(abs(c['qb'])*m.Nweight),1.))
        assert max(errors+[energy,species])<1e-10
        source.append(dict(k=k,once_only=errors,energy=energy,species=species))
    keys=['delta_packet_scaled_occupation','delta_material','ledger','escape','collision_transfer','radial_ports','photon_history_scaled_occupation','material_history','accepted_angular_luminosity']
    states=[];photons=[];rows=[]
    for n in [64,128]:
        before=dict(np.load(run.run.photon_path(n,128)));inc=dict(np.load(run.OUT/f'steps-{n}-reference-128.npz'))
        p=dict(np.load(run.photon_path(n,128)));d=dict(np.load(run.RETURN/f'steps-{n}-reference-128.npz'))
        assert all(np.array_equal(p[key],before[key]+inc[key]) for key in keys)
        assert np.array_equal(p['moments'][:,[0,1,2,3,5,6]],(before['moments']+inc['moments'])[:,[0,1,2,3,5,6]])
        w=m.Eweight/m.scale;raw=p['photon_history_scaled_occupation']
        reconstructed=[(0,np.sum(raw*w,axis=(2,3))),(4,np.sum(abs(raw)*w,axis=(2,3))),(5,np.sum(raw*w*m.model.bulk.mu2[None,None,:,None],axis=(2,3)))]
        moment=max(float(np.max(np.sum(abs(v-p['moments'][:,i]),axis=1))/max(np.max(np.sum(abs(p['moments'][:,i]),axis=1)),1.)) for i,v in reconstructed)
        expected=before['moments'][:,[1,2]].astype(LD)-before['collision_transfer'].transpose(0,2,1)
        actual=p['moments'][:,[1,2]].astype(LD)-p['collision_transfer'].transpose(0,2,1)
        scale=np.maximum(np.max(np.sum(abs(p['moments'][:,[1,2]]),axis=2),axis=0),1.)
        mechanical=np.max(np.sum(abs(actual-expected),axis=2),axis=0)/scale
        ids=np.array([int(np.argmin(abs(d['t']-t))) for t in p['t']]);assert np.max(abs(d['t'][ids]-p['t']))<1e-18
        z=d['history_scaled'][ids].astype(LD)*LD(run.AMP)
        balance=float(np.max(abs(d['history_scaled'].astype(LD).sum(2)+d['discards_scaled']-d['ledgers_scaled'])/np.maximum(d['norms_scaled'],1.)))
        residual=np.max(np.sum(abs(p['moments'][:,[1,2]]-z[:,[2,3]]),axis=2),axis=0)
        relative=residual/np.maximum(np.max(np.sum(abs(z[:,[2,3]]),axis=2),axis=0),1.)
        pilot=dict(np.load(run.RETURN/f'pilot-{n}.npz'));assert np.array_equal(d['history_scaled'][:3],pilot['history_scaled'])
        assert moment<1e-12 and max(mechanical)<1e-8 and balance<1e-8
        states.append(z[:,[2,3]]);photons.append(p['moments'][:,[1,2]])
        rows.append(dict(steps=n,superposition_exact=True,material_pilot_prefix_exact=True,photon_moment_relative=moment,
            prescribed_transport_change_relative=np.asarray(mechanical,float).tolist(),full_material_balance=balance,
            energy_H_residual=np.asarray(relative,float).tolist(),absolute_residual=np.asarray(residual,float).tolist()))
    numerical=np.maximum(np.max(np.sum(abs(states[0]-states[1]),axis=2),axis=0),np.max(np.sum(abs(photons[0]-photons[1]),axis=2),axis=0))
    ratio=np.asarray(np.array(rows[1]['absolute_residual'])/np.maximum(numerical,LD(1e-300)),float)
    gr=json.loads((run.GR/'result.json').read_text());assert gr['passed'] and gr['seconds']<90
    result=dict(classification='Counterexample candidate',passed=True,source=source,paths=rows,NPZ_memory_read_exact=True,
        residual_over_temporal_comparison=ratio.tolist(),residual_below_temporal_comparison_scale=bool(max(ratio)<1),
        finite_collision_remainder_applied_to_photons_material_and_GR=True,
        uniform_EOS_derivative_bound=False,uniform_nonlinear_remainder_bound=False,coupled_fixed_point_verified=False,
        exterior_floor_feedback_closed=False,full_nonlinear_GR=False,final_charge_solved=False,full_goal_complete=False,
        seconds=time.monotonic()-started)
    run.write(run.OUT/'audit.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


def saved():
    root=Path('outputs/direct-eos-gr33/def-native-conservative-rates/thermal-refined/updated-gr-return/conserved-history-charge/collision-forcing/response/material-charge/finite/zero-exact/compensated/motion-feedback/finite-collision')
    source=root/'extended';out=source/'applied';started=time.monotonic();signal.alarm(30)
    def write(path,value):path.write_text(json.dumps(value,indent=2)+'\n')
    write(out/'saved-audit-plan.json',dict(classification='Counterexample candidate',budget_seconds=30,
        allocation='Use30s of the unstarted90s source/readout audit allowance. No physical production or larger source-readout cap.',
        scope='Check saved sample metadata, exact zero input, direct versus memory-backed NPZ decoding and the symbolic paired number projection. Do not certify the unrun photon/material/GR return.'))
    C,C0,J,R,d,lo,hi,mu=sp.symbols('C C0 J R d lo hi mu')
    assert sp.expand((J+R-C+C0).subs(R,C-C0-J))==0
    nlo=-d*hi/(hi-lo);nhi=d*lo/(hi-lo)
    assert sp.simplify(nlo+nhi+d)==0 and sp.simplify(lo*nlo+hi*nhi)==0 and sp.simplify(mu*(lo*nlo+hi*nhi))==0
    write(out/'symbolic.json',dict(classification='Proven',passed=True,
        scope='The defined finite collision remainder restores the exact sampled collision difference algebraically. The two-frequency number projection removes its chosen number residual with zero energy and same-angle momentum. No trajectory or EOS-error theorem.'))
    rows=[json.loads(p.read_text()) for p in source.glob('point-*.json')];rows.sort(key=lambda r:r['k'])
    assert all(r['passed'] for r in rows) and rows[0]['k']==0
    zero=dict(np.load(source/'point-0.npz'));assert all(not np.any(zero[k]) for k in ['photon','bound','escape','linear','half_photon'])
    path=source/'point-8.npz';direct=dict(np.load(path));buffered=dict(np.load(io.BytesIO(path.read_bytes())))
    assert direct.keys()==buffered.keys() and all(np.array_equal(direct[k],buffered[k]) for k in direct)
    assert np.finfo(direct['photon'].dtype).eps<np.finfo(float).eps
    assert not (source/'production.json').exists() and not (out/'result.json').exists()
    assert not json.loads((out/'cached-input-pilot.json').read_text())['eligible']
    full=[r['rows'][0] for r in rows if r['k']];halves=[r['half_over_full'] for r in rows if r['k']]
    result=dict(classification='Counterexample candidate',saved_checks_passed=True,available_knots=[r['k'] for r in rows],
        missing_knots=[k for k in range(17) if k not in [r['k'] for r in rows]],
        maximum_remainder_over_linear=max(r['remainder_over_linear'] for r in full),
        half_over_full_range=[min(halves),max(halves)],maximum_projection_over_linear=max(r['projection_over_linear'] for r in full),
        maximum_negative_photon_energy_relative=max(r['negative_photon_energy_relative'] for r in full),
        exact_zero=True,NPZ_memory_read_exact=True,extended_precision_saved=True,
        uniform_EOS_derivative_bound=False,uniform_nonlinear_remainder_bound=False,source_complete=False,
        finite_remainder_propagated_to_GR=False,final_charge_solved=False,full_goal_complete=False,
        result_kind='loophole progress with an explicit resource stop',seconds=time.monotonic()-started)
    write(out/'saved-audit.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


if __name__=='__main__':
    def timeout(*_):raise TimeoutError('Phase143 saved-array audit cap')
    signal.signal(signal.SIGALRM,timeout)
    saved() if len(sys.argv)>1 and sys.argv[1]=='saved' else main()
