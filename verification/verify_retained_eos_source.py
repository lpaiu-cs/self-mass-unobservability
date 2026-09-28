"""Check the recorded native replacement and its same-time infinity mapping."""
from pathlib import Path
import hashlib,io,json,signal,time

ROOT=Path('outputs/direct-eos-gr33');OUT=ROOT/'retained-eos-source'
CURRENT=ROOT/'native-retained-completion';read=lambda p:json.loads(p.read_text())


def main():
    start=time.monotonic()
    receipt=read(OUT/'receipt.json');result=read(OUT/'result.json');assert receipt['returncode']==0 and result['passed']
    plan=read(OUT/'plan.json');assert not (OUT/'verification.json').exists()
    remaining=min(120-result['seconds'],1140-receipt['seconds']);cap=min(60.,remaining);assert cap>0
    def timeout(*_):raise TimeoutError('Remaining Phase149 verification allocation')
    signal.signal(signal.SIGALRM,timeout);signal.setitimer(signal.ITIMER_REAL,max(.001,cap-(time.monotonic()-start)))
    import numpy as np
    import sympy as sp
    load=lambda p:np.load(io.BytesIO(p.read_bytes()))
    rows=read(OUT/'production-samples.json');production=read(OUT/'production.json');pilot=read(OUT/'pilot.json')
    assert production['passed'] and read(OUT/'audit.json')['passed'] and pilot['eligible']
    assert len(rows)==production['states']==len({(r['kind'],r['it'],r['cell']) for r in rows})
    assert sum(r['reused'] for r in rows)==sum(production['exact_reused'].values())
    assert {r['native_owner'] for r in rows if not r['reused']}=={'deep','high','low'}
    for region in pilot['regions'].values():
        assert region['reseed'] and all(r['relative']<plan['gates']['reseed'] for r in region['reseed'])
    fine=load(OUT/'wave-fine.npz');applied=load(OUT/'applied-charge.npz');normalized=load(OUT/'normalized-charge.npz')
    original=load(CURRENT/'gr/wave-128-g8.npz');infinity=load(CURRENT/'infinity/completed/retained-128-a8-r8.npz')
    ids=np.array([int(np.argmin(abs(fine['t']-t))) for t in infinity['t']])
    assert np.max(abs(fine['t'][ids]-infinity['t']))<1e-18
    assert np.array_equal(ids,np.arange(0,129,8)) and np.array_equal(normalized['t'],infinity['t'])
    corrected_scalar=applied['free_scalar'][ids]+infinity['exterior']
    rebuilt=(corrected_scalar+infinity['mass_term'])/(1-infinity['epsilon'])
    norm=float(np.max(abs(rebuilt)));error=float(np.max(abs(rebuilt-normalized['normalized'])))/norm
    assert error<1e-12
    linear=float(np.max(abs(applied['free_scalar']-original['free_scalar']-fine['free_scalar'])))/float(np.max(abs(original['free_scalar'])))
    assert linear<plan['gates']['linear']
    A,s,e,d=sp.symbols('A s e d')
    assert sp.simplify((s+d+A*e)/(1-e)-(s+A*e)/(1-e)-d/(1-e))==0
    flags=['native_correction_time_convergence_tested','uniform_EOS_derivative_bound','full_floor_feedback_enclosed',
        'coupled_fixed_point_verified','nonlinear_GR','final_charge_solved','full_goal_complete']
    assert all(result[k] is False for k in flags)
    assert production['seconds']<900 and production['worker_seconds']<2400 and production['native_calls']<80000
    assert pilot['seconds']<120 and result['seconds']<120
    files=[Path(__file__),OUT/'plan.json',OUT/'receipt.json',OUT/'result.json',OUT/'production.json',
        OUT/'production-samples.json',OUT/'normalized-charge.npz',OUT/'applied-charge.npz']
    bindings={p.as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    seconds=time.monotonic()-start;assert seconds<cap
    evidence=dict(classification='Counterexample candidate',passed=True,states=len(rows),
        independent_time_mapping=True,independent_normalization_relative=error,linear_application_relative=linear,
        symbolic_normalization_checked=True,prior_photon_emission_held_fixed=True,physical_time_paths_replayed=0,
        full_goal_complete=False,seconds=seconds,cap_seconds=cap,
        action_seconds_with_verification=receipt['seconds']+seconds,
        bindings=bindings)
    (OUT/'verification.json').write_text(json.dumps(evidence,indent=2)+'\n')
    assert time.monotonic()-start<cap
    signal.setitimer(signal.ITIMER_REAL,0.);print(json.dumps(evidence),flush=True)


if __name__=='__main__':main()
