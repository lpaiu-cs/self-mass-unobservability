"""Postprocess completed native implicit paths; never alter a running attempt."""
import argparse
import json
from pathlib import Path
import numpy as np

import gr_implicit_coupled_evolution as evolution

e, OUT = evolution.e, evolution.OUT


def cones(z):
    """Floating check of the existing full heat-law principal matrices.

    The local rest-frame matrices are those derived in gr_heat_entropy_closure.
    This is a physical validity check on saved evolving states, not an interval
    certificate or a modification of the original time-refinement criterion.
    """
    aux, w, rho, P = z['aux'], z['w'], z['rho'], z['P']
    b, r = rho*aux[:, 10]/w, P*aux[:, 5]/w
    d, therm = P*aux[:, 6]/w, (P-rho*aux[:, 9])/w
    a = z['K']*z['T']/(e.C**2*w*e.TAU)
    j, kr, kt = z['Q']/w, -aux[:, 22], 5-aux[:, 23]
    M = np.zeros((len(w), 4, 4))
    N = M.copy()
    M[:, 0, 0] = 1
    M[:, 1, 1], M[:, 1, 2], M[:, 2, 2:] = b, 2*j, 1
    M[:, 3, 0], M[:, 3, 1], M[:, 3, 2], M[:, 3, 3] = -j*kr/2, -j*kt/2, a, 1
    N[:, 0, 2], N[:, 1, 2], N[:, 1, 3] = 1, therm, 1
    N[:, 2, 0], N[:, 2, 1], N[:, 2, 2], N[:, 3, 1] = r, d, 2*j, a
    roots = np.linalg.eigvals(np.linalg.solve(M, N))
    speed = np.max(abs(roots.real), axis=1)
    return dict(maximum_local_rest_characteristic_speed_over_c=float(speed.max()),
        maximum_characteristic_imaginary_part=float(np.max(abs(roots.imag))),
        sampled_cone_inside_light_cone=bool(np.all(roots.imag == 0) and np.all(speed < 1)))


def resolved(refinement):
    folder = OUT/f'path-{refinement}'
    link = folder/'continuation.json'
    if not link.exists():
        return folder, OUT/'plan.json'
    record = json.loads(link.read_text())
    continued = e.ROOT/record['path']
    assert continued.resolve().is_relative_to(e.g.OUT.resolve())
    assert e.digest(continued/'manifest.json') == record['manifest_sha256']
    resume = json.loads((continued/'resume.json').read_text())
    assert (e.ROOT/resume['parent']).resolve() == folder.resolve()
    for name, digest in resume['parent_files_sha256'].items():
        assert e.digest(folder/name) == digest, name
    for step in range(resume['start']+1):
        name = f'step-{step:04d}.npz'
        assert e.digest(continued/name) == e.digest(folder/name)
    replay = json.loads((continued/'checkpoint-replay.json').read_text())
    assert replay['passed'] and replay['fresh_native']
    path = continued.parent/'plan.json'
    current, previous = json.loads(path.read_text()), json.loads((OUT/'plan.json').read_text())
    assert current['parent_plan_sha256'] == e.digest(OUT/'plan.json')
    for key in ['cells', 'time_edges_tau', 'duration_tau', 'duration_seconds', 'refinements',
                'nonlinear_absolute_tolerances', 'nonlinear_relative_tolerance', 'time_refinement_gate', 'boundary']:
        assert current[key] == previous[key], key
    return continued, path


def completed(refinement):
    folder, plan_path = resolved(refinement)
    result = json.loads((folder/'result.json').read_text())
    assert result['completed'] and not (folder/'failure.json').exists()
    manifest = json.loads((folder/'manifest.json').read_text())
    assert manifest['plan_sha256'] == e.digest(plan_path)
    for rel, digest in manifest['sha256'].items():
        assert e.digest(folder/rel) == digest, rel
    plan = json.loads(plan_path.read_text())
    for rel, digest in plan['bindings'].items():
        assert e.digest(e.ROOT/rel) == digest, rel
    assert result['steps'] == refinement*(len(plan['time_edges_tau'])-1)
    logs = [json.loads(line) for line in (folder/'iterations.jsonl').read_text().splitlines()]
    stages = {}
    for row in logs:
        stages.setdefault((row['step'], row['stage']), []).append(row)
    assert set(stages) == {(i, j) for i in range(1, result['steps']+1) for j in [1, 2]}
    for rows in stages.values():
        assert [r['iteration'] for r in rows] == list(range(len(rows)))
        assert rows[-1]['residual_norm'] <= 1
        assert len(rows) <= plan['maximum_stage_iterations']
    final = np.load(folder/f"step-{result['steps']:04d}.npz")
    initial = np.load(folder/'step-0000.npz')
    star = evolution.initialize(None)
    assert np.all(initial['delta'] == 0)
    y = star.base+final['delta']
    star.material_cache = {e.material_key(row): aux for row, aux in zip(
        zip(y[:, 0], y[:, 1], y[:, 4:]), final['aux'])}
    z = star.state(y)
    for key in ['m', 'mf', 'a', 'N', 'Q']:
        assert np.array_equal(z[key], final[key]), key
    defect, budget = evolution.energy_budget(star, final['delta'], z, final['integrated_face_mass_flux'])
    assert np.array_equal(defect, final['normalized_cell_energy_defect'])
    cone = cones(z)
    maximum_path_energy_defect = float(np.max(abs(defect)))
    # Inspect every saved state, since an admissible endpoint cannot rescue an
    # earlier sampled excursion outside the declared heat model's light cone.
    for step in range(result['steps']):
        with np.load(folder/f'step-{step:04d}.npz') as saved:
            sample = star.base+saved['delta']
            star.material_cache = {e.material_key(row): aux for row, aux in zip(
                zip(sample[:, 0], sample[:, 1], sample[:, 4:]), saved['aux'])}
            sample_state = star.state(sample)
            for key in ['m', 'mf', 'a', 'N', 'Q']:
                assert np.array_equal(sample_state[key], saved[key]), (step, key)
            sampled = cones(sample_state)
            cone['sampled_cone_inside_light_cone'] &= sampled['sampled_cone_inside_light_cone']
            for key in ['maximum_local_rest_characteristic_speed_over_c', 'maximum_characteristic_imaginary_part']:
                cone[key] = max(cone[key], sampled[key])
            sample_defect, _ = evolution.energy_budget(star, saved['delta'], sample_state, saved['integrated_face_mass_flux'])
            assert np.array_equal(sample_defect, saved['normalized_cell_energy_defect'])
            maximum_path_energy_defect = max(maximum_path_energy_defect, float(np.max(abs(sample_defect))))
    start_b = np.sum(initial['a']*np.exp(star.base[:, 0])*star.volume)
    end_b = np.sum(z['a']*z['D']*star.volume)
    row = dict(classification='Counterexample candidate', refinement=refinement,
        steps=result['steps'], duration_seconds=float(final['time_seconds']),
        maximum_logged_accepted_stage_norm=max(rows[-1]['residual_norm'] for rows in stages.values()),
        native_stage_evaluations=sum(rows[-1].get('native_evaluations', len(rows)) for rows in stages.values()),
        maximum_changes=np.max(abs(final['delta'][:, :4]), axis=0).astype(float).tolist(),
        maximum_composition_change=float(np.max(abs(final['delta'][:, 4:]))),
        relative_total_baryon_budget_defect=float((end_b-start_b-final['boundary_exchange'][0])/start_b),
        maximum_metric_changes={key:float(np.max(abs(final[key]-initial[key]))) for key in ['m', 'a', 'N']},
        maximum_saved_cell_energy_defect_over_initial_heat_capacity=maximum_path_energy_defect,
        **budget, **cone)
    assert row['duration_seconds'] == plan['duration_seconds']
    return final['delta'].copy(), row


def compare():
    target = OUT/'time-refinement.json'
    assert not target.exists(), 'Keep every original refinement verdict.'
    paths = [completed(r) for r in [1, 2, 4]]
    errors = np.array([np.max(abs(right[0][:, :4]-left[0][:, :4]), axis=0)
                       for left, right in zip(paths[:-1], paths[1:])])
    assert np.all(errors > 0)
    orders = np.log2(errors[0]/errors[1])
    passed = bool(np.all(errors[1] < errors[0]) and np.min(orders[[1, 3]]) >= 1.5)
    result = dict(classification='Counterexample candidate', passed=passed,
        endpoint_maximum_differences=errors.astype(float).tolist(), observed_orders=orders.tolist(),
        paths=[p[1] for p in paths], native_EOS_error_certified=False,
        rigorous_time_error_bound=False, physical_exterior_match=False, observational_closure=False)
    e.write(target, result)
    files = [Path(__file__), target, OUT/'plan.json']
    for refinement in [1, 2, 4]:
        folder, plan_path = resolved(refinement)
        files.extend([folder/'manifest.json', plan_path])
        link = OUT/f'path-{refinement}/continuation.json'
        if link.exists():
            files.append(link)
    e.write(OUT/'time-refinement-manifest.json', dict(sha256={p.relative_to(e.ROOT).as_posix():e.digest(p) for p in files}))
    print('NATIVE IMPLICIT TIME REFINEMENT', json.dumps(result), flush=True)
    assert passed, 'The original finite refinement gate failed; preserve this result.'


def verify():
    for rel, digest in json.loads((OUT/'time-refinement-manifest.json').read_text())['sha256'].items():
        assert e.digest(e.ROOT/rel) == digest, rel
    result = json.loads((OUT/'time-refinement.json').read_text())
    assert result['passed']
    for refinement, row in zip([1, 2, 4], result['paths']):
        assert completed(refinement)[1] == row
    assert evolution.symbolic()['passed']
    print('PASS completed native implicit paths and declared finite refinement', flush=True)


def report():
    from verify_direct_eos_gr import append_scoped_material_progress
    verify()
    result = json.loads((OUT/'time-refinement.json').read_text())
    final = result['paths'][-1]
    assert all(row['sampled_cone_inside_light_cone'] for row in result['paths']), 'Preserve the computed trajectory, but resolve the physical cone failure before issuing this progress report.'
    body = f'''분류: Counterexample candidate. 동일5735구역·보존 초기 자료·물리식을 유지하고 강직한 열·유체 결합 적분에 SDIRK2를 적용했다. 근사 행렬에서만 계량과 EOS 미분을 고정하고 수락 잔차에는 실제 native EOS·불투명도·조성 이류·유체·열식과 재계산 질량/lapse를 사용했다. 수락된 단계의 기록된 최대 정규화 잔차는 {final['maximum_logged_accepted_stage_norm']:.9g}이다. 이는 선언한 비선형 대수 잔차 기준이며 외부 독립 EOS 인증은 아니다.

분류: Proven. 지정 SDIRK2 공식의2차 조건과 L 안정성을 기호 검산했다. 이 선형 안정성 결과는 비선형 수렴·공간 오차·물리 모델의 정확성을 보증하지 않는다.

분류: Counterexample candidate. 세 경로는 같은 좌표 시간 {final['duration_seconds']:.12g}초까지 {[p['steps'] for p in result['paths']]}단계로 완료했다. 초기 과도구간을 포함하는 같은 비균일 시간 격자를 순서대로 이분했다. 종료점 [lnrho,lnT,v/c,Q/초기 엔탈피] 차이의 관측 차수는 {result['observed_orders']}이고 원래의 네 차이 감소 및 온도·열유속 차수1.5 이상 기준을 통과했다. 가장 세밀한 경로의 최대 변화는 {final['maximum_changes']}이다. 최대 구역 에너지 결함/초기 열용량은 {final['maximum_cell_energy_defect_over_initial_heat_capacity']:.9g}, 전체 바리온 수지 상대 결함은 {final['relative_total_baryon_budget_defect']:.9g}이다. 같은 native 내부에너지 증분과 공통 면 교환으로 에너지 수지를 계산했고 질량/열량을 재조정하지 않았다.

분류: Conjectural. 유한 시간 경로의 진전이며 물리 외부 경계·핵반응 결합·공간 및 native EOS/연속 오차 인증·실제 구동/관측 폐쇄는 남는다. 지정 완화시간의 계산 모형을 실제 항성의 전체 진화나 투고 가능성 인증으로 확대하지 않는다.
'''
    brief = ('분류: Counterexample candidate. 실제 native 결합 잔차를 푸는 SDIRK2로 같은5735구역의0.421초 경로와 세 시간 해상도 대조를 완료했다. 분류: Conjectural. 물리 외부 경계·반응·EOS/연속 오차·관측 폐쇄는 남는다.')
    append_scoped_material_progress('implicit-coupled-evolution', '강직한 실제 결합 진화의 암시적 적분', body, brief,
        'Native 결합 암시적 경로·유한 시간 대조. 물리 경계·반응·EOS/연속 오차·관측은 미완료.', verify)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['completed', 'compare', 'verify', 'report'])
    parser.add_argument('--refinement', type=int, choices=[1, 2, 4], default=1)
    args = parser.parse_args()
    if args.command == 'completed':
        print(json.dumps(completed(args.refinement)[1], indent=2))
    else:
        globals()[args.command]()
