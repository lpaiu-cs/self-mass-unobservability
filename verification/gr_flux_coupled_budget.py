"""Measure actual endpoint energy budgets without subtracting large mass prefixes."""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from types import FunctionType
import argparse
import json

import numpy as np
import gr_coupled_evolution as e
import gr_conservative_evolution as restriction
import gr_flux_coupled_evolution as flux

OUT = flux.OUT


def plan4():
    flux.verify_initial()
    path = OUT/'time-refinement-plan.json'
    assert not path.exists()
    e.write(path, dict(classification='Counterexample candidate', steps=[4, 8, 16],
        minimum_thermal_heat_order=1.5,
        gate='All four maximum endpoint differences decrease; lnT and Q observed order at least 1.5. Same gate as the earlier coupled-evolution time comparison. No continuous error certificate.',
        timing='Additional four-step path specified while the original planned 8/16 corrected paths are still running; endpoint differences have not been evaluated.',
        bindings={str(p.relative_to(e.ROOT)):e.digest(p) for p in [Path(__file__), OUT/'initial-manifest.json']}))


def verify_four():
    plan = flux.verify_initial()
    for rel, h in json.loads((OUT/'time-refinement-plan.json').read_text())['bindings'].items():
        assert e.digest(e.ROOT/rel) == h, rel
    return dict(plan, steps=[4])


def refine4(workers):
    runner = FunctionType(flux.run.__code__, dict(flux.run.__globals__, verify_initial=verify_four))
    runner(4, workers)


def compare3():
    verify_four()
    plan = json.loads((OUT/'time-refinement-plan.json').read_text())
    path = OUT/'time-refinement.json'
    assert not path.exists()
    arrays, files = [], [OUT/'time-refinement-plan.json']
    for steps in plan['steps']:
        folder = OUT/f'path-{steps}'
        assert json.loads((folder/'result.json').read_text())['completed']
        manifest = folder/'manifest.json'
        for rel, h in json.loads(manifest.read_text())['sha256'].items(): assert e.digest(folder/rel) == h, rel
        arrays.append(np.load(folder/f'step-{steps:04d}.npz')['delta'][:, :4])
        files.append(manifest)
    for k in ['m', 'mf', 'a', 'N', 'Q', 'totals']:
        baseline = np.load(OUT/'path-4/step-0000.npz')[k]
        assert all(np.array_equal(np.load(OUT/f'path-{n}/step-0000.npz')[k], baseline) for n in [8, 16])
    errors = np.array([np.max(abs(b-a), axis=0) for a,b in zip(arrays[:-1], arrays[1:])])
    orders = np.log2(errors[0]/errors[1])
    passed = bool(np.all(errors[1] < errors[0]) and np.min(orders[[1, 3]]) >= plan['minimum_thermal_heat_order'])
    e.write(path, dict(classification='Counterexample candidate', passed=passed,
        endpoint_maximum_differences=errors.astype(float).tolist(), observed_orders=orders.astype(float).tolist(),
        rigorous_time_error_bound=False))
    files.append(path)
    e.write(OUT/'time-refinement-manifest.json', dict(sha256={str(p.relative_to(e.ROOT)):e.digest(p) for p in files}))
    print('CORRECTED COUPLED TIME REFINEMENT', passed, orders, flush=True)
    assert passed, 'Preserve the failed refinement; do not relax its criterion.'


def calculate(steps, aux):
    initial = np.load(OUT/'initial.npz')
    final = np.load(OUT/f'path-{steps}/step-{steps:04d}.npz')
    base, delta = initial['base'], final['delta']
    aux0, volume = initial['aux'], initial['volume']
    rho0 = np.exp(base[:, 0])
    rest0 = (base[:, 4:]/e.g.c.A)@e.g.c.W*e.C**2
    drest = (delta[:, 4:]/e.g.c.A)@e.g.c.W*e.C**2
    eps0 = rho0*(rest0+aux0[:, 2])
    deps = rho0*((rest0+aux0[:, 2])*np.expm1(delta[:, 0])
                +np.exp(delta[:, 0])*(drest+(aux[:, 2]-aux0[:, 2])))
    v = base[:, 2]+delta[:, 2]
    Q = (base[:, 3]+delta[:, 3])*initial['qscale']
    dE = (deps+(aux[:, 1]+eps0)*v*v+2*Q*v)/(1-v*v)
    heat = rho0*aux0[:, 10]*volume
    exchange = np.diff(final['integrated_face_mass_flux'])/e.GRAV
    defect = dE*volume+exchange
    normalized = defect/heat
    mass_defect = np.r_[e.ld(0), np.cumsum(e.GRAV*defect)]
    peak = int(np.argmax(abs(normalized)))
    report = dict(classification='Counterexample candidate', steps=steps,
        maximum_cell_energy_defect_over_initial_heat_capacity=float(abs(normalized[peak])),
        worst_original_cell=int(initial['indices'][peak]),
        total_energy_defect_over_initial_heat_capacity=float(defect.sum()/heat.sum()),
        maximum_face_mass_defect_over_initial_total=float(np.max(abs(mass_defect))/
                    (e.GRAV*np.sum(eps0*volume))),
        maximum_composition_change=float(np.max(abs(delta[:, 4:]))),
        maximum_metric_changes={k:float(np.max(abs(final[k]-np.load(OUT/f'path-{steps}/step-0000.npz')[k])))
                                for k in ['m', 'a', 'N']},
        raw_mass_prefix_diagnostic_superseded=True, native_or_continuum_error_certified=False,
        exact_discrete_energy_conservation=False)
    return report, dict(delta_energy_density=dE, cell_energy_defect=defect,
                        normalized_cell_energy_defect=normalized, face_mass_defect=mass_defect)


def budget(steps, workers):
    flux.verify_initial()
    folder = OUT/f'path-{steps}'
    manifest = json.loads((folder/'manifest.json').read_text())
    assert json.loads((folder/'result.json').read_text())['completed']
    for rel, h in manifest['sha256'].items(): assert e.digest(folder/rel) == h, rel
    assert not (OUT/f'budget-{steps}.json').exists()
    initial = np.load(OUT/'initial.npz')
    y = initial['base']+np.load(folder/f'step-{steps:04d}.npz')['delta']
    with ProcessPoolExecutor(max_workers=workers, initializer=e.worker_init) as pool:
        aux = np.asarray(list(pool.map(e.material, zip(y[:, 0], y[:, 1], y[:, 4:]), chunksize=4)), dtype=e.ld)
    assert np.all(np.isfinite(aux)) and np.all(aux[:, 10] > 0)
    report, data = calculate(steps, aux)
    np.savez_compressed(OUT/f'budget-{steps}.npz', aux=aux, **data)
    e.write(OUT/f'budget-{steps}.json', report)
    files = [Path(__file__), OUT/'initial-manifest.json', folder/'manifest.json',
             folder/f'step-{steps:04d}.npz', OUT/f'budget-{steps}.npz', OUT/f'budget-{steps}.json']
    e.write(OUT/f'budget-{steps}-manifest.json', dict(sha256={p.relative_to(e.ROOT).as_posix():e.digest(p) for p in files}))
    print('STABLE COUPLED ENERGY BUDGET', json.dumps(report), flush=True)


def verify():
    flux.verify_initial()
    assert flux.symbolic()['passed']
    for steps in [8, 16]:
        for rel, h in json.loads((OUT/f'budget-{steps}-manifest.json').read_text())['sha256'].items():
            assert e.digest(e.ROOT/rel) == h, rel
        data = np.load(OUT/f'budget-{steps}.npz')
        report, calculated = calculate(steps, data['aux'])
        assert report == json.loads((OUT/f'budget-{steps}.json').read_text())
        for k, v in calculated.items(): assert np.array_equal(v, data[k]), k
        folder = OUT/f'path-{steps}'
        for rel, h in json.loads((folder/'manifest.json').read_text())['sha256'].items():
            assert e.digest(folder/rel) == h, rel
    identity = json.loads((OUT/'initial-energy-identity.json').read_text())
    assert identity['passed']
    for rel, h in identity['bindings'].items(): assert e.digest(e.ROOT/rel) == h, rel
    stopped = json.loads((restriction.OUT/'heat-divergence-stop.json').read_text())
    assert stopped['stopped'] and all(not row['completed'] for row in stopped['paths'])
    for rel, h in stopped['sha256'].items(): assert e.digest(restriction.OUT/rel) == h, rel
    comparison = json.loads((OUT/'comparison.json').read_text())
    assert comparison['completed'] and comparison['same_initial_state']
    for rel, h in json.loads((OUT/'comparison-manifest.json').read_text())['sha256'].items():
        assert e.digest(e.ROOT/rel) == h, rel
    for rel, h in json.loads((OUT/'time-refinement-manifest.json').read_text())['sha256'].items():
        assert e.digest(e.ROOT/rel) == h, rel
    assert json.loads((OUT/'time-refinement.json').read_text())['passed']
    print('PASS conserved initial moments, actual corrected paths and reproducible stable energy budgets', flush=True)


def report():
    from verify_direct_eos_gr import append_scoped_material_progress
    verify()
    r = json.loads((OUT/'restriction.json').read_text())
    c = json.loads((OUT/'comparison.json').read_text())
    b = [json.loads((OUT/f'budget-{n}.json').read_text()) for n in [8, 16]]
    refinement = json.loads((OUT/'time-refinement.json').read_text())
    last = c['runs'][-1]['last']
    body = f'''분류: Counterexample candidate. 이번 작업은 실제 적분기의 초기 보존량과 열유속 차분 불일치를 수정한 loophole progress다. 기존 16절점 자료의 각 구역 바리온·좌표 에너지 적분을 목표로 고정했다. 에너지에서 공통 질량 면과 한 절점 계량을 계산하고 rho=B/(aV)를 정한 뒤 같은 native EOS의 내부에너지를 역산했다. 모든 5735구역이 사전 문턱을 통과했다. 원 총 바리온·질량에 대한 절대 상대 편차는 각각 {r['maximum_original_baryon_relative_defect']:.10g}, {r['maximum_original_mass_relative_defect']:.10g}이다. 이전 한 절점 초기 질량 편차 2.483851531e-7을 줄였으며, 목표 총질량이나 구적 가중치를 재조정하지 않았다. 구역별 최대 바리온/에너지 상대 복원 잔차는 {r['maximum_cell_baryon_relative_defect']:.8g}/{r['maximum_cell_energy_relative_defect']:.8g}이다. 이는 보존적 초기 재표현이지 고차 공간 진화나 엔트로피 보존 사상은 아니다.

분류: Proven. Q_r=(a/N) div_r(NQ/a)-Q(nu_r-a_r/a+2/r)인 연속 곱 미분 항등식을 검산했다. 방사 Einstein 제약과 초기 v=0을 적용하면 물질 내부에너지 식과 바리온 압축으로부터 E_t=-c div_r(NQ/a)를 얻는다. 수정 적분기는 Q의 공간 미분에 실제 공통 면 유속 차분을 사용한다. 유한 속도·조성 이류·다른 공간 연쇄 법칙의 이산 보존을 모두 증명한 것은 아니다.

분류: Counterexample candidate. 이전 적분기는 내부에너지 Q 미분에 절점 차분, 총에너지에는 면 유속 차분을 써서 표면에서 큰 불일치를 만들었다. 동일 초기 자료의 열용량 기준 에너지 변화율 잔차 최대가 261.3455949/s에서 9.139682031e-15/s로 줄었으며 초기 정지 상태에 대한 수치 검사다. 오류가 확인된 초기 재표현 경로는 8단계 계획의 6단계,16단계 계획의 8단계에서 모든 결과를 보존하고 중단했다. 이전 4/8/16 시간 수렴 결과는 원 공간 연산자에 대한 결과로 남지만, 큰 표면 온도 변화의 물리적 해석을 뒷받침하지 않는다.

분류: Counterexample candidate. 수정된 실제 EOS·열·유체·조성 이류·질량·계량 경로를 동일 0.000411300008893초까지 4·8·16단계로 완료했다. 추가4단계는 원8/16 경로 실행 도중 종료점 차이를 계산하기 전에 계획했으며, 예전과 같은 유한 수렴 문턱을 썼다. 초기 보존 자료와 초기 계량은 세 경로에서 동일하다.16단계 최대 변화 [lnrho,lnT,v/c,Q/초기 엔탈피]는 {last['maximum_changes']}이다.8/16 종료점 최대 차이는 {c['endpoint_maximum_8_16_differences']}이고 4/8/16 관측 차수는 {refinement['observed_orders']}이다. 네 차이의 감소와 온도/열유속 차수1.5 이상 기준을 통과했지만 엄밀한 시간 오차 상계는 아니다.

분류: Counterexample candidate. 표면 셀의 작은 열에너지는 큰 누적 질량을 뺀 값으로 분해할 수 없으므로, 별도 EOS 종료점 평가와 expm1 밀도·조성 정지에너지 증분으로 Delta E를 계산했다. 같은 면의 누적 에너지 교환을 더한 8/16단계 구역 에너지 결함/초기 열용량 최댓값은 각각 {b[0]['maximum_cell_energy_defect_over_initial_heat_capacity']:.10g}/{b[1]['maximum_cell_energy_defect_over_initial_heat_capacity']:.10g}이다.16단계 전체 열용량 기준 에너지 결함은 {b[1]['total_energy_defect_over_initial_heat_capacity']:.10g}이다. 원 누적 질량 차분 진단은 그대로 보존하되 이 국소 증분 결과로 해석을 정정한다. 정확한 이산 에너지 보존이나 native EOS 오차 인증으로 확대하지 않는다.

분류: Conjectural. 남은 실제 결합 병목은 유한 속도·조성 수송까지 일관된 국소 에너지 보존과 열 완화보다 긴 시간 구간을 감당하는 적분이다. 물리 대기/외부 경계, 핵반응의 장시간 결합, 물리 EOS·연속 미분 오차 및 실제 구동·관측 폐쇄는 남는다. 계산 경계와 짧은 지정 시간을 실제 항성 진화의 완성으로 세지 않는다.

재현: verification/gr_flux_coupled_budget.py verify. 원 경로 중단 상태·수정된 두 전체 시간 경로·보존 초기 자료·원시 종료점 EOS 및 면 에너지 교환을 함께 고정했다.
'''
    brief = ('분류: Counterexample candidate. 실제 적분기의 초기 질량 편차를 약1.04e-12로 줄이고 공통 면 열유속과 내부에너지 미분을 일치시켰다. 수정된5735구역 4/8/16 시간 경로와 원 유한 수렴 기준을 통과했다. 이전 표면 온도 변화의 물리적 해석은 지지되지 않으며 원 실패/시간 경로를 보존했다. 분류: Conjectural. 유한 속도·조성의 국소 에너지 보존, 긴 시간 적분·물리 경계·EOS·관측 폐쇄는 남는다.')
    append_scoped_material_progress('flux-coupled-evolution', '초기 질량 보존과 실제 열유속 결합 진화 수정', body, brief,
        '보존 초기 재표현과 면 열유속 차분 수정.5735구역 동일 기간 4/8/16 실제 진화·유한 수렴; 국소 보존·장시간·물리 경계·EOS·관측 미완료.', verify)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['budget', 'verify', 'report', 'plan4', 'refine4', 'compare3'])
    parser.add_argument('--steps', type=int, choices=[8, 16], default=16)
    parser.add_argument('--workers', type=int, default=8)
    args = parser.parse_args()
    if args.command == 'budget': budget(args.steps, args.workers)
    elif args.command == 'refine4': refine4(args.workers)
    else: globals()[args.command]()
