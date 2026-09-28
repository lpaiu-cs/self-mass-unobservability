"""Actual coupled evolution with matching heat flux and baryon chain rules."""
from pathlib import Path
from types import FunctionType, SimpleNamespace
import argparse
import json
import subprocess

import numpy as np
import sympy as sp
import gr_coupled_evolution as e
import gr_conservative_evolution as initial_driver
import gr_flux_coupled_evolution as heat
import gr_flux_coupled_budget as analysis

OUT = e.g.OUT/'gr-baryon-flux-evolution'


class BaryonStar(heat.FluxStar):
    def state(self, y):
        self.primitives = y
        return super().state(y)

    def gradient(self, value, odd=False):
        z = getattr(self, 'current', None)
        if z is not None and value is z['D']:
            # Use exactly the primitive gradients subtracted on converting
            # material rates back to coordinate rates. D has no X dependence.
            return (z['D']*super().gradient(self.primitives[:, 0])
                    +z['rho']*z['W']**3*z['v']*super().gradient(self.primitives[:, 2], odd=True))
        return super().gradient(value, odd=odd)


def symbolic():
    eta, v = sp.symbols('eta v', real=True)
    D = sp.exp(eta)/sp.sqrt(1-v*v)
    assert sp.simplify(sp.diff(D, eta)-D) == 0
    assert sp.simplify(sp.diff(D, v)-sp.exp(eta)*v/(1-v*v)**sp.Rational(3, 2)) == 0
    assert heat.symbolic()['passed']
    return dict(classification='Proven', passed=True,
        scope='Exact baryon primitive chain and the previous heat product identity. Matching gradients cancel the material/coordinate conversion terms algebraically; finite time, evolving metric and species-energy conservation remain numerical questions.')


def prepare():
    previous = heat.verify_initial()
    assert not OUT.exists()
    OUT.mkdir()
    files = [Path(__file__), Path(heat.__file__), Path(analysis.__file__),
             heat.OUT/'initial-manifest.json', heat.OUT/'budget-8-manifest.json', heat.OUT/'baryon-chain-stop.json']
    plan = dict(previous,
        checkpoint=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        steps=[4, 8, 16],
        bindings=dict(previous['bindings'], **{p.relative_to(e.ROOT).as_posix():e.digest(p) for p in files}),
        evolution='Same conservative reference moments, EOS, proper tau, computational boundary, heat-flux derivative and SSPRK2. For the baryon material/coordinate conversion use D_r=D*lnrho_r+rho*W^3*v*v_r with exactly the primitive gradients later subtracted. Do not separately difference D and then assume the discrete chain rule. Momentum and passive composition discretizations remain unchanged.',
        symbolic=symbolic())
    for name in ['initial.npz', 'restriction.json']:
        (OUT/name).write_bytes((heat.OUT/name).read_bytes())
    e.write(OUT/'plan.json', plan)
    e.write(OUT/'time-refinement-plan.json', dict(classification='Counterexample candidate',
        steps=[4, 8, 16], minimum_thermal_heat_order=1.5,
        gate='All four endpoint differences decrease; lnT/Q observed order at least 1.5. Same finite gate as earlier paths. Fixed before these corrected paths run.',
        rigorous_time_error_bound=False))
    e.write(OUT/'initial-manifest.json', dict(sha256={p.name:e.digest(p) for p in
        [OUT/'initial.npz', OUT/'restriction.json', OUT/'plan.json', OUT/'time-refinement-plan.json']}))
    verify_initial()


engine = SimpleNamespace(**vars(e))
engine.Star = BaryonStar
namespace = dict(vars(initial_driver), OUT=OUT, evolution=engine)
verify_initial = FunctionType(initial_driver.verify_initial.__code__, namespace)
namespace['verify_initial'] = verify_initial
run = FunctionType(initial_driver.run.__code__, namespace)
budget_namespace = dict(vars(analysis), OUT=OUT, flux=SimpleNamespace(verify_initial=verify_initial), verify_four=verify_initial)
calculate = FunctionType(analysis.calculate.__code__, budget_namespace)
budget_namespace['calculate'] = calculate
budget = FunctionType(analysis.budget.__code__, budget_namespace)
compare = FunctionType(analysis.compare3.__code__, budget_namespace)


def verify():
    verify_initial()
    assert symbolic()['passed']
    for steps in [4, 8, 16]:
        folder = OUT/f'path-{steps}'
        manifest = json.loads((folder/'manifest.json').read_text())
        assert manifest['initial_manifest_sha256'] == e.digest(OUT/'initial-manifest.json')
        assert json.loads((folder/'result.json').read_text())['completed']
        for rel, h in manifest['sha256'].items(): assert e.digest(folder/rel) == h, rel
    for steps in [8, 16]:
        for rel, h in json.loads((OUT/f'budget-{steps}-manifest.json').read_text())['sha256'].items():
            assert e.digest(e.ROOT/rel) == h, rel
        data = np.load(OUT/f'budget-{steps}.npz')
        result, arrays = calculate(steps, data['aux'])
        assert result == json.loads((OUT/f'budget-{steps}.json').read_text())
        for k, value in arrays.items(): assert np.array_equal(value, data[k]), k
    for rel, h in json.loads((OUT/'time-refinement-manifest.json').read_text())['sha256'].items():
        assert e.digest(e.ROOT/rel) == h, rel
    assert json.loads((OUT/'time-refinement.json').read_text())['passed']
    for folder, name in [(initial_driver.OUT, 'heat-divergence-stop.json'), (heat.OUT, 'baryon-chain-stop.json')]:
        stopped = json.loads((folder/name).read_text())
        assert stopped['stopped']
        for rel, h in stopped['sha256'].items(): assert e.digest(folder/rel) == h, rel
    print('PASS actual corrected coupled paths, finite time refinement and stable endpoint energy budgets', flush=True)


def report():
    from verify_direct_eos_gr import append_scoped_material_progress
    verify()
    r = json.loads((OUT/'restriction.json').read_text())
    t = json.loads((OUT/'time-refinement.json').read_text())
    b = [json.loads((OUT/f'budget-{n}.json').read_text()) for n in [8, 16]]
    last = json.loads((OUT/'path-16/result.json').read_text())['last']
    body = f'''분류: Counterexample candidate. 이번 작업은 실제 결합 적분기의 보존 오차 두 원인을 수정한 loophole progress다. 기존16절점 참조의 각 구역 바리온·좌표 에너지를 보존적으로 초기 상태에 옮겼다. 에너지에서 공통 질량/계량을 계산하고 rho=B/(aV)를 정한 뒤 native 내부에너지를 역산했다. 총질량이나 구적 가중치를 재조정하지 않았다.5735구역 모두 사전 문턱을 통과했고 원 바리온/질량에 대한 절대 상대 편차는 {r['maximum_original_baryon_relative_defect']:.10g}/{r['maximum_original_mass_relative_defect']:.10g}이다. 이전 초기 질량 편차2.483851531e-7을 줄였으나 고차 공간 진화나 엔트로피 보존 사상은 아니다.

분류: Proven. Q_r=(a/N)div_r(NQ/a)-Q(nu_r-a_r/a+2/r)와 D_r=D(lnrho)_r+rho W^3 v v_r를 검산했다. 열식에 실제 공통 면 유속 차분을 쓰고, 바리온의 물질/좌표 시간 변환에는 실제로 나중에 빼는 것과 같은 원시 변수 기울기를 쓴다. 연속 곱/연쇄 법칙을 서로 다른 차분에 자동 적용했던 두 불일치를 제거한다. 기호 항등식은 유한 시간·계량·조성의 완전한 이산 보존을 보증하지 않는다.

분류: Counterexample candidate. 열유속 수정은 같은 초기 정지 상태의 에너지 변화율 잔차/열용량을 최대261.3455949/s에서9.139682031e-15/s로 줄였다. 그러나 열 수정만 한8단계 종료점에는 구역 에너지 결함/초기 열용량 최대0.03146190002가 남았다. 같은 종료 상태의 표면에서는 조성 변화율이 정확히0이고, 바리온 차분까지 맞추면 에너지 변화율 잔차가 약168.4039/s에서9.2564e-6/s로 줄었다. 원 부분 경로와 열 수정만 한8단계 전체 경로/16단계 중단 상태를 보존했다. 그 경로용 추가4단계 계획은 실행하지 않았으며 최종 세 경로는 별도 소스·계획으로 다시 고정했다. 이전 큰 표면 온도 변화와 시간 수렴만으로 물리적 결과를 주장하지 않는다.

분류: Counterexample candidate. 두 수정이 들어간 실제 EOS·불투명도·열·유체·조성 이류·질량·계량 경로를 원5735구역과 같은0.000411300008893초까지4·8·16단계로 완료했다.16단계 최대 변화 [lnrho,lnT,v/c,Q/초기 엔탈피]는 {last['maximum_changes']}이다.4/8 및8/16 종료점 최대 차이는 {t['endpoint_maximum_differences']}이고 관측 차수는 {t['observed_orders']}이다. 초기 자료는 같으며 기존의 네 차이 감소/온도·열유속 차수1.5 이상 기준을 통과했다. 이는 이 공간 근사의 유한 시간 수렴 근거이며 연속 시간 오차 인증은 아니다.

분류: Counterexample candidate. 큰 누적 질량 차분이 작은 표면 열에너지를 분해하지 못하므로, 원시 EOS 종료점과 expm1 밀도·조성 정지에너지 증분을 사용해 Delta E와 각 공통 면의 누적 교환을 직접 비교했다.8/16단계 구역 에너지 결함/초기 열용량 최댓값은 {b[0]['maximum_cell_energy_defect_over_initial_heat_capacity']:.10g}/{b[1]['maximum_cell_energy_defect_over_initial_heat_capacity']:.10g}이고,16단계 전체 열용량 기준 결함은 {b[1]['total_energy_defect_over_initial_heat_capacity']:.10g}이다.16단계의 저장 질량·방사 계량·lapse 최대 변화는 {b[1]['maximum_metric_changes']}이다. 원 누적 질량 차분 진단을 보존하되 이 증분 계산으로 국소 수지를 판단한다. 정확한 이산 에너지 보존 또는 native EOS 오차 인증으로 확대하지 않는다.

분류: Conjectural. 조성/운동량까지 일관된 유한 공간 보존, 긴 시간 구간의 강직한 열·유체 적분, 물리 대기/외부 경계와 핵반응 결합은 남는다. 계산 경계·지정 완화시간의 짧은 경로를 실제 항성 matching·물리 EOS·연속 미분 인증 또는 관측 추론 완성으로 세지 않는다.

재현: verification/gr_baryon_flux_evolution.py verify. 보존 초기 자료, 원 오류 경로, 수정4/8/16 전체 경로와 원시 EOS 에너지 증분을 함께 고정했다.
'''
    brief = ('분류: Counterexample candidate. 초기 질량 편차를 약1.04e-12로 줄이고 열 면 유속과 바리온 원시 변수 차분의 두 불일치를 실제 적분기에서 수정했다. 수정5735구역 4/8/16 전체 경로와 원 유한 시간 수렴 기준을 통과했다. 원 오류 경로·국소 에너지 결함은 보존했다. 분류: Conjectural. 공간 보존·장시간·물리 경계·EOS·관측 폐쇄는 남는다.')
    append_scoped_material_progress('baryon-flux-evolution', '실제 결합 진화의 초기 보존과 두 공간 차분 수정', body, brief,
        '초기 질량 및 열·바리온 차분 수정.5735구역 4/8/16 실제 경로·유한 수렴; 공간/장시간/물리 EOS·경계·관측은 미완료.', verify)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['prepare', 'run', 'budget', 'compare', 'verify', 'report'])
    parser.add_argument('--steps', type=int, choices=[4, 8, 16], default=16)
    parser.add_argument('--workers', type=int, default=8)
    args = parser.parse_args()
    if args.command == 'run': run(args.steps, args.workers)
    elif args.command == 'budget': budget(args.steps, args.workers)
    else: globals()[args.command]()
