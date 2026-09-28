"""Archive the completed corrected pair and its own accepted charge/mass test."""
from pathlib import Path
from types import FunctionType
import datetime, importlib.util, sys

spec=importlib.util.spec_from_file_location('b',Path('.phase259-publish.py'))
b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master
read,write,sha,copy=b.read,b.write,b.sha,b.copy
out=root/'outputs/direct-eos-gr33/native-corrected-charge'
manifest=out.parent/'native-corrected-charge-manifest.json'
note=root/'notes/REQUEST260_CORRECTED_SAME_SOLUTION_CHARGE_KO.md'


def package():
    actual=runtime/'native-short-return259-work'
    charge=runtime/'native-short-return-charge259-work/full'
    exterior=runtime/'native-short-return-exterior259-work'
    status=read(actual/'pipeline-status.json')
    physical=read(actual/'result.json');compact=read(charge/'charge/result.json')
    ext=read(exterior/'full/result.json');audit=read(exterior/'full/audit.json')
    assert status['state']=='completed' and all(x['returncode']==0 for x in status['completed'])
    assert physical['passed'] and compact['charge_comparison_admitted'] and ext['passed'] and audit['passed']
    assert compact['same_actual_returned_solution_read'] and ext['same_actual_returned_solution_read']
    assert status['compact_sign_survives'] and status['frozen_exterior_sign_survives']
    assert [v['actual_completed_steps'] for v in physical['rows']]==[119,231]
    assert not physical['scientific_gates_changed']
    old=read(runtime/'native-material-exterior257-work/full/audit.json')
    mass_time=audit['controls']['low']['time']['homogeneous_mass_cm']
    old_time=old['controls']['low']['time']['homogeneous_mass_cm']
    assert not old['passed'] and old_time>.02 and mass_time<.02
    plans=[actual/'plan.json',actual/'controller-start.json',charge/'endpoint/plan.json',
           charge/'dense/plan.json',charge/'charge/plan.json',exterior/'full/plan.json',exterior/'full/audit.json']
    bindings={}
    for plan in plans:
        for p,h in read(plan)['bindings'].items():
            assert sha(runtime/p)==h,(str(plan),p)
            assert p not in bindings or bindings[p]==h,p
            bindings[p]=h
    dense={str(n):read(charge/f'dense/source-{n}-check.json') for n in [64,128]}
    assert all(x['passed'] and x['dense_stage_max']<1e-12 for x in dense.values())
    assert read(charge/'charge/polynomial-audit.json')['passed']
    assert read(exterior/'full/symbolic.json')['passed']
    for folder in [actual,charge,exterior]:
        for p in folder.rglob('*receipt.json'):assert read(p).get('error') is None,p
    preserved={p:h for p,h in read(b.manifest)['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    assert not out.exists();out.mkdir()
    for folder,label in [(actual,'actual259'),(charge,'charge259'),(exterior,'exterior259')]:
        for p in folder.rglob('*'):
            if p.is_file() and p.suffix in ['.json','.log','.py'] and not any(k in p.parts for k in ['sweep-0','initialization','__pycache__']):
                copy(p,out/label/p.relative_to(folder))
    for n in [64,128]:
        for name in [f'recovered-{n}.npz',f'sweep-1/photons/return-{n}.npz']:
            copy(actual/name,out/'actual259'/name)
        for name in [f'endpoint/gr/endpoint-{n}.npz',f'dense/gr/source-{n}.npz']:
            copy(charge/name,out/'charge259'/name)
    copy(charge/'dense/geometry.npz',out/'charge259/dense/geometry.npz')
    for n,q in [(64,8),(128,4),(128,8)]:
        copy(charge/f'charge/gr/fields-{n}-g{q}.npz',out/f'charge259/charge/gr/fields-{n}-g{q}.npz')
    for p in (exterior/'full').glob('*.npz'):copy(p,out/'exterior259/full'/p.name)
    copy(runtime/'.phase259-followthrough.py',out/'followthrough.py')
    copy(Path(__file__),out/'publication-producer.py')
    c=next(v for v in compact['components'] if v['clock']==128)['values']['endpoint_compact_with_metric']
    e=next(v for v in ext['components'] if v['clock']==128)
    now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    final=dict(classification='Counterexample candidate',passed=True,
        verdict='CORRECTED_SAME_SOLUTION_PASSES_MASS_TIME_AND_RETAINS_CONDITIONAL_NEGATIVE_CHARGE',
        actual_same_joint_solution_GR_computed=True,actual_steps=[119,231],
        full_declared_period=True,same_horizon_seconds=physical['same_horizon_seconds'],
        original_mass_time_relative=old_time,corrected_mass_time_relative=mass_time,mass_time_gate=.02,
        corrected_charge_time_relative=audit['controls']['low']['time']['normalized_standalone'],
        compact_sign_survives=True,frozen_exterior_sign_survives=True,
        conditional_fine_compact=c,conditional_fine_exterior=e,
        original_failures_preserved=True,scientific_gates_changed=False,
        physical_final_charge_solved=False,final_charge_conclusion='unadjudicated',full_goal_complete=False,
        snapshot_KST=now)
    note.write_text(f'''# 修正된 동일 결합 해의 질량·전하 판정

분류: Counterexample candidate. **압력일 구간 끝 미분을 수정한 같은 결합 해에서, 기존의 조건부 음의 전하 결론이 유지됐다.** 원 119/231단계가 모두 끝났으며 전체 기간은 {physical['same_horizon_seconds']:.16g}초다. 원 실패의 질량 시간 대조 {100*old_time:.8f}%는 {100*mass_time:.8f}%로 줄어, 완화하지 않은 2% 기준을 통과했다. 조건부 compact 판독과 고정 외부·질량 독립 감사가 모두 통과했다. 이번은 사용자가 정한 동일 해의 최종 전하 판정에 도달한 loophole progress다. 전체 물리 전하 폐쇄까지 완료했다는 뜻은 아니다. 게시={now}.

분류: Counterexample candidate. 단계258의 왼쪽 구간 끝 미분을 실제 적분에 적용했고, 단계259의 짧은 내부 반복과 참 잔차 갱신으로 비용을 줄였다. 저장 coarse 8단계는 정확한 재시작 검산 후 재사용했다. 이후 coarse/fine의 기존 시계, 공간·주파수·각도 격자, 구동, 기간, 선형·물질·비선형·보존·구성·전하 수락 기준을 유지했다. 단계257의 원 실패와 단계258의 의도적 작업 교체 기록을 변경하지 않았다. 기존 방정식에서 예측한 약 1/3의 질량 시간 오차 감소가 실제 판독에서도 나타났지만, 남은 1.32954%를 영으로 선언하지 않는다.

분류: Counterexample candidate. 두 경로의 참 단계 잔차 최댓값은 coarse={physical['rows'][0]['maximum_true_stage']:.12e}, fine={physical['rows'][1]['maximum_true_stage']:.12e}다. 에너지 보존 차이는 coarse={physical['rows'][0]['energy_balance_relative']:.12e}, fine={physical['rows'][1]['energy_balance_relative']:.12e}다. 전체 시간 상태 대조 최댓값={max(physical['time_relative']):.12e}; 같은 해의 dense 단계 잔차는 coarse={dense['64']['dense_stage_max']:.12e}, fine={dense['128']['dense_stage_max']:.12e}로 원 1e−12 기준을 통과했다.

분류: Counterexample candidate. 미세 시계의 compact high={float(c['high']):.12e}, low={float(c['low']):.12e}, low/high={c['low_over_high']:.12e}다. 동일 해의 고정 외부·질량 정규화 후 high={float(e['high']):.12e}, low 증분={float(e['same_solution_low_increment']):.12e}, 합={float(e['total']):.12e}이며 음의 부호가 유지됐다. 작은 low는 배경에 더해 소실시키지 않고 별도 성분 및 100자리 유리식으로 판독했다. low의 정규화 전하 시간 대조는 {100*final['corrected_charge_time_relative']:.8f}%다. 단일 출력의 부호만 보고 질량 실패를 무시하지 않았고, 이번에는 질량 기준도 실제 통과했다.

분류: Proven. 저장된 Radau 방출의 1차·2차 primitive와 질량 정규화 유리식에 대한 기호 항등식이 통과했다. 증명 범위는 해당 대수식이며, 동적 외부 전파나 균일 EOS·시간 오차의 증명이 아니다.

분류: Counterexample candidate. coarse 속행={read(actual/'coarse-receipt.json')['seconds']:.2f}초, fine={read(actual/'fine-receipt.json')['seconds']:.2f}초, 전체 자동 체인={status['elapsed_seconds']:.2f}초였다. 후속 판독은 저장된 같은 해를 사용했고 새 물리 단계가 없다. 제어기와 생산자 해시, 계획·실제 입력 해시, 판독·감사 바인딩을 다시 확인했다. 실제 계산 산출물은 WSL의 longdouble 상태를 그대로 복사했으며 Windows에서 재해석하지 않았다.

분류: Conjectural. 다음 병목은 실제 좌표 에너지 변환, 시간 의존 외부 광자의 에너지·방향·도착 변화, 그와 일치하는 질량·스칼라 경계 및 배경 질량 정규화를 같은 해에 연결하는 일이다. 기존 고정 외부 판독에 과거 다른 해의 진단값을 가산하지 않는다. 현재 해는 한 번의 GR 되먹임이며 자기 GR 수렴, 연속 EOS/균일 미분·공간 오차, 완전 비선형성, 정적 EFT와 관측 연결은 열린다. 이번 수락을 전하의 전역적·관측적 확정으로 확대하지 않는다.
'''.replace('修正','수정'),encoding='utf-8',newline='\n')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        p=root/f'docs/{name}.md';prefixes[p.relative_to(root).as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
        text=f'\n\n## 단계260 — 수정된 동일 해의 조건부 전하 수락\n\n분류: Counterexample candidate. 원119/231단계를 완료한 수정 결합 해의 질량 시간 차이는 {100*old_time:.8f}%에서 {100*mass_time:.8f}%로 줄어 원2%기준을 통과했다. 같은 해의 compact 및 고정 외부·질량 판독에서 기존 음의 전하 부호가 유지됐다. 원 실패·수락 기준은 보존했다. 전체 물리 전하의 미판정 사유는 이제 이 수치 수락 실패가 아니라 시간 의존 물리 외부·배경 정규화·자기GR 및 원 EOS/관측 폐쇄다. 이번은 동일 해의 전하까지 도달한 loophole progress다. [근거](../notes/{note.name}).\n'
        with p.open('ab') as f:f.write(text.encode())
    write(out/'result.json',final)
    write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes,verified_runtime_bindings=bindings))
    files=[note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes,sha256={p.relative_to(root).as_posix():sha(p) for p in files})
    write(manifest,final);m=read(master);m['sha256'].update(final['sha256'])
    m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_corrected_charge']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
