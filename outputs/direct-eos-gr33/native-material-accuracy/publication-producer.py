"""Publish an actual material-accuracy repair through its same-solution charge."""
from pathlib import Path
from types import FunctionType
import datetime,importlib.util,sys
spec=importlib.util.spec_from_file_location('b',Path('.phase256-publish.py'));b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha,copy=b.read,b.write,b.sha,b.copy
out=root/'outputs/direct-eos-gr33/native-material-accuracy'
manifest=out.parent/'native-material-accuracy-manifest.json'
note=root/'notes/REQUEST257_MATERIAL_ACCURACY_CHARGE_KO.md'


def package():
    actual=runtime/'native-material-accuracy257-work'
    charge=runtime/'native-material-charge257-work/full'
    exterior=runtime/'native-material-exterior257-work'
    for p in [actual/'controller-status.json',charge/'controller-status.json']:
        assert read(p)['state']=='completed',p
    assert read(exterior/'full-controller-status.json')['state']=='failed'
    physical=read(actual/'result.json');compact=read(charge/'charge/result.json')
    ext=read(exterior/'full/result.json');audit=read(exterior/'full/audit.json')
    assert physical['passed'] and compact['charge_comparison_admitted'] and ext['passed'] and not audit['passed']
    mass_time=audit['controls']['low']['time']['homogeneous_mass_cm']
    assert mass_time>.02 and audit['conditional_compact_sign_survives_frozen_exterior'] is None
    ledger=read(exterior/'mass-ledger.json');assert not ledger['original_audit_passed']
    assert read(actual/'restart-regression.json')['every_saved_array_exact']
    assert read(runtime/'phase257-symbolic-adapter-check.json')['passed']
    for p,h in read(actual/'plan.json')['bindings'].items():assert sha(runtime/p)==h,p
    for plan in [charge/'endpoint/plan.json',charge/'dense/plan.json',charge/'charge/plan.json',exterior/'full/plan.json']:
        for p,h in read(plan)['bindings'].items():assert sha(runtime/p)==h,(plan,p)
    dense={str(n):read(charge/f'dense/source-{n}-check.json') for n in [64,128]}
    assert all(v['passed'] and v['dense_stage_max']<1e-12 for v in dense.values())
    modules=[root/'verification'/n for n in ['finish_returned_material_accuracy.py','read_material_return_charge.py','read_material_return_exterior.py']]
    assert all(sha(p)==sha(runtime/'verification'/p.name) for p in modules)
    preserved={p:h for p,h in read(b.manifest)['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    assert not out.exists();out.mkdir()
    for folder,label in [(actual,'actual257'),(charge,'charge257'),(exterior,'exterior257')]:
        for p in folder.rglob('*'):
            if p.is_file() and p.suffix in ['.json','.log','.py'] and not any(k in p.parts for k in ['sweep-0','initialization','__pycache__']):
                copy(p,out/label/p.relative_to(folder))
    for n in [64,128]:
        for name in [f'recovered-{n}.npz',f'sweep-1/photons/return-{n}.npz']:
            copy(actual/name,out/'actual257'/name)
        for name in [f'endpoint/gr/endpoint-{n}.npz',f'dense/gr/source-{n}.npz']:
            copy(charge/name,out/'charge257'/name)
    copy(charge/'dense/geometry.npz',out/'charge257/dense/geometry.npz')
    for n,q in [(64,8),(128,4),(128,8)]:copy(charge/f'charge/gr/fields-{n}-g{q}.npz',out/f'charge257/charge/gr/fields-{n}-g{q}.npz')
    for p in (exterior/'full').glob('*.npz'):copy(p,out/'exterior257/full'/p.name)
    for p in exterior.glob('mass-ledger-*.npz'):copy(p,out/'exterior257'/p.name)
    failed=runtime/'native-right-charge256-work/full'
    for name in ['controller-status.json','controller-start.json','dense_source.stderr.log','dense/source-receipt.json','dense/source-64-check.json','dense/source-128-check.json']:
        copy(failed/name,out/'failed256-charge'/name)
    for name in ['phase257-neutral-precision.json','phase257-symbolic-adapter-check.json']:
        copy(runtime/name,out/name)
    for name in ['.phase257-followthrough.py','.phase257-charge-followthrough.py','.phase257-exterior-followthrough.py','.phase257-neutral-precision.py','.phase257-readout-check.py','.phase257-mass-ledger.py']:
        copy(root/name,out/name.lstrip('.'))
    copy(Path(__file__),out/'publication-producer.py')
    row=read(actual/'run-128.json');cost=read(actual/'fine-receipt.json')
    c=next(v for v in compact['components'] if v['clock']==128)['values']['endpoint_compact_with_metric']
    e=next(v for v in ext['components'] if v['clock']==128)
    old=read(failed/'dense/source-128-check.json')
    now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    note.write_text(f'''# 물질 성분 정확도 수정 후 동일 결합 해의 전하

분류: Counterexample candidate. **최종 물리 전하는 아직 미판정이다.** 전체 기간의 같은 119/231단계 해에서 물질 성분의 실제 단계 잔차를 수정하고, 원 compact 전하 판독까지 통과했다. 조건부 compact 부호 유지={compact['conditional_compact_sign_survives_one_return']}. 그러나 고정 초기 외부·질량 판독은 낮은 성분의 동질 질량 시간 대조가 {100*mass_time:.8f}%로 원 2% 기준을 넘어 독립 감사에서 탈락했다. 따라서 외부·질량을 포함한 부호는 수락하지 않았다. 이전 256의 전체 해는 dense 기준에 실패했으므로 수락된 전체 전하 기준값이 없었다. 게시={now}.

분류: Counterexample candidate. 256의 실제 전체 적분·시간 대조는 통과했지만, 같은 해의 fine 220단계에서 두 번째 내부 단계 H의 상대 잔차 {old['dense_stage_max']:.12e}가 원 dense 기준 1e−12를 넘었다. 저장 보존 변수에서 직접 계산한 80자리 잔차는 1.719880873310006e−12였고, 원 방식과의 차이는 2.19e−20 규모였다. 따라서 읽기 과정의 덧셈 오차를 고쳐 해결할 문제가 아니었다. 이 진단의 floor 처리 한계와 원 실패는 원본에 보존한다.

분류: Proven. Radau 연속식의 적분 가중치를 단계 1/3과 1에 평가하면 각각 (5/12,−1/12), (3/4,1/4)이며, 고정 대각 보존 단위변환은 단계 잔차식과 교환한다. 기호 검사는 이 항등식과 메타데이터 키 회귀를 확인하며, 균일 미분·수치 오차나 물리 결과를 증명하지 않는다.

분류: Counterexample candidate. 기존 오른쪽 전처리 풀이를 실제 결합 단계에 적용하고, 각 물질 Etilde/H/B/S 성분이 자기 크기에 대해 1e−13 미만이어야 수락하도록 내부 기준을 강화했다. 원 선형 1e−14/합성 물리 1e−13, 실제 비선형 1e−12, 구성 0.2%, 보존 1e−8, port 1e−12, 시간 2%와 dense 1e−12는 완화하지 않았다. 완료된 coarse 119단계와 fine 215단계를 재사용하고, 수정의 인과적 영향을 전파할 fine 216~231단계만 다시 풀었다. 저장 재시작의 모든 배열이 정확히 일치했다. 후반 실제 물질 잔차 최댓값={row['maximum_new_true_material_stage']:.12e}, 비선형={row['maximum_true_stage']:.12e}, 합성 물리={row['maximum_true_physical_stage']:.12e}; 원 시간 대조 최댓값={max(physical['time_relative']):.12e}.

분류: Counterexample candidate. 수정된 동일 해의 광자·물질·출구·적용 계량을 전하 원천에 연결했다. 문제였던 fine 220단계의 H 잔차는 {dense['128']['dense_stage'][439][1]:.12e}로 감소했다. 전체 dense 최댓값도 coarse={dense['64']['dense_stage_max']:.12e}, fine={dense['128']['dense_stage_max']:.12e}로 원 1e−12를 통과했다. 다항 원천과 세 GR 대조 경로는 통과했다. 고정 외부와 질량 접합의 독립 감사는 실행했으나 아래 한 항에서 실패했다. 판독에는 물리 적분을 추가하지 않았다. 시작 전 SHA256 키가 sha257로 바뀐 연결 오류는 실패 소스·로그로 보존하고 두 adapter의 키를 고쳤다. 이 연결 오류에서는 물리 계산과 판독 산출물이 변경되지 않았다.

분류: Counterexample candidate. 미세 시계의 compact high={float(c['high']):.12e}, low={float(c['low']):.12e}, low/high={c['low_over_high']:.12e}이다. 세 GR 대조에서 시간 최댓값은 {max(compact['time'].values()):.12e}, 구적 차이는 {compact['controls']['quadrature']:.12e}, 독립 GR 차이는 {compact['controls']['independent_GR']:.12e}였다. 이 범위에서 앞선 조건부 음의 부호가 유지됐다. 같은 해의 고정 외부·질량 정규화 후 명목값은 high={float(e['high']):.12e}, low 증분={float(e['same_solution_low_increment']):.12e}, 합계={float(e['total']):.12e}이지만, 감사 실패로 수락된 결론이 아니다. 서로 다른 크기의 성분을 분리해 읽고 실제 동질 질량 port를 포함했으며, 다른 해의 진단값을 가산하지 않았다.

분류: Counterexample candidate. 외부 감사의 유일한 초과 항은 low 동질 질량의 시간 차이다. coarse/fine 질량은 각각 9.80204685784675e−67/9.609309017999381e−67cm였다. 나머지 정규화 low 전하의 시간 차이 {100*audit['controls']['low']['time']['normalized_standalone']:.8f}%가 통과하더라도, 이를 이유로 독립 질량 기준을 없애지 않는다. 동일 저장 해의 질량 제약을 열에너지·정지에너지·광자·내외부 port로 분해해 원 값을 6.81e−14 이내에서 재현했다. 두 시계의 에너지 차이 {ledger['clock_difference_erg']:.12e}erg 중 기체 비정지 에너지가 {ledger['clock_difference_terms_erg']['gas_nonrest']:.12e}erg를 차지한다. 80자리 입력 합산과의 상대 차이는 최대 {max(v['summation_relative'] for v in ledger['rows']):.12e}로, 2.00574% 차이는 이 합산 반올림으로 설명되지 않는다. 이 분해는 적분·원천 지도의 차이를 분리할 다음 작업의 근거이며, 유일 원인 증명이나 연속 해 오차 상계가 아니다. 기존 경로를 추가 적분하지 않았고 탈락 판정은 보존한다.

분류: Conjectural. 새 fine 계산의 비용은 실행 전 5~60분으로 가정했고 6시간/16GiB/CPU3의 여유를 두었다. 실제 fine receipt는 {cost['seconds']:.2f}초, {cost['peak_RSS_bytes']/1024**3:.3f}GiB였다. 가정 범위는 소폭 넘었지만 승인 상한 내에서 정확도를 유지하여 완료했다. 완료 구간을 반복 중단하거나 원 물리 기간·격자를 확대하지 않았다.

분류: Conjectural. 다음 수락 병목은 같은 해의 기체 비정지 에너지 차이가 지배하는 동질 질량 시간 기준이다. 추가 경로를 자동 확대하거나 2%를 반올림해 통과시키지 않는다. 동시에 완전한 물리 전하를 위해서는 255의 물리 경계 에너지 변환을 실제 외부 전파·전파 중 에너지 변화·광선/도착시간 및 질량/스칼라 접합에 연결하고, 그 계량을 같은 결합 해에 적용해야 한다. 구동 계량에서 경계 에너지를 무한대까지 보존량으로 간주하거나 옛 진단을 이번 전하에 더하지 않는다. 자기 GR 고정점, 배경 재정규화, EOS 인증·균일 미분/시간·공간/경계 오차, 완전 비선형·정적 EFT 비교·관측 연결과 188의 초기 국소 실패도 남는다. 이번 결과는 실제 물질 단계 오차를 수정한 동일 해에서 compact 결론까지 연결한 loophole progress이며, 외부·질량 감사나 전체 연구의 완료가 아니다.
''',encoding='utf-8')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        p=root/f'docs/{name}.md';prefixes[p.relative_to(root).as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
        with p.open('ab') as f:f.write(('\n\n## 단계257 — 물질 성분 정확도 수정과 동일 해 전하\n\n분류: Counterexample candidate. 최종 물리 전하는 미판정이다. 전체 119/231단계 실제 해의 fine 220단계 물질 성분 오차를 더 엄격한 실제 단계 수락으로 수정했다. coarse 119단계와 fine 215단계를 재사용하고 후반 16단계만 다시 풀어, 원 dense 1e−12와 실제 시간·compact 전하 기준을 통과했다. 조건부 compact 음의 부호는 유지됐다. 그러나 고정 외부·질량 감사는 low 동질 질량의 시간 차이 2.005740886%가 원 2%를 넘어 탈락했다. 차이는 주로 기체 비정지 에너지 항에 있으며 합산 반올림으로 설명되지 않는다. 원 256 실패와 연결 키 오류, 새 감사 탈락을 보존한다. 물리 외부 에너지/전파 중 에너지 변화·질량 접합과 실제 계량 반환, 자기 GR 및 EOS/균일 오차/비선형/정적/관측 조건도 남는다. [근거](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',passed=False,verdict='COMPACT_PASS_EXTERIOR_MASS_TIME_FAIL',actual_same_joint_solution_GR_computed=True,
        actual_returned_steps=[119,231],full_period_actual_return_completed=True,
        actual_material_stage_repaired=True,original256_dense_failure_preserved=True,
        original_scientific_gates_preserved=True,internal_material_acceptance_strengthened=True,
        same_solution_dense_compact_charge_completed=True,frozen_exterior_evaluated=True,
        frozen_exterior_mass_audit_passed=False,failed_mass_time_relative=mass_time,
        conditional_compact_sign_survives=compact['conditional_compact_sign_survives_one_return'],
        conditional_frozen_exterior_sign_survives=audit['conditional_compact_sign_survives_frozen_exterior'],
        physical_final_charge_solved=False,final_charge_conclusion='unadjudicated',full_goal_complete=False,snapshot_KST=now)
    write(out/'result.json',final);write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes))
    files=modules+[note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes,sha256={p.relative_to(root).as_posix():sha(p) for p in files});write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_material_accuracy']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
