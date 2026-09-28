"""Publish the actual short same-return charge and its gated long continuation."""
from pathlib import Path
from types import FunctionType
import datetime,importlib.util,sys
spec=importlib.util.spec_from_file_location('b',Path('.phase245-publish.py'));b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha,copy=b.read,b.write,b.sha,b.copy
out=root/'outputs/direct-eos-gr33/native-compensated-return-charge'
manifest=out.parent/'native-compensated-return-charge-manifest.json'
note=root/'notes/REQUEST247_SAME_RETURN_CHARGE_KO.md'

def package():
    dense=runtime/'native-dense-returned246-work';charge=runtime/'native-compensated-charge247-work'
    d=read(dense/'check/result.json');r=read(charge/'check/result.json');poly=read(charge/'check/polynomial-audit.json')
    assert d['representation_controls_passed'] and d['source_time_passed']
    assert read(dense/'check/geometry-result.json')['passed'] and read(dense/'check/symbolic.json')['passed']
    assert r['charge_comparison_admitted'] and r['conditional_compact_sign_survives_one_return'] and poly['passed']
    assert read(charge/'check/controller-status.json')['state']=='completed'
    for folder,actions in [(dense/'check',['prepare','geometry','source']),
        (charge/'check',['prepare','field1288','field648','field1284','collect','audit','compare'])]:
        for action in actions:assert read(folder/f'{action}-receipt.json')['error'] is None
    assert not out.exists();out.mkdir();preserved={p:h for p,h in read(b.manifest)['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    bindings=[]
    for plan in [dense/'check/plan.json',charge/'check/plan.json',charge/'check/polynomial-audit.json',charge/'full-controller-start.json']:
        for p,h in read(plan)['bindings'].items():
            actual=runtime/p.removeprefix('/home/lpaiu/work/native-retained-tail-runtime/')
            assert sha(actual)==h,(str(plan),p)
            bindings.append(dict(plan=str(plan.relative_to(runtime)),source=p,sha256=h))
    for folder,label in [(dense/'check','dense'),(charge/'check','charge')]:
        for p in folder.iterdir():
            if p.is_file() and p.suffix in ['.json','.py','.npz']:copy(p,out/label/p.name)
        for p in (folder/'gr').iterdir():
            if p.is_file() and p.suffix in ['.npz','.json']:copy(p,out/label/'gr'/p.name)
    high=runtime/'native-stage-metric227-work'
    for n in [64,128]:
        copy(high/f'run-{n}.json',out/f'high/run-{n}.json')
        for q in ([8,4] if n==128 else [8]):
            for ext in ['.json','.npz']:copy(high/f'gr/fields-{n}-g{q}{ext}',out/f'high/fields-{n}-g{q}{ext}')
    for name in ['.phase246-followthrough.py','.phase246-launch.ps1','.phase247-followthrough.py',
        '.phase247-polynomial-audit.py','.phase247-full-followthrough.py','.phase247-full-launch.ps1']:
        copy(root/name,out/name.lstrip('.'))
    copy(Path(__file__),out/'publication-producer.py')
    for folder,names in [(dense,['controller-start.json']),(charge,['full-controller-start.json'])]:
        for name in names:copy(folder/name,out/name if folder==charge else out/('dense-'+name))
    snapshots={}
    for number,folder in [(236,'native-complete-return236-work'),(239,'native-common-arithmetic239-work'),
        (244,'native-full-captured244-work'),(245,'native-returned-source245-work'),(246,'native-dense-returned246-work'),(247,'native-compensated-charge247-work')]:
        for name in ['controller-status.json','stage-progress-128.json','capture-128.json','full-controller-status.json']:
            p=runtime/folder/name
            if p.exists():
                value=read(p);snapshots[f'{number}/{name}']=value
                (out/str(number)).mkdir(exist_ok=True);write(out/str(number)/('snapshot-'+name),value)
    now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    fine=r['components'][1]['values']['endpoint_compact_with_metric'];geometry=read(dense/'check/geometry-result.json')
    costs={a:read(charge/f'check/{a}-receipt.json') for a in ['field1288','field648','field1284']}
    note.write_text(f'''# 같은 실제 반환 해의 compact 전하와 긴 구간 후속 실행

분류: Counterexample candidate. **최종 물리 전하는 여전히 미판정이다. 짧은 실제 결합 해에서는 한 번의 GR 반환을 반영한 뒤 기존 compact 전하의 음의 부호가 유지됐다.** 이는 앞서 따로 진화한 응답을 사후 가산한 결과가 아니라, 실제227의 같은15/29단계 보상 풀이에서 높은 성분과 낮은 성분을 각각 읽은 결과다. 시각={now}. 이번 작업은 loophole progress이며 전체 완료가 아니다.

분류: Counterexample candidate. 적용 구간은0.429303889741ms다. fine128의 저장된 compact 전하 정규화에서 높은 성분은{float(fine['high']):.12e}, 같은 해의 낮은 성분은{float(fine['low']):.12e}, 비는{fine['low_over_high']:.12e}다. 거친 경로도 음의 부호를 유지했다. 이 작은 낮은 성분을 잃지 않도록 성분을 분리해서 계산했다. 최종 binary 입력을100자리 Decimal로 합산한 것은 기록 보존 수단이며100자리 물리 정확도를 뜻하지 않는다.

분류: Counterexample candidate. 단계246은 실제 반환 기하의 u와lambda에 독립적인 원천 계수 지도를 사용한다. 최초 입사장의 한 가지 모드 지도를 그대로 쓰지 않았다. 같은 실제 Radau 물질·광자·출구 이력을 사용하고, 같은 질량 제약에서 유도한 floor의 한쪽 극한과 점프를 보존했다. 적용 노드의 값은 정확히 유지했고, 보류한 실제 기하 표본과의 공간L1 차이는 최대{max(geometry['held_out_spatialL1'].values()):.12e}<0.002였다. 이는 표본 검사이며 균일한 시간 또는 미분 오차 보장이 아니다.

분류: Proven. 한쪽 극한을 사용하는 Hermite 값·미분 항등식의 기호 검사를 통과했다. 이는 보간식의 대수적 항등식이며 물리적 시간 오차 정리가 아니다.

분류: Counterexample candidate. 단계247은 높은 성분과 정확히 같은73개 GR 출력 시각, 같은 배경·반경 연산자·potential 표현으로 낮은 원천을 지연 전파했다. 모든 구간의 내부3개 점에서 독립적인 두 기하장 원천 평가와 합성 다항식을 비교했다. 거친228개/미세219개 표본의 최대 상대차는{max(v for row in poly['rows'] for v in row['relative'].values()):.12e}<1e-12였다. 내부 점 검사도 연속 구간 전체의 균일 오차 증명으로 해석하지 않는다.

분류: Counterexample candidate. 원 수락 기준을 그대로 적용했다. 전하 장의 시간 대조 U={r['time']['U']:.12e}<0.02, U_t={r['time']['U_t']:.12e}, U_x={r['time']['U_x']:.12e}, 반경 적분 차수 대조={r['controls']['quadrature']:.12e}<0.002, 독립 Jordan 반경 적분={r['controls']['independent_GR']:.12e}<1e-9였다. U의1.918%는2%기준에 가깝기 때문에 충분한 오차 여유나 최종 부호의 물리 보장을 주장하지 않는다. 기존188의 초기 국소 실패와 이전 산술·표현 실패는 보존한다.

분류: Conjectural. 긴 구간은 실제236→245→246→247의 같은 이력으로 이어진다.236의111/215단계 짝 대조가 통과하면 같은 해의 끝점·연속 원천을 만들고, 독립 다항식 검사 후 세 전하 장을 병렬 계산한다. 짧은 장 계산 실측은{costs['field1288']['seconds']:.3f}/{costs['field648']['seconds']:.3f}/{costs['field1284']['seconds']:.3f}초, 최대RSS는{max(v['peak_RSS_bytes'] for v in costs.values())}bytes였다. 후속 계산에는 CPU0/8/12, 각16GiB·최대2시간, 의존 계산 대기12시간을 배정했다. 실제 긴 입력의 출력 수와 원천 구간 수로 작업량을 추정하고2배 시간 여유+5분 및1.5배 메모리 여유가 이 예산에 드는지 확인한다. 이는 미측정 규모 가정이므로 확정 종료시각은 아니다. 한 작업이 실패하면 남은 동료 작업을 종료하고 원 실패를 보존한다. 기준을 낮추거나 물리 격자를 자동 확대하지 않는다.

분류: Counterexample candidate. 게시 시점에는 짧은 적용만 완료됐고 긴247은 같은246연속 원천을 기다린다. 기존239와236의 물리 진화 및244–246후속 계산을 재시작하지 않았다. 이번 판독에서 새 물질·광자 단계는 없다.

분류: Conjectural. 남은 판정 범위는 전체 선언 기간, EOS/미분, 시간 표현의 균일 오차, 공간·외부 경계, 자기GR 고정점과 완전 비선형 진화, 정적 비교 모형·관측 연결·무한대 전하 정규화다. 짧은 한 번의 반환에서 얻은 부호와 낮은 성분 비를 자기GR 수축 증명이나 전체 목표의 완료로 대체하지 않는다. 연구의 가치는 지배 오차를 해결한 동일 결합 해에서 최종 전하 결론이 유지되는지로 계속 평가한다.
''',encoding='utf-8')
    tails={
        'model-definition':'실제 반환 기하의 독립 u/lambda 지도와 같은 Radau 물질·광자·출구를 사용해 짧은 해의 두 성분 전하를 읽었다. 한 모드 지도와 다른 해의 진단값 가산을 사용하지 않았다.',
        'observable-targets':'최종 전하 미판정. 짧은 같은 반환 해의 compact 전하는 음의 부호를 유지했고 낮은/높은 성분 비는1.63043e-27이었다. 시간 대조1.918%는 원2%기준에 가까우며 무한대 전하나 물리 부호 보장은 아니다.',
        'adiabatic-limit':'한쪽 Hermite 값·미분의 기호 항등식을 확인했다. 같은 질량 제약의 floor 점프와 실제 끝점 값을 유지했다. 표본 보간 검사와 계수 항등식을 균일 EOS/시간 미분 정리로 대체하지 않는다.',
        'nonadiabatic-regime':'같은73개 높은 성분 GR 출력 시각에서 실제 낮은 성분을 전파했다. 원 시간·반경 차수·독립 적분 기준을 통과했고 긴 구간에도 같은 판독을 자동 연결했다.',
        'failure-ledger-dynamic-chi':'최초 구동의 한 모드 기하 계수 지도를 반환 기하에 그대로 적용할 수 없는 병목을 독립 u/lambda 원천 지도와 실제 floor 점프 보존으로 해결했다. 원188국소 실패와 이전 거절 기록 및 수락 기준은 보존한다. 짧은 시간 대조가 기준에 가까운 한계도 유지한다.',
        'dynamic-charge-completion':'짧은 같은 실제 반환 해의 compact 전하 부호가 유지됐지만 최종 물리 전하와 전체 목표는 미판정이다. 긴236→245→246→247자동 연결을 시작했고 실제 입력 크기와 실측 비용으로2시간/16GiB의 세 병렬 판독을 허용한다. 전체 기간·균일 오차·자기GR/비선형·정적/관측·무한대 범위는 그대로 남는다.'}
    prefixes={}
    for name,body in tails.items():
        p=root/f'docs/{name}.md';prefixes[p.relative_to(root).as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
        with p.open('ab') as f:f.write(('\n\n## 단계246–247 — 같은 실제 반환 해의 compact 전하 판독\n\n분류: Counterexample candidate. '+body+' [근거](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',actual_same_joint_solution_GR_computed=True,
        same_actual_returned_solution_read=True,conditional_compact_sign_survives_one_return=True,
        short_charge_result=r,short_dense_result=d,independent_polynomial_audit=poly,
        full_return_charge_completed=False,live_snapshots=snapshots,snapshot_KST=now,
        uniform_temporal_error_certificate=False,self_GR_return_closed=False,physical_final_charge_solved=False,
        final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(out/'result.json',final);write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes,plan_bindings=bindings))
    modules=[root/'verification'/v for v in ['read_dense_returned_source.py','read_compensated_return_charge.py']]
    for p in modules:assert sha(p)==sha(runtime/'verification'/p.name)
    files=modules+[note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes,sha256={p.relative_to(root).as_posix():sha(p) for p in files});write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_compensated_return_charge']={k:v for k,v in final.items() if k!='sha256'};write(master,m)

check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
