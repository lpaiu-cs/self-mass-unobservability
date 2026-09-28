"""Publish actual failures and their bounded repairs, without inferred success."""
from pathlib import Path
from types import FunctionType
import importlib.util,shutil,sys
previous=Path(__file__).with_name('.phase186-publish.py')
if not previous.exists():previous=Path(__file__).with_name('previous-publication-producer.py')
spec=importlib.util.spec_from_file_location('previous',previous);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha=b.read,b.write,b.sha
out=root/'outputs/direct-eos-gr33/native-linear-repair';manifest=out.parent/'native-linear-repair-manifest.json'
folders=['native-return-time189-work','native-full-finish190-work','native-right-precondition191-work','native-linear-floor192-work','native-stable-residual193-work']
names=['resolve_returned_native_time.py','finish_full_incident_horizon.py','repair_full_incident_linear_budget.py','right_precondition_full_interval.py','diagnose_linear_residual_floor.py','stabilize_full_interval_residual.py']
modules=[root/'verification'/n for n in names]
note=root/'notes/REQUEST189_193_FINAL_INTERVAL_REPAIR_KO.md'


def package():
    assert not manifest.exists()
    for p in modules:assert sha(p)==sha(runtime/'verification'/p.name)
    w=[runtime/n for n in folders];diagnosis=read(w[0]/'diagnosis.json');branch=read(w[0]/'branch-result.json')
    assert diagnosis['passed'] and branch['passed'] and read(w[0]/'symbolic.json')['passed']
    original=runtime/'native-full-horizon185-work';assert read(original/'status.json')['state']=='failed'
    original15=read(original/'comparison-15.json');assert original15['passed']
    for folder,receipt in [(w[1],'probe-receipt.json'),(w[1],'repair-coarse-receipt.json'),(w[2],'check-receipt.json'),(w[4],'check-receipt.json')]:
        assert read(folder/receipt)['error'] is not None
    assert read(w[2]/'system-reconstruction.json')['passed']
    arithmetic=read(w[3]/'arithmetic.json');assert read(w[3]/'arithmetic-receipt.json')['error'] is None
    consistent_receipt=read(w[4]/'consistent_check-receipt.json')
    passed=consistent_receipt['error'] is None
    consistent=read(w[4]/'consistent/stable-check-result.json') if passed else dict(passed=False,error=consistent_receipt['error'])
    old=read(out.parent/'native-return-evolution-manifest.json');preserved={p:h for p,h in old['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    aliases={sha(p):p for folder in w for p in folder.rglob('*producer.py')};bindings=[]
    for folder in w:
        for plan in folder.glob('*plan.json'):
            for name,h in read(plan).get('bindings',{}).items():
                actual=runtime/name.removeprefix('/home/lpaiu/work/native-retained-tail-runtime/')
                if sha(actual)!=h:actual=aliases[h]
                assert sha(actual)==h,(plan,name)
                bindings.append(dict(plan=str(plan.relative_to(runtime)),source=name,resolved=str(actual),sha256=h))
    out.mkdir();reused={};known=read(out.parent/'native-return-evolution/publication.json')['reused']
    for folder in w:
        for src in sorted(folder.rglob('*')):
            if not src.is_file():continue
            rel=src.relative_to(folder);key=folder.name+'/'+rel.as_posix()
            if rel.as_posix() in known and sha(src)==known[rel.as_posix()]['sha256']:
                reused[key]=known[rel.as_posix()];continue
            dst=out/folder.name/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);assert sha(src)==sha(dst)
    for src in [original/n for n in ['status.json','plan.json','failure.json','run-receipt.json','comparison-15.json']]:
        dst=out/'original-terminal'/src.name;dst.parent.mkdir(exist_ok=True);shutil.copyfile(src,dst)
    for src,name in [(Path(__file__),'publication-producer.py'),(previous,'previous-publication-producer.py'),(b.helper,'publication-helper.py')]:shutil.copyfile(src,out/name)
    final=dict(classification='Counterexample candidate',linear_system_repair_passed=passed,
        accepted_full_horizon_fraction='15/16',accepted_horizon_seconds=original15['same_horizon_seconds'],
        full_horizon_completed=False,actual_same_joint_solution_GR_computed=False,
        GR_return_time_error_resolved=False,physical_final_charge_solved=False,final_charge_conclusion='unadjudicated',full_goal_complete=False,
        returned_branch_diagnosis=branch,arithmetic_diagnosis=arithmetic,consistent_arithmetic_check=consistent)
    status=('저장된 실패 선형 시스템은 일관된 고정밀 행 평가 뒤 원 선형 기준을 통과했다. 실제 남은 결합 구간에 적용한 수락 결과는 아직 없다.' if passed else '일관된 고정밀 행 평가까지 시험했지만 원 선형 기준을 통과하지 못했다. 실제 남은 결합 구간의 진화는 재개하지 않았다.')
    reduction=1-branch['frozen_branch_quadrature_over_original_difference']/branch['original_quadrature_over_original_difference']
    note.write_text(f'''# 반환 분기 오차와 마지막 결합 구간의 선형 풀이

분류: Counterexample candidate. **최종 전하는 미판정이다.** {status} 앞서 수락된15/16구간,3.219779173ms는 보존되며 원3.434431118ms 전체 기간은 미완료다. 원 계산은 종료됐고 실행 중으로 보고하지 않는다. 이 기록은 진단·풀이의 중간 진전이며 사용자 기준의 최종 과학적 성과를 대신하지 않는다.

분류: Counterexample candidate. GR 반환의 바리온 시간 오차5.211566%는 유지된다. 저장 Radau 단계의 복원 차이 최대{diagnosis['dense_stage_relative_max']:.9g}, 끝점 차이의 분해 잔차{diagnosis['decomposition_relative']:.9g}였다. native 원천의 시간 적분·원 상태 분기 항의 L1은 전체 차이의{diagnosis['baryon_quadrature_and_base_history_L1_over_difference']:.9g}배이고 진화 상태 차이 항은{diagnosis['baryon_evolved_state_L1_over_difference']:.9g}배다. 항들은 상쇄하므로 원인 확률로 읽지 않는다. 마지막 native 분기를 고정한 읽기 전용 대조에서 해당 적분 차이가{100*reduction:.9g}%감소했다. 이는 지배적인 분기 변화의 진단이며 고정 분기로 진화한 물리 해를 수락한 것이 아니다.

분류: Counterexample candidate. 원 전체 입사 해의 마지막 실패를 재현하고 실제 RHS·제안값과 Krylov 이력을 저장했다. 수락된 앞 구간은 다시 계산하지 않았다. 선형 보정만4회에서7회로 늘린 실제 남은 구간 시도는224.924초 뒤 원 벡터 잔차2.455638e-14에서 실패했다. fine 경로는 실행하지 않았다. 7회 시도의 실패 위치는 별도로 저장되지 않아 원 첫 실패와 같은 단계라고 단정하지 않는다.

분류: Counterexample candidate. 저장 실패 시스템의 RHS·제안값을 비트 단위로 복원했다. 같은 전처리기를 우측에 적용하고 원4회 보정을 유지한 시험은49.241초 뒤5.002190e-14로 실패했다. 네 물리 모멘트는 통과했으나 원 벡터 기준1e-14를 완화하지 않았다. GMRES 내부의 작은 전처리 잔차와 실제 방정식 잔차를 구분한다. 이 차이는 [SciPy의 GMRES 정의](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.gmres.html)와 실제 저장 이력을 함께 확인했다.

분류: Counterexample candidate. 실패한 전체 반복값을 저장하기 위한 동일 시스템 재현49.711초와 무반복 고정밀 대조24.595초를 각각 사전 상한75/45초 안에 수행했다. 정규화된 해의 노름1.32977e14와 RHS4.20827e8의 큰 차이에서 바리온 잔차가 지배했다. 원 평가의 바리온 상대 잔차{arithmetic['original_baryon_relative']:.9g}, 같은 저장 값과 정확한 이진 계수의80자리 대조{arithmetic['exact_binary_baryon_relative']:.9g}, 두 평가의 차이{arithmetic['original_evaluation_error']:.9g}였다. 잔차 평가 오차 자체가 기준을 넘으며 저장 해의 정확한 잔차도 아직 기준을 넘는다. 단순 long-double 행렬 선조립도 해결하지 못했다.

분류: Counterexample candidate. 바리온 잔차만80자리로 평가한 첫 수리는73.547초 뒤 실패했다. float 입력의 내부 반복과 long-double 입력의 검사에서 서로 다른 산술 경로를 사용하는 문제가 드러났다. 모든 입력에 같은 바리온 다항식 평가를 적용하는 별도90초 시험을 시행했다. {status} 그 마지막 실행은{consistent_receipt['seconds']:.3f}초였고 원4회 보정,1e-14벡터·1e-13물리 모멘트 기준을 유지했다. 해당 행에 충돌 바리온 원천은 없으며 다른 성분은 기존 연산자를 쓴다. 이 선택적 정밀화는 물리 EOS 인증이나 균일 오차 보장이 아니다.

분류: Proven. Radau 누적 가중치와 끝점 차이의 분해 항등식, 바리온 행의 선조립 다항식 항등식을 검산했다. 기호 항등식은 같은 해의 최종 전하나 시간 수렴을 증명하지 않는다.

분류: Conjectural. 다음 실제 수락 목표는 남은 결합 구간을 원 단계식·시간·에너지·출구 기준으로 완료하는 것이다. 동시에 작은 GR 반환의 시간 분기 오차와 자기 결합 폐쇄가 남아 있다. 이후 같은 해의 에너지·경계 이력으로 전하를 다시 읽어야 한다. 이전 부호를 계승하거나 다른 해의 진단량을 사후 가산하지 않는다. 원 실패·불변 앞부분·수락 기준을 보존한다.
''',encoding='utf-8')
    prefixes={}
    tails={
        'model-definition':'원 전체 입사 해는15/16구간 뒤 선형 풀이에서 종료했다. 저장 실패 시스템과 동일한 RHS·제안값에서 정밀도·전처리를 검사했다. '+status,
        'observable-targets':'최종 전하는 미판정이다. 선형 시스템의 검사나 반환 분기 진단을 실제 전체 결합 해의 전하 수락으로 대체하지 않는다.',
        'adiabatic-limit':'반환 native 분기 고정은 원천 차이를 줄이는 진단이었다. 실제 시간별 분기로 진화한 해의 수렴이나 물리적 정적 흡수를 증명하지 않는다.',
        'nonadiabatic-regime':'반환 바리온 시간 오차5.211566퍼센트와 전체 입사 해 마지막 구간 실패가 모두 남는다. 입력·시간격자·물리 수락 기준을 유지했다.',
        'failure-ledger-dynamic-chi':'4회·7회 보정, 우측 전처리, 선택적 고정밀 잔차의 실패를 보존했다. 바리온의 큰 항 상쇄로 잔차 평가 오차 자체가 원 기준을 넘는 것을 고정밀 대조로 확인했다. '+status,
        'dynamic-charge-completion':'최종 전하 결론은 미판정이다. 원 전체 입사 해는15/16수락 후 종료했다. '+status+' 반환 시간 오차·자기 GR 폐쇄·같은 해의 전하 판독은 계속 남는다.'}
    for name,body in tails.items():
        f=root/f'docs/{name}.md';prefixes[f.relative_to(root).as_posix()]=dict(bytes=f.stat().st_size,sha256=sha(f))
        with f.open('ab') as h:h.write(('\n\n## 단계189–193 — 마지막 구간의 잔차 병목\n\n분류: Counterexample candidate. '+body+' [근거](../notes/'+note.name+').\n').encode())
    write(out/'final-result.json',final);write(out/'publication.json',dict(reused=reused,previous_nondoc=preserved,document_prefixes=prefixes,plan_bindings=bindings))
    files=modules+[note]+[root/k for k in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(sha256={p.relative_to(root).as_posix():sha(p) for p in files},document_prefixes=prefixes);write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_final_interval_linear_repair']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
