"""Publish209 and immutable210 entry evidence; never snapshot a live result."""
from pathlib import Path
from types import FunctionType
import importlib.util,shutil,sys
p=Path(__file__).with_name('.phase186-publish.py')
if not p.exists():p=Path(__file__).with_name('previous-publication-producer.py')
spec=importlib.util.spec_from_file_location('p',p);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha=b.read,b.write,b.sha
out=root/'outputs/direct-eos-gr33/native-momentum-entry';manifest=out.parent/'native-momentum-entry-manifest.json'
folders=[runtime/n for n in ['native-momentum209-work','native-momentum-continuation210-work']]
modules=[root/'verification'/n for n in ['repair_native_momentum_residual.py','continue_precise_momentum.py']]
note=root/'notes/REQUEST209_210_ACTUAL_MOMENTUM_CONTINUATION_KO.md'


def package():
    assert not manifest.exists()
    audit,entry=folders;r=read(audit/'result.json');restart=read(entry/'restart-check.json')
    assert r['exact_saved_rhs_and_residual'] and r['precision_difference']==0 and not r['linear_gates_passed']
    assert restart['passed'] and restart['accepted_actual_steps']==115 and restart['new_physical_steps']==0
    assert read(entry/'linear-seed-identity.json')['passed']
    assert all(sha(m)==sha(runtime/'verification'/m.name) for m in modules)
    controller=root/'.phase210-followthrough.py';assert sha(controller)==read(entry/'controller-start.json')['source_sha256']
    old=read(out.parent/'native-flux-precision-manifest.json')
    preserved={name:h for name,h in old['sha256'].items() if not name.startswith('docs/')}
    for name,h in preserved.items():assert sha(root/name)==h,name
    known={h:dict(path=name,sha256=h) for name,h in read(master)['sha256'].items() if not name.startswith('docs/')}
    known.update({v['sha256']:v for v in read(out.parent/'native-flux-precision/publication.json')['reused'].values()})
    bindings=[]
    for folder in folders:
        for name,h in read(folder/'plan.json')['bindings'].items():
            actual=runtime/name.removeprefix('/home/lpaiu/work/native-retained-tail-runtime/')
            assert sha(actual)==h,name
            bindings.append(dict(plan=folder.name+'/plan.json',source=name,sha256=h))
    out.mkdir();reused={}
    # Only files that the running continuation never rewrites are copied from210.
    frozen=['plan.json','symbolic.json','prepare-receipt.json','check-receipt.json','restart-check.json','check-result.json','controller-start.json','linear-seed-identity.json']
    sources=[s for s in audit.rglob('*') if s.is_file()]+[entry/n for n in frozen]
    for src in sorted(sources):
        key=src.relative_to(runtime).as_posix();digest=sha(src)
        if digest in known and (root/known[digest]['path']).exists():
            assert sha(root/known[digest]['path'])==digest;reused[key]=known[digest];continue
        dst=out/key;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);assert sha(dst)==sha(src)==digest
        known[digest]=dict(path=dst.relative_to(root).as_posix(),sha256=digest)
    for src,name in [(Path(__file__),'publication-producer.py'),(p,'previous-publication-producer.py'),(b.helper,'publication-helper.py'),(controller,'followthrough-producer.py')]:shutil.copyfile(src,out/name)
    note.write_text('''# 실제 운동량 산술 수정과 결합 적분 속행

분류: Counterexample candidate. **최종 전하는 미판정이다.** 209의 저장 실패 시스템 분석을 마치고, 그 수정을210의 실제 광자·물질 결합 풀이에 연결했다. 이 문서의210 근거는 재시작 검사와 실행 입력을 동결한 기록이다. 실행 종료·116번째 단계 수락·전체 기간 통과를 선언하는 문서가 아니다. 연구 분류는 loophole progress다.

분류: Counterexample candidate. 202는115번째 실제 단계를 수락한 뒤, 다음 선형 풀이에서1.794402994e-14로 원1e-14기준을 실패했다. 4181.741초로7200초 상한 전에 종료됐으므로 예산 부족에 의한 중단이 아니다. 지배 잔차는 반경 운동량 셀261이며, 기존 바리온 전용 보정이 처리하는 성분과 다르다. 기존 상태·실패·기준을 유지했다.

분류: Counterexample candidate. 209는 저장 우변·분기 제안·기존 잔차를 비트 단위로 재현했다. native 운동량 행의 I-hAJ를 이진 계수 그대로 고정밀로 조립하고 같은 충돌 항과 결합했다. 40/80자리 평가 차이는 저장 longdouble에서0이었다. 기존 산술과의 차이는 상대5.587227623e-14로 원 선형 기준보다 크다. 고정밀 평가에서 저장 해의 선형 잔차는4.874143907e-14, 실제 비선형 단계 잔차는9.950524852e-12여서 둘 다 실패다. 물리 모멘트가 통과한 것을 전체 벡터 기준의 대체로 삼지 않았다.

분류: Proven. 운동량 행의 I-hAJ 사전 조립과 같은 충돌 항을 합한 식은 원 선형식과 대수적으로 같다. 기호 검산을 통과했다. 이 항등식은 EOS·미분·시간 오차의 물리 인증이 아니다.

분류: Counterexample candidate. 210은 기존60자리 실제 native 질량 유속과 에너지 기준 변환을 유지하고, 모든 Krylov 및 잔차 입력에 동일한80자리 native 운동량 행 평가를 적용한다. 기존 제한된 잔차 보정의 대상도 B 또는 S 성분으로 확장했다. 충돌 자체의 모멘트 계산은 기존 longdouble이므로 모든 결합 연산을 임의 정밀도로 바꿨다고 주장하지 않는다. 최종 수락에는 원 실제 단계식·벡터·물리·보존·출구·시간 기준을 그대로 요구한다.

분류: Counterexample candidate. 115개 수락 상태와 전체 이력·수지를 정확히 복원하는 검사를 통과했고, 재시작 검사에서 새 물리 적분은0회였다. 남은 coarse4/fine16단계만 속행한다. 실패한202분기 제안과 누적 선형 제안은 우변·분기 제안의 정확한 일치를 확인한 뒤 초기 제안으로만 재사용했다. 새 수락 상태로 간주하지 않았다.

계산 정책: coarse7200초,fine10800초,각 가상 메모리6GiB·CPU1스레드,새 Newton8회·선형 보정12회다. 사용자 지시에 따라 충분한 속행 예산을 유지한다. 실행 중 소스·계획은 고정한다. 승인된 같은 두 경로가 원 대조를 통과하면 같은 에너지·경계 이력의 GR 읽기로 자동 연결한다. 이것은 종료시각 예측이나 무제한 격자·기간·경로 확대가 아니다.

분류: Conjectural. 수정 적용이 남은 실제 단계를 수락하게 하는지, 지배적인 GR 시간 표현·상태 오차를 고친 같은 자기 결합 해에서 최종 전하의 결론이 유지되는지는 실행 결과로 판정해야 한다. 초기 원천·GR 시간 실패와 원 작은 반환5.211566% 실패는 유지한다. 무한대 전하·전체 물리 EOS·균일 미분·공간·경계·완전 비선형·정적/관측 연결은 여전히 미완료다.
''',encoding='utf-8')
    tails={
        'model-definition':'native 운동량 I-hAJ와 원 충돌 항의 일관된 고정밀 평가를 실제210결합 풀이에 연결했다. B/Eref 수정은 유지한다. 충돌 모멘트 자체까지 고정밀화한 것은 아니다.',
        'observable-targets':'최종 전하는 미판정이다. 115개 수락 상태에서 실제 속행을 시작했으며 여기의210근거는 실행 입력 동결이다. 다음 단계나 최종 전하의 수락 결과로 세지 않는다.',
        'adiabatic-limit':'운동량 행 사전 조립 항등식과40/80자리 대조를 확인했다. 물리 EOS·균일 미분·시간·경계 오차 상계가 아니다.',
        'nonadiabatic-regime':'저장 해는 수정 산술에서도 선형4.87414e-14와 실제 단계9.95052e-12로 실패한다. 우변·분기 제안의 정확한 일치 뒤 실제210재풀이의 초기 제안으로만 재사용한다.',
        'failure-ledger-dynamic-chi':'202후속 실패는 시간 상한 전의 운동량 벡터 기준 미달이다. 물리 모멘트 통과로 대체하지 않는다. 일관된 S산술과 B/S제한 보정을 적용한210의 재시작 검사는 통과했지만 다음 단계 수락은 이 기록에서 판정하지 않는다.',
        'dynamic-charge-completion':'최종 전하는 미판정이다. 115단계의 같은 상태·전체 이력을 재사용하여 실제 운동량 수정 풀이를 시작했다. 이 동결 기록을 완료로 세지 않으며 GR시간오차·자기결합·무한대전하·전체EOS/미분/공간/경계/비선형/정적/관측 요건을 유지한다.'}
    prefixes={}
    for name,body in tails.items():
        f=root/f'docs/{name}.md';prefixes[f.relative_to(root).as_posix()]=dict(bytes=f.stat().st_size,sha256=sha(f))
        with f.open('ab') as h:h.write(('\n\n## 단계209–210 — 운동량 수정과 실제 속행 입력 동결\n\n분류: Counterexample candidate. '+body+' [근거](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',momentum_audit=r,restart=restart,entry_evidence_only=True,terminal_evolution_result_included=False,original_failures_preserved=True,
        actual_same_joint_solution_GR_computed=False,physical_final_charge_solved=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(out/'entry-result.json',final);write(out/'publication.json',dict(reused=reused,previous_nondoc=preserved,document_prefixes=prefixes,plan_bindings=bindings))
    files=modules+[note]+[root/k for k in prefixes]+[f for f in out.rglob('*') if f.is_file()]
    final.update(sha256={f.relative_to(root).as_posix():sha(f) for f in files},document_prefixes=prefixes);write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest);m['native_momentum_entry']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
