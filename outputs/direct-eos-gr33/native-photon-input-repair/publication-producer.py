"""Preserve both failed photon repairs without altering the running actual solver."""
from pathlib import Path
from types import FunctionType
import importlib.util,shutil,sys
helper=Path('outputs/direct-eos-gr33/native-radau-pulse-gr/previous-publication-producer.py')
spec=importlib.util.spec_from_file_location('b',helper);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha=b.read,b.write,b.sha
out=root/'outputs/direct-eos-gr33/native-photon-input-repair';manifest=out.parent/'native-photon-input-repair-manifest.json'
module=root/'verification/restart_photon_recovery_from_archive.py';note=root/'notes/REQUEST220_PHOTON_INPUT_AND_REFINEMENT_KO.md'


def package():
    assert not out.exists();out.mkdir();work=runtime/'native-archival-photon220-work'
    assert sha(module)==sha(runtime/'verification'/module.name)==read(work/'refine-receipt.json')['source_sha256']
    seed=read(work/'seed-identity.json');row=read(work/'clock-128/snapshot-02.json');audit=read(work/'refinement-audit/result.json')
    assert seed['exact_original_checkpoint_photon_input'] and not row['passed'] and not audit['original_port_gate_passed']
    assert read(work/'symbolic.json')['passed'] and read(work/'refinement-audit/receipt.json')['error'] is None
    assert read(work/'fine-receipt.json')['error'] and read(work/'refine-receipt.json')['error']
    previous=read(out.parent/'native-rounding-continuation-manifest.json');preserved={k:h for k,h in previous['sha256'].items() if not k.startswith('docs/')}
    for k,h in preserved.items():assert sha(root/k)==h,k
    for src in work.rglob('*'):
        if not src.is_file() or any(p in ['sweep-0','sweep-1','__pycache__'] for p in src.relative_to(work).parts):continue
        dst=out/'completed'/src.relative_to(work);dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);assert sha(dst)==sha(src)
    for src,name in [(runtime/'native-front-continuation184-work/sweep-1/photons/input-128.npz','original-input-128.npz'),
        (runtime/'native-equation-recovery208-work/recovered-128.npz','reused-prefix-128.npz'),
        (root/'.phase220-refinement-audit.py','refinement-audit-producer.py'),(Path(__file__),'publication-producer.py')]:
        shutil.copyfile(src,out/name);assert sha(src)==sha(out/name)
    note.write_text('''# 광자 구간 입력과 추가 선형 보정의 실제 대조

분류: Counterexample candidate. **최종 전하의 기존 결론은 미판정이다.** 실제 물질·광자 결합218은 별도 고정 소스로 속행하고 있으며, 이 기록은 그 해의 GR 판독에 필요한 광자 이력 복원에서 남은 경계 실패를 다룬다. 본220의 수정은 출구 수락에 실패했으므로 전체 복원 완료나 최종 과학적 성과로 세지 않는다. 연구 분류는 loophole progress의 실패 원인 축소다.

분류: Counterexample candidate.219에서 원 누적 순서와 경계 시각 차이는 원인이 아니었다.220은 원184생산자가T/32에서 실제 사용한restart_x를 같은 구간의 입력으로 복구했다. 원 체크포인트의 물질·시각 접두부와 저장 광자 스냅샷은 정확히 일치했다. 이전 복원 광자 입력과의 전체 상대 차이는5.310380180e-17이었다. 이미 통과한8개 앞 광자 블록의 모멘트·경계 이력은 유지하고, 실패 구간의 나머지8개 조건부 광자 블록만 다시 풀었다. 물질 적분은0개이며 격자·기간·경로·물리식·원 기준을 바꾸지 않았다.

분류: Counterexample candidate.250.96초·최대RSS3.254GB에 같은T/16까지 계산했으나 반경 출구4.547206623e-12>1e-12로 실패했다. 이전4.740709463e-12에서 조금 줄었지만 이 입력 수정으로 해결됐다고 볼 수 없다. 원 전체 결합식·native 비트 동일성은 각 새 쌍에서 통과했다. 최종 광자 끝점1.29074e-16, 물질 수지 최대1.16686e-13, 각도 출구1.70699e-16도 통과했지만 이를 반경 출구 실패의 대체로 쓰지 않았다. 기존214의15단계 입력과 이번220의15단계 입력은 별도 이력이며 서로 섞지 않는다.

분류: Counterexample candidate.실패한 마지막 광자 쌍과 우변을 저장했으므로 추가 대조는 그 쌍에서 시작했다. 원 초기값·우변·잔차를 정확하게 재현하고 내부 선형 목표만1e-18/물리1e-17로 강화했다. 원 물리 수락 기준을 완화한 것이 아니다.11회 추가Krylov보정,123.43초 후 잔차1.992275903e-17에서 강화한 목표에 실패했다. 물리 모멘트는2.92704e-18이하였다. 이 강화 실패도 원220실패와 함께 보존했다.

분류: Counterexample candidate.추가 반복을 더하지 않고 저장된 마지막 보정의 실제 영향을 별도로 재평가했다. 해의 상대 변화는5.07318e-19이며 원 전체 결합식9.55195e-17와native동일성은 통과했다. 그러나 반경 출구는4.545628262e-12로 여전히 실패하고 바깥 광자 수 출구의 변화는3.81917e-14에 그쳤다. 같은 마지막 블록에 반복 수만 더 주는 것은 이번 경계 차이를 실질적으로 줄이지 못했다. 이 결과를 원 전체 광자 풀이의 최소 달성 오차 증명이나 모든 앞 구간의 오차 배제로 확대하지 않는다.

분류: Proven. 저장 기체를 알고 있는 항으로 옮겨 같은 선형 광자 블록을 푸는 정확 실수 항등식을 기호 검사로 다시 확인했다. 이는 부동소수점 복원의 정확 일치, EOS 균일 오차 또는 물리 해의 정확도 증명이 아니다.

실행은각6GiB·CPU1스레드,8블록900초·저장쌍 보정900초·재평가300초의 사전 예산 안에서 종료했다. 두 실패 모두 벽시간/메모리 부족으로 발생한 것이 아니다. 실제 결합218의2/3시간 예산과 실행 소스는 손대지 않았다. 원220완료 직전 쌍, 강화 실패 쌍, 각각의 마지막 수락 입력·이력과 모든 원 기준을 보존했다.

분류: Conjectural. 남은 경계 차이는 앞 구간에서 누적된 광자 근사, 저장 기체 좌표의 역변환과 충돌 평가, 원 결합 풀이가 남긴 광자 오차를 분리해야 한다. 현재 데이터로 유일 원인을 선언하지 않는다. 원 기체·시각·출구를 사후 맞춰서 통과시키거나 다른 해의 진단을 최종 전하에 더하지 않는다. 최종 목표는 지배 오차를 해결한 같은 결합 해와 그 에너지·경계 이력에서 최종 전하의 결론을 판정하는 것이며, 전체 시간 대조·자기GR·EOS/미분·공간·경계·비선형·정적/관측·무한대 전하 조건은 그대로 남는다.
''',encoding='utf-8')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        f=root/f'docs/{name}.md';prefixes[f.relative_to(root).as_posix()]=dict(bytes=f.stat().st_size,sha256=sha(f))
        with f.open('ab') as h:h.write(('\n\n## 단계220 — 원 광자 입력과 추가 반복의 한계\n\n분류: Counterexample candidate. 최종 전하는 미판정이다. 원 체크포인트의 광자 입력으로 실패 구간8블록을 다시 풀어도 반경 출구4.5472e-12는 원1e-12기준을 실패했다. 저장된 마지막 선형쌍만11회 더 보정한 뒤에도 출구4.5456e-12로 거의 바뀌지 않았다. 원 전체 결합식·native동일성·에너지/물질수지 통과로 이 경계 실패를 대체하지 않는다. 원 실패와 강화목표 실패를 모두 보존하고 동일 반복을 자동 확대하지 않았다. 실제 결합218의 고정 실행은 별개로 유지한다. [근거](../notes/'+note.name+').\n').encode())
    result=dict(classification='Counterexample candidate',actual_same_joint_solution_GR_computed=False,
        exact_original_photon_input_replayed=True,conditional_blocks_recomputed=8,original_port_failure_preserved=True,
        original_equation_after_extra_corrections_passed=audit['original_equation']['passed'],last_refinement_did_not_resolve_port=True,
        numerical_result=audit,new_physical_steps=0,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(out/'result.json',result);write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes,reused={}))
    files=[module,note]+[root/k for k in prefixes]+[f for f in out.rglob('*') if f.is_file()]
    result.update(document_prefixes=prefixes,sha256={f.relative_to(root).as_posix():sha(f) for f in files});write(manifest,result)
    m=read(master);m['sha256'].update(result['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest);m['native_photon_input_repair']={k:v for k,v in result.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
