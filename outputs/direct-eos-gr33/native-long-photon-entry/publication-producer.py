"""Freeze recovery admission without publishing mutable running outputs."""
from pathlib import Path
from types import FunctionType
import importlib.util,shutil,sys
helper=Path('outputs/direct-eos-gr33/native-radau-pulse-gr/previous-publication-producer.py')
spec=importlib.util.spec_from_file_location('b',helper);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha=b.read,b.write,b.sha
work=runtime/'native-long-photon213-work';out=root/'outputs/direct-eos-gr33/native-long-photon-entry';manifest=out.parent/'native-long-photon-entry-manifest.json'
note=root/'notes/REQUEST213_LONG_SAME_SOLUTION_PHOTON_HISTORY_KO.md';module=root/'verification/recover_remaining_joint_photons.py'


def package():
    assert not manifest.exists();assert read(work/'restart-check.json')['passed']
    assert sha(module)==sha(runtime/'verification'/module.name)==read(work/'controller-start.json')['producer_sha256']
    controller=root/'.phase213-followthrough.py';assert sha(controller)==read(work/'controller-start.json')['source_sha256']
    old=read(out.parent/'native-radau-pulse-gr-manifest.json');preserved={k:h for k,h in old['sha256'].items() if not k.startswith('docs/')}
    for k,h in preserved.items():assert sha(root/k)==h,k
    for k,h in dict(read(work/'plan.json')['bindings'],**read(work/'plan.json')['reused']).items():assert sha(runtime/k.removeprefix('/home/lpaiu/work/native-retained-tail-runtime/'))==h,k
    out.mkdir()
    for name in ['plan.json','restart-check.json','symbolic.json','check-receipt.json','prepare-receipt.json','prepare_retry-receipt.json','preparation-fix.json','failed-prepare-producer.py','controller-start.json']:shutil.copyfile(work/name,out/name)
    shutil.copyfile(controller,out/'controller-producer.py');shutil.copyfile(__file__,out/'publication-producer.py')
    note.write_text('''# 같은 장기 결합 해의 누락 광자 이력 속행

분류: Counterexample candidate. **최종 전하는 미판정이다.** 이 문서는213의 재사용·재시작 검사와 실행 입력을 동결한다. 이후 실행의 종료·수락 단계·전체 기간 또는GR통과를 선언하지 않는다. 연구 분류는 loophole progress다.

분류: Counterexample candidate. 전체 선언 기간3.434431ms중 이미 두 경로에서 수락된15/16=3.219779ms의 동일 물질 단계·광자 끝점·에너지·경계 이력을 사용한다. coarse111/fine215개 기존 단계 중 앞4/8개의205/208광자 복원을 재사용하고 빠진107/207개의 조건부 광자 블록을 복원한다. 원 물질 재적분이나 새 물리 경로가 아니다. 본210의 마지막 기간 속행과 소스·계획은 유지한다.

분류: Counterexample candidate. 저장 단계 시각·가중치·물질 상태·native 및 충돌 변화율의 접두부가 기존 복원 입력과 비트 단위로 같음을 확인했다. 복원된 앞부분의 광자 끝점은 현재 원 저장 끝점과5.179e-17/6.250e-17이내로 일치하며 원1e-12기준을 통과했다. 재시작 체크포인트의 광자 상태·모멘트·충돌·경계 배열 값은 정확히 일치했다. 초기 준비의 출력시각과 누적단계끝점 차이±4.337e-19는 기존1e-18시각 일치 문턱으로 처리했고 실패 근거를 보존했다. 물리 단계 시각·가중치는 바꾸지 않았다.

분류: Conjectural. 각 새 광자 쌍은 원 조건부1e-14/물리1e-13과 전체 결합식1e-12/물리1e-13기준을 통과해야 한다. 실제 native변화율은 원 저장값과 정확히 같아야 한다. 모든 원 저장 출력에서 광자끝점1e-12,동일 충돌의 물질수지1e-8,반경·각도 출구1e-12를 검사한다.205의 추가적인 엄격한H교환률 실패는 유지하고 모멘트를 맞춰 변형하지 않는다. 기존 수락 앞부분을 새 독립 전체식 검사로 세지 않는다.

실행 예산: coarse2시간/fine3시간,각CPU1스레드·가상메모리6GiB로 두 경로를 병렬 실행한다. 실행 전 WSL가용메모리는약50.6GiB였고 기존210계산을 포함한 세 프로세스 상한은18GiB다. 초기 구간 실측으로부터 남은 복원은 coarse36–81분/fine69–155분으로 외삽했으며 후반 속도는 미측정이다. 예산은 완료시각 보증이 아니다. 단계 직전 원자적 체크포인트를 저장하고, 한 경로가 실패하면 다른 소유 복원도 중단한다. 임의 재시작·시간경로·물리격자·기간 확대는 하지 않는다.

분류: Conjectural. 복원이 끝나도 마지막1/16기간은 본210수락 뒤 같은 방식으로 연결해야 한다. 최종 전하는 완성한 동일 해에서 판독하고, 초기T/64국소상대오차40.86%를 전체 최종 전하 오차로 곧바로 확대하지 않는다. 원 시간 실패, 자기GR·무한대전하·전체EOS/균일미분/공간/경계/완전비선형/정적/관측 요건은 모두 유지한다.
''',encoding='utf-8')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        f=root/f'docs/{name}.md';prefixes[f.relative_to(root).as_posix()]=dict(bytes=f.stat().st_size,sha256=sha(f))
        with f.open('ab') as h:h.write(('\n\n## 단계213 — 같은 장기 해의 광자 이력 속행 입력\n\n분류: Counterexample candidate. 최종 전하는 미판정이다. 이미 수락된공통15/16기간의 물질 이력을 재사용해 누락 광자107/207단계 복원을 시작했다. 앞4/8단계 재사용·시각/가중치/상태 접두부·재시작 검사는 통과했다. 이 기록은 실행 입력 동결이며 종료·GR·최종 전하 수락이 아니다. 원 실패와 전체 완료 요건을 유지한다. [근거](../notes/'+note.name+').\n').encode())
    result=dict(classification='Counterexample candidate',entry_evidence_only=True,terminal_result_included=False,actual_same_joint_solution_GR_computed=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(out/'entry-result.json',result);write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes,reused={}))
    files=[module,note]+[root/k for k in prefixes]+[f for f in out.rglob('*') if f.is_file()]
    result.update(document_prefixes=prefixes,sha256={f.relative_to(root).as_posix():sha(f) for f in files});write(manifest,result)
    m=read(master);m['sha256'].update(result['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest);m['native_long_photon_entry']={k:v for k,v in result.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
