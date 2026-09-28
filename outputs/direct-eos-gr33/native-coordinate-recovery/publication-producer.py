"""Freeze terminal evolution evidence and the corrected recovery admission."""
from pathlib import Path
from types import FunctionType
import importlib.util,shutil,sys
helper=Path('outputs/direct-eos-gr33/native-radau-pulse-gr/previous-publication-producer.py')
spec=importlib.util.spec_from_file_location('b',helper);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha=b.read,b.write,b.sha
out=root/'outputs/direct-eos-gr33/native-coordinate-recovery';manifest=out.parent/'native-coordinate-recovery-manifest.json'
note=root/'notes/REQUEST213_215_COORDINATE_RECOVERY_AND_STAGE_LIMIT_KO.md'
module=root/'verification/resume_coordinate_exact_photons.py'


def package():
    assert not manifest.exists();out.mkdir()
    w210=runtime/'native-momentum-continuation210-work';w213=runtime/'native-long-photon213-work';w214=runtime/'native-coordinate-photon214-work';w215=runtime/'native-broad-polish215-work'
    assert read(w210/'controller-status.json')['state']=='failed' and read(w210/'failure-64.json')['actual_accepted_steps']==116
    check=read(w214/'restart-check.json');assert check['passed'] and sum(r['stages'] for r in check['rows'])==652
    assert read(w215/'result.json')['exact_RHS_and_residual'] and not read(w215/'basis/result.json')['linear_passed']
    assert sha(module)==sha(runtime/'verification'/module.name)==read(w214/'resume-plan.json')['producer_sha256']==read(w214/'resume-controller-start.json')['producer_sha256']
    old=read(out.parent/'native-long-photon-entry-manifest.json');preserved={k:h for k,h in old['sha256'].items() if not k.startswith('docs/')}
    for k,h in preserved.items():assert sha(root/k)==h,k
    def copy(src,dst):
        dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);assert sha(src)==sha(dst)
    for folder,names in [
        (w210,['controller-status.json','coarse-receipt.json','failure-64.json','stage-progress-64.json','linear-64.json','polish.json','last-accepted-64.npz','failed-linear-64.npz','controller.stderr.log']),
        (w213,['controller-status.json','fine-receipt.json','progress-64.json','progress-128.json','accepted-64.npz','accepted-128.npz','rejected-original-128.npz','clock-128/original-equation-9.json','native-replay.json','inverse-probe.json']),
        (w214,['plan.json','prepare-receipt.json','check-receipt.json','check_interval-receipt.json','failed-check-producer.py','array-guard-failed-producer.py','time-probe.json','owner-probe.json','cache-probe.json','interval-amendment.json','symbolic.json','restart-check.json','controller-start.json','controller-status.json','coarse-receipt.json','coarse.stderr.log','resume-plan.json','resume-input-64.npz','resume-input-128.npz','resume-controller-start.json']),
        (w215,['plan.json','receipt.json','result.json','polish.json','proposal.npz','setup-failed-producer.py','setup-failed-plan.json','setup-failed-receipt.json','closure-failed-producer.py','closure-failed-plan.json','closure-failed-receipt.json','basis/plan.json','basis/receipt.json','basis/result.json','basis/polish.json','basis/proposal.npz'])]:
        for name in names:copy(folder/name,out/folder.name/name)
    for p in (w214/'first-dispatch').rglob('*'):
        if p.is_file():copy(p,out/w214.name/p.relative_to(w214))
    producers=['.phase213-native-replay.py','.phase213-inverse-probe.py','.phase214-time-probe.py','.phase214-owner-probe.py','.phase214-cache-probe.py','.phase214-followthrough.py','.phase214-resume-followthrough.py','.phase214-resume-prepare.py','.phase215-polish-probe.py','.phase215-basis-probe.py']
    for name in producers:copy(root/name,out/'producers'/name.removeprefix('.'))
    copy(Path(__file__),out/'publication-producer.py')
    note.write_text('''# 동일 해의 좌표 복원 수정과 후속 실제 단계 경계

분류: Counterexample candidate. **최종 전하의 기존 결론은 미판정이다.** 실제 결합 계산210은 종전 실패 지점인116번째 단계를 원 기준으로 통과했다. 이후117번째 선형 풀이에서 실패했다. 별도의 광자 이력 복원은 저장 좌표와 원 구간별 모델 수명을 복구해 같은 해에서 속행한다. 연구 분류는 loophole progress이며 최종 과학적 완료가 아니다.

분류: Counterexample candidate. 실제116단계의 결합식 잔차는9.4399304441e-14<1e-12, 물리 모멘트 최대1.23351e-17<1e-13이다. 수락 시각은3.353936639ms이며 coarse 전체119단계 중116개다. 약3889.54초 후117단계의 선형 잔차4.093584129e-10>1e-14로 종료했다.2시간 한도보다 일찍 끝났으므로 자원 제한이 원인은 아니다.116수락 상태·모든 이력과 실패117제안을 보존했고 fine 전체 기간 속행은 시작하지 않았다.

분류: Counterexample candidate. 저장117제안의 우변과 잔차를 비트 단위로 재현했다. 기존 잔차 보정의 적용 문턱만1e-10에서1e-8로 넓히자 선형 잔차가3.301275091e-13까지 감소했지만 원1e-14기준은 실패했다. 실제 비선형 잔차도3.788254248e-7로 실패했다. 보정 후보에서 바리온 열의 일괄 제외를 없애고 기존 계수×ULP 문턱으로 각각 걸러도 추가 개선은 없었다. 어떤 제안도 수락하지 않았고 새로운 물리 단계·GMRES 반복은0개다. 작은 물리 모멘트 잔차로 원 벡터 기준을 대체하지 않는다.

분류: Conjectural. 남은 선형 잔차는 큰 정규화 변수의 표현 간격과 관련된 것으로 의심된다. 바리온 보정 뒤 운동량 성분이 제한하고 제한된 열 보정은 정체했다. 더 많은 동일 반복만으로 해결된다고 가정하지 않는다. 다음 실제 진화 수정은 저장 제안에서 작은 보정분을 잃지 않는 표현과 원 식 평가를 검증한 뒤 적용해야 하며, 실패한 선형·비선형 제안을 그대로 통과시키면 안 된다.

분류: Counterexample candidate. 광자 복원213은 coarse6/fine9단계까지 진행한 뒤 fine10번째 쌍에서 native 저장값 비트 일치가 실패했다. 그 쌍의 전체 결합식 잔차5.36494e-17과 물리 모멘트 최대1.70506e-16은 원 기준을 통과했지만 별도의 비트 일치 실패를 보존했다. 한 셀의 저장 Eref 역변환이 원 보존량과0.5만큼 달랐고 이웃 셀 운동량 변화율 차이1024는 해당 값4ULP, 전체 운동량 변화율 L1의1.46579e-28이었다. Jacobian 호출 전후 차이는 없었다. 저장 보존량을 정확히 재현하는 이웃 부동소수점 좌표를 찾자 native 변화율도 정확히 복원됐다. 변화율이나 광자 모멘트에 맞춰 좌표를 적합하지 않았다.

분류: Counterexample candidate. 단일 모델을 전체 이력에 계속 쓰는 추가 오류도 확인했다. coarse17번째 저장 단계에서 원 native 값과 달랐지만 원 생산자처럼 해당 구간 시작에서 모델을 새로 구성하면 정확히 재현됐다. 시각의 자료형 변경이나 열역학 미분 캐시만 지우는 방법은 해결하지 못했다. 원185생산자의 정규 출력마다 모델을 재생성하는 수명을 그대로 복원했다. 보존량 역변환과 이 수명을 함께 적용한 결과652개 저장 단계(222coarse/430fine)의 모든 보존량과 native 값이 비트 단위로 일치했다. 바뀐 정규화 좌표는8/13개이며 물리 이력 자체는 바꾸지 않았다. 이 검사는449.79초, 최대RSS약1.572GB였고 물질 재적분·광자 풀이0개다.

분류: Proven. B의 정규화 역변환과Eref=Etilde+kappa*B의 정확 실수 좌표 항등식은 기호 계산으로 확인했다. 유한 정밀도의 정확 복원은 위652개 저장 상태에 대한 별도 수치 확인이며 균일 EOS·미분 오차 정리가 아니다.

분류: Counterexample candidate. 수정214를 실제 광자 복원에 적용해 이전 fine 실패 쌍과 다음 쌍, coarse의 다음 두 단계를 원 결합식·native 비트 일치 기준으로 통과했다. coarse T/16출력에서 광자 끝점1.67808e-17, 같은 해의 물질 수지 최대3.14468e-15, 반경 출구1.74092e-13도 원 기준을 통과했다. 이후 원 구간 전환의 배열 비교를 scalar처럼 처리한 구현 오류로 멈췄으며 이 실패도 보존했다. 배열 동일성 비교로 고치고 수락 직후에도 원자적 체크포인트를 저장하도록 수정했다. 역변환·모델 수명 함수의 AST가652단계 검사 때와 같음을 확인했다. 마지막으로 저장된coarse7/fine11복원 단계에서 재개한다. coarse8번째 광자 끝점은 오류 직전에 저장되지 않아 그 조건부 블록 하나만 다시 계산한다. 수락된 물질 단계를 다시 적분하지 않는다.

실행 예산: 사용자 지시에 따라 coarse2시간·fine3시간, 각CPU1스레드·가상메모리6GiB로 병렬 속행한다. 원 물리 격자·시간 경로·기간과 조건부1e-14/물리1e-13, 원 결합식1e-12/물리1e-13, native비트 동일성, 끝점·출구1e-12, 물질 수지1e-8은 그대로다. 구간당 생성자 약14초의 비용을 추가로 허용했다. 이전 중단 기록과 재개 입력의 해시를 고정했다. 이 문서의214부분은 수락된 앞부분·재개 입력 증거이며 이후 전체 복원의 완료나 수락을 선언하지 않는다.

분류: Conjectural. 전체 복원과 마지막 결합 기간, 두 시간 경로 비교를 통과한 동일 해에 한해 Radau 물질·광자·실제 입사 펄스·경계 이력을 GR 최종 판독으로 연결한다. 초기T/64의 국소 GR 상대오차40.86%는 전체 최종 전하 오차로 확대하지 않는다. 자기GR·무한대전하·EOS/균일미분/공간/경계/완전비선형/정적/관측 완료 조건과 원 실패들은 유지한다.
''',encoding='utf-8')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        f=root/f'docs/{name}.md';prefixes[f.relative_to(root).as_posix()]=dict(bytes=f.stat().st_size,sha256=sha(f))
        with f.open('ab') as h:h.write(('\n\n## 단계213–215 — 동일 좌표·연산자 복원과 실제 단계 경계\n\n분류: Counterexample candidate. 최종 전하는 미판정이다. 실제116단계는 원 기준을 통과했지만117선형 풀이와 저장 제안의 실제 비선형식은 실패했다. 넓힌 잔차 보정만으로 수락하지 않았다. 저장652단계의 정확 보존좌표와 원 구간별 모델 수명을 복구했고 실제 광자 재개에 적용해 기존 실패 쌍을 통과했다. 배열 비교 오류를 고친 뒤 기존 체크포인트에서 속행한다. 이 기록은 전체 복원·자기GR·최종 전하 통과가 아니다. 원 실패·기준·전체 완료 조건을 유지한다. [근거](../notes/'+note.name+').\n').encode())
    result=dict(classification='Counterexample candidate',actual_same_joint_solution_GR_computed=False,terminal210accepted_steps=116,terminal210failed_next_linear=True,all652conserved_and_native_exact=True,recovery214_entry_and_first_dispatch_only=True,recovery214_terminal_resume_result_included=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(out/'result.json',result);write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes,reused={}))
    files=[module,note]+[root/k for k in prefixes]+[f for f in out.rglob('*') if f.is_file()]
    result.update(document_prefixes=prefixes,sha256={f.relative_to(root).as_posix():sha(f) for f in files});write(manifest,result)
    m=read(master);m['sha256'].update(result['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest);m['native_coordinate_recovery']={k:v for k,v in result.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
