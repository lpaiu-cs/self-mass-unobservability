"""Preserve the original coupled capture and its actual continuation entry."""
from pathlib import Path
from types import FunctionType
import importlib.util,shutil,sys
helper=Path('outputs/direct-eos-gr33/native-radau-pulse-gr/previous-publication-producer.py')
spec=importlib.util.spec_from_file_location('b',helper);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha=b.read,b.write,b.sha
out=root/'outputs/direct-eos-gr33/native-original-photon-capture';manifest=out.parent/'native-original-photon-capture-manifest.json'
note=root/'notes/REQUEST221_223_ORIGINAL_PHOTON_CAPTURE_KO.md'
modules=[root/'verification'/name for name in ['capture_original_joint_photons.py','continue_captured_photon_history.py']]


def copy(src,dst):
    before=sha(src);dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);assert before==sha(dst)==sha(src),src


def package():
    assert not out.exists();out.mkdir()
    control=runtime/'native-original-state221-work';capture=runtime/'native-joint-capture222-work';entry=runtime/'native-captured-photon223-work'
    result=read(capture/'result.json');replay=read(capture/'physical-replay.json')
    assert result['passed'] and replay['passed'] and result['all_physical_array_values_exact']
    assert result['snapshot']['radial_port_relative']<1e-12 and result['snapshot']['endpoint_relative']==0
    assert read(capture/'symbolic.json')['passed'] and read(capture/'fine-receipt.json')['error'] is None
    assert read(control/'receipt.json')['error'].startswith("AssertionError(('Original endpoint collision identity'")
    assert read(entry/'prepare-receipt.json')['error'] is None and read(entry/'restart-check.json')['passed']
    assert sha(modules[0])==sha(runtime/'verification'/modules[0].name)==read(capture/'fine-receipt.json')['source_sha256']
    assert sha(modules[1])==sha(runtime/'verification'/modules[1].name)==read(entry/'controller-start.json')['producer_sha256']
    previous=read(out.parent/'native-photon-input-repair-manifest.json')
    preserved={p:h for p,h in previous['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    reused={}
    for folder,name in [(control,'endpoint-control'),(capture,'original-capture')]:
        for src in folder.rglob('*'):
            if not src.is_file() or '__pycache__' in src.parts:continue
            rel=src.relative_to(folder)
            if rel.parts[0]=='sweep-0' or (rel.parts[0]=='sweep-1' and not src.name.startswith('capture-128')):
                reused[f'{name}/{rel.as_posix()}']=dict(runtime_path=str(src),sha256=sha(src));continue
            copy(src,out/name/rel)
    for name in ['plan.json','restart-check.json','prepare-receipt.json','controller-start.json','input-64.npz','input-128.npz']:
        copy(entry/name,out/'continuation-entry'/name)
    # Freeze actual accepted progress, not a live controller file in the manifest.
    for n in [64,128]:
        p=entry/f'progress-{n}.json';row=read(p);assert row['completed']>{64:11,128:16}[n]
        copy(p,out/'continuation-entry'/p.name)
    for src,name in [(root/'.phase221-original-state.py','endpoint-control-producer.py'),
                     (root/'.phase223-followthrough.py','continuation-controller.py'),(Path(__file__),'publication-producer.py')]:copy(src,out/name)
    note.write_text('''# 원 결합 광자 단계의 직접 회수와 후속 적용

분류: Counterexample candidate. **최종 전하의 기존 결론은 미판정이다.** 이번에는 기존 실패 구간의 실제 결합 단계를 그대로 재현하고, 누락됐던 광자 상태와 경계 이력을 직접 저장하여 원 출구 기준을 통과했다. 연구 분류는 loophole progress의 연결 병목 해소다. 전체 기간·GR 시간 오차·자기 결합·최종 전하의 수락은 아직 아니다.

분류: Counterexample candidate.221 끝점 대조는 원 충돌 비트 동일성 검사에서 실패했다. 19.34초, 최대RSS2.461GB를 사용했다. 보존좌표 동일성과 선택된 guide의 재역변환 동일성은 확인했지만, 입력 guide 자체가184에서 보존량을 역변환해 만든 값이었다. 따라서 독립적으로 저장한 원 기체 좌표와의 대조로 해석할 수 없으며, 기체 역변환을 원인에서 배제하지 않는다. 원 충돌 비트 실패와 잘못된 출처 가정의 한계를 보존한다.

분류: Counterexample candidate.222는 동일한 원184 초기화·restored-128 입력·단계 함수·방정식·3회Newton 규칙으로 원 구간8개 실제 결합 단계를 재생했다. 수락된 광자·기체 쌍을 반환 직후 복사하고, 기존 경계 함수가 실제 반환한 값을 한 번 저장했다. 추가 경계 호출이나 원천·풀이 변경은 없었다. 처음8개 광자 복원 이력은 재사용했다. 새 물리 기간이나 격자를 계산하지 않았다.

분류: Counterexample candidate.전체213.20초, 최대RSS2.568GB에서 원184의 물질·광자·native·충돌·바닥 보정·에너지·각도/반경 출구 이력과 내부 끝점 배열을 모두 정확히 재현했다. 원 단계 생성 소스도 정확히 같았다. 같은T/16 끝점에서 광자 상대 차이는0, 물질 수지 최대1.17724e-13, 각도 출구6.65773e-18, 반경 출구1.61324e-13<1e-12로 통과했다. 기존 조건부 복원214/220의4.74e-12/4.55e-12 실패는 그대로 유지한다. 사후 출구 정합이나 다른 해의 진단 가산으로 통과시킨 것이 아니다. 실제 원 근사해의 단계 기록을 직접 확보했을 때 해당 실패를 해결했다는 제한된 결과다.

분류: Proven. 기존 Radau 가중치의0·1·2차 모멘트 항등식을 유리수로 검사했다. 이는 원 단계 관찰 코드의 물리 수렴이나 부동소수점 균일 오차를 증명하지 않는다.

222 준비 후 실행 전에 생성함수의inspect 조회를 이미 저장된 원 단계 소스의 정확 비교로 고쳤다. 준비 당시 소스·계획은 별도 보존하고 수정 계획에 양쪽SHA와 이유를 기록했다. 실제 재생 실행 중 소스와 계획은 바꾸지 않았다.

분류: Counterexample candidate.223은 수정된 미세 경로16단계 체크포인트와 기존 거친 경로11단계 체크포인트를 실제 후속 복원에 적용했다. 추가 조건부 광자 단계를 두 경로 모두 수락한 시점의 진행 기록을 동결했다. 남은100/199블록은 같은 원 결합식·native 비트 일치·끝점·수지·각도/반경 출구 기준으로 검증한다. 원652개 native/보존좌표 대조는214의 검증을 계승하며 이번에 새로 수행했다고 세지 않는다. 실패한 스냅샷의 전체 광자·기체 쌍과 경계 이력을 즉시 저장하도록 하여 같은 누락 때문에 재계산하지 않게 했다.

자원은각CPU1스레드·6GiB, 거친 경로7200초·미세 경로10800초이며 독립 두 경로를 병렬 실행한다. 이전250.96초/8조건부 블록의 실측을 바탕으로20~45초/블록을 가정하면33~75분/66~149분이다. 후반 속도는 미측정이며 종료 보장은 아니다. 원 기준 실패 시 다른 소유 복원만 중단하고 마지막 수락 상태·실패 제안을 보존한다. 별도 실제 결합218은 고정 소스·계획으로 계속 실행한다. 이번 배포는223 시작과 실제 추가 단계 적용의 증거이며 그 후의 완료를 주장하지 않는다.

분류: Conjectural. 완성된 동일 이력을 실제 GR 원천·시간 판독에 연결하고 지배적인 시간 오차를 해결한 뒤, 같은 해의 에너지·경계 이력에서 최종 전하를 판정해야 한다. 전 기간, 전체EOS·균일 미분, 공간·경계, 자기GR·비선형, 정적 비교·관측·무한대 전하 조건은 그대로 남는다. 이력 기록의 수와 개별 단계 통과를 최종 연구 성과로 대체하지 않는다.
''',encoding='utf-8')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        f=root/f'docs/{name}.md';prefixes[f.relative_to(root).as_posix()]=dict(bytes=f.stat().st_size,sha256=sha(f))
        with f.open('ab') as h:h.write(('\n\n## 단계221–223 — 원 결합 광자 단계의 회수와 적용\n\n분류: Counterexample candidate. 최종 전하는 미판정이다. 짧은 원 결합 구간8단계를 재생하며 누락된 실제 광자 단계를 저장했고, 모든 원 물리 배열·내부 끝점을 정확히 재현했다. 같은 이력의 반경 출구1.61324e-13은 원1e-12기준을 통과했다. 수정된16단계 체크포인트를 후속 복원에 실제 적용하여 두 경로가 추가 단계를 수락했다. 원 조건부 실패와 출처가 독립적이지 않았던221끝점 대조 실패를 보존한다. 전체 이력·GR 시간 오차·최종 전하 수락은 아직 아니다. [근거](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',actual_same_joint_solution_GR_computed=False,
        original_coupled_capture=result,original_physical_replay=replay,endpoint_control_failure_preserved=True,
        actual_continuation_started=True,continuation_full_history_completed=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(out/'result.json',final);write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes,reused=reused))
    files=modules+[note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes,sha256={p.relative_to(root).as_posix():sha(p) for p in files});write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest);m['native_original_photon_capture']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
