"""Preserve the bounded restart and actual full-period return continuation."""
from pathlib import Path
from types import FunctionType
import datetime,importlib.util,sys
spec=importlib.util.spec_from_file_location('b',Path('.phase248-publish.py'));b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha,copy=b.read,b.write,b.sha,b.copy
out=root/'outputs/direct-eos-gr33/native-full-return-continuation'
manifest=out.parent/'native-full-return-continuation-manifest.json'
note=root/'notes/REQUEST249_FULL_RETURN_CONTINUATION_KO.md'

def package():
    w=runtime/'native-full-return249-work';r=read(w/'restart-regression.json');anchor=read(w/'anchor/result.json');boundary=read(w/'boundary-extension.json')
    module=root/'verification/complete_returned_period.py';assert r['passed'] and anchor['passed']
    assert sha(module)==sha(runtime/'verification'/module.name)==read(w/'check-receipt.json')['source_sha256']
    for p,h in r['bindings'].items():assert sha(runtime/p.removeprefix('/home/lpaiu/work/native-retained-tail-runtime/'))==h,p
    assert not out.exists();out.mkdir();preserved={p:h for p,h in read(b.manifest)['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    for name in ['restart-regression.json','check-receipt.json','boundary-extension.json','controller-start.json','controller-status.json']:
        copy(w/name,out/name)
    copy(w/'anchor/result.json',out/'native-anchor.json');copy(w/'check/sweep-1/photons/zero-64.json',out/'zero-step-readout.json')
    for name in ['.phase249-native-anchor.py','.phase249-boundary.py','.phase249-followthrough.py','.phase249-launch.ps1']:
        copy(root/name,out/name.lstrip('.'))
    copy(Path(__file__),out/'publication-producer.py');snapshots={}
    for folder,names in [('native-common-arithmetic239-work',['controller-status.json','stage-progress-128.json']),
        ('native-complete-return236-work',['controller-status.json','capture-128.json']),
        ('native-retarded-extension248-work',['controller-status.json']),('native-compensated-charge247-work',['full-controller-status.json'])]:
        for name in names:
            value=read(runtime/folder/name);snapshots[f'{folder}/{name}']=value;write(out/('snapshot-'+folder+'-'+name),value)
    now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    change=next(v for v in boundary['rows'] if v['key']=='actual_delta_lambda_rate')
    note.write_text(f'''# 저장된 실제 결합 해에서 전체 기간으로 이어가기

분류: Counterexample candidate. **최종 물리 전하는 미판정이다.** 단계247의 짧은 같은 반환 해에서 유지된 compact 부호는 그 범위에 한정한다. 이번에는 완성되는248전체 기간 GR을 실제 물질·광자 방정식에 적용하도록 연결했다. 전체를 처음부터 다시 풀지 않고14/16저장 상태에서15/16구간을 겹쳐 재계산한 뒤 마지막16/16구간을 푼다. 게시={now}. 후속 실제 진화는 아직 의존 결과를 기다리며, 재시작 검사만으로 완료를 선언하지 않는다.

분류: Counterexample candidate. 기존 짧은227과 긴236의 같은 시각들을 대조하니 actual_delta_lambda_rate는 종전 마지막 시각에서{change['old_final']:.12e}만큼 달랐다. 척도는 짧은 이력 전체에서의 최대 절댓값이며, 그 마지막 값 자체를 분모로 한 비율이 아니다. 그 이전 시각의 차이는 같은 척도에서{change['prefix']:.12e}였다. 이전 종료점에서 과거쪽 마지막 원천 다항식 미분을 읽던 연산이 연장 후 내부 시각의 오른쪽 다항식 미분을 읽는다. 끝점의 기존 입력을 조용히 유지하면서 전체 기간의 같은 입력이라고 주장하면 안 된다. 이 때문에15/16끝점을 포함하는 마지막 저장 구간 하나를 다시 푼다. 이전 수락 해와 실패를 덮어쓰지 않는다.

분류: Counterexample candidate. coarse14/16저장의103실제 단계를0새 단계로 복원하여44개 모든 NPZ배열이 정확히 일치했다. 상태·광자·보존량·floor·Radau 이력·출구·적분 수지·반복 진단을 재사용한다. 별도의 최종 coarse native anchor 확인에서 추가 구간 첫/마지막 원래 단계의 네 보존 채널 오차 최대는{max(v for row in anchor['rows'] for v in row['native_relative']):.12e}<1e-12였다. 이는 이어 풀 준비가 된 근거이며 실제 추가 단계 수락을 대신하지 않는다.

분류: Conjectural. 실행은236의 실제 짝 반환 수락과248의 전체 기간 원천/GR 수락을 모두 요구한다. 동일 full-period 원천·실제 angular 방출·해석적 원천 미분으로 계량을 만든다.14/16까지의 새/저장 계량 차이가 원천 척도1e-12미만인지 확인하고 그 이전에 실제 적용했던 값을 그대로 유지한다. 접합 후 metric 시간2%·구적0.2%기준을 다시 확인한다. 종전 종료점15/16은 새 내부 미분을 사용하며, 그 값으로 실제15/16과16/16물질·광자 단계를 풀어 같은 해의 출구를 캡처한다. 다른 해의 전하나 진단 수치를 사후 가산하지 않는다.

분류: Conjectural. coarse103/fine199개 실제 단계를 그대로 재사용하고, coarse16/fine32개 단계만 새로 계산한다. 같은 셀·두 시계·전체 선언 기간을 유지한다.236coarse111단계 실측2258초를 기준으로 대략6–12분/12–25분을 예상하되 후반 분기 비용은 미측정이다. CPU3·16GiB,계량1시간·coarse2시간·fine4시간·의존 대기12시간으로 충분한 여유를 둔다. 원 단계/물리 모멘트/native anchor/분기/구성식/보존/출구 및10채널2%짝 대조 기준은 유지한다. 실패 상태와 수락 체크포인트를 보존하고 자동 격자 확대나 수락 기준 완화를 하지 않는다.

분류: Conjectural. 이 실행이 통과해도 같은 전체 기간 반환 해의 연속 원천·전하 판독과 자기GR 고정점·완전 비선형 폐쇄는 별도다. EOS·균일 미분/시간·공간/경계 오차, 정적 비교·관측 및 무한대 전하 정규화도 그대로 남는다. 연구 가치 기준을 유지하며 이번 단계는 loophole progress다.
''',encoding='utf-8')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        p=root/f'docs/{name}.md';prefixes[p.relative_to(root).as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
        with p.open('ab') as f:f.write(('\n\n## 단계249 — 종료점 미분을 반영한 실제 결합 해의 연장\n\n분류: Counterexample candidate. 최종 전하 미판정. 기존 종료점이 내부 시각이 되면 원천 다항식의 미분이 바뀌므로14/16저장 상태에서15/16을 겹쳐 다시 풀고16/16까지 연결한다. coarse103단계·44개 저장 배열의0단계 복원은 정확했고 late native anchor도 원 기준을 통과했다.236/248실제 수락 뒤 같은 전체 입력으로 물질·광자를 이어 풀도록 연결했다. 원 해·실패·모든 수락 기준 및 전체 최종 전하 범위를 유지한다. [근거](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',actual_same_joint_solution_GR_computed=False,
        exact_saved_restart_verified=True,restart_regression=r,native_anchor=anchor,terminal_derivative_change=boundary,
        actual_full_period_return_completed=False,live_snapshots=snapshots,snapshot_KST=now,
        final_charge_conclusion='unadjudicated',physical_final_charge_solved=False,full_goal_complete=False)
    write(out/'result.json',final);write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes,reused=r['bindings']))
    files=[module,note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes,sha256={p.relative_to(root).as_posix():sha(p) for p in files});write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_full_return_continuation']={k:v for k,v in final.items() if k!='sha256'};write(master,m)

check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
