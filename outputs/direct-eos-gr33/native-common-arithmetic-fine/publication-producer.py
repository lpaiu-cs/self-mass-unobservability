"""Freeze the old fine failure and the actual common-arithmetic continuation."""
from pathlib import Path
from types import FunctionType
import datetime,importlib.util,sys
spec=importlib.util.spec_from_file_location('b',Path('.phase238-publish.py'));b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha,copy=b.read,b.write,b.sha,b.copy
out=root/'outputs/direct-eos-gr33/native-common-arithmetic-fine'
manifest=out.parent/'native-common-arithmetic-fine-manifest.json'
note=root/'notes/REQUEST239_COMMON_ARITHMETIC_FINE_KO.md'


def package():
    assert not out.exists();out.mkdir();old=read(b.manifest)
    preserved={p:h for p,h in old['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    previous=runtime/'native-independent-fine232-work';assert read(previous/'controller-status.json')['state']=='failed'
    for p in previous.iterdir():
        if p.is_file() and p.suffix in ['.json','.log','.py']:copy(p,out/'232'/p.name)
    copy(previous/'rejected-joint-stage.npz',out/'232/rejected-joint-stage.npz')
    reused={str(previous/'last-accepted-128.npz'):sha(previous/'last-accepted-128.npz')}
    w=runtime/'native-common-arithmetic239-work'
    for name in ['plan.json','prepare-receipt.json','check-receipt.json','restart-check.json','symbolic.json','controller-start.json']:
        copy(w/name,out/'239'/name)
    snapshots={}
    for number,name in [(239,'native-common-arithmetic239-work'),(238,'native-true-momentum238-work'),(236,'native-complete-return236-work')]:
        base=runtime/name
        for name in ['controller-status.json','stage-progress-128.json','prefix-64.json','prefix-128.json','prefix-result.json','prefix-receipt.json','fine-receipt.json','result.json']:
            p=base/name
            if p.exists():
                r=read(p);write(out/str(number)/('snapshot-'+name),r);snapshots[f'{number}/{name}']=r
        if read(base/'controller-status.json')['state'] in ['completed','failed']:
            for p in base.iterdir():
                if p.is_file() and p.suffix in ['.json','.log','.py']:copy(p,out/str(number)/p.name)
    copy(root/'.phase239-followthrough.py',out/'controller-239.py');copy(Path(__file__),out/'publication-producer.py')
    failure=read(previous/'rejected-joint-stage.json');receipt=read(previous/'fine-receipt.json');prefix=snapshots['238/prefix-64.json']
    now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat();state=snapshots['239/controller-status.json']
    note.write_text(f'''# 동일 산술의 원 미세 경로 속행

분류: Counterexample candidate. **최종 전하는 아직 미판정이다.** 거친 경로119단계는 완료됐으나, 이전 산술로 병행하던232미세 경로는220단계 수락 후221단계 원 비선형 잔차에서 실패했다. 따라서 그 경로를 수정된 거친 해와 단순히 짝지어 수락할 수 없다. 검증된215단계 원 체크포인트에서 거친 경로와 같은 열 좌표·바리온·운동량 수정으로 나머지16단계만 실제 속행하는239를 시작했다.

분류: Counterexample candidate.232의 실제221단계 마지막 잔차는{failure['equations'][-1]['relative']:.12e}>1e-12, 벽시간은{receipt['seconds']:.3f}초로3시간 상한 이내다. 원8회 Newton 제안·선형/물리 기준과 실패 제안을 보존했다. 이 실패만으로 원인의 유일성을 단정하지 않으며239의 실제 풀이 결과로 수정 효과를 판정한다. 기존232의 추가5단계는 다른 산술의 증거로 보존하고 수정 해에 섞지 않는다.

분류: Counterexample candidate.239의 재시작은215단계 체크포인트와 저장 NPZ의 모든 배열을 원소별로 정확히 재현했다. 앞215단계나 완료된 거친119단계를 재적분하지 않는다. 원128시계·531셀·기간·방정식·수락 기준을 유지한다. 거친 완료는 실제238결과를 재사용하고 원 단계의 광자 모멘트와 경계 반환값을 수락 순간에 저장한다. 두 경로에 같은 수정 산술을 적용한 뒤 원10채널·2%시간 대조로 연결한다.

분류: Counterexample candidate.238의 거친 저장prefix236개 Radau 순간은 수정된 native 물질률로 재평가하여 원 수지 기준을 통과했다. 최대 구성률 상대 변화={max(prefix['native_rate_change']):.12e}, 같은 이력 물질 수지 최대={max(prefix['same_prefix_material_balance']):.12e}다. 미세 저장prefix430순간은 별도 검증 중이다. 이 수지 대조는 이전 모든 전체 벡터 방정식의 균일 인증이 아니다.239의 짝 비교는 그 원 수지 검증이 성공한 뒤에만 실행하도록 연결했다.

계산은CPU10·12GiB·6시간 상한으로 실행한다. 마지막 거친 단계의116.38초 실측을 기준으로16개 미세 단계가 약32–64분일 수 있다고 계획했으나 후반 분기 비용은 미측정이며 완료시각 보장이 아니다. 실제 실행 중인236의 독립GR3코어와238의 이력 검증을 변경하지 않았다.232를 기다린 뒤에야 수정 경로를 시작하는 직렬 지연을 피했으며232는 자체 수치 실패로 종료됐다. 격자·시계·정확도나 물리 경로를 자동 확대하지 않았다.

분류: Counterexample candidate. 스냅샷 시각={now},239상태={state['state']},동작={state.get('action')}. 실행 상태는 성공 판정이 아니다. 제한 구간과 공통15/16GR수락·거친 전체 기간 완료를 전체 자기GR·완전 비선형·최종 전하 성과로 확대하지 않는다.

분류: Proven. 거친 경로에서 검증한 원Radau source/collision 묶음 항등식을 재사용한다. 재시작 배열의 정확 일치와 이 항등식은 물리 폐쇄 정리가 아니다.

분류: Conjectural. 동일 산술의 실제 미세 경로와 짝 시간 판정을 통과하면 그 같은 완전 이력의 에너지·경계·광자를 GR 및 최종 전하에 연결한다. EOS/미분·공간·경계·자기GR 고정점·완전 비선형·정적 비교·관측·무한대의 원 미해결 조건은 유지한다. 이번 작업은 loophole progress다.
''',encoding='utf-8')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        p=root/f'docs/{name}.md';prefixes[p.relative_to(root).as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
        with p.open('ab') as f:f.write(('\n\n## 단계239 — 동일 산술의 원 미세 경로\n\n분류: Counterexample candidate. 최종 전하 미판정. 이전 산술의232미세 경로는220단계 수락 후221비선형 잔차에서 실패했다. 원 실패와 수락 상태를 보존한다.239는215단계 모든 저장 배열을 정확히 재시작하고, 거친119단계를 통과한 같은 열 좌표·B/S 산술로 원16개 미세 단계와 짝 시간 비교를 진행한다.238의 거친 저장prefix 물질 수지는 통과했으나 전체 벡터 균일 인증이나 최종 전하 완료가 아니다. [근거](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',actual_same_joint_solution_GR_computed=True,old_fine_failure_preserved=True,
        exact_fine_restart=True,live_snapshots=snapshots,snapshot_KST=now,paired_time_admitted=False,
        final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(out/'result.json',final);write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes,reused=reused))
    module=root/'verification/complete_common_arithmetic_fine.py';assert sha(module)==sha(runtime/'verification'/module.name)
    files=[module,note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes,sha256={p.relative_to(root).as_posix():sha(p) for p in files});write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_common_arithmetic_fine']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
