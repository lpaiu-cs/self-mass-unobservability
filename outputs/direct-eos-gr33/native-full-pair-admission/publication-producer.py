"""Preserve actual completed physical/GR-return pairs and the audit-only repair."""
from pathlib import Path
from types import FunctionType
import datetime,importlib.util,sys
spec=importlib.util.spec_from_file_location('b',Path('.phase249-publish.py'));b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha,copy=b.read,b.write,b.sha,b.copy
out=root/'outputs/direct-eos-gr33/native-full-pair-admission'
manifest=out.parent/'native-full-pair-admission-manifest.json'
note=root/'notes/REQUEST250_COMPLETED_PHYSICAL_PAIRS_KO.md'

def package():
    w=runtime/'native-full-admission250-work';repair=read(w/'result.json')
    primary=runtime/'native-common-arithmetic239-work';returned=runtime/'native-complete-return236-work'
    a=read(primary/'result.json');r=read(returned/'result.json')
    assert repair['passed'] and a['passed'] and r['passed']
    assert read(primary/'controller-status.json')['state']==read(returned/'controller-status.json')['state']=='completed'
    assert repair['completed_physical_histories_unchanged'] and repair['original_producer_unchanged']
    assert not out.exists();out.mkdir();preserved={p:h for p,h in read(b.manifest)['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    for p in w.rglob('*'):
        if p.is_file():copy(p,out/'repair'/p.relative_to(w))
    for folder,label in [(primary,'primary'),(returned,'returned')]:
        for name in ['result.json','controller-status.json','fine-receipt.json','audit-receipt.json']:
            assert read(folder/'fine-receipt.json')['error'] is None and read(folder/'audit-receipt.json')['error'] is None
            copy(folder/name,out/label/name)
    copy(returned/'coarse-receipt.json',out/'returned/coarse-receipt.json')
    copy(primary/'prefix-result.json',out/'primary/prefix-result.json')
    copy(root/'.phase250-admit.py',out/'admission-repair-producer.py');copy(Path(__file__),out/'publication-producer.py')
    inputs=[primary/f'sweep-1/photons/complete-{n}.npz' for n in [64,128]]
    inputs += [returned/f'sweep-1/photons/return-{n}.npz' for n in [64,128]]+[returned/f'recovered-{n}.npz' for n in [64,128]]
    bindings={str(p.relative_to(runtime)):sha(p) for p in inputs}
    for p,h in read(w/'plan.json')['physical_bindings'].items():assert sha(runtime/p)==h,p
    snapshots={}
    for folder in ['native-full-captured244-work','native-returned-source245-work','native-dense-returned246-work','native-compensated-charge247-work','native-retarded-extension248-work','native-full-return249-work']:
        name='full-controller-status.json' if '247' in folder else 'controller-status.json'
        value=read(runtime/folder/name);snapshots[folder]=value;write(out/('snapshot-'+folder+'.json'),value)
        start='full-controller-start.json' if '247' in folder else 'controller-start.json'
        copy(runtime/folder/start,out/(folder+'-'+start))
    now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    note.write_text(f'''# 전체 원 경로와 긴 실제 GR 반환 쌍의 완료

분류: Counterexample candidate. **최종 물리 전하는 여전히 미판정이다.** 실제 원 물질·광자 경로119/231단계가 전체 선언 기간3.434431117929ms에서 완료됐고, 같은 공통15/16이력의 GR을 실제 방정식에 반환한111/215단계도3.219779173058ms에서 완료됐다. 원 경로의 전체 기간과 GR 반환의 공통 기간을 구분한다. 게시={now}. 이번은 실제 두 물리 경로의 수락 진전이며, 원천 검사만으로 전체 완료를 대체하지 않는다.

분류: Counterexample candidate. 원 전체119/231경로의 기존10채널 시간 대조 최대는{max(a['time_relative']):.12e}, 즉{100*max(a['time_relative']):.8f}%로 원2%기준을 통과했다. 높은 정밀도의 B/열/S산술을 사용하는 마지막 원 단계들도 실제 수락됐고 같은 저장 앞부분의 물질 수지를 다시 확인했다. prefix 수지 확인은 모든 과거 전체 벡터 방정식이나 균일 오차 인증은 아니다. fine실행은{read(primary/'fine-receipt.json')['seconds']:.3f}초였으며 이전 서로 다른 산술의 결과를 이어 붙이지 않았다.

분류: Counterexample candidate. 실제 GR 반환111/215쌍의 기존10채널 시간 대조 최대는{max(r['time_relative']):.12e}, 즉{100*max(r['time_relative']):.8f}%로 통과했다. 동일한 원천의 실제 단계 GR과 계량률을 그 해의 물질·광자에 적용했고, 수락 순간의 광자 모멘트·에너지·반경 출구를 함께 저장했다. coarse/fine실행은{read(returned/'coarse-receipt.json')['seconds']:.3f}/{read(returned/'fine-receipt.json')['seconds']:.3f}초다. 한 번의 보상 GR 반환이며 자기GR 고정점이나 완전 비선형 해로 세지 않는다.188초기 국소 실패 Eref약2.42%,H약18.8%와 이전 실패들은 그대로다. 더 긴 구간의 전역 대조 성공으로 국소 실패를 지우지 않는다.

분류: Counterexample candidate.239fine은 실제231단계 완료 후 짝 검사에서 coarse의interval-15-64비교 파일을 찾지 못했다. 물리 단계나 시간 기준이 실패한 것이 아니라 준비 시 파일 연결이 누락돼 판정에 도달하지 못한 오류였다. 이미 수락된238의 같은 prefix NPZ/JSON을 원 SHA와 대조해 연결했다. 원239producer와 수락 기준을 바꾸지 않고 audit만 재실행했다. 기존 실패 로그·receipt·controller·전 결과를 보존했고, 완료된 두 물리 NPZ의 SHA가 전후 같음을 확인했다. 복구와 검사 비용={repair['seconds']:.3f}초, 새 물질/광자 단계=0이다.

분류: Counterexample candidate. 복구 도구의 첫 입력 확인은 같은 파일의 원본/복사본 두 SHA별칭을 한 개로 잘못 가정해 중단됐다. 같은 해시의 모든 별칭을 검증하도록 수정했고, 해당 거절 producer를 보존했다. 과거 파일이나 과학적 기준을 바꾸지 않았다. 상류 audit중단을 받아 아직 아무 작업도 시작하지 않은244/248/249대기 컨트롤러만 실제 종료를 확인한 뒤 실패 상태를 보존하고 다시 시작했다. 원 물리 적분과 진행 중이던245–247은 재시작하지 않았다.

분류: Conjectural. 다음 경로는 두 갈래다.236의 같은 반환 해는245–246연속 원천을 거쳐247전하 판독으로 이어진다.239의 전체 원 해는244원천→248기존 GR장 재사용 연장→249끝점 미분 변경을 반영한 실제 전체 기간 반환으로 이어진다.249는14/16이후만 이어 풀며, 최종 판독에는 그 동일 해의 원천·출구·계량을 사용해야 한다. 현재 실행과 준비 범위는 게시 스냅샷에 보존했다.

분류: Conjectural. 전체 반환 뒤 최종 물리 전하의 유지/변경/제한을 확정하려면 시간·공간·EOS·미분·경계 및 비선형/자기GR 오차가 결론을 지배하는지 평가하고, 같은 재고의 정적 비교·관측량과 무한대 전하 정규화를 연결해야 한다. 짧은247compact부호 결과를 전체 결론으로 확대하지 않는다. 전체 완료 조건과 연구 가치 기준을 유지하며 이번은 loophole progress다.
''',encoding='utf-8')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        p=root/f'docs/{name}.md';prefixes[p.relative_to(root).as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
        with p.open('ab') as f:f.write(('\n\n## 단계250 — 전체 원 경로와 긴 실제 GR 반환 쌍 수락\n\n분류: Counterexample candidate. 최종 전하 미판정. 원 물질·광자의 전체119/231단계는 시간 대조 최대0.00267%, 실제 GR 반환 공통111/215단계는0.01733%로 원2%기준을 통과했다. 원 전체 경로는 이미 완성됐고 누락된 비교 파일만 연결해 동일 audit를 다시 실행했으며 물리 이력 SHA는 그대로다.188초기 국소 실패·균일 오차 및 전체 최종 전하 범위를 유지한다. 같은 반환 해의 전하 판독과 전체 기간의 실제 반환으로 이어간다. [근거](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',actual_same_joint_solution_GR_computed=True,
        full_primary_photon_material_pair_completed=True,primary_time_relative=a['time_relative'],
        common_15of16_actual_GR_return_pair_completed=True,returned_time_relative=r['time_relative'],
        full_period_GR_return_completed=False,immutable_runtime_physical_bindings=bindings,
        actual_history_directory=str(runtime),live_snapshots=snapshots,snapshot_KST=now,
        final_charge_conclusion='unadjudicated',physical_final_charge_solved=False,full_goal_complete=False)
    write(out/'result.json',final);write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes,reused=bindings))
    files=[note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes,sha256={p.relative_to(root).as_posix():sha(p) for p in files});write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_full_pair_admission']={k:v for k,v in final.items() if k!='sha256'};write(master,m)

check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
