"""Preserve the tested causal reuse and full-period GR continuation."""
from pathlib import Path
from types import FunctionType
import datetime,importlib.util,sys
spec=importlib.util.spec_from_file_location('b',Path('.phase247-publish.py'));b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha,copy=b.read,b.write,b.sha,b.copy
out=root/'outputs/direct-eos-gr33/native-retarded-extension'
manifest=out.parent/'native-retarded-extension-manifest.json'
note=root/'notes/REQUEST248_CAUSAL_GR_EXTENSION_KO.md'

def package():
    w=runtime/'native-retarded-extension248-work';r=read(w/'regression.json');receipt=read(w/'check-receipt.json')
    assert r['passed'] and receipt['error'] is None and read(w/'symbolic-extension.json')['passed']
    module=root/'verification/extend_retarded_history.py';assert sha(module)==sha(runtime/'verification'/module.name)==receipt['source_sha256']
    for p,h in r['bindings'].items():assert sha(runtime/p)==h,p
    assert not out.exists();out.mkdir();preserved={p:h for p,h in read(b.manifest)['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    for p in w.iterdir():
        if p.is_file() and p.suffix in ['.py','.json']:copy(p,out/p.name)
    for name in ['.phase248-followthrough.py','.phase248-launch.ps1']:copy(root/name,out/name.lstrip('.'))
    copy(Path(__file__),out/'publication-producer.py')
    snapshots={}
    for folder,names in [('native-common-arithmetic239-work',['controller-status.json','stage-progress-128.json']),
        ('native-complete-return236-work',['controller-status.json','capture-128.json']),
        ('native-full-captured244-work',['controller-status.json']),('native-compensated-charge247-work',['full-controller-status.json'])]:
        for name in names:
            value=read(runtime/folder/name);snapshots[f'{folder}/{name}']=value
            write(out/('snapshot-'+folder+'-'+name),value)
    now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    maxerr=max(v for row in r['rows'] for v in row['independent_full_field_reproduction'].values())
    note.write_text(f'''# 같은 인과적 GR 이력을 재사용하는 전체 기간 연장

분류: Counterexample candidate. **최종 물리 전하의 결론은 미판정이다.** 앞서 짧은 동일 반환 해의 compact 전하 부호가 유지된 범위는 단계247에 보존했다. 이번에는 전체 선언 기간의 마지막1/16을 GR에 연결할 때 수락된535개 출력 시각을 다시 계산하는 비용을 제거했다. 완성된244원천이 통과하면 자동으로 실제 전체 기간의 최초 GR 장에 적용한다. 게시={now}. 아직 전체 기간의 GR 반환 해나 최종 전하 완료를 주장하지 않는다.

분류: Proven. 고정된 지연 연산자 G에서 t이전 원천이 같으면 G[S](t)의 적분 구간에 들어가는 원천이 같다. 같은 과거 시계의 선형 potential 표현과 같은 반복 차수를 사용하면 potential 기여의 과거도 같다. 순간 질량 제약은 같은 시각의 원천·장과 같은 공간 연산자로 정해진다. 이 유한 표현의 인과적 재사용 논증은 균일 연속체 오차나 EOS 물리 인증을 포함하지 않는다. 코드에서는 과거 원천·배경·상태 계수·기하 지도·실제 구동을 비트 단위로 비교하고, 두 접합 시각에서 자유장/potential/U/U_t/U_x를 독립 재계산한다.

분류: Counterexample candidate. 기존236의 독립 세 전체 GR 장을 기준으로535시각 중532개를 재사용하고, 나머지3개와 접합2개를 계산했다. 모든 저장 장과 질량·부피·압력 등 후속 배열의 최대 상대차는{maxerr:.12e}였다. 원천 계수 하나를 한 ULP 바꾼 입력도 거절했다. 검사는{receipt['seconds']:.3f}초, 최대RSS={receipt['peak_RSS_bytes']}bytes였으며 새 물질·광자 단계는 없다. 실제 마지막 미완료 구간의 결과는 이 회귀 확인과 구분한다.

분류: Proven. 기존 Radau 보간 항등식과0–9차 특성선 다항식 검사를 다시 실행했고 최대 오차는{max(read(w/'symbolic-extension.json')['box_U_Ut_over_c_Ux_errors']):.12e}였다. 적용 범위는 해당 다항식 연산자의 검사다.

분류: Conjectural. 후속 순서는239의 실제 전체 기간 짝 수락→244의 같은 물질·광자·출구 원천→248의 전체 기간 GR 장이다. 새 시각과 접합2개만 계산하되, 과거 전체를 적분에 포함하는 동일 원천 다항식을 사용한다. 미래 입력을 생략하거나 이전 전하에 무관한 진단을 더하지 않는다. 원 시간2%, 반경 구적0.2%, 독립 적분1e-9, potential 반복1e-8기준을 유지한다. 기존535시각의 직접·질량 응력장과 potential장을 각각 보존하며 실제 source prefix가 한 비트라도 바뀌면 재사용을 중단한다.

분류: Conjectural. 세 독립 작업은CPU1/5/9,각16GiB·2시간으로 실행하고 의존 계산은 최대12시간 기다린다. 회귀에서 다섯 시각 계산의 실측 {', '.join(f"{row['seconds']:.3f}" for row in r['rows'])}초를 사용한다. 실제 추가 시각과 원천 구간 수에 비례한 보수적 추정에2배 시간 여유+5분 및1.5배 메모리 여유를 적용해 예산에 드는지 확인한다. 이 추정에는 고정 setup 비용도 반복 곱해지며 실제 후반 속도는 미측정이다. 추가 물리 격자나 경로를 실행하지 않고, 한 작업 실패 시 다른 작업을 종료한다.

분류: Counterexample candidate. 초기 회귀에서 작업 디렉터리의 sweep-1이 없어 초기화가 실패했고, 원 설정·배경을 hardlink하는 초기화로 수정했다. 독립 회귀 결과와 기존 초기 설정 파일의 이름 충돌도 실행 전에 분리했다. 원 실패 producer/receipt와 이름 변경 전 성공 결과를 함께 보존했으며, 수정본의 세 장 회귀를 다시 통과했다. 물리 방정식이나 수락 기준은 바꾸지 않았다.

분류: Conjectural. 전체 선언 기간의 최초 GR 장 이후에는 그 장을 같은 물질·광자 해에 실제 반환하고, 동일 해의 전하를 판독해야 한다. 진행 중인236–247의 공통15/16반환과 이 전체 기간 최초 장을 혼동하지 않는다. 균일 시간/미분·EOS·공간·경계 오차, 자기GR 고정점·완전 비선형·정적 비교·관측·무한대 전하 요구사항은 유지한다. 이번 작업은 loophole progress이며 전체 목표는 진행 중이다.
''',encoding='utf-8')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        p=root/f'docs/{name}.md';prefixes[p.relative_to(root).as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
        with p.open('ab') as f:f.write(('\n\n## 단계248 — 같은 인과적 GR 이력의 전체 기간 연장\n\n분류: Counterexample candidate. 최종 전하 미판정. 기존535시각 중532개를 재사용하고 마지막3개와 접합2개만 계산해 독립 세 전체 GR 장을 모든 저장 배열에서 정확히 재현했다. 실제 원천 계수 변화는 거절한다. 전체244원천 수락 뒤 원535시각과 같은 과거 원천을 그대로 유지하며 추가 시각만 계산하도록 연결했다. 원 수락 기준과 초기 실패를 보존하고 전체 기간 실제 반환·자기GR·물리 오차·최종 전하 범위는 축소하지 않는다. [근거](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',actual_same_joint_solution_GR_computed=False,
        causal_prefix_regression_passed=True,regression=r,full_period_extension_completed=False,
        live_snapshots=snapshots,snapshot_KST=now,actual_full_period_GR_return_closed=False,
        final_charge_conclusion='unadjudicated',physical_final_charge_solved=False,full_goal_complete=False)
    write(out/'result.json',final);write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes,reused=r['bindings']))
    files=[module,note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes,sha256={p.relative_to(root).as_posix():sha(p) for p in files});write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_retarded_extension']={k:v for k,v in final.items() if k!='sha256'};write(master,m)

check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
