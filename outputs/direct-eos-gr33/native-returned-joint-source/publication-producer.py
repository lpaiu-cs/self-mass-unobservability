"""Preserve actual returned-source application and the gated long consumer."""
from pathlib import Path
from types import FunctionType
import datetime,importlib.util,sys
spec=importlib.util.spec_from_file_location('b',Path('.phase244-publish.py'));b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha,copy=b.read,b.write,b.sha,b.copy
out=root/'outputs/direct-eos-gr33/native-returned-joint-source'
manifest=out.parent/'native-returned-joint-source-manifest.json'
note=root/'notes/REQUEST245_RETURNED_JOINT_SOURCE_KO.md'


def package():
    w=runtime/'native-returned-source245-work';r=read(w/'check/result.json');assert r['representation_controls_passed'] and r['source_time_passed']
    assert not out.exists();out.mkdir();preserved={p:h for p,h in read(b.manifest)['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    for p in (w/'check').iterdir():
        if p.is_file() and p.suffix in ['.json','.py']:copy(p,out/'check'/p.name)
    for n in [64,128]:copy(w/f'check/gr/endpoint-{n}.npz',out/f'check/gr/endpoint-{n}.npz')
    copy(w/'controller-start.json',out/'controller-start.json')
    for name in ['.phase245-followthrough.py','.phase245-launch.ps1','.phase245-source-scale.py']:
        copy(root/name,out/name.lstrip('.'))
    copy(Path(__file__),out/'publication-producer.py');snapshots={}
    for number,folder in [(236,'native-complete-return236-work'),(239,'native-common-arithmetic239-work'),(244,'native-full-captured244-work'),(245,'native-returned-source245-work')]:
        for name in ['controller-status.json','stage-progress-128.json','capture-64.json','capture-128.json','coarse-receipt.json','run-64.json']:
            p=runtime/folder/name
            if p.exists():
                value=read(p);snapshots[f'{number}/{name}']=value;write(out/str(number)/('snapshot-'+name),value)
    now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat();cost=read(w/'check/source-receipt.json')
    scale=read(w/'check/component-scale.json');rows=r['rows']
    note.write_text(f'''# 같은 실제 GR 반환 해의 원천 판독

분류: Counterexample candidate. **최종 전하 결론은 아직 미판정이다.** 짧은 실제227반환 해의 물질·광자·자체 출구와, 그 해를 구동한 반환 기하를 전하 원천 판독에 연결했다. 진행 중인 긴236의 실제 수락 이후 같은 판독을 자동 적용한다. 게시 시각={now}. 짧은 확인을 긴 해나 최종 전하의 성공으로 세지 않는다.

분류: Counterexample candidate.224소비자는 최초 구동 기하를 읽고212의 dense 표현은 한 가지 구동 모드를 전제로 한다.245는 실제227 StageDriver와 그 해에 적용한 동일한 fine metric을 사용한다. 모든 실제 끝점에서 u·lambda·실제 lambda 시간률이 저장 입력과 정확히 일치하는지 확인했다. 원 안정 에너지 변환·보존량 역변환·native 압력·같은 광자 모멘트·floor·Radau 출구 이력을 재사용하며 새 물질·광자·GR 적분은 없다.

분류: Counterexample candidate. 첫 적용은 원64/128시계의15/29단계,0.429303889741ms에 한정된다. 원천 시간 차이 최대={max(r['source_time'].values()):.12e}<0.02, 압력 probe 최대={max(v['pressure_probe'] for v in rows):.12e}<0.002, 압력 mapping 최대={max(v['pressure_mapping'] for v in rows):.12e}<1e-12, 물질 수지 최대={max(x for v in rows for row in v['local_material_ledger'] for x in row):.12e}<1e-8로 통과했다. 실행={cost['seconds']:.3f}초,peakRSS={cost['peak_RSS_bytes']}bytes다. 실패한188초기 국소 대조와 다른 원 실패는 유지한다.

분류: Counterexample candidate. 같은 실제 보상 풀이의 낮은 성분을 높은 성분과 분리해 읽었다. 짧은 fine 경로의 시간별 공간L1 최대 비는 바리온8.90319e-27, 비정지 에너지/trace약1.5666e-24, 광자 에너지3.93990e-23이다. 이는 실제 원천 성분의 크기이며 전하 비율·균일 오차·고정점 수축 증명이 아니다. 단순 합산은 이 작은 성분을 반올림으로 잃을 수 있으므로 이후 같은 선형 readout도 성분을 유지해야 한다. 다른 해의 진단 전하를 사후 가산하는 방식은 허용하지 않는다.

분류: Proven. 고정된 선형 원천 연산자의 높은/낮은 성분 항등식을 확인했다. 물리 EOS·분기 전체·완전 비선형 또는 자기GR 수축 정리는 아니다.

분류: Conjectural. 긴236의111/215실제 단계가 원 짝 대조까지 통과하면 이력과 applied metric에 SHA를 결속해 같은 끝점 판독을 수행한다. CPU6·16GiB·45분으로 배정했다.224의326단계 끝점185.92초와 이번44단계56.86초가 실측 근거이며, 전체 readout은 대략4–8분으로 예상하지만 후반 비용은 미측정이다. 의존 계산에는 기존 거친3/미세6시간과 검사 시간을 감안해 최대9시간을 기다린다. live PID·시작시각·부팅ID를 확인하고 실패 시 중단한다. 기존236/239/244의 소스·계획·실행은 바꾸지 않았다.

분류: Conjectural. 이 결과는 같은 해의 끝점 원천이다. 실제 반환 기하의 연속 시간 표현과 같은 Radau 상태의 retarded 전하, 전체 기간, 자기GR 고정점, EOS/미분·공간/경계·완전 비선형·정적 비교·관측·무한대 정규화는 남아 있다. 최초 구동 모드의 dense 다항식을 반환 기하에 그대로 적용하지 않는다. 연구 가치 기준과 완료 범위는 그대로이며 이번 작업은 loophole progress다.
''',encoding='utf-8')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        p=root/f'docs/{name}.md';prefixes[p.relative_to(root).as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
        with p.open('ab') as f:f.write(('\n\n## 단계245 — 실제 반환 기하와 같은 해의 원천 판독\n\n분류: Counterexample candidate. 최종 전하 미판정. 짧은 실제227반환 해의15/29단계 끝점을 그 해의 반환 기하·물질·광자·출구와 함께 판독했고, 원천 시간 최대0.635%와 압력/수지가 원 기준을 통과했다. 진행 중인 긴236이 실제 수락된 뒤 동일 판독을 적용한다. 높은/낮은 성분을 유지하며 최초 구동의 한 모드 표현을 반환 기하로 대체하지 않는다. 끝점 검사는 연속 시간·최종 전하·물리 폐쇄의 완료가 아니다. 원 실패와 전체 완료 범위를 유지한다. [근거](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',actual_same_joint_solution_GR_computed=False,
        actual_returned_solution_endpoint_source_applied=True,short_source_result=r,
        live_snapshots=snapshots,snapshot_KST=now,full_return_readout_completed=False,
        dense_source_completed=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(out/'result.json',final);write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes,
        reused=scale['bindings']))
    module=root/'verification/read_returned_joint_source.py';assert sha(module)==sha(runtime/'verification'/module.name)
    files=[module,note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes,sha256={p.relative_to(root).as_posix():sha(p) for p in files});write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_returned_joint_source']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
