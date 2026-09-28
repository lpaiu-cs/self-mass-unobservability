"""Preserve the bounded completion of the accepted coarse photon history."""
from pathlib import Path
from types import FunctionType
import datetime,importlib.util,sys
spec=importlib.util.spec_from_file_location('b',Path('.phase239-publish.py'));b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha,copy=b.read,b.write,b.sha,b.copy
out=root/'outputs/direct-eos-gr33/native-complete-coarse-photons'
manifest=out.parent/'native-complete-coarse-photons-manifest.json'
note=root/'notes/REQUEST240_243_COMPLETE_COARSE_PHOTONS_KO.md'


def package():
    w=runtime/'native-preimage-photon243-work';state=read(w/'controller-status.json')
    assert state['state'] in ['completed','failed'];assert not out.exists();out.mkdir()
    preserved={p:h for p,h in read(b.manifest)['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    for p in w.rglob('*'):
        if p.is_file() and (p.suffix in ['.json','.log','.py'] or p.name in ['recovered-64.npz','conditional-recovered-64.npz','rejected-original-64.npz','rejected-linear-64.npz'] or p.match('preimage-*.npz')):
            copy(p,out/'243'/p.relative_to(w))
    for number,folder in [(240,'native-coarse-photon240-work'),(241,'native-versioned-photon241-work'),(242,'native-stored-time-photon242-work')]:
        failed=runtime/folder
        for p in failed.rglob('*'):
            if p.is_file() and (p.suffix in ['.json','.log','.py'] or p.name=='rejected-original-64.npz'):
                copy(p,out/str(number)/p.relative_to(failed))
    prefix=runtime/'native-true-momentum238-work'
    for name in ['prefix-64.json','prefix-128.json','prefix-result.json','prefix-receipt.json','controller-status.json']:
        copy(prefix/name,out/'238'/name)
    pfx=read(prefix/'prefix-result.json');assert pfx['passed']
    snapshots={};reused={}
    for number,folder in [(236,'native-complete-return236-work'),(239,'native-common-arithmetic239-work')]:
        base=runtime/folder
        for name in ['controller-status.json','stage-progress-128.json','capture-64.json','capture-128.json','metric-result.json','coarse-receipt.json','fine-receipt.json','result.json']:
            p=base/name
            if p.exists():
                value=read(p);write(out/str(number)/('snapshot-'+name),value);snapshots[f'{number}/{name}']=value
        if number==236:
            for action,n,q in [('field1288',128,8),('field648',64,8),('field1284',128,4)]:
                p=base/f'{action}-receipt.json'
                if not p.exists():continue
                copy(p,out/'236'/p.name)
                p=base/action/f'field-{n}-g{q}-check.json'
                if p.exists():copy(p,out/'236'/p.name)
                p=base/action/f'gr/fields-{n}-g{q}.npz'
                if p.exists():reused[str(p)]=sha(p)
    for name in ['accepted-64.npz']:
        p=w/name
        if p.exists():reused[str(p)]=sha(p)
    for number in [240,241,242,243]:copy(root/f'.phase{number}-followthrough.py',out/f'controller-{number}.py')
    for name in ['.phase242-defect.py','.phase242-times.py']:
        copy(root/name,out/'242'/name.lstrip('.'))
    p=runtime/'native-stored-time-photon242-work/actual-defect.npz';reused[str(p)]=sha(p)
    copy(Path(__file__),out/'publication-producer.py')
    result=read(w/'result.json') if (w/'result.json').exists() else {}
    completed=state['state']=='completed' and result.get('passed',False)
    now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    outcome='거친 원119단계의 광자·물질·출구 이력을 모두 연결했다.' if completed else '원 기준 미달로 복구가 중단됐으며 실패 자료를 보존했다.'
    metric=snapshots.get('236/metric-result.json',{});fine=snapshots.get('239/stage-progress-128.json',{})
    note.write_text(f'''# 거친 전체 기간의 실제 광자 이력 연결

분류: Counterexample candidate. **최종 전하는 아직 미판정이다.** {outcome} 시각={now}, 실행 상태={state['state']}. 이 작업은 이미 수락된 물질 해의 GR 원천을 완성하기 위한 연결이며, 추가 물질 진화나 최종 전하 판정을 대신하지 않는다.

분류: Counterexample candidate. 기존223의111단계와235/238의 실제118/119단계 광자 캡처를 재사용했다. 빠진112–117의6단계에만 원래 광자 Radau 블록을 저장된 물질 단계에 조건부로 풀었다.114부터는 당시199의 확장 정밀도 물질 복원,115–117에서는 당시202의60자리 B 산술을 우변·실제 잔차·에너지 좌표 공변환에 함께 적용한다. 모든 새 쌍은 원 전체 벡터/물리 잔차, native 저장률의 정확 일치, 각도 출구 기준을 통과해야 한다. 물질 재적분이나 저장률을 맞추는 보정은 없다.

분류: Counterexample candidate.240은112/113복구를 통과했으나114에서199물질 복원 변경을 누락하여 native 저장률 정확 일치와 원 벡터 기준에 실패했다. 실패 자료는 보존한다.241은 실제199함수를 다시 사용하며 수락113끝점과 저장된114광자 제안을 재사용한다. 먼저 같은 저장 쌍의 정확 native 일치·전체 방정식을 검사하고 실제 남은 복구와 전체 수지로 연결한다. 성공한 진단만으로 광자 이력 완료를 선언하지 않는다.

분류: Counterexample candidate.241은 실제114/115쌍의 원 전체 방정식과 native 저장률 정확 일치를 통과했다.116의 광자 쌍을 저장한 뒤 재계산한t+c*h와 실제 저장 시각의1ULP차이로 고정밀B잔차 캐시 조회가 실패했다.242는 원 저장 시각을B우변과 잔차에도 동일하게 사용하고115수락 끝점·116저장 제안을 재사용한다. 원 시간 일치1e-18·모든 잔차/수지 기준과 두 실패를 보존한다. 시간 가중치·실제 입력·물리 방정식을 바꾸지 않았다.

분류: Counterexample candidate.242는 native 저장률이 정확히 같아졌지만116의 전체 벡터3.79152e-11로 실패했다. 지배 잔차는 둘째Radau상태의셀261 B에서0.0312499996이었다.243은 독립적으로 재계산한 원202 B율과 원Radau식에 더 가까운 정규화B의 인접 표현값을 최대8개까지 비교하되, 모든 저장 보존량이 비트 단위로 동일한 경우만 허용한다. 물리 보존 상태를 수정하지 않으며, 실제 광자 우변·원 native 정확 일치·전체 벡터/물리·끝점·수지를 다시 판정한다. 복원된 모든 내부 좌표가 원 실행의 비트와 같다는 별도 주장은 하지 않는다. 실패와 수정 전후 좌표를 보존한다.

분류: Counterexample candidate. 복구 끝점117은231의 실제 수락 체크포인트와 비교한다. 초기 준비에서235의 체크포인트를117로 잘못 지정했으나 그 파일은118상태여서 단계 수 검사가 실행 전에 차단했다. 실패 producer/receipt를 보존하고 실제117체크포인트에 SHA를 결속한 뒤 진행했다. 원 기준을 낮추지 않았다.

분류: Counterexample candidate. 최종 조립 결과={result}. 조립 시 앞111단계의 시간·가중치·광자·충돌·출구 배열을 정확히 보존하고, 실제 캡처4개 Radau 순간의 시간·가중치·충돌률을 원119단계 해와 정확히 대조한다. 마지막 광자 끝점은 같은 실제 결합 체크포인트에서 가져오며, 전체 기간의 물질 수지와 반경·각도 출구도 원 기준으로 판정한다.

분류: Counterexample candidate.238의 저장prefix 검증은 두 경로 모두 완료됐다. 거친236순간·미세430순간의 물질 수지 최대는 각각{max(pfx['rows'][0]['same_prefix_material_balance']):.12e}, {max(pfx['rows'][1]['same_prefix_material_balance']):.12e}다. 이는 저장 수지의 재검증이며 모든 과거 벡터 방정식의 균일 인증은 아니다.239원 미세 경로와 짝2%판정,236공통15/16실제GR반환은 별도로 계속된다. 스냅샷은 게시 시점 자료이며 실행 상태를 성공으로 세지 않는다.

분류: Counterexample candidate. 게시 시점239수락 단계={fine.get('solved_actual_steps')}, 실제 마지막 잔차={fine.get('original_stage_relative')}다. 수정 산술에서 이전232실패221단계를 실제로 통과했지만 미세 전체 기간이나 짝 시간 수락으로 확대하지 않는다.236의535실제시각 GR장3경로와 반환 계량 기준 통과={metric.get('passed')}. 최대 계량 시간 차이={max(metric.get('controls',{}).get('time',{'unavailable':float('nan')}).values())}, 최대 구적 차이={max(metric.get('controls',{}).get('quadrature',{'unavailable':float('nan')}).values())}이며 각각 원0.02/0.002기준으로 판정했다. 이 입력이 실제 결합 반환 적분으로 전달됐고 그 same-solution결과가 다음 판정 대상이다. 공통15/16한 번의 반환이며 전체 기간·자기GR고정점·무한대 전하 완료가 아니다.

복구는 비어 있는CPU2에8GiB·1시간을 배정했다.240은95.21초,241은검사 포함112.51초였다. 각 실패에서 수락된 광자 단계와 실패 쌍을 저장하여 재풀이하지 않았다.243에서도 새 조건부 풀이가 필요한 것은117단계 하나뿐이었다. 기존236/239의 소스·계획·프로세스를 변경하지 않았다. 새 격자·기간·물리 경로를 추가하지 않고 전체 기간 연결에 필요한 누락분만 계산한다.

분류: Proven. 원Radau0–2차 모멘트 항등식을 유리수로 검사했다. 이 대수 항등식은 물리 EOS·미분·공간·경계 오차 보장이 아니다.

분류: Conjectural. 전체 미세 경로와 원 시간 대조가 통과한 뒤 이 같은 해의 광자·물질·경계를 전체 기간 GR 원천으로 사용한다. 실제 반환 해에서 전하 결론이 유지되는지가 연구 가치 기준이며, 자기GR·완전 비선형·EOS/미분·공간/경계·정적 비교·관측·무한대 정규화는 여전히 완료 조건이다. 이번 작업은 loophole progress다.
''',encoding='utf-8')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        p=root/f'docs/{name}.md';prefixes[p.relative_to(root).as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
        with p.open('ab') as f:f.write(('\n\n## 단계240–243 — 같은 전체 기간의 누락 광자 연결\n\n분류: Counterexample candidate. 최종 전하 미판정. '+outcome+' 기존111단계와 실제118/119광자 캡처를 재사용하고 누락6단계만 원 방정식·정확 native 이력·끝점·수지·출구 기준으로 판정한다.240의물질 복원 버전,241의시각 캐시,242의정규화B잔차 실패를 보존한다. 실제199복원·원 저장 시각·모든 보존량이 정확히 같은B좌표 복원을 적용하고 수락 상태·실패 광자 쌍을 재사용한다.238의 두 저장prefix 수지는 통과했으나 균일 벡터 인증은 아니다. 원 미세 경로/짝 시간 판정과 실제GR반환, 전체 전하 완료 조건은 유지한다. [근거](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',actual_same_joint_solution_GR_computed=True,
        complete_coarse_photon_history=completed,coarse_history_result=result,prefix_material_ledger_passed=True,
        live_snapshots=snapshots,snapshot_KST=now,paired_time_admitted=False,
        final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(out/'result.json',final);write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes,reused=reused))
    modules=[root/f'verification/{name}.py' for name in ['complete_coarse_photon_history','resume_versioned_photon_history','finish_stored_time_photons','finish_exact_preimage_photons']]
    for module in modules:assert sha(module)==sha(runtime/'verification'/module.name)
    files=modules+[note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes,sha256={p.relative_to(root).as_posix():sha(p) for p in files});write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_complete_coarse_photons']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
