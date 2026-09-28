"""Preserve a measured cost repair applied to the unchanged coupled equations."""
from pathlib import Path
from types import FunctionType
import datetime,importlib.util,json,sys
spec=importlib.util.spec_from_file_location('b',Path('.phase258-publish.py'));b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,sha,copy=b.read,b.sha,b.copy
out=root/'outputs/direct-eos-gr33/native-short-krylov'
manifest=out.parent/'native-short-krylov-manifest.json'
note=root/'notes/REQUEST259_ACTUAL_KRYLOV_COST_KO.md'


def write(p,value):
    Path(p).write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n',encoding='utf-8',newline='\n')


def package():
    test=runtime/'native-short-krylov259-work';actual=runtime/'native-short-return259-work'
    result=read(test/'result.json');restart=read(actual/'restart-regression.json');takeover=read(actual/'supersession.json')
    assert result['passed'] and result['speedup']>2 and restart['passed'] and restart['every_saved_array_exact']
    assert read(actual/'symbolic-and-adapter-check.json')['passed']
    assert takeover['new_canonical_steps']>=takeover['old_actual_steps']
    modules=[root/'verification'/n for n in ['check_returned_krylov_cost.py','continue_returned_short_krylov.py']]
    assert all(sha(p)==sha(runtime/'verification'/p.name) for p in modules)
    for p,h in read(actual/'plan.json')['bindings'].items():assert sha(runtime/p)==h,p
    old=runtime/'native-left-metric258-work';old_status=read(old/'pipeline-status.json')
    assert old_status['state']=='failed' and old_status['completed'][0]['returncode']==-15
    preserved={p:h for p,h in read(b.manifest)['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    assert not out.exists();out.mkdir()
    for p in test.rglob('*'):
        if p.is_file() and p.suffix in ['.json','.npz','.py'] and not any(k in p.relative_to(test).parts for k in ['sweep-0','sweep-1','metric','gr','__pycache__']):copy(p,out/'saved-system'/p.relative_to(test))
    for name in ['plan.json','prepare-receipt.json','check-receipt.json','restart-regression.json','symbolic-and-adapter-check.json','prefix-seed.json','controller-start.json','controller-status.json','pipeline-status.json','capture-64.json','supersession-intent.json','supersession.json','charge-reader.py','exterior-reader.py']:
        copy(actual/name,out/'actual-start'/name)
    for name in ['controller-start.json','controller-status.json','pipeline-status.json','capture-64.json','coarse.stderr.log']:
        copy(old/name,out/'superseded258'/name)
    rows=[p for p in (actual/'sweep-1/photons').glob('interval-*-64.json') if read(p)['actual_completed_steps']==takeover['new_canonical_steps']]
    assert len(rows)==1;accepted=rows[0];row=read(accepted);assert row['passed']
    copy(accepted,out/'actual-start'/accepted.name);copy(accepted.with_suffix('.npz'),out/'actual-start'/accepted.with_suffix('.npz').name)
    copy(actual/'prefix-recovered-64.npz',out/'actual-start/prefix-recovered-64.npz')
    copy(root/'.phase259-followthrough.py',out/'followthrough.py');copy(Path(__file__),out/'publication-producer.py')
    for p in root.glob('.phase259*.log'):copy(p,out/'logs'/p.name.lstrip('.'))
    now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat();count=read(actual/'capture-64.json')['actual_steps']
    note.write_text(f'''# 기존 수락 기준을 유지한 실제 결합 계산 속행

분류: Counterexample candidate. **최종 물리 전하의 결론은 아직 미판정이다.** 단계258의 압력일 구간 끝 미분 수정은 유지하고, 비싼 내부 선형 반복을 줄인 풀이를 실제 광자·물질 결합 진화에 적용했다. 저장된 coarse 8단계를 재적분 없이 정확히 재사용했으며 새 방식으로 {takeover['new_canonical_steps']}단계까지 수락한 뒤 기존 진행량 {takeover['old_actual_steps']}단계를 따라잡아 기존 작업을 종료했다. 이 문서는 {now} 스냅샷이며 당시 새 수락 단계는 {count}/119다. 전체 두 경로와 같은 해의 질량·전하 판독은 실행 중이다.

분류: Counterexample candidate. 기존 오른쪽 전처리의 내부 GMRES는 작은 내부 잔차를 얻고도 이중 정밀도 재시작을 오래 반복했다. 같은 저장 단계의 실제 초기 상태와 수락 물질 방향에서 선형화한 동일 행렬·우변·초기 추정값을 비교했다. 원래 첫 Newton 행렬을 재현했다는 주장은 하지 않는다. 저장된 단계의 참 비선형 잔차를 먼저 검산했으며 벡터 3.4790e−16, 물질 최댓값 1.4723e−14로 기존 기준을 통과했다.

분류: Counterexample candidate. 내부 한 번의 재시작 후 long-double 참 잔차를 갱신하면 63회/2번의 선형 보정으로 {result['rows'][0]['seconds']:.3f}초가 걸렸다. 기존 최대20회 재시작 방식은 1172회/1번 보정으로 {result['rows'][1]['seconds']:.3f}초였다. 이 저장 계의 비용 비율은 {result['speedup']:.3f}배이며 두 해의 상대 차이는 {result['solution_relative_difference']:.12e}였다. 두 방식 모두 원 선형 벡터1e−14·물리량1e−13·개별 물질1e−13 기준을 통과했다. 이 비율을 전체 실제 진화의 속도나 물리 정확도 증명으로 확대하지 않는다.

분류: Counterexample candidate. 비교 도구의 최초 재구성은 실제 단계와 달리 collision(source=True)를 빠뜨려 거절됐다. 원 실패 소스·영수증·수치를 보존하고 도구를 원 방정식과 일치시켰다. 중복된 절대 소스 해시 경로를 갱신하지 않아 초기화 전에 종료된 실행도 보존했다. 실제258생산자와 그 수락 기준은 변경하지 않았다. 유효 비교는 {read(test/'check-receipt.json')['seconds']:.2f}초였으며 격자나 추가 물리 경로를 만들지 않았다.

분류: Proven. 오른쪽 보정의 잔차 항등식 r_new=(b−A x)−A P y는 내부 반복 횟수와 독립이다. 이를 기호 검산했다. 따라서 내부 한도는 제안값의 비용을 바꾸지만 수락된 해의 정확도는 별도의 참 잔차 검사로 판정해야 한다. 이 항등식은 수렴·조건수·전역 오차 보장이 아니다.

분류: Counterexample candidate. 기존 canonical 체크포인트의 모든 저장 배열은 0단계 재시작에서 정확히 일치했다. 저장 native 수송률도 오차 0으로 재현했다. 상태·물질 및 광자 이력·보존 장부·floor·가이드·기존 잔차 기록을 그대로 복원했고 원259이전 입력·계량·119/231단계 시계를 유지했다. 새 실제 canonical 구간도 원 단계·물리·개별 물질·native·구성·보존 기준을 통과했다. 기존258진행보다 뒤처진 상태에서 기존 작업을 종료하지 않았다. 원258프로세스의 -15종료는 PID·시작 tick·boot·명령을 확인한 의도적 비용 교체이며 수치 수락 실패와 구별한다.

분류: Conjectural. 실제 풀이는 최대12번의 짧은 내부 제안 후에도 기존 참 잔차 기준을 충족하지 못하면 마지막 제안에서 기존20회 재시작 풀이를 이어간다. 원3회 Newton·모든 물리/수치 수락 기준은 유지한다. 원 coarse/fine 및 해당 해의 연속 원천·compact 전하·고정 외부·질량 독립 감사를 자동으로 연결했다. CPU4·16GiB, coarse6시간/fine8시간의 여유를 유지하며 신규 격자·기간·시계를 추가하지 않는다. 저장 선형계의 속도는 실측됐으나 후반 전체 단계 속도는 미측정이므로 계획상 남은 coarse30~120분/fine60~240분 및 판독50~90분은 가정 범위다.

분류: Conjectural. 최종 가치는 지배 오차를 수정한 동일 결합 해의 질량·전하 기준으로 판정한다. 원257질량 시간2.005740886%실패와2%기준, 단계258수정으로 해결하지 않은 나머지 coarse 구적 오차를 유지한다. 다른 해의 진단 에너지를 전하에 가산하지 않는다. 물리 외부 전파·접합·같은 해의 계량 재적용, 자기GR, EOS/균일 오차, 완전 비선형·정적 비교·관측 연결은 여전히 열린다. 이번은 실제 진화를 속행하는 loophole progress이며 전체 연구 완료가 아니다.
''',encoding='utf-8',newline='\n')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        p=root/f'docs/{name}.md';prefixes[p.relative_to(root).as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
        with p.open('ab') as f:f.write(('\n\n## 단계259 — 같은 결합 해의 반복 비용 수정과 속행\n\n분류: Counterexample candidate. 최종 전하의 결론은 미판정이다. 원258미분 수정과 수락 기준을 유지한 동일 저장 선형계에서 잦은 참 잔차 갱신이 비용을15.677배 줄였다. 그 방식을 실제 결합 진화에 적용했고, 저장8단계의 모든 배열·native률을 정확히 복원한 뒤 새 canonical 구간이 기존 진행량을 따라잡았음을 확인해 기존 프로세스를 교체했다. 비용 비교를 전체 물리 정확도나 최종 성과로 확대하지 않는다. 같은 새 해의 질량·전하 자동 판독과 전체 물리 폐쇄의 미완료 범위를 유지한다. [근거](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',passed=False,verdict='ACTUAL_UNCHANGED_EQUATION_COST_REPAIR_PENDING_FINAL_CHARGE',
        saved_system_passed=True,saved_system_speedup=result['speedup'],exact_restart_passed=True,
        actual_new_canonical_passed=True,actual_same_joint_solution_GR_computed=True,
        current_coarse_saved_steps=count,new_actual_return_completed=False,new_final_charge_read=False,
        original257_mass_time_failure_preserved=True,original_scientific_gates_preserved=True,
        physical_final_charge_solved=False,final_charge_conclusion='unadjudicated',full_goal_complete=False,snapshot_KST=now)
    write(out/'result.json',final);write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes))
    files=modules+[note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes,sha256={p.relative_to(root).as_posix():sha(p) for p in files});write(manifest,final)
    masterdata=read(master);masterdata['sha256'].update(final['sha256']);masterdata['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    masterdata['native_short_krylov']={k:v for k,v in final.items() if k!='sha256'};write(master,masterdata)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
