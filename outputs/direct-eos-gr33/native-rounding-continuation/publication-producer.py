"""Freeze the actual117th acceptance and the unchanged remaining failure boundaries."""
from pathlib import Path
from types import FunctionType
import importlib.util,shutil,sys
helper=Path('outputs/direct-eos-gr33/native-radau-pulse-gr/previous-publication-producer.py')
spec=importlib.util.spec_from_file_location('b',helper);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha=b.read,b.write,b.sha
out=root/'outputs/direct-eos-gr33/native-rounding-continuation';manifest=out.parent/'native-rounding-continuation-manifest.json'
note=root/'notes/REQUEST216_219_ACTUAL_ROUNDING_CONTINUATION_KO.md'


def package():
    assert not manifest.exists()
    w216=runtime/'native-full-momentum216-work';w217=runtime/'native-rounding217-work';w218=runtime/'native-rounding-continuation218-work';w219=runtime/'native-port-order219-work'
    assert read(w216/'result.json')['full_arithmetic_change']==0 and not read(w216/'repair-result.json')['linear_passed']
    assert not read(w217/'result.json')['linear_passed'] and read(w217/'refined-result.json')['linear_passed']
    assert not read(w217/'refined-result.json')['actual_stage_passed']
    r=read(out/'first-actual-acceptance/stage-progress-64.json')
    assert r['solved_actual_steps']==117 and r['original_stage_relative']<1e-12 and max(r['physical_stage_relative'])<1e-13
    assert read(w218/'restart-check.json')['passed'] and read(w218/'linear-seed-identity.json')['passed']
    assert read(w219/'result.json')['original_order_relative']>1e-12 and read(w219/'result.json')['timestamp_change_count']==0
    old=read(out.parent/'native-coordinate-recovery-manifest.json');preserved={k:h for k,h in old['sha256'].items() if not k.startswith('docs/')}
    for k,h in preserved.items():assert sha(root/k)==h,k
    for work,names in [
        (w216,['plan.json','prepare-receipt.json','check-receipt.json','result.json','comparison.npz','collision-producer.py','symbolic.json','repair-plan.json','repair-receipt.json','repair-result.json','polish.json']),
        (w217,['plan.json','prepare-receipt.json','check-receipt.json','result.json','search.json','proposal.npz','pair-producer.py','refinement-plan.json','refine-receipt.json','refined-result.json','refined-search.json','refined-proposal.npz']),
        (w218,['plan.json','prepare-receipt.json','check-receipt.json','restart-check.json','symbolic.json','linear-seed-identity.json','controller-start.json','expanded-balanced-linear.py','expanded-postfloor-evolve.py']),
        (w219,['plan.json','receipt.json','result.json','snapshot-input.npz','linear-128.json','clock-128/snapshot-02.json','clock-128/original-equation-15.json'])]:
        for name in names:
            src=work/name;dst=out/work.name/name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);assert sha(dst)==sha(src)
    for name in ['.phase218-followthrough.py','.phase219-port-replay.py']:
        shutil.copyfile(root/name,out/name.removeprefix('.'))
    shutil.copyfile(__file__,out/'publication-producer.py')
    note.write_text('''# 표현 간격 보정을 실제 결합 진화에 적용

분류: Counterexample candidate. **최종 전하의 기존 결론은 미판정이다.** 다만 저장 제안의 보정에 그치지 않고 같은 결합 진화에 적용해 종전에 실패한117번째 실제 단계를 원 기준으로 통과했다. 이번 결과는 loophole progress의 중간 수치 병목 해소이며, 최종 전하·자기GR·전체 연구 완료가 아니다.

분류: Counterexample candidate. 기존210의117선형 실패4.093584129e-10을 보존했다.215의 작은 표현 간격 열 보정은3.301275091e-13에서 멈췄다.216에서는 광자 충돌·탈출·기체 운동량 행을40/80자리로 조립하고 독립 직접 합과 대조했지만 잔차가 바뀌지 않았다. 광자 열을 추가 후보로 허용해도 제한 행에서 쓸 수 있는 방향이 추가되지 않았고 기준은 여전히 실패했다. 따라서 충돌 산술을 이번 잔차의 원인이나 해결책으로 선언하지 않았다.

분류: Counterexample candidate.217은 같은 원 선형 연산자에 대해 큰 표현 간격 변수의 인접 실수와 작은 표현 간격 변수의 최소제곱 보정을 함께 사용했다. 두 변수 검색은68.66초 만에1.668182294e-14까지 줄였지만 원1e-14기준을 실패했다. 세 변수 공동 보정은 추가39.43초에7.980254442e-15, 물리 모멘트 최대1.76470e-19로 원 선형 기준을 통과했다. 원 전체 연산자를 실제 저장 배열에 다시 적용해 판정했으며, 최종 전하·출구값·물리 매개변수에 맞춘 보정이 아니다. 기존 실제 비선형 잔차3.788254250e-7은 실패한 채 보존했고 이 제안을 물리 상태로 수락하지 않았다.

분류: Proven. 선택한 열의 보정은 원 선형식의 잔차 항등식을 보존한다. 기호 검사를 실행했다. 이 항등식은 수렴 보장, 전체 EOS 오차 상계나 비선형 물리 결과의 증명이 아니다.

분류: Counterexample candidate.218에서116단계의 저장 상태·모든 물리 이력을 정확하게 복구했고, 원 우변·Newton 분기 제안·수정 제안의 잔차가 비트 단위로 같은지 검사한 뒤 실제 풀이를 재개했다. 새 보정기를 실제 Krylov/잔차 경로에 넣은 결과117단계가3.380768132ms까지 진전했다. 독립적인 실제 비선형 잔차5.271475431e-13<1e-12, 물리 모멘트 최대1.711251101e-17<1e-13로 통과했다. 첫 수락 기록과 그 뒤 재시작에 쓸 전체 상태를 복사 전후SHA 일치 검사와 함께 보존했다. 이것은 선형 예비 검사만 통과한 기록과 구분되는 실제 결합 단계 수락이다.

실행 예산은 사용자 지시에 따라 coarse2시간·fine3시간, 각CPU1스레드·가상메모리6GiB, Newton8회·선형 보정12회·GMRES restart80/maxiter10을 유지했다. 실제 재개 입력에서coarse3개·fine16개 단계가 남아 있었다. 본 묶음은 첫117수락과 실행 입력을 고정한 기록이며 실행 후반의 완료를 선언하지 않는다. 실행 소스와 계획은 고정했고 수락 기준·공간 격자·물리 기간·경로 수는 바꾸지 않았다. 전체 물리 상태와 실패 제안을 저장하여 다음 재개에서 이미 수락한 앞부분을 반복 적분하지 않는다.

분류: Counterexample candidate. 별도의219에서는214가 저장하지 못한 fine16번째 조건부 광자 블록 하나만 재계산해 이전 반경 출구 실패를 정확히 재현했다.45.95초, 최대RSS2.751GB였고 새 물질 적분은0개다. 원 식·native동일성·광자 끝점·물질수지·각도 출구는 통과하지만 누적 반경 출구4.740709463e-12>1e-12는 실패했다. 원 생산자와 같은 비정규화binary64 Radau 누적 순서로 계산해도4.740513315e-12로 실패하며 저장 경계 시각32개는 모두 기존 복원 시각과 정확히 같았다. 따라서 누적 순서 또는 시각 복원만으로 해결된다는 가설은 이 사례에서 기각했다. 차이는 바깥 경계의 광자 수 출구24.19705608463339 대24.197056084518685에 남았고 에너지 출구는 표시 정밀도에서 같았다. 작은 양자수 출구를 흔드는 광자 상태/복원 오차의 위치는 아직 확정하지 않았으며 원 기준을 완화하지 않는다.

분류: Conjectural. 다음 결정 기준은 남은 실제 결합 단계와 두 시간 경로 비교를 마친 동일 해에서, 아직 실패한 광자 경계 이력과 GR 시간 오차를 해결한 뒤 최종 전하의 기존 결론이 유지되는지다. 첫T/64국소 GR 시간 차이40.86퍼센트를 전 기간 전하 오차로 확대하지 않는다. EOS 자료·균일 미분·공간·경계·완전 비선형·자기GR·정적 비교·관측·무한대 전하의 완료 조건과 기존 실패는 유지한다.
''',encoding='utf-8')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        f=root/f'docs/{name}.md';prefixes[f.relative_to(root).as_posix()]=dict(bytes=f.stat().st_size,sha256=sha(f))
        with f.open('ab') as h:h.write(('\n\n## 단계216–219 — 보정기를 실제 결합 진화에 적용\n\n분류: Counterexample candidate. 최종 전하는 미판정이다. 같은117선형식의 인접 표현 좌표 보정이7.9803e-15로 통과한 뒤, 이를 실제 진화에 적용해117단계를 원 비선형5.2715e-13와 물리 모멘트 기준으로 수락했다.116저장 상태·이력과 선형 우변·분기·잔차 일치를 검사했고 원 실패들은 보존했다. 광자 이력의 누적 반경 출구 실패는 원 누적 순서와 정확히 같은 시각으로도 유지되므로 전체 복원을 수락하지 않는다. 본 기록은 첫 실제 진전과 계속 실행할 입력의 증거이며 전체 기간·자기GR·전하의 완료가 아니다. [근거](../notes/'+note.name+').\n').encode())
    result=dict(classification='Counterexample candidate',actual_same_joint_solution_GR_computed=False,actual_stage117_accepted=True,
        linear_gate=read(w217/'refined-result.json'),actual117=r,continuation_entry_only=True,full_remaining_pair_completed=False,
        photon_recovery_port_failure_preserved=True,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(out/'result.json',result);write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes,reused={}))
    modules=[root/'verification'/n for n in ['repair_full_momentum_collision.py','solve_native_rounding_neighbours.py','continue_rounding_aware_native.py']]
    for p in modules:assert sha(p)==sha(runtime/'verification'/p.name)
    files=modules+[note]+[root/k for k in prefixes]+[f for f in out.rglob('*') if f.is_file()]
    result.update(document_prefixes=prefixes,sha256={f.relative_to(root).as_posix():sha(f) for f in files});write(manifest,result)
    m=read(master);m['sha256'].update(result['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest);m['native_rounding_actual_continuation']={k:v for k,v in result.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
