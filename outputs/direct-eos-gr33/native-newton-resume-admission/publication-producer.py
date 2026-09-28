"""Bind the actual stage failure and the user's wider continuation budget."""
from pathlib import Path
from types import FunctionType
import importlib.util,shutil,sys
previous=Path(__file__).with_name('.phase186-publish.py')
if not previous.exists():previous=Path(__file__).with_name('previous-publication-producer.py')
spec=importlib.util.spec_from_file_location('previous',previous);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha=b.read,b.write,b.sha
out=root/'outputs/direct-eos-gr33/native-newton-resume-admission';manifest=out.parent/'native-newton-resume-admission-manifest.json'
old=runtime/'native-stable-full194-work';work=runtime/'native-resumed-newton195-work'
modules=[root/'verification'/n for n in ['complete_stable_full_interval.py','resume_full_native_newton.py']]
note=root/'notes/REQUEST194_195_NATIVE_RESIDUAL_CONTINUATION_KO.md'


def package():
    assert not manifest.exists()
    for p in modules:assert sha(p)==sha(runtime/'verification'/p.name)
    fail=read(old/'failure-64.json');receipt=read(old/'coarse-receipt.json');check=read(work/'restart-check.json');plan=read(work/'plan.json')
    assert fail['actual_accepted_steps']==112 and receipt['error'] and check['passed'] and check['new_physical_steps']==0
    prior=read(out.parent/'native-linear-repair-manifest.json');preserved={p:h for p,h in prior['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    known={h:dict(path=p,sha256=h) for p,h in preserved.items()}
    known.update({r['sha256']:r for r in read(out.parent/'native-linear-repair/publication.json')['reused'].values()})
    aliases={sha(old/'prepared-producer.py'):old/'prepared-producer.py'};bindings=[]
    for folder in [old,work]:
        for name,h in read(folder/'plan.json')['bindings'].items():
            p=runtime/name.removeprefix('/home/lpaiu/work/native-retained-tail-runtime/')
            if sha(p)!=h:p=aliases[h]
            assert sha(p)==h,name;bindings.append(dict(source=name,resolved=str(p),sha256=h))
    out.mkdir();reused={}
    sources=[p for p in old.rglob('*') if p.is_file()]
    sources += [work/n for n in ['plan.json','prepare-receipt.json','check-receipt.json','restart-check.json','expanded-resumed-run.py','expanded-eight-newton-stage.py','controller-start.json','sweep-1/photons/resume-check-64.npz','sweep-1/photons/resume-check-64.json']]
    for src in sources:
        folder=old if src.is_relative_to(old) else work;key=folder.name+'/'+src.relative_to(folder).as_posix();digest=sha(src)
        if digest in known:reused[key]=known[digest];continue
        dst=out/key;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);assert sha(dst)==digest
    for src,name in [(Path(__file__),'publication-producer.py'),(previous,'previous-publication-producer.py'),(b.helper,'publication-helper.py'),(root/'.phase195-followthrough.py','followthrough-producer.py')]:shutil.copyfile(src,out/name)
    note.write_text(f'''# 실제 native 단계 실패와 저장 상태의 속행

분류: Counterexample candidate. **최종 전하는 미판정이다.** 선형 상쇄 수정을 실제 결합 구간에 적용한194에서 첫 하위 단계를 수락했지만 다음 단계의 실제 native 잔차가2.15484855e-12로 원1e-12기준을 실패했다. 물리 모멘트는4.71e-17이하였으나 이를 벡터 잔차 실패의 면제로 쓰지 않았다. 선형 풀이 병목은 넘겼고 원 비선형 단계 판정은 아직 통과하지 않았다. fine 경로는 이 실패 뒤 시작하지 않았다.

분류: Counterexample candidate. 194는 최대3회 Newton 제안을 모두 사용했다. 해당 실패는 시간 초과가 아니다. 실제 재개 실행은{receipt['seconds']:.3f}초였고 마지막 수락 상태·단계 이력·floor·충돌·광자 출구와 실패한 Radau 쌍을 저장했다. 그 전 별도 프로세스·WSL 중단에는 원인을 단정하지 않고 보수적180초 비용을 기록했다. 중단 전 새 수락 단계는 저장되지 않았으며 이전15/16지점에서 재시작했다.

사용자 지시: 2026-09-24 계산 리소스를 더 넉넉히 주고 실험을 속행해 전체 비용을 줄이라는 요청을 반영했다. coarse30분, fine45분, CPU1스레드·가상 메모리6GiB로 여유를 확보했다. Newton 제안은 단계별 최대8회로 늘렸다. 이미3회 수행한 실패 단계에는 남은5회만 추가한다. 선형 보정은4회로 유지하고 모든 원 정확도·구성식·보존·출구·2%시간 기준과 기간·시간격자는 유지한다. 자원 한도 변경을 물리 수락 기준 완화로 설명하지 않는다.

분류: Counterexample candidate. 195의 무적분 왕복에서 마지막 수락 하위 단계의 내부 변수·장부·각도 출구·실제 시각·물질 단계를 정확히 복원했다. 새 물리 단계는0개다. 따라서 앞15/16구간과 추가 수락1단계를 재계산하지 않고 coarse7개, fine16개 하위 단계만 남긴다. 실패한 세 번째 Radau 쌍은 다음 Newton의 시작 제안으로만 쓰며 수락 상태로 승격하지 않는다. 이 저장점에서 누락된 첫 새 선형 잔차의 단일 최댓값은 통과 조건1e-14의 상계로 명시하며 실측값을 꾸며내지 않는다.

분류: Proven. 재시작에서 사용한 Radau 가중치의0..2차 모멘트는 정확한 유리수 항등식을 만족한다. 이는 향후 진화의 수렴이나 전하의 증명이 아니다.

분류: Conjectural. 실행 순서는 coarse실제 물리 수락, fine실제 물리 수락, 동일 해의 열 가지 시간 대조다. 실패하면 그 상태를 보존하고 멈추며 더 촘촘한 시간격자나 기간을 자동으로 추가하지 않는다. 이 문서는 실행 수락 계획과 복원 검사를 고정한 기록이며, 이후 생성될 생산 receipt와 판정이 실제 완료 여부를 결정한다. GR 반환 시간 분기 오차·자기 결합·같은 해의 최종 전하·EOS/균일 미분/공간/경계/비선형/관측 폐쇄는 계속 남는다.
''',encoding='utf-8')
    tails={
        'model-definition':'선형 상쇄 수정을 실제 결합 진화에 적용했으나 다음 native 단계는3회 Newton 뒤 원 잔차 기준을 실패했다. 저장된 마지막 수락 상태와 반복 제안을 재사용하는8회 예산으로 속행한다.',
        'observable-targets':'최종 전하는 미판정이다. 실제 단계의 벡터 잔차 실패는 작은 물리 모멘트 오차로 면제하지 않았다. 수락된 같은 해의 전체 이력이 준비된 뒤 전하를 읽는다.',
        'adiabatic-limit':'무적분 재시작의 내부 변수·floor·충돌·실제 출구·단계 시각 복원을 확인했다. 복원 항등식은 동역학이나 균일 미분의 증명이 아니다.',
        'nonadiabatic-regime':'자원 예산을30/45분과8회 Newton으로 넓혔다. 원 시간격자·기간·2퍼센트 시간 대조와 실제 단계 정확도는 유지한다. 기존 실패3회를 포함해 최대8회만 허용한다.',
        'failure-ledger-dynamic-chi':'194의 실제 native 잔차2.15485e-12>1e-12실패를 보존했다. 195는 추가 수락1단계와 실패한3회 제안을 재사용하며, 이전 상태 복원을 새 물리 단계로 세지 않는다.',
        'dynamic-charge-completion':'최종 전하 결론은 미판정이다. 사용자의 총비용 기준으로 속행 자원에 여유를 두었으며 실제 native 단계·전체 기간·반환 시간 오차·자기 GR·동일 해 전하 판독은 결과로 판정한다.'}
    prefixes={}
    for name,body in tails.items():
        p=root/f'docs/{name}.md';prefixes[p.relative_to(root).as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
        with p.open('ab') as f:f.write(('\n\n## 단계194–195 — 실제 단계 실패와 넉넉한 속행 예산\n\n분류: Counterexample candidate. '+body+' [근거](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',admission_and_exact_restart_passed=True,actual_native_stage_rejection=fail,
        new_physical_steps_in_restart_check=0,coarse_remaining_substeps=7,fine_remaining_substeps=16,
        source_budgets_seconds=plan['budgets'],maximum_Newton_proposals=8,accuracy_gates_unchanged=True,
        actual_same_joint_solution_GR_computed=False,full_horizon_completed=False,physical_final_charge_solved=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(out/'final-result.json',final);write(out/'publication.json',dict(reused=reused,previous_nondoc=preserved,document_prefixes=prefixes,plan_bindings=bindings))
    files=modules+[note,root/'AGENTS.md']+[root/k for k in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(sha256={p.relative_to(root).as_posix():sha(p) for p in files},document_prefixes=prefixes);write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_wider_Newton_resume_admission']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
