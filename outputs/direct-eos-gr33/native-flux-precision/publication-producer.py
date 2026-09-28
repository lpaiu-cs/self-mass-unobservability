"""Publish terminal actual evolution, preserving both failed hypotheses."""
from pathlib import Path
from types import FunctionType
import importlib.util,shutil,sys
previous=Path(__file__).with_name('.phase186-publish.py')
if not previous.exists():previous=Path(__file__).with_name('previous-publication-producer.py')
spec=importlib.util.spec_from_file_location('previous',previous);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha=b.read,b.write,b.sha
out=root/'outputs/direct-eos-gr33/native-flux-precision';manifest=out.parent/'native-flux-precision-manifest.json'
folders=[runtime/n for n in ['native-postfloor-continuation200-work','native-energy-coordinate201-work','native-flux-precision202-work']]
modules=[root/'verification'/n for n in ['continue_postfloor_native.py','precise_native_baryon.py','continue_native_flux_precision.py']]
note=root/'notes/REQUEST200_202_ACTUAL_PRECISE_NATIVE_FLUX_KO.md'


def package():
    assert not manifest.exists()
    failed,probe,w=folders;control=read(w/'controller-status.json')
    assert control['state'] in ['completed','failed']
    assert read(failed/'controller-status.json')['state']=='failed'
    assert read(failed/'failure-64.json')['actual_accepted_steps']==114
    assert read(probe/'fast-result.json')['exact_match_to_full_bank_conversion']
    assert read(w/'restart-check.json')['passed'] and read(w/'controls.json')['passed']
    gr=None;fields=None
    if control['state']=='completed':
        gr=read(w/'gr-controller-status.json');assert gr['state'] in ['completed','failed']
        folder=runtime/'native-complete-joint-gr-work';folders.append(folder)
        modules.append(root/'verification/read_completed_joint_gr.py')
        if (folder/'fields.json').exists():fields=read(folder/'fields.json')
    for p in modules:assert sha(p)==sha(runtime/'verification'/p.name),p
    old=read(out.parent/'native-primitive-precision-manifest.json')
    preserved={p:h for p,h in old['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    known={h:dict(path=p,sha256=h) for p,h in preserved.items()}
    known.update({r['sha256']:r for r in read(out.parent/'native-primitive-precision/publication.json')['reused'].values()})
    aliases={sha(p):p for folder in folders for p in folder.glob('*.py')};bindings=[]
    for folder in folders:
        for plan in folder.glob('*plan.json'):
            for name,h in read(plan).get('bindings',{}).items():
                p=runtime/name.removeprefix('/home/lpaiu/work/native-retained-tail-runtime/')
                if sha(p)!=h:p=aliases[h]
                assert sha(p)==h,name;bindings.append(dict(plan=str(plan),source=name,resolved=str(p),sha256=h))
    for n in [200,202]:
        path=root/f'.phase{n}-followthrough.py';folder=folders[0 if n==200 else 2]
        assert sha(path)==read(folder/'controller-start.json')['source_sha256']
    out.mkdir();reused={}
    for folder in folders:
        for src in sorted(folder.rglob('*')):
            if not src.is_file():continue
            key=folder.name+'/'+src.relative_to(folder).as_posix();digest=sha(src)
            if digest in known:reused[key]=known[digest];continue
            dst=out/key;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);assert sha(dst)==digest
            known[digest]=dict(path=dst.relative_to(root).as_posix(),sha256=digest)
    helpers=[(Path(__file__),'publication-producer.py'),(previous,'previous-publication-producer.py'),(b.helper,'publication-helper.py')]
    helpers += [(root/f'.phase{n}-followthrough.py',f'followthrough-{n}-producer.py') for n in [200,202]]
    for src,name in helpers:shutil.copyfile(src,out/name)
    result=read(w/'result.json') if (w/'result.json').exists() else None
    failure=next((read(w/f'failure-{n}.json') for n in [128,64] if (w/f'failure-{n}.json').exists()),None)
    counts={}
    for n,initial in [(64,114),(128,215)]:
        path=w/f'path-{n}.json';progress=w/f'stage-progress-{n}.json'
        counts[n]=read(path)['audit']['actual_completed_steps'] if path.exists() else read(progress)['solved_actual_steps'] if progress.exists() else initial
    passed=bool(result and result['passed'] and control['state']=='completed')
    gr_passed=bool(gr and gr['state']=='completed' and fields and fields['passed'])
    event=(f'두 실제 경로의 전체 기간과 원 시간 대조를 통과했다. 수락 단계는 {counts[64]}/{counts[128]}개다.' if passed else
           f'속행은 {control["action"]}에서 종료됐고 두 경로의 수락 단계는 {counts[64]}/{counts[128]}개다. 전체 기간 통과로 판정하지 않는다.')
    gr_event=('같은 전체 해의 에너지·경계 이력을 사용한 GR 원천·장 판독도 통과했다.' if gr_passed else '같은 전체 해의 GR 원천·장 판독은 아직 수락하지 않았다.')
    note.write_text(f'''# 고정밀 질량 유속을 실제 결합 해에 적용

분류: Counterexample candidate. **최종 전하는 미판정이다.** {event} {gr_event} 같은 해의 자기 GR·무한대 전하·물리 EOS/균일 미분/공간/경계/비선형/정적/관측 폐쇄는 유지한다. 다른 해의 진단을 더하거나 이전 전하 부호를 상속하지 않는다. 분류상 loophole progress다.

분류: Counterexample candidate. 200은 실제 post-floor 수락 상태로 초기 Newton 제안을 맞춰 속행했다. 1659.174초 동안 여덟 선형 풀이는 모두 원 기준을 통과했으나 실제 비선형 단계의 잔차는 마지막5.96166e-12로1e-12 기준을 실패했다. 수락 단계는114개로 유지됐다. 초기 제안 변경만으로 문제가 해결된다는 가설을 채택하지 않는다. 시간 상한 중단과 정확도 실패를 구분한다.

분류: Counterexample candidate. 201은 저장된199실패 쌍에서 에너지 좌표 왕복 변환 차이0, 실제 native와 선택 분기 평가의 일치, affine 예측과 실제 단계의 차이1.83663e-11을 확인했다. 같은 원시 변수 복원·재구성·HLL 식을 독립 고정밀 산술로 평가한40/70자리 대조는1.88117e-41차이로 일치했지만 저장 쌍의 실제 바리온 잔차1.83046e-11은 여전히 실패였다. 고정밀 진단으로 이전 상태를 수락하지 않았다. 미사용 광자 계수의 변환을 제외하여 동일 평가값을 유지하면서 두 단계 평가를약38초에서1.28초로 줄였다.

분류: Counterexample candidate. 202는60자리 질량 유속을 실제 Newton 우변의 native-Jg와 독립적인 실제 Radau 잔차에 모두 연결했다. 같은 공유 면 유속으로 raw B와 경계 이력을 계산하고, Eref=Etilde+kappa B가 유지되도록 Etilde 변화율을 함께 수정했다. 다른 물리식·EOS 은행·광자 계수·floor·기간·격자·수락 기준은 유지했다. 이것은 모든 EOS/깊은 셀 연산자를 고정밀화한 것도, 균일 오차 상계를 얻은 것도 아니다.

분류: Counterexample candidate. 200의114단계 수락 상태·전체 이력·수지를 정확히 복원했다. 두 경로의658개 저장 물질 단계를 수정 산술로 재평가하여 구성식 변화 최대5.20802e-12와 동일 이력 물질 수지 최대1.32926e-13으로 원0.002/1e-8기준을 통과했다. 과거 수락 구간을 재적분하지 않았다. 저장 쌍에 대한40/70자리 대조·half/double 구성식 탐침은 통과했지만 그 쌍 자체는 여전히 원 단계 기준을 실패했다.

분류: Counterexample candidate. 실제 수정 방정식으로 기존 거부 제안에서 다시 풀어115번째 단계를 잔차1.270943e-14<1e-12, 물리 모멘트1.231e-17이하로 수락했다. 최종 속행 결과는 {event} 원 실패·실제 단계 이력·두 경로 대조·중단 사유는 생산 결과와 함께 보존한다. 작은 반환 성분의 시간 분기 오차5.211566%는 별도 미해결이며 이 단계 수락으로 해소되지 않는다.

분류: Proven. 질량 변화율의 산술 수정을 Etilde_dot에 -kappa delta B_dot로 반영하면 Eref_dot를 보존한다는 기호 항등식을 검산했다. 이는 물리 EOS 인증·시간 수렴·최종 전하 증명이 아니다.

계산 정책: prefix20분,coarse2시간,fine3시간,CPU1스레드·가상 메모리6GiB,단계당 새 Newton8회·선형 보정12회다. 기존 여덟 실패 제안은 보존하고 고친 산술의 새 풀이를 구분한다. 저장된 유효 상태와 제안을 재사용하며, 실행 중 소스·계획은 바꾸지 않는다. 기준 미달에 따른 자동 격자·기간·경로 확대는 하지 않는다.
''',encoding='utf-8')
    tails={
        'model-definition':'고정밀 native 질량 유속을 실제 Newton 우변·독립 단계 잔차·같은 공유 면과 Etilde 변환에 함께 연결했다. 원 물리식·EOS 은행·구동·격자는 유지한다.',
        'observable-targets':'최종 전하는 미판정이다. '+event+' '+gr_event+' 같은 해의 자기 GR·무한대 전하 및 지배 오차를 함께 판정해야 한다.',
        'adiabatic-limit':'658개 저장 물질 단계의 수정 구성식·동일 이력 수지를 확인했다. 전체 EOS 정밀도·균일 미분 상계 또는 과거 모든 국소 단계의 재인증이 아니다.',
        'nonadiabatic-regime':'이전 실패 실제 단계를 원 잔차1.270943e-14로 수락하여115단계에 도달했다. '+event,
        'failure-ledger-dynamic-chi':'post-floor 초기 제안 변경만으로 실제 비선형 실패가 해소되지 않았다. 에너지 좌표 왕복 가설도 저장 쌍에서 기각했다. 고정밀 평가에서도 옛 쌍은 실패하여 실제 재풀이 전에는 수락하지 않았다. 원 실패와 작은 GR 반환의5.211566퍼센트 시간 실패를 보존한다.',
        'dynamic-charge-completion':'최종 전하는 미판정이다. '+event+' '+gr_event+' 전체 자기 결합·무한대 전하·정적/관측·오차 요건을 축소하지 않는다.'}
    prefixes={}
    for name,body in tails.items():
        p=root/f'docs/{name}.md';prefixes[p.relative_to(root).as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
        with p.open('ab') as f:f.write(('\n\n## 단계200–202 — 고정밀 실제 유속과 결합 진화\n\n분류: Counterexample candidate. '+body+' [근거](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',controller=control,accepted_actual_steps=counts,
        full_horizon_completed=passed,original_failures_preserved=True,high_precision_actual_native_flux=True,
        actual_same_joint_solution_GR_computed=gr_passed,gr_controller=gr,gr_fields=fields,
        physical_final_charge_solved=False,self_GR_return_closed=False,final_charge_conclusion='unadjudicated',
        full_goal_complete=False,actual_result=result,actual_failure=failure)
    write(out/'final-result.json',final);write(out/'publication.json',dict(reused=reused,previous_nondoc=preserved,document_prefixes=prefixes,plan_bindings=bindings))
    files=modules+[note]+[root/k for k in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(sha256={p.relative_to(root).as_posix():sha(p) for p in files},document_prefixes=prefixes);write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_flux_precision']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
