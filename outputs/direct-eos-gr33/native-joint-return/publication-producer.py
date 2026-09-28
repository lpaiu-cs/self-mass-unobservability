"""Bind the same-solution metric and actual compensated stage application."""
from pathlib import Path
from types import FunctionType
import importlib.util,shutil,sys
previous=Path(__file__).with_name('.phase186-publish.py')
if not previous.exists():previous=Path(__file__).with_name('previous-publication-producer.py')
spec=importlib.util.spec_from_file_location('previous',previous);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha=b.read,b.write,b.sha
work=runtime/'native-joint-return187-work';out=root/'outputs/direct-eos-gr33/native-joint-return'
manifest=out.parent/'native-joint-return-manifest.json';module=root/'verification/return_joint_gr_geometry.py'
note=root/'notes/REQUEST187_SAME_SOLUTION_GR_RETURN_KO.md'


def package():
    assert not manifest.exists() and sha(module)==sha(runtime/'verification'/module.name)
    metric=read(work/'metric-result.json');applied=read(work/'compensated-result.json')
    assert metric['passed'] and applied['passed'] and not read(work/'application-admission.json')['passed']
    assert read(work/'symbolic.json')['passed'] and read(work/'odd-source-symbolic.json')['passed']
    assert read(work/'compensated_current-receipt.json')['error'] is None
    total=sum(read(p)['seconds'] for p in work.glob('compensated*-receipt.json'));assert total<90
    old=read(out.parent/'native-joint-gr-manifest.json');preserved={p:h for p,h in old['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    aliases={sha(p):p for p in work.glob('*producer*.py')};aliases[sha(module)]=module;bindings=[]
    for plan in work.glob('*plan*.json'):
        for name,h in read(plan).get('bindings',{}).items():
            actual=runtime/name.removeprefix('/home/lpaiu/work/native-retained-tail-runtime/')
            if Path(name).name==module.name and sha(actual)!=h:actual=aliases[h]
            assert sha(actual)==h,(plan,name)
            bindings.append(dict(plan=plan.name,source=name,resolved=str(actual),sha256=h))
    out.mkdir();reused={};inputs=read(work/'plan.json')['reused']
    old_reuse=read(out.parent/'native-joint-gr/publication.json')['reused']
    for src in sorted(work.rglob('*')):
        if not src.is_file():continue
        rel=src.relative_to(work);key='native-joint-return187-work/'+rel.as_posix()
        if key in inputs:
            prev=root/old_reuse[rel.as_posix()]['path'];assert sha(src)==sha(prev)==inputs[key]
            reused[rel.as_posix()]=dict(path=prev.relative_to(root).as_posix(),sha256=sha(src));continue
        dst=out/'completed'/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);assert sha(src)==sha(dst)
    for src,name in [(Path(__file__),'publication-producer.py'),(previous,'previous-publication-producer.py'),(b.helper,'publication-helper.py')]:shutil.copyfile(src,out/name)
    tails={
        'model-definition':'같은 해의 정확한 중심 질량 제약과 실제 Radau 각도 출구로 생성 GR의 lapse를 구성했다. 반환 성분을 별도로 유지하여 현재 물질 분기와 실제 광자 주파수·충돌 원천에 적용했다. 새 시간 적분은 아직 없다.',
        'observable-targets':'최종 전하는 여전히 미판정이다. 작은 반환 입력의 단계식 적용은 통과했지만 그 원천을 포함한 실제 결합 진화와 동일 해 전하 판독이 남는다. 입력 크기만으로 최종 전하 상계를 선언하지 않는다.',
        'adiabatic-limit':'현재 입사 광자 원천은 부호 대칭화한 주파수 drift 연산자를 사용한다. 그 부호 선형 항등식을 기호 검증했다. 물질은 현재 상태의 분기를 기준으로 하며 전역 가산성이나 이후 모든 상태의 분기 안정성을 인증하지 않는다.',
        'nonadiabatic-regime':'정확한 중심·lapse 입력의 시간 오차 최대0.0604457퍼센트로 통과했다. 두 실제 상태에서 반환을 단계 원천에 적용했고 전체 식 대조의 물질 최대1.101e-8·광자 최대1.263e-10으로 원0.002기준을 통과했다.',
        'failure-ledger-dynamic-chi':'직접 계량 덧셈은 작은 반환 성분을 최대100퍼센트 소실해 거절했다. 과거 한쪽 upwind 보상식의 재사용은 현재 부호 대칭 원천과37퍼센트 불일치했고 분기 이동만으로 설명되지 않았다. 실제 호출 함수를 재사용해 고쳤으며 두 실패와 개별80자리 대조를 보존했다. 검증한 개별 공식을 실제 연산자와 혼동하지 않는다.',
        'dynamic-charge-completion':'동일 해의 생성 GR을 정확한 중심·실제 출구 lapse로 만들고 보상 방식으로 현재 물질·광자 단계 원천에 적용했다. 이 통과는 새 결합 시간 적분이 아니며 전 기간·실제 반환 진화·동일 해 최종 전하는 계속 미완료다.'}
    prefixes={}
    for name,body in tails.items():
        f=root/f'docs/{name}.md';prefixes[f.relative_to(root).as_posix()]=dict(bytes=f.stat().st_size,sha256=sha(f))
        with f.open('ab') as h:h.write(('\n\n## 단계187 — 같은 해의 생성 GR 반환 입력\n\n분류: Counterexample candidate. '+body+' [근거](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',passed=True,joint_prefix_accepted=True,actual_same_joint_solution_GR_computed=True,
        exact_center_and_actual_packet_lapse=True,actual_current_state_stage_return_applied=True,
        metric_controls=metric['controls'],compensated_application=applied,compensated_total_seconds=total,
        direct_addition_rejected=True,wrong_frequency_owner_rejected=True,full_horizon_completed=False,
        GR_return_time_evolved=False,final_charge_conclusion='unadjudicated',physical_final_charge_solved=False,full_goal_complete=False)
    write(out/'final-result.json',final)
    write(out/'publication.json',dict(reused=reused,previous_nondoc=preserved,document_prefixes=prefixes,plan_bindings=bindings))
    files=[module,note]+[root/k for k in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(sha256={p.relative_to(root).as_posix():sha(p) for p in files},document_prefixes=prefixes)
    write(manifest,final);m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_same_joint_solution_GR_return_input']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
