"""Bind actual joint-stage evidence without promoting it to a charge verdict."""
from pathlib import Path
import hashlib,importlib.util,json,shutil,sys
helper=Path(__file__).with_name('.phase174-publish.py')
if not helper.exists():helper=Path(__file__).with_name('publication-helper.py')
spec=importlib.util.spec_from_file_location('helper',helper);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha,git=b.read,b.write,b.sha,b.git
work=runtime/'native-fluid-radau179-work';out=root/'outputs/direct-eos-gr33/native-fluid-radau'
manifest=out.parent/'native-fluid-radau-manifest.json';module=root/'verification/couple_native_fluid_radau.py'
note=root/'notes/REQUEST179_JOINT_NATIVE_FLUID_RADAU_KO.md'


def package():
    assert not manifest.exists() and (work/'branch64-receipt.json').exists()
    assert sha(module)==sha(runtime/'verification'/module.name)
    old=read(out.parent/'native-stage-collisions-manifest.json');preserved={p:h for p,h in old['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    aliases={sha(p):p for p in work.glob('*producer.py')};aliases[sha(module)]=module;bindings=[]
    for plan in work.glob('*plan.json'):
        for name,h in read(plan).get('bindings',{}).items():
            rel=name.removeprefix('/home/lpaiu/work/native-retained-tail-runtime/')
            actual=runtime/rel
            if Path(name).name==module.name and (not actual.exists() or sha(actual)!=h):actual=aliases[h]
            assert sha(actual)==h,(plan,name)
            bindings.append(dict(plan=plan.name,source=name,resolved=str(actual),sha256=h))
    out.mkdir();reused={}
    for src in sorted(work.rglob('*')):
        if not src.is_file():continue
        rel=src.relative_to(work)
        if rel.parts[:1]==('sweep-0',) and src.suffix=='.npz':
            prior=out.parent/'native-pressure-reciprocal/completed/sweep-1'/rel.parts[1]/src.name
            assert sha(src)==sha(prior);reused[rel.as_posix()]=dict(path=prior.relative_to(root).as_posix(),sha256=sha(src));continue
        dst=out/'completed'/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);assert sha(src)==sha(dst)
    for src,name in [(Path(__file__),'publication-producer.py'),(Path(b.__file__),'publication-helper.py')]:shutil.copyfile(src,out/name)
    coarse=read(work/'sweep-1/photons/pilot-64.json');fine=read(work/'fine-admission.json')
    assert coarse['passed'] and coarse['actual_B_S_E_H_joint_unknowns'] and not coarse['lagged_material_input_used']
    assert not fine['eligible'] and not (work/'sweep-1/photons/pilot-128.npz').exists()
    assert read(work/'symbolic.json')['passed'] and read(work/'branch64-receipt.json')['error'] is None
    final=dict(classification='Counterexample candidate',passed=False,coarse_joint_prefix_accepted=True,joint_prefix_accepted=False,
        original_fixed_support_rejected=True,original_finite_Jacobian_stage_rejected=True,
        selected_branch_Jacobian_applied_to_actual_evolution=True,actual_B_S_E_H_joint_unknowns=True,
        lagged_material_input_used=False,coarse_completed_macro_steps=coarse['completed_steps'],coarse_completed_substeps=coarse['actual_completed_steps'],
        true_stage_relative=coarse['native_true_equation'],physical_stage_relative=coarse['native_true_moments'],
        joint_material_ledger_relative=coarse['joint_material_ledger_relative'],
        energy_balance_relative=coarse['energy_balance_relative'],species_balance_relative=coarse['species_balance_relative'],
        fine_admission=fine,time_comparison_completed=False,full_horizon_authorized=False,full_horizon_completed=False,
        GR_charge_readout_executed=False,final_charge_conclusion='unadjudicated',physical_final_charge_solved=False,full_goal_complete=False)
    tails={
        'model-definition':'광자와 B/S/Etilde/H를 같은 실제 Radau 단계에서 푼64시계의 T/16 해가 단계·미분·수지 기준을 통과했다. B/S와 재고 변위를 현재 미지수에서 원천·압력에 전달하며 이전 물질 이력을 고정 입력으로 사용하지 않는다.128시간 대조와 전체 GR·전하는 미완료다.',
        'observable-targets':'최종 전하 결론은 미판정이다. 실제6개 하위 단계에서 물질과 광자의 동시 식이 통과한 것을 전하 결론의 유지나 새로운 관측량 검출로 대신하지 않는다.',
        'adiabatic-limit':'분기별 직접 열 작용은 실제 방향 minmod·HLL·donor의 현재 선택을 고정한 Newton 제안이다. 원 native 방정식으로 수락하며 전역 선형성이나 균일 미분 정리로 확대하지 않는다. 활동영역 floor 투영은 같은 해의 별도 제거 수지에 남겼다.',
        'nonadiabatic-regime':'수정된 연산자를 실제64경로에 적용해6개 하위 단계를 완료했다. 최종 물리 단계 잔차 최대4.52e-15, 같은 해의 네 물질 수지 잔차 최대2.45e-15였다.128경로의 실측 기반 추정421초가 원 예산의 남은210초를 넘어 시간 대조는 실행하지 않았다.',
        'failure-ledger-dynamic-chi':'고정 물질 활동영역 가정은 기각했고 실제 끝점 floor 투영을 복원했다. 유한 차분 Jacobian의 첫 실제 단계는 운동량 잔차6.22e-10으로 원1e-13을 실패했다. 현재 분기를 고정한 직접 열 방식은 같은 원 식의 실제64경로6개 하위 단계에서 최대4.52e-15로 통과했다. 원 실패·최대3회 반복·정밀도·수락 기준을 보존했다.128경로는 예산 입장에서 거절됐으며 물리적 시간 실패나 통과로 분류하지 않는다.',
        'dynamic-charge-completion':'지배 오차를 고친 동일 결합 해의 최종 전하라는 사용자 기준은 아직 미충족이다. 실제 광자·네 물질 보존량의 동시 예비 해는64경로에서 통과했지만 시간 대조·전체 기간·실제 계량 반환·전하 판독은 미완료다. 이전 선택 전하의 부호를 계승하지 않는다.'}
    prefixes={}
    assert final['final_charge_conclusion']=='unadjudicated' and not final['full_goal_complete']
    for name,text in tails.items():
        f=root/f'docs/{name}.md';prefixes[f.relative_to(root).as_posix()]=dict(bytes=f.stat().st_size,sha256=sha(f))
        with f.open('ab') as h:h.write(('\n\n## 단계179 — 동일 광자·물질 단계의 실제 판정\n\n분류: Counterexample candidate. '+text+' [실행과 판정](../notes/'+note.name+').\n').encode())
    assert len(prefixes)==6
    write(out/'publication.json',dict(reused=reused,previous_nondoc=preserved,document_prefixes=prefixes,plan_bindings=bindings))
    write(out/'final-result.json',final)
    files=[module,note]+[root/k for k in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(sha256={p.relative_to(root).as_posix():sha(p) for p in files},document_prefixes=prefixes)
    write(manifest,final);m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_fluid_joint_radau']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


def check(head=False):
    m=read(manifest);a=read(master)
    for p,h in m['sha256'].items():assert sha(root/p)==h and a['sha256'][p]==h,p
    for p,h in read(out/'publication.json')['previous_nondoc'].items():assert sha(root/p)==h,p
    for p,v in m['document_prefixes'].items():assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest()==v['sha256'],p
    p='docs/dynamic-charge-completion.md';assert (root/p).read_bytes().splitlines()[:20]==git('show','d5744f766719f0a29c7acb913c7984bdfd3624e2:'+p).splitlines()[:20]
    assert a['sha256'][manifest.relative_to(root).as_posix()]==sha(manifest)
    paths=list(m['sha256'])+[manifest.relative_to(root).as_posix(),'paper/revision-manifest.json']
    Path('C:/Users/lpaiu/AppData/Local/Temp/phase179-paths').write_bytes(b'\0'.join(p.encode() for p in paths)+b'\0')
    if head:
        for p in paths:assert hashlib.sha256(git('show','HEAD:'+p)).hexdigest()==sha(root/p),p
    print(json.dumps(dict(bound_files=len(m['sha256']),paths=len(paths),prefixes_preserved=6,final_charge_conclusion=m['final_charge_conclusion'],joint_prefix_accepted=m['joint_prefix_accepted'])))


if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
