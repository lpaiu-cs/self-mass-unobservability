"""Publish the actual stage-return improvement and the surviving rejection."""
from pathlib import Path
import hashlib,importlib.util,json,shutil,sys
helper=Path(__file__).with_name('.phase174-publish.py')
if not helper.exists():helper=Path(__file__).with_name('publication-helper.py')
spec=importlib.util.spec_from_file_location('helper',helper);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha,git=b.read,b.write,b.sha,b.git
work=runtime/'native-stage-collisions178-work';out=root/'outputs/direct-eos-gr33/native-stage-collisions'
manifest=out.parent/'native-stage-collisions-manifest.json';module=root/'verification/return_native_stage_collisions.py'
note=root/'notes/REQUEST178_ACTUAL_STAGE_COLLISION_RETURN_KO.md'

def package():
    assert not manifest.exists();capture=read(work/'capture-result.json');material=read(work/'material-result.json');state=read(work/'saved-state-check.json')
    assert capture['passed'] and not material['passed'] and all(r['passed'] for r in material['rows'])
    assert material['time_comparison'][0]>.02 and not state['lagged_B_S_accepted']
    assert sha(module)==sha(runtime/'verification'/module.name)
    old=read(out.parent/'native-neutral-coupled-manifest.json');preserved={p:h for p,h in old['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    out.mkdir();reused={}
    for src in sorted(work.rglob('*')):
        if not src.is_file():continue
        rel=src.relative_to(work)
        if src.name.startswith('face-') and src.suffix=='.npz':prior=out.parent/'native-neutral-coupled/completed'/src.name
        elif rel.parts[:1]==('sweep-0',) and src.suffix=='.npz':prior=out.parent/'native-pressure-reciprocal/completed/sweep-1'/rel.parts[1]/src.name
        else:prior=None
        if prior is not None:
            assert sha(src)==sha(prior);reused[rel.as_posix()]=dict(path=prior.relative_to(root).as_posix(),sha256=sha(src));continue
        dst=out/'completed'/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);assert sha(src)==sha(dst)
    for src,name in [(Path(__file__),'publication-producer.py'),(Path(b.__file__),'publication-helper.py'),(root/'.phase178-state-check.py','saved-state-check-producer.py')]:shutil.copyfile(src,out/name)
    tails={
        'model-definition':'원 광자 해의 실제 Radau 충돌 단계를 자유 물질에 전달해 H 수송을 원0.2% 이내로 맞췄다. 그러나 바리온 시간 대조2.2826%와 B/S 상호 입력이 원 기준을 실패했다. 전체 결합·GR·최종 전하는 미판정이다.',
        'observable-targets':'최종 전하 유지 여부는 미판정이다. 실제 단계 충돌을 같은 예비 해에 적용한 수소 수송 개선은 확인했지만 전체 물질·광자 상호 일치의 수락으로 대체하지 않는다.',
        'adiabatic-limit':'알려진 충돌 누적량 C(t)에 대해 z=U+C, Udot=F(U+C)는 원 강제 물질식의 대수 변환이다. 이 적용은 stiff 광자·열·수소 Radau 방정식을 바꾸지 않는다. 변환 항등식만으로 실제 시간 정확도를 인증하지 않으며 B 시간 실패를 유지한다.',
        'nonadiabatic-regime':'같은 T/16 광자 해의 실제 충돌·gas·충격량 단계를 저장해 반환했다. 원12종 물리 배열은 정확히 재생됐고 실제 H 수송 차이는0.1349%/0.01590%로 줄었다. B 시간 대조는2.2826%로 실패했으며 장기 경로를 확대하지 않았다.',
        'failure-ledger-dynamic-chi':'177의약66% H-C 수송 불일치는 실제 충돌 단계 전달과 고정 SSPRK3 물질 적분에서 원0.2% 이내로 줄었다. 이전 실패는 보존한다. 그러나 B 시간 대조2.2826% 및 기존 광자 입력과 새 B/S의 큰 차이는 실패로 남았다. 읽기 검사 정규화 오류도 원 성분 L1으로 고쳤다. H 통과를 전체 결합 수락으로 승격하지 않는다.',
        'dynamic-charge-completion':'지배 오차를 고친 동일 해의 최종 전하라는 기준은 아직 미충족이다. 수소 수송의 실제 예비 연결은 개선됐으나 바리온 시간 정확도와 물질 운동의 상호 입력이 원 기준에 미달했다. 전체 기간·GR·전하는 미실행이다.'}
    prefixes={}
    for name,text in tails.items():
        f=root/f'docs/{name}.md';prefixes[f.relative_to(root).as_posix()]=dict(bytes=f.stat().st_size,sha256=sha(f))
        with f.open('ab') as h:h.write(('\n\n## 단계178 — 실제 충돌 이력 반환과 남은 물질 실패\n\n분류: Counterexample candidate. '+text+' [실제 실행과 남은 판정](../notes/'+note.name+').\n').encode())
    write(out/'publication.json',dict(reused=reused,previous_nondoc=preserved,document_prefixes=prefixes))
    final=dict(classification='Counterexample candidate',passed=False,original_physical_arrays_exactly_replayed=True,
        actual_stage_collisions_returned=True,neutral_transport_original_gate_passed=True,
        mechanical_H_relative=[v['mechanical_H_relative'] for v in material['rows']],
        native_H_rate_relative=[v['native_H_rate_relative'] for v in material['rows']],
        material_time_comparison=material['time_comparison'],material_time_accepted=False,
        lagged_B_S_relative=[v['lagged_B_S_relative'] for v in state['rows']],lagged_material_input_accepted=False,
        full_horizon_authorized=False,full_horizon_completed=False,GR_charge_readout_executed=False,
        final_charge_conclusion='unadjudicated',physical_final_charge_solved=False,full_goal_complete=False)
    write(out/'final-result.json',final)
    files=[module,note]+[root/k for k in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(sha256={p.relative_to(root).as_posix():sha(p) for p in files},document_prefixes=prefixes)
    write(manifest,final);m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_stage_collision_return']={k:v for k,v in final.items() if k!='sha256'};write(master,m)

def check(head=False):
    m=read(manifest);a=read(master)
    for p,h in m['sha256'].items():assert sha(root/p)==h and a['sha256'][p]==h,p
    for p,h in read(out/'publication.json')['previous_nondoc'].items():assert sha(root/p)==h,p
    for p,v in m['document_prefixes'].items():assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest()==v['sha256'],p
    p='docs/dynamic-charge-completion.md';assert (root/p).read_bytes().splitlines()[:20]==git('show','d5744f766719f0a29c7acb913c7984bdfd3624e2:'+p).splitlines()[:20]
    assert a['sha256'][manifest.relative_to(root).as_posix()]==sha(manifest)
    paths=list(m['sha256'])+[manifest.relative_to(root).as_posix(),'paper/revision-manifest.json']
    Path('C:/Users/lpaiu/AppData/Local/Temp/phase178-paths').write_bytes(b'\0'.join(p.encode() for p in paths)+b'\0')
    if head:
        for p in paths:assert hashlib.sha256(git('show','HEAD:'+p)).hexdigest()==sha(root/p),p
    print(json.dumps(dict(bound_files=len(m['sha256']),paths=len(paths),prefixes_preserved=6,neutral_transport_accepted=True,full_coupling_accepted=False,final_charge_conclusion='unadjudicated')))

if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
