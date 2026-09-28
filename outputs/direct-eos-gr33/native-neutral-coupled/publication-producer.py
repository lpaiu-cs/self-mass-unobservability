"""Bind the actual177 prefix and rejected free-material transport without promotion."""
from pathlib import Path
import hashlib,importlib.util,json,shutil,sys
helper=Path(__file__).with_name('.phase174-publish.py')
if not helper.exists():helper=Path(__file__).with_name('publication-helper.py')
spec=importlib.util.spec_from_file_location('helper',helper);base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
root,runtime,master=base.root,base.runtime,base.master
read,write,sha,git=base.read,base.write,base.sha,base.git
work=runtime/'native-neutral-coupled177-work';out=root/'outputs/direct-eos-gr33/native-neutral-coupled'
manifest=out.parent/'native-neutral-coupled-manifest.json'
modules=['couple_native_neutral_transport.py','complete_native_neutral_transport.py']
note=root/'notes/REQUEST177_NATIVE_NEUTRAL_COUPLING_KO.md'

def package():
    assert not manifest.exists();p=read(work/'pilot-result.json');m=read(work/'material-prefix-result.json')
    assert p['passed'] and not m['passed'] and not p['full_horizon_completed']
    for name in modules:assert sha(root/'verification'/name)==sha(runtime/'verification'/name)
    old=read(out.parent/'native-pressure-reciprocal-manifest.json');preserved={k:v for k,v in old['sha256'].items() if not k.startswith('docs/')}
    for k,v in preserved.items():assert sha(root/k)==v,k
    out.mkdir();reused={}
    for src in sorted(work.rglob('*')):
        if not src.is_file():continue
        rel=src.relative_to(work)
        if rel.parts[:1]==('sweep-0',) and src.suffix=='.npz':
            prior=out.parent/'native-pressure-reciprocal/completed/sweep-1'/rel.parts[1]/src.name
            assert sha(src)==sha(prior);reused[rel.as_posix()]=dict(path=prior.relative_to(root).as_posix(),sha256=sha(src));continue
        dst=out/'completed'/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);assert sha(src)==sha(dst)
    shutil.copyfile(__file__,out/'publication-producer.py');shutil.copyfile(base.__file__,out/'publication-helper.py')
    prefixes={};tails={
        'model-definition':'원 native H 수송을 실제 광자·열·수소 Radau 방정식에 넣은 예비 두 경로는 통과했다. 그러나 같은 광자 이력의 자유 물질 반환에서 H-C 수송이 원 기준을 실패했으므로 동일 결합 해 수락과 최종 전하는 미판정이다.',
        'observable-targets':'최종 전하 결론의 유지 여부는 미판정이다. 단계 방정식의 수치 수락을 최종 관측량 성과로 대체하지 않는다. 실제 자유 물질 반환의 수송 불일치 때문에 전체 기간·GR 판독은 실행하지 않았다.',
        'adiabatic-limit':'Hdot=Cdot+T_H이면 (H-C)dot=T_H인 대수 항등식을 실제 수송 적분 장부로 확인했다. 그러나 다른 시간 보간을 소비하는 자유 물질 해와 수송이 같다는 보장은 없으며, 그 실패는 정적 흡수나 no-go 판정이 아니다.',
        'nonadiabatic-regime':'동일T/16 구간의 실제6/12 Radau 하위 단계는 원2% 시간·단계·보존 기준을 통과했다. 실제 반환 물질의 H 수송은 약66% 불일치해 원0.2%를 실패했다. 실제 단계별 충돌 전달이 다음 수리 대상이다.',
        'failure-ledger-dynamic-chi':'전역 선형 유속 인증의 첫 절점 부호 반전 실패를 보존했다. 행렬은 제안 연산자로만 쓰고 원 native 단계 잔차를 검사한 예비 해는 통과했으나 실제 자유 물질 반환에서 H-C 수송 차이66.05%/65.97%로 실패했다. 입력 초기화·guard 실패도 보존했다. 전체 기간을 실행하지 않았다.',
        'dynamic-charge-completion':'지배 오차를 수정한 동일 결합 해의 최종 전하라는 기준은 아직 미충족이다. 실제 수송을 결합 단계에 넣었으나 자유 물질과의 수송 일치는 실패했다. 전체 기간·GR·전하는 미실행이며 이전 선택 전하 부호를 계승하지 않는다.'}
    for name,body in tails.items():
        f=root/f'docs/{name}.md';prefixes[f.relative_to(root).as_posix()]=dict(bytes=f.stat().st_size,sha256=sha(f))
        with f.open('ab') as h:h.write(('\n\n## 단계177 — 단계 내 H 수송과 실제 반환 실패\n\n분류: Counterexample candidate. '+body+' [실행과 판정](../notes/'+note.name+').\n').encode())
    write(out/'publication.json',dict(reused=reused,previous_nondoc=preserved,document_prefixes=prefixes))
    final=dict(classification='Counterexample candidate',passed=False,actual_native_H_Radau_prefix_passed=True,
        photon_time_comparison=p['time_comparison'],native_equation_relative=max(r['native_equation_relative'] for r in p['rows']),
        native_partition_relative=max(r['neutral_partition_relative'] for r in p['rows']),
        actual_free_material_prefix_completed=True,free_material_transport_accepted=False,
        free_material_mechanical_H_relative=[r['mechanical_H_relative'] for r in m['rows']],
        free_material_native_rate_relative=[r['native_H_rate_relative'] for r in m['rows']],
        full_horizon_authorized=False,GR_charge_readout_executed=False,
        final_charge_conclusion='unadjudicated',physical_final_charge_solved=False,full_goal_complete=False)
    write(out/'final-result.json',final)
    files=[root/'verification'/n for n in modules]+[note]+[root/k for k in prefixes]+[f for f in out.rglob('*') if f.is_file()]
    final.update(sha256={f.relative_to(root).as_posix():sha(f) for f in files},document_prefixes=prefixes)
    write(manifest,final);a=read(master);a['sha256'].update(final['sha256']);a['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    a['native_neutral_coupled']={k:v for k,v in final.items() if k!='sha256'};write(master,a)

def check(head=False):
    a=read(manifest);b=read(master)
    for p,h in a['sha256'].items():assert sha(root/p)==h and b['sha256'][p]==h,p
    for p,h in read(out/'publication.json')['previous_nondoc'].items():assert sha(root/p)==h,p
    for p,v in a['document_prefixes'].items():assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest()==v['sha256'],p
    p='docs/dynamic-charge-completion.md';assert (root/p).read_bytes().splitlines()[:20]==git('show','d5744f766719f0a29c7acb913c7984bdfd3624e2:'+p).splitlines()[:20]
    assert b['sha256'][manifest.relative_to(root).as_posix()]==sha(manifest)
    paths=list(a['sha256'])+[manifest.relative_to(root).as_posix(),'paper/revision-manifest.json']
    Path('C:/Users/lpaiu/AppData/Local/Temp/phase177-paths').write_bytes(b'\0'.join(p.encode() for p in paths)+b'\0')
    if head:
        for p in paths:assert hashlib.sha256(git('show','HEAD:'+p)).hexdigest()==sha(root/p),p
    print(json.dumps(dict(bound_files=len(a['sha256']),paths=len(paths),preserved_prefixes=6,final_charge_conclusion='unadjudicated')))

if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
