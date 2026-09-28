"""Publish the actual same-solution GR connection and preserve its failures."""
from pathlib import Path
import hashlib,importlib.util,json,shutil,sys
helper=Path(__file__).with_name('.phase174-publish.py')
if not helper.exists():helper=Path(__file__).with_name('publication-helper.py')
spec=importlib.util.spec_from_file_location('helper',helper);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha,git=b.read,b.write,b.sha,b.git
work=runtime/'native-joint-gr186-work';out=root/'outputs/direct-eos-gr33/native-joint-gr'
manifest=out.parent/'native-joint-gr-manifest.json';module=root/'verification/read_full_incident_joint_gr.py'
note=root/'notes/REQUEST186_SAME_JOINT_SOLUTION_GR_KO.md'


def package():
    assert not manifest.exists() and sha(module)==sha(runtime/'verification'/module.name)
    source=read(work/'sources.json');field=read(work/'fields.json')
    assert source['passed'] and field['passed'] and read(work/'symbolic.json')['passed']
    assert all(read(work/f'{a}-receipt.json')['error'] is None for a in ['source_finish','fields'])
    assert sum(read(work/f'{a}-receipt.json')['seconds'] for a in ['source','source_retry','source_finish'])<90
    old=read(out.parent/'native-full-horizon-admission-manifest.json')
    preserved={p:h for p,h in old['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    aliases={sha(p):p for p in work.glob('*producer.py')};aliases[sha(module)]=module;bindings=[]
    for plan in work.glob('*plan.json'):
        for name,h in read(plan).get('bindings',{}).items():
            actual=runtime/name.removeprefix('/home/lpaiu/work/native-retained-tail-runtime/')
            if Path(name).name==module.name and sha(actual)!=h:actual=aliases[h]
            assert sha(actual)==h,(plan,name)
            bindings.append(dict(plan=plan.name,source=name,resolved=str(actual),sha256=h))
    out.mkdir();reused={};inputs=read(work/'plan.json')['reused']
    previous_reuse=read(out.parent/'native-full-horizon-admission/publication.json')['reused']
    for src in sorted(work.rglob('*')):
        if not src.is_file():continue
        rel=src.relative_to(work);key='native-joint-gr186-work/'+rel.as_posix()
        if key in inputs:
            prev=root/previous_reuse[rel.as_posix()]['path'];assert sha(src)==sha(prev)==inputs[key]
            reused[rel.as_posix()]=dict(path=prev.relative_to(root).as_posix(),sha256=sha(src));continue
        dst=out/'completed'/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);assert sha(src)==sha(dst)
    for n in [64,128]:
        src=runtime/f'native-full-horizon185-work/sweep-1/photons/interval-02-{n}.npz'
        dst=out/f'input/interval-02-{n}.npz';dst.parent.mkdir(exist_ok=True);shutil.copyfile(src,dst);assert sha(src)==sha(dst)
    src=runtime/'native-full-horizon185-work/comparison-02.json';shutil.copyfile(src,out/'input/comparison-02.json')
    for src,name in [(Path(__file__),'publication-producer.py'),(helper,'publication-helper.py')]:shutil.copyfile(src,out/name)
    tails={
        'model-definition':'전체 입사장으로 함께 진화한 실제 네 물질 변수·광자·출구 이력을 GR 원천으로 연결했다. Etilde 직접 복원, 현재 재고 변위와 반경 응력, 실제 비영 기하 입력을 사용한다. 다른 해의 진단·전하는 더하지 않았다.',
        'observable-targets':'최종 전하는 미판정이다. 수락된0.4293ms의 동일 해에서 실제 compact GR 응답을 계산했지만 입사파 전체의 무한대 전하나 자기 결합 완료로 해석하지 않는다.',
        'adiabatic-limit':'기존 안정적 primitive 역변환과 현재 바리온 재고 변위를 사용했다. 비정지 에너지의 기호 변환과 압력 탐침은 통과했으나 균일 EOS 미분 정리나 완전 비선형 진화의 증명은 아니다.',
        'nonadiabatic-regime':'같은 해의 짧은 구간에서 생성 GR 장까지 계산했다. GR 시간 대조0.0366408퍼센트,4/8적분 차수 대조7.827e-16,독립 직접 적분2.371e-14로 원 기준을 통과했다. 새 물리 적분은 없었다.',
        'failure-ledger-dynamic-chi':'GR 연결 구현의 namespace 누락·배열 목록 abs 오류와 SHA 별칭 불일치를 실패 소스·계획·receipt와 함께 보존했다. 수정 뒤 원천·GR 기준을 통과했고 원 예산 안에서 끝났다. 기존 물리 실패는 그대로이며 중심 J 보간과 lapse 경계는 실제 반환 전에 재구성해야 한다.',
        'dynamic-charge-completion':'동일 결합 해의 실제 GR 원천·장 연결은 통과했다. 전체 기간185는 별도 고정 실행 중이다. 생성 GR을 실제 물질·광자에 반환하고 같은 해의 질량 정규화·최종 전하를 판정하는 작업은 아직 남아 있으며 전체 목표를 완료 처리하지 않는다.'}
    prefixes={}
    for name,body in tails.items():
        f=root/f'docs/{name}.md';prefixes[f.relative_to(root).as_posix()]=dict(bytes=f.stat().st_size,sha256=sha(f))
        with f.open('ab') as h:h.write(('\n\n## 단계186 — 동일 결합 해에서 생성 GR로 연결\n\n분류: Counterexample candidate. '+body+' [근거](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',passed=True,joint_prefix_accepted=True,
        actual_same_joint_solution_GR_computed=True,physical_horizon_seconds=field['physical_horizon_seconds'],
        source_controls=source,field_controls=field['controls'],same_solution_energy_and_ports=True,
        full_horizon_completed=False,GR_return_to_matter_applied=False,final_charge_conclusion='unadjudicated',
        physical_final_charge_solved=False,full_goal_complete=False)
    write(out/'final-result.json',final)
    write(out/'publication.json',dict(reused=reused,previous_nondoc=preserved,document_prefixes=prefixes,plan_bindings=bindings))
    files=[module,note]+[root/k for k in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(sha256={p.relative_to(root).as_posix():sha(p) for p in files},document_prefixes=prefixes)
    write(manifest,final);m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_same_joint_solution_GR']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


def check(head=False):
    m=read(manifest);a=read(master)
    for p,h in m['sha256'].items():assert sha(root/p)==h and a['sha256'][p]==h,p
    for p,h in read(out/'publication.json')['previous_nondoc'].items():assert sha(root/p)==h,p
    for p,v in m['document_prefixes'].items():assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest()==v['sha256'],p
    p='docs/dynamic-charge-completion.md';assert (root/p).read_bytes().splitlines()[:20]==git('show','d5744f766719f0a29c7acb913c7984bdfd3624e2:'+p).splitlines()[:20]
    assert a['sha256'][manifest.relative_to(root).as_posix()]==sha(manifest)
    paths=list(m['sha256'])+[manifest.relative_to(root).as_posix(),'paper/revision-manifest.json']
    Path('C:/Users/lpaiu/AppData/Local/Temp/phase186-paths').write_bytes(b'\0'.join(p.encode() for p in paths)+b'\0')
    if head:
        for p in paths:assert hashlib.sha256(git('show','HEAD:'+p)).hexdigest()==sha(root/p),p
    print(json.dumps(dict(bound_files=len(m['sha256']),paths=len(paths),prefixes_preserved=len(m['document_prefixes']),same_solution_GR=m['actual_same_joint_solution_GR_computed'],final_charge_conclusion=m['final_charge_conclusion'])))


if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
