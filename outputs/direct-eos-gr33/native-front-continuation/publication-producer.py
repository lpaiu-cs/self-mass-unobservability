"""Publish actual continuation and same-solution history preservation."""
from pathlib import Path
from types import FunctionType
import importlib.util,shutil,sys
source=Path(__file__).with_name('.phase179-publish.py')
if not source.exists():source=Path(__file__).with_name('previous-publication-producer.py')
spec=importlib.util.spec_from_file_location('previous',source);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha=b.read,b.write,b.sha
work=runtime/'native-front-continuation184-work';out=root/'outputs/direct-eos-gr33/native-front-continuation'
manifest=out.parent/'native-front-continuation-manifest.json';module=root/'verification/continue_full_material_front.py'
note=root/'notes/REQUEST184_SAME_SOLUTION_CONTINUATION_KO.md'


def package():
    assert not manifest.exists() and sha(module)==sha(runtime/'verification'/module.name)
    result=read(work/'time-result.json');audit=read(work/'audit-result.json');prefix=read(work/'prefix-audit.json')
    assert audit['passed'] and prefix['passed'] and read(work/'restart-check.json')['passed']
    assert read(work/'dispatch-admission.json')['source_sha256']==sha(module)
    bindings=[]
    for name in ['plan.json','restore-repair.json']:
        for p,h in read(work/name)['bindings'].items():
            target=runtime/p.removeprefix('/home/lpaiu/work/native-retained-tail-runtime/')
            if target.name==module.name:target=work/'initial-producer.py'
            assert sha(target)==h,p;bindings.append(dict(plan=name,path=p,resolved=str(target),sha256=h))
    for p,h in read(work/'alias-repair.json')['restored_runtime_files'].items():assert sha(runtime/p.replace('\\','/'))==h,p
    old=read(out.parent/'native-material-front-manifest.json');preserved={p:h for p,h in old['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    old_reuse=read(out.parent/'native-material-front/publication.json')['reused'];out.mkdir();reused={};omitted={}
    for src in sorted(work.rglob('*')):
        if not src.is_file():continue
        rel=src.relative_to(work);key=rel.as_posix()
        if key in read(work/'reuse.json'):
            row=read(work/'reuse.json')[key];source_rel=Path(row['path']).relative_to('native-material-front183-work').as_posix()
            prior=root/old_reuse[source_rel]['path'] if source_rel in old_reuse else out.parent/'native-material-front/completed'/source_rel
            assert sha(src)==sha(prior)==row['sha256'];reused[key]=dict(path=prior.relative_to(root).as_posix(),sha256=sha(src));continue
        # Reproducible zero-step serialization copies are not another result.
        # Keep their receipts and hashes; original states and both continued
        # physical paths are already retained. Leave runtime files untouched.
        if src.suffix=='.npz' and src.name.startswith(('input-','roundtrip-','restored-','identity-')):
            omitted[key]=dict(sha256=sha(src),bytes=src.stat().st_size,reason='Zero-step serialization copy; reconstruct from immutable parent checkpoints with this producer.');continue
        dst=out/'completed'/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);assert sha(src)==sha(dst)
    for src,name in [(Path(__file__),'publication-producer.py'),(source,'previous-publication-producer.py'),(Path(b.b.__file__),'publication-helper.py')]:shutil.copyfile(src,out/name)
    final=dict(result,joint_prefix_accepted=result['passed'],same_solution_audit=audit,prefix_audit=prefix,
        original_physical_prefix_recomputed=False,all_prior_state_and_history_values_preserved=True,
        zero_step_serialization_copies_omitted=True,old_failures_preserved=True,
        full_horizon_authorized=False,physical_final_charge_solved=False,gross_recovery_derivative_control_accepted=False)
    final['coarse_bytes_reused']=False;final['accepted_prefix_bytes_reused']=True
    prefixes={}
    for name,body in read(work/'publication-text.json').items():
        f=root/f'docs/{name}.md';prefixes[f.relative_to(root).as_posix()]=dict(bytes=f.stat().st_size,sha256=sha(f))
        with f.open('ab') as h:h.write(('\n\n## 단계184 — 동일 해를 보존한 실제 연속 적분\n\n분류: Counterexample candidate. '+body+' [실행과 판정](../notes/'+note.name+').\n').encode())
    assert len(prefixes)==6
    write(out/'publication.json',dict(reused=reused,omitted_reproducible_serializations=omitted,previous_nondoc=preserved,document_prefixes=prefixes,plan_bindings=bindings))
    write(out/'final-result.json',final)
    files=[module,note]+[root/k for k in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(sha256={p.relative_to(root).as_posix():sha(p) for p in files},document_prefixes=prefixes)
    write(manifest,final);m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_front_continuation']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
