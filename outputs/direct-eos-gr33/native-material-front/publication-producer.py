"""Publish the material-front repair, preserving original time-pair failure."""
from pathlib import Path
from types import FunctionType
import importlib.util,shutil,sys
source=Path(__file__).with_name('.phase179-publish.py')
if not source.exists():source=Path(__file__).with_name('previous-publication-producer.py')
spec=importlib.util.spec_from_file_location('previous',source);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha=b.read,b.write,b.sha
work=runtime/'native-material-front183-work';out=root/'outputs/direct-eos-gr33/native-material-front'
manifest=out.parent/'native-material-front-manifest.json';module=root/'verification/resolve_full_material_front.py'
note=root/'notes/REQUEST183_NATIVE_MATERIAL_FRONT_KO.md'


def package():
    assert not manifest.exists() and sha(module)==sha(runtime/'verification'/module.name)
    result=read(work/'time-result.json');audit=read(work/'audit-result.json');assert audit['passed']
    assert read(work/'equivalence.json')['passed']
    bindings=[]
    for p,h in read(work/'plan.json')['bindings'].items():
        target=runtime/p.removeprefix('/home/lpaiu/work/native-retained-tail-runtime/')
        assert sha(target)==h,p;bindings.append(dict(path=p,sha256=h))
    old=read(out.parent/'native-full-time-manifest.json');preserved={p:h for p,h in old['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    old_reuse=read(out.parent/'native-full-time/publication.json')['reused'];out.mkdir();reused={}
    for src in sorted(work.rglob('*')):
        if not src.is_file():continue
        rel=src.relative_to(work);key=rel.as_posix()
        if key in read(work/'reuse.json'):
            row=read(work/'reuse.json')[key];source_rel=Path(row['path']).relative_to('native-full-time182-work').as_posix()
            prior=root/old_reuse[source_rel]['path'] if source_rel in old_reuse else out.parent/'native-full-time/completed'/source_rel
            assert sha(src)==sha(prior)==row['sha256'];reused[key]=dict(path=prior.relative_to(root).as_posix(),sha256=sha(src));continue
        dst=out/'completed'/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);assert sha(src)==sha(dst)
    for src,name in [(Path(__file__),'publication-producer.py'),(source,'previous-publication-producer.py'),(Path(b.b.__file__),'publication-helper.py')]:shutil.copyfile(src,out/name)
    final=dict(result,joint_prefix_accepted=result['passed'],same_solution_audit=audit,
        coarse_saved_run_base_clock=128,coarse_equivalent_new_base_clock=64,coarse_actual_substeps=4,
        fine_base_clock=128,fine_actual_substeps=8,old_time_failure_preserved=True,
        full_horizon_authorized=False,physical_final_charge_solved=False,gross_recovery_derivative_control_accepted=False)
    prefixes={}
    for name,body in read(work/'publication-text.json').items():
        f=root/f'docs/{name}.md';prefixes[f.relative_to(root).as_posix()]=dict(bytes=f.stat().st_size,sha256=sha(f))
        with f.open('ab') as h:h.write(('\n\n## 단계183 — 실제 물질 전선의 시간 분해\n\n분류: Counterexample candidate. '+body+' [실행과 판정](../notes/'+note.name+').\n').encode())
    assert len(prefixes)==6
    write(out/'publication.json',dict(reused=reused,previous_nondoc=preserved,document_prefixes=prefixes,plan_bindings=bindings))
    write(out/'final-result.json',final)
    files=[module,note]+[root/k for k in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(sha256={p.relative_to(root).as_posix():sha(p) for p in files},document_prefixes=prefixes)
    write(manifest,final);m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_material_front']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
