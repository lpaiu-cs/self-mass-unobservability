"""Publish the full-input joint trial with its failed controls preserved."""
from pathlib import Path
from types import FunctionType
import importlib.util,shutil,sys
source=Path(__file__).with_name('.phase179-publish.py')
if not source.exists():source=Path(__file__).with_name('previous-publication-producer.py')
spec=importlib.util.spec_from_file_location('previous',source);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha=b.read,b.write,b.sha
work=runtime/'native-full-incident181-work';out=root/'outputs/direct-eos-gr33/native-full-incident'
manifest=out.parent/'native-full-incident-manifest.json';module=root/'verification/couple_full_incident_fluid.py'
note=root/'notes/REQUEST181_FULL_INCIDENT_JOINT_KO.md'


def package():
    assert not manifest.exists() and sha(module)==sha(runtime/'verification'/module.name)
    receipt=read(work/'pilot-receipt.json');control=read(work/'check-result.json')
    assert control['passed'] and sha(module)==read(work/'chain-control-plan.json')['source_sha256']
    bindings=[]
    for p,h in read(work/'plan.json')['bindings'].items():
        target=work/'initial-producer.py' if Path(p).name==module.name else runtime/p.removeprefix('/home/lpaiu/work/native-retained-tail-runtime/')
        assert sha(target)==h,p;bindings.append(dict(path=str(p),sha256=h))
    old=read(out.parent/'native-fluid-time-manifest.json');preserved={p:h for p,h in old['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    audit=read(work/'saved-audit.json');assert audit['passed']
    out.mkdir();reused={};old_reuse=read(out.parent/'native-fluid-radau/publication.json')['reused']
    for src in sorted(work.rglob('*')):
        if not src.is_file():continue
        rel=src.relative_to(work);key=rel.as_posix()
        if key in read(work/'reuse.json'):
            row=old_reuse[key];assert sha(src)==sha(root/row['path'])==row['sha256'];reused[key]=row;continue
        dst=out/'completed'/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);assert sha(src)==sha(dst)
    for src,name in [(Path(__file__),'publication-producer.py'),(source,'previous-publication-producer.py'),(Path(b.b.__file__),'publication-helper.py')]:shutil.copyfile(src,out/name)
    result=read(work/'pilot-result.json') if (work/'pilot-result.json').exists() else dict(passed=False,error=receipt['error'])
    final=dict(result,classification='Counterexample candidate',joint_prefix_accepted=bool(result['passed'] and receipt['error'] is None),
        full_horizon_authorized=False,physical_final_charge_solved=False,final_charge_conclusion='unadjudicated',full_goal_complete=False,
        drive_scope='full_declared_incident_primary_plus_saved_first_Born_with_conservative_redshift',
        self_GR_return_closed=False,gross_recovery_derivative_control_accepted=False,
        original180_baryon_time_failure_preserved=True,check=control,receipt=receipt,same_solution_audit=audit)
    text=read(work/'publication-text.json')
    prefixes={}
    for name,body in text.items():
        f=root/f'docs/{name}.md';prefixes[f.relative_to(root).as_posix()]=dict(bytes=f.stat().st_size,sha256=sha(f))
        with f.open('ab') as h:h.write(('\n\n## 단계181 — 전체 입사 구동의 동일 단계 연결\n\n분류: Counterexample candidate. '+body+' [실행과 판정](../notes/'+note.name+').\n').encode())
    assert len(prefixes)==6
    write(out/'publication.json',dict(reused=reused,previous_nondoc=preserved,document_prefixes=prefixes,plan_bindings=bindings))
    write(out/'final-result.json',final)
    files=[module,note]+[root/k for k in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(sha256={p.relative_to(root).as_posix():sha(p) for p in files},document_prefixes=prefixes)
    write(manifest,final);m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_full_incident_joint']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
