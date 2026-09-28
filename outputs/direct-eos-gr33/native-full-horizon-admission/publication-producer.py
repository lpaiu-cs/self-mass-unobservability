"""Freeze admission evidence before the bounded full-horizon producer starts."""
from pathlib import Path
from types import FunctionType
import importlib.util,shutil,sys
source=Path(__file__).with_name('.phase179-publish.py')
if not source.exists():source=Path(__file__).with_name('previous-publication-producer.py')
spec=importlib.util.spec_from_file_location('previous',source);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha=b.read,b.write,b.sha
work=runtime/'native-full-horizon185-work';out=root/'outputs/direct-eos-gr33/native-full-horizon-admission'
manifest=out.parent/'native-full-horizon-admission-manifest.json';module=root/'verification/complete_full_incident_horizon.py'
note=root/'notes/REQUEST185_FULL_INCIDENT_HORIZON_KO.md'


def package():
    assert not manifest.exists() and sha(module)==sha(runtime/'verification'/module.name)
    plan=read(work/'plan.json');gate=read(work/'admission-check.json')
    assert gate['passed'] and read(work/'check-receipt.json')['error'] is None
    assert read(work/'status.json')['state']=='prepared' and not (work/'run-receipt.json').exists()
    bindings=[]
    for p,h in plan['bindings'].items():
        target=runtime/p.removeprefix('/home/lpaiu/work/native-retained-tail-runtime/')
        assert sha(target)==h,p;bindings.append(dict(path=p,sha256=h))
    old=read(out.parent/'native-front-continuation-manifest.json');preserved={p:h for p,h in old['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    old_reuse=read(out.parent/'native-front-continuation/publication.json')['reused'];out.mkdir();reused={}
    for src in sorted(work.rglob('*')):
        if not src.is_file():continue
        rel=src.relative_to(work);key=rel.as_posix()
        if key in read(work/'reuse.json'):
            row=read(work/'reuse.json')[key];source_rel=Path(row['path']).relative_to('native-front-continuation184-work').as_posix()
            prev=root/old_reuse[source_rel]['path'] if source_rel in old_reuse else out.parent/'native-front-continuation/completed'/source_rel
            assert sha(src)==sha(prev)==row['sha256'];reused[key]=dict(path=prev.relative_to(root).as_posix(),sha256=sha(src));continue
        dst=out/'prepared'/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);assert sha(src)==sha(dst)
    for src,name in [(Path(__file__),'publication-producer.py'),(source,'previous-publication-producer.py'),(Path(b.b.__file__),'publication-helper.py')]:shutil.copyfile(src,out/name)
    tails={
        'model-definition':'동일한184의 전체 외부 입사 결합 해를 원3.43443ms까지 이어갈 계획을 고정했다. 원 EOS·배경·대기·방정식·파형·시계·한 번의 국소 분할을 유지하며 새로운 해나 과거 전하를 더하지 않는다.',
        'observable-targets':'최종 전하 결론은 미판정이다. 이번 전체 기간의 물질·광자 응답은 그 해의 생성 GR과 질량·전하 판독을 위한 필수 입력이다. 전체 적분 통과를 최종 관측량 판정으로 바꾸지 않는다.',
        'adiabatic-limit':'선택된 분기에서 푸는 Newton 제안과 원 native 단계식 판정을 유지한다. 현재 경로는 처방된 입사장에 대한 retained 응답이며 보편적 미분 정리나 완전 비선형 항성 해가 아니다.',
        'nonadiabatic-regime':'이미 수락한0.21465ms와 모든 누적 이력을 재사용한다. 새 실제 하위 단계는111개와215개다. 실측 예상 약2시간16분, 가정상 여유 비용 약4시간, 총 실행 상한4시간30분으로 등록했다.',
        'failure-ledger-dynamic-chi':'각 정규 출력에서 원 여섯 광자/열·네 물질량2퍼센트 시간 기준과 같은 해의 물리 단계·구성식·보존·출구·앞부분 보존을 검사한다. 실패나 비용 상한 초과 즉시 멈추며 자동 재시도·추가 분할·해상도·기간·sweep 확대를 하지 않는다.',
        'dynamic-charge-completion':'사용자 기준은 지배 오차를 수정한 동일 해의 최종 전하다. 전체 기간을 향한 유한 예산을 등록했으며 시작 전 저장 해의 단계·수지·출구 대조가 통과했다. 전체 기간·생성 GR 반환·동일 해 전하 판독은 아직 완료되지 않았다.'}
    prefixes={}
    for name,body in tails.items():
        f=root/f'docs/{name}.md';prefixes[f.relative_to(root).as_posix()]=dict(bytes=f.stat().st_size,sha256=sha(f))
        with f.open('ab') as h:h.write(('\n\n## 단계185 — 동일 해의 전체 기간 계산 등록\n\n분류: Conjectural. '+body+' [고정 계획](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',passed=True,scope='Admission and saved-prefix verification only; no185physical steps yet.',
        joint_prefix_accepted=True,full_horizon_admitted=True,full_horizon_completed=False,GR_return_closed=False,
        final_charge_conclusion='unadjudicated',physical_final_charge_solved=False,full_goal_complete=False,
        forecast=plan['forecast'],budgets=plan['budgets'],actual_steps=plan['actual_steps'],admission_check=gate)
    write(out/'publication.json',dict(reused=reused,previous_nondoc=preserved,document_prefixes=prefixes,plan_bindings=bindings))
    write(out/'admission-result.json',final)
    files=[module,note]+[root/k for k in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(sha256={p.relative_to(root).as_posix():sha(p) for p in files},document_prefixes=prefixes)
    write(manifest,final);m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_full_horizon_admission']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
