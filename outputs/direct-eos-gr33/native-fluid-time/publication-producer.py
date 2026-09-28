"""Publish the original time pair; reuse the prior binding verifier."""
from pathlib import Path
from types import FunctionType
import importlib.util,shutil,sys
source=Path(__file__).with_name('.phase179-publish.py')
if not source.exists():source=Path(__file__).with_name('previous-publication-producer.py')
spec=importlib.util.spec_from_file_location('previous',source);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha=b.read,b.write,b.sha
work=runtime/'native-fluid-time180-work';out=root/'outputs/direct-eos-gr33/native-fluid-time'
manifest=out.parent/'native-fluid-time-manifest.json';module=root/'verification/complete_native_fluid_time.py'
note=root/'notes/REQUEST180_JOINT_FLUID_TIME_KO.md'


def package():
    assert not manifest.exists() and sha(module)==sha(runtime/'verification'/module.name)
    pair=read(work/'pilot-result.json');audit=read(work/'audit-result.json');assert audit['passed']
    for p,h in read(work/'plan.json')['bindings'].items():assert sha(runtime/p.removeprefix('/home/lpaiu/work/native-retained-tail-runtime/'))==h,p
    old=read(out.parent/'native-fluid-radau-manifest.json');preserved={p:h for p,h in old['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    out.mkdir();reused={};old_reuse=read(out.parent/'native-fluid-radau/publication.json')['reused']
    for src in sorted(work.rglob('*')):
        if not src.is_file():continue
        rel=src.relative_to(work);key=rel.as_posix()
        if key in read(work/'reuse.json'):
            prior=root/old_reuse[key]['path'] if key in old_reuse else out.parent/'native-fluid-radau/completed'/rel
            assert sha(src)==sha(prior);reused[key]=dict(path=prior.relative_to(root).as_posix(),sha256=sha(src));continue
        dst=out/'completed'/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);assert sha(src)==sha(dst)
    for src,name in [(Path(__file__),'publication-producer.py'),(source,'previous-publication-producer.py'),(Path(b.b.__file__),'publication-helper.py'),(root/'.phase180-scope.py','superposition-boundary-producer.py'),(root/'.phase180-baryon.py','baryon-localization-producer.py')]:shutil.copyfile(src,out/name)
    final=dict(pair,joint_prefix_accepted=pair['passed'],same_solution_audit_passed=True,
        original_coarse_bytes_unchanged=True,full_horizon_authorized=False,physical_final_charge_solved=False,
        drive_scope='redshift_correction_only_with_zero_additional_metric',full_physical_drive_evolved=False,
        solution_superposition_certified=False,next_mainline='Apply the full physical driver and corrected source to one joint equation before final charge readout; do not prioritize a full correction-only path.')
    write(out/'final-result.json',final)
    assert not pair['passed'] and pair['material_time'][0]>.02 and max(pair['photon_time'])<.02
    tails={
        'model-definition':'원128경로를 실제 완료했다. 광자와 네 물질량의 동시 단계식·보존은 두 시계에서 통과했으나 B 시간 차이2.8667%가 원2%를 실패했다. 구동 범위는 추가 계량0의 적색편이 보정 원천이며 전체 물리 구동은 아직 같은 식에 들어 있지 않다.',
        'observable-targets':'최종 전하는 미판정이다. 바리온 시간 실패를 유지하고 전체 기간을 확대하지 않았다. 보정 응답을 과거 선택 전하에 더하는 방식은 원 물질 연산자의 분기에서 중첩이 보장되지 않아 최종 판정 경로가 아니다.',
        'adiabatic-limit':'현재 방향의 직접 Jacobian이 원 단계식을 통과해도 물질 방향 연산자의 전역 중첩은 따라오지 않는다. 원 minmod의 영 기울기에서 비가산 반례를 실행 확인했다. 전체 구동과 수정 원천을 같은 식에 적용하는 경로가 필요하다.',
        'nonadiabatic-regime':'원64/128시계의 동일 T/16 보정 응답은 모든 여섯 광자/열 채널 시간 차이가2% 안이지만 B가2.8667%로 실패했다. 운동량0.5572%,에너지0.4461%,H0.3805%다. 실제 단계·출구·제거량 장부를 저장해 같은 해 감사는 통과했다.',
        'failure-ledger-dynamic-chi':'128실제12개 하위 단계와64저장6개 하위 단계는 잔차·보존 기준을 통과했지만 B 시간 대조2.8667%는 실패했다. B 차이의90.82%가 대기143–149셀에 있고 floor 제거 차이의 L1은 전체 B 차이 L1의1.39e-7배다. 유속 누적 차이를 원인 항목으로 좁혔으나 그 시간 오차의 유일 원인은 미확정이다. 추가 시계·전체 기간은 실행하지 않는다.',
        'dynamic-charge-completion':'최종 전하 결론은 아직 미판정이다. 현재까지 연결한 해는 적색편이 보정 원천의 응답이며 전체 물리 구동과 자기 GR 반환을 닫은 해가 아니다. 이 보정 해의 B 시간 실패도 유지한다. 전체 구동과 수정 수송을 같은 방정식에서 풀고 그 해의 에너지·경계·현재 계량으로 전하를 판독하는 경로를 우선한다.'}
    prefixes={}
    for name,text in tails.items():
        f=root/f'docs/{name}.md';prefixes[f.relative_to(root).as_posix()]=dict(bytes=f.stat().st_size,sha256=sha(f))
        with f.open('ab') as h:h.write(('\n\n## 단계180 — 실제 동시 해의 원 시간 대조\n\n분류: Counterexample candidate. '+text+' [실행과 판정](../notes/'+note.name+').\n').encode())
    assert len(prefixes)==6
    write(out/'publication.json',dict(reused=reused,previous_nondoc=preserved,document_prefixes=prefixes))
    files=[module,note]+[root/k for k in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(sha256={p.relative_to(root).as_posix():sha(p) for p in files},document_prefixes=prefixes)
    write(manifest,final);m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_fluid_joint_time']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
