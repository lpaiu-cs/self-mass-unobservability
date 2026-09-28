"""Preserve strict recovery failures and the failed original GR-source gate."""
from pathlib import Path
from types import FunctionType
import importlib.util,shutil,sys
p=Path(__file__).with_name('.phase186-publish.py')
if not p.exists():p=Path(__file__).with_name('previous-publication-producer.py')
spec=importlib.util.spec_from_file_location('previous',p);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha=b.read,b.write,b.sha
out=root/'outputs/direct-eos-gr33/native-resolved-return';manifest=out.parent/'native-resolved-return-manifest.json'
folders=[runtime/n for n in ['native-stage-recovery204-work','native-exact-stage205-work','native-resolved-return206-work']]
modules=[root/'verification'/n for n in ['recover_joint_photon_stages.py','recover_exact_joint_stages.py','return_resolved_joint_history.py']]
note=root/'notes/REQUEST204_206_STAGE_HISTORY_GR_SOURCE_KO.md'


def package():
    assert not manifest.exists();a,c,w=folders
    assert read(c/'recovered-64.json')['passed'] and read(c/'coarse-receipt.json')['error'] is None
    assert 'collision_relative' in read(a/'fine-receipt.json')['error'] and 'collision_relative' in read(c/'fine-receipt.json')['error']
    audit=read(c/'proposal-audit.json');norm=read(w/'original-source-norm.json');source=read(w/'sources.json')
    assert audit['original_joint_stage_passed'] and not audit['exact_archival_rate_identity_passed']
    assert not norm['passed'] and not source['passed'] and read(w/'source-receipt.json')['error']
    assert not (w/'coarse-receipt.json').exists() and not (w/'fine-receipt.json').exists()
    assert all(sha(module)==sha(runtime/'verification'/module.name) for module in modules)
    old=read(out.parent/'native-primitive-precision-manifest.json');preserved={p:h for p,h in old['sha256'].items() if not p.startswith('docs/')}
    for name,h in preserved.items():assert sha(root/name)==h,name
    known={h:dict(path=name,sha256=h) for name,h in read(master)['sha256'].items() if not name.startswith('docs/')}
    known.update({r['sha256']:r for r in read(out.parent/'native-primitive-precision/publication.json')['reused'].values()})
    aliases={sha(f):f for folder in folders for f in folder.glob('*.py')};bindings=[]
    for folder in folders:
        for plan in folder.glob('*plan.json'):
            for name,h in read(plan).get('bindings',{}).items():
                actual=runtime/name.removeprefix('/home/lpaiu/work/native-retained-tail-runtime/')
                if sha(actual)!=h:actual=aliases[h]
                assert sha(actual)==h,name;bindings.append(dict(plan=str(plan),source=name,resolved=str(actual),sha256=h))
    out.mkdir();reused={}
    for folder in folders:
        for src in sorted(folder.rglob('*')):
            if not src.is_file():continue
            key=folder.name+'/'+src.relative_to(folder).as_posix();digest=sha(src)
            if digest in known and (root/known[digest]['path']).exists():
                assert sha(root/known[digest]['path'])==digest;reused[key]=known[digest];continue
            dst=out/key;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);assert sha(dst)==digest
            known[digest]=dict(path=dst.relative_to(root).as_posix(),sha256=digest)
    helpers=[(Path(__file__),'publication-producer.py'),(p,'previous-publication-producer.py'),(b.helper,'publication-helper.py'),(root/'.phase205-stage-audit.py','proposal-audit-producer.py'),(root/'.phase206-source-norm.py','original-source-norm-producer.py')]
    for src,name in helpers:shutil.copyfile(src,out/name)
    assert sha(out/'proposal-audit-producer.py')==read(c/'proposal-audit-retry-receipt.json')['source_sha256']
    assert sha(out/'original-source-norm-producer.py')==norm['source_sha256']
    note.write_text('''# 실제 단계 이력의 GR 연결과 초기 시간 오차

분류: Counterexample candidate. **최종 전하는 미판정이다.** 저장된 실제 광자·물질 단계의 GR 원천 연결을 실행했지만 초기 원천의 두 경로 시간 차이가 원2%기준을 넘었다. 새 GR 반환 진화·자기 결합·최종 전하의 해결로 판정하지 않는다. 분류상 loophole progress는 수락 경계와 재사용 가능한 같은 해 이력의 확정까지다.

분류: Counterexample candidate. 204는 저장된 모든 실제 기체 단계에 조건부인 원 광자 Radau 블록만 풀었다. 물질 Newton 이력을 재적분하지 않았다. 최초 sparse LU의float128입력 오류는 실제 단계 이전에 발생했으며, 근사 전처리 입출력만float64로 맞췄다. 실제 연산자와 잔차 정밀도는 유지했다. coarse4단계 끝점은 저장값과1.1414e-16차이로 통과했다. fine여섯째 단계의 작은H교환률 일치 오차1.56354e-12가 추가1e-12기준을 넘었다.

분류: Counterexample candidate. 205는 반올림된 두 시각 차이 대신 원 저장 Radau 가중치로 공통h를 정확히 복원하고, 실제 저장 단계 시각을 사용했다. coarse4단계는87.399초에 끝나고 끝점 상대차1.0922e-17로 통과했다. fine은 같은 여섯째 단계에서H교환률 상대차1.61034e-12로 다시 실패했다. 시각 차이만이 유일 원인은 아니다. 첫5수락 단계와 거부 쌍의 완전한 광자 상태·기체·모멘트·반경 및 각도 출구를 보존했다.

분류: Counterexample candidate. 저장된 fine거부 쌍을 원 전체 광자·4기체 Radau 방정식에 독립적으로 넣었다. 벡터 잔차1.11368e-15<1e-12, 물리 잔차 최대3.83412e-15<1e-13이고 native변화율은 저장값과 정확히 같았다. 작은 교환률 상대 일치와 원 전체 단계식은 서로 다른 검사다. 이 결과로204/205의 엄격한 저장값 일치 실패를 통과로 바꾸지 않았으며, 모멘트 적합도 하지 않았다.

분류: Counterexample candidate. 206은 원188반환 실험의T/64=53.662986마이크로초,coarse2/fine4단계를 그대로 사용했다. 이 범위는205fine실패보다 앞이며, 필요한 모든 광자 복원 단계가 원 잔차·저장 교환률·각도 출구 기준을 통과했다. 더 긴T/32복원 실패를 짧게 줄여 성공이라 부르지 않는다. 실제 저장 기체의 닫는 단계를 원floor규칙으로 읽고, 같은 광자E/Pr/N와Radau가중 반경 출구를 사용했다. 실제 시각의 압력·배경은 기존 두 canonical소유자를 원 시간 가중치로 결합했다.

분류: Counterexample candidate. 원천 수집은38.017초에 끝났다. 압력 독립 탐침 최대1.330e-10<0.002, 압력 대응 최대1.533e-16<1e-12, 같은 물질 수지 최대1.076e-13<1e-8이다. 그러나 기존186공간L1/시간최대 비교에서B6.983765%,광자E7.498528%,광자Pr7.627175%,metric stress6.983761%로 원2%기준을 실패했다. 206최초 비교가 셀최대 norm을 사용한 구현 불일치도 보존하고, 저장 배열만으로 원공간L1 norm을 재검산했다. 양쪽 모두 실패이므로 통과 판정은 없었다. 후기 전체구간 신호로 초기 오차를 다시 정규화하지 않았다.

분류: Conjectural. 초기 원 해의 시간 오차가 작은 GR 반환의 시간 실패에 기여할 수 있다. 지금의 결과는 GR입력의 선형 보간만 고치면 문제가 해결된다는 가설을 지지하지 않는다. 다음 판단은 같은 해의 초기 원천 오차가 반환·최종 전하에 미치는 영향을 추적하고 실제 지배 항을 수정하는 것이다. 격자·기간을 자동 확대하거나 물리 분기를 임의로 고정하지 않는다.

분류: Proven. 기존 특성선 분할의 다항식 자체 검산을 통과했다. 이는 전체 EOS, 균일 미분, 연속 시간 오차 또는 최종 전하의 증명이 아니다.

계산 정책: 실제 반환 예산은coarse900초/fine1500초,원천·장·계량 각600초,6GiB/1스레드로 여유를 뒀다. 원천 실패 뒤 그 하류 계산은 실행하지 않았다. 진행 중인202본 결합 계산은별도coarse7200초/fine10800초 예산과 고정 소스를 유지한다. 본 패키지는202의 진행 중 결과를 종료된 결과로 보존하지 않는다.
''',encoding='utf-8')
    tails={
        'model-definition':'저장된 실제 기체 단계에 조건부인 원 광자 블록으로 누락 이력을 복원했다. 엄격한 검사 통과 앞부분만 원188기간의 실제 GR 원천에 사용했다. 원 긴 복원 실패는 유지한다.',
        'observable-targets':'최종 전하는 미판정이다. 초기 원천의B·광자E·Pr가 원시간 기준을 실패해 새 실제 GR 반환 진화의 수락으로 연결하지 않았다.',
        'adiabatic-limit':'광자 조건부 블록 항등식과 특성선 다항식 검산은 물리 EOS·균일 미분·시간 오차 상계가 아니다. 작은 교환률 상대 일치 검사와 원 전체 결합 단계식의 통과를 구분한다.',
        'nonadiabatic-regime':'원T/64구간의 실제 단계별 원천을 연결했으나 공간L1시간차B6.983765%,광자E7.498528%,광자Pr7.627175%로2%기준을 실패했다. 원 압력 대응·같은 이력 수지는 통과했다.',
        'failure-ledger-dynamic-chi':'204/205fine H저장 교환률 실패를 보존했다. exacth만 고쳐도 해소되지 않았다. 원전체식에 통과한 그 제안으로 복원 실패를 재분류하지 않았다. 206셀최대norm 불일치를 원L1norm으로 독립 재검산해도 원천시간 기준을 실패했다.',
        'dynamic-charge-completion':'최종 전하는 미판정이다. 같은해초기원천의시간오차가 새 GR 반환 적분의 수락을 막는다. 늘린 예산으로 본 결합 속행을 유지하며 자기 GR·동일 해 전하·EOS/미분/공간/경계/정적/관측 요건을 축소하지 않는다.'}
    prefixes={}
    for name,body in tails.items():
        f=root/f'docs/{name}.md';prefixes[f.relative_to(root).as_posix()]=dict(bytes=f.stat().st_size,sha256=sha(f))
        with f.open('ab') as h:h.write(('\n\n## 단계204–206 — 실제 단계 이력과 초기 GR 원천 오차\n\n분류: Counterexample candidate. '+body+' [근거](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',strict_coarse_recovery_passed=True,strict_fine_recovery_passed=False,original_joint_proposal_audit=audit,source=source,original_source_norm=norm,
        GR_return_evolution_executed=False,actual_same_joint_solution_GR_computed=False,original_failures_preserved=True,physical_final_charge_solved=False,self_GR_return_closed=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(out/'final-result.json',final);write(out/'publication.json',dict(reused=reused,previous_nondoc=preserved,document_prefixes=prefixes,plan_bindings=bindings))
    files=modules+[note]+[root/k for k in prefixes]+[f for f in out.rglob('*') if f.is_file()]
    final.update(sha256={f.relative_to(root).as_posix():sha(f) for f in files},document_prefixes=prefixes);write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest);m['native_resolved_return']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
