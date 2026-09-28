"""Preserve actual compensated time evolution and its same-stage anchors."""
from pathlib import Path
from types import FunctionType
import importlib.util,shutil,sys
previous=Path(__file__).with_name('.phase186-publish.py')
if not previous.exists():previous=Path(__file__).with_name('previous-publication-producer.py')
spec=importlib.util.spec_from_file_location('previous',previous);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha=b.read,b.write,b.sha
work=runtime/'native-return-evolution188-work';out=root/'outputs/direct-eos-gr33/native-return-evolution'
manifest=out.parent/'native-return-evolution-manifest.json';module=root/'verification/evolve_same_solution_gr_return.py'
note=root/'notes/REQUEST188_SAME_SOLUTION_GR_RETURN_EVOLUTION_KO.md'


def package():
    assert not manifest.exists() and sha(module)==sha(runtime/'verification'/module.name)
    result=read(work/'result.json');assert read(work/'symbolic.json')['passed']
    for action in ['prepare','coarse','fine']:
        r=read(work/f'{action}-receipt.json');assert r['error'] is None
        assert r['seconds']<read(work/'plan.json')['budgets'][action]
    old=read(out.parent/'native-joint-return-manifest.json')
    preserved={p:h for p,h in old['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    bindings=[]
    for name,h in read(work/'plan.json')['bindings'].items():
        actual=runtime/name.removeprefix('/home/lpaiu/work/native-retained-tail-runtime/')
        assert sha(actual)==h,name;bindings.append(dict(source=name,sha256=h))
    out.mkdir();reused={};inputs=read(work/'plan.json')['reused']
    old_reuse=read(out.parent/'native-joint-return/publication.json')['reused']
    for src in sorted(work.rglob('*')):
        if not src.is_file():continue
        rel=src.relative_to(work);key='native-return-evolution188-work/'+rel.as_posix()
        if key in inputs:
            prev=root/old_reuse[rel.as_posix()]['path'];assert sha(src)==sha(prev)==inputs[key]
            reused[rel.as_posix()]=dict(path=prev.relative_to(root).as_posix(),sha256=sha(src));continue
        dst=out/'completed'/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);assert sha(src)==sha(dst)
    for src,name in [(Path(__file__),'publication-producer.py'),(previous,'previous-publication-producer.py'),(b.helper,'publication-helper.py')]:shutil.copyfile(src,out/name)
    rows=result['rows'];err=max(result['time_relative']);stage=max(r['maximum_true_stage'] for r in rows)
    location=read(work/'baryon-localization.json');assert location['ledger_relative']<1e-8
    moments=max(r['maximum_true_physical_stage'] for r in rows)
    balance=max(max(r['audit']['local_material_balance']) for r in rows)
    port=max(r['audit']['angular_port_relative'] for r in rows)
    anchor=max(max(c['native_relative']) for r in rows for c in r['anchor_checks'])
    costs={a:read(work/f'{a}-receipt.json')['seconds'] for a in ['coarse','fine']}
    verdict='통과했다' if result['passed'] else '실패했다'
    body=f'''# 동일 결합 해를 기준으로 한 GR 반환의 실제 시간 적분

분류: Counterexample candidate. **최종 전하의 기존 결론은 여전히 미판정이다.** 같은 광자·네 물질 변수 해에서 생성한 GR을 반환하는 첫 동시 시간 적분을 실행했다. 원 두 시계의 첫 T/64={result['same_horizon_seconds']*1e6:.12g}마이크로초, 실제2/4단계를 완료했다. 전 기간 또는 자기 GR 반복의 완료가 아니다.

분류: Imported from prior work. 입력은185의 동일 결합 해와187에서 구성한 정확한 중심 질량·실제 각도 출구 lapse다. 원64/128시계, 기존 전선 분할,531셀,EOS·배경·입사 진폭과 단계 방정식을 유지했다. 원래 단계 상태와 반환 증분을 함께 읽으며, 별도 영 상태의 비선형 응답이나 과거 전하에 보정을 더하지 않았다. 실행 중인185와 그 의존 소스·계획은 바꾸지 않았다.

분류: Counterexample candidate. 반환이 입사장보다 작아 직접 덧셈에서 소실되는 문제를 피하도록 상태·입력을 큰 성분과 작은 성분으로 유지했다. 실제 native 연산자의 비교 피연산자를 기록하고 두 성분으로 분기를 선택했다. 큰 성분에서 정확히 같은 값인 분기는 작은 성분이 결정할 수 있다. 큰 성분의 비영 분기 경계를 넘거나8회 안에 분기가 안정되지 않으면 거절한다. 해당 경로의 단계 검사가 다른 모든 상태의 분기 안정성이나 균일 EOS 미분을 증명하지는 않는다.

분류: Proven. 같은 affine 분기에서 F(z+d)-F(z)=Jd라는 항등식과, 합하면 사라지는 작은 성분이 같은 값 사이의 min/argmin을 결정하는 산술 예를 검사했다. 이는 독립 비선형 해의 전역 중첩 정리가 아니다.

분류: Counterexample candidate. 저장된 실제 Eref/B/S/H 단계로 복원한 원천과 원 저장 원천의 상대 차이는 최대{anchor:.9g}였다. 반환 증분의 원 단계식 잔차 최대{stage:.9g}<1e-12, 네 물리 모멘트 잔차 최대{moments:.9g}<1e-13, 같은 증분 해의 물질 수지 최대{balance:.9g}<1e-8, 실제 가중치 광자 출구 차이 최대{port:.9g}<1e-12였다. 원천·충돌·floor·출구는 이 실제 단계들에서 함께 저장했다.

분류: Counterexample candidate. 여섯 광자/열·네 물질량의 시간 비교는 최대{100*err:.9g}%이며 원2%기준을 {verdict}. 이 결과는 작은 반환 증분의 시간 대조다. 큰 원 해의 이산화·저장 오차가 이 작은 증분보다 작다는 주장이나 전체 전하의 오차 상계로 바꾸지 않는다. 전체 입력 대비 작은 오차만으로 반환 해상도를 인정하지 않는다.

분류: Counterexample candidate. 실패한 성분은 바리온이다. 저장 이력에서 그 끝점 차이는 native 유속 적분 차이로 설명되며, 두 수지의 차이는 상대{location['ledger_relative']:.9g}다. 바닥 처리 차이의 L1은 끝점 차이의 {location['floor_difference_over_final_difference']:.9g}에 불과하다. 셀144–147이 차이의 약97.6%를 차지한다. 아직 유일한 시간 오차 원인을 확정한 것은 아니다. 처음 진단에서는 물리 입사 진폭과 계산 정규화 AMP를 혼동했으며 원 AMP를 가져와 바로잡았다. 잘못된 진단과 코드는 rejected 이름으로 보존한다.

분류: Counterexample candidate. 실제 벽시간은 첫 경로{costs['coarse']:.3f}초, 둘째 경로{costs['fine']:.3f}초로 등록한150/250초 상한 안이다. CPU1스레드·가상 메모리4GiB의 경로별 제한을 유지했고, 자동 시간격자·기간·GR반복 확대는 실행하지 않았다.

분류: Conjectural. 다음은 이 시간 실패의 원인을 수정한 실제 반환 적분이다. 실패한 해의 기간을 늘려 상대 차이를 희석하거나 격자를 자동 세분화하지 않는다. 먼저 세 절점 GR의 시간 표현과 원 상태 분기의 시간 변화가 초기 native 유속에 주는 영향을 구분하고, 그 결과로 필요한 수정과 계산 예산을 정해야 한다. 이후 같은 해의 에너지·경계·생성 GR에서 반복 잔차와 최종 전하를 판정한다. 전체 기간·자기 결합·EOS 균일 미분·공간·경계·완전 비선형·정적 비교·관측 폐쇄는 계속 남는다.

분류: Counterexample candidate. [판정](../outputs/direct-eos-gr33/native-return-evolution/final-result.json), [실행 코드](../verification/evolve_same_solution_gr_return.py), [결합 기록](../outputs/direct-eos-gr33/native-return-evolution-manifest.json)에 근거를 보존한다. 실제 반환 시간 적분까지 연결했지만 원 시간 기준 실패로 확장을 중단했다. 전체 목표는 미완료다.
'''
    note.write_bytes(body.encode())
    tails={
        'model-definition':'같은 저장 단계의 원 상태와 작은 반환 증분을 함께 사용하는 동시 Radau 방정식을 실제 적분했다. native 분기는 두 성분으로 결정하며 독립 비선형 응답의 중첩으로 바꾸지 않는다.',
        'observable-targets':'최종 전하는 미판정이다. 첫T/64의 실제 GR 반환 증분과 자체 에너지·경계 이력이 생겼지만 전체 기간·자기 결합·동일 해의 최종 전하 판독은 남는다.',
        'adiabatic-limit':'같은 affine 분기의 증분 항등식과 작은 성분의 정확한 동률 선택을 검사했다. 비영 큰 성분 경계의 교차는 거절한다. 저장 단계 검사는 전역 선형성이나 균일 EOS 미분 정리가 아니다.',
        'nonadiabatic-regime':f'같은T/64에서 실제2/4단계 GR 반환 적분을 완료했다. 열 가지 시간 대조 최대{100*err:.9g}퍼센트로 원2퍼센트 기준을 {verdict}. 현재 계량은 같은 원 해에서 생성한 첫 고정 반환이며 자기 GR 반복은 미완료다.',
        'failure-ledger-dynamic-chi':'반환 적분은 실행됐지만 바리온 시간 차이5.211566퍼센트로 원2퍼센트 기준을 실패했다. 차이는 native 유속 누적에 있고 floor는 지배적이지 않다. 원인 수정 전 기간·격자 확대를 하지 않는다. 직접 덧셈 소실 실패와 원 해/작은 증분의 오차 경계를 유지한다.',
        'dynamic-charge-completion':'GR 반환 성분을 같은 원 단계 상태에 연결한 실제 동시 적분2/4단계를 완료했다. 최종 전하 결론은 미판정이다. 전 기간·계량 시간 표현·자기 결합 잔차와 같은 해 전하·물리 폐쇄가 계속 필요하다.'}
    prefixes={}
    for name,text in tails.items():
        f=root/f'docs/{name}.md';prefixes[f.relative_to(root).as_posix()]=dict(bytes=f.stat().st_size,sha256=sha(f))
        with f.open('ab') as h:h.write(('\n\n## 단계188 — 같은 단계 해의 GR 반환 시간 적분\n\n분류: Counterexample candidate. '+text+' [근거](../notes/'+note.name+').\n').encode())
    final=dict(result,actual_same_joint_solution_GR_computed=True,GR_return_time_evolved=True,
        full_horizon_completed=False,physical_final_charge_solved=False,cost_seconds=costs)
    write(out/'final-result.json',final)
    write(out/'publication.json',dict(reused=reused,previous_nondoc=preserved,document_prefixes=prefixes,plan_bindings=bindings))
    files=[module,note]+[root/k for k in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(sha256={p.relative_to(root).as_posix():sha(p) for p in files},document_prefixes=prefixes)
    write(manifest,final);m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_same_solution_GR_return_evolution']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
