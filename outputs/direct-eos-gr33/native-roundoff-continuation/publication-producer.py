"""Publish terminal actual continuations, preserving their original failures."""
from pathlib import Path
from types import FunctionType
import importlib.util,shutil,sys
previous=Path(__file__).with_name('.phase186-publish.py')
if not previous.exists():previous=Path(__file__).with_name('previous-publication-producer.py')
spec=importlib.util.spec_from_file_location('previous',previous);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha=b.read,b.write,b.sha
out=root/'outputs/direct-eos-gr33/native-roundoff-continuation';manifest=out.parent/'native-roundoff-continuation-manifest.json'
folders=[runtime/n for n in ['native-resumed-newton195-work','native-balanced-krylov196-work','native-roundoff-polish197-work']]
modules=[root/'verification'/n for n in ['balance_full_native_krylov.py','polish_full_native_roundoff.py']]
note=root/'notes/REQUEST195_197_ACTUAL_NATIVE_CONTINUATION_KO.md'


def package():
    assert not manifest.exists()
    terminal=[read(w/'controller-status.json') for w in folders]
    assert all(s['state'] in ['completed','failed'] for s in terminal)
    for p in modules:assert sha(p)==sha(runtime/'verification'/p.name)
    a,c,w=folders
    previous_manifest=read(out.parent/'native-newton-resume-admission-manifest.json')
    preserved={p:h for p,h in previous_manifest['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    known={h:dict(path=p,sha256=h) for p,h in preserved.items()}
    known.update({r['sha256']:r for r in read(out.parent/'native-newton-resume-admission/publication.json')['reused'].values()})
    aliases={sha(p):p for folder in folders for p in folder.glob('*producer.py')};bindings=[]
    for folder in folders:
        for name,h in read(folder/'plan.json')['bindings'].items():
            p=runtime/name.removeprefix('/home/lpaiu/work/native-retained-tail-runtime/')
            if sha(p)!=h:p=aliases[h]
            assert sha(p)==h,name;bindings.append(dict(source=name,resolved=str(p),sha256=h))
    out.mkdir();reused={}
    for folder in folders:
        for src in sorted(folder.rglob('*')):
            if not src.is_file():continue
            key=folder.name+'/'+src.relative_to(folder).as_posix();digest=sha(src)
            if digest in known:reused[key]=known[digest];continue
            dst=out/key;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);assert sha(dst)==digest
            known[digest]=dict(path=dst.relative_to(root).as_posix(),sha256=digest)
    for src,name in [(Path(__file__),'publication-producer.py'),(previous,'previous-publication-producer.py'),(b.helper,'publication-helper.py'),(root/'.phase197-followthrough.py','followthrough-producer.py'),(root/'.phase197-inspect.py','saved-array-inspection.py')]:shutil.copyfile(src,out/name)
    result=read(w/'result.json') if (w/'result.json').exists() else None
    failure=read(w/'failure-64.json') if (w/'failure-64.json').exists() else None
    coarse=read(w/'path-64.json') if (w/'path-64.json').exists() else None
    passed=bool(result and result['passed'] and terminal[-1]['state']=='completed')
    actual=coarse['audit']['actual_completed_steps'] if coarse else failure['actual_accepted_steps'] if failure else 113
    event=('두 실제 경로의 전체 기간과 원 시간 대조를 통과했다.' if passed else
           f'실제 속행은 {terminal[-1]["state"]} 상태이며, coarse 수락 단계는 {actual}개다. 원 수락 기준을 넘긴 결과를 완료로 승격하지 않았다.')
    note.write_text(f'''# 저장된 실제 결합 해의 선형 정밀도와 속행

분류: Counterexample candidate. **최종 전하는 미판정이다.** {event} 같은 해의 자기 GR 반환, 무한대 전하, EOS/균일 미분/공간/경계/비선형/관측 폐쇄는 여전히 별도 요건이다. 기존 전하 부호를 상속하거나 다른 해의 진단을 가산하지 않는다.

분류: Counterexample candidate. 195에서 기존 실패 native 단계를 원 단계 잔차5.4311478e-13<1e-12로 수락했고 총113개 실제 하위 단계에 도달했다. 다음 단계는4회 선형 보정 뒤 벡터1.56758e-9, 물리 모멘트0.11624/0.01502로 실패했다. 시간 초과가 아니었으며 마지막 수락 상태와 실패 반복 전체를 저장했다.

분류: Counterexample candidate. 196은 동일 RHS·분기 제안을 정확히 확인하고 그 반복값을 재사용했다. 가역 대각 스케일과 오른쪽 전처리로 전체 벡터 및 네 물리 성분을 함께 반영하고, Krylov restart80/maxiter10 및12회 선형 보정을 적용했다. 실제 계산{read(c/'coarse-receipt.json')['seconds']:.3f}초 뒤 물리 오차는1.81e-19이하였으나 벡터7.2300250e-14>1e-14로 실패했다. 새 물리 단계는 수락하지 않았다.

분류: Counterexample candidate. 저장된196배열에서 잔차 노름3.7854e-5는 거의 모두 셀261의 두 B 방정식에 있었다. 해당 B 값은 약2.1864e15와-5.5040e15이고 최소 표현 간격은1.2207e-4와4.8828e-4였다. 작은 보정이 반올림으로 사라질 수 있음을 보여주지만, 이것만으로 원 기준이 달성 불가능하다고 증명한 것은 아니다.

분류: Counterexample candidate. 197은 같은 전체 연산자에 연결된 표현 정밀도가 더 나은 물질 열을 최대24개 사용해 잔차 최소화 보정을 적용했다. 갱신한 실제 배열의 전체 벡터·네 물리 잔차를 다시 계산하며, 그 뒤 원 native 비선형 단계·구성식·보존·출구·2퍼센트 시간 기준이 최종 수락을 결정한다. 이전195와196의 마지막 수락 상태 전체가 정확히 동일함을 확인해 수락된113단계를 재계산하지 않았다.

분류: Proven. 가역 스케일 변환과 선택 열 보정은 원 선형 방정식을 보존하는 대수 항등식이다. 해당 검사는 수렴이나 물리 전하의 증명이 아니다.

자원 정책: 사용자 승인에 따라 coarse30분·fine45분, CPU1스레드·가상 메모리6GiB, 단계별8회 Newton과12회 선형 보정을 허용한다. 정확도·기간·격자는 유지했다. 새 실행을 반복 중단해서 수락 이력과 유효 반복을 버리는 비용을 피한다. 종료 receipt·실패 배열·원 생산 코드와 계획은 동봉된 해시로 고정한다.
''',encoding='utf-8')
    tails={
        'model-definition':'같은 결합 해의113수락 단계를 재사용했다. 물리 성분별 Krylov 스케일과 표현 가능한 물질 좌표의 잔차 보정을 실제 native 풀이에 적용했다. 방정식·물리 판정 기준은 유지했다.',
        'observable-targets':'최종 전하는 미판정이다. '+event+' 반환 GR과 같은 해의 전하 판독까지 통과해야 전하 결론을 내릴 수 있다.',
        'adiabatic-limit':'가역 스케일·선택 열 보정의 대수 동일성을 확인했다. 동일성은 동역학 수렴·미분 오차 상계·물리 전하의 증명이 아니다.',
        'nonadiabatic-regime':'원 시간격자와 기간을 유지한 실제 속행에30/45분 및8회 Newton·12회 선형 반복을 배정했다. '+event,
        'failure-ledger-dynamic-chi':'195의 다음 단계 물리·선형 실패와196의7.230025e-14벡터 잔차 실패를 원 배열·코드·receipt로 보존했다. 작은 물리 모멘트 잔차로 벡터 기준을 면제하지 않았다.',
        'dynamic-charge-completion':'최종 전하 결론은 미판정이다. '+event+' 완료 기준은 지배 오차를 해결한 동일 결합 해에서 자기 GR·최종 전하와 나머지 물리 오차를 함께 판정하는 것이다.'}
    prefixes={}
    for name,body in tails.items():
        p=root/f'docs/{name}.md';prefixes[p.relative_to(root).as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
        with p.open('ab') as f:f.write(('\n\n## 단계195–197 — 실제 잔차 보정과 결합 속행\n\n분류: Counterexample candidate. '+body+' [근거](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',controller=terminal[-1],
        actual_coarse_accepted_steps=actual,full_horizon_completed=passed,original_failures_preserved=True,
        actual_same_joint_solution_GR_computed=False,physical_final_charge_solved=False,
        final_charge_conclusion='unadjudicated',full_goal_complete=False,
        actual_result=result,actual_failure=failure)
    write(out/'final-result.json',final);write(out/'publication.json',dict(reused=reused,previous_nondoc=preserved,document_prefixes=prefixes,plan_bindings=bindings))
    files=modules+[note]+[root/k for k in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(sha256={p.relative_to(root).as_posix():sha(p) for p in files},document_prefixes=prefixes);write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_roundoff_continuation']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
