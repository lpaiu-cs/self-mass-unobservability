"""Bind the native precision repair and its actual terminal evolution together."""
from pathlib import Path
from types import FunctionType
import importlib.util,shutil,sys
previous=Path(__file__).with_name('.phase186-publish.py')
if not previous.exists():previous=Path(__file__).with_name('previous-publication-producer.py')
spec=importlib.util.spec_from_file_location('previous',previous);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha=b.read,b.write,b.sha
out=root/'outputs/direct-eos-gr33/native-primitive-precision';manifest=out.parent/'native-primitive-precision-manifest.json'
folders=[runtime/n for n in ['native-affine-precision198-work','native-precise-continuation199-work']]
modules=[root/'verification'/n for n in ['repair_native_affine_precision.py','continue_precise_native.py']]
note=root/'notes/REQUEST198_199_NATIVE_PRECISION_CONTINUATION_KO.md'


def package():
    assert not manifest.exists()
    control=read(folders[-1]/'controller-status.json');assert control['state'] in ['completed','failed']
    a,w=folders
    assert read(a/'primitive_check-receipt.json')['error'] is None
    assert (w/'coarse-receipt.json').exists()
    for p in modules:assert sha(p)==sha(runtime/'verification'/p.name)
    old=read(out.parent/'native-roundoff-continuation-manifest.json')
    preserved={p:h for p,h in old['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    known={h:dict(path=p,sha256=h) for p,h in preserved.items()}
    known.update({r['sha256']:r for r in read(out.parent/'native-roundoff-continuation/publication.json')['reused'].values()})
    aliases={sha(p):p for folder in folders for p in folder.glob('*producer.py')};bindings=[]
    for folder in folders:
        for plan in folder.glob('*plan.json'):
            for name,h in read(plan).get('bindings',{}).items():
                p=runtime/name.removeprefix('/home/lpaiu/work/native-retained-tail-runtime/')
                if sha(p)!=h:p=aliases[h]
                assert sha(p)==h,name;bindings.append(dict(plan=str(plan),source=name,resolved=str(p),sha256=h))
    out.mkdir();reused={}
    for folder in folders:
        for src in sorted(folder.rglob('*')):
            if not src.is_file():continue
            key=folder.name+'/'+src.relative_to(folder).as_posix();digest=sha(src)
            if digest in known:reused[key]=known[digest];continue
            dst=out/key;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);assert sha(dst)==digest
            known[digest]=dict(path=dst.relative_to(root).as_posix(),sha256=digest)
    for src,name in [(Path(__file__),'publication-producer.py'),(previous,'previous-publication-producer.py'),(b.helper,'publication-helper.py'),(root/'.phase199-followthrough.py','followthrough-producer.py')]:shutil.copyfile(src,out/name)
    result=read(w/'result.json') if (w/'result.json').exists() else None
    failure=read(w/'failure-64.json') if (w/'failure-64.json').exists() else None
    coarse=read(w/'path-64.json') if (w/'path-64.json').exists() else None
    progress=read(w/'stage-progress-64.json')
    passed=bool(result and result['passed'] and control['state']=='completed')
    actual=coarse['audit']['actual_completed_steps'] if coarse else failure['actual_accepted_steps'] if failure else progress['solved_actual_steps']
    event=('두 실제 경로의 전체 기간과 원 시간 대조를 통과했다.' if passed else
           f'속행은 {control["action"]}에서 종료됐고 coarse 수락 단계는 {actual}개다. 전체 기간 및 최종 전하를 완료로 판정하지 않았다.')
    note.write_text(f'''# 물질 복원의 정밀도 수정과 실제 결합 진화

분류: Counterexample candidate. **최종 전하는 미판정이다.** {event} 같은 해의 자기 GR 반환·무한대 전하·EOS/균일 미분/공간/경계/비선형/관측 요건을 유지한다. 다른 해의 진단을 가산하거나 기존 전하 부호를 상속하지 않는다.

분류: Counterexample candidate. 198에서 197의 마지막 실패 쌍을 정확히 재현했다. affine 우변의 산술 수정 효과는1.46816e-13으로 실제 native/선택 affine 차이5.18366e-11을 설명하지 못했다. 이 가설과 미통과 결과를 보존하고 affine 전용 변경을 실제 해법으로 채택하지 않았다.

분류: Counterexample candidate. 실제 native 물질 복원에서 다섯 기본 float 배열과 다섯 float 변환이 확장 정밀도 상태를 binary64로 내리고 있었다. 실제 및 선택 분기의 같은 복원 함수에서 이 변환만 확장 정밀도로 유지했다. EOS 은행·광자 계수 맵·물리식·분기 규칙은 유지했다. 저장 쌍의 native/선택 affine 차이가7.56658e-13으로 감소했고 원 구성식 탐침은 통과했으나, 저장된 옛 제안의 실제 단계 잔차2.69905e-10은 여전히 실패였다. 진단 성공으로 그 상태를 수락하지 않았다.

분류: Counterexample candidate. 199는 기존 두 경로의 저장된656개 물질 Radau 단계를 수정 복원으로 다시 읽었다. native 변화는 최대5.2081e-12, 동일 이력 물질 수지는 최대1.3314e-13으로 원 구성식0.002·수지1e-8 기준을 통과했다. 수락 이력을 재적분하지 않았다. 이 검사는 균일 오차 상계나 과거 모든 국소 벡터 방정식의 새 인증이 아니다.

분류: Counterexample candidate. 수정된 실제 방정식으로 마지막 거부 쌍에서 다시 풀어, 기존 실패 단계를 원 native 벡터9.6130523e-13<1e-12와 물리 잔차1.9767e-17이하로 수락했다. 실제 수락 이력이113에서114단계로 진전됐다. 그 뒤 같은 상태의 남은 구간을 속행한 최종 결과와 원 실패 배열은 동봉된 생산 결과가 소유한다. 최종 coarse 단계 수는{actual}개이며, 전체 두 경로 통과 여부는{passed}다.

계산 정책: 저장 상태와 유효 제안을 재사용하고 coarse30분·fine45분, CPU1스레드·가상 메모리6GiB, 새 Newton8회·선형 보정12회를 허용했다. 정확도·기간·격자는 유지했다. 197의 실패한8회 제안과199의 수정 산술 제안을 구분한다. 실행 중 생산 코드·계획은 수정하지 않았다.
''',encoding='utf-8')
    tails={
        'model-definition':'native 물질 복원의 중간 배열·명시 변환에서 확장 정밀도를 유지하고 실제 결합 진화에 적용했다. EOS 은행·광자 계수·물리식·분기는 그대로다.',
        'observable-targets':'최종 전하는 미판정이다. '+event+' 같은 해의 자기 GR·무한대 전하와 지배 오차를 함께 판정해야 한다.',
        'adiabatic-limit':'저장된656개 물질 단계의 수정 구성식·동일 이력 수지를 확인했다. 과거 모든 국소 벡터 방정식의 새 인증이나 균일 EOS 미분 보장으로 해석하지 않는다.',
        'nonadiabatic-regime':'기존 실패 실제 단계를 원 native 잔차9.6130523e-13으로 수락해114단계에 도달했고 동일 해를 속행했다. '+event,
        'failure-ledger-dynamic-chi':'affine 우변 산술이 지배 원인이라는 가설은 대조로 기각했다. 정밀도 수정만 적용한 옛 쌍도 원 단계 잔차를 실패하여 실제 재풀이 전에는 수락하지 않았다. 이후 속행의 실패와 수락 기준을 보존한다.',
        'dynamic-charge-completion':'최종 전하 결론은 미판정이다. '+event+' 수락된 원천·단계·기간을 최종 전하 완료로 대체하지 않는다.'}
    prefixes={}
    for name,body in tails.items():
        p=root/f'docs/{name}.md';prefixes[p.relative_to(root).as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
        with p.open('ab') as f:f.write(('\n\n## 단계198–199 — native 정밀도와 실제 진화\n\n분류: Counterexample candidate. '+body+' [근거](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',controller=control,actual_coarse_accepted_steps=actual,
        full_horizon_completed=passed,original_failures_preserved=True,native_precision_repair_applied_to_actual_evolution=True,
        actual_same_joint_solution_GR_computed=False,physical_final_charge_solved=False,
        final_charge_conclusion='unadjudicated',full_goal_complete=False,actual_result=result,actual_failure=failure)
    write(out/'final-result.json',final);write(out/'publication.json',dict(reused=reused,previous_nondoc=preserved,document_prefixes=prefixes,plan_bindings=bindings))
    files=modules+[note]+[root/k for k in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(sha256={p.relative_to(root).as_posix():sha(p) for p in files},document_prefixes=prefixes);write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_primitive_precision']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
