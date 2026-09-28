"""Publish same-solution exterior readout and the actual full-return status."""
from pathlib import Path
from types import FunctionType
import datetime,importlib.util,sys
spec=importlib.util.spec_from_file_location('b',Path('.phase251-publish.py'));b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha,copy=b.read,b.write,b.sha,b.copy
out=root/'outputs/direct-eos-gr33/native-returned-exterior'
manifest=out.parent/'native-returned-exterior-manifest.json'
note=root/'notes/REQUEST252_SAME_RETURN_EXTERIOR_KO.md'


def package():
    w=runtime/'native-returned-exterior252-work';r=read(w/'common/result.json');audit=read(w/'common/audit.json')
    assert r['passed'] and audit['passed'] and read(w/'common-receipt.json')['error'] is None
    modules=[root/'verification'/p for p in ['read_returned_exterior.py','verify_returned_exterior.py']]
    assert sha(modules[0])==read(w/'common-receipt.json')['source_sha256']
    for p,h in read(w/'common/plan.json')['bindings'].items():assert sha(runtime/p)==h,p
    for p,h in audit['bindings'].items():assert sha(runtime/p.removeprefix('/home/lpaiu/work/native-retained-tail-runtime/'))==h,p
    assert not out.exists();out.mkdir()
    preserved={p:h for p,h in read(b.manifest)['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    for p in (w/'common').iterdir():
        if p.is_file():copy(p,out/'common'/p.name)
    for name in ['check.json','check-receipt.json','common-start.json','common-receipt.json','full-controller-start.json','full-controller-status.json']:copy(w/name,out/name)
    for name in ['.phase252-common.stdout.log','.phase252-common.stderr.log','.phase252-followthrough.py']:copy(root/name,out/name.lstrip('.'))
    copy(Path(__file__),out/'publication-producer.py')
    returned=runtime/'native-full-return249-work';metric=read(returned/'full/metric-result.json')
    assert metric['passed'] and read(returned/'metric-receipt.json')['error'] is None
    for name in ['controller-start.json','controller-status.json','metric-receipt.json','full/metric-result.json','full/capture-64.json']:
        if (returned/name).exists():copy(returned/name,out/'full-return'/name)
    copy(runtime/'native-full-charge251-work/full/controller-status.json',out/'full-charge-status.json')
    now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    fine=next(v for v in r['components'] if v['clock']==128)
    note.write_text(f'''# 같은 반환 해의 외부·질량 판독

분류: Counterexample candidate. **최종 물리 전하는 미판정이다.** 공통15/16기간의 같은 실제 물질·광자·GR 반환 해에 고정된 초기 외부 계량의 연산자를 적용했다. 기존 compact 전하의 음의 부호는 외부 scalar, 도착 에너지, 같은 원천의 homogeneous 질량을 포함한 조건부 판독에서도 유지됐다. 전체 기간249반환과251전하 판독 뒤 이 연산자를 자동 적용하도록 연결했다. 게시={now}.

분류: Counterexample candidate. fine의 compact 높은 성분은−1.368965112963e−51, 외부 scalar는+6.041864583481e−66, 도착 에너지는−5.933275462481e−4erg다. 음의 도착 에너지는 기준에 대한 방출 감소다. 같은 원천의 homogeneous 질량3.303079699849e−49cm를0으로 가정하거나 방출에 흡수하지 않았다. 이 질량은 실제 적용했던 GR 계량의 질량 제약과 정확히 일치했다. 초기 질량 기준 정규화의 높은 성분은{fine['high']}이고 compact 대비 변화는0.00330957%다. 동일 해의 낮은 성분은{fine['same_solution_low_increment']}로 별도 보존했다. 작은 성분의 질량 정규화는 그 compact 값의 약46.6%를 바꾸므로 그 항을 단순 생략하지 않는다.

분류: Proven. 원 실제 Radau 두 순간의 luminosity를 잇는 일차 다항식의 적분이 원3/4·1/4가중치와 같음을 기호 검사했다. 첫째·둘째 원시함수를 각 구간 직접 적분과 대조했고, 0 입력 및 signed 입력을 검사했다. 유리식 전하 차분도 독립160자리 산술과 대조했다. 이들은 저장 표현의 항등식·산술 검사이며 물리 정확도 인증이 아니다.

분류: Counterexample candidate. 실제 각도 출구와 원천 끝점의 차이는9.36e−16이하<1e−12다. 외부 scalar의 각도 구적 최대0.05498%, 반경 구적0.0001047%, 도착 에너지까지 포함한 시간 대조 최대0.005585%로 원0.2%/2%기준을 통과했다. 질량·정규화 성분을 따로 검사해 낮은 질량의 시간 차이1.24491%, 낮은 정규화 전하0.387767%, 높은 정규화 전하0.00515001%도 원2%기준으로 수락했다. 전체 합이 큰 성분에 가려 작은 성분의 실패를 숨기지 않도록 성분별로 판정했다.

분류: Conjectural. 이 판독은 초기 고정 광선과 기준 에너지 방출에 조건부다. 시간 의존 계량의 물리 좌표 에너지 변환·외부 광선의 일과 기존 배경 질량 재정규화는 미완료다.166에서 이미 구분한 기준 에너지와−p_t를 혼동하지 않으며, 본 결과를 전체 물리 무한대 전하나 완전 비선형 진화로 표시하지 않는다. EOS·균일 미분/시간·공간/경계·자기GR·정적 EFT 비교·관측 연결의 전체 완료 조건과188초기 국소 실패도 유지한다.

분류: Counterexample candidate. 실제 전체 기간의249계량 계산은{read(returned/'metric-receipt.json')['seconds']:.3f}초에 완료됐고, 원 metric 시간·구적·광선 및 과거 구간 접합 기준을 통과했다. 그 입력을 실제 결합 coarse 진화에 적용하기 시작했다. 게시 스냅샷의 단계 수는 진행 상황이며 전체119/231반환 쌍의 완료를 뜻하지 않는다. 앞서 완료된 원 해를 다시 적분하지 않고 마지막 두 구간만 이어 푼다.

실행 기록: 이번 공통 기간 외부·질량 판독 본문{read(w/'common-receipt.json')['seconds']:.3f}초, 최대RSS{read(w/'common-receipt.json')['peak_RSS_bytes']}바이트, 새 물리 단계0. 기존4/8구적과17개 외부 출력 시각을 사용했다.30분·16GiB상한과 기존 수락 기준을 유지하며 전체 기간 후속도 같은 범위다. 분류상 조건부 loophole progress이며 전체 목표는 진행 중이다.
''',encoding='utf-8')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        p=root/f'docs/{name}.md';prefixes[p.relative_to(root).as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
        with p.open('ab') as f:f.write(('\n\n## 단계252 — 같은 반환 해의 외부·질량 판독\n\n분류: Counterexample candidate. 최종 물리 전하 미판정. 긴 같은 반환 해의 실제 Radau 방출과 동일 원천 질량을 고정 외부 연산자로 읽어 음의 부호를 유지했다. 작은 낮은 성분의 질량 정규화 변화도 따로 적용하고 시간·구적을 성분별로 판정했다. 물리 에너지 변환·시간 의존 외부·배경 질량 재정규화와 전체 EOS/미분/공간/비선형/관측 범위는 남는다. 전체249계량은 완료되어 실제 마지막 두 구간의 결합 진화에 적용됐으며,251후속에 같은 외부 판독을 연결했다. [근거](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',actual_same_joint_solution_GR_computed=True,
        conditional_same_return_exterior_admitted=True,components=r['components'],audit_controls=audit['controls'],
        full_period_metric_applied_to_actual_return=True,full_period_actual_return_completed=False,
        final_charge_conclusion='unadjudicated',physical_final_charge_solved=False,full_goal_complete=False,snapshot_KST=now)
    write(out/'result.json',final);write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes,reused=read(w/'common/plan.json')['bindings']))
    files=modules+[note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes,sha256={p.relative_to(root).as_posix():sha(p) for p in files});write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_returned_exterior']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
