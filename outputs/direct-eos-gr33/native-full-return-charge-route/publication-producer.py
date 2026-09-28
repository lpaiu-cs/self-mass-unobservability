"""Preserve the same-return charge result and the resumed full-period chain."""
from pathlib import Path
from types import FunctionType
import datetime,importlib.util,sys
spec=importlib.util.spec_from_file_location('b',Path('.phase250-publish.py'));b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha,copy=b.read,b.write,b.sha,b.copy
out=root/'outputs/direct-eos-gr33/native-full-return-charge-route'
manifest=out.parent/'native-full-return-charge-route-manifest.json'
note=root/'notes/REQUEST251_FULL_RETURN_CHARGE_ROUTE_KO.md'


def package():
    w=runtime/'native-full-charge251-work';extension=runtime/'native-retarded-extension248-work'
    charge=runtime/'native-compensated-charge247-work/full';r=read(charge/'result.json')
    assert read(charge.parent/'full-controller-status.json')['state']=='completed'
    reg=read(w/'regression.json');assert reg['passed'] and read(w/'check-retry/controller-status.json')['state']=='completed'
    assert read(extension/'prefix-check.json')['passed']
    assert read(extension/'controller-status.json')['state']=='completed' and read(extension/'result.json')['GR_return_admitted']
    assert read(w/'full/controller-status.json')['state'] in ['waiting','running','completed']
    modules=[root/'verification'/v for v in ['read_full_return_charge.py','extend_unchanged_source_prefix.py','align_return_restart_clocks.py']]
    assert read(w/'regression-receipt.json')['source_sha256']==sha(modules[0])
    assert read(extension/'prefix-check-receipt.json')['source_sha256']==sha(modules[1])
    assert read(runtime/'native-full-return249-work/clock-check-receipt.json')['source_sha256']==sha(modules[2])
    assert not out.exists();out.mkdir();preserved={p:h for p,h in read(b.manifest)['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    for folder,label in [(w/'extension-failure','extension-failure'),(w/'check-initialization-failure','initialization-failure')]:
        for p in folder.rglob('*'):
            if p.is_file():copy(p,out/label/p.relative_to(folder))
    for p in (w/'interruption-0605').rglob('*'):
        if p.is_file() and p.suffix in ['.json','.log','.py']:copy(p,out/'interruption-0605'/p.relative_to(w/'interruption-0605'))
    for p in (w/'clock-failure').rglob('*'):
        if p.is_file() and p.suffix in ['.json','.log','.py']:copy(p,out/'clock-failure'/p.relative_to(w/'clock-failure'))
    for name in ['clock-alignment.json','clock-check-receipt.json']:
        copy(runtime/'native-full-return249-work'/name,out/name)
    for folder,label in [(w/'check','failed-check'),(w/'check-retry','check')]:
        for p in folder.rglob('*'):
            if p.is_file() and p.suffix in ['.json','.log']:copy(p,out/label/p.relative_to(folder))
    for name in ['regression.json','regression-receipt.json','source-prefix-difference.json']:copy(w/name,out/name)
    for name in ['prefix-check.json','prefix-check-receipt.json','plan.json','prepare-receipt.json','dispatch-forecast.json','result.json','audit-receipt.json','field1288-receipt.json','field648-receipt.json','field1284-receipt.json']:
        copy(extension/name,out/'extension'/name)
    for p in (extension/'gr').glob('fields-*.json'):copy(p,out/'extension/gr'/p.name)
    for folder,label in [(runtime/'native-full-captured244-work','primary-source'),(runtime/'native-returned-source245-work/full','returned-endpoint'),(runtime/'native-dense-returned246-work/full','returned-dense'),(charge,'returned-charge')]:
        for p in folder.iterdir():
            if p.is_file() and p.suffix=='.json':copy(p,out/label/p.name)
        if label=='returned-charge':
            for p in (folder/'gr').glob('fields-*.json'):copy(p,out/label/'gr'/p.name)
    for name in ['.phase251-followthrough.py','.phase251-extension-followthrough.py','.phase251-rearm-extension.py','.phase251-source-prefix.py','.phase251-reader-launch.py','.phase251-resume.py','.phase251-resume-return.py','.phase251-aligned-return.py','.phase251-rearm-clocks.py','.phase251-retire-partial-seed.py']:
        copy(root/name,out/name.removeprefix('.'))
    copy(Path(__file__),out/'publication-producer.py')
    bindings={}
    for folder in [runtime/'native-full-captured244-work',runtime/'native-dense-returned246-work/full',charge,runtime/'native-complete-return236-work',extension]:
        for p in (folder/'gr').glob('*.npz'):bindings[str(p.relative_to(runtime))]=sha(p)
    for row in reg['rows']:
        for key in ['file','reference']:
            p=runtime/row[key];bindings[row[key]]=sha(p)
    snapshots={}
    for folder,status in [(runtime/'native-compensated-charge247-work','full-controller-status.json'),(extension,'controller-status.json'),(runtime/'native-full-return249-work','controller-status.json'),(w/'full','controller-status.json')]:
        label=folder.relative_to(runtime).as_posix().replace('/','-');value=read(folder/status);snapshots[label]=value;write(out/(label+'-status.json'),value)
        start='full-controller-start.json' if '247' in label else 'controller-start.json';copy(folder/start,out/(label+'-start.json'))
    fine=next(v for v in r['components'] if v['clock']==128)['values']['endpoint_compact_with_metric']
    now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    note.write_text(f'''# 같은 반환 해의 긴 전하 판독과 전체 기간 연결

분류: Counterexample candidate. **최종 물리 전하는 여전히 미판정이다.** 실제236의 공통15/16기간111/215단계 해를 그 자체의 물질·광자·반경 출구·적용 계량으로 읽어 GR장과 compact 전하를 계산했다. 전하 대조 수락={r['charge_comparison_admitted']}, 이 제한된 한 번의 반환에서 기존 부호 유지={r['conditional_compact_sign_survives_one_return']}. 전체 기간249반환 뒤 판독은 별도251연결로 이어진다. 게시={now}.

분류: Counterexample candidate. 긴246연속 원천의 원 시간 대조 최대는1.232083616174%이며, 같은 소스의247전하 장 시간 차이는U={r['time']['U']:.12e},U_t={r['time']['U_t']:.12e},U_x={r['time']['U_x']:.12e}다. 반경 구적={r['controls']['quadrature']:.12e}, 독립 GR 적분={r['controls']['independent_GR']:.12e}. 원2%시간·0.2%구적·1e-9독립 적분 기준을 그대로 적용했다. 초기188국소 실패를 지우거나 전역 통과로 대체하지 않는다.

분류: Counterexample candidate. fine compact 전하의 high={float(fine['high']):.12e},low={float(fine['low']):.12e},low/high={fine['low_over_high']:.12e}. 두 성분은 동일한 실제 반환 해와 GR 시각·반경 연산자의 성분이며, 서로 다른 진단 해를 더한 값이 아니다. 작은 성분을 분리하고 최종 이진 입력을100자리로 합산하는 방식은100자리 물리 정확도를 뜻하지 않는다. 부호 판정은 저장된 유한 연산자와 한 번의 보상 반환에 조건부다.

분류: Counterexample candidate.244전체 원천이 완료된 뒤248의 완전 prefix동일성 검사에서 기하 계수 차이가 드러났다. 이전 종료 시각과 내부 canonical시각 사이1ULP차이 때문에, coarse의 마지막 기하 구간 기울기가 달라지고 fine에는 매우 짧은 종단 구간이 따로 있었다. 전체 원천이나 과거 파일을 수정하지 않았다. canonical14/16이전 마지막 출력까지 모든 필요한 원천 계수·끝점·배경이 비트 단위로 같음을 검증하고,494개 GR출력만 재사용하며81개를 다시 계산하도록 변경했다. 원래535개 모두를 재사용하려던 실패는 보존한다. 두 겹침 출력과 원1e-9/1e-12기준으로 캐시 일관성을 다시 판정한다.

분류: Counterexample candidate. 이 수정 방식으로 세 전체기간 GR장을 실제 계산해 수락했다. 세 경로의 겹침 출력에서 free/potential/U/U_t/U_x차이는 모두0이었으며, 전체 시간 대조 최대3.507970401735e-5, 구적8.107722479311e-16, 독립 적분1.278517805621e-14로 기존 기준을 통과했다. 원천 수정 없이 성립하는 범위로 재사용을 줄인 실제 계산 결과다. 이 입력의 실제 전체기간 반환은249에서 이어진다.

분류: Counterexample candidate. 새 전체기간 판독 경로는249실제 반환,244원천의 질량 제약,248high장을 사용한다. 기존245–247소스/전하 연산자를 재사용하며 새 물질·광자 적분은 없다. 짧은227해를 이 연결로 끝까지 읽어 끝점·기하·연속 원천·세 GR장의 모든 저장 배열 및 두 성분 전하가 기존 결과와 정확히 같음을 확인했다. 첫 병렬 확인에서는 legacy초기화가 공용 generated코드를 동시에 써 읽기 실패가 발생했다. 원 실패·producer를 보존했고 초기화에만 파일 잠금을 걸어 전체 재현을 통과했다. 본 GR장은 계속 병렬 계산한다.

분류: Conjectural.248수정 경로는 실측5query비용에 실제 출력 수와 원천 구간 수를 적용한 nominal약16.4분, 계획상 상단약37.8분이다. 각 field2시간·16GiB를 유지한다.251후속 판독에는 endpoint90분·geometry60분·dense source120분·각 GR장3시간·16GiB를 둔다. 실제 dispatch직전247실측과 입력 크기로 예측을 갱신한다. 원 정확도 기준·전체 기간·두 물리 경로를 유지하며 자동 격자 확대는 없다.

06:05KST무렵 WSL부팅 식별자가 변경돼 실행 중인 reader와 대기 controller가 종료됐다. 사용자가 직접 재시작했음을 확인했으며 과학적 수락 실패로 세지 않는다. 기존 부팅 식별자·상태·로그·중단 시 생성물과 새 시작 정보를 보존했다. 이미 완료된 물리 해, source/geometry와247의580초 구적장을 재사용하고 receipt가 없는 미완료 action만 다시 실행했다. 계산 producer와 기준을 바꾸지 않았으며, 공용 legacy초기화만 전체 reader에 걸친 파일 잠금으로 직렬화했다. 판독과 실제 반환의 부팅별 시작 증거는 게시 스냅샷과 interruption기록에 둔다.

분류: Counterexample candidate.249재시작 준비에는 두 경로의 시각을 비트 단위로 같다고 요구하는 과도한 전제가 있었다. 저장 coarse/fine종료 시각은4.336808689942018e-19초 차이이며, 각각은 자신의 원 경로와 비트 단위로 같다. 두 시각과 상태를 바꾸지 않고 기존 StageDriver의1e-18초 정렬 규칙으로 연결했다.1e-17초 어긋난 부정 대조는 거절했다. 원 물리 단계·paired시간·구적·보존 기준은 그대로다. 첫 재개 때 남은 불완전 seed디렉터리의 존재 검사에도 걸렸으며, 해당 실패와 모든 링크를 별도 위치에 보존한 뒤 새 작업 디렉터리로 재개했다. 물리 단계 실행 전의 준비 오류였다.

분류: Conjectural. 남은 전체 판정은249전체 기간 반환 및251같은 해 전하 판독, 시간 표현의 균일 오차·EOS/미분·공간/경계·자기GR 고정점과 완전 비선형·같은 재고의 정적 비교·관측·무한대 전하 정규화다. 본 단계는 loophole progress이며, 국소 수치 재현이나 compact부호를 최종 물리 전하의 완료로 세지 않는다.
''',encoding='utf-8')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        p=root/f'docs/{name}.md';prefixes[p.relative_to(root).as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
        with p.open('ab') as f:f.write((f'\n\n## 단계251 — 같은 반환 해의 긴 전하 판독과 전체기간 연결\n\n분류: Counterexample candidate. 최종 물리 전하 미판정. 공통15/16실제 해의 compact전하 대조 수락={r["charge_comparison_admitted"]}, 조건부 부호 유지={r["conditional_compact_sign_survives_one_return"]}. 전체 원천의 종단 시간 표현 때문에494개 정확히 같은 과거 GR출력만 재사용하고81개를 재계산한다. 원 source-prefix실패를 보존했다. 전체249반환을 동일 해의 전하 판독에 연결했으며, 짧은 해의 모든 원천·GR배열과 전하를 정확히 재현했다. 전체 EOS·미분·공간·경계·비선형/자기GR·정적/관측·무한대 범위와 기존 실패는 유지한다. [근거](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',actual_same_joint_solution_GR_computed=True,
        common_15of16_charge_admitted=r['charge_comparison_admitted'],conditional_compact_sign_survives_one_return=r['conditional_compact_sign_survives_one_return'],
        same_solution_compact_components=r['components'],full_period_readout_regression_passed=True,
        exact_prefix_extension_repaired=True,immutable_runtime_bindings=bindings,live_snapshots=snapshots,snapshot_KST=now,
        full_goal_complete=False,physical_final_charge_solved=False,final_charge_conclusion='unadjudicated')
    write(out/'result.json',final);write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes,reused=bindings))
    files=modules+[note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes,sha256={p.relative_to(root).as_posix():sha(p) for p in files});write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_full_return_charge_route']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
