"""Publish the applied dense GR bridge, including its original arithmetic failure."""
from pathlib import Path
from types import FunctionType
import importlib.util,shutil,sys
helper=Path('outputs/direct-eos-gr33/native-radau-pulse-gr/previous-publication-producer.py')
spec=importlib.util.spec_from_file_location('b',helper);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha=b.read,b.write,b.sha
out=root/'outputs/direct-eos-gr33/native-complete-radau-bridge';manifest=out.parent/'native-complete-radau-bridge-manifest.json'
module=root/'verification/read_complete_radau_history.py';note=root/'notes/REQUEST224_DENSE_GR_BRIDGE_KO.md'


def copy(src,dst):
    h=sha(src);dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);assert h==sha(dst)==sha(src)


def package():
    assert not out.exists();out.mkdir();work=runtime/'native-complete-radau224-check-work';control=runtime/'native-complete-radau224-control'
    r=read(work/'result.json');s=read(work/'sources.json');audit=read(work/'driver-polynomial-audit.json')
    assert r['numerical_controls_passed'] and s['representation_controls_passed'] and audit['passed']
    assert read(work/'arithmetic-regression.json')['passed'] and read(work/'symbolic.json')['passed']
    for action in ['endpoints','source','fields']:assert read(work/f'{action}-receipt.json')['error'] is None
    tested=work/'tested-consumer.py';assert sha(tested)==read(work/'fields-receipt.json')['consumer_sha256']
    assert sha(module)==sha(runtime/'verification'/module.name)==read(control/'start.json')['producer_sha256']
    old=read(out.parent/'native-original-photon-capture-manifest.json');preserved={p:h for p,h in old['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    reused={}
    for src in work.rglob('*'):
        if not src.is_file() or '__pycache__' in src.parts:continue
        rel=src.relative_to(work)
        if rel.parts[0] in ['sweep-0','sweep-1']:
            reused[rel.as_posix()]=dict(runtime_path=str(src),sha256=sha(src));continue
        copy(src,out/'completed'/rel)
    for src,name in [(root/'.phase224-prefix-check.py','prefix-producer.py'),(root/'.phase224-energy-probe.py','energy-probe-producer.py'),
                     (root/'.phase224-followthrough.py','followthrough-producer.py'),(control/'start.json','followthrough-entry.json'),(Path(__file__),'publication-producer.py')]:copy(src,out/name)
    maximum=max(max(row['errors'].values()) for row in audit['rows'])
    note.write_text(f'''# 수락된 같은 이력의 배경 구간을 연결한 GR 판독

분류: Counterexample candidate. **최종 전하의 결론은 미판정이다.** 실제 수락된T/8 이력에서 배경 구간 전환·Radau 단계·원 입사 파형·에너지/출구를 같은 GR 판독에 연결했다. 이 연결은 원 수치 대조를 통과했다. 전체 기간, 자기GR·비선형 반환과 무한대 최종 전하를 완료한 것은 아니다. 연구 분류는 loophole progress의 같은 해 판독 연결이다.

분류: Counterexample candidate. 전체 공통 이력을 읽기 전에223의 기존 수락 체크포인트에서 거친15/미세29개 물리 단계를 동결했다. 두 경로의 원T/8 스냅샷은 끝점·물질 수지·각도/반경 출구를 통과했다. 실제 물리 적분이나 새 시간 경로는0개다.224는 기존206 끝점 판독과212 밀집 시간 판독을 재사용하면서, 첫 배경 구간을 전체에 외삽하지 않고 각 실제 배경 구간의 계수를 연결한다. 원8차 펄스와 저장 Born 항을 그대로 사용하며 실제 파형을 부분구간 중점마다 독립 대조한다.

분류: Counterexample candidate. 첫 실제 적용은 배경 에너지 상쇄 때문에 원천 다항식 검사가1.415335158e-12>1e-12로 실패했다. 원 끝점·물질 수지·단계 방정식·압력 대응은 통과했지만 이 실패를 대체하지 않았다. 실패한 단계14의5/6 지점에서 동일 이진 입력을50·80자리로 계산한 값은 long-double로 정확히 같았다. 기존 배경 에너지 대비항의 오차는 전체 판독량 기준8.22193e-13이었다. 단순히 끝점별로 계산 순서를 바꾼 결과도6.43892e-13의 차이가 남았으므로, 해당 대비항을80자리로 계산한 뒤 원 선형 시간 보간에 넣었다. 물리식·EOS 입력·해·수락 기준은 바꾸지 않았다. 원 실패·계획·소스·receipt와 SHA 바인딩 수정 중 발생한 실행 전 실패도 보존했다.

분류: Proven. 배경의 선형 보간과 에너지 대비항 계산을 교환하는 정확 실수 항등식을 기호 검사했다. 기존 한 구간 소비자의 계수는 정확히 재현했고, 같은 선형 배경을 두 구간으로 분할한 대조도 통과했다. 이 항등식과 고정 입력의 산술 대조는 전체EOS 정확도나 균일 미분 오차의 증명이 아니다.

분류: Counterexample candidate. 수정된 동일 이력의 두 경로에서 원천 표현 최대 잔차는9.77514e-16, 원 끝점과의 차이는2.50491e-16 이하였다. 독립 실제 파형 대조 최대 상대 차이는{maximum:.9e}이다. 원 공간L1 시간 대조의 최대 차이는0.068799퍼센트로2퍼센트 기준을 통과했다. 이전T/64의 실패는 그대로이며, 더 긴 구간의 전역 척도로 통과한 것을 이전 국소 오류의 소멸이나 균일 시간 오차 보증으로 해석하지 않는다.

분류: Counterexample candidate. 이 수정된 원천을 실제 지연 GR에 적용한 시간 경로 차이는U {r['time']['U']:.9e}, U_t {r['time']['U_t']:.9e}, U_x {r['time']['U_x']:.9e}이다.4/8차 공간 적분 대조는{r['controls']['quadrature']:.9e}, 독립 직접 GR 적분 대조는{r['controls']['independent_GR']:.9e}이다. 원천 시간 통과={r['source_time_passed']}, 장 시간 통과={r['field_time_passed']}, 이 제한 구간의 기존GR 반환 수락 조건={r['GR_return_admitted']}. 실제 물질에 GR을 다시 반환한 새 진화는 아직0개이며, 출력의 compact 전하나 부호를 무한대 최종 전하로 계승하지 않는다.

같은 코드의 후속 소비자는223이 원 기준으로 전체 공통15/16 이력을 수락한 경우에만 시작한다. 원천/장 시간 대조가 실패하면 그 판정을 보존하며 물리 반환을 자동 수락하지 않는다. 해상도·물리 기간·경로를 추가하지 않는다. 실행은CPU1스레드·12GiB, 끝점40분·원천90분·GR장120분의 여유 예산을 둔다. 이번T/8의실측을 전체 구간에 외삽하면 출력 시각 수와 파형 경계 절단 수에 따라 시간이 늘어나며, 후반 비용은 보장이 아니다. 기존218 실제 결합 진화와223 광자 복원은 고정 소스·계획을 유지한다.

분류: Conjectural. 전 기간의 수락된 같은 해에서 지배 오차를 해결한 뒤 자기GR·비선형 반환, EOS/미분·공간/경계 오차, 동일 재고 정적 비교와 관측량·무한대 전하를 연결해야 한다. 작은 구간의 수치 수락은 최종 연구 성과의 대체가 아니다.
''',encoding='utf-8')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        p=root/f'docs/{name}.md';prefixes[p.relative_to(root).as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
        with p.open('ab') as f:f.write(('\n\n## 단계224 — 배경 구간을 연결한 동일 이력 GR 판독\n\n분류: Counterexample candidate. 최종 전하는 미판정이다. 수락된T/8 물질·광자·출구 이력을 실제 배경 구간과 원 입사 파형으로 GR에 연결했다. 배경 에너지 상쇄로 실패한 원 다항식 기준을80자리 동일식 계산으로 해결하고 원 수치 대조를 통과했다. 원 실패와 이전T/64 시간 실패는 유지한다. 전체 공통 이력이 수락되면 같은 판독기를 적용하며, 전 기간·자기GR·비선형·최종 전하의 수락은 아직 아니다. [근거](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',actual_same_joint_solution_GR_computed=True,
        same_solution_prefix_result=r,source_representation_passed=True,old_arithmetic_failure_preserved=True,
        full_common_history_followthrough_armed=True,full_declared_period_completed=False,self_GR_return_completed=False,
        final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(out/'result.json',final);write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes,reused=reused))
    files=[module,note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes,sha256={p.relative_to(root).as_posix():sha(p) for p in files});write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest);m['native_complete_radau_bridge']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
