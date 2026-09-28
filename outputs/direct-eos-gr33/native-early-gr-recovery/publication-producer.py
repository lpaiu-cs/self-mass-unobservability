"""Publish207/208; neither original source failure nor final-charge scope changes."""
from pathlib import Path
from types import FunctionType
import importlib.util,shutil,sys
p=Path(__file__).with_name('.phase186-publish.py')
if not p.exists():p=Path(__file__).with_name('previous-publication-producer.py')
spec=importlib.util.spec_from_file_location('p',p);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha=b.read,b.write,b.sha
out=root/'outputs/direct-eos-gr33/native-early-gr-recovery';manifest=out.parent/'native-early-gr-recovery-manifest.json'
folders=[runtime/n for n in ['native-early-gr207-work','native-equation-recovery208-work']]
modules=[root/'verification'/n for n in ['measure_early_gr_time_error.py','complete_joint_photon_recovery.py']]
note=root/'notes/REQUEST207_208_GR_TIME_ERROR_AND_PHOTON_HISTORY_KO.md'


def package():
    assert not manifest.exists();g,r=[read(w/'result.json') for w in folders]
    assert g['diagnostic_controls_passed'] and not g['source_admission_passed']
    assert r['original_endpoint_and_ledger_passed'] and not r['strict_archival_rate_identity_passed']
    assert read(folders[0]/'measure-receipt.json')['error'] is None and read(folders[1]/'fine-receipt.json')['error'] is None and read(folders[1]/'audit-receipt.json')['error'] is None
    assert all(sha(m)==sha(runtime/'verification'/m.name) for m in modules)
    old=read(out.parent/'native-resolved-return-manifest.json');preserved={p:h for p,h in old['sha256'].items() if not p.startswith('docs/')}
    for name,h in preserved.items():assert sha(root/name)==h,name
    known={h:dict(path=name,sha256=h) for name,h in read(master)['sha256'].items() if not name.startswith('docs/')}
    known.update({v['sha256']:v for v in read(out.parent/'native-resolved-return/publication.json')['reused'].values()})
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
    for src,name in [(Path(__file__),'publication-producer.py'),(p,'previous-publication-producer.py'),(b.helper,'publication-helper.py')]:shutil.copyfile(src,out/name)
    note.write_text('''# 초기 GR 시간 오차와 같은 해의 광자 이력 복원

분류: Counterexample candidate. **최종 전하는 미판정이다.** 207은 이미 실패한 초기 원천의 오차가 GR 읽기에서 작아지는지 측정했고, 오히려 큰 차이가 남았다. 208은 누락된 fine광자 단계 이력을 원 방정식·끝점·수지 기준으로 완결했다. 어느 결과도 자기 GR·최종 전하의 해결이나206원천 시간 실패의 통과가 아니다.

분류: Counterexample candidate. 207은 같은206원천의coarse/fine 및 fine을coarse저장 시각으로 줄인 자료를 동일 특성선 GR 연산자로 읽었다. 세 번째는 새 물리 경로가 아닌 저장 자료의 투영이다. 공통 시각의 coarse-fine 차이를 원천 상태 차이와 저장 시각 사이 표현 차이로 정확히 분해했다. U의 최대 차이25.8719%,상태 항11.9050%,표현 항35.2599%다. 각 항의 최대 위치가 다르고 서로 상쇄하므로 이 크기를 단순히 더하지 않는다. 분해 잔차는1.13e-16이하,4/8차수 차이2.695e-14,독립 직접 GR적분 차이1.740e-12다. 원2%시간 기준이 해소된 것이 아니다.

분류: Counterexample candidate. 해당 첫T/64compact끝점은 fine=-4.068961503e-63,coarse=-5.675710782e-63,투영fine=-6.258745133e-63다. 이는 각 저장 원천의 compact 진단 읽기이며 무한대 최종 전하, 실제 자기 결합 해의 전하, 검출 가능한 관측량으로 해석하지 않는다. 새로운 물리 진화는 없었다.

분류: Counterexample candidate. 208은205의 추가적인 저장H교환률 일치 실패를 그대로 유지했다. 원 전체 결합식·물리 모멘트 기준으로 저장 거부 제안을 독립 검사했던 근거에 따라, 그 판정과 구분되는 복원 실험을 선언했다. 기존5단계의 엄격한 통과 이력과 여섯째 단계 제안을 비트 그대로 재사용하고 실제 새 조건부 광자 풀이는 마지막2단계에만 수행했다. 기체 이력·원 물리식·시간 격자·끝점은 바꾸지 않았다. 모멘트를 저장값에 맞춰 적합하지 않았다.

분류: Counterexample candidate. fine8단계 복원은88.036초에 완료했다. 남은 단계5/6/7의 전체 결합식 상대 잔차는1.11368e-15/1.76012e-16/1.31497e-16으로 원1e-12기준, 물리 잔차 최대3.83412e-15로 원1e-13기준을 통과했다. 실제 native변화율은 각각 원 저장값과 정확히 같았다. 앞5단계는 기존 엄격한 원천·광자 검사 근거를 재사용했으며 새 전체 상태 재검사를 했다고 주장하지 않는다.

분류: Counterexample candidate. 원T/32=107.325972마이크로초 저장 끝점의 광자 상태 차이3.17034e-17<1e-12,전체 재구성 충돌을 사용한 같은 물질 수지 최대9.87996e-15<1e-8,각도 출구1.85150e-16<1e-12,반경 출구7.37552e-13<1e-12로 통과했다. 이로써 coarse4/fine8의 모든 실제 광자 모멘트·물질 단계·반경/각도 출구가 더 촘촘한 GR 입력 재구성에 이용 가능하다. 정확한 저장 교환률 일치와 전체 기간 복원은 여전히 수락하지 않았다.

분류: Conjectural. 다음 수정은 드문 끝점만 잇는 GR 원천 표현에 실제 Radau 단계 및 연속 다항식 정보를 반영하는 것이다. 기존 초기 원 해의 시간 차이도 별도로 남아 있어, 표현 개선만으로 최종 전하가 안정된다고 예측하지 않는다. 수정은 같은 해와 경계 이력에 적용하고 원 실패·원 수락 기준을 보존해야 한다. 본 결합202는 고정된 별도 프로세스에서 속행 중이며 이 묶음은 그 종료 판정이 아니다.

분류: Proven. 기존 특성선 다항식 자체 검산을 통과했다. 선형 연산자 출력의 차이 분해는 대수 항등식이다. 이 사실들은 물리 EOS·균일 미분·연속 시간·공간·경계·비선형 폐쇄의 증명이 아니다.

실행 예산:207기존 자료 읽기600초,208fine1800초·audit180초,각6GiB/1스레드. 전체 물질 이력을 반복 적분하거나 새 물리 격자·기간을 늘리지 않았다. 실제 생산 코드·계획·수치·의존성 SHA와 이전 실패를 함께 보존했다.
''',encoding='utf-8')
    tails={
        'model-definition':'원 실제 기체 이력을 고정한 광자 조건부 복원을T/32끝점까지 마쳤다. 앞5단계의 엄격한 검사와 남은3단계의 원 전체식 검사를 구분하고 같은 충돌·끝점·경계 수지를 연결했다.',
        'observable-targets':'최종 전하는 미판정이다. 초기U시간차25.8719%가 남고 같은 짧은구간compact읽기를 무한대 전하로 승격하지 않았다. 원 방정식에 맞는 누락 광자 이력을 확보했다.',
        'adiabatic-limit':'시간 표현과 상태 차이의 선형 읽기 분해 및 특성선 검산은 전체 균일 오차 증명이 아니다. 엄격한 저장H교환률 실패는 원 전체 단계식 통과와 구분해 유지한다.',
        'nonadiabatic-regime':'초기GR차이의 시간 표현 항이35.2599%,상태 항이11.9050%의 최대 크기를 보였다. 두 항은 상쇄하며 원시간 기준은 실패다. 같은fine8단계광자복원은 끝점3.17e-17와 원수지·출구 기준을 통과했다.',
        'failure-ledger-dynamic-chi':'206초기원천시간 실패는 GR연산자에서도 해소되지 않았다. 205추가교환률일치 실패를 통과로 바꾸지 않고, 구분된208원 전체식/끝점/수지 복원을 완료했다. 앞5단계를 새 전체식 검사로 세지 않는다.',
        'dynamic-charge-completion':'최종 전하는 미판정이다. 실제 단계 이력을 확보해 GR시간 표현을 고칠 입력이 준비됐지만 수정된 자기GR진화·최종전하는 아직 없다. 본 결합202와 물리EOS/미분/공간/경계/비선형/정적/관측 요건은 계속 유지한다.'}
    prefixes={}
    for name,body in tails.items():
        f=root/f'docs/{name}.md';prefixes[f.relative_to(root).as_posix()]=dict(bytes=f.stat().st_size,sha256=sha(f))
        with f.open('ab') as h:h.write(('\n\n## 단계207–208 — GR 시간 오차 전달과 광자 이력 완결\n\n분류: Counterexample candidate. '+body+' [근거](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',GR_error_transfer=g,photon_recovery=r,original_failures_preserved=True,actual_same_joint_solution_GR_computed=False,diagnostic_compact_GR_computed=True,
        GR_return_evolution_executed=False,physical_final_charge_solved=False,self_GR_return_closed=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(out/'final-result.json',final);write(out/'publication.json',dict(reused=reused,previous_nondoc=preserved,document_prefixes=prefixes,plan_bindings=bindings))
    files=modules+[note]+[root/k for k in prefixes]+[f for f in out.rglob('*') if f.is_file()]
    final.update(sha256={f.relative_to(root).as_posix():sha(f) for f in files},document_prefixes=prefixes);write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest);m['native_early_gr_recovery']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
