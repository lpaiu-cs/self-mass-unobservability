"""Publish the applied GR time repair and its unchanged physical rejection."""
from pathlib import Path
from types import FunctionType
import importlib.util,shutil,sys
p=Path(__file__).with_name('.phase186-publish.py')
if not p.exists():p=Path(__file__).with_name('previous-publication-producer.py')
spec=importlib.util.spec_from_file_location('p',p);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha=b.read,b.write,b.sha
out=root/'outputs/direct-eos-gr33/native-radau-pulse-gr';manifest=out.parent/'native-radau-pulse-gr-manifest.json'
folders=[runtime/n for n in ['native-radau-gr211-work','native-driver-radau212-work']]
modules=[root/'verification'/n for n in ['apply_radau_history_gr.py','apply_driver_aware_radau_gr.py']]
note=root/'notes/REQUEST211_212_RADAU_PULSE_GR_KO.md'


def package():
    assert not manifest.exists();failed,w=folders;r=read(w/'result.json')
    assert read(failed/'source-receipt.json')['error'] and read(failed/'source-64-check.json')['polynomial_max']>.04
    assert r['numerical_controls_passed'] and not r['source_time_passed'] and not r['field_time_passed'] and not r['GR_return_admitted']
    assert read(w/'source_retry-receipt.json')['error'] is None and read(w/'fields_finish-receipt.json')['error'] is None
    assert read(w/'driver-polynomial-fixed-audit.json')['passed'] and not read(w/'driver-polynomial-audit.json')['passed']
    assert all(sha(m)==sha(runtime/'verification'/m.name) for m in modules)
    assert sha(w/'gr-dispatch-producer.py')==read(w/'fields-receipt.json')['source_sha256']
    for name,h in read(w/'fields-reuse.json')['bindings'].items():assert sha(runtime/name)==h,name
    old=read(out.parent/'native-momentum-entry-manifest.json');preserved={name:h for name,h in old['sha256'].items() if not name.startswith('docs/')}
    for name,h in preserved.items():assert sha(root/name)==h,name
    known={h:dict(path=name,sha256=h) for name,h in read(master)['sha256'].items() if not name.startswith('docs/')}
    known.update({v['sha256']:v for v in read(out.parent/'native-momentum-entry/publication.json')['reused'].values()})
    aliases={sha(f):f for folder in folders for f in folder.glob('*.py')};bindings=[]
    for folder in folders:
        for name,h in read(folder/'plan.json')['bindings'].items():
            actual=runtime/name.removeprefix('/home/lpaiu/work/native-retained-tail-runtime/')
            if sha(actual)!=h:actual=aliases[h]
            assert sha(actual)==h,name;bindings.append(dict(plan=folder.name+'/plan.json',source=name,resolved=str(actual),sha256=h))
    out.mkdir();reused={}
    for folder in folders:
        for src in sorted(folder.rglob('*')):
            if not src.is_file():continue
            digest=sha(src);key=src.relative_to(runtime).as_posix()
            if digest in known and (root/known[digest]['path']).exists():
                assert sha(root/known[digest]['path'])==digest;reused[key]=known[digest];continue
            dst=out/key;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);assert sha(src)==sha(dst)==digest
            known[digest]=dict(path=dst.relative_to(root).as_posix(),sha256=digest)
    for src,name in [(Path(__file__),'publication-producer.py'),(p,'previous-publication-producer.py'),(b.helper,'publication-helper.py'),(root/'.phase212-driver-audit.py','driver-audit-producer.py'),(root/'.phase212-amplitude-fix.py','amplitude-fix-producer.py')]:shutil.copyfile(src,out/name)
    note.write_text('''# 실제 Radau·입사 펄스를 GR 시간 원천에 적용

분류: Counterexample candidate. **최종 전하는 미판정이다.** 같은 저장 결합 해의 단계 내부 다항식과 실제 입사 펄스를 GR 연산자에 적용했다. 수치 연산 대조는 통과했지만 두 시간 경로의 U차이는40.8632%로 원2%기준을 실패했다. 보간 수정만으로 초기 결합 해의 시간 오차가 해소된다는 가설을 채택하지 않는다. 실제 자기GR 반환은 수락하지 않았다. 연구 분류는 loophole progress다.

분류: Counterexample candidate. 기존188첫T/64=53.662986마이크로초의 coarse2/fine4단계만 재사용했다. 물리적 새 적분·격자·기간·경로는 없다. 물질의 실제 초기 상태와 두 Radau 단계, 같은 광자의 재구성 모멘트 및 실제 반경 유속을 사용한다. 각 구간의 시작은 post-floor, 닫힘은 pre-floor로 두어 floor 점프를 임의로 평탄화하지 않는다. 실제 Radau 변화율과 단계 복원 차이는 coarse3.79204e-14/fine2.14988e-13로 원1e-12기준을 통과했다.

분류: Counterexample candidate. 211은 전체 원천을3차식으로 놓았다가 중간 시각 읽기에서4.74956%차이를 보여 실패했다. 물질 상태는2차 Radau식이지만 외부 입력은 움직이는8차 compact pulse이므로 전체 원천을 같은 낮은 차수로 놓을 수 없다. 실패 생산 코드·원천·receipt를 보존했다.

분류: Counterexample candidate. 212는 실제 물질·광자 상태의 zero-geometry3차 읽기와, 시간에 선형인 배경 기하 대응×정확한8차 입사 펄스·저장 Born성분을 분리해 동일 원천을 구성했다. 기존 단계/Born시각 및 각 셀의 실제 pulse 도달·종료 시각에서만 다항식을 나누었다.517/519개 구간은 GR 원천의 대수적 지지집합 분할이며 새 물리 적분 경로가 아니다. 독립 중간 시각 전체 읽기의 차이는3.26081e-13이하, 기존 원천 끝점 차이는6.20811e-17이하, 압력 탐침1.32931e-10이하, 압력 대응1.962e-16이하로 원 기준을 통과했다.

분류: Counterexample candidate. 전달 직전 live Driver대조에서 응답 좌표의 정규화AMP를 실제 구동ETA로 잘못 저장한 메타데이터 오류를 발견했다. 실제 구동 진폭으로 정정하고 나머지 모든 물리 배열이 비트 단위로 같음을 확인했다. 잘못된 원천 파일·실패 대조·수정 근거를 보존했다. 정정 뒤 모든517/519구간의 중간점에서 실제 Driver와 최종 소비 다항식의 차이는3.796e-16이하로 통과했다. 물리 구동 진폭을 바꾼 실험이 아니다.

분류: Proven. Radau보간의 단계 항등식과, 지지 경계·특성선 절단을 포함한0–9차 단항 원천의 U/U_t/U_x대조를 통과했다. 이것은 선언한 유한 표현의 대수 검산이며 EOS·균일 미분·공간·시간 연속체 오차 증명이 아니다.

분류: Counterexample candidate. 수정 원천의 실제 GR3경로 계산은21.234/14.180/13.074초에 끝났다. 독립 반경 적분의 numpy.interp자료형 오류는 완료한 장을 재계산하지 않고 고쳤다. 원4/8차수 대조3.79457e-14와 독립 Jordan반경 적분1.98335e-12는 원0.002/1e-9기준을 통과했다. U시간 차이는40.8632%,U_t15.1682%,U_x10.2553%다. 기존 직선 시간 표현과 fine U의 차이는22.8402%이며 표현 변경이 실제 GR값에 반영됐다. 원천 끝점의 기존 B6.98376%·광자E7.49853%등 시간 실패도 그대로 남는다.

분류: Counterexample candidate. 같은 짧은 구간의 compact끝점은 fine=-2.987085482e-63,coarse=-8.975356624e-64다. 같은 부호만으로 안정된 무한대 전하나 관측량을 주장하지 않는다. 저장 시각에서 potential반복 성분은 전체U대비1.90231e-15였고 그 항의 시간 표현은 여전히 기존 직선이다. 원 physical self-GR·최종 전하·전체 EOS/균일미분/공간/경계/완전비선형/정적/관측 요건은 미완료다.

분류: Conjectural. 지정된 초기 반환 실험을 수락하려면 초기 구동의 도달·시간별 native분기를 포함한 원 해의 시간 정확도를 해결해야 한다. 이 첫T/64의 국소 상대오차를 전체 기간 최종 전하의 지배 오차로 곧바로 해석하지 않는다. 전체 결합 해의 최종 판독에서는 같은 에너지·경계 이력으로 오차를 다시 평가해야 한다. 본210후반 속행은 별도 고정 프로세스이며 이 GR읽기 결과로 그 종료나 다음 단계 수락을 대신하지 않는다.

계산 예산: 원천·장 각각20분,6GiB/CPU1스레드. 원천은45.812초,완료 장 재사용 후 독립 판독은약7초였다. 실제 생산 코드·실패·수정·입력SHA를 함께 보존했다.
''',encoding='utf-8')
    tails={
        'model-definition':'실제 Radau상태와 정확한 입사 pulse/Born시간식을 같은 GR원천에 연결했다. 기하 응답은 선형 배경 대응과 원 펄스의 곱으로 유지하며 floor점프·실제 경계 유속을 보존한다.',
        'observable-targets':'최종 전하는 미판정이다. 시간 표현 수정 후 실제U대조40.8632%로 원2%를 실패했다. 짧은compact끝점 부호를 무한대 전하나 관측량으로 해석하지 않는다.',
        'adiabatic-limit':'0–9차 특성선 대수 대조와 실제 Driver다항식 일치를 확인했다. 물리 EOS·균일 미분·공간·경계·연속 시간의 오차 상계가 아니다.',
        'nonadiabatic-regime':'정확한 입사 지지 시각과 실제 Radau이력을 GR에 적용했다.4/8차수3.79457e-14·독립적분1.98335e-12는 통과했지만U시간40.8632%를 실패했다. 표현 수정만으로 원 해의 시간 오차가 해소되지 않았다.',
        'failure-ledger-dynamic-chi':'전체 원천3차 가정, 진폭 메타데이터 및 독립 읽기 자료형 실패를 보존하고 수정했다. 정정된 실제GR판독에서도 원천과U시간 기준은 실패다. 원 실패나2%문턱을 완화하지 않는다.',
        'dynamic-charge-completion':'최종 전하는 미판정이다. 실제Radau·입사펄스 시간 표현을 같은GR연산자에 적용했으나 원 해의 시간오차가 남았다. 실제자기GR반환·무한대전하·전체EOS/미분/공간/경계/비선형/정적/관측 요건을 계속 유지한다.'}
    prefixes={}
    for name,body in tails.items():
        f=root/f'docs/{name}.md';prefixes[f.relative_to(root).as_posix()]=dict(bytes=f.stat().st_size,sha256=sha(f))
        with f.open('ab') as h:h.write(('\n\n## 단계211–212 — 실제 Radau·입사 펄스의 GR 적용\n\n분류: Counterexample candidate. '+body+' [근거](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',result=r,actual_same_joint_solution_GR_computed=True,GR_readout_accepted=False,GR_return_evolution_executed=False,
        original_failures_preserved=True,physical_final_charge_solved=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(out/'final-result.json',final);write(out/'publication.json',dict(reused=reused,previous_nondoc=preserved,document_prefixes=prefixes,plan_bindings=bindings))
    files=modules+[note]+[root/k for k in prefixes]+[f for f in out.rglob('*') if f.is_file()]
    final.update(sha256={f.relative_to(root).as_posix():sha(f) for f in files},document_prefixes=prefixes);write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest);m['native_radau_pulse_GR']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
