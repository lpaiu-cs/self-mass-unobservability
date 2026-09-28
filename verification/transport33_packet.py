"""Freeze the completed transport/source evidence without closing Request33."""
import json, shutil, sys
from pathlib import Path
import numpy as np
import sympy as sp
import direct_eos_gr as g

OP=g.OUT/'opacity';OUT=g.OUT/'transport-checkpoint'
CODE=['native_opacity','opacity_electrons','opacity_derivatives','opacity_internal','opacity_cubic',
    'opacity_cubic_certificate','opacity_cubic_propagation','opacity_tables','opacity_table_diagnostics',
    'opacity_mixed_term','source_retry','gr_opacity_check','transport33_packet']


def read(path): return json.loads(path.read_text())


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def audit():
    OUT.mkdir(exist_ok=True)
    aliases={
        'verification/native_opacity.py':[OP/n for n in ['before-type1-native_opacity.py',
            'before-internal-native_opacity.py','before-cubic-native_opacity.py','before-tables-native_opacity.py']]+
            [OP/'cubic/before-context-native_opacity.py',OP/'tables/before-enrichment-native_opacity.py'],
        'verification/opacity_tables.py':[OP/'tables'/n for n in ['before-scalar-compat-opacity_tables.py',
            'before-enrichment-opacity_tables.py','before-mixed-term-opacity_tables.py']],
        'verification/opacity_cubic_certificate.py':[OP/'cubic/certificate/before-piecewise-certificate.py']}
    bindings=[]
    for plan in OP.rglob('*plan.json'):
        item=read(plan)
        for key in ['inputs_sha256','input_sha256']:
            if not isinstance(item.get(key),dict): continue
            for rel,digest in item[key].items():
                path=g.ROOT/rel
                candidates=[path,*aliases.get(rel,[])]
                matches=[p for p in candidates if p.exists() and g.c.sha(p)==digest]
                assert matches,(str(plan),rel,digest)
                bindings.append(dict(plan=plan.relative_to(g.ROOT).as_posix(),input=rel,
                    verified_source=matches[0].relative_to(g.ROOT).as_posix(),sha256=digest))
    source=read(g.OUT/'reference-source.json');assert source['completed']
    raw=dict(np.load(g.OUT/'reference-cached-native.npz'));aux=np.load(g.OUT/'reference-auxiliary.npy')
    assert np.array_equal(raw['aux_used'],aux)
    assert g.c.sha(g.OUT/'reference-cached-native.npz')==source['native_source_sha256']
    assert g.c.sha(g.OUT/'reference-cached-corrected.npz')==source['corrected_source_sha256']
    corrected=g.s.v.corrected_values(raw);saved=dict(np.load(g.OUT/'reference-cached-corrected.npz'))
    assert all(np.array_equal(value,saved[key]) for key,value in corrected.items())
    old=dict(np.load(g.OLD/'initial-native.npz'));identity=dict(np.load(g.OUT/'cached-identity-v2-native.npz'))
    assert all(np.array_equal(value,identity[key]) for key,value in old.items())
    for label in ['cached-identity-v2','reference-cached']:
        trace=read(g.OUT/(label+'-trace.json'));assert trace['calls']==5735 and trace['elapsed_s']<180
    original=(g.OUT/'source-retry/original-native_eos_bridge.py').read_text()
    expected=original.replace("inp=np.load(OUT/(label+'-input.npz'));","inp=dict(np.load(OUT/(label+'-input.npz')));")
    expected=expected.replace('except (ValueError,IndexError,KeyError,OSError): pass','except (ValueError,IndexError,KeyError,OSError,StopIteration): pass')
    assert expected==(g.OUT/'source-retry/materialized_native_eos_bridge_v2.py').read_text()
    assert expected.count('time.monotonic()-begin<180')==2
    electron=np.load(OP/'electrons/replacement.npy');new=dict(np.load(OP/'new-electrons-captured.npz'))
    baseline=dict(np.load(OP/'baseline-type1-captured.npz'))
    assert np.array_equal(new['used'],electron) and np.array_equal(new['outputs'],baseline['outputs'])
    assert not read(OP/'derivatives/result.json')['passed']
    assert read(OP/'cubic/certificate/piecewise-result.json')['all_12828_saved_affine_paths_resolved']
    restored=read(OP/'tables/result.json');legacy=read(OP/'tables/mixed-term/result.json')
    assert not restored['passed'] and not legacy['passed']
    assert restored['maximum_derivative_score']<1e-3 and restored['maximum_two_step_difference']<1e-3
    assert legacy['maximum_native_value_relative_change']<1e-4
    bound=read(OP/'tables/mixed-term/continuous-correction-bound.json');assert bound['rectangles']==29400
    save('audit.json',dict(classification='Proven',passed=True,verified_input_bindings=bindings,
        source_arrays_and_corrections_verified=True,prior_failed_gates_preserved=True,
        physical_EOS_or_opacity_certified=False,full_GR_evolution=False))
    x=sp.symbols('x');full=x*x/4;left=full-x*(1-x)/4;right=full-(x-1)*(2-x)/4
    assert left.subs(x,1)==right.subs(x,1)==sp.Rational(1,4)
    dl=sp.diff(left,x).subs(x,1);dr=sp.diff(right,x).subs(x,1)
    assert dl==sp.Rational(3,4) and dr==sp.Rational(1,4)
    save('mixed-term-C1-counterexample.json',dict(classification='Proven',passed=True,
        polynomial='x^2*y^2 on x grids [0,1,2], y grids [0,1], evaluated at y=1/2',
        omitted_value_continuous_at_x1=True,left_derivative=str(dl),right_derivative=str(dr),
        restored_derivative='1/2',
        conclusion='Omitting the mixed spline term can destroy C1 continuity even when the full spline exactly represents a smooth polynomial.'))
    print('TRANSPORT CHECKPOINT AUDIT',len(bindings),'input bindings; source replay and failed gates preserved',flush=True)


def notes():
    assert read(OUT/'audit.json')['passed'] and not (g.OUT/'transport-note-history.json').exists()
    prior=read(g.OUT/'progress-bindings.json');history={}
    for rel,digest in prior['sha256'].items():
        path=g.ROOT/rel;assert g.c.sha(path)==digest
        target=g.OUT/'transport-previous-notes'/path.name;target.parent.mkdir(exist_ok=True);shutil.copy2(path,target)
        history[rel]=dict(snapshot=target.relative_to(g.ROOT).as_posix(),sha256=digest)
    g.save('transport-note-history.json',history)
    shutil.copy2(g.OUT/'progress-bindings.json',g.OUT/'before-transport-progress-bindings.json')
    text='''

## Request 33 원천 연결과 불투명도 미분 검증

분류: Counterexample candidate. 새 EOS의 5,735구역 반응 입력을 실제 native 함수에 연결했다. 첫 계측은 180초 제한으로 실패했다. 불변 압축 입력을 한 번만 읽도록 바꾸고, 작성 중인 프로파일의 불완전 헤더를 기다리도록 수정한 별도 계측 코드에서 같은 제한을 유지했다. 이전 원천의 비트 일치를 확인한 뒤 새 입력 계측을 약 15.3초에 마쳤다. 새 전자 입력의 불투명도 개입은 반환값 변화 0이었다. 실제 온도 범위에서 이 입력을 사용하는 Compton 경로가 비활성이었다.

분류: Counterexample candidate. 원래 불투명도 미분은 사전 0.1% 기준을 실패했다. 실제 단조 3차 보간 19,242회를 읽어, 값에 적용된 선택 분기의 미분과 미분값을 별도로 보간한 결과의 불일치를 확인했다. 선택된 식을 직접 미분해 복사·전도 결합에 전달했다. 이어 원본 계수를 읽은 이중정밀도 계산에서 전도 혼합항 누락을 확인했다. 정수 z36th에 1/36을 대입하는 보관 코드가 그 원인이었다.

분류: Proven. 저장된 네 보간값의 선형 변화 경로 12,828개를 정확한 유리수 연산으로 판정했다. 12,826개는 단일 식으로 보증했고 나머지 두 경로는 분할했으며 두 비영 미분 점프를 보존했다. 누락된 전도 혼합항은 전체 29,400개 저장 격자 구간에서 값·두 미분의 조건부 상계를 얻었다. 실제 표 값의 물리 오차나 전체 EOS 미분의 보증으로 확대하지 않는다.

분류: Counterexample candidate. 혼합항을 복구한 계산의 유한 미분 최대 점수는 5.404e-4로 통과했지만, 원래 native 값과의 최대 차이 6.323e-4는 원래 값 대조 기준 1e-4를 실패했다. 혼합항을 0으로 둔 별도 재현은 native 값 차이 6.251e-6으로 값 대조를 통과했으나 미분 최대 점수 1.0624e-3으로 전체 판정은 실패했다. 이 원래 판정들을 유지한다. 아직 생성되지 않은 새 GR 상태에 적용할 별도 검증 기준과 모델 해시를 고정했다.

분류: Conjectural. 물리 EOS·표 데이터 오차·전체 연속 미분, 자체 수송·유체·계량·대기의 시간 진화, 실제 비영 scalar 구동·전하 및 완전한 관측 추론은 남는다. 동적 chi 후보의 계산 기반을 보완했으며 신규 관측량 판정과 최종 투고본 갱신은 보류한다. 근거: [단계 33 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).
'''
    for rel in prior['sha256']:
        if rel.startswith('docs/'):
            with (g.ROOT/rel).open('ab') as stream: stream.write(text.encode('utf-8'))
    report=g.ROOT/'notes/REQUEST33_DIRECT_EOS_GR_KO.md'
    extra='''

분류: Counterexample candidate. 새 EOS 고정 상태의 반응 원천에서 이전 모델과의 최대 상대 변화는 반응 벡터 9.998e-6, 가열 6.028e-9, 중성미자 원천 7.656e-13이었다. 입력과 상태를 고정한 EOS 교체 비교이며 시간 진화나 관측 검출이 아니다. 근거는 `reference-source.json`과 `source-retry/`다.

분류: Proven. 전도 혼합항의 상계는 `delta log10 conductivity`, 그 밀도·온도 로그 미분에 각각 0.0474504, 0.949008, 0.949008이다. 이는 29,400개 저장 격자와 볼록 전하 보간을 덮는 보수적 전역 상계다. 정확한 유리수는 `opacity/tables/mixed-term/continuous-correction-bound.json`에 보존한다. 이 값들을 실제 표 데이터나 전체 항성의 오차 막대로 사용하지 않는다.

분류: Proven. 매끄러운 다항식 x^2 y^2조차 혼합항을 버리면 지정 격자 경계에서 좌·우 미분이 3/4, 1/4로 달라질 수 있다. 복구한 다항식의 미분은 양쪽 모두 1/2이다. 따라서 작은 유한차분 간격만으로 원래 구현의 전구간 C1 미분을 보증할 수 없다. 근거: `transport-checkpoint/mixed-term-C1-counterexample.json`.

분류: Counterexample candidate. 최초 이상 쌍입자 상계 영역 beta<=0.005는 실제 표본의 최대 beta=0.00544555를 포함하지 않았다. 이 영역 실패를 보존했다. 별도로 선언한 beta<=0.006, eta>=-50 영역에는 5,735구역의 구역당 아홉 평가점이 포함됐다.

분류: Proven. 확대된 이상기체 영역의 총·순 자유입자수 상대 차이는 1e-100 미만이다. 상호작용 플라스마 및 미분의 물리 오차 인증은 아니다. 근거: `opacity/electrons/expanded-pair-bound.json`.

분류: Counterexample candidate. 원래 항성 호출은 필요한 표만 늦게 읽으므로, 이중정밀도 조성 좌표가 요구하는 인접 표가 누락되어 첫 재계산이 중단됐다. 별도 합성 호출로 원래 X/Z 격자의 124개 복사 표를 모두 읽었고, 앞서 저장한 모든 계수·격자는 비트 단위로 같았다. 합성 호출은 표를 읽기 위한 대조이며 새로운 항성 조성이나 진화가 아니다. 이전 미완료 로그와 수용 기준은 보존했다.
'''
    with report.open('ab') as stream: stream.write((text+extra).encode('utf-8'))
    paper_path=g.ROOT/'paper/revision-manifest.json';paper=read(paper_path)
    paper['request33_transport_checkpoint']=dict(report='notes/REQUEST33_DIRECT_EOS_GR_KO.md',
        status='새 EOS 원천 입력 연결과 원래 원천 비트 재현 완료. 불투명도 보간·혼합항 결함을 분리하고 조건부 구간 보증을 추가했다. 원래 실패 판정과 전체 물리 폐쇄 미완료를 유지한다.')
    paper['sha256'].update({rel:g.c.sha(g.ROOT/rel) for rel in prior['sha256'] if rel.startswith('docs/')})
    paper_path.write_text(json.dumps(paper,ensure_ascii=False,indent=2)+'\n')
    prior['sha256']={rel:g.c.sha(g.ROOT/rel) for rel in prior['sha256']};g.save('progress-bindings.json',prior)


def freeze():
    assert read(OUT/'audit.json')['passed'] and not (OUT/'manifest.json').exists()
    paths=[p for p in OP.rglob('*') if p.is_file()]
    paths += [p for p in (g.OUT/'source-retry').rglob('*') if p.is_file()]
    paths += [g.ROOT/'verification'/(name+'.py') for name in CODE]
    paths += [g.OUT/n for n in ['reference-state.npz','reference-auxiliary.npy','reference-auxiliary.json','reference-source.json']]
    for prefix in ['reference-cached-','cached-identity-v2-']: paths += [p for p in g.OUT.glob(prefix+'*') if p.is_file()]
    paths += [OUT/'audit.json',OUT/'mixed-term-C1-counterexample.json',g.OUT/'gr-opacity/plan.json']
    save('manifest.json',dict(classification='Proven',scope='Completed source/opacity checkpoint only; Request33 GR and the full physical/observation goal remain open.',
        sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in sorted(set(paths))},
        physical_EOS_or_opacity_certified=False,full_GR_evolution=False,submission_updated=False))
    verify()


def verify():
    manifest=read(OUT/'manifest.json')
    for rel,digest in manifest['sha256'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    print('PASS TRANSPORT CHECKPOINT',len(manifest['sha256']),'artifact/code SHA; original failures preserved; physical closure open',flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
