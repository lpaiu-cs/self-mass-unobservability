"""Request30 evidence packaging and byte-level history checks."""
from pathlib import Path
import hashlib, json, subprocess, sys
import closure_precision as v


def finalize():
    assert not (v.OUT/'manifest.json').exists()
    inverse=json.loads((v.OUT/'cold-inverse.json').read_text())
    unified=json.loads((v.OUT/'unified-source-audit.json').read_text())
    interval=json.loads((v.OUT/'weak-interval-audit.json').read_text())
    v.save('gates.json',dict(classification='Counterexample candidate',
        exact_forcing_audit_mismatch_explained=True,
        seeded_EOS_repeatability=json.loads((v.OUT/'cold-EOS-test.json').read_text())['all_repeatable'],
        saved_endpoint_equations_resolved=all(r['passed'] for r in inverse['rows']),
        original_failed_local_cell_replayed=json.loads((v.OUT/'failed-cell-replay.json').read_text())['outcome']['passed'],
        finite_initial_EOS_checks=json.loads((v.OUT/'full-EOS-audit.json').read_text())['passed'],
        native_weak_derivative_zero_route_confirmed=True,
        declared_corrected_source_finite_derivative_checks=unified['passed'],
        fixed_auxiliary_four_weak_derivative_interval=interval['passed'],
        physical_common_EOS_certified=False,continuous_EOS_derivative_error_certified=False,
        original_native_derivative_certified=False,full_composition_Jacobian_certified=False,
        updated_GR_trajectory_recomputed=False,composition_time_gate_resolved=False,
        whole_star_physical_transport=False,full_GR_thermal_fluid_metric_evolution=False,
        scalar_actual_drive_charge_map=False,complete_nonlinear_observation=False,
        final_submission_package_updated=False))
    additions={
        'model-definition':'분류: Counterexample candidate. 같은 광도 입력과 고정 초기값 FreeEOS로 57,350개 저장 원천 끝점 방정식을 기존 허용량 안에서 재역산했다. 미지원 원소의 물리 EOS 인증과 전체 GR 재진화는 별도다. 실제 약반응 표의 선형 보간·열·중성미자 기여를 독립 재구성한 배정밀도 원천 함수를 정의했다.',
        'observable-targets':'분류: Counterexample candidate. 실제 반응 계산기의 weaklib 경로가 온도·밀도 미분을 0으로 지정함을 끝점 대조로 확인했다. 보정된 원천 함수의 전체 초기 상태 유한 미분 대조가 통과했다.\n\n분류: Conjectural. 이 원천 함수에서 새 GR 시간 경로·질량 정규화 scalar 읽기·실제 구동·관측 모형을 다시 계산해야 한다. 이전 조건부 scalar 수치를 새 미시물리 모형의 결과로 소급하지 않는다.',
        'adiabatic-limit':'분류: Proven. 선형 혼합 반응률의 온도 미분에는 가중치 미분 항이 필요하다. 고정 표·혼합 분기에서 지수 일차식의 3차 미분으로 중심 차분 오차를 제한할 수 있다.\n\n분류: Conjectural. 반응 미분을 수정해도 전체 항성의 빠른 모드 제거·orbital 완화시간·비영 읽기 잔여량이 자동으로 보장되지 않는다.',
        'nonadiabatic-regime':'분류: Counterexample candidate. 같은 배정밀도 약반응 함수와 EOS 보조 입력 대조에서 896개 핵 가열 구역의 벡터·열 미분 기준 1e-3을 모두 통과했다.\n\n분류: Proven. 지정된 네 약반응의 고정 밀도·자유전자·조성 모형에 한해 h=5e-5의 미분 계산 및 차분 절단 오차를 구간 보증했다. 정규화 상계 최대는 4.59802e-11 미만이다. 전체 EOS·반응·관측의 연속 미분 보증은 포함하지 않는다.',
        'failure-ledger-dynamic-chi':'분류: Counterexample candidate. Request29의 큰 독립 역산 잔차에는 exp(nu)**2와 exp(2nu)의 광도 차감 차이가 주로 기여했다. 동일 강제 입력에서는 최대 잔차/허용량이 1.71875 이하로 줄었고, 별도 고정 초기값 역산은 최대 0.5625로 통과했다. 기존 실패 배열은 보존한다.\n\n분류: Counterexample candidate. 약반응 혼합 가중치만 수정한 벡터 기준은 실패했다. 실제 경로는 약반응 온도·밀도 미분을 모두 0으로 덮으며 Q 및 Qnu의 미분도 전달하지 않는다. 단정밀도 선형 보간을 실제 표로 재현하고, 별도 배정밀도 기여와 완전한 부분 미분으로 초기 유한 기준을 통과했다. 실패한 설정 끝점 대조와 중간 미분 재구성을 보존했다.\n\n분류: Conjectural. 전체 조성 Jacobian, 물리 EOS·연속 오차, 새 GR 시간 경로·자체 수송·유체·계량·대기, 비영 scalar 구동과 전체 관측 추론은 남는다.'}
    history=json.loads((v.OUT/'historical-note-bindings.json').read_text())
    for name,paragraph in additions.items():
        rel='docs/'+name+'.md';p=v.ROOT/rel;snap=v.ROOT/history[rel]['snapshot']
        assert p.read_bytes()==snap.read_bytes(),rel
        with p.open('ab') as f:
            f.write(('\n\n## Request 30 EOS 역산과 약반응 미분 수정\n\n'+paragraph+'\n\n근거: [한글 보고서](../notes/REQUEST30_CLOSURE_PRECISION_KO.md).\n').encode('utf-8'))
    p=v.ROOT/'paper/revision-manifest.json';paper=json.loads(p.read_text())
    paper['request30_supporting_note_update']=dict(evidence_manifest='outputs/closure-precision30/manifest.json',
        historical_notes='outputs/closure-precision30/historical-note-bindings.json',
        status='EOS 국소 역산 재현성·약반응 미분 경로 수정과 지정 약반응 구간 보증; 물리 EOS·새 GR 진화·관측 폐쇄 미완료',
        artifact_status='Request12 원고 PDF/ZIP은 역사 산출물로 보존한다.')
    for name in additions:
        rel='docs/'+name+'.md';paper['sha256'][rel]=v.c.sha(v.ROOT/rel)
    p.write_text(json.dumps(paper,ensure_ascii=False,indent=2)+'\n')
    prov=json.loads((v.OLD/'provenance.json').read_text())
    import mpmath
    deps=Path(mpmath.__file__).parent
    v.save('provenance.json',dict(classification='Imported from prior work',
        previous_provenance_sha256=v.c.sha(v.OLD/'provenance.json'),
        executable_and_libraries=prov['executable_and_libraries'],
        interval_source_root=str(deps),interval_sources={p.relative_to(deps).as_posix():v.c.sha(p) for p in deps.rglob('*.py')},
        original_native_outputs_modified=False,declared_new_numerical_function=True))
    paths=[p for p in v.OUT.rglob('*') if p.is_file() and p.name!='manifest.json']
    paths += [v.ROOT/n for n in history]
    paths += [v.ROOT/n for n in ['verification/closure_precision.py','verification/verify_closure_precision.py',
        'notes/REQUEST30_CLOSURE_PRECISION_KO.md','.gitattributes']]
    v.save('manifest.json',dict(classification='Proven',sha256={p.relative_to(v.ROOT).as_posix():v.c.sha(p) for p in sorted(set(paths))}))


def verify():
    stages=['validated-variational','remaining-levers15','nbody-readout16','nonzero-drive17','thermal-wd18',
        'thermal-restart19','thermal-robustness20','gr-mass21','thermal-closure22','baryon-entropy23',
        'reactive-energy24','fresh-microphysics25','remaining-closure26','native-closure27','conservative-cell28','common-eos29','closure-precision30']
    history={f'outputs/{a}/manifest.json':f'outputs/{b}/historical-note-bindings.json' for a,b in zip(stages,stages[1:])};count=0
    for label in ['outputs/research-remediation/manifest.json',*history,'paper/revision-manifest.json','outputs/closure-precision30/manifest.json']:
        old=json.loads((v.ROOT/history[label]).read_text()) if label in history else {}
        for name,digest in json.loads((v.ROOT/label).read_text())['sha256'].items():
            path=v.ROOT/name
            if name in old:
                bind=old[name];path=v.ROOT/bind['snapshot'];assert digest==bind['sha256']
                if bind.get('historical_manifest'):
                    before=json.loads(path.read_text());after=json.loads((v.ROOT/name).read_text())
                    for k,value in before.items():
                        if k!='sha256': assert after[k]==value,k
                    for k,value in before['sha256'].items():
                        if k not in old: assert after['sha256'][k]==value,k
                elif not bind.get('code_revision'): assert (v.ROOT/name).read_bytes().startswith(path.read_bytes())
            assert v.c.sha(path)==digest,(label,name);count+=1
    prior=json.loads((v.OLD/'provenance.json').read_text());prov=json.loads((v.OUT/'provenance.json').read_text())
    for path,digest in prov['executable_and_libraries'].items(): assert v.c.sha(Path(path))==digest,path
    for path,digest in prior['FreeEOS_source_sha256'].items(): assert v.c.sha(v.c.gr.SOURCE/path)==digest,path
    data=json.loads((v.ROOT/prior['unchanged_native_data_binding']).read_text())
    for path,digest in data['sha256'].items(): assert v.c.sha(Path(data['root'])/path)==digest,path
    for path,digest in prov['interval_sources'].items(): assert v.c.sha(Path(prov['interval_source_root'])/path)==digest,path
    gates=json.loads((v.OUT/'gates.json').read_text())
    for key in ['seeded_EOS_repeatability','saved_endpoint_equations_resolved','original_failed_local_cell_replayed',
        'finite_initial_EOS_checks','native_weak_derivative_zero_route_confirmed',
        'declared_corrected_source_finite_derivative_checks','fixed_auxiliary_four_weak_derivative_interval']: assert gates[key] is True,key
    for key in ['physical_common_EOS_certified','continuous_EOS_derivative_error_certified','original_native_derivative_certified',
        'full_composition_Jacobian_certified','updated_GR_trajectory_recomputed','composition_time_gate_resolved',
        'whole_star_physical_transport','full_GR_thermal_fluid_metric_evolution','scalar_actual_drive_charge_map',
        'complete_nonlinear_observation','final_submission_package_updated']: assert gates[key] is False,key
    print('PASS',count,'artifact/history SHA;',len(data['sha256']),'native data SHA;',len(prior['FreeEOS_source_sha256']),
        'FreeEOS source SHA;',len(prov['interval_sources']),'interval source SHA',flush=True)


def git_blobs():
    expected=dict(json.loads((v.OUT/'manifest.json').read_text())['sha256'])
    expected['outputs/closure-precision30/manifest.json']=v.c.sha(v.OUT/'manifest.json')
    with subprocess.Popen(['git','cat-file','--batch'],cwd=v.ROOT,stdin=subprocess.PIPE,stdout=subprocess.PIPE) as proc:
        for rel,digest in expected.items():
            proc.stdin.write(('HEAD:'+rel+'\n').encode());proc.stdin.flush()
            header=proc.stdout.readline().split();assert len(header)==3 and header[1]==b'blob',(rel,header)
            remaining=int(header[2]);actual=hashlib.sha256()
            while remaining:
                chunk=proc.stdout.read(min(remaining,1048576));assert chunk
                actual.update(chunk);remaining-=len(chunk)
            assert proc.stdout.read(1)==b'\n' and actual.hexdigest()==digest,rel
        proc.stdin.close();assert proc.wait(timeout=10)==0
    print('PASS',len(expected),'raw Git blobs at',subprocess.check_output(['git','rev-parse','HEAD'],cwd=v.ROOT,text=True).strip(),flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
