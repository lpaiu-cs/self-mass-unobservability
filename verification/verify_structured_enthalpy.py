"""Request32 verdicts and immutable history; numerical failure is retained."""
from pathlib import Path
import hashlib, json, subprocess, sys
import structured_enthalpy as e
import audit_structured_enthalpy as audit
import direct_ion_eos as direct


def progress():
    assert not (e.OUT/'progress-note-bindings.json').exists()
    text='\n\n## Request 32 직접 원소 EOS 중간 검증\n\n분류: Counterexample candidate. Li·Be·B·F를 직접 포함한 24원소 EOS 후보를 별도로 만들었다. 5,735개 기존 원소 호환성 대조는 비트 단위로 일치하고 6개 핵종의 희박 고온 전자수 대조도 통과했다. 실제 전체 조성의 미분 검사는 한 구역에서 실패했으며, 기존 EOS에서도 같은 실패를 재현했다. 실제 전자 교환 인자의 Cody–Thacher η=1 근사 경계 통과를 확인했다. 작은 차분의 국소 통과를 연속 오차 보증으로 취급하지 않는다.\n\n분류: Conjectural. 직접 Fermi 적분 대조, 동위원소 및 물리 EOS 오차, 자체 GR 수송·유체·계량 진화와 실제 구동·관측 연결은 계속 수행해야 한다.\n\n근거: [단계 32 보고서](../notes/REQUEST32_STRUCTURED_ENTHALPY_KO.md).\n'
    history=json.loads((e.OUT/'historical-note-bindings.json').read_text());bindings={}
    for rel,row in history.items():
        if not rel.startswith('docs/'): continue
        p=e.ROOT/rel;assert p.read_bytes()==(e.ROOT/row['snapshot']).read_bytes()
        with p.open('ab') as f: f.write(text.encode('utf-8'))
        bindings[rel]=e.c.sha(p)
    p=e.ROOT/'paper/revision-manifest.json';paper=json.loads(p.read_text())
    paper['request32_in_progress']=dict(note='notes/REQUEST32_STRUCTURED_ENTHALPY_KO.md',
        status='직접 원소 EOS 대조와 근사 미분 경계 실패를 기록했다. 4단계 경로 및 최종 검증은 미완료다.')
    paper['sha256'].update(bindings);p.write_text(json.dumps(paper,ensure_ascii=False,indent=2)+'\n')
    e.save('progress-note-bindings.json',dict(classification='Proven',sha256=bindings))


def finalize():
    assert not (e.OUT/'manifest.json').exists()
    e.collect();audit.stages();audit.embedding_moment_control();audit.inverse_energy_audit()
    evo=json.loads((e.OUT/'evolution.json').read_text())
    checks=json.loads((e.OUT/'stage-audit.json').read_text())
    assert evo['complete'] and checks['complete']
    for name in ['control','initial-control','reference-control','transport-control',
                 'implicit-transport-control','embedding-moment-control','strict-entropy-control']:
        assert json.loads((e.OUT/(name+'.json')).read_text())['passed'],name
    strict=json.loads((e.OUT/'strict-entropy-projection.json').read_text());assert strict['completed']
    direct.provenance_control()
    for name in ['supported-control','pure-control','ct-interval','isotope-no-go','fermi-polylog-control','fermi-tail-control']:
        assert json.loads((direct.OUT/(name+'.json')).read_text())['passed'],name
    eos_checks={name:json.loads((direct.OUT/(name+'mixture-control.json')).read_text())
        for name in ['', 'integral-', 'tight-integral-','full-integral-']}
    primitive=json.loads((direct.OUT/'fermi-interval/result.json').read_text())
    assert primitive['complete'] and primitive['passed'] and len(primitive['rows'])==12
    assert json.loads((direct.OUT/'fermi-interval/rule-control.json').read_text())['passed']
    rows=evo['rows'];energy=all(r['energy_passed'] for r in rows)
    finest=rows[-1];refined=finest['refinement']['passed']
    closed=['physical_common_EOS_certified','continuous_EOS_derivative_error_certified',
        'all_zone_continuous_source_derivative_certified','whole_star_physical_transport',
        'full_GR_thermal_fluid_metric_evolution','scalar_actual_drive_charge_map',
        'complete_nonlinear_observation','final_submission_package_updated']
    e.save('gates.json',dict(classification='Counterexample candidate',
        direct_predictor_and_corrector_GR_paths_completed=True,
        all_path_energy_gates_passed=energy,finest_path_energy_gate_passed=finest['energy_passed'],
        finest_composition_temperature_time_gate_passed=refined,
        strict_entropy_projection_completed=True,strict_entropy_projection_energy_passed=strict['energy_passed'],
        direct_element_compatibility_passed=True,direct_element_dilute_electron_count_passed=True,
        direct_element_CT_derivative_gate_passed=eos_checks['']['passed'],
        direct_fermi_default_derivative_gate_passed=eos_checks['integral-']['passed'],
        direct_fermi_tight_derivative_gate_passed=eos_checks['tight-integral-']['passed'],
        full_fermi_exchange_derivative_gate_passed=eos_checks['full-integral-']['passed'],
        nonrelativistic_Fermi_upper_tail_bound_proved=True,
        registered_Fermi_primitive_point_errors_certified=True,
        CT_piece_noncontinuity_proved=True,isotope_group_exactness_obstruction_proved=True,
        strict_root_full_time_refinement_completed=False,
        fixed_flux=True,fixed_initial_preconditioner=True,**{k:False for k in closed}))
    verdict='; '.join(f"{r['steps']}단계 에너지 {r['energy_score']:.8g} ({'통과' if r['energy_passed'] else '실패'})" for r in rows)
    verdict+=f". 마지막 조성 시간 점수 {finest['refinement']['composition_score']:.8g}, 로그 온도 차이 {finest['refinement']['lnT_difference']:.8g}로 시간 대조는 {'통과' if refined else '실패'}다."
    verdict+=f" 별도 에너지 단위 엔트로피 역산 GR 재투영의 점수는 {strict['energy_score']:.8g}로 {'통과' if strict['energy_passed'] else '실패'}다. 이 별도 재투영은 강화한 역산법의 전체 시간 재적분이 아니다."
    additions={
        'model-definition':'분류: Counterexample candidate. 고정 기준 압력의 엔탈피를 열 좌표로 사용하고, 매 단계의 예측 상태와 잔차 보정 상태에 각각 새 EOS 역산·GR 접합을 수행했다. 이전 실패 종점을 새 예측 상태로 대체하지 않았다. '+verdict,
        'observable-targets':'분류: Counterexample candidate. 같은 최종 시점의 보정 GR 상태와 누적 조성 증분을 독립 대조했다. '+verdict+'\n\n분류: Conjectural. 자체 수송과 방사 대기, 비영 scalar 구동·전하, 전체 비선형 관측 추론은 계속 연결해야 한다.',
        'adiabatic-limit':'분류: Proven. 양의 원소 치환 가중치가 이온 수·평균 전하·전하 제곱합을 모두 보존하면 치환 전하의 분산이 0이어야 한다. 따라서 지원하지 않는 전하를 다른 전하들로 치환하는 현재 사상은 완전 이온화 한계의 Coulomb 입력까지 정확히 보존할 수 없다. 조성만의 에너지 기준 이동으로 이 한계를 제거할 수 없다.',
        'nonadiabatic-regime':'분류: Counterexample candidate. 구조 피드백을 명시적으로 다시 평가한 엔탈피 보정법으로 1·2·4단계를 적분했다. '+verdict+' 고정된 광도 면 입력과 초기 행렬을 사용했으며, 이 결과는 물리 수송이나 관측 완화 pole의 검출이 아니다.',
        'failure-ledger-dynamic-chi':'분류: Counterexample candidate. '+verdict+' 이전 단계의 시간·에너지 실패는 그대로 보존한다. 내부 확산의 고정 계수 대조에서 명시적 양성 시간 간격은 약 4.64e-5초다. 별도 암시적 밴드 풀이 대조가 통과해도 실제 표면 열손실·대류·대기·비선형 EOS와 GR 결합의 통과를 뜻하지 않는다.\n\n분류: Proven. 수·전하 보존 원소 치환의 Coulomb 전하 제곱합 불일치를 명시했다.\n\n분류: Conjectural. 물리 EOS 및 연속 오차, 전체 GR 진화와 실제 구동·관측 폐쇄는 아직 완료되지 않았다.'}
    eos_verdict='분류: Counterexample candidate. 직접 원소 후보의 CT 미분 대조는 '+('통과' if eos_checks['']['passed'] else '실패')+'다. 직접 Fermi 적분의 원래 정밀도 및 강화 정밀도 전구역 대조는 각각 '+('통과' if eos_checks['integral-']['passed'] else '실패')+', '+('통과' if eos_checks['tight-integral-']['passed'] else '실패')+'다. 이 유한 대조는 물리·연속 EOS 오차 보증이 아니다.'
    eos_verdict+=' 교환 항의 CT 호출까지 직접 적분으로 바꾼 후보의 전구역 대조는 '+('통과' if eos_checks['full-integral-']['passed'] else '실패')+'다.'
    eos_proof='분류: Proven. 보관 CT 근사식의 η=1 함수 및 일차 미분 점프가 0을 배제함을 계수 구간으로 보였다. 또한 서로 다른 동위원소 질량의 부분 이온화 평형은 조성만의 엔트로피 이동으로 원소 단위 EOS에서 정확히 복원되지 않는다.'
    eos_proof+=' 새 공통 Fermi 평가의 η=1·4 지정 원시 값·미분 12개는 구간 구적으로 반폭 1e-11 이하에 감쌌고, 해당 정규화 출력의 절대오차 상계는 모두 1e-9 이하다. 이는 지정 원시 함수의 점별 보증이며 전체 EOS나 GR 오차 보증이 아니다.'
    additions={name:paragraph+'\n\n'+eos_verdict+'\n\n'+eos_proof for name,paragraph in additions.items()}
    history=json.loads((e.OUT/'historical-note-bindings.json').read_text())
    for name,paragraph in additions.items():
        rel='docs/'+name+'.md';p=e.ROOT/rel;snapshot=e.ROOT/history[rel]['snapshot']
        if (e.OUT/'progress-note-bindings.json').exists():
            expected=json.loads((e.OUT/'progress-note-bindings.json').read_text())['sha256'][rel]
            assert e.c.sha(p)==expected and p.read_bytes().startswith(snapshot.read_bytes()),rel
        else: assert p.read_bytes()==snapshot.read_bytes(),rel
        with p.open('ab') as f:
            f.write(('\n\n## Request 32 구조 피드백 엔탈피 시간 대조\n\n'+paragraph+
                '\n\n근거: [한글 보고서](../notes/REQUEST32_STRUCTURED_ENTHALPY_KO.md).\n').encode('utf-8'))
    p=e.ROOT/'paper/revision-manifest.json';paper=json.loads(p.read_text())
    paper.pop('request32_in_progress',None)
    paper['request32_supporting_note_update']=dict(evidence_manifest='outputs/structured-enthalpy32/manifest.json',
        historical_notes='outputs/structured-enthalpy32/historical-note-bindings.json',status=verdict,
        artifact_status='Request12 원고 PDF/ZIP은 역사 산출물로 보존한다. 물리 EOS·전체 GR·관측 폐쇄 및 최종 투고 상태는 미완료다.')
    for name in additions: paper['sha256']['docs/'+name+'.md']=e.c.sha(e.ROOT/'docs'/f'{name}.md')
    p.write_text(json.dumps(paper,ensure_ascii=False,indent=2)+'\n')
    note=e.ROOT/'notes/REQUEST32_STRUCTURED_ENTHALPY_KO.md'
    with note.open('a') as f: f.write('\n\n분류: Counterexample candidate. 지정 보정 경로를 모두 계산하고 독립 검증했다. '+verdict+'\n\n'+eos_verdict+'\n\n'+eos_proof+'\n')
    prior=json.loads((e.OLD/'provenance.json').read_text())
    e.save('provenance.json',dict(classification='Imported from prior work',
        previous_provenance_sha256=e.c.sha(e.OLD/'provenance.json'),
        executable_and_libraries=prior['executable_and_libraries'],original_native_outputs_modified=False,
        direct_ion_provenance='outputs/structured-enthalpy32/direct-ions/provenance.json',
        direct_ion_provenance_sha256=e.c.sha(direct.OUT/'provenance.json'),
        source_function='Request30 four-weak double-table correction and seeded FreeEOS auxiliaries; new direct GR predictor and corrector stages.'))
    paths=[p for p in e.OUT.rglob('*') if p.is_file() and p.name!='manifest.json']
    paths += [e.ROOT/n for n in history]
    paths += [e.ROOT/n for n in ['verification/structured_enthalpy.py','verification/audit_structured_enthalpy.py',
        'verification/verify_structured_enthalpy.py','verification/direct_ion_eos.py','verification/fermi_interval.py',
        'notes/REQUEST32_STRUCTURED_ENTHALPY_KO.md','.gitattributes']]
    e.save('manifest.json',dict(classification='Proven',sha256={p.relative_to(e.ROOT).as_posix():e.c.sha(p) for p in sorted(set(paths))}))


def verify():
    stages=['validated-variational','remaining-levers15','nbody-readout16','nonzero-drive17','thermal-wd18',
        'thermal-restart19','thermal-robustness20','gr-mass21','thermal-closure22','baryon-entropy23',
        'reactive-energy24','fresh-microphysics25','remaining-closure26','native-closure27','conservative-cell28',
        'common-eos29','closure-precision30','conservative-star31','structured-enthalpy32']
    history={f'outputs/{a}/manifest.json':f'outputs/{b}/historical-note-bindings.json' for a,b in zip(stages,stages[1:])};count=0
    for label in ['outputs/research-remediation/manifest.json',*history,'paper/revision-manifest.json','outputs/structured-enthalpy32/manifest.json']:
        old=json.loads((e.ROOT/history[label]).read_text()) if label in history else {}
        for name,digest in json.loads((e.ROOT/label).read_text())['sha256'].items():
            path=e.ROOT/name
            if name in old:
                bind=old[name];path=e.ROOT/bind['snapshot'];assert digest==bind['sha256']
                if bind.get('historical_manifest'):
                    before=json.loads(path.read_text());after=json.loads((e.ROOT/name).read_text())
                    for k,value in before.items():
                        if k!='sha256': assert after[k]==value,k
                    for k,value in before['sha256'].items():
                        if k not in old: assert after['sha256'][k]==value,k
                elif not bind.get('code_revision'): assert (e.ROOT/name).read_bytes().startswith(path.read_bytes())
            assert e.c.sha(path)==digest,(label,name);count+=1
    p29=json.loads((e.s.v.OLD/'provenance.json').read_text());p30=json.loads((e.s.OLD/'provenance.json').read_text())
    for path,digest in p30['executable_and_libraries'].items(): assert e.c.sha(Path(path))==digest,path
    for path,digest in p29['FreeEOS_source_sha256'].items(): assert e.c.sha(e.c.gr.SOURCE/path)==digest,path
    native=json.loads((e.ROOT/p29['unchanged_native_data_binding']).read_text())
    for path,digest in native['sha256'].items(): assert e.c.sha(Path(native['root'])/path)==digest,path
    for path,digest in p30['interval_sources'].items(): assert e.c.sha(Path(p30['interval_source_root'])/path)==digest,path
    gates=json.loads((e.OUT/'gates.json').read_text());evo=json.loads((e.OUT/'evolution.json').read_text())
    assert evo['complete'] and [r['steps'] for r in evo['rows']]==[1,2,4]
    assert gates['all_path_energy_gates_passed']==all(r['energy_passed'] for r in evo['rows'])
    assert gates['finest_path_energy_gate_passed']==evo['rows'][-1]['energy_passed']
    assert gates['finest_composition_temperature_time_gate_passed']==evo['rows'][-1]['refinement']['passed']
    strict=json.loads((e.OUT/'strict-entropy-projection.json').read_text())
    assert strict['completed'] and gates['strict_entropy_projection_energy_passed']==strict['energy_passed']
    assert not gates['strict_root_full_time_refinement_completed'] and not strict['full_strict_root_time_path_recomputed']
    dp=json.loads((direct.OUT/'provenance.json').read_text())
    for path,digest in dp['libraries'].items(): assert e.c.sha(Path(path))==digest,path
    for root,files in dp['sources'].items():
        for rel,digest in files.items(): assert e.c.sha(Path(root)/rel)==digest,rel
    for prefix,key in [('', 'direct_element_CT_derivative_gate_passed'),
            ('integral-','direct_fermi_default_derivative_gate_passed'),
            ('tight-integral-','direct_fermi_tight_derivative_gate_passed'),
            ('full-integral-','full_fermi_exchange_derivative_gate_passed')]:
        assert gates[key]==json.loads((direct.OUT/(prefix+'mixture-control.json')).read_text())['passed']
    assert gates['CT_piece_noncontinuity_proved'] and json.loads((direct.OUT/'ct-interval.json').read_text())['passed']
    assert gates['isotope_group_exactness_obstruction_proved'] and json.loads((direct.OUT/'isotope-no-go.json').read_text())['passed']
    assert gates['nonrelativistic_Fermi_upper_tail_bound_proved'] and json.loads((direct.OUT/'fermi-tail-control.json').read_text())['passed']
    primitive=json.loads((direct.OUT/'fermi-interval/result.json').read_text())
    assert gates['registered_Fermi_primitive_point_errors_certified'] and primitive['complete'] and primitive['passed']
    assert len(primitive['rows'])==12 and json.loads((direct.OUT/'fermi-interval/rule-control.json').read_text())['passed']
    for key in ['physical_common_EOS_certified','continuous_EOS_derivative_error_certified',
        'all_zone_continuous_source_derivative_certified','whole_star_physical_transport','full_GR_thermal_fluid_metric_evolution',
        'scalar_actual_drive_charge_map','complete_nonlinear_observation','final_submission_package_updated']:
        assert gates[key] is False,key
    print('PASS',count,'artifact/history SHA;',len(native['sha256']),'native data SHA;',len(p29['FreeEOS_source_sha256']),
        'FreeEOS source SHA;',len(p30['interval_sources']),'interval source SHA',flush=True)
    print('PASS',len(dp['libraries']),'new EOS library SHA;',sum(map(len,dp['sources'].values())),'candidate source SHA',flush=True)


def git_blobs():
    expected=dict(json.loads((e.OUT/'manifest.json').read_text())['sha256'])
    expected['outputs/structured-enthalpy32/manifest.json']=e.c.sha(e.OUT/'manifest.json')
    with subprocess.Popen(['git','cat-file','--batch'],cwd=e.ROOT,stdin=subprocess.PIPE,stdout=subprocess.PIPE) as proc:
        for rel,digest in expected.items():
            proc.stdin.write(('HEAD:'+rel+'\n').encode());proc.stdin.flush()
            header=proc.stdout.readline().split();assert len(header)==3 and header[1]==b'blob',(rel,header)
            remaining=int(header[2]);actual=hashlib.sha256()
            while remaining:
                chunk=proc.stdout.read(min(remaining,1048576));assert chunk
                actual.update(chunk);remaining-=len(chunk)
            assert proc.stdout.read(1)==b'\n' and actual.hexdigest()==digest,rel
        proc.stdin.close();assert proc.wait(timeout=10)==0
    print('PASS',len(expected),'raw Git blobs at',subprocess.check_output(['git','rev-parse','HEAD'],cwd=e.ROOT,text=True).strip(),flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
