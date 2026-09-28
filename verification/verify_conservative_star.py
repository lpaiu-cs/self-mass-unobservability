"""Request31 packaging: current numerical verdicts and immutable prior evidence."""
from pathlib import Path
import hashlib, json, subprocess, sys
import conservative_star as s


def finalize():
    assert not (s.OUT/'manifest.json').exists()
    comp=json.loads((s.OUT/'composition-derivatives.json').read_text())
    evo=json.loads((s.OUT/'evolution.json').read_text())
    tangent=json.loads((s.OUT/'source-tangent.json').read_text())
    assert comp['complete'] and len(comp['rows'])==25 and evo['complete']
    assert json.loads((s.OUT/'tail-input-complete.json').read_text())['complete']
    partial=json.loads((s.OUT/'composition-partial-control.json').read_text())
    assert partial['complete'] and partial['passed']
    chain=json.loads((s.OUT/'h2-chain-tangent.json').read_text())
    assert chain['all_zone_finite_chain_passed'] and chain['burning_direct_passed']
    assert tangent['H2_chain_candidate_used'] and not tangent['H2_small_step_candidate_used']
    assert [r['steps'] for r in evo['rows']]==[1,2,4]
    for name in ['compensation-control','source-control','refined-he3','refined-lnT','structural-split-control',
                 'enthalpy-chart-control','enthalpy-coordinate-predictor','zero-source-projection']:
        result=json.loads((s.OUT/(name+'.json')).read_text())
        assert result.get('passed',result.get('finite_step_passed')),name
    chart=json.loads((s.OUT/'residual-chart-audit.json').read_text())
    assert chart['complete'] and chart['passed'] and len(chart['rows'])==5735
    assert json.loads((s.OUT/'structure-arithmetic-audit.json').read_text())['baseline_bitwise_equal']
    trial=json.loads((s.OUT/'residual-corrected-trial.json').read_text())
    assert trial['new_forced_quasistatic_trial_completed'] and not trial['time_refinement_completed']
    energy=all(r['energy_passed'] for r in evo['rows'])
    finest_energy=evo['rows'][-1]['energy_passed']
    refined=evo['rows'][-1]['refinement']['passed']
    s.save('gates.json',dict(classification='Counterexample candidate',
        all_25_composition_directions_recomputed=True,independent_initial_tangent_used=True,
        actual_new_forced_GR_paths_recomputed=True,all_path_energy_gates_passed=energy,
        finest_path_energy_gate_passed=finest_energy,
        finest_composition_temperature_time_gate_passed=refined,
        residual_corrected_single_trial_completed=True,
        residual_corrected_single_trial_energy_passed=trial['energy_passed'],
        residual_corrected_time_refinement_completed=False,
        fixed_flux=True,fixed_initial_tangent=True,
        physical_common_EOS_certified=False,continuous_EOS_derivative_error_certified=False,
        all_zone_continuous_source_derivative_certified=False,whole_star_physical_transport=False,
        full_GR_thermal_fluid_metric_evolution=False,scalar_actual_drive_charge_map=False,
        complete_nonlinear_observation=False,final_submission_package_updated=False))
    wording=f"새 1·2·4단계 고정 광도 GR 경로의 에너지 기준 일괄 판정은 {'통과' if energy else '실패'}다. 가장 미세한 4단계의 에너지 기준은 {'통과' if finest_energy else '실패'}, 마지막 조성·온도 시간 대조는 {'통과' if refined else '실패'}다. 거친 경로의 실패는 미세 경로의 판정과 구분하여 보존한다."
    wording+=f" 별도 구조 피드백 잔차 보정 1단계의 에너지 점수는 {trial['energy_score']:.8g}로 {'통과' if trial['energy_passed'] else '실패'}했으며, 그 보정 방법의 시간 수렴은 아직 검증하지 않았다."
    additions={
        'model-definition': '분류: Counterexample candidate. 같은 EOS·반응 함수를 조성 25방향과 온도·밀도에서 재평가했다. 직접 미분한 초기 연산자는 명시적인 고정 적분 연산자로 사용하며 각 시간 단계의 비선형 원천과 EOS는 새로 계산한다. '+wording,
        'observable-targets': '분류: Counterexample candidate. 반환된 조성 편미분과 EOS·조성 보조량을 포함한 전체 미분을 구분했다. 실제 값과 대조하지 않은 Jacobian을 관측 오차 보증으로 사용하지 않는다.\n\n분류: Conjectural. 새 상태의 자체 수송·대기·유체·계량, 실제 scalar 구동과 전체 비선형 관측 추론은 여전히 연결해야 한다.',
        'adiabatic-limit': '분류: Proven. 바리온 질량 좌표에서 구대칭 준정적 Fourier 식은 면적 제곱·적색편이·온도 기울기로 표현되며 지정 면 이산식은 일정한 적색편이 온도에서 영유속을 준다.\n\n분류: Conjectural. 이 항등식이나 유한 시간 대조는 orbital 완화시간, 빠른 모드 제거와 물리 오차의 연속 보증을 대신하지 않는다.',
        'nonadiabatic-regime': '분류: Counterexample candidate. 누적 조성 증분, 정지질량 변화, EOS 엔탈피, 직접 적분한 총 중성미자 손실과 광도 발산을 연결한 별도 원천 적분을 실행했다. 최대 Be7 시간 차이는 구조 재조정의 반응 피드백으로 약 99.995% 재현됐다. 기준 압력 엔탈피 좌표와 전체 함수의 잔차를 사용한 별도 보정도 실행했다. '+wording+' 물리적 수송과 외부 구동의 응답 pole을 이 결과에서 추정하지 않는다.',
        'failure-ledger-dynamic-chi': '분류: Counterexample candidate. He3의 원래 간격 실패와 온도 혼합 경계를 넘은 중성미자 미분 실패를 보존하고 별도 작은 간격을 대조했다. 가열 구역 마스크 밖의 H2 실패도 별도 기록한다. 음의 광도에서 1 erg/s로 분모가 잘리는 native 출력은 그대로 대류 비율로 해석할 수 없다. 중심 질량 간격을 생략한 수송 출력 재현 실패와 실제 경계를 반영한 대조를 모두 보존한다.\n\n분류: Conjectural. 물리 EOS·전구간 미분, 자체 열수송·대류·대기·완전 GR 진화, 비영 scalar 구동·전하와 전체 관측 추론의 완료는 아직 주장하지 않는다.'}
    history=json.loads((s.OUT/'historical-note-bindings.json').read_text())
    for name,paragraph in additions.items():
        rel='docs/'+name+'.md';p=s.ROOT/rel;snapshot=s.ROOT/history[rel]['snapshot']
        assert p.read_bytes()==snapshot.read_bytes(),rel
        with p.open('ab') as f:
            f.write(('\n\n## Request 31 조성 미분과 보존형 GR 재적분\n\n'+paragraph+
                '\n\n근거: [한글 보고서](../notes/REQUEST31_CONSERVATIVE_STAR_KO.md).\n').encode('utf-8'))
    p=s.ROOT/'paper/revision-manifest.json';paper=json.loads(p.read_text())
    paper['request31_supporting_note_update']=dict(evidence_manifest='outputs/conservative-star31/manifest.json',
        historical_notes='outputs/conservative-star31/historical-note-bindings.json',
        status=wording+' 물리 EOS·전체 GR·관측 폐쇄는 미완료다.',
        artifact_status='Request12 원고 PDF/ZIP은 역사 산출물로 보존한다.')
    for name in additions: paper['sha256']['docs/'+name+'.md']=s.c.sha(s.ROOT/'docs'/f'{name}.md')
    p.write_text(json.dumps(paper,ensure_ascii=False,indent=2)+'\n')
    note=s.ROOT/'notes/REQUEST31_CONSERVATIVE_STAR_KO.md';contents=note.read_text()
    contents=contents.replace('분류: Conjectural. 단계 31은 진행 중이다.',
        '분류: Counterexample candidate. 단계 31의 지정 경로와 별도 보정 1단계 계산을 마쳤다.')
    contents=contents.replace('별도 구조 피드백 보정 경로를 검사하고 있다.',
        '별도 구조 피드백 보정 1단계를 실행했으며, 이 방법의 시간 정밀화는 후속 작업으로 남아 있다.')
    contents=contents.replace('이 값을 이용한 별도 잔차 보정 경로를 실행 중이다.',
        '이 값을 이용한 별도 잔차 보정 1단계를 실행했다.')
    contents=contents.replace('전체 계산과 감사가 끝나기 전에는 단계 31의 최종 manifest나 최종 투고 상태를 선언하지 않는다.',
        '단계 31의 지정 결과와 역사 산출물은 manifest 및 검증 코드로 결합한다. 최종 투고 상태나 전체 연구의 완료는 선언하지 않는다.')
    contents+='\n\n분류: Counterexample candidate. '+wording+' 전체 물리 EOS·미분 보증·자체 GR 진화·비영 구동·관측 연결의 목표는 계속 열려 있다.\n'
    note.write_text(contents)
    prior=json.loads((s.OLD/'provenance.json').read_text())
    s.save('provenance.json',dict(classification='Imported from prior work',
        previous_provenance_sha256=s.c.sha(s.OLD/'provenance.json'),
        executable_and_libraries=prior['executable_and_libraries'],
        original_native_outputs_modified=False,source_function='Request30 double-table correction, fresh seeded FreeEOS auxiliary inputs',
        optional_seven_input_partial_control='Composition moments and EOS inputs are recorded and restored in a disposable child. This is a separate partial derivative control, never a modification of returned results.'))
    paths=[p for p in s.OUT.rglob('*') if p.is_file() and p.name!='manifest.json']
    paths += [s.ROOT/n for n in history]
    paths += [s.ROOT/n for n in ['verification/conservative_star.py','verification/verify_conservative_star.py',
        'notes/REQUEST31_CONSERVATIVE_STAR_KO.md','.gitattributes']]
    s.save('manifest.json',dict(classification='Proven',sha256={p.relative_to(s.ROOT).as_posix():s.c.sha(p) for p in sorted(set(paths))}))


def verify():
    stages=['validated-variational','remaining-levers15','nbody-readout16','nonzero-drive17','thermal-wd18',
        'thermal-restart19','thermal-robustness20','gr-mass21','thermal-closure22','baryon-entropy23',
        'reactive-energy24','fresh-microphysics25','remaining-closure26','native-closure27','conservative-cell28',
        'common-eos29','closure-precision30','conservative-star31']
    history={f'outputs/{a}/manifest.json':f'outputs/{b}/historical-note-bindings.json' for a,b in zip(stages,stages[1:])};count=0
    for label in ['outputs/research-remediation/manifest.json',*history,'paper/revision-manifest.json','outputs/conservative-star31/manifest.json']:
        old=json.loads((s.ROOT/history[label]).read_text()) if label in history else {}
        for name,digest in json.loads((s.ROOT/label).read_text())['sha256'].items():
            path=s.ROOT/name
            if name in old:
                bind=old[name];path=s.ROOT/bind['snapshot'];assert digest==bind['sha256']
                if bind.get('historical_manifest'):
                    before=json.loads(path.read_text());after=json.loads((s.ROOT/name).read_text())
                    for k,value in before.items():
                        if k!='sha256': assert after[k]==value,k
                    for k,value in before['sha256'].items():
                        if k not in old: assert after['sha256'][k]==value,k
                elif not bind.get('code_revision'): assert (s.ROOT/name).read_bytes().startswith(path.read_bytes())
            assert s.c.sha(path)==digest,(label,name);count+=1
    p29=json.loads((s.v.OLD/'provenance.json').read_text());p30=json.loads((s.OLD/'provenance.json').read_text())
    for path,digest in p30['executable_and_libraries'].items(): assert s.c.sha(Path(path))==digest,path
    for path,digest in p29['FreeEOS_source_sha256'].items(): assert s.c.sha(s.c.gr.SOURCE/path)==digest,path
    native=json.loads((s.ROOT/p29['unchanged_native_data_binding']).read_text())
    for path,digest in native['sha256'].items(): assert s.c.sha(Path(native['root'])/path)==digest,path
    for path,digest in p30['interval_sources'].items(): assert s.c.sha(Path(p30['interval_source_root'])/path)==digest,path
    gates=json.loads((s.OUT/'gates.json').read_text());evo=json.loads((s.OUT/'evolution.json').read_text())
    assert gates['all_path_energy_gates_passed']==all(r['energy_passed'] for r in evo['rows'])
    assert gates['finest_path_energy_gate_passed']==evo['rows'][-1]['energy_passed']
    assert gates['finest_composition_temperature_time_gate_passed']==evo['rows'][-1]['refinement']['passed']
    trial=json.loads((s.OUT/'residual-corrected-trial.json').read_text())
    assert gates['residual_corrected_single_trial_energy_passed']==trial['energy_passed']
    assert not gates['residual_corrected_time_refinement_completed'] and not trial['time_refinement_completed']
    for key in ['physical_common_EOS_certified','continuous_EOS_derivative_error_certified',
        'all_zone_continuous_source_derivative_certified','whole_star_physical_transport','full_GR_thermal_fluid_metric_evolution',
        'scalar_actual_drive_charge_map','complete_nonlinear_observation','final_submission_package_updated']:
        assert gates[key] is False,key
    print('PASS',count,'artifact/history SHA;',len(native['sha256']),'native data SHA;',len(p29['FreeEOS_source_sha256']),
        'FreeEOS source SHA;',len(p30['interval_sources']),'interval source SHA',flush=True)


def git_blobs():
    expected=dict(json.loads((s.OUT/'manifest.json').read_text())['sha256'])
    expected['outputs/conservative-star31/manifest.json']=s.c.sha(s.OUT/'manifest.json')
    with subprocess.Popen(['git','cat-file','--batch'],cwd=s.ROOT,stdin=subprocess.PIPE,stdout=subprocess.PIPE) as proc:
        for rel,digest in expected.items():
            proc.stdin.write(('HEAD:'+rel+'\n').encode());proc.stdin.flush()
            header=proc.stdout.readline().split();assert len(header)==3 and header[1]==b'blob',(rel,header)
            remaining=int(header[2]);actual=hashlib.sha256()
            while remaining:
                chunk=proc.stdout.read(min(remaining,1048576));assert chunk
                actual.update(chunk);remaining-=len(chunk)
            assert proc.stdout.read(1)==b'\n' and actual.hexdigest()==digest,rel
        proc.stdin.close();assert proc.wait(timeout=10)==0
    print('PASS',len(expected),'raw Git blobs at',subprocess.check_output(['git','rev-parse','HEAD'],cwd=s.ROOT,text=True).strip(),flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
