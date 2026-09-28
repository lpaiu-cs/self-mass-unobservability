"""Close Phase113 evidence without repeating any coupled integration."""
from pathlib import Path
import hashlib
import json
import shutil
import sys
import numpy as np
ROOT=Path('/mnt/e/lab/self-mass-unobservability');OUT=ROOT/'outputs/direct-eos-gr33/def-native-cold-coupling'
sys.path.insert(0,str(ROOT/'verification'))
import def_native_cold_coupling as t
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n')


def archive():
    # Recover historical producer bytes only when their original digest agrees.
    src=(OUT/'scaled-producer.py').read_text();small=src[:src.index('def diagnostic_build():')].rstrip()+'\n\n\n'+src[src.index("if __name__=='__main__'"):]
    for line in ['import shutil\n','import subprocess\n','from types import FunctionType\n',"COLD=prior.optical.ex.old.cold\n","CACHE=COLD.CACHE.parent/'native-cold-coupling'\n"]:small=small.replace(line,'')
    candidates=[small,small.replace('\n\n\n','\n\n')]
    wanted=json.loads((OUT/'plan.json').read_text())['bindings'][str(ROOT/'verification/def_native_cold_coupling.py')]
    found=False
    for s in candidates:
        for newline in ['\n','\r\n']:
            b=s.replace('\n',newline).encode()
            if hashlib.sha256(b).hexdigest()==wanted:(OUT/'initial-producer.py').write_bytes(b);found=True
    write(OUT/'initial-producer-recovery.json',dict(matched=found,original_sha256=wanted))
    for name in ['operands113.py','lowest113.py']:
        p=Path('/mnt/c/Users/lpaiu/AppData/Local/Temp')/name
        if p.exists():shutil.copyfile(p,OUT/name)
    # Diagnostic source was overwritten by the separate repair candidate.
    for name in ['mod_excitation.f90','excitation_pi.f90','excitation_sum.f90']:
        s=(t.COLD.CACHE/name).read_text()
        if name=='excitation_pi.f90':
            anchor='                 mu(ion) = exparg*qstar(nmin_s,izqstar)'
            s=s.replace(anchor,"                 if(tl.lt.log(240._fp_kind)) write(*,*) 'COLD_EXP',iz,nmin_s,izqstar,exparg,qstar(nmin_s,izqstar),c2t*bion(ion),plop(ion),qh2plus\n"+anchor)
        (OUT/('diagnostic-'+name)).write_text(s)
    mapping={}
    for phase in ['diagnostic','scaled','logarithmic']:
        receipt=json.loads((OUT/(phase+'-build.json')).read_text())
        for name,digest in receipt['bindings'].items():
            p=Path(name)
            if p.suffix=='.f90':p=OUT/(phase+'-'+p.name)
            assert sha(p)==digest,(phase,p,digest,sha(p));mapping[name+'@'+phase]=dict(path=str(p),sha256=digest)
    write(OUT/'build-binding-audit.json',dict(classification='Counterexample candidate',passed=True,bindings=mapping))
    symbolic=json.loads((OUT/'symbolic.json').read_text());assert symbolic['passed']
    checks={p.stem:json.loads(p.read_text()) for p in OUT.glob('*.json') if p.name in ['diagnostic.json','partition-operands.json','partition-diagnostic.json','scaled-controls.json','lowest-probe.json','logarithmic-controls.json','support.json','support-controls.json','audit.json']}
    calls=sum(v['native_calls'] for v in checks.values())
    write(OUT/'resource-accounting.json',dict(classification='Counterexample candidate',measured_native_calls=calls,
        measured_breakdown={k:v['native_calls'] for k,v in checks.items()},additional_constructor_call_upper_bound=36,
        reason='Conservative upper bound for table/model constructors in support controls, one actual continuation, readout attempts and endpoint audit; no full trajectory repetition.',
        actual_continuation_seconds=json.loads((OUT/'resumed-896-128.json').read_text())['seconds'],continuation_budget_seconds=190,
        fluid_macroscopic_steps=24,reused_complete_steps=104,reused_coarse_paths=2,CPU_threads=1,peak_memory_measured=False,
        note='The original signed Infinity/NaN diagnostic values and wrong-slot partition inspection are preserved as failed diagnostic output, not accepted finite constitutive data.'))
    # The retained native/EOS and same support are tested. Compile revised
    # sources without running a second integration merely for reporting.
    compile((ROOT/'verification/def_native_cold_coupling.py').read_text(),str(ROOT/'verification/def_native_cold_coupling.py'),'exec')
    compile((ROOT/'verification/verify_native_cold_coupling.py').read_text(),str(ROOT/'verification/verify_native_cold_coupling.py'),'exec')
    write(OUT/'completion-boundary.json',dict(classification='Counterexample candidate',
        native_cold_overflow_repaired=True,original_failed_state_recovered=True,original_coupled_fine_horizon_completed=True,
        segment_conservation_passed=True,actual_native_endpoint_audit_passed=True,
        nominal_original_space_gates_passed=True,strict_restart_readout_audit_passed=False,overall_phase_verdict=False,
        remaining='The corrected nominal full-history spatial differences are0.46253percent trace and0.67978percent mass. The newly registered restart overlap1e-8 and independent readout1e-10 checks fail narrowly/at recovery roundoff. Retain their failures. Recovering internal energy a second time introduces a subtractive trace error; next use conserved-energy trace identity and an explicit recovery error budget, without repeating evolution.',
        full_interior_mechanics=False,frequency_convergence=False,full_microphysics=False,full_GR_scalar_feedback=False,final_charge_solved=False,full_goal_complete=False))
    print(json.dumps(dict(archive=True,initial_producer_found=found,native_calls=calls)))


def documents():
    paragraphs={
      'model-definition':'분류: Counterexample candidate. 차가운 H2+ 분배함수의 중간 overflow와 극미량 H 광학 판독을 공통 지수·로그 표현으로 수정했다. 같은 native 재고·준위로 실패 셀의170.318K 근을 얻고 필요한 cold 상태136개를 연결해 실제 양방향 광자·움직이는 대기의 원 fine 경로를 완주했다. 온도·화학 재고를 절단하지 않았고 기존 warm 보간기와 두 coarse 경로를 재사용했다. 전체 내부 밀도·계량 고정과 미시 반응 범위는 그대로다.',
      'observable-targets':'분류: Counterexample candidate. 수정 입력으로 원128단계 결합 진화를 완주했다. 기준값 오류 수정 뒤 전체 구간의 명목 공간 차이는 대기 trace0.462530%, 외부 질량0.679781%로 원2% 안이다. 별도로 등록한 엄격한 재시작 중첩·독립 읽기 검사는 실패다. 이 결과를 최종 scalar 전하나 관측 완화 검출로 승격하지 않는다. 다음 원천 읽기는 보존 변수와 회복 오차를 사용한다.',
      'adiabatic-limit':'분류: Proven. 공통 스케일을 원 분배함수와 그 원래 편미분에 적용하면 정확한 실수 산술에서 물리적 곱과 미분 비가 보존된다. 같은 물질 trace는 보존 변수로 tau−S*v−3p라 쓸 수 있다. 분류: Counterexample candidate. 이번 로그 표현과 실제 결합 완주는 저온 계산 경계를 해결했으나 정적 계수 흡수에 관한 단열 no-go 경계를 바꾸지 않는다. 지정 초기 과도응답의 체적 적분은 관측 독립성을 보증하지 않는다.',
      'nonadiabatic-regime':'분류: Counterexample candidate. 원fine 경로가 저온 native EOS 수리 뒤 실제 광자·유체 방정식으로3.4344311179ms 전체 구간을 완주했다. 추가24단계106.16초이며 새 구간 수지·native 끝점 대조를 통과했다. 명목 공간 차이는 원 기준 안이지만 엄격한 읽기 실패를 보존한다. 내부 반경·주파수·전체 기계 진화·GR/scalar 되먹임과 최종 질량 정규화는 미완료다.',
      'failure-ledger-dynamic-chi':'분류: Counterexample candidate. 기존109단계 실패는 H2+ 분배함수를 보상 Boltzmann 인자보다 먼저 계산한 중간 overflow로 추적했다. 분배함수·원 미분의 공통 스케일과 H 광학 로그값으로 수정해 실제 원경로를 완주했다. 최초 컴파일 인자 누락·80K 광학 판독 실패·잘못된 진단 배열 행·원 EOS 실패는 보존한다. 재시작 보고 기준값 오류는 실제 상태를 바꾸지 않고 수정했다. 중첩 trace1.3335e-8은 등록1e-8을, 독립 끝점6.4195e-7은1e-10을 넘어 전체 형식 판정은false다. 보존 변수 trace와 회복 오차의 직접 전달이 다음 최소 작업이며 진화를 자동 반복하지 않는다.',
      'dynamic-charge-completion':'분류: Counterexample candidate. native_cold_overflow_repaired, original_failed_state_recovered, original_fine_coupled_horizon_completed, segment_conservation_passed, actual_native_endpoint_audit_passed, nominal_original_space_gates_passed는true다. strict_restart_readout_audit_passed, overall_phase_verdict, full_interior_mechanics, frequency_convergence, full_microphysics, full_GR_scalar_feedback, final_charge_solved, full_goal_complete는false다. 실제 EOS 병목 해결은loophole progress, 공통 스케일·보존 trace 항등식은theorem progress다. 다음은 저장 원천의 보존 표현과 회복 오차를 실제 전하에 전달하는 일이다.'}
    prefixes={}
    for name,paragraph in paragraphs.items():
        p=ROOT/'docs'/f'{name}.md';old=p.read_bytes();assert '## 단계113' not in old.decode()
        prefixes[str(p.relative_to(ROOT))]=dict(bytes=len(old),sha256=hashlib.sha256(old).hexdigest())
        addition='\n\n## 단계113 — 저온 EOS 병목을 고친 실제 결합 완주\n\n'+paragraph+'\n\n상세: [단계113 보고](../notes/REQUEST113_NATIVE_COLD_COUPLING_KO.md).\n'
        p.write_bytes(old+addition.encode())
    write(OUT/'documentation-prefixes.json',prefixes)
    write(OUT/'osk-write-status.json',dict(detail='sha256:5785317ba70da2a975ef1542337f2109547bb2adf9258c137df942faa10271bf',root='sha256:bed1cadb6d592442060b00a50d66a5723f2fc871ee58f2712d2ad70f94ce94fe',scope='sha256:b6ac316e5eb174415f088b4aa17403ecab02c5fe95097ea6c9f9f2d363bdc621',scope_chars=1488,body_and_scope_written=True))


def verify():
    import re
    for p,record in json.loads((OUT/'documentation-prefixes.json').read_text()).items():assert hashlib.sha256((ROOT/p).read_bytes()[:record['bytes']]).hexdigest()==record['sha256']
    a=np.load(OUT/'first-readout-resumed-896-128.npz');b=np.load(OUT/'resumed-896-128.npz')
    for k in a.files:
        if k not in ['bulk_trace','atmosphere_trace']:assert np.array_equal(a[k],b[k]),k
    result=json.loads((OUT/'result.json').read_text());assert result['continuation']['passed'] and not result['passed']
    assert result['continuation']['completed_steps']==128 and json.loads((OUT/'audit.json').read_text())['passed']
    assert max(result['space_comparison'].values())<.02 and max(result['overlap'].values())>1e-8
    digests={sha(p):str(p.relative_to(ROOT)) for p in OUT.glob('*.py')}
    wanted={data['source_sha256']:p.stem for p in OUT.glob('*.json') if (data:=json.loads(p.read_text())).get('source_sha256')}
    for source in list(OUT.glob('*producer.py'))+[ROOT/'verification/def_native_cold_coupling.py']:
        s=source.read_text()
        for match in re.finditer(r'^def \w+\(',s,re.M):
            prefix=s[:match.start()].rstrip()
            for blanks in [1,2,3]:
                candidate=prefix+'\n'*blanks+"if __name__=='__main__':globals()[sys.argv[1]]()\n"
                for newline in ['\n','\r\n']:
                    raw=candidate.replace('\n',newline).encode();h=hashlib.sha256(raw).hexdigest()
                    if h in wanted and h not in digests:
                        archive=OUT/(wanted[h]+'-producer.py');archive.write_bytes(raw);digests[h]=str(archive.relative_to(ROOT))
    provenance={}
    for p in OUT.glob('*.json'):
        data=json.loads(p.read_text());h=data.get('source_sha256')
        if h:provenance[p.name]=dict(sha256=h,archive=digests.get(h),inherited_prior_source=h==sha(ROOT/'verification/def_native_two_way_atmosphere.py'))
    write(OUT/'producer-binding-audit.json',dict(classification='Counterexample candidate',bindings=provenance,
        limitation='Nested resumed result source_sha256 originally reported inherited Phase112 source; resume-plan and resume-producer bind the actual wrapper. Reporting correction retains original bytes and no dynamics were rerun.'))
    write(OUT/'closure-audit.json',dict(classification='Counterexample candidate',passed=True,original_document_prefixes_preserved=True,
        state_arrays_unchanged_by_readout_correction=True,cold_bottleneck_fixed_in_original_trajectory=True,failed_strict_readout_verdict_preserved=True,
        nominal_space_comparison=result['space_comparison'],full_goal_complete=False))
    print(json.dumps(dict(closure_checks=True,strict_phase_verdict=result['passed'])))


def manifest():
    paths=[ROOT/'verification/def_native_cold_coupling.py',ROOT/'verification/verify_native_cold_coupling.py',ROOT/'notes/REQUEST113_NATIVE_COLD_COUPLING_KO.md']
    paths += [ROOT/p for p in json.loads((OUT/'documentation-prefixes.json').read_text())]
    paths += [p for p in OUT.rglob('*') if p.is_file() and '__pycache__' not in str(p)]
    files={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in sorted(paths)}
    target=ROOT/'outputs/direct-eos-gr33/native-cold-coupling-manifest.json'
    write(target,dict(classification='Counterexample candidate',checkpoint='1297eba71',actual_original_coupled_horizon_completed=True,strict_readout_passed=False,final_charge_solved=False,full_goal_complete=False,files=files))
    master=ROOT/'paper/revision-manifest.json';d=json.loads(master.read_text());d['sha256'].update(files);d['sha256'][str(target.relative_to(ROOT))]=sha(target)
    d['native_cold_coupling']=dict(classification='Counterexample candidate',report='notes/REQUEST113_NATIVE_COLD_COUPLING_KO.md',manifest=str(target.relative_to(ROOT)),cold_failure_fixed_in_actual_evolution=True,strict_readout_passed=False,full_goal_complete=False)
    write(master,d)
    for p,h in d['sha256'].items():assert sha(ROOT/p)==h,(p,'Master hash mismatch')
    print(json.dumps(dict(phase_files=len(files),master_hashes=len(d['sha256']),all_current_hashes_passed=True)))


if __name__=='__main__':globals()[sys.argv[1]]()
