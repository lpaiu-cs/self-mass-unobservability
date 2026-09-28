"""Evidence closure; never executes another fluid/photon trajectory."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
ROOT=Path('/mnt/e/lab/self-mass-unobservability');OUT=ROOT/'outputs/direct-eos-gr33/def-native-coupled-charge'
sys.path.insert(0,str(ROOT/'verification'))
import def_native_coupled_charge as task
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):Path(p).write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n')


def prepare():
    # Match both response components to exactly the same observer reference.
    d=dict(np.load(OUT/'first-exterior-energy.npz'));source=np.load(OUT/'source-896-128.npz');q=np.load(OUT/'wave-896-128.npz')['normalized_direct'];M=float(source['M_cm']);alpha0=-float(source['K_cm'])/M
    eps=task.G*d['arrived_energy_erg']/(task.C**4*M);eps0=task.G*d['frozen_initial_energy_erg']/(task.C**4*M)
    direct=q-q[0];d['direct_plus_mass_normalization']=(direct+alpha0*eps)/(1-eps)
    d['excess_over_frozen_photon_normalization']=d['direct_plus_mass_normalization']-alpha0*eps0/(1-eps0)
    np.savez_compressed(OUT/'exterior-energy.npz',**d);result=json.loads((OUT/'first-exterior.json').read_text())
    result.update(endpoint_direct_plus_mass_normalization=float(d['direct_plus_mass_normalization'][-1]),endpoint_excess_over_frozen_photon_normalization=float(d['excess_over_frozen_photon_normalization'][-1]),observer_reference_aligned=True)
    write(OUT/'exterior.json',result)
    write(OUT/'observer-reference.json',dict(classification='Counterexample candidate',
        correction='The photon mass term was already relative to observer u=0; subtract the material direct value at the same u=0 as well. This is a stored-array algebraic correction, not a rerun.',
        direct_at_observer_zero=float(q[0]),relative_to_endpoint=float(abs(q[0]/q[-1])),gates_unchanged=True,first_result_preserved=True))
    write(OUT/'first-audit-failure.json',dict(classification='Counterexample candidate',passed=False,
        failure='Exact equality of two differently padded pairwise floating-point baryon sums failed.',
        correction='The intended invariant is that fixed deep-density columns contain zero baryon change. Verify those columns directly. Manufactured and source-seed acceptance thresholds are unchanged.'))
    paragraphs={
      'model-definition':'분류: Counterexample candidate. 실제 결합 이력의 보존 변수에서 물질 trace와 반경 응력·광자 에너지/압력을 읽었다.1,100km 내부와 움직이는 대기의 저장 원천을 동일 지연 적분에 넣고, 출사 광자를 실제 진공 null geodesic으로 외부 질량 분모에 연결했다. 내부 밀도·계량 고정, 안쪽 바리온 점 보상, 각도 빈 내부 상수 재구성은 명시적 조건이다. 전체 물리 전하로 승격하지 않는다.',
      'observable-targets':'분류: Counterexample candidate. 실제 결합 원천의 지연 직접 성분2.2368e-27과 출사 광자의 질량 정규화4.1530e-27의 합은 선언 모형에서6.3897e-27로 양수다. 초기 출사 점유를 고정한 질량 성분만 뺀 값도5.7626e-27이다. 완전한 정적 비교 항성, nuisance 흡수, 전체 GR 전하 또는 궤도 관측 검출이 아니다. 내부 기계·계량·반경/주파수·각도 오차를 닫아야 한다.',
      'adiabatic-limit':'분류: Proven. 보존 변수의 비정지 trace는 tau−S*v−3p다. 정규화 응답은 전하 변화와 복사 질량 감소를 같은 분모에 넣은(delta_alpha_direct+alpha0*epsilon)/(1-epsilon)로 쓸 수 있다. 분류: Counterexample candidate. 이번 양의 조건부 과도응답은 정적 EFT 흡수에 대한 단열 no-go를 바꾸지 않는다. 광행 지연을 궤도 완화시간 검출로 해석하지 않는다.',
      'nonadiabatic-regime':'분류: Counterexample candidate. 보존 변수 판독의 초기 추정 민감도1.367e-11, 지연 직접 성분의 시간 대조0.19965percent·저장17/9시각 대조0.71036percent가 등록 기준 안이다. 외부 광자의 각도별 도달 시간을 적용했다. 마지막 각도 빈의 거의 반경 방향 부분만 이 짧은 구간에 도착한다. 구적 수렴을 물리적 빈 내부 분포의 인증으로 세지 않는다.',
      'failure-ledger-dynamic-chi':'분류: Counterexample candidate. 단계113의 엄격한 읽기 실패는 보존한다. 단계114의 새 보존 변수 판독과 지정 성분 대조는 통과했다. 첫 읽기의 구적점 의존 관측 시계, 감사의 서로 다른 부동소수 합산 순서 비트 동일성 검사, 두 성분의 미세한 관측 기준값 차이는 원 출력과 함께 보존하고 읽기에서 수정했다. 유체·열·광자 이력을 재실행하지 않았다. 깊은 내부 기계적 바리온 재배치와 계량/scalar 원천이 다음 직접 병목이다.',
      'dynamic-charge-completion':'분류: Counterexample candidate. conservative_saved_source_readout_passed, actual_coupled_retarded_direct_component, vacuum_photon_arrival_mass_component, declared_component_controls_passed는true다. full_interior_mechanics, radial_frequency_certification, full_angular_shape_error_bound, full_GR_scalar_feedback, complete_exterior_mass_bookkeeping, final_charge_solved, full_goal_complete는false다. 지정 성분 합은양수이며loophole progress이나 목표는active다. 다음은 깊은 내부의 실제 압력·광자 힘과 바리온/운동량 응답을 원천에 연결한다.'}
    prefixes={}
    for name,body in paragraphs.items():
        p=ROOT/'docs'/f'{name}.md';data=p.read_bytes();assert '## 단계114' not in data.decode()
        prefixes[str(p.relative_to(ROOT))]=dict(bytes=len(data),sha256=hashlib.sha256(data).hexdigest())
        p.write_bytes(data+('\n\n## 단계114 — 실제 결합의 지연 전하와 광자 질량 분모\n\n'+body+'\n\n상세: [단계114 보고](../notes/REQUEST114_NATIVE_COUPLED_CHARGE_KO.md).\n').encode())
    write(OUT/'documentation-prefixes.json',prefixes)
    write(OUT/'osk-write-status.json',dict(detail='sha256:e2e13e2d3c62d5825e5ba823affab1be5657506f10959ef47cda8513cdcb9dea',root='sha256:bc9ebbdd2dad7a7497614a76683beb8d3ff49195fbb36ae94694dfe7b09febac',scope='sha256:1401a4a487feb2d9cefde9a489d86114a2123bbf67112d87656276841aa7bd73',scope_chars=1472))
    write(OUT/'resource-accounting.json',dict(classification='Counterexample candidate',new_fluid_steps=0,new_native_states=0,
        native_constructor_calls_upper_bound=10,source_seconds=9.095874516999999,first_readout_seconds=4.732050233999999,
        corrected_clock_readout_seconds=4.908551183,exterior_seconds=4.7365411800000015,independent_exterior_audit_seconds=.163538089,
        note='Source45s/readout45s/exterior30s preregistered caps. Constructor bound also includes initial read-only geometry inspection; no new native EOS bank. WSL startup and hashing excluded; peak memory not measured.'))


def verify():
    for p,row in json.loads((OUT/'documentation-prefixes.json').read_text()).items():assert hashlib.sha256((ROOT/p).read_bytes()[:row['bytes']]).hexdigest()==row['sha256']
    for name in ['sources','result','audit','exterior','exterior-audit','symbolic']:assert json.loads((OUT/(name+'.json')).read_text())['passed'],name
    p=json.loads((OUT/'plan.json').read_text());bindings=p['bindings'].copy();original=str(ROOT/'verification/def_native_coupled_charge.py')
    assert sha(OUT/'first-producer.py')==bindings.pop(original)
    for name,h in bindings.items():assert sha(name)==h,name
    write(OUT/'closure-audit.json',dict(classification='Counterexample candidate',passed=True,old_document_prefixes_preserved=True,
        actual_input_bindings_verified=True,original_failed_strict_readout_preserved=not json.loads((task.cold.OUT/'result.json').read_text())['passed'],
        new_sources_and_declared_direct_ray_controls_passed=True,full_goal_complete=False))
    for name in ['def_native_coupled_charge.py','verify_native_coupled_charge.py']:compile((ROOT/'verification'/name).read_text(),name,'exec')
    print(json.dumps(dict(closure_passed=True,full_goal_complete=False)))


def manifest():
    paths=[ROOT/'verification/def_native_coupled_charge.py',ROOT/'verification/verify_native_coupled_charge.py',ROOT/'notes/REQUEST114_NATIVE_COUPLED_CHARGE_KO.md']
    paths += [ROOT/p for p in json.loads((OUT/'documentation-prefixes.json').read_text())]
    paths += [p for p in OUT.rglob('*') if p.is_file() and '__pycache__' not in str(p)]
    files={str(p.relative_to(ROOT)):sha(p) for p in sorted(paths)};target=ROOT/'outputs/direct-eos-gr33/native-coupled-charge-manifest.json'
    write(target,dict(classification='Counterexample candidate',checkpoint='a19c21e12',declared_component_checks_passed=True,final_charge_solved=False,full_goal_complete=False,files=files))
    master=ROOT/'paper/revision-manifest.json';d=json.loads(master.read_text());d['sha256'].update(files);d['sha256'][str(target.relative_to(ROOT))]=sha(target)
    d['native_coupled_charge']=dict(classification='Counterexample candidate',report='notes/REQUEST114_NATIVE_COUPLED_CHARGE_KO.md',manifest=str(target.relative_to(ROOT)),declared_direct_and_mass_components_computed=True,final_charge_solved=False,full_goal_complete=False)
    write(master,d)
    for p,h in d['sha256'].items():assert sha(ROOT/p)==h,p
    print(json.dumps(dict(phase_files=len(files),master_hashes=len(d['sha256']),all_current_hashes_passed=True)))


if __name__=='__main__':globals()[sys.argv[1]]()
