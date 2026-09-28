"""Bind the actual two-way GR response evidence and append project verdicts."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
ROOT=Path('/mnt/e/lab/self-mass-unobservability');OUT=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,v):Path(p).write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n')


def prepare():
    texts={
        'model-definition':'분류: Counterexample candidate. 실제 양방향 물질·광자의 에너지와 반경 압력을 선형 Einstein 질량 제약과 scalar 원천에 넣었다. 안쪽 광자 유속은 더 깊은 내부의 실제 에너지 감소로 반영했다. 추가 계량 변화의 좌표 재고·엔트로피 고정 체적 응답이라는 명시한 tangent 폐쇄를 사용하며, 물질·광자를 동적 계량에서 재진화한 것은 아니다.',
        'observable-targets':'분류: Counterexample candidate. 영역 내 계량 응력/질량 기여는 각각−1.4640e-32,−1.3986e-34다. 임의 비음수 외향 외부 광자 분포와 scalar 퍼텐셜 모든 반복을 허용한 추가 GR 상계 뒤에도, 고정된 실제 원천의 관측 직접 변화 하한은+1.14170e-27다. 이 조건부 하한은 원천 공간/주파수/내부 각도·초기 GR 및 전체 재진화 오차를 포함하지 않는다.',
        'adiabatic-limit':'분류: Proven. 추가 계량 변화에 좌표 재고·엔트로피·조성을 고정한 체적 응답을 포함하면 J_prime+(nu_prime+lambda_prime)J=4pi*r²*A⁴*E_forced다. 분류: Counterexample candidate. 이 보존 체적 tangent 식으로 실제 저장된 비단열 원천의 GR 응답을 계산했다. 실제 광자·물질 궤적을 계량과 함께 재계산하거나 전체 비단열 미시 폐쇄를 보인 것은 아니다.',
        'nonadiabatic-regime':'분류: Counterexample candidate. trace가0인 광자도 실제 E−Pr와 에너지 질량 제약을 통해 scalar를 구동한다. 해당 항을 포함한64/128 시간 차이는0.943percent다. 이번 짧은 과도응답에서 이 추가 선형 GR 항은 부호를 뒤집지 못하지만 궤도 구동·완화 식별 검출은 아니다.',
        'failure-ledger-dynamic-chi':'분류: Counterexample candidate. 저장 에너지의 질량 제약에 내부 경계 광자 공급의 반대 감소를 누락하지 않았다. 두 경계 유속과 실제 영역 에너지의 차이를 맞추는 상수 없이 측정하고 그 영향을 상계에 남겼다. 첫 상계 보고의 NumPy bool 직렬화 실패를 보존했다. 새 GR 구간은 원천 자체의 반경/내부 각도·초기 완전 Einstein 제약·자유 물질 경계·동적 GR 미완료를 해소하지 않는다.',
        'dynamic-charge-completion':'분류: Counterexample candidate. actual_matter_and_photon_metric_source_applied, actual_inner_energy_debit_applied, declared_all_orders_scalar_potential_bound, external_outward_angle_independent_GR_bound는true다. initial_full_Einstein_constraints_matched, free_mechanical_interface, source_continuum_error_certified, full_dynamic_GR_feedback, final_charge_solved, full_goal_complete는false다. 다음은 고정 원천의 작은 GR 항을 더 정밀하게 반복하기보다 실제 물질 접합과 원천의 공간/초기 제약을 닫는 일이다.'}
    prefixes={}
    for name,body in texts.items():
        p=ROOT/'docs'/(name+'.md');data=p.read_bytes();assert '## 단계117' not in data.decode()
        prefixes[str(p.relative_to(ROOT))]=dict(bytes=len(data),sha256=hashlib.sha256(data).hexdigest())
        p.write_bytes(data+('\n\n## 단계117 — 실제 물질·광자 원천의 GR 기여\n\n'+body+'\n\n상세: [단계117 보고](../notes/REQUEST117_NATIVE_FEEDBACK_GR_KO.md).\n').encode())
    write(OUT/'documentation-prefixes.json',prefixes)


def verify():
    for name,row in json.loads((OUT/'documentation-prefixes.json').read_text()).items():assert hashlib.sha256((ROOT/name).read_bytes()[:row['bytes']]).hexdigest()==row['sha256']
    for name in ['sources','result','bound','audit']:assert json.loads((OUT/(name+'.json')).read_text())['passed'],name
    bound=json.loads((OUT/'bound.json').read_text());audit=json.loads((OUT/'audit.json').read_text());z=np.load(OUT/'wave-128.npz')
    direct=z['components'][0]-z['components'][0,0];metric=z['components'][1:].sum(0);metric-=metric[0]
    relative=(max(abs(metric))+2*bound['uncomputed_GR_normalized_bound'])/max(abs(direct));assert relative<.02
    assert audit['arbitrary_outward_angular_mass_conditional_lower']>0
    for name in ['def_native_feedback_gr.py','verify_native_feedback_gr.py']:compile((ROOT/'verification'/name).read_text(),name,'exec')
    write(OUT/'closure-audit.json',dict(classification='Counterexample candidate',passed=True,
        observer_difference_all_GR_change_bound_over_direct=float(relative),two_event_uncertainty_not_single_endpoint=True,
        old_document_prefixes_unchanged=True,original_two_way_trajectories_unchanged=True,
        new_fluid_steps=0,new_photon_steps=0,full_physical_error_enclosed=False,final_charge_solved=False,full_goal_complete=False))
    print('CLOSURE PASSED',float(relative))


def manifest():
    paths=[ROOT/'verification'/name for name in ['def_native_feedback_gr.py','verify_native_feedback_gr.py']]
    paths += [ROOT/'notes/REQUEST117_NATIVE_FEEDBACK_GR_KO.md']
    paths += [ROOT/name for name in json.loads((OUT/'documentation-prefixes.json').read_text())]
    paths += [p for p in OUT.rglob('*') if p.is_file() and '__pycache__' not in str(p)]
    files={str(p.relative_to(ROOT)):sha(p) for p in sorted(paths)};target=ROOT/'outputs/direct-eos-gr33/native-feedback-gr-manifest.json'
    write(target,dict(classification='Counterexample candidate',checkpoint='e1862c3be',actual_GR_source_and_conditional_enclosure_passed=True,final_charge_solved=False,full_goal_complete=False,files=files))
    p=ROOT/'paper/revision-manifest.json';d=json.loads(p.read_text());d['sha256'].update(files);d['sha256'][str(target.relative_to(ROOT))]=sha(target)
    d['native_feedback_GR']=dict(classification='Counterexample candidate',report='notes/REQUEST117_NATIVE_FEEDBACK_GR_KO.md',manifest=str(target.relative_to(ROOT)),
        actual_GR_forcing_and_conditional_bound=True,conditional_positive_remnant=True,full_dynamic_GR_feedback=False,final_charge_solved=False,full_goal_complete=False)
    write(p,d)
    for name,h in d['sha256'].items():assert sha(ROOT/name)==h,name
    print(json.dumps(dict(phase_files=len(files),master_hashes=len(d['sha256']),all_current_hashes_verified=True)))


if __name__=='__main__':globals()[sys.argv[1]]()
