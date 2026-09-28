"""Freeze the actual two-way path evidence without rerunning its evolution."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
ROOT=Path('/mnt/e/lab/self-mass-unobservability');OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'verification'))
import def_native_interior_feedback as task

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,v):Path(p).write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n')


def prepare():
    texts={
        'model-definition':'분류: Counterexample candidate. 내부 면 바리온·운동량, 실제 압축/이류 총에너지, 비선형 수소·열·각도별 광자를 같은 시간 경로에 양방향 연결했다. 실제 밀도·유속은 EOS와 이동 광자 충돌에 되돌아가고 반대 운동량·일은 물질에 들어간다. 전체 GR 계량과 기계 경계의 자유 접합은 아직 없고, 내부 밀도·배경 재고 보정은 일차다.',
        'observable-targets':'분류: Counterexample candidate. 실제 양방향64/128 경로에서 직접 전하 끝점은+1.14392e-27, 광자 질량 정규화까지 더한 합은+5.29688e-27다. 시간 차이는 각각0.943percent와0.207percent다. 단계115의 약99percent 상쇄는 새 모형에서 유지되지 않는다. 이동 산란 끝점 규칙도 동시에 고쳤으므로 차이 전부를 단일 되먹임 항의 효과로 식별하지 않는다. 최종 물리 전하·정적 EFT 이탈 검출은 아니다.',
        'adiabatic-limit':'분류: Proven. 보존 에너지와 운동량을 사용하면 E−S*v−3p=D*sqrt(1−v²)*(cx*c²+u)−3p다. 고정 GR 총에너지 유속을 진화시킨 후 정지·운동 에너지를 빼서 내부에너지를 복원할 때 독립 p*dV 항을 다시 더하면 중복이다. 분류: Counterexample candidate. 이번 비단열 실제 결합은 이 표현을 사용하며 정적 계수 흡수의 기존 경계를 뒤집지는 않는다.',
        'nonadiabatic-regime':'분류: Counterexample candidate. 실제 밀도 변화와 유속을 광자·열·수소 교환에 되돌린 두 경로가 원 에너지·바리온·시간 문턱을 통과했다. 깊은 내부 비정지 trace 시간 차이는0.0918percent다. 영속도에서 유한 산란 누출이 남던 공통 주파수 경계를 양의 보간으로 고쳐야 실제 결합 양수성을 유지할 수 있었다. 이는 짧은 초기 과도응답이며 궤도 주기의 물리적 구동과 완화 검출은 아니다.',
        'failure-ledger-dynamic-chi':'분류: Counterexample candidate. 초기 Newton 부호, Wien 꼬리의 상대 미분 및 영속도 불연속 산란 누출 때문에 첫 결합 단계가 실패했다. 각각 원 상태와 실행 소스를 보존하고 근본 항을 수정한 뒤 실제64/128 경로를 완주했다. 원 산란 커널로 생산한 과거 이력은 재인증하지 않는다. 단계115의 거의 완전한 상쇄는 양방향 모형의 결론이 아니다. 새 경로의 전체 중성수소 누적 원장·자유 기계 접합·공간/각도·전체 GR·최종 전하 미완료를 남긴다.',
        'dynamic-charge-completion':'분류: Counterexample candidate. actual_two_way_interior_material_photon_feedback, original64_and128_horizons_completed, energy_and_baryon_gates_passed, direct_and_total_time_gates_passed는true다. complete_neutral_species_trajectory_audit, free_mechanical_interface, spatial_frequency_continuum_certified, full_angular_error_bound, full_GR_scalar_feedback, final_charge_solved, full_goal_complete는false다. 실제 연결 병목을 해소한loophole progress이며, 유한 모형의 성공과 전체 물리 폐쇄를 구분한다.'}
    prefix={}
    for name,body in texts.items():
        p=ROOT/'docs'/(name+'.md');data=p.read_bytes();assert '## 단계116' not in data.decode()
        prefix[str(p.relative_to(ROOT))]=dict(bytes=len(data),sha256=hashlib.sha256(data).hexdigest())
        p.write_bytes(data+('\n\n## 단계116 — 내부 물질과 광자의 실제 양방향 진화\n\n'+body+'\n\n상세: [단계116 보고](../notes/REQUEST116_NATIVE_INTERIOR_FEEDBACK_KO.md).\n').encode())
    write(OUT/'documentation-prefixes.json',prefix)
    production=json.loads((OUT/'production.json').read_text());bank=json.loads((OUT/'bank.json').read_text());audit=json.loads((OUT/'audit.json').read_text())
    write(OUT/'resources.json',dict(classification='Counterexample candidate',new_spectral_native_calls=bank['native_calls'],new_spectral_native_seconds=bank['seconds'],
        actual_endpoint_native_audit_calls=audit['native_calls'],independent_audit_seconds=audit['seconds'],
        completed_coarse_seconds=production['paths'][0]['seconds'],completed_fine_recovery_seconds=production['paths'][1]['seconds'],recovery_dispatch_seconds=production['recovery_seconds'],
        pilots_seconds=[json.loads((OUT/f'pilot-{n}.json').read_text())['seconds'] for n in [64,128]],readout_seconds=json.loads((OUT/'result.json').read_text())['seconds'],
        original_dispatch_cap_seconds=420,recovery_dispatch_cap_seconds=400,unsaved_fine_runtime='Unknown; retained as lost work under the original dispatch cap, not subtracted from its cost.',
        CPU_threads=1,peak_memory_measured=False,all_constructor_native_calls_counted=False,
        note='Action timings omit import, WSL startup, hashing, persistence verification and documentation. The original64 history and both2step pilot prefixes were reused. No extra fluid resolution, frequency, angle or time path was run.'))
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        trace='For D=rho*W, E=rho*h*W^2-p and S=rho*h*W^2*v, E-S*v-3p=rho*(cx*c^2+u)-3p.',
        packet='Two positive interpolation weights sum to one and reproduce the shifted packet energy; the same linear weights preserve radial momentum for a fixed outgoing angle.',
        boundary='Summing minus adjacent face-flux differences cancels the internal baryon faces. Opposite total-energy face fluxes cancel at the shared interface.',
        scope='Algebraic identities only; they are not a continuous EOS, full mechanical-interface or dynamical metric certificate.'))


def verify():
    for p,row in json.loads((OUT/'documentation-prefixes.json').read_text()).items():assert hashlib.sha256((ROOT/p).read_bytes()[:row['bytes']]).hexdigest()==row['sha256']
    aliases={
        'plan.json':{'verification/def_native_interior_feedback.py':'bank-producer.py','verification/def_native_two_way_atmosphere.py':'scattering-owner-before.py'},
        'production-budget.json':{'verification/def_native_interior_feedback.py':'production-producer.py'},
        'readout-plan.json':{'verification/verify_native_interior_feedback.py':'readout-producer.py'}}
    for name in ['plan.json','production-budget.json','recovery-plan.json','readout-plan.json']:
        for path,value in json.loads((OUT/name).read_text())['bindings'].items():
            relative=str(Path(path).relative_to(ROOT));target=OUT/aliases.get(name,{}).get(relative,'') if relative in aliases.get(name,{}) else Path(path)
            assert sha(target)==value,(name,path,target)
    for name in ['bank','scattering-audit','production','result','audit']:assert json.loads((OUT/(name+'.json')).read_text())['passed'],name
    for prefix in ['first-','second-','third-']:assert not json.loads((OUT/(prefix+'pilot-64.json')).read_text())['passed']
    # Check the existing atmospheric accumulated species ledger separately.
    model=task.Coupled();f=model.flow;m=model.m;rows=[]
    for steps in [64,128]:
        z=np.load(OUT/f'coupled-{steps}.npz');den=np.sum(f.initial[3]*m.vol)
        error=float(abs(np.sum((z['U'][3]-f.initial[3])*m.vol)+z['discard'][3]-z['ledger'][2]-z['ledger'][4])/den)
        assert error<1e-9
        assert int(z['completed_steps'])==steps
        rows.append(dict(steps=steps,atmosphere_accumulated_neutral_species_relative=error))
    final=np.load(OUT/'coupled-128.npz');checkpoint=np.load(OUT/'coupled-128-checkpoint.npz')
    for key in final.files:assert np.array_equal(final[key],checkpoint[key]),key
    for name in ['def_native_interior_feedback.py','verify_native_interior_feedback.py','def_native_two_way_atmosphere.py']:compile((ROOT/'verification'/name).read_text(),name,'exec')
    write(OUT/'closure-audit.json',dict(classification='Counterexample candidate',passed=True,declared_subset_gates_passed=True,
        old_document_prefixes_unchanged=True,original_inputs_or_exact_producer_archives_bound=True,original_failed_pilots_preserved=True,
        final_checkpoint_bitwise_equal=True,checkpoint_sha256=sha(OUT/'coupled-128-checkpoint.npz'),species=rows,
        full_deep_plus_atmosphere_neutral_trajectory_ledger=False,all_original_physical_gates_closed=False,
        final_charge_solved=False,full_goal_complete=False))
    print('CLOSURE PASSED; complete species ledger and full goal remain false')


def manifest():
    paths=[ROOT/'verification'/name for name in ['def_native_interior_feedback.py','verify_native_interior_feedback.py','def_native_two_way_atmosphere.py']]
    paths += [ROOT/'notes/REQUEST116_NATIVE_INTERIOR_FEEDBACK_KO.md']
    paths += [ROOT/p for p in json.loads((OUT/'documentation-prefixes.json').read_text())]
    # The final compressed trajectory already contains every checkpoint field.
    # Do not duplicate the192MB operational checkpoint or transient process ID.
    excluded={'coupled-128-checkpoint.npz','recovery-wsl.pid'}
    paths += [p for p in OUT.rglob('*') if p.is_file() and p.name not in excluded and '__pycache__' not in str(p) and p.suffix!='.tmp']
    files={str(p.relative_to(ROOT)):sha(p) for p in sorted(paths)};target=ROOT/'outputs/direct-eos-gr33/native-interior-feedback-manifest.json'
    write(target,dict(classification='Counterexample candidate',checkpoint='39e0bae7e',actual_two_way_paths_and_charge_time_comparison_passed=True,
        all_original_physical_gates_closed=False,full_goal_complete=False,final_charge_solved=False,files=files))
    p=ROOT/'paper/revision-manifest.json';d=json.loads(p.read_text());d['sha256'].update(files);d['sha256'][str(target.relative_to(ROOT))]=sha(target)
    d['native_interior_feedback']=dict(classification='Counterexample candidate',report='notes/REQUEST116_NATIVE_INTERIOR_FEEDBACK_KO.md',manifest=str(target.relative_to(ROOT)),
        actual_two_way_feedback=True,previous_near_total_direct_cancellation_not_retained=True,old_kernel_histories_recertified=False,full_goal_complete=False,final_charge_solved=False)
    write(p,d)
    for name,h in d['sha256'].items():assert sha(ROOT/name)==h,name
    print(json.dumps(dict(phase_files=len(files),master_hashes=len(d['sha256']),all_current_hashes_verified=True)))


if __name__=='__main__':globals()[sys.argv[1]]()
