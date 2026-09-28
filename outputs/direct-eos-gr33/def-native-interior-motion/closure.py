"""Freeze the native interior-motion evidence and update repository pointers."""
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
        'model-definition':'분류: Counterexample candidate. 실제 열·광자 저장 구동에 native 고정 조성 단열 미분을 연결해 기존16셀 내부의 면 바리온·운동량과 압축 압력을 진화시켰다. 원 대기의 실제 질량 이동과 내부 질량 감소를 일치시켰다. 고정 계량의 선형 기계 응답이며 광자·열로의 완전한 양방향 되먹임은 아직 없다.',
        'observable-targets':'분류: Counterexample candidate. 내부 질량 재배치를 지연 원천에 넣자 기존 직접 성분+2.2368e-27에−2.2569e-27이 더해져 약99percent 상쇄되고 직접 합은−2.0151e-29다. 기존 출사 광자 질량 정규화까지 더한 합은+4.1328e-27이다. 작은 잔여의 부호·크기는 내부 반경 오차가 닫히기 전 물리적 검출로 해석하지 않는다. 이전의 양수 성분 해석은 이 새 조건부 결과와 함께 읽는다.',
        'adiabatic-limit':'분류: Proven. 고정 조성 로그 미분의 단열 압축 계수는 K=P_rho+P_T*(P/rho-u_rho)/u_T다. 분류: Counterexample candidate. 이를 실제 native 상태에서 구해 물질 이동에 따른 배경 재고·엔트로피 이류와 압축 압력에 적용했다. 국소 단열 기계 응답을 포함했다고 전체 열평형이나 정적 EFT 흡수 경계를 뒤집지는 않는다.',
        'nonadiabatic-regime':'분류: Counterexample candidate. 내부 최대 변위0.9005mm와 상대 밀도 변화2.372e-8에도 정지질량 재배치가 직접 전하와 같은 크기였다. 작은 유체 변화만으로 전하 원천을 생략할 수 없다. 새 기계64/128 시간 대조와 독립 강제 진동자·충돌 첫 모멘트 대조는 통과했으나 원16셀 공간 오차·전체 복사 되먹임은 미폐쇄다.',
        'failure-ledger-dynamic-chi':'분류: Counterexample candidate. 첫 내부 포트의 누적 질량 부호를 잘못 적용해 내부 단독 수지는 통과했으나 내부+실제 대기 수지가 실패했다. 원 소스·명목 판정·실패 감사를 보존하고 공통 포트 부호를 고쳤다. 수정 뒤 실제 포트 대비 합산 질량 오차는1.586e-9다. 단계114의 고정 내부 직접 성분은 실제 내부 이동을 포함하면 약99percent 상쇄되므로 최종 물리 신호로 사용할 수 없다. 완전한 양방향 복사/GR·공간 오차가 남는다.',
        'dynamic-charge-completion':'분류: Counterexample candidate. actual_internal_momentum_and_baryons_evolved, native_adiabatic_pressure_feedback, joint_baryon_audit_passed는true다. full_two_way_interior_radiation, spatial_continuum_certified, full_GR_scalar_feedback, final_charge_solved, full_goal_complete는false다. 직접 성분의 큰 상쇄를 실제 계산한loophole progress다. 다음은 이 물질 압축·유속을 실제 광자·열 교환에 되돌려 상쇄 이후의 전하 합을 다시 판정하는 일이다.'}
    prefix={}
    for name,body in texts.items():
        p=ROOT/'docs'/(name+'.md');data=p.read_bytes();assert '## 단계115' not in data.decode()
        prefix[str(p.relative_to(ROOT))]=dict(bytes=len(data),sha256=hashlib.sha256(data).hexdigest())
        p.write_bytes(data+('\n\n## 단계115 — 내부 운동과 직접 전하 상쇄\n\n'+body+'\n\n상세: [단계115 보고](../notes/REQUEST115_NATIVE_INTERIOR_MOTION_KO.md).\n').encode())
    write(OUT/'documentation-prefixes.json',prefix)
    write(OUT/'operational-failures.json',dict(classification='Counterexample candidate',preserved=True,
        preparation='Referenced exterior.npz instead of existing exterior-energy.npz; prepare failed before writing plan or starting science. Corrected filename.',
        forcing='The896 forcing NPZ was written, then abs(list) failed in report serialization. Reused that NPZ; replaced with numpy.abs for the448 report. No missing saved force or repeated native bank.',
        physics='First outer baryon-port sign failed independent joint-domain audit; see first-joint-baryon-audit.json and port-repair.json. First raw result is not accepted despite its old component-only passed flag.'))
    native=json.loads((OUT/'native.json').read_text());motion=json.loads((OUT/'motion.json').read_text());first=json.loads((OUT/'first-motion.json').read_text())
    write(OUT/'resources.json',dict(classification='Counterexample candidate',native_bank_calls=native['native_calls'],native_bank_seconds=native['seconds'],
        constructor_native_calls_additional_upper_bound=32,first_motion_path_seconds=sum(x['seconds'] for x in first['paths']),corrected_motion_path_seconds=sum(x['seconds'] for x in motion['paths']),
        first_readout_seconds=json.loads((OUT/'first-result.json').read_text())['seconds'],corrected_readout_seconds=json.loads((OUT/'result.json').read_text())['seconds'],
        independent_audit_seconds=json.loads((OUT/'audit.json').read_text())['seconds'],new_heat_photon_or_atmosphere_steps=0,CPU_threads=1,peak_memory_measured=False,
        wall_clock_note='Action measurements exclude Python imports, WSL startup, hashing and documentation. Native bank cap60s, each force action30s, combined mechanical paths below45s, combined readouts below45s. Failed report forcing action elapsed was not serialized; tool wall time6.138s includes startup.'))


def verify():
    for p,row in json.loads((OUT/'documentation-prefixes.json').read_text()).items():assert hashlib.sha256((ROOT/p).read_bytes()[:row['bytes']]).hexdigest()==row['sha256']
    bindings=json.loads((OUT/'plan.json').read_text())['bindings'];producer=str(ROOT/'verification/def_native_interior_motion.py')
    assert sha(OUT/'first-producer.py')==bindings[producer]
    for p,h in bindings.items():
        if p!=producer:assert sha(p)==h,p
    for name in ['forcing','native','motion','result','audit']:assert json.loads((OUT/(name+'.json')).read_text())['passed'],name
    assert not json.loads((OUT/'first-joint-baryon-audit.json').read_text())['passed']
    for cells in [448,896]:
        old=np.load(OUT/f'first-forcing-{cells}.npz');new=np.load(OUT/f'forcing-{cells}.npz')
        for key in old.files:
            if key!='outer_mass':assert np.array_equal(old[key],new[key]),key
        source=np.load(ROOT/f'outputs/direct-eos-gr33/def-native-coupled-charge/source-{cells}-128.npz')
        assert np.array_equal(new['outer_mass'],np.asarray(source['baryon_g'].sum(1),float))
    for name in ['def_native_interior_motion.py','verify_native_interior_motion.py']:compile((ROOT/'verification'/name).read_text(),name,'exec')
    write(OUT/'closure-audit.json',dict(classification='Counterexample candidate',passed=True,actual_original_inputs_unchanged=True,
        old_document_prefixes_unchanged=True,first_port_failure_preserved=True,only_saved_forcing_port_sign_changed=True,final_charge_solved=False,full_goal_complete=False))
    print('CLOSURE PASSED; full_goal_complete=false')


def manifest():
    paths=[ROOT/'verification/def_native_interior_motion.py',ROOT/'verification/verify_native_interior_motion.py',ROOT/'notes/REQUEST115_NATIVE_INTERIOR_MOTION_KO.md']
    paths += [ROOT/p for p in json.loads((OUT/'documentation-prefixes.json').read_text())]
    paths += [p for p in OUT.rglob('*') if p.is_file() and '__pycache__' not in str(p)]
    files={str(p.relative_to(ROOT)):sha(p) for p in sorted(paths)};target=ROOT/'outputs/direct-eos-gr33/native-interior-motion-manifest.json'
    write(target,dict(classification='Counterexample candidate',checkpoint='d7708d9a2',native_forced_interior_and_joint_audit_passed=True,final_charge_solved=False,full_goal_complete=False,files=files))
    p=ROOT/'paper/revision-manifest.json';d=json.loads(p.read_text());d['sha256'].update(files);d['sha256'][str(target.relative_to(ROOT))]=sha(target)
    d['native_interior_motion']=dict(classification='Counterexample candidate',report='notes/REQUEST115_NATIVE_INTERIOR_MOTION_KO.md',manifest=str(target.relative_to(ROOT)),direct_component_almost_cancelled=True,final_charge_solved=False,full_goal_complete=False)
    write(p,d)
    for name,h in d['sha256'].items():assert sha(ROOT/name)==h,name
    print(json.dumps(dict(phase_files=len(files),master_hashes=len(d['sha256']),all_current_hashes_verified=True)))


if __name__=='__main__':globals()[sys.argv[1]]()
