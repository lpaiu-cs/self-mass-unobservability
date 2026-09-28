"""Bind the repaired initial projection without promoting it to evolution."""
from pathlib import Path
import hashlib,json

root=Path(__file__).resolve().parents[3];out=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text())
def write(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
result=read(out/'finite-volume/balanced-result.json');audit=read(out/'finite-volume/audit.json')
assert result['passed'] and audit['mathematical_projection_passed'] and audit['constraint_input_matched']
assert not audit['production_ready'] and not result['final_charge_solved']
provenance=[]
for plan,archive,field in [('plan.json','cauchy-producer.py','bindings'),('balanced-plan.json','balanced-producer.py','bindings'),('audit-plan.json','audit-producer.py','source_sha256'),('finite-volume/first-audit-plan.json','finite-volume/first-audit-producer.py','source_sha256')]:
    p=read(out/plan);h=sha(out/archive)
    if field=='bindings':
        values=[v for k,v in p[field].items() if k.endswith('def_native_initial_constraints.py')]
        assert values and all(v==h for v in values),(plan,archive)
    else:assert p[field]==h,(plan,archive)
    provenance.append(dict(plan=plan,executed_producer_archive=archive,sha256=h))
write(out/'executed-provenance.json',dict(passed=True,rows=provenance,
    explanation='Earlier source-bound plans retain their executed hashes. Their code is preserved in these archives; the current implementation includes subsequent explicit repairs.'))
docs=[root/'docs'/n for n in ['model-definition.md','observable-targets.md','adiabatic-limit.md','nonadiabatic-regime.md','failure-ledger-dynamic-chi.md','dynamic-charge-completion.md']]
assert not (out/'documentation-prefixes.json').exists()
write(out/'documentation-prefixes.json',{str(p.relative_to(root)):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in docs})
paragraphs={
'model-definition.md':'분류: Counterexample candidate. 실제19+512셀의 저장 질량·광자 에너지/반경 압력을 보존하는 구간별 유한체적 원천으로 초기 질량·lapse·순간 스칼라 균형을 풀었다. 중심 정칙성·phi_infinity=.001과Pi_phi=0을 유지하고 실제 광자 순유속의 외재곡률을 저장했다. 가스에는 native로 고정한 국소 등엔트로피 Gamma 연장을 쓴다. 매끈한 중심값 보간은 실제 셀 재고와 달라 기각했다. 이 새 초기 상태를 기존 진화기에 아직 설치하지 않았다.',
'observable-targets.md':'분류: Counterexample candidate. 초기 배경의 정적 질량·K 보정은 동적 관측 전하가 아니다. 새 초기 입력의 셀별 재고·native 대조는 통과했지만 새 실제 결합 궤적과 전하 판정은 없다. 단계119의 조건부 양수 하한을 새 배경에 승계하지 않는다. 작은 계량 상대 보정만으로 훨씬 작은 동적 신호의 영향이 작다고 추정하지 않는다.',
'adiabatic-limit.md':'분류: Proven. 선언한 Gamma 연장 E(w)=E0*w+p0*(w^Gamma-w)/(Gamma-1), P(w)=p0*w^Gamma는dE/dlnw=E+P를 만족한다.\n\n분류: Counterexample candidate. 이 구성 관계와 순간 스칼라 균형은 열적 정상상태나 전체 native EOS의 증명이 아니며 정적 EFT 흡수 경계를 바꾸지 않는다.',
'nonadiabatic-regime.md':'분류: Counterexample candidate. 새 초기 스칼라 가속 균형은 원 저장 Hermite Cauchy장을 유지한 해와 다른 명시적 초기 조건이다. 원 Cauchy 해의 잔차를 보존했다. 광자 유속·외재곡률이 비영이므로 순간 scalar 균형을 정확한 정적 복사 항성으로 부르지 않는다. 새 상태의 실제 비단열 응답은 재진화 후 판정해야 한다.',
'failure-ledger-dynamic-chi.md':'분류: Counterexample candidate. 첫 초기 제약 해의 매끈한 중심값 보간은 실제 셀 질량과 최대15.8312percent 불일치했다. 저장 셀 내용물의 유한체적 원천으로 바꾸어2.9420e-16까지 일치시켰다. 첫 native 역산은 고정 재고 경로에 평형 열미분을 사용해 실패했고, 저장된 constrained cvT로 고쳐 원 문턱을 통과했다. 원 실패·초기 Hermite scalar 가속 잔차·20호출 진단 종료를 보존한다.',
'dynamic-charge-completion.md':'분류: Counterexample candidate. actual_finite_volume_initial_source_matched, momentarily_scalar_balanced_initial_projection, known_region_momentum_constraint, native_initial_anchor_checks는true다. smooth_anchor_initial_source_accepted, new_initial_state_installed_in_transport, new_coupled_trajectory_completed, full_native_continuum_EOS, complete_core_momentum_source, full_dynamic_GR_feedback, final_charge_solved, full_goal_complete는false다. 기존initial_full_Einstein_constraints_matched의 물리적 전 범위 판정도false로 유지한다. 다음은 새 계량·체적·EOS 기준·광자 주파수 표현·중력 힘을 함께 실제 진화에 설치하는 일이다.'}
for p in docs:
    added='\n\n## 단계120 — 실제 셀 재고와 초기 제약 연결\n\n'+paragraphs[p.name]+'\n\n상세: [단계120 보고](../notes/REQUEST120_INITIAL_CONSTRAINTS_KO.md).\n'
    p.write_bytes(p.read_bytes()+added.encode('utf-8'))
paths=docs+[root/'notes/REQUEST120_INITIAL_CONSTRAINTS_KO.md',root/'verification/def_native_initial_constraints.py',root/'verification/verify_native_initial_constraints.py']
paths+=sorted(p for p in out.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc')
files={str(p.relative_to(root)):sha(p) for p in paths};manifest=root/'outputs/direct-eos-gr33/native-initial-constraints-manifest.json'
write(manifest,dict(classification='Counterexample candidate',checkpoint='537d57c',actual_finite_volume_initial_source_matched=True,
    momentarily_scalar_balanced=True,new_initial_state_installed_in_transport=False,new_fluid_steps=0,full_dynamic_GR_feedback=False,final_charge_solved=False,full_goal_complete=False,files=files))
master=root/'paper/revision-manifest.json';d=read(master);d['sha256'].update(files);d['sha256'][str(manifest.relative_to(root))]=sha(manifest)
d['native_initial_constraints']=dict(classification='Counterexample candidate',report='notes/REQUEST120_INITIAL_CONSTRAINTS_KO.md',manifest=str(manifest.relative_to(root)),
    actual_finite_volume_initial_source_matched=True,momentarily_scalar_balanced=True,smooth_anchor_initial_source_accepted=False,
    new_initial_state_installed_in_transport=False,full_native_continuum_EOS=False,full_dynamic_GR_feedback=False,final_charge_solved=False,full_goal_complete=False)
write(master,d)
for path,h in d['sha256'].items():assert sha(root/path)==h,path
for path,p in read(out/'documentation-prefixes.json').items():assert hashlib.sha256((root/path).read_bytes()[:p['bytes']]).hexdigest()==p['sha256'],path
print(json.dumps(dict(manifest_files=len(files),master_files=len(d['sha256']),passed=True)))
