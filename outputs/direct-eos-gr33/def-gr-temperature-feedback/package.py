"""Preserve the temperature-feedback result without changing Phase80's verdict."""
from pathlib import Path
import hashlib
import json
import subprocess

OUT=Path('outputs/direct-eos-gr33/def-gr-temperature-feedback')
read=lambda p:json.loads(p.read_text(encoding='utf-8'))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
write=lambda p,d:p.write_bytes((json.dumps(d,ensure_ascii=False,indent=2)+'\n').encode('utf-8'))
manifest=OUT.parent/'gr-temperature-feedback-manifest.json';assert not manifest.exists()
r=read(OUT/'result.json');audit=read(OUT/'replay.json');check=read(OUT/'temperature-check.json')
assert r['numerical_target_passed'] and r['temperature_only_feedback_below_target'] and audit['saved_replay_passed']
assert audit['accounted_compute_seconds']<720
names=['전체 속도','scalar','기존 절단면 속도','새 지지 경계 속도']
table='\n'.join(f'| {name} | {row["correction_relative"]:.8g} | {row["time_absolute_to_baseline"]:.8g} | {row["contour_absolute_to_baseline"]:.8g} |' for name,row in zip(names,r['comparisons'].values()))
report=Path('notes/REQUEST81_GR_TEMPERATURE_FEEDBACK_KO.md')
text=f'''# 단계81 — 온도 되먹임을 실제 GR 열수송 경로에 연결

분류: Counterexample candidate. **GR에서 복원한 온도를 4,012개 전도 면의 구동에 되돌려 넣고, 바뀐 열유속이 다시 GR에 작용하는 결합 응답을 계산했다.** 기존과 같은0.23080495568542375초·65시각에서 네 응답 중 가장 큰 보정은 원 응답 척도의{audit['largest_response_correction_relative']:.8g}다. 이번 온도 항에 대해 사전 판정 문턱1e−4보다 작다는 결론을 수락했다. 열 입력이 고정돼 온도 변화가 전도에 전혀 돌아가지 않던 부분은 이 제한된 모형에서 해소됐다.

분류: Imported from prior work. 단계80에서 수락한4차 GR 공간, 물질 EOS 접선, 열 pole, 보존 debit과 운동량 lift를 재사용했다. 광자·유한 이온 점유 표는 한 외곽 상태의 결과이므로 이를 전 반경 입력으로 확대하지 않았다. 새 EOS 호출은0이다. 단계80의 원 네 성분 수렴 판정과 그 이전 실패 기록은 변경하지 않는다.

## 연결한 식과 남긴 범위

분류: Counterexample candidate. 같은 물리적 압력·바리온·질량 제약으로 온도를 복원하고, 고정 Eulerian 반경에서의 구동으로 변환했다. 아래 A0,N0,T0와 수송 계수·완화율·면 위치는 배경값이다.

```text
Delta lnT = (d lnT / d ln rho)_s Delta ln rho - loss/(rho*geo*cvT)
delta lnT_E = Delta lnT - xi*d lnT0/dr
delta theta = A0*N0*T0*delta lnT_E

G(z) = K + z^2 M
u = G(z)^(-1) load E,  q = u - H E
E = E0 + R(z) [Tq*u + (TE-Tq*H)*E]
```

분류: Proven. 물리적 q=u−HE를 온도식에 대입하면 피드백 맵에도 −Tq*H가 필요하다. 같은 내부 면의 에너지를 한쪽에서 빼고 이웃에 넣으므로 총 열 debit은 상쇄된다. 두 셀의 양의 전도 행렬에 대해 합계 보존과 온도차 감쇠를 기호 검사했다. 이 항등식이 전체 GR 연산자의 균일 수축성을 증명하는 것은 아니다.

분류: Counterexample candidate. 모든 native 물질 표본에서 canonical 압력으로 별도 복원한 온도와 새 맵의 상대 차이는{check['canonical_pressure_temperature_relative']:.8g}다. 저장 GR 끝점의 최대 |delta lnT_E|는{check['maximum_Eulerian_delta_lnT']:.8g}였다. 압력·밀도·계량의 GR 응답은 유지하되, **수송 법칙의 계량 인자·전도 계수·pole 자체의 변화는 연결하지 않았다.** 외층·광자·전체 비선형 되먹임 완료로 읽지 않는다.

## 실제 전파 판정

분류: Counterexample candidate. 같은 전체 행렬을 사용하여 각 전달점의 결합식을 반복으로 풀었다. 가장 많이 사용한 보정은{r['maximum_corrections']}회이며, 표본에서 두 번째/첫 번째 보정의 최대 비는{r['sampled_second_to_first']:.8g}다. 결합식 잔차는{r['feedback_residual']:.8g}, GR 원식 잔차는{r['linear_residual']:.8g}였다. 이는4096개 역변환 contour의 표본 검사이지 모든 복소 주파수에 대한 엄밀한 수축 상계가 아니다.

분류: Counterexample candidate. 아래 값은 **보정 자체의 가중 노름을 원 응답의 최대 가중 노름으로 나눈 비율**이다. 거의 같은 두 전체 해를 빼지 않고 보정항을 직접 전파했다. 시간은1024–2048계수 차이, contour는2048–4096점 차이다. 작은 보정의 정확한 유효숫자보다, 원 응답의1e−4를 바꿀 수 있는지를 판정하도록 두 오차에 원 응답 대비1e−6 기준을 미리 정했다.

| 성분 | 온도 되먹임 보정 | 시간 대조 | contour 대조 |
|---|---:|---:|---:|
{table}

분류: Counterexample candidate. 네 성분이 같은 경로의 절대 판정 기준을 통과했다. 원 끝점의 독립 재생 차이는{r['baseline_endpoint_replay']:.8g}, 내부 열 수지 결함은{r['heat_balance']:.8g}다. 이번에는4차 보정 한 경로만 계산했으며, 그 보정에 대한 새 공간 수렴이나 물리 오차 인증을 주장하지 않는다. 기존 선형 GR의 수락 기준을 완화한 것이 아니라 별도 온도 항의 영향 판정이다.

## 비용과 다음 결정

분류: Counterexample candidate. 최초16전달점의 보수적 예상639.05초는 경로당600초 한도를 넘어 전체 실행을 보류했다. 원 풀이와 계획을 보존한 뒤 작은 보정항에만 잔차 개선 횟수를 줄였다. 같은16점의 native 보정 차이는1.6669128e−14, 기본 GR 응답은 비트 일치했고 원 잔차 문턱도 유지했다. 재측정 예상518.39초 후에 전체 계산을 실행했다. 실제 전체 경로는{r['seconds']:.2f}초, 두 pilot·독립 온도 검사·저장 재생을 합한 계산은{audit['accounted_compute_seconds']:.2f}초로 원 총720초 안이다. 단일 CPU/BLAS 스레드, 관측 최대 메모리는{r['memory_GB']:.3f}GB였다.

분류: Conjectural. 이 구간의 **온도만의 전도 되먹임**을 위해 기간이나 해상도를 확대할 근거는 얻지 못했다. 다음 가치가 큰 연결은 누락된 외층 광자·물질 에너지와 운동량 교환을 실제 반경 GR식에 넣는 것이다. 이를 위해 한 상태의 표가 실제로 적용되는 반경·상태 범위와 경계 조건을 먼저 고정해야 한다. 수송 계수·계량 되먹임, 완전한 물리 EOS 오차, 동적 전하·구동·정적 비교·관측 추론도 미완료다. 이번 성과는 **loophole progress: 온도 되먹임의 실제 결합 및 지정 구간 영향 판정**이다.
'''
report.write_bytes(text.replace('\n','\r\n').encode('utf-8'))
paragraphs={
 'model-definition.md':'분류: Counterexample candidate. 같은 GR 압력·바리온·질량 제약으로 복원한 Eulerian 온도를4012면 전도 pole에 되돌려 넣고 결합식을 실제로 풀었다. q=u−HE이므로 온도 맵도 −Tq*H를 포함한다. 수송 법칙의 계량 인자·전도 계수·pole·면 위치는 고정한 부분 폐쇄다.',
 'observable-targets.md':f'분류: Counterexample candidate. 같은65시각의 네 native 응답에 대한 온도 되먹임 보정은 최대{audit["largest_response_correction_relative"]:.6g}이며 원 응답 대비1e−4 판정 문턱 아래다. 이것은 GR 상태 읽기의 제한된 영향 판정이고 새로운 관측량·위상 지연의 비흡수성·전하 검출을 뜻하지 않는다.',
 'adiabatic-limit.md':'분류: Counterexample candidate. 원 유한 열 pole을 유지한 온도 되먹임 보정이 지정 구간에서 작았다. 이를 전체 수송의 단열 극한, 미지 외층 영유속, 정적 EFT 붕괴 또는 모든 완화 효과의 부재로 승격하지 않는다.',
 'nonadiabatic-regime.md':'분류: Counterexample candidate. 해상도·기간을 늘리지 않고 전체 행렬의 온도 되먹임을 같은 역변환으로 전파했다. 보정의 시간·contour 차이는 원 응답 대비 사전1e−6 기준을 통과했다. 새 공간 인증이나 모든 복소 주파수의 엄밀 수축 상계는 아니다.',
 'failure-ledger-dynamic-chi.md':'분류: Counterexample candidate. 온도만의 고정 수송 모형을 전체 비선형 열·광자 폐쇄로 읽는 단계는 허용되지 않는다. 빠진 항은 수송 계량 인자·전도 계수·pole·기하 변화와 외층 광자 교환이다. 최초 비용 예상은 경로당 한도 초과로 보류됐으며 원 계획·pilot을 보존하고 수치 동등성 검사 후 원 예산 안에서 실행했다.',
 'dynamic-charge-completion.md':'분류: Counterexample candidate. GR→온도→전도→GR의 실제 결합을 지정 구간에서 완결하고 원 응답 대비 작은 영향을 확인했다. 단계80의 수치 완료를 유지하되 full_temperature_metric_coefficient_feedback, full_GR_photon_feedback_evolved, full_nonlinear_evolution, full_dynamic_charge_solved는 모두 false다.'}
prefixes={}
for name,paragraph in paragraphs.items():
    p=Path('docs')/name;before=p.read_bytes();prefixes[p.as_posix()]=dict(bytes=len(before),sha256=sha(p))
    addition='\r\n\r\n## Phase81 — 온도 되먹임의 실제 결합\r\n\r\n'+paragraph+'\r\n\r\n분류: Imported from prior work. 식·원 계획·비용·독립 검사·범위는 [단계81 보고서](../notes/REQUEST81_GR_TEMPERATURE_FEEDBACK_KO.md)에 둔다.\r\n'
    p.write_bytes(before+addition.encode('utf-8'));assert p.read_bytes()[:len(before)]==before
write(OUT/'documentation-prefixes.json',prefixes)
files=[Path('verification/def_gr_temperature_feedback.py'),Path('verification/def_gr_temperature_feedback_fast.py'),report,
       *[Path('docs')/name for name in paragraphs],*[p for p in OUT.iterdir() if p.is_file()]]
write(manifest,dict(classification='Counterexample candidate',checkpoint='ad93f3b77',
    progress_class='loophole progress; actual temperature-only conductive GR feedback and finite-interval impact decision',
    numerical_target_passed=True,full_nonlinear_evolution=False,
    files={p.as_posix():sha(p) for p in sorted(files)}))
paper=Path('paper/revision-manifest.json');prior=read(paper);assert len(prior)==189 and 'request81_temperature_feedback' not in prior
entry=dict(classification='Counterexample candidate',actual_temperature_feedback_evolved=True,numerical_target_passed=True,
    temperature_only_feedback_below_target=True,full_temperature_metric_coefficient_feedback=False,
    full_GR_photon_feedback_evolved=False,full_nonlinear_evolution=False,full_dynamic_charge_solved=False,
    new_EOS_calls=0,report=report.as_posix(),report_sha256=sha(report),
    evidence_manifest=manifest.as_posix(),evidence_manifest_sha256=sha(manifest))
write(paper,dict(prior,request81_temperature_feedback=entry))
now=read(paper);assert all(now[k]==v for k,v in prior.items())
stage=[p.as_posix() for p in files]+[manifest.as_posix(),paper.as_posix()]
subprocess.run(['git','-c','core.autocrlf=false','add','-f','--',*stage],check=True)
print(json.dumps(dict(manifest_files=len(files),staged_files=len(stage),paper_entries=len(now),report=report.as_posix()),indent=2))
