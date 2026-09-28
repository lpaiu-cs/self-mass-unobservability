# 통합 연구 후속 검증 완료 기록

이 문서는 1–3번 레버의 완료 당시 기록이다. 이후 4–5번의 추가 결과와 현재 원고는 [최신 완료 기록](levers-4-5-completion-2026-09-09.md)에 정리했다.

사용자 승인 순서인 OSK 정리 → 1 nuisance → 2 coverage → 3 물리적 matching으로 진행했다. 입력 통합 원고는 4897038이다. 이 기록은 이전 `unified-revision-2026-09-09.md`의 남은 과제 목록을 아래 범위에서 갱신한다. 과거 REQUEST10 원시 결과와 검출 판정은 보존한다.

| 단계 | 완료한 결과 | 판정과 범위 |
| --- | --- | --- |
| 1. nuisance | 90차원 전체 공간, 19차원 누락 공간과 약한 모드의 혼합을 재구성했다. 기존 표를 재현하고 cutoff/편향 민감도를 계산했다. | Imported from prior work: 누락 방향을 물리적 0으로 둘 근거가 없다. 전체 공간을 기본으로 채택했다. |
| 2. coverage | 조건별 8192회 생성·추정 검증, 전체 grid의 별도 512회 검사, 알려진 공분산 양성대조군, 공분산을 REML로 추정하는 후속 검증을 실행했다. | Imported from prior work: 상관잡음 무시 시 최소 58.25%, 지정한 공분산 추정 시 최소 94.59%의 점별 포함률. 이 값들은 시험 사례의 기술적 최솟값이다. |
| 3. 물리적 matching | 독립 scalar charge의 작용·힘·에너지 수지, 단일극 근사 오차, 계수 대응, 동반성 되먹임 및 드라이브 위상 조건을 유도하고 14개 검사를 통과했다. | Counterexample candidate: 조건부 EFT 구현을 얻었다. Proven: 기존 보조 physical-drive의 위상 closure가 유도와 불일치하며 공통 시간 이동으로 수정할 수 없다. |

Status: Proven. 안정 분지와 동일한 고정 동반성 전하/질량비에서 tau=Gamma/kappa, beta=Ustar*a_w^2/(kappa*m_p)이다. 기존 tau→입자 질량 및 궤도 준공명 해석은 성립하지 않는다. 단일극 근사는 관성항·비선형성·동반성 되먹임·초기상태에 대한 별도 조건을 요구한다.

Status: Proven. 유도한 선도차수 공면 potential drive의 위상 closure는 저장 매개변수에서 3.11837 rad이며, 기존 보조 분석은 0이다. 그 보조 beta_phys 수치를 이 물리 모형의 상한으로 해석하는 주장은 철회했다. 기존 단위-drive beta 표는 정의된 현상론적 벤치마크로만 남긴다.

Status: Imported from prior work. K=10 점별 구간의 실패 사례와 전체 grid envelope 검사를 구별했다. 실제 시험한 K=10 envelope의 최소 포함률은 511/512였다. 따라서 원래 envelope가 검증에서 실패했다고 주장하지 않으며, 유한한 성공 사례를 보편 보장으로 확대하지도 않는다.

Status: Proven. 알려진 공분산·불편 Gaussian 추정량에서는 사용한 대칭 Gaussian-mass 구간의 빈도주의 포함률이 최소 95%임을 증명했다. 공분산 추정의 성능은 별도 Monte Carlo 결과이고, 이 증명의 자동 귀결이 아니다.

## 원고와 검증

통합 원고에 4.3절 조건부 scalar-charge 구현, 4.4절 물리 위상 gate, 5.5절 nuisance/coverage 감사, 표 2와 그림 2, 등록·재현 부록을 추가했다. 초록·서론·논의·정준 문서 다섯 개를 실제 결과에 맞췄다. 상세 설계와 결과는 `notes/REQUEST11_*`에 있다.

완성 원고는 `paper/manuscript.md`, 생성 TeX는 `paper/main.tex`, PDF는 `output/pdf/free-fall-identifiability.pdf`, 이동 가능한 소스 묶음은 `output/submission/free-fall-identifiability-source.zip`이다. `paper/package_revision.py`가 새 분석·입력·결과·문서의 SHA-256 manifest를 작성한다.

작업 분류: theorem progress 및 loophole progress. 이론 경계와 조건부 검증을 갖춘 하나의 방법론 논문으로 정리했다. 실패한 물리 해석을 성공한 경험적 결과로 바꾸지 않았다.

검증 기록: symbolic smoke 통과, exact-character 8/8 sector 일치, 물리 matching 14개 검사 통과, 통합 검증의 이론·기존 표·coverage 표·결과 counts·manifest 검사 통과. PDF 16페이지를 전부 렌더링하여 확인했고, 소스 ZIP을 별도 디렉터리에 풀어 재빌드한 PDF의 텍스트와 16페이지 렌더링이 원본과 동일했다. TeX 로그에 미해결 인용·참조, missing character, overfull/underfull box가 없다. 기존 REQUEST10 파일의 diff는 없다.

OSK는 기존 프로젝트 개요·A/B·물리 후보·허브 노드를 정정하고 scope 요약을 갱신했다. 기존 완화시간→입자질량 및 준공명 해석도 철회했다. OSK 검증기는 PASS, 해당 연구 scope의 도달 불가능 노드는 0이었다. 다른 scope의 기존 dangling-reference 경고는 이 작업의 범위 밖이다.

## 여전히 주장하지 않는 결과

Status: Conjectural. 실제 중성자별 EOS와 백색왜성의 전하를 계산한 수치적 body matching, 일 단위 완화의 실현, 수정된 물리 드라이브에 대한 상한, 모든 천체물리 잡음을 포괄하는 timing likelihood는 이번 결과로 확정되지 않았다. 필요한 물리 입력이 없는 부분을 임의로 정해 완료했다고 하지 않는다.

현재 제출 주장은 조건부 이론·방법론이다. 보편적 SEP 개선이나 특정 scalar-tensor 이론 배제를 주결론으로 삼으려면 위 추가 물리·관측 작업이 필요하다. 저널 투고·공개 push·DOI 등록은 수행하지 않았다.
