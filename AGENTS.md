# AGENTS.md

## Mission

This repository is dynamic-loophole-first.

Status: Imported from prior work. The earlier theorem work left A4 as the active assumption to attack: no orbital-timescale internal state variable in the free-fall sector.

Status: Conjectural. This repository asks whether the smallest such state, `chi_A`, creates a real observable beyond a static finite-dimensional sensitivity-manifold EFT.

## Primary Scope

Work on the free-fall-style response model first.
Do not start with the clock sector.
Do not reopen LLR / MLRS / PEP / Nutimo / runtime / build-environment work.
Do not make static primitive-family audits the mainline unless they are strictly needed to define the drive basis `Y`.

## Core Target

Status: Counterexample candidate. The first target is

```text
tau_chi * d chi_A / dt + chi_A = alpha * F(Y)
m_A(Y, chi_A) = m_A^(0) * [1 + c_Y F(Y) + c_chi chi_A]
```

Status: Conjectural. The target observable must be one of:

- phase-lagged quadrature relative to a static drive,
- sidebands or mixed-frequency response,
- frequency-dependent transfer not absorbable into static coefficients,
- a sharply stated no-go boundary showing collapse of the minimal model.

## Allowed Outputs

Codex may produce:

- theorem or no-go statements,
- assumption ledgers,
- short analytic derivations,
- symbolic response checks,
- explicit counterexample candidates,
- failure ledgers with exact collapse conditions.

Codex must not:

- do more empirical LLR / MLRS / PEP work,
- do more pulsar / Nutimo runtime work,
- invent empirical claims without derivation or source,
- silently strengthen assumptions,
- count a static parameter redefinition as novelty.

## Claim Labels

Every substantive scientific claim must be tagged as exactly one of:

- Proven
- Imported from prior work
- Conjectural
- Counterexample candidate

## Maintained Files

Maintain these files during dynamic-chi work:

- `docs/model-definition.md`
- `docs/observable-targets.md`
- `docs/adiabatic-limit.md`
- `docs/nonadiabatic-regime.md`
- `docs/failure-ledger-dynamic-chi.md`

If a model fails to escape collapse, record the exact failing step and the minimal missing assumption in `docs/failure-ledger-dynamic-chi.md` (the theorem track owns `docs/failure-ledger.md`).

## MVP Done Rule

A task is done only if:

- the relevant markdown notes are updated,
- symbolic checks run without error,
- and the result is classified as theorem progress or loophole progress.

For this MVP, done means:

1. the one-state `chi_A` model is defined,
2. the monochromatic response is solved,
3. the adiabatic collapse boundary is written,
4. the non-adiabatic observable classification is written,
5. either a genuine observable candidate is isolated or the exact no-go boundary is stated.

## Git Discipline

Before each major task:

- create a checkpoint commit.

After each major task:

- create a checkpoint commit with a one-line scientific summary.

## Long Computation Budget

2026-09-16 사용자 합의: 오래 걸리는 계산을 함부로 남발하지 않는다. 연구 레버를 모두 진행하라는 포괄적 지시를 무제한 계산 예산으로 해석하지 않는다.

- 장기 계산을 시작하거나 크게 확대하기 전에 해결할 주장, 결과에 따라 달라질 결정, 필요한 정확도와 수락 기준을 명시한다. 수치 검증의 성공이 어떤 물리적 주장까지 뒷받침하는지 구분한다.
- 저장 결과 재사용, 해석적 판단, 작은 대표 구간 등 더 저렴한 대안부터 확인한다. 이미 수락한 동일 방정식의 이력을 이유 없이 재계산하지 않는다.
- 대표 구간의 실측 속도로 예상 벽시간, CPU/GPU·메모리 사용, 실행·반복 예산을 정한다. 서로 다른 경로나 후반 구간의 미측정 속도는 가정으로 표시하고, 점 추정만 제시하지 않는다.
- 종료·실패·중단 기준과 추가 계산의 범위를 실행 전에 정한다. 목표상 이득이 불분명하거나 예산을 넘으면 새 장기 계산을 시작하기 전에 계획을 재평가한다.
- 기준 미달을 이유로 해상도·기간·경로 수를 자동 확대하지 않는다. 과거 실패를 보존하고 수락 기준을 완화하지 않는다. 검증 기록의 개수보다 실제 병목 해결을 우선한다.
- 현재 승인된 GR 계산은 원 70/140/280 세 경로와 최종 판정까지 완료한다. 그 이후 더 촘촘한 경로를 자동 실행하지 않는다. 실행 중인 소스·계획은 바꾸지 않는다.
- 기존 승인 범위의 자율 진행과 반복 허락 생략은 유지한다. 실행 세부와 측정치는 저장소 문서에, 재사용할 교훈은 기존 OSK 상세 노드에 기록한다.

2026-09-24 후속 사용자 지시: 계산 예산을 지나치게 촘촘히 잡아 실험을 반복 중단하지 않는다. 저장 상태에서 실험을 속행하는 편이 전체 비용을 줄이면 벽시간과 반복 예산에 충분한 여유를 둔다. 정확도·물리 수락 기준은 유지하고, 기존 실패와 예산 변경 이유를 기록한다. 이 지시는 무제한 격자·기간·경로 확대가 아니다.

## Research Value Criterion

2026-09-24 사용자 지시: 연구의 가치는 **현재 지배적인 오차를 해결한 동일 결합 해에서 최종 전하의 결론이 유지되는가**로 평가한다.

- 원천·보존·예비 적분·개별 수렴 통과는 최종 판정을 위한 중간 검증이다. 이를 최종 과학적 성과나 전체 완료로 대체하지 않는다.
- 지배 오차를 수정한 실제 광자·물질·GR 결합 해와 같은 해의 에너지·경계 이력을 최종 전하 판독에 연결한다. 다른 해의 진단값을 사후 가산해 결론을 만들지 않는다.
- 보고에는 최종 전하의 기존 결론이 유지되는지, 변경되는지, 아직 판정 불가인지 먼저 밝히고 그 근거와 남은 오차를 구분한다. 원 실패와 수락 기준은 유지한다.
