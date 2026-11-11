# Physica Scripta 투고 준비본 - 30 September 2026

상태: **저자 검토·직접 제출용 패키지. 아직 제출하지 않았음.**

제목: *Finite-frequency identifiability of relaxation responses in free-fall timing*.
대상: **Physica Scripta / Research Paper / 구독형 출판(선택적 OA 신청 안 함)**.

## 제출 파일

| 파일 | 용도 |
| --- | --- |
| `manuscript.pdf` | 심사용 본문. 방법과 증명·정확 연산 대조 예제·조건부 J0337 사례를 포함한다. |
| `supplement.pdf` | Supplementary Information. 상세 감사·정적 분류·항성 계산과 한계를 보존한다. |
| `supplementary-code.zip` | Supplementary Code. 정확 연산 예제, 저장 comparator 기록, 공개 JSON 추론 재생 및 필요한 소스 해시 대조 파일. |
| `cover-letter.pdf` | 편집자 서한. 이전 GRG 편집 단계 거절과 이번 기여·한계를 명시한다. |
| `manuscript-source.zip` | 요청 시 제출할 본문 LaTeX 소스와 bibliography. 첫 제출 PDF와 대응한다. |
| `abstract-plain.txt`, `submission-fields.txt` | 웹 양식에 옮길 초록과 메타데이터. |
| `manifest.json` | 실제 전달 파일의 SHA-256, 소스 바인딩, 준비 상태. 접수증이 아니다. |

`manuscript.md`와 `supplement.md`가 이 저널 준비본의 편집 원본이다. `main.tex`와 `supplement.tex`는 생성물이다. 이전 `paper/manuscript.md` 및 `output/submission/grg/`는 보존된 이전 판본이며 덮어쓰지 않는다. Supplementary Information의 과학 본문과 부록은 이전 GRG 판본을 보존한다. 새 정확 연산 예제는 본문과 코드에서 별도로 식별한다.

## 이번 수정의 범위

**Proven.** 비교 모형의 유한 주파수 보간 가능성, nuisance 제거 후 분리 가능성, 측정 정보량을 연결하는 절차를 본문에 명시한다. 표준 투영·보간 결과를 새로운 일반 수학으로 주장하지 않는다. 무차원 대조 예제는 동일한 응답이 그대로 남는 경우, 투영으로 사라지는 경우, full rank인데 임의로 약해지는 경우를 정확 연산으로 확인한다.

**Imported from prior work.** J0337 수치는 기존 저장 기록을 재사용한다. 미분 오차, 2일 지연의 근소한 거절, 지정 공분산 밖의 미검증 범위와 과거 잘못된 구동 위상 해석을 유지한다. 새로운 timing/GR/항성 진화 또는 Monte Carlo 계산을 하지 않는다.

**Conjectural.** 원래 최종 전하의 물리·관측 폐쇄 여부는 이번 편집으로 개선되거나 확정되지 않는다. 항성 계산과 열 지연 꼬리의 미보장은 SI에 남으며 본문 방법을 검증하는 증거로 쓰지 않는다. 이 패키지의 검수 통과는 저널 수락 가능성이나 새 empirical SEP 제약의 인증이 아니다.

## 재현과 빌드

저장소 루트에서 구성된 Windows 호스트는 명령 앞에 `rtk proxy`를 붙인다.

```text
python verification/projected_identifiability.py
python verification/replay_public_inference.py
python symbolic/state_identifiability.py
```

앞의 두 스크립트는 실제 timing engine을 실행하지 않는다. 세 번째는 기존 14개 정확 항등식을 확인하고 그 파생 JSON만 다시 쓴다. Supplementary Code ZIP을 독립 폴더에 풀어도 같은 세 명령이 동작한다. 재생의 성공은 미분이나 관측 잡음의 검증과 구분한다.

빌드는 기존 `paper/build_manuscript.py`, Pandoc, Tectonic을 사용한다. PowerShell에서 `TECTONIC` 환경변수에 설치된 실행 파일을 지정한 후 `python paper/build_physica_scripta.py`를 실행한다. 현재 호스트 실행 파일은 `C:/Users/lpaiu/AppData/Local/Programs/tectonic/tectonic.exe`이다. 이 빌더는 GRG manifest의 모든 기존 파일이 그대로인지 먼저 확인하고, 현재 폴더의 Markdown·참고문헌·그림과 cover-letter.tex만 빌드한다. 원래 연구 manifest를 덮어쓰지 않는다.

## 저자가 직접 마무리할 절차

1. 본문과 커버레터의 과학적 해석·AI 사용 공개·책임 진술을 검토한다.
2. [공식 저널 페이지](https://publishingsupport.iopscience.iop.org/journals/physica-scripta/)의 Submit an article 경로를 열고 Research Paper로 제출한다. 준비한 파일과 `submission-fields.txt`를 사용한다. ORCID 등 확인되지 않은 식별자는 임의로 추가하지 않았다.
3. 일반 구독형 출판을 선택하면 출판료는 없다. 유료 OA 선택이나 결제는 이 패키지에 포함하지 않는다. [공식 비용 안내](https://publishingsupport.iopscience.iop.org/journals/physica-scripta/about-physica-scripta/)
4. 저자 정책 선언 및 최종 제출은 본인이 수행하고 접수번호를 보존한다. IOP는 제출 저자 본인이 제출해야 한다고 명시한다. [Identity fraud and impersonation 정책](https://publishingsupport.iopscience.iop.org/ethical-policy-journals/)

IOP의 실질적 AI 사용 공개 요건에 따라 Methods와 Acknowledgements에 도구·사용 범위·검증과 한계를 기술했다. GRG 제출은 2026-09-29 거절로 끝났으며 이전 패키지와 접수·결정 증거는 별도 보존한다.
