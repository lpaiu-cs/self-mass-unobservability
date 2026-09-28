# 통합 원고에 백색왜성 전하 절 반영

분류: Conjectural. **목표는 단계276–280의 영문 절 초안(`docs/white-dwarf-free-fall-charge-section.md`)을 통합 원고 `paper/manuscript.md`의 §4.6으로 넣고, `main.tex`를 같은 도구로 다시 생성해 논문 검증을 통과시키는 것이다.** 사용자는 2026-09-27에 Pandoc을 설치하고 원고 통합을 진행하라고 지시했다. 기록=2026-09-27 KST.

## 준비

분류: Counterexample candidate. winget으로 Pandoc 3.11을 사용자 범위에 설치했다. 설치 파일은 github.com/jgm/pandoc/releases/download/3.11/pandoc-3.11-windows-x86_64.msi이고, winget이 설치 파일 해시를 확인했다. 요청한 설치에 필요한 winget 원본·패키지 약관 동의를 포함한다. 바꾸지 않은 `manuscript.md`로 `build_manuscript.py`를 돌리자 커밋된 `main.tex`가 바이트 단위로 재현됐다(sha256 29bb4eaf…). 따라서 원래 빌드와 같은 변환 결과를 준다.

분류: Imported from prior work. 새 참고문헌 세 편의 서지를 원 기록에서 확인했다(2026-09-27).
- Kaplan et al. 2014, "Spectroscopy of the inner companion of the pulsar PSR J0337+1715", ApJL 783, L23, doi:10.1088/2041-8205/783/1/L23, arXiv:1402.0407. 초록은 T_eff 15,800±100 K, log g 5.82±0.05, R 0.091±0.005 R☉의 수소 대기(DA) 백색왜성을 보고한다.
- Ransom et al. 2014, "A millisecond pulsar in a stellar triple system", Nature 505, 520–524 (Crossref).
- Bertotti, Iess & Tortora 2003, "A test of general relativity using radio links with the Cassini spacecraft", Nature 425, 374–376 (Crossref).

## 방법

분류: Conjectural.
1. 초안을 §4.5 뒤, §5 앞에 §4.6으로 넣는다. 원고 관례에 맞춰 §표기는 "Section X"로 바꾼다. 선언 봉투의 H·He 혼합 조성과 관측된 DA 대기의 차이를 남는 가정에 더한다.
2. 초록에 한 문장, 논의(§6)에 한 문단을 더한다. 데이터·코드 가용성 절에 새 산출물 위치를 적고, 날짜를 2026-09-27로 바꾼다.
3. `references.bib`에 세 항목을 더한다.
4. Pandoc 3.11로 `main.tex`를 다시 생성하고 `verification/verify_unified_paper.py`를 돌린다. 원고의 표 대조가 바뀌지 않았는지도 여기서 확인한다.
5. `paper/revision-manifest.json`의 원고·`main.tex`·참고문헌·README 해시를 갱신한다. PDF와 제출 zip은 TeX가 없어 다시 만들지 않는다. 이 둘이 §4.6 이전 판이라는 사실을 README와 manifest에 적는다.

분류: Conjectural. 판정 규칙: 논문 검증이 통과하고 `main.tex` 차이가 추가한 절·문장과 날짜에만 있으면 통합을 수락한다. 검증이 실패하면 원인을 기록하고 원고 변경을 되돌린다.
