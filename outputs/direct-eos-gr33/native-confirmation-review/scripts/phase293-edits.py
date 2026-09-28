"""Phase 293: minor fixes from the confirmation reviews (all three: minor revision then accept).

Scope of the J0337 conclusion (six evaluated lags, instantaneous plus single relaxation), marginal 2-day rejection and its
robust form for beta >= 0, archived-derivative condition, ratio range, diagonal weighting of the K=1 columns, relation of the two
p-values, unit amplitudes in the comparator phase sets, REML citations, internal wording, and SM alignment. CRLF kept; every
anchor must occur exactly once.
"""
from pathlib import Path
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')

def apply(rel, edits):
    p = root/rel; t = p.read_bytes().decode('utf-8')
    for old, new in edits:
        assert t.count(old) == 1, (rel, old[:70])
        t = t.replace(old, new)
    p.write_bytes(t.encode('utf-8')); print(rel, len(edits), 'edits')

apply('paper/manuscript.md', [
    ("had narrowed the intervals on the relaxation coefficient by factors of 2.1--27.5.",
     "had narrowed the intervals on the relaxation coefficient by factors of 2.1 to 27.5, depending on lag and interval width."),
    ("The recorded carrier coefficients exceed its null threshold, but no combination of instantaneous and relaxing response to the leading physical drive at any evaluated lag reproduces them; no relaxation is detected, and the excess is not attributed.",
     "With the archived timing derivatives, the recorded carrier coefficients exceed its null threshold, but an instantaneous plus single-relaxation "
     "response to the leading physical drive does not reproduce them at any of six evaluated lags (marginally at 2 days, robustly for the sign that the "
     "scalar-charge model requires); no relaxation is detected, and the excess is not attributed."),
    ("and a phase gate on the physical drive in J0337 (Section 4).\r\n",
     "and a phase gate on the physical drive in J0337 (Section 4). These boundaries specialize classical results, recalled below; what is new is their "
     "formulation for a nuisance-projected free-fall timing measurement, the comparator classes fixed by conservative or reciprocal dissipative physics, "
     "the drive phase gate, and their application to a stored J0337 analysis. They do not establish a new theory of gravity or a general unobservability theorem.\r\n"),
    ("**Imported from prior work.** Each boundary specializes a classical result. Theorem 1 is", "**Imported from prior work.** Theorem 1 is"),
    (" and the transient obstruction is a statement of insufficient excitation. What is new is their formulation for a nuisance-projected free-fall timing measurement, the comparator classes fixed by conservative or reciprocal dissipative physics, the drive phase gate, and their application to a stored J0337 analysis. These results do not establish a new theory of gravity or a general unobservability theorem.",
     " and the transient obstruction is a statement of insufficient excitation."),
    ("and \\(t_{\\rm asc}=t_{\\rm pericenter}-P\\varpi/(2\\pi)\\) \\cite{voisin2025planet}. **Proven.** The frozen parameters then give \\(\\varpi_p=1.69406454\\) and \\(\\varpi_b=1.67084424\\) radians (97.06 and 95.73 degrees, close to the 97.62 and 95.62 degrees of Ref. \\cite{ransom2014triple}) and",
     "and \\(t_{\\rm asc}=t_{\\rm pericenter}-P\\varpi/(2\\pi)\\) \\cite{voisin2025planet}; the published pericenters are 97.62 and 95.62 degrees \\cite{ransom2014triple}. "
     "**Proven.** The frozen parameters then give \\(\\varpi_p=1.69406454\\) and \\(\\varpi_b=1.67084424\\) radians (97.06 and 95.73 degrees) and"),
    ("**Counterexample candidate.** We therefore withdraw the physical-beta numbers of that dictionary as constraints on this realization.",
     "**Counterexample candidate.** Coupling limits that our earlier analysis derived from that dictionary are therefore withdrawn as constraints on this realization."),
    ("The registered-grid causal maximum of \\(|\\widehat\\beta|/\\sigma_F\\) is 2.2851,",
     "For the causal (lagging) template, the registered-grid maximum of \\(|\\widehat\\beta|/\\sigma_F\\) is 2.2851,"),
    ("the narrowest stored construction, fails under omitted means (Section 5.3).",
     "the narrowest stored construction, fails under omitted means (Section 5.3). Both K=1 columns use the stored diagonal weighting, which also fails "
     "under correlated extra Fourier power (Table 2); the full-space K=1 values come from the registered nuisance audit."),
    ("The last column is the full/truncated ratio at K=10.}", "The last column is the full/truncated ratio at K=10; the Full, K=1 column is from the registered nuisance audit.}"),
    ("by restricted maximum likelihood (REML), where the columns of L",
     "by restricted maximum likelihood (REML) \\cite{patterson1971reml}, a simple form of the correlated-noise likelihoods used in pulsar timing "
     "\\cite{vanhaasteren2013noise}, where the columns of L"),
    ("and the leading physical phases of Section 4)", "and the leading physical phases of Section 4, all with unit carrier amplitudes)"),
    ("so the data reject vanishing carrier coefficients. The best fit within the physical plane leaves a statistic of 12.84 at 2 days, just above the threshold, and 14.39--14.59 at the other five evaluated lags. Every physical lag section is therefore empty: no combination of instantaneous and relaxing response to the leading physical drive reproduces the carrier coefficients at the calibrated 95 percent level, and no interval on \\(\\beta\\) or EOS-matched bound follows. **Conjectural.** Every relaxing response to this drive lies in the rejected plane, so the excess is not a relaxation signal; it is not attributed.",
     "so the data reject vanishing carrier coefficients. This test differs from the unit-drive scan of Section 5.2 (p=0.26) in statistic, template and "
     "covariance: it tests all six carrier coefficients at once. The best fit within the physical plane leaves a statistic of 12.84 at 2 days and "
     "14.39--14.59 at the other five evaluated lags. Every evaluated physical lag section is therefore empty: at these six lags, no instantaneous term plus "
     "single relaxation driven by the leading physical drive reproduces the carrier coefficients at the calibrated 95 percent level, and no interval on "
     "\\(\\beta\\) or EOS-matched bound follows. The 2-day rejection is marginal: the frozen threshold is the largest of four calibration order statistics "
     "(12.27--12.82), each with a Monte Carlo standard error near 0.13, although 12.84 exceeds all four and the nominal chi-square value 12.59. "
     "**Proven.** The rejection is robust for the equal-charge realization of Section 4, which has \\(\\beta\\ge0\\). The statistic is a convex quadratic in b, "
     "the \\(\\beta=0\\) line lies in every section, and the unconstrained \\(\\widehat\\beta\\) is negative at five lags and \\(7.4\\times10^{-12}\\) at "
     "18 days; so with \\(\\beta\\ge0\\) every section minimum is at least the 18-day value, 14.59. **Conjectural.** These results use the archived timing "
     "derivatives, which Section 5.6 shows to be uncertain. At the evaluated lags the excess is not explained by a response to this drive, so it gives no "
     "evidence for relaxation; most of it lies in directions that the drive does not span, and it is not attributed."),
    ("and calibrated inference is available within specified covariance families. The recorded carrier coefficients carry an excess that the leading physical drive, with or without relaxation, does not reproduce; no relaxation is detected. The measured derivative sensitivity and the unphysical full pulse compensation prevent empirical promotion.",
     "and calibrated inference is available within specified covariance families. At six evaluated lags between 2 and 500 days, the recorded carrier "
     "coefficients carry an excess that the leading physical drive, with or without a single relaxation, does not reproduce; no relaxation is detected. "
     "The physical-drive fits have standard errors of \\(3.5\\times10^{-10}\\) to \\(1.7\\times10^{-9}\\) in \\(\\beta\\), about 2 to 10 in \\(B=\\beta/U_*\\). "
     "The measured derivative sensitivity and the unphysical full pulse compensation keep the analysis from serving as an empirical constraint."),
    ("and no mechanism that produces a relaxation at the stored 2--500-day lags with a coupling these data could detect is identified here.",
     "and no mechanism that produces a relaxation at the stored 2--500-day lags with a coupling B of that size is identified here."),
    ("the verification programs named below and these reviews.", "the verification program named below and these reviews."),
])

apply('paper/supplement.md', [
    ("**Repository:** lpaiu-cs/self-mass-unobservability\r\n", ""),
    ("unchanged except for corrections found in the final review: the periastron convention of the J0337 drive (Sections 4.4, 5.6, 5.9 and 5.10, Figures 3 and 4, and the scale list of Section 4.6) and the reproduction notes of the data-availability section.",
     "unchanged except for corrections made in this version: the periastron convention of the J0337 drive (Sections 4.4, 5.6, 5.9 and 5.10, Figures 3 "
     "and 4, and the scale list of Section 4.6), wording aligned with the main text in Sections 3.3, 3.5 and 5.3 and in the Use of AI tools section, "
     "and the reproduction notes of the data-availability section."),
    ("An earlier version of this analysis read eta as \\(e\\cos\\varpi\\); the final review found the error, and the values below and in Sections 5.6, 5.9 and 5.10 are corrected.",
     "An earlier version of this analysis read eta as \\(e\\cos\\varpi\\); the values below and in Sections 5.6, 5.9 and 5.10 are corrected in this version."),
    ("are distinct positive frequencies with known nonzero drive and deprojected readout.",
     "are distinct positive frequencies at which the drive is known and nonzero, so that the values \\(G(i\\omega_k)\\) are available."),
    ("Instantaneous mass at zero time is inseparable from \\(c_0\\).", "A relaxation-time atom at \\(\\tau=0\\) would be indistinguishable from \\(c_0\\)."),
    ("This does not address the factor 2.1--17.4 full/truncated dependence,", "This does not address the factor 2.1--17.4 full/truncated dependence over the 65 stored lags,"),
    ("all above the threshold, so every evaluated physical lag section is empty: no combination of instantaneous and relaxing response to the leading physical drive reproduces the carrier coefficients at the calibrated level.",
     "all above the threshold, so every evaluated physical lag section is empty: at these six lags, no instantaneous term plus single relaxation driven by "
     "the leading physical drive reproduces the carrier coefficients at the calibrated level. The 2-day rejection is marginal; with the sign \\(\\beta\\ge0\\) "
     "of the equal-charge realization every section minimum is at least 14.59. The region, its threshold and the omnibus statistic do not depend on the "
     "drive phases, so only the sections were recomputed."),
    ("GPT-6-Astra, Claude Opus 5.5 and Claude Fable 5.1 independently reviewed drafts of Section 4.6 in several rounds;",
     "GPT-6-Astra, Claude Opus 5.5 and Claude Fable 5.1 (Anthropic, run as a Claude Code subagent) independently reviewed drafts of Section 4.6 in several "
     "rounds and the complete manuscript before submission;"),
])

apply('output/submission/cover-letter-prd.tex', [
    ("nuisance truncation had narrowed its intervals by factors of 2.1--27.5, a calibrated six-coefficient region is validated on independent simulations, and the recorded carrier coefficients show an excess that the leading physical drive, with or without relaxation, does not reproduce. No relaxation is detected, and failed empirical-promotion gates are reported.",
     "nuisance truncation had narrowed its intervals by factors of 2.1 to 27.5, depending on lag and interval width, a calibrated six-coefficient region "
     "is validated on independent simulations, and at six evaluated lags the recorded carrier coefficients show an excess that the leading physical drive, "
     "with or without a single relaxation, does not reproduce. No relaxation is detected, and failed validation gates are reported."),
])

apply('output/submission/submission-checklist-prd.md', [
    ("이를 바로잡자 J0337 물리 구동 결론이 바뀌었다. 여섯 계수 초과를 물리 구동이 어느 지연에서도 재현하지 못한다.",
     "이를 바로잡자 J0337 물리 구동 결론이 바뀌었다. 평가한 여섯 지연에서 여섯 계수 초과를 물리 구동이 재현하지 못한다(2일은 근소하고, β≥0이면 확실하다). "
     "세 심사자의 확인 심사는 모두 \"경미 수정 후 수락\"이었고, 그 지적도 반영했다(`notes/REQUEST293_CONFIRMATION_REVIEW_KO.md`)."),
])
