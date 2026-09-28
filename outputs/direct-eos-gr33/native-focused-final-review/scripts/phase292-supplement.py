"""Phase 292: correct the periastron convention and its consequences in the Supplemental Material (CRLF kept).

Each replacement anchor must occur exactly once. Values come from the recomputed outputs (compare292.json).
"""
from pathlib import Path
p = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672/paper/supplement.md')
t = p.read_bytes().decode('utf-8')
edits = [
    # SM abstract: state the corrections and map the main text to the SM.
    ("It is the full earlier version of the work; apart from this title and abstract, its text is unchanged. ",
     "It is the full earlier version of the work. Apart from this title and abstract, its text is unchanged except for corrections found in the final review: "
     "the periastron convention of the J0337 drive (Sections 4.4, 5.6, 5.9 and 5.10, Figures 3 and 4, and the scale list of Section 4.6) and the "
     "reproduction notes of the data-availability section. Main-text Sections 2 and 3.1--3.3 correspond to Sections 3.1--3.5 here, main-text "
     "Section 3.4 to Section 5.8, main-text Section 4 to Sections 4.3--4.5, main-text Section 5 to Sections 5.1--5.10, and the white-dwarf paragraph "
     "of the main-text discussion to Section 4.6. "),
    # Section 4.4: convention.
    ("Its table and the frozen parameter values identify this release's eta as \\(e\\cos\\varpi\\) and kappa as \\(e\\sin\\varpi\\); their names must not be interpreted using a different implementation's convention.",
     "The released timing code uses \\(\\eta=e\\sin\\varpi\\) and \\(\\kappa=e\\cos\\varpi\\), the ELL1 convention \\cite{lange2001ell1}; the frozen parameters "
     "then place the pericenters close to those of Ref. \\cite{ransom2014triple}. An earlier version of this analysis read eta as \\(e\\cos\\varpi\\); the final "
     "review found the error, and the values below and in Sections 5.6, 5.9 and 5.10 are corrected."),
    ("The frozen parameters give \\(\\varpi_p=-0.12326821\\), \\(\\varpi_b=-0.10004792\\) and \\(\\mathcal C=3.11837236\\) radians.",
     "The frozen parameters give \\(\\varpi_p=1.69406454\\), \\(\\varpi_b=1.67084424\\) and \\(\\mathcal C=3.16481295\\) radians."),
    ("differ from Equation (\\ref{eq:physical-drive}) by 1.69406454, 1.67084424 and \\(\\pi\\) radians.",
     "differ from Equation (\\ref{eq:physical-drive}) by 0.12326821, 0.10004792 and \\(\\pi\\) radians."),
    # Section 4.6: the physical-drive endpoint no longer exists.
    ("- \\(4.1\\times10^{-10}\\) (physical-drive endpoint)\r\n", ""),
    # Section 5.6: comparator minima with the corrected physical phase set.
    ("is 0.08289 at N=1, \\(8.36\\times10^{-6}\\) at N=2, \\(3.19\\times10^{-6}\\) at N=3, and \\(2.88\\times10^{-6}\\) at N=4. The last case widens the unit-noise standard error by about 589 times. The even-only comparator retains at least 0.003035 in the tested cases.",
     "is 0.005389 at N=1, \\(6.29\\times10^{-5}\\) at N=2, \\(1.17\\times10^{-5}\\) at N=3, and \\(4.63\\times10^{-6}\\) at N=4. The last case widens the unit-noise standard error by about 465 times. The even-only comparator retains at least 0.0311 in the tested cases."),
    ("between \\(9.3\\times10^{-30}\\) and \\(1.4\\times10^{-24}\\) is rounding error",
     "between \\(5.4\\times10^{-30}\\) and \\(5.3\\times10^{-25}\\) is rounding error"),
    # Section 5.9: every physical lag section is empty.
    ("The recorded data's six-coefficient omnibus statistic is 16.3525, above that threshold. Its null is that all six carrier coefficients vanish, not \\(\\beta=0\\) with a physical instantaneous term fitted. Every evaluated physical lag section includes \\(\\beta=0\\). At 2 and 500 days, their absolute-beta endpoints are \\(4.10217\\times10^{-10}\\) and \\(8.31672\\times10^{-9}\\). These conditional intersections also account for mismatch in other coefficient directions, so smaller endpoints than a pointwise interval are not automatically stronger physical constraints. No EOS-matched bound is restored.",
     "The recorded data's six-coefficient omnibus statistic is 16.3525, above that threshold (nominal chi-square p of about 0.012 for six degrees of freedom). "
     "Its null is that all six carrier coefficients vanish, not \\(\\beta=0\\) with a physical instantaneous term fitted. The minimum statistic over the physical "
     "plane is 12.84 at 2 days and 14.39--14.59 at the other five evaluated lags, all above the threshold, so every evaluated physical lag section is empty: "
     "no combination of instantaneous and relaxing response to the leading physical drive reproduces the carrier coefficients at the calibrated level. "
     "No beta interval or EOS-matched bound follows, and the excess is not attributed. With the earlier periastron error, every section contained \\(\\beta=0\\); "
     "that result is withdrawn."),
    # Section 5.10: live factors with the corrected drive.
    ("increases beta standard errors by factors 1.00096--1.00523 across the nine tested lag/covariance combinations.",
     "changes beta standard errors by factors 0.999998--1.00742 across the nine tested lag/covariance combinations."),
    ("standard errors by factors 0.8680--0.8785 and shifts estimates.", "standard errors by factors 0.7858--0.8722 and shifts estimates."),
    # Section 6.
    ("Corrected leading-drive and joint-region calculations address defined conditional questions; the live transient responses cover three lags.",
     "Corrected leading-drive and joint-region calculations address defined conditional questions and reject the leading physical-drive plane at every evaluated lag; the live transient responses cover three lags."),
    # Data availability: withdrawn outputs and reproduction of the registered simulations.
    ("the separately hashed Request 12 returns contain the new live evaluations described in Section 5.10.",
     "the separately hashed Request 12 returns contain the new live evaluations described in Section 5.10. Outputs superseded by the periastron-convention "
     "correction are kept in \\path{outputs/research-completion/withdrawn-periastron-convention/}. Rerunning the registered six-coefficient validation with the "
     "same seeds on another linear-algebra backend reproduces its inclusion rates only within Monte Carlo error, so the registered simulation rows are retained."),
]
for old, new in edits:
    assert t.count(old) == 1, old[:80]
    t = t.replace(old, new)
p.write_bytes(t.encode('utf-8'))
print('applied', len(edits), 'edits; bytes', len(t.encode('utf-8')))
