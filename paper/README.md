# Unified manuscript

The current paper is **Identifying a relaxing internal state in free-fall timing: finite-frequency boundaries and an application to PSR J0337+1715**.

- `manuscript.md` is the focused paper (28 September 2026); `main.tex` is generated.
- `supplement.md` is its Supplemental Material: the previous full manuscript (*Static response and dynamical identifiability in free-fall tests*), with the periastron correction, aligned claim wording, thermal-cutoff qualification and public-replay instructions. It keeps its own numbering, and the paper cites it as "SM Section 4.6" and so on. `supplement.tex` is generated with `python paper/build_manuscript.py supplement.md supplement.tex`.
- `references.bib` supplies numbered BibTeX citations.
- `figures/conditional-intervals.pdf` plots the stored table without refitting (paper and SM).
- `figures/coverage-validation.pdf` plots the registered linear-model validation (SM).
- `figures/comparator-phase-validation.pdf` plots comparator information and phase envelopes (SM).
- `figures/remaining-levers-validation.pdf` plots independent joint-region inclusion and measured transient effects (SM).
- `revision-manifest.json` records input SHA-256 digests. `submission-manifest.json` binds the current submission files; older phase manifests bind their historical revisions.
- `../output/pdf/free-fall-identifiability.pdf` and `../output/pdf/free-fall-identifiability-supplement.pdf` are the rendered paper and SM.
- SM Section 4.6 (the worked white-dwarf calculation) was integrated on 28 September 2026 after independent review, with Pandoc 3.11. The PDFs and the journal source archive below are built with Tectonic 0.17.0 for the Physical Review D submission; `../output/submission/` also holds the cover letter, a plain-text abstract and the submission checklist. Do not run `package_revision.py` as is: it rewrites `revision-manifest.json` from its fixed 9 September list and would drop the later bindings. Package with `../outputs/direct-eos-gr33/native-focused-manuscript/scripts/phase291-package.py` instead.
- The final independent review (28 September 2026) found that the physical drive had used the wrong periastron convention. The released Nutimo code uses the ELL1 convention, eta = e sin(varpi) and kappa = e cos(varpi), and `symbolic/physical_matching.py` and `verification/physical_drive_completion.py` now follow it. The affected outputs were regenerated, and the superseded versions are in `../outputs/research-completion/withdrawn-periastron-convention/`. `verification/simultaneous_inference.py validate` recomputes the registered simulations together with the data section. Its simulated noise is drawn in an eigenbasis that is not unique for repeated eigenvalues, so under a different linear-algebra configuration (for example `OPENBLAS_NUM_THREADS=1` instead of 4) its counts agree with the registered run only within Monte Carlo error; a fixed configuration reproduces itself exactly. The committed file therefore keeps the registered simulation rows (see its `data_revision` field).
- `../output/submission/free-fall-identifiability-source.zip` contains the journal source of the paper: `main.tex`, `main.bbl`, `references.bib`, its one figure and a README; it compiles on its own.
- `../docs/unified-revision-2026-09-09.md` maps the review issues to corrections and remaining limitations.

The separate Paper A/B sources and builders are historical snapshots, superseded for submission. Do not combine their previous headlines with this revision's conclusions. The root Makefile contains symbolic targets; use the explicit commands below to build both submission documents.

`../docs/remaining-levers-2026-09-09.md` is a historical lever ledger. The current submission scope and qualifications are in the paper, SM and `../notes/REQUEST295_SUBMISSION_CORRECTIONS_KO.md`. The conditional inference and worked stellar model do not establish complete nonlinear pulse/noise inference or overall physical closure.

## Request 12 reproduction

From the repository root, prefix these commands with `rtk proxy` on the configured Windows host:

```bash
python symbolic/state_identifiability.py
python verification/physical_drive_completion.py
python verification/simultaneous_inference.py calibrate
python verification/simultaneous_inference.py validate
python verification/gap_pair_audit.py
python verification/runtime12_analysis.py
```

The original calibration was committed before validation. Reproducing the fixed seeds is not a new independent validation. These analyses require NumPy/SymPy and reuse the existing fixed-array helpers. Runtime adjudication reads the committed full return vectors; it does not launch the engine.

Live reproduction uses the existing external WSL release in `~/work/nutimo_pilot`. `verification/prepare_runtime12.py` creates new source/run copies and refuses to replace an existing build. Its exact compile command is in `outputs/research-completion/runtime12/build.json`; the source patch is alongside it. Set `LD_LIBRARY_PATH` to the new `nutimo_request12/src`, `TEMPO2` to the existing `install/third_party/tempo2`, and `OMP_NUM_THREADS=2`. Run `runtime_completion.py transient`, `runtime_completion.py jac W 4` for W=0,1,2,3, then `runtime_completion.py nonlinear` after the gap inputs exist. Each producer writes to `outputs/research-completion/runtime12` and records hashes. Jacobian workers use separate run directories. Existing return files are checkpoints: use a fresh output directory for an independent reproduction, preserving the committed records. The nonlinear producer explicitly rejects the originally proposed e>=1 fractions and runs only the registered local follow-through.

## Build

Python with the repository's NumPy/SymPy dependencies and Matplotlib, Pandoc, and a TeX distribution are needed. From the repository root:

```bash
python paper/build_figures.py
python paper/build_manuscript.py
python paper/build_manuscript.py supplement.md supplement.tex
cd paper
latexmk -pdf -interaction=nonstopmode -halt-on-error -output-directory=build main.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error -output-directory=build supplement.tex
```

Alternatively, from the root use `tectonic --keep-logs --keep-intermediates --outdir paper/build paper/main.tex` (and the same for `paper/supplement.tex`) after creating `paper/build`. The builder accepts an absolute Pandoc executable path through the `PANDOC` environment variable. `pypandoc-binary` supplies one when Pandoc is not installed. The committed TeX and figure can be compiled directly without Pandoc or Python.

Also compile the cover letter with `tectonic --outdir output/submission output/submission/cover-letter-prd.tex`. After compiling, run `python outputs/direct-eos-gr33/native-focused-manuscript/scripts/phase291-package.py` to copy both PDFs, write the plain-text abstract, assemble the portable source ZIP and refresh the submission hashes. It resolves its own checkout, independently of the working directory. Existing research-input bindings are retained. The archive omits the bulky timing arrays; the manifest identifies them in the repository. `paper/package_revision.py` is kept only as a historical record (see above).

## Public replay

The submission snapshot is tagged `prd-submission-2026-09-28`. A public clone and NumPy suffice for:

```bash
python verification/replay_public_inference.py
```

Status: Imported from prior work. The small `outputs/research-completion/public-inference.json` record holds the six carrier coefficients, precision matrix, physical-drive parameters, evaluated lags and threshold. The script independently forms the two drive columns and solves the whitened least-squares problem, reproducing the omnibus statistic and all six sections, including the nonnegative-beta constraint. This is a replay conditional on the exported fit, not a timing-engine run, derivative certificate or new coverage calibration. `--export` regenerates the record only when the private frozen arrays are available; fixed seeds and the calibration remain unchanged.

## Bounded verification

The full first command below requires omitted binary inputs (including the baseline, 35 Jacobians and all manifest-bound arrays). It cannot finish on the public snapshot alone. Use the public replay above for the compact conditional inference. Analytic-only scripts can be run separately without those arrays.

```bash
python verification/verify_unified_paper.py
python symbolic/checks/test_symbolic_smoke.py
python verification/tier1_survivor_exact.py
python symbolic/chi_relaxation_response.py
python symbolic/chi_two_frequency_response.py
python verification/verify_ce_a3a4.py
python symbolic/frequency_sweep_distinguishability.py
python verification/verify_identities.py
python verification/verify_survivors.py
```

Status: Proven. These checks cover finite-order regularity, exact interpolation and its obstruction, the declared static catalog, ODE and linear-sideband identities, and table arithmetic. The manifest check verifies bytes when the manifest is present; it does not establish statistical coverage or provenance before those bytes were recorded.

Status: Imported from prior work. Timing outputs and failed/corrected gates in `request10_external/` are reused records. None of these commands performs a new Nutimo fit or simulation. The frequency-sweep script regenerates its analytic JSON/TSV summaries only.

## Ordered research completion

The 1/2/3 follow-through is recorded in `../notes/REQUEST11_1_NUISANCE_AUDIT_RESULT.md`, `../notes/REQUEST11_2_COVERAGE_RESULT.md` and `../notes/REQUEST11_3_MATCHING_RESULT.md`. Designs precede their respective runs; the estimated-covariance extension was explicitly registered after the first coverage outcomes.

```bash
OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 python verification/nuisance_audit.py
OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 python verification/coverage_audit.py
OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 python verification/estimated_covariance_audit.py
python symbolic/physical_matching.py
OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 python verification/comparator_audit.py
OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 python verification/phase_state_audit.py
OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 python verification/phase_refinement.py
```

In PowerShell, set `$env:OPENBLAS_NUM_THREADS='4'` and `$env:OMP_NUM_THREADS='4'` before the Python commands. Outputs live in `../outputs/research-completion/`; simulations use fixed seeds and stored linear arrays. These commands do not launch the timing engine. The nuisance audit checks only frozen REQUEST10 entries in the manifest so it can precede manuscript regeneration. Refresh the revision manifest after deliberate source/output changes before running the complete paper check.

Status: Imported from prior work. The audit retains the full 90-direction nuisance span. Estimated-covariance simulations give minimum U coverage 94.59% within the prescribed family; this is not an astrophysical likelihood certification.

Status: Proven. Force-level matching now gives tau=Gamma/kappa and beta=Ustar*a_w^2/(kappa*m_p) under explicit stable-branch/equal-charge/reduction assumptions. The leading physical-drive phase closure rejects the historical auxiliary mapping, so its physical-beta interpretation is withdrawn.

Status: Conjectural. Numerical EOS-to-body matching, corrected physical-drive inference and a validated astrophysical likelihood remain outside the conditional paper. This revision does not claim a universal empirical SEP exclusion.

The final two requested levers are recorded in `../notes/REQUEST11_4_COMPARATOR_RESULT.md` and `../notes/REQUEST11_5_PHASE_STATE_RESULT.md`; the current completion report is `../docs/levers-4-5-completion-2026-09-09.md`. Comparator testing includes a physically specified reciprocal fast-relaxation class. Expanded phases include an analytic upper envelope and a separately registered local refinement. Initial-state insufficiency and the weak long-gap stress are reported as boundaries, not successful universal validation.
