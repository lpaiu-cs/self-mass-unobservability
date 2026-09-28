from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import zipfile

ROOT = Path(__file__).resolve().parents[1]
files = [
    "paper/manuscript.md", "paper/main.tex", "paper/references.bib",
    "paper/build_manuscript.py", "paper/build_figures.py", "paper/README.md",
    "paper/figures/conditional-intervals.pdf",
    "paper/figures/coverage-validation.pdf", "paper/package_revision.py",
    "paper/figures/comparator-phase-validation.pdf",
    "verification/comparator_audit.py", "verification/phase_state_audit.py", "verification/phase_refinement.py",
    "outputs/research-completion/comparator-audit.json", "outputs/research-completion/phase-state-audit.json",
    "outputs/research-completion/phase-refinement.json",
    "notes/REQUEST11_4_COMPARATOR_PLAN.md", "notes/REQUEST11_4_COMPARATOR_RESULT.md",
    "notes/REQUEST11_5_PHASE_STATE_PLAN.md", "notes/REQUEST11_5B_PHASE_REFINEMENT_PLAN.md",
    "notes/REQUEST11_5_PHASE_STATE_RESULT.md", "docs/levers-4-5-completion-2026-09-09.md",
    "verification/nuisance_audit.py", "verification/coverage_audit.py",
    "verification/estimated_covariance_audit.py", "symbolic/physical_matching.py",
    "outputs/research-completion/nuisance-audit.json",
    "outputs/research-completion/nuisance-intervals.csv",
    "outputs/research-completion/coverage-audit.json",
    "outputs/research-completion/estimated-covariance-audit.json",
    "outputs/research-completion/physical-matching.json",
    "notes/REQUEST11_1_NUISANCE_AUDIT_PLAN.md", "notes/REQUEST11_1_NUISANCE_AUDIT_RESULT.md",
    "notes/REQUEST11_2_COVERAGE_PLAN.md", "notes/REQUEST11_2B_ESTIMATED_COVARIANCE_PLAN.md",
    "notes/REQUEST11_2_COVERAGE_RESULT.md", "notes/REQUEST11_3_MATCHING_PLAN.md",
    "notes/REQUEST11_3_MATCHING_RESULT.md", "docs/research-completion-2026-09-09.md",
    "README.md", "docs/paper-claims-vs-nonclaims.md", "docs/theorem-package.md",
    "docs/boundary-escape-map.md", "docs/release-note.md",
    "docs/model-definition.md", "docs/observable-targets.md", "docs/adiabatic-limit.md",
    "docs/nonadiabatic-regime.md", "docs/failure-ledger-dynamic-chi.md",
    "verification/verify_unified_paper.py", "verification/tier1_survivor_exact.py",
    "symbolic/nonanalytic_jet_demo.py", "symbolic/chi_relaxation_response.py",
    "symbolic/chi_two_frequency_response.py", "symbolic/frequency_sweep_distinguishability.py",
    "request10_external/scripts/sep_common.py",
    "request10_external/scripts/sep_phase_marg_10_8e.py",
    "request10_external/baseline_planetGR.npz",
    "request10_external/finite_jacobian_v2.npy",
    "request10_external/finite_jacobian_v2_meta.json",
    "request10_external/finite_jacobian.npy", "request10_external/finite_jacobian_meta.json",
    "request10_external/carrier_projection_rank.json",
    "request10_external/sep_dynamic/col_SEP_D.npz",
    "request10_external/sep_dynamic/sep_dynamic_columns.npz",
    "request10_external/sep_dynamic/sep_phase_marg_10_8e.json",
    "request10_external/sep_dynamic/sep_rn_robustness_10_8h.json",
    "request10_external/sep_dynamic/sep_gateG2.json",
    "request10_external/sep_dynamic/sep_gateG2wp.json",
    "request10_external/sep_dynamic/sep_quadrature_overlap_10_8g.json",
    "request10_external/sep_dynamic/turn_search_10_8d.json",
]
pdf = ROOT / "output/pdf/free-fall-identifiability.pdf"
files += [
    'verification/physical_drive_completion.py', 'verification/simultaneous_inference.py',
    'verification/prepare_runtime12.py', 'verification/runtime_completion.py',
    'verification/runtime12_analysis.py', 'verification/gap_pair_audit.py',
    'symbolic/state_identifiability.py', 'docs/remaining-levers-2026-09-09.md',
    'paper/figures/remaining-levers-validation.pdf',
    'outputs/research-completion/corrected-physical-drive.json',
    'outputs/research-completion/simultaneous-calibration.json',
    'outputs/research-completion/simultaneous-validation.json',
    'outputs/research-completion/state-identifiability.json',
    'outputs/research-completion/gap-pair-audit.json',
    'outputs/research-completion/runtime12-analysis.json',
]
files += [p.relative_to(ROOT).as_posix() for p in sorted((ROOT/'notes').glob('REQUEST12_*.md'))]
files += [p.relative_to(ROOT).as_posix() for p in sorted((ROOT/'outputs/research-completion/runtime12').iterdir()) if p.is_file()]
pdf.parent.mkdir(parents=True, exist_ok=True)
shutil.copyfile(ROOT / "paper/build/main.pdf", pdf)
files.append(pdf.relative_to(ROOT).as_posix())
manifest = {
    "revision": "Unified manuscript with remaining-lever inference, live validation and state identifiability, 2026-09-09",
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "reviewed_pre_revision_commit": "3bc2fce",
    "before_task_checkpoint": "14de460",
    "registration_and_result_commits": {"nuisance": ["e790643", "b8d931b", "3506522"], "coverage": ["9412142", "12d9f9a", "6217da2"], "matching": ["1e1cd54", "8edbe84"], "comparator": ["fb269b9", "65f4b65"], "phase_state": ["fb269b9", "5b9dba2", "7636d72"]},
    "remaining_levers_commits": ["abb973e", "6d7f3f3", "150fa43", "d8eb72d", "a08c36a", "05d7af9"],
    "scope": "Conditional theory and inference, plus isolated live transient/derivative and local nonlinear checks. Numerical EOS matching, derivative-error certification, complete pulse reconnection and universal SEP exclusion remain incomplete.",
    "paths_relative_to": "repository root",
    "sha256": {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
}
(ROOT / "paper/revision-manifest.json").write_text(json.dumps(manifest,indent=2)+"\n", encoding="utf-8")
archive = ROOT / "output/submission/free-fall-identifiability-source.zip"
archive.parent.mkdir(parents=True, exist_ok=True)
members = {
    "paper/main.tex":"main.tex", "paper/build/main.bbl":"main.bbl",
    "paper/references.bib":"references.bib", "paper/manuscript.md":"manuscript.md",
    "paper/figures/conditional-intervals.pdf":"figures/conditional-intervals.pdf",
    "paper/figures/coverage-validation.pdf":"figures/coverage-validation.pdf",
    "paper/figures/comparator-phase-validation.pdf":"figures/comparator-phase-validation.pdf",
    "paper/figures/remaining-levers-validation.pdf":"figures/remaining-levers-validation.pdf",
    "paper/revision-manifest.json":"revision-manifest.json",
}
with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as package:
    for source, target in members.items():
        package.write(ROOT/source, target)
    package.writestr("README.txt", "Unified revised manuscript, 9 September 2026.\n\nCompile main.tex with latexmk -pdf main.tex, or tectonic main.tex.\nThe bibliography, generated .bbl and figure are included.\nNo timing engine or Pandoc is needed to compile this TeX source.\n\nrevision-manifest.json uses repository-root paths and identifies the full\nresearch inputs; those timing arrays are not duplicated in this source zip.\nSee the repository paper/README.md for bounded symbolic reproduction.\n\nThis is a local source package, not evidence of journal submission or a DOI deposit.\nAuthor affiliation/contact and the selected journal's formatting must be supplied at submission.\n")
print(f"Packaged {len(manifest['sha256'])} SHA-256 inputs; {pdf}; {archive}")
