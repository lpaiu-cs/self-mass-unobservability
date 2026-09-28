"""Bind the completed, rejected repair attempts without altering prior records."""
from pathlib import Path
import hashlib
import json
import subprocess

root=Path.cwd();out=root/'outputs/direct-eos-gr33/def-gr-spatial-repair'
digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
sources=[root/'verification'/name for name in ['def_gr_spatial_repair.py','def_gr_energy_factor.py',
    'def_gr_balanced_modes.py','def_gr_spatial_direct.py','def_gr_spatial_banded.py',
    'def_gr_spatial_banded_refined.py','def_gr_hierarchical.py','def_gr_direct_time.py']]
docs=[root/'docs'/name for name in ['model-definition.md','observable-targets.md','adiabatic-limit.md',
    'nonadiabatic-regime.md','failure-ledger-dynamic-chi.md','dynamic-charge-completion.md']]
report=root/'notes/REQUEST76_GR_SPATIAL_REPAIR_KO.md'
manifest=root/'outputs/direct-eos-gr33/gr-spatial-repair-milestone-manifest.json'
paper=root/'paper/revision-manifest.json'
assert not manifest.exists()
for p in docs:
    old=subprocess.check_output(['git','show','HEAD:'+p.relative_to(root).as_posix()])
    assert p.read_bytes().startswith(old),p
previous=json.loads(subprocess.check_output(['git','show','HEAD:paper/revision-manifest.json']))
current=json.loads(paper.read_text());assert current==previous and len(current)==184
paths=sources+docs+[report]+[p for p in out.rglob('*') if p.is_file() and '__pycache__' not in p.parts]
evidence={p.relative_to(root).as_posix():digest(p) for p in sorted(paths)}
result=json.loads((out/'result.json').read_text());assert not result['passed'] and not result['original_failure_resolved']
payload=dict(classification='Counterexample candidate',checkpoint='95d99987',
    progress_class='loophole progress; rejected actual repair candidates',decision=result['decision'],
    prior_milestone_manifest_sha256=digest(root/'outputs/direct-eos-gr33/gr-canonical-evolution-milestone-manifest.json'),
    actual_same_input_applied=True,scalar_direct_time_passed=True,all_four_time_gates_passed=False,
    new_space_comparison_completed=False,original_failure_resolved=False,full_dynamic_charge_solved=False,sha256=evidence)
manifest.write_bytes((json.dumps(payload,ensure_ascii=False,indent=2)+'\n').encode())
current['request76_same_input_spatial_repair']=dict(classification='Counterexample candidate',
    progress_class=payload['progress_class'],actual_same_input_applied=True,scalar_direct_time_passed=True,
    all_four_time_gates_passed=False,new_space_comparison_completed=False,passed=False,original_failure_resolved=False,
    full_dynamic_charge_solved=False,new_EOS_calls=0,
    next_bottleneck='Obtain simultaneous four-component propagation acceptance in the same polynomial GR space, preserving weak scalar and local velocity, then finish the fixed1/2/4 spatial comparison. No mixed-method splicing or automatic refinement.',
    report=report.relative_to(root).as_posix(),report_sha256=digest(report),
    evidence_manifest=manifest.relative_to(root).as_posix(),evidence_manifest_sha256=digest(manifest))
paper.write_bytes((json.dumps(current,ensure_ascii=False,indent=2)+'\n').encode())
assert all(current[k]==v for k,v in previous.items()) and len(current)==185
print(json.dumps(dict(files_bound=len(evidence),paper_entries=185,prior_entries_preserved=184,document_prefixes_preserved=6)))
