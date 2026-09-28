"""Bind Phase77 evidence and append its unaccepted verdict, preserving history."""
from pathlib import Path
import json
import hashlib
import subprocess

root=Path.cwd();out=root/'outputs/direct-eos-gr33/def-gr-inverse-projection'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
sources=[root/'verification'/name for name in ['def_gr_inverse_projection.py','def_gr_field_projection.py',
    'def_gr_gauss_time.py','def_gr_multishift.py','def_gr_field_inverse.py','def_gr_inverse_gram.py']]
docs=[root/'docs'/name for name in ['model-definition.md','observable-targets.md','adiabatic-limit.md',
    'nonadiabatic-regime.md','failure-ledger-dynamic-chi.md','dynamic-charge-completion.md']]
report=root/'notes/REQUEST77_GR_COUPLED_REPAIR_KO.md'
manifest=root/'outputs/direct-eos-gr33/gr-coupled-repair-milestone-manifest.json'
paper=root/'paper/revision-manifest.json';assert not manifest.exists()
for p in docs:
    prior=subprocess.check_output(['git','show','HEAD:'+p.relative_to(root).as_posix()])
    assert p.read_bytes().startswith(prior),p
prior=json.loads(subprocess.check_output(['git','show','HEAD:paper/revision-manifest.json']))
current=json.loads(paper.read_text());assert current==prior and len(current)==185
paths=sources+docs+[report]+[p for p in out.rglob('*') if p.is_file() and '__pycache__' not in p.parts]
result=json.loads((out/'result.json').read_text());assert not result['propagation_passed'] and not result['original_failure_resolved']
payload=dict(classification='Counterexample candidate',checkpoint='1260c415',progress_class=result['progress_class'],
    decision=result['decision'],prior_milestone_manifest_sha256=sha(root/'outputs/direct-eos-gr33/gr-spatial-repair-milestone-manifest.json'),
    actual_same_input_applied=True,propagation_passed=False,spatial_comparison_started=False,
    original_failure_resolved=False,full_dynamic_charge_solved=False,
    sha256={p.relative_to(root).as_posix():sha(p) for p in sorted(paths)})
manifest.write_bytes((json.dumps(payload,ensure_ascii=False,indent=2)+'\n').encode())
current['request77_coupled_GR_propagation_repair']=dict(classification='Counterexample candidate',progress_class=payload['progress_class'],
    actual_same_input_applied=True,propagation_passed=False,spatial_comparison_started=False,passed=False,
    original_failure_resolved=False,full_dynamic_charge_solved=False,new_EOS_calls=0,
    direct_velocity_difference_reduction=result['direct_velocity_difference_reduction'],
    next_bottleneck='Separate fast-wave remainder and arithmetic effects using component error estimates on the saved failed histories; retain the common propagation and spatial acceptance boundary without component splicing or automatic refinement.',
    report=report.relative_to(root).as_posix(),report_sha256=sha(report),
    evidence_manifest=manifest.relative_to(root).as_posix(),evidence_manifest_sha256=sha(manifest))
assert all(current[k]==v for k,v in prior.items()) and len(current)==186
paper.write_bytes((json.dumps(current,ensure_ascii=False,indent=2)+'\n').encode())
print(json.dumps(dict(files_bound=len(paths),prior_paper_entries_preserved=185,paper_entries=186,document_prefixes_preserved=6)))
