"""Freeze Phase59 component decisions without promoting a failed photon input."""
from pathlib import Path
import json
import subprocess
import def_photon_partition_audit as audit

ex=audit.ex;h=audit.h;ROOT=h.ROOT;OUT=audit.OUT
NAMES=['native-groups','loss-bounds','line-resolved','explicit-grid','edge-repair',
       'display-control','display-fixed','matter-split','partition-audit']


def read(path):return json.loads(path.read_text())


def main():
    assert not (OUT/'milestone-result.json').exists()
    dirs=[OUT.parent/('def-photon-'+name) for name in NAMES]
    bindings=0
    for folder in dirs:
        for path,sha in read(folder/'plan.json')['bindings'].items():
            assert h.digest(ROOT/path)==sha,path
            bindings+=1
    assert audit.proof()['passed']
    control=OUT.parent/'def-photon-display-control'
    gray=[]
    for label in ['original','outward']:
        lines=(control/(label+'.txt')).read_text().splitlines()
        start=next(i for i,s in enumerate(lines) if s.startswith('Rosseland and Planck'))
        gray.append(lines[start:start+3])
    assert gray[0]==gray[1]
    ex.write(control/'gray-replay.json',dict(classification='Counterexample candidate',passed=True,gray_tokens_identical=True))
    # These are retrospective records of observed failures, not preregistrations.
    failures={
        'explicit-grid':dict(action='pilot submit',error='HTTP Error 403: Forbidden',submitted_queries=1,
            unrun_queries=15,raw_response_body_saved=False,interpretation='Long explicit payload rejected; no reason for rejection established.'),
        'edge-repair':dict(action='pilot reader',expected_groups=2,returned_groups=1,submitted_queries=1,
            unrun_queries=15,interpretation='Output-window filtering established by the separate same-cookie display intervention.'),
        'line-resolved':dict(action='reader interpretation correction',expected_groups=999,returned_groups=998,
            interpretation='Generation-loss hypothesis withdrawn. Display-fixed replay recovered every group and preserved all old rows; original records and reader-failure.json retained.')}
    for name,record in failures.items():
        ex.write(OUT.parent/('def-photon-'+name)/'failure-disposition.json',
            dict(classification='Counterexample candidate',retrospective=True,**record))
    # Preserve the two actual alarm wrappers used before the retrieval code
    # acquired its own alarm. These add no new network request.
    for name,module,owner in [('loss-bounds','def_photon_loss_bounds','m'),('line-resolved','def_photon_line_resolved','m.previous')]:
        (OUT.parent/('def-photon-'+name)/'bounded-fetch-launcher.py').write_text(
            f'import signal\nimport {module} as m\noriginal={owner}.FunctionType\n'
            'def bounded(code,env,argdefs=None):\n    fn=original(code,env,argdefs=argdefs)\n'
            '    def invoke(row):\n        signal.alarm(20)\n        try:return fn(row)\n        finally:signal.alarm(0)\n'
            f'    return invoke\n{owner}.FunctionType=bounded\nm.fetch()\n')
    native=read(OUT.parent/'def-photon-native-groups/result.json')
    matter=read(OUT.parent/'def-photon-matter-split/result.json')
    partitions=read(OUT/'result.json');fine=read(OUT.parent/'def-photon-display-fixed/result.json')
    assert native['passed'] and matter['passed'] and not partitions['passed'] and not fine['passed']
    times=[read(p)['seconds'] for p in (OUT.parent/'def-photon-native-groups').glob('*-retrieval.json')]
    times += [read(OUT.parent/('def-photon-'+name)/'retrieval.json')['seconds'] for name in ['loss-bounds','line-resolved','display-fixed','partition-audit']]
    times += [read(control/'result.json')['seconds']]
    result=dict(classification='Counterexample candidate',progress_class='loophole progress; conditional theorem progress',
        decision='PARTIAL_INPUT_PROGRESS_PHOTON_SOURCE_REJECTED',
        current_background_matter_radiation_split_passed=True,native_gray_recombination_passed=True,
        native_group_partition_additivity_passed=False,local_photon_dynamic_input_accepted=False,
        conditional_moment_theorem_passed=True,display_filter_repair_passed=True,
        angular_gain_kernel_complete=False,material_energy_exchange_evolved=False,
        physical_surface_photon_flux_closed=False,heat_velocity_time_order_passed=False,
        whole_star_heat_closed=False,full_dynamic_charge_solved=False,
        frozen_plan_bindings_verified=bindings,recorded_network_seconds=sum(times),
        network_time_excludes='Two pilots without stored duration, preparation, offline analysis and documentation.',
        total_submit_attempts=73,native_EOS_calls=6,new_stellar_steps=0,
        source_limitation='The fine partitions fail a necessary common-measure condition. Their conditional loss envelopes are not accepted bounds for actual photon transport.',
        next_bottleneck='Obtain a partition-consistent absorption/scattering representation and reproduce both means without normalization; then close gain, material exchange and physical atmosphere flux before long GR evolution.')
    ex.write(OUT/'milestone-result.json',result)
    sources=[ROOT/'verification'/('def_photon_'+name.replace('-','_')+'.py') for name in NAMES]+[Path(__file__)]
    artifacts=sorted(sources+[p for folder in dirs for p in folder.iterdir() if p.is_file()])
    ex.write(OUT/'manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():h.digest(p) for p in artifacts}))
    docs=[ROOT/'docs'/name for name in ['model-definition.md','observable-targets.md','adiabatic-limit.md','nonadiabatic-regime.md','failure-ledger-dynamic-chi.md','dynamic-charge-completion.md']]
    report=ROOT/'notes/REQUEST59_PHOTON_INPUTS_KO.md'
    global_manifest=OUT.parent/'gr-photon-inputs-milestone-manifest.json'
    ex.write(global_manifest,dict(classification='Counterexample candidate',checkpoint='0f6b09ed',
        prior_milestone_manifest_sha256=h.digest(OUT.parent/'gr-radial-conduction-milestone-manifest.json'),
        decision=result['decision'],full_dynamic_charge_solved=False,
        sha256={p.relative_to(ROOT).as_posix():h.digest(p) for p in docs+[report,OUT/'manifest.json']}))
    paper=ROOT/'paper/revision-manifest.json';before=paper.read_bytes()
    assert before==subprocess.check_output(['git','show','HEAD:paper/revision-manifest.json'],cwd=ROOT)
    old=json.loads(before);assert len(old)==167
    entry=dict(classification='Counterexample candidate',progress_class=result['progress_class'],
        matter_radiation_split_passed=True,native_gray_recombination_passed=True,
        fine_photon_dynamic_input_accepted=False,source_partition_additivity_failed=True,
        physical_surface_flux_closed=False,whole_star_heat_closed=False,full_dynamic_charge_solved=False,
        report=report.relative_to(ROOT).as_posix(),report_sha256=h.digest(report),
        evidence_manifest=global_manifest.relative_to(ROOT).as_posix(),evidence_manifest_sha256=h.digest(global_manifest),
        next_bottleneck=result['next_bottleneck'])
    pos=before.rfind(b'}');prefix=before[:pos].rstrip()
    addition=json.dumps({'request59_photon_inputs':entry},ensure_ascii=False,indent=2)[1:-1].strip('\n')
    after=prefix+b',\n'+addition.encode()+b'\n'+before[pos:]
    assert after.startswith(prefix) and all(json.loads(after)[k]==v for k,v in old.items())
    paper.write_bytes(after)
    print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)
    print('BOUND ARTIFACTS',len(artifacts),'PRESERVED PAPER ENTRIES',len(old),flush=True)


if __name__=='__main__':main()
