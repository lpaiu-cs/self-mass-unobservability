"""Append the actual terminal recovery outcome without erasing its admission."""
from pathlib import Path
import importlib.util,shutil
spec=importlib.util.spec_from_file_location('p',Path('.phase213215-publish.py'));p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p)
work=p.runtime/'native-coordinate-photon214-work';out=p.out;manifest=p.read(p.manifest)
assert not (out/'terminal-amendment.json').exists()
status=p.read(work/'resume-controller-status.json');assert status['state']=='failed'
row=p.read(work/'clock-128/snapshot-02.json');assert not row['passed'] and row['radial_port_relative']>1e-12
shutil.copyfile(p.manifest,out/'published-entry-manifest.json');shutil.copyfile(p.note,out/'published-entry-note.md')
for name in ['resume-controller-status.json','fine_resume-receipt.json','fine_resume.stderr.log','progress-64.json','progress-128.json','clock-128/snapshot-02.json','clock-128/original-equation-15.json','accepted-64.npz','accepted-128.npz']:
    dst=out/'terminal'/name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(work/name,dst);assert p.sha(dst)==p.sha(work/name)
shutil.copyfile(__file__,out/'terminal-publication-producer.py')
message='분류: Counterexample candidate. 재개214는coarse11/fine15개 광자 복원 단계를 수락한 뒤 fine16번째 쌍을 포함한T/16출력의 누적 반경 출구 상대 차이4.740709463e-12>1e-12로 실패했다. 해당 쌍의 원 결합식·native비트 일치는 통과했고 광자 끝점1.29139e-16,물질수지 최대1.17006e-13,각도출구1.72126e-16도 원 기준을 통과했다. 실패 원인은 벽시간·메모리 한도가 아니며, 반경 출구 기준을 완화하거나 전체 복원을 수락하지 않는다. 원 출구 누적 순서와 작은 복원 오차의 영향을 아직 분리하지 않았으므로 원인을 단정하지 않는다. 두 자식 프로세스는 종료됐고 마지막 수락 체크포인트를 보존했다. 본210의116단계 수락·117실패, 최종 전하 미판정과 전체 완료 조건도 유지한다.'
with p.note.open('a',encoding='utf-8') as f:f.write('\n\n'+message+'\n')
for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
    with (p.root/f'docs/{name}.md').open('a',encoding='utf-8') as f:f.write('\n\n분류: Counterexample candidate. 214재개는원T/16누적반경출구 기준4.741e-12>1e-12로 종료됐다. 광자끝점·물질수지·각도출구와 개별 결합식/native동일성 통과로 이를 대체하지 않는다. 마지막coarse11/fine15복원 체크포인트를 보존했고 최종 전하는 미판정이다. [종료 근거](../notes/'+p.note.name+').\n')
amendment=dict(classification='Counterexample candidate',resume214_terminal_failed=True,radial_port_gate_failed=True,snapshot=row,workers_finished=True,original_entry_manifest_sha256=p.sha(out/'published-entry-manifest.json'),original_gates_unchanged=True,final_charge_conclusion='unadjudicated',full_goal_complete=False)
p.write(out/'terminal-amendment.json',amendment)
result=p.read(out/'result.json');result.update(recovery214_terminal_resume_result_included=True,resume214_terminal_failed=True);p.write(out/'result.json',result)
manifest.update(recovery214_terminal_resume_result_included=True,resume214_terminal_failed=True)
for f in [p.note,p.module,*[p.root/k for k in manifest['document_prefixes']],*[f for f in out.rglob('*') if f.is_file()]]:manifest['sha256'][f.relative_to(p.root).as_posix()]=p.sha(f)
p.write(p.manifest,manifest);master=p.read(p.master);master['sha256'].update(manifest['sha256']);master['sha256'][p.manifest.relative_to(p.root).as_posix()]=p.sha(p.manifest);master['native_coordinate_recovery']={k:v for k,v in manifest.items() if k!='sha256'};p.write(p.master,master);p.check()
