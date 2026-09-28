"""Publish the completed175 trajectories and their failed reciprocal gate."""
from pathlib import Path
import hashlib,importlib.util,json,shutil,sys

helper=Path(__file__).with_name('.phase174-publish.py')
if not helper.exists():helper=Path(__file__).with_name('previous-publication-helper.py')
spec=importlib.util.spec_from_file_location('previous_publication',helper)
previous=importlib.util.module_from_spec(spec);spec.loader.exec_module(previous)
root,runtime,master=previous.root,previous.runtime,previous.master
read,write,sha,git=previous.read,previous.write,previous.sha,previous.git
out=root/'outputs/direct-eos-gr33/native-pressure-reciprocal'
manifest=out.parent/'native-pressure-reciprocal-manifest.json'
prefix_path=root/'.phase175176-doc-prefixes.json'
work=runtime/'native-pressure-reciprocal175-work'
modules=['return_native_pressure_reciprocity.py','complete_native_pressure_reciprocity.py','read_native_pressure_charge.py']
notes=['REQUEST175_ACTUAL_MATERIAL_RECIPROCITY_KO.md','REQUEST176_SAME_SOLUTION_CHARGE_READOUT_KO.md']


def package():
    assert not manifest.exists() and not (out/'final-result.json').exists()
    r=read(work/'result.json');b=read(work/'block-result.json');photon=read(work/'photon-result.json')
    material=read(work/'sweep-1/material/production.json');local=read(work/'mechanical-localization.json')
    assert not r['passed'] and not b['passed'] and not r['reciprocal_block_accepted']
    assert photon['passed'] and photon['full_photon_horizon_completed'] and material['passed']
    assert 'Actual finite reciprocal block' in read(work/'closure-block_resume-receipt.json')['error']
    assert not (runtime/'native-pressure-charge176-work').exists()
    assert read(work/'charge-collector-code-check.json')['passed']
    for name in modules:assert sha(root/'verification'/name)==sha(runtime/'verification'/name)
    prefixes=read(prefix_path);old=read(out.parent/'native-pressure-matter-manifest.json')
    preserved={p:h for p,h in old['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    out.mkdir(exist_ok=True);copied={};reused={};intermediate={}
    for src in sorted(work.rglob('*')):
        if not src.is_file():continue
        rel=src.relative_to(work);key=rel.as_posix()
        if rel.parts[:1]==('sweep-0',) and src.suffix=='.npz':
            parent='native-pressure-return/completed/sweep-1/photons' if rel.parts[1]=='photons' else 'native-pressure-matter/completed/sweep-1/material-flux-tangent'
            oldfile=out.parent/parent/src.name;assert sha(src)==sha(oldfile),src
            reused[key]=dict(path=oldfile.relative_to(root).as_posix(),sha256=sha(src));continue
        if src.name.startswith('interval-') and src.suffix=='.npz':intermediate[key]=sha(src);continue
        dst=out/'completed'/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst)
        assert sha(dst)==sha(src);copied[key]=sha(dst)
    sources=[(Path(__file__),'publication-producer.py'),(Path(previous.__file__),'previous-publication-helper.py'),
        (root/'.phase175-followthrough.py','stopped-followthrough-producer.py'),
        (root/'.phase175-mechanical-inspect.py','mechanical-localization-producer.py'),
        (root/'.phase176-codecheck.py','charge-collector-check-producer.py'),
        (runtime/'.phase176-metric-preflight.json','metric-preflight.json'),
        (runtime/'.phase176-pre-metric-producer.py','pre-metric-charge-producer.py')]
    for src,name in sources:shutil.copyfile(src,out/name)
    write(out/'publication.json',dict(copies=copied,reused_published_files=reused,
        runtime_only_accepted_intermediate_restart_archives=intermediate,
        scope='Full175 photon/material trajectories, prefixes/final checkpoints, failed gates and original producers are byte-identical. Intermediate photon restart archives remain runtime-local and hash-listed. The176 charge producer is prepared but NOT executed; only its actual-object preflight and code regression exist.'))
    used=sum(read(f)['seconds'] for f in work.glob('*-receipt.json'))
    final=dict(classification='Counterexample candidate',passed=False,
        actual_previous_material_returned_to_photons=True,full_photon_horizon_completed=True,
        full_free_material_return_completed=True,reciprocal_block_accepted=False,
        maximum_block_defect=b['maximum_block_defect'],block_rows=b['rows'],
        photon_time_comparison=photon['time_comparison'],material_time_comparison=material['time_comparison'],
        mechanical_localization=local['rows'],original_total_action_budget_seconds=5950,action_seconds_used=used,
        charge_collector_cache_defect_identified=True,charge_collector_code_repaired_and_checked=True,
        corrected_GR_charge_readout_executed=False,previous_selected_source_geometry_revalidated=False,
        final_charge_conclusion='unadjudicated',physical_final_charge_solved=False,full_goal_complete=False)
    write(out/'final-result.json',final)
    tails={
        'model-definition':'수정 자유 물질을 실제 광자 방정식에 반환하고, 새 광자 이력을 같은 자유 물질에 적용해 두 전체 경로를 완료했다. 그러나 비충돌 수소 수송 M_H와 그 변화율의 상호 입력 잔차가 원0.2% 기준을 넘었다. 결합 해 수락과 최종 전하 유지 여부는 미판정이다.',
        'observable-targets':'단계175의 실제 상호 반환 뒤 M_H 입력 기준이 실패했으므로 예정한176 GR·최종 전하 판독은 실행하지 않았다. 작은 광자/물질 E/H 차이나 단계 원천 잔차를 전체 상호 결합 수락으로 대체하지 않는다.',
        'adiabatic-limit':'M=H-C의 충돌 상쇄 항등식은 유지되지만 원 비충돌 수소 입력이 반복 후 바뀌었다. 이 입력 실패를 정적 흡수·no-go 결과로 해석하지 않으며 최종 전하 결론은 미판정이다.',
        'nonadiabatic-regime':'동일 전체 기간의 실제117/227 Radau 하위 단계와 반환 물질 경로를 완료했다. 모든688 단계의 실제 원천 검사를 수행했지만 M_H 자체와 시간 변화 입력이 원 기준에 미달했다. 새 시간 격자나 자동 추가 반복을 시작하지 않았다.',
        'failure-ledger-dynamic-chi':'원0.2% 상호 입력 기준에서 M_H 약20.71%, dM_H 약42.14%로 실패했다. 원 단계 원천 차이는 최대0.03242%이나 이것으로 별도의 실패 기준을 지우지 않는다. 저장 배열의 부동소수점 규모 대조로 단순 차감 반올림 설명을 배제했다. GR 수집기의 생성자 계량 캐시 결함도 발견·코드 수정했으나 실제 GR 실행은 이 실패 때문에 보류했다.',
        'dynamic-charge-completion':'현재 지배 오차를 해결한 동일 결합 해의 최종 전하 결론이라는 사용자 기준은 아직 충족하지 못했다. 광자·물질 전체 반환을 실행했지만 비충돌 수소 상호 입력 기준이 실패했다. physical_final_charge_solved와full_goal_complete는false이며176 물리 판독은 미실행이다.'}
    for name,text in tails.items():
        file=root/f'docs/{name}.md';prior=prefixes[file.relative_to(root).as_posix()]
        assert file.stat().st_size==prior['bytes'] and sha(file)==prior['sha256']
        with file.open('ab') as f:f.write(('\n\n## 단계175 — 실제 상호 반환 완료와 수소 입력 실패\n\n분류: Counterexample candidate. '+text+' [실행 및 원 판정](../notes/'+notes[0]+'), [미실행 GR 판독](../notes/'+notes[1]+').\n').encode())
    report='\n\n## 전체 경로 결과와 수락 실패\n\n분류: Counterexample candidate. **최종 전하의 기존 결론 유지 여부는 미판정이다.** 실제 이전 자유 물질 입력을 반환한 광자 전체 기간과 그 광자를 소비한 자유 물질 전체 기간을 완료했다. 그러나 상호 결합의 원0.2% 입력 기준이 실패했다.\n\n| 항목 | 64 시계 | 128 시계 |\n|---|---:|---:|\n'
    for key,label in [('M_H','비충돌 수소 M_H 차이'),('dM_H','비충돌 수소 변화율 차이'),('B','바리온 입력 차이')]:
        report+=f"| {label} | {b['rows'][0]['lagged_inputs'][key]*100:.8g}% | {b['rows'][1]['lagged_inputs'][key]*100:.8g}% |\n"
    report+=f"\n분류: Counterexample candidate. 광자 여섯 시간 차이는{photon['time_comparison']}, 물질 B/S/E/H 시간 차이는{material['time_comparison']}이다. 모든688 실제 Radau 단계의 원천 차이는 최대{max(v for row in b['rows'] for v in row['stage_forcing_defect'])*100:.8g}%이고 광자/물질 E/H 잔차는 최대{max(v for row in b['rows'] for v in row['paired_E_H']):.8g}다. 이 작은 값으로 별도로 등록한 M_H 실패를 취소하지 않는다.\n"
    x=local['rows'][1]
    report+=f"\n분류: Counterexample candidate. 저장 배열만의 대조에서128시계의 M_H 최대 L1은{x['new_mechanical_H_L1']:.9g}, 변화 L1은{x['change_L1']:.9g} 중성수소 수였다. 전체 H 최대 L1은{x['total_H_L1']:.9g}이며 binary64 반올림 척도/관측 변화는{x['float64_scale_over_change']:.9g}다. 약20.7% 차이는 단순 차감 반올림 규모보다 훨씬 크다. 이 검사는 원인 전체나 수렴률의 증명은 아니다.\n"
    report+='\n분류: Counterexample candidate. 실제688 단계 검사의2배 여유 비용 예상428.832초가 원300초 상한을 넘어 먼저 중단됐다. 거절16.890초를 보존하고 완료된 물질 생산에서 미사용150초를 옮겨 같은 총5950초 내에서430초 상한으로 재개했다. 물리 수락 기준·진화 해·단계 수는 바꾸지 않았다. 그 결과 위 상호 입력 실패가 확정됐고 자동 추가 sweep은 시작하지 않았다.\n'
    report+=f'\n기록된175 action 벽시간 합계는{used:.3f}초다. 176 GR 실행은 없다. 모든 원 실패·코드·실행 내역은 [산출물](../outputs/direct-eos-gr33/native-pressure-reciprocal/final-result.json)에 결속했다.\n'
    report+='\n분류: Conjectural. 다음 결정 대상은 비충돌 수소 수송의 상호 결합이다. 현재 전체 해를 재사용해 해당 수송 블록을 직접 결합하거나, 원 기준을 만족할 비용이 산정된 재설계만 진행해야 한다. 측정되지 않은 추가 반복을 연속 실행하거나 시간 간격을 자동 세분화하지 않는다. 수락 뒤에만176의 동일 해 원천·에너지·경계 전하 판독으로 넘어간다.\n'
    with (root/'notes'/notes[0]).open('ab') as f:f.write(report.encode())
    report='\n\n## 선행 조건의 실제 실패 — 판독 미실행\n\n분류: Counterexample candidate.175 전체 경로는 완료됐지만 M_H·dM_H 상호 입력이 원 기준에 미달했다. 따라서 이 노트의176 물리 판독은 실행하지 않았으며 새 전하값이나 부호 판정은 없다. 생성자 계량 불일치의 실제 객체 확인과 시간별 필드 소비 코드 회귀 검사만 완료했다. 기존 선택 전하의 원천 계량도 아직 재검증되지 않았다.\n'
    with (root/'notes'/notes[1]).open('ab') as f:f.write(report.encode())
    old_master=sha(master);write(out/'preservation.json',dict(previous_phase174_nondoc=preserved,document_prefixes=prefixes,previous_master_sha256=old_master))
    files=[root/'verification'/name for name in modules]+[root/'notes'/name for name in notes]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(previous_master_sha256=old_master,preserved_document_prefixes=prefixes,sha256={p.relative_to(root).as_posix():sha(p) for p in files})
    write(manifest,final);m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_pressure_reciprocal_charge']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


def check(mode):
    m=read(manifest);master_data=read(master)
    for p,h in m['sha256'].items():assert sha(root/p)==h and master_data['sha256'][p]==h,p
    for p,h in read(out/'preservation.json')['previous_phase174_nondoc'].items():assert sha(root/p)==h,p
    for p,v in m['preserved_document_prefixes'].items():assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest()==v['sha256'],p
    p='docs/dynamic-charge-completion.md';assert (root/p).read_bytes().splitlines()[:20]==git('show','d5744f766719f0a29c7acb913c7984bdfd3624e2:'+p).splitlines()[:20]
    assert master_data['sha256'][manifest.relative_to(root).as_posix()]==sha(manifest)
    paths=list(m['sha256'])+[manifest.relative_to(root).as_posix(),'paper/revision-manifest.json']
    Path('C:/Users/lpaiu/AppData/Local/Temp/native-pressure-reciprocal-175176-paths').write_bytes(b'\0'.join(p.encode() for p in paths)+b'\0')
    if mode=='head':
        for p in paths:assert hashlib.sha256(git('show','HEAD:'+p)).hexdigest()==sha(root/p),p
    print(json.dumps(dict(bound_files=len(m['sha256']),published_paths=len(paths),prefixes_preserved=6,reciprocal_block_accepted=False,charge_readout_executed=False,final_charge_conclusion='unadjudicated')))


if __name__=='__main__':
    action=sys.argv[1]
    if action=='package':package()
    check(action)
