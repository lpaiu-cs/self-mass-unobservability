from pathlib import Path
import hashlib,json,shutil,subprocess,sys,time

root=Path('E:/lab/self-mass-unobservability')
runtime=Path('//wsl.localhost/Ubuntu-22.04/home/lpaiu/work/native-retained-tail-runtime')
work=runtime/'native-pressure-matter174-work'
out=root/'outputs/direct-eos-gr33/native-pressure-matter'
manifest=out.parent/'native-pressure-matter-manifest.json';master=root/'paper/revision-manifest.json'
read=lambda p:json.loads(Path(p).read_text(encoding='utf-8-sig'))
modules=['resolve_native_material_inverse.py','resolve_native_material_face_response.py','resolve_native_material_flux_tangent.py','continue_native_material_flux_tangent.py']
docs=['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for v in iter(lambda:f.read(1024**2),b''):h.update(v)
    return h.hexdigest()


def write(p,v):Path(p).write_bytes((json.dumps(v,indent=2,ensure_ascii=False)+'\n').encode())
def git(*args):return subprocess.check_output(['rtk','proxy','git',*args],cwd=root)


def audit():
    import numpy as np
    import continue_native_material_flux_tangent as producer
    start=time.monotonic();w=producer.OUT;folder=producer.prior.paths(1)[1]
    final=read(w/'result.json');production=read(folder/'production.json')
    assert final['passed'] and production['passed']
    for name in ['material-continuation-production','material-continuation-residual']:
        assert read(w/f'{name}-receipt.json')['error'] is None
    for name in ['plan','flux-tangent-plan','material-continuation-plan']:
        for p,h in read(w/f'{name}.json')['bindings'].items():assert sha(p)==h,p
    rows=[];hist=[]
    for n in [64,128]:
        with np.load(folder/f'steps-{n}-reference-128.npz') as d,np.load(folder/f'pilot-{n}.npz') as p:
            assert int(d['completed'])==n and d['time']==d['t'][-1]
            for key in ['t','history_scaled','ledgers_scaled','norms_scaled','discards_scaled']:
                assert np.array_equal(d[key][:len(p[key])],p[key]),key
            balance=float(np.max(abs(np.sum(d['history_scaled'],axis=2,dtype=np.longdouble)+d['discards_scaled']-d['ledgers_scaled'])/np.maximum(d['norms_scaled'],1.)))
            hist.append(d['history_scaled'][::n//16].copy());assert balance<1e-8
        ph=producer.prior.paths(1)[0]/f'steps-{n}-reference-128.npz'
        old=Path('native-pressure-return173-work/sweep-1/photons')/ph.name
        assert sha(ph)==sha(old)
        _,_,port=producer.run.packets(ph)
        rows.append(dict(clock=n,all_history_balance=balance,photon_input_byte_equal=True,angular_port=port))
    errors=producer.run.c.relative(*hist)
    assert max(errors)<.02 and max(abs(a-b) for a,b in zip(errors,production['time_comparison']))<1e-14
    cost=sum(read(p)['seconds'] for p in w.glob('*-receipt.json'))+read(w/'probe-localization.json')['seconds']
    assert cost<640
    write(w/'full-material-audit.json',dict(classification='Counterexample candidate',passed=True,rows=rows,
        time_comparison=errors,accepted_prefix_arrays_preserved=True,all_same_photon_inputs_preserved=True,
        action_and_localization_seconds=cost,original_total_cap_seconds=640,
        physical_final_charge_solved=False,full_goal_complete=False,seconds=time.monotonic()-start))
    print(json.dumps(read(w/'full-material-audit.json')),flush=True)


def snapshot():
    p=root/'.phase174-doc-prefixes.json';assert not p.exists()
    write(p,{f'docs/{name}.md':dict(bytes=(root/f'docs/{name}.md').stat().st_size,sha256=sha(root/f'docs/{name}.md')) for name in docs})


def package():
    assert not manifest.exists();a=read(work/'full-material-audit.json');r=read(work/'result.json')
    prod=read(work/'sweep-1/material-flux-tangent/production.json')
    assert a['passed'] and prod['passed'] and r['passed'] and not r['physical_final_charge_solved']
    prefixes=read(root/'.phase174-doc-prefixes.json')
    previous=read(out.parent/'native-pressure-return-manifest.json')
    preserved={p:h for p,h in previous['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    for name in modules:assert sha(root/'verification'/name)==sha(runtime/'verification'/name)
    out.mkdir();copies={};reused={}
    for p in work.rglob('*'):
        if not p.is_file():continue
        rel=p.relative_to(work)
        if rel.parts[:2]==('sweep-1','photons') and p.suffix=='.npz':
            old=out.parent/'native-pressure-return/completed/sweep-1/photons'/p.name
            assert old.is_file() and sha(old)==sha(p)
            reused[rel.as_posix()]=dict(path=old.relative_to(root).as_posix(),sha256=sha(old));continue
        q=out/'completed'/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        assert sha(p)==sha(q);copies[rel.as_posix()]=sha(q)
    shutil.copyfile(__file__,out/'publication-producer.py')
    write(out/'publication.json',dict(copies=copies,reused_published_photon_arrays=reused,
        metadata='Failed/pilot rows retain their original inactive nominal_arithmetic_probe16label. The accepted analytic operator uses local EOS log steps5e-6/1e-5/2e-5; full path metadata names these explicitly. Forward-probe indicators are not uniform error bounds.'))
    summary=dict(classification='Counterexample candidate',same_photon_history_applied_to_free_material=True,
        full_horizon_free_material_completed=True,original_time_gate_passed=True,time_comparison=prod['time_comparison'],
        physical_derivative_maximum=max(row['directional_relative'] for row in prod['rows']+read(work/'sweep-1/material-flux-tangent/pilot.json')['rows']),
        photon_material_energy_H_residual=r['rows'],
        reciprocal_block_accepted=False,GR_return_completed=False,final_charge_conclusion='unadjudicated',
        physical_final_charge_solved=False,full_goal_complete=False,cost=a['action_and_localization_seconds'])
    write(out/'final-result.json',summary)
    tails={
        'model-definition':'동일 수정 광자 이력을 자유 물질에 적용한 전체64/128경로가 원 기준을 통과했다. 대기는 보존량의 직접 원시변수 미분, 방향 minmod와 HLL 대수 미분을 사용한다. 현재 추가 계량이0인 반환에만 적용하며 동일 EOS·접합 면·보존 장부를 유지한다.',
        'observable-targets':'수정된 광자와 전체 자유 물질 해를 연결했지만 상호 반환과 GR 판독은 미완료다. 최종 전하의 기존 부호·크기가 유지되는지는 여전히 미판정이다.',
        'adiabatic-limit':'HLL 대수의 미분 항등식은 선택한 파속 분기에서 증명했다. 이 항등식과 유한 경로 대조는 전체 EOS·미분 오차나 정적 비흡수성의 증명이 아니다.',
        'nonadiabatic-regime':'실제4/8단계에서 통과한 물질 해를 보존하여3.434431ms 전체 경로까지 이어갔다. 같은 광자 입력을 재계산하지 않았으며 다음에는 실제 물질 입력을 광자에 되돌려야 한다.',
        'failure-ledger-dynamic-chi':'원 전역 섭동, 역산 허용오차 강화, 국소 유속 차분의 실패를 모두 보존한다. 마지막 방식은 끝점 검사만 통과하고 실제 중간 궤적에서 실패했다. HLL 대수의 직접 미분을 실제 경로에 적용하여 원 미분·보존·시간 기준을 통과했다. 원 호출수 단위 예산 거절과 물리 하위 단계 단위 재평가도 보존한다.',
        'dynamic-charge-completion':'full_horizon_free_material_completed=true이며 same_photon_history_applied_to_free_material=true다. reciprocal_block_accepted,GR_return_completed,physical_final_charge_solved,full_goal_complete는false다. 연구 가치 기준인 수정된 동일 결합 해의 최종 전하 유지 여부는 미판정이다.'}
    for name,text in tails.items():
        p=root/f'docs/{name}.md';v=prefixes[p.relative_to(root).as_posix()]
        assert p.stat().st_size==v['bytes'] and sha(p)==v['sha256']
        with p.open('ab') as f:f.write(('\n\n## 단계174 — 동일 수정 해의 실제 자유 물질 반환\n\n분류: Counterexample candidate. '+text+' [단계174 보고](../notes/REQUEST174_SAME_SOLUTION_MATERIAL_PLAN_KO.md).\n').encode())
    note=root/'notes/REQUEST174_SAME_SOLUTION_MATERIAL_PLAN_KO.md'
    report='\n\n## 실제 실행과 판정\n\n분류: Counterexample candidate. **최종 전하의 유지 여부는 미판정이다.** 동일 단계173광자 해를 자유 물질에 적용하고, 수락한4/8단계에서 전체64/128경로까지 이어갔다. '
    report+=f"전 구간 B·S·E·H 시간 차이는 각각{[100*x for x in prod['time_comparison']]}%이고 원2%기준을 통과했다. 실제 물질 하위 단계는{[row['substeps'] for row in prod['rows']]}개다. 미분 차이 최대는{100*summary['physical_derivative_maximum']}%로 원0.2%기준 안이다. 모든 저장 절점의 보존과 실제 각도 출구를 재검사했다.\n\n"
    report+='분류: Counterexample candidate. 원 전역 보존량 차분은 실제 예비 경로에서 실패했다. 대기 역산의 실제 허용오차는 이미2e-14였고,1e-14로 강화해도 실패했다. 보존량의 직접 원시변수 미분과 면별 크기 조절은 끝점에서 통과했으나 실제 중간 시각에서 운동량0.927%로 실패했다. 마지막으로 방향 minmod와 HLL의 대수 미분을 적용하고 EOS의 국소 p,u,gamma 편미분만 반·기본·두 배 간격으로 비교했다. 실패를 없애거나 기준을 바꾸지 않았다.\n\n'
    report+='분류: Proven. 선택한 파속 분기에서 HLL 몫 미분 항등식을 SymPy로 확인했다. 정확한 동률은 방향 min/max, minmod 모서리는 한쪽 방향 미분으로 처리한다. 이는 전체 EOS 균일 미분 보장이 아니다.\n\n'
    report+=f"운영 집계: 모든 실행 action과 국소화 검사 합계{a['action_and_localization_seconds']:.3f}초로640초 안이다. 호출당 비용 예측537.55초는 거절했다. 전역 차분에서 원 연산자를 두 번, 대수 미분에서는 한 번 부르므로 실제 SSP 하위 단계를 공통 단위로 다시 계산했다. 이전 전체 경로의 늦은 CFL 단계 수와 큰 쪽의 실측 단가,2배 여유를 유지한379.29초가 원450초 상한에 들어가 생산을 수락했다. 광자 전체 경로는 재실행하지 않았다.\n\n"
    report+=f"분류: Counterexample candidate. 전체 물질 해와 현재 광자/열/H 해의 E/H 잔차는{[x['photon_material_energy_H_residual'] for x in r['rows']]}다. 이것은 상호 결합 수락이 아니다. 현재 광자 경로가 사용한 이전 자유 물질 입력과 새 실제 물질 입력의 차이를 다음에 광자 방정식으로 반환해야 한다. GR·에너지·경계 장부와 최종 전하 판독은 남아 있다.\n\n"
    report+='분류: Conjectural. 이 수리는 저장된 배경 위의 추가 계량0인 물질 응답에만 적용된다. 전체 비선형 별, 균일 EOS/미분·공간·보간·외부 경계 오차 및 정적 비교/관측 폐쇄는 미완료다. 기존 선택 전하의 부호를 이 해로 이전하지 않는다.\n'
    with note.open('ab') as f:f.write(report.encode())
    before=sha(master);write(out/'preservation.json',dict(previous_phase173_nondoc=preserved,document_prefixes=prefixes,previous_master_sha256=before))
    files=[root/'verification'/name for name in modules]+[root/'verification/return_native_pressure_matter.py',note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    summary.update(previous_master_sha256=before,preserved_document_prefixes=prefixes,sha256={p.relative_to(root).as_posix():sha(p) for p in files})
    write(manifest,summary);data=read(master);data['sha256'].update(summary['sha256']);data['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    data['native_pressure_same_solution_material_return']={k:v for k,v in summary.items() if k!='sha256'};write(master,data)


def check(mode):
    m=read(manifest);data=read(master)
    for p,h in m['sha256'].items():assert sha(root/p)==h and data['sha256'][p]==h,p
    for p,h in read(out/'preservation.json')['previous_phase173_nondoc'].items():assert sha(root/p)==h,p
    for p,v in m['preserved_document_prefixes'].items():assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest()==v['sha256'],p
    p='docs/dynamic-charge-completion.md';assert (root/p).read_bytes().splitlines()[:20]==git('show','d5744f766719f0a29c7acb913c7984bdfd3624e2:'+p).splitlines()[:20]
    assert data['sha256'][manifest.relative_to(root).as_posix()]==sha(manifest)
    paths=list(m['sha256'])+[manifest.relative_to(root).as_posix(),'paper/revision-manifest.json']
    Path('C:/Users/lpaiu/AppData/Local/Temp/native-pressure-matter-174-paths').write_bytes(b'\0'.join(p.encode() for p in paths)+b'\0')
    if mode in ['staged','head']:
        if mode=='staged':
            tracked=set(filter(None,git('ls-tree','-r','--name-only','-z','HEAD','--',*paths).decode().split('\0')))
            changed={p for p in paths if p not in tracked or hashlib.sha256(git('show','HEAD:'+p)).hexdigest()!=sha(root/p)}
            assert set(filter(None,git('diff','--cached','--name-only','-z').decode().split('\0')))==changed
        for p in paths:assert hashlib.sha256(git('show',(':' if mode=='staged' else 'HEAD:')+p)).hexdigest()==sha(root/p),p
    print(json.dumps(dict(bound_files=len(m['sha256']),published_paths=len(paths),prefixes_preserved=6,full_material_completed=True,final_charge_conclusion='unadjudicated')))


if __name__=='__main__':
    action=sys.argv[1]
    if action in ['audit','snapshot']:globals()[action]()
    else:
        if action=='package':package()
        check(action)
