from pathlib import Path
import hashlib,json,shutil,subprocess,sys,time

root=Path('E:/lab/self-mass-unobservability')
runtime=Path('//wsl.localhost/Ubuntu-22.04/home/lpaiu/work/native-retained-tail-runtime')
work=runtime/'native-pressure-return173-work'
out=root/'outputs/direct-eos-gr33/native-pressure-return'
manifest=out.parent/'native-pressure-return-manifest.json';master=root/'paper/revision-manifest.json'
read=lambda p:json.loads(p.read_text(encoding='utf-8-sig'))

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for chunk in iter(lambda:f.read(1024**2),b''):h.update(chunk)
    return h.hexdigest()

def write(p,v):p.write_bytes((json.dumps(v,ensure_ascii=False,indent=2)+'\n').encode())
def git(*args):return subprocess.check_output(['rtk','proxy','git',*args],cwd=root)


def audit():
    import numpy as np
    from fractions import Fraction as F
    import continue_native_pressure_response as run
    start=time.monotonic();w=run.OUT;r=read(w/'result.json')
    assert r['passed'] and read(w/'status.json')['state']=='completed'
    assert read(w/'run-receipt.json')['error'] is None
    assert not (w/'full-history-audit.json').exists()
    for p,h in read(w/'plan.json')['bindings'].items():assert sha(p)==h,p
    values=[];ports=[];files=[Path(__file__),Path(run.__file__)]
    for n in [64,128]:
        p=w/f'sweep-1/photons/steps-{n}-reference-128.npz';initial=w/f'sweep-1/photons/pilot-{n}.npz'
        error,v,count=run.verify_packet(p,initial);values.append(v);ports.append(error);files.extend([p,initial])
        assert count==r['actual_steps'][str(n)]
    rows=[]
    for k in range(2,17):
        a,b=[v[:k+1] for v in values]
        errors=np.max(np.sum(abs(a-b),axis=2),axis=0)/np.maximum(np.max(np.sum(abs(b),axis=2),axis=0),1e-290)
        record=w/f'comparison-{k:02d}.json';saved=read(record);files.append(record)
        assert saved['passed'] and max(errors)<.02
        assert np.max(abs(errors-np.asarray(saved['time_comparison'])))<1e-14
        rows.append(dict(interval=k,time_comparison=errors.astype(float).tolist()))
    # Unequal accepted substeps still carry the exact Radau moments.
    total=[]
    for power in range(3):
        value=F(0)
        for a,b in [(F(0),F(1,4)),(F(1,4),F(1))]:
            h=b-a;value+=h*(F(3,4)*(a+h/3)**power+F(1,4)*b**power)
        assert value==F(1,power+1);total.append(str(value))
    write(w/'full-history-audit.json',dict(classification='Counterexample candidate',passed=True,
        final_history_reproduces_every_interval_verdict=True,rows=rows,angular_port_relative=ports,
        accepted_prefix_arrays_preserved=True,final_charge_conclusion='unadjudicated',
        bindings={str(p):sha(p) for p in files},seconds=time.monotonic()-start))
    write(w/'symbolic.json',dict(classification='Proven',passed=True,moments=total,
        scope='Exact degree0..2 Radau moments on unequal intervals; no full coupled error certificate.'))
    print(json.dumps(dict(full_history_audit=True,angular_port_relative=ports,seconds=time.monotonic()-start)))


def package():
    assert not manifest.exists();final=read(work/'result.json');audit_result=read(work/'full-history-audit.json')
    assert final['passed'] and audit_result['passed'] and not final['full_goal_complete']
    assert read(work/'status.json')['state']=='completed'
    prefixes=read(root/'.phase173-doc-prefixes.json');prior=read(out.parent/'native-pressure-front-manifest.json')
    preserved={p:h for p,h in prior['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    module=root/'verification/continue_native_pressure_response.py';assert sha(module)==sha(runtime/'verification'/module.name)
    receipts={p.name:read(p) for p in work.glob('*-receipt.json')}
    assert set(receipts)=={'prepare-receipt.json','run-receipt.json'}
    assert all(v['error'] is None and v['source_sha256']==sha(module) for v in receipts.values())
    assert receipts['run-receipt.json']['seconds']<4800
    dest=out/'completed';dest.mkdir(parents=True);copies={};omitted=[]
    for p in work.rglob('*'):
        if not p.is_file():continue
        rel=p.relative_to(work)
        keep=p.suffix in ['.json','.py'] or (p.suffix=='.npz' and (rel.parts[0]=='sweep-0' or p.stem in ['pilot-64','pilot-128','steps-64-reference-128','steps-128-reference-128']))
        if not keep:
            omitted.append(dict(path=rel.as_posix(),bytes=p.stat().st_size));continue
        q=dest/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        assert sha(q)==sha(p);copies[rel.as_posix()]=sha(q)
    shutil.copyfile(module,dest/'final-producer.py');copies['final-producer.py']=sha(module)
    shutil.copyfile(__file__,out/'publication-producer.py')
    cost=dict(action_seconds=sum(v['seconds'] for v in receipts.values()),CPU_seconds=sum(v['CPU_seconds'] for v in receipts.values()),
        maximum_RSS_bytes=max(v['peak_RSS_bytes'] for v in receipts.values()),run_cap_seconds=4800,
        scope='Prepare and continuation bodies. Read-only audit, imports, publication and Git are separate; no free-material or GR steps in this phase.')
    final.update(cost=cost,full_history_audit=True,final_charge_conclusion='unadjudicated on the corrected coupled photon/free-material/GR solution')
    write(out/'final-result.json',final)
    write(out/'publication.json',dict(copies=copies,all_final_states_and_accepted_prefixes_copied=True,
        intermediate_runtime_files_not_duplicated=omitted,
        rationale='Final arrays reproduce every interval verdict. Intermediate/checkpoint NPZ files stay in the runtime; do not duplicate their growing histories in Git.'))
    before=sha(master);data=read(master)
    tails={
      'model-definition':'분류: Counterexample candidate. 동일 수정 원천과 국소 이분 규칙으로 광자·열·수소 결합 해를 지정 전체 기간까지 이어갔다. 실제117/227단계의 출구 시각·가중치와 충돌·일 이력을 저장했다.',
      'observable-targets':'분류: Counterexample candidate. 전체 기간 광자·열·수소 시간 대조와 출구 이력이 원 기준을 통과했다. 자유 물질·GR 반환이 아직 없으므로 최종 전하의 기존 결론 유지 여부는 미판정이다.',
      'adiabatic-limit':'분류: Proven. 서로 다른 하위 구간에서도 실제 Radau 가중치는2차 이하 다항식의 적분 모멘트를 보존한다. 이 항등식과 수치 경로 통과는 정적 no-go나 전체 물리 오차의 증명이 아니다.',
      'nonadiabatic-regime':'분류: Counterexample candidate. 수락한6/12단계 예비 해에서 이어간3.434431ms 경로가 매 비교 절점의 원 여섯 시간 기준을 통과했다. 같은 해의 물질·GR 상호 반환을 완료해야 최종 전하를 읽을 수 있다.',
      'failure-ledger-dynamic-chi':'분류: Counterexample candidate. 원 실패를 유지한 채 입력 도달 구간의 한 번 이분으로 통과한 수정 해가 전체 기간에서도 시간 기준을 유지했다. 이를 서로 다른 해의 통과 성분 조합이나 진단값의 질량 가산으로 대체하지 않았다.',
      'dynamic-charge-completion':'분류: Counterexample candidate. full_horizon_photon_thermal_H_completed=true이다. free_material_response_completed,physical_final_charge_solved,full_goal_complete는false이며 최종 전하 결론은 미판정이다.'}
    for name,line in tails.items():
        p=root/f'docs/{name}.md';v=prefixes[p.relative_to(root).as_posix()]
        assert p.stat().st_size==v['bytes'] and sha(p)==v['sha256']
        with p.open('ab') as f:f.write(('\n\n## 단계173 — 수정 결합 해의 전체 기간 계속\n\n'+line+' [단계173 보고](../notes/REQUEST173_COUPLED_CONTINUATION_PLAN_KO.md).\n').encode())
    note=root/'notes/REQUEST173_COUPLED_CONTINUATION_PLAN_KO.md'
    percentages=[100*v for v in final['time_comparison']]
    report='\n\n## 실행 결과\n\n분류: Counterexample candidate. 최종 전하 결론은 미판정이다. 지정 전체 기간의 광자·열·수소 진화와15개 후속 비교 구간을 완료했다. '
    report+=f"실제 단계 수는{final['actual_steps']}이고, 마지막 누적 시간 대조의 광자E·물질E·H·충격량·광자압력·물질압력 차이는 각각{percentages}%로 원2% 기준 안이다. 초기 수락 해의 배열을 보존했고, 최종 이력으로 모든 중간 판정을 재현했다.\n\n"
    report+=f"운영 집계: 준비·진화 action {cost['action_seconds']:.3f}초, CPU {cost['CPU_seconds']:.3f}초, 최대 RSS {cost['maximum_RSS_bytes']}바이트. 전체 경로 재실행 없이326개 남은 실제 단계를 진행했다.\n\n"
    report+='분류: Conjectural. 다음 필수 단계는 같은 충돌 전달의 자유 물질 반환, 상호 결합 잔차와 GR·경계 폐쇄다. 이 결과 자체로 최종 전하 안정성이나 전체 EOS·비선형·관측 폐쇄를 주장하지 않는다.\n'
    with note.open('ab') as f:f.write(report.encode())
    write(out/'preservation.json',dict(previous_phase172_nondoc=preserved,document_prefixes=prefixes,previous_master_sha256=before))
    files=[module,note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    summary=dict(final,previous_master_sha256=before,preserved_document_prefixes=prefixes,
        sha256={p.relative_to(root).as_posix():sha(p) for p in files})
    write(manifest,summary);data['sha256'].update(summary['sha256']);data['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    data['native_pressure_full_horizon_continuation']={k:v for k,v in summary.items() if k!='sha256'};write(master,data)


def check(mode):
    m=read(manifest);data=read(master)
    for p,h in m['sha256'].items():assert sha(root/p)==h and data['sha256'][p]==h,p
    for p,h in read(out/'preservation.json')['previous_phase172_nondoc'].items():assert sha(root/p)==h,p
    for p,v in m['preserved_document_prefixes'].items():assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest()==v['sha256'],p
    completion='docs/dynamic-charge-completion.md'
    assert (root/completion).read_bytes().splitlines()[:20]==git('show','d5744f766719f0a29c7acb913c7984bdfd3624e2:'+completion).splitlines()[:20]
    assert data['sha256'][manifest.relative_to(root).as_posix()]==sha(manifest)
    paths=list(m['sha256'])+[manifest.relative_to(root).as_posix(),'paper/revision-manifest.json']
    Path('C:/Users/lpaiu/AppData/Local/Temp/native-pressure-return-173-paths').write_bytes(b'\0'.join(p.encode() for p in paths)+b'\0')
    if mode in ['staged','head']:
        if mode=='staged':
            tracked=set(filter(None,git('ls-tree','-r','--name-only','-z','HEAD','--',*paths).decode().split('\0')))
            changed={p for p in paths if p not in tracked or hashlib.sha256(git('show','HEAD:'+p)).hexdigest()!=sha(root/p)}
            assert set(filter(None,git('diff','--cached','--name-only','-z').decode().split('\0')))==changed
        for p in paths:assert hashlib.sha256(git('show',(':' if mode=='staged' else 'HEAD:')+p)).hexdigest()==sha(root/p),p
    print(json.dumps(dict(bound_files=len(m['sha256']),published_paths=len(paths),preserved_prefixes=6,full_photon_horizon_passed=True,final_charge_conclusion='unadjudicated')))


if __name__=='__main__':
    if sys.argv[1]=='audit':audit()
    else:
        if sys.argv[1]=='package':package()
        check(sys.argv[1])
