"""Archive the phase-266 depth decomposition of the accepted primary charge (diagnostic, no new physical step).

root is this session's worktree; the main checkout is read only (publication-chain helpers).
"""
from pathlib import Path
import datetime,hashlib,importlib.util,json,os,subprocess,sys
legacy=Path('E:/lab/self-mass-unobservability')
here=Path(__file__).resolve();os.chdir(legacy)
spec=importlib.util.spec_from_file_location('prior',legacy/'.phase264-publish.py');b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root=Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
runtime=b.runtime;master=root/'paper/revision-manifest.json';read,write,sha,copy=b.read,b.write,b.sha,b.copy
previous_manifest=root/'outputs/direct-eos-gr33/native-photon-boundary-manifest.json'
out=root/'outputs/direct-eos-gr33/native-depth-decomposition'
manifest=out.parent/'native-depth-decomposition-manifest.json'
note=root/'notes/REQUEST266_PRIMARY_DEPTH_DECOMPOSITION_KO.md'
FILES=['.phase266-depth-probe.py','.phase266-depth-probe-ineffective.py','.phase266-depth-b128.json','.phase266-depth-b64.json',
       '.phase266-depth-b128-ineffective.json','.phase266-depth-b64-ineffective.json','.phase266-depth-v128.json','.phase266-depth-v128.stderr.log',
       '.phase266-depth-v2.json','.phase266-charge-history.py','.phase266-charge-history.json','.phase266-bands-summary.py','.phase266-bands-summary.txt']


def package():
    for f in FILES:assert (runtime/f).is_file(),f
    fine,coarse=read(runtime/'.phase266-depth-b128.json'),read(runtime/'.phase266-depth-b64.json')
    assert fine.get('error') is None and coarse.get('error') is None
    rows={r['band']:r for r in fine['rows']};rows64={r['band']:r for r in coarse['rows']}
    total=rows['all']['charge'];assert total==fine['saved_free_endpoint']
    parts=[k for k in rows if k!='all'];s=sum(rows[k]['charge'] for k in parts);s64=sum(rows64[k]['charge'] for k in parts)
    closure=(s-total)/total;closure64=(s64-rows64['all']['charge'])/rows64['all']['charge'];assert abs(closure)<1e-12 and abs(closure64)<1e-12
    v2=read(runtime/'.phase266-depth-v2.json');vr={r['band']:r['charge'] for r in v2['rows']}
    assert vr['copy']==vr['all']==v2['saved_free_endpoint']
    bad=[read(runtime/f'.phase266-depth-b{n}-ineffective.json') for n in [128,64]]
    assert all(len({r['charge'] for r in x['rows']})==1 for x in bad)
    hist=read(runtime/'.phase266-charge-history.json')
    preserved={p:h for p,h in read(previous_manifest)['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    assert not out.exists();out.mkdir()
    for f in FILES:copy(runtime/f,out/f.lstrip('.'))
    copy(here,out/'publication-producer.py')
    now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    share=lambda k:rows[k]['charge']/total
    neg=sum(rows[k]['charge'] for k in parts if rows[k]['charge']<0);pos=sum(rows[k]['charge'] for k in parts if rows[k]['charge']>0)
    rel=lambda k:(rows64[k]['charge']-rows[k]['charge'])/rows[k]['charge']
    deep=max(abs(rows[f'cells:{i}-{i}']['charge']) for i in range(8))
    shallow=sum(abs(rows[k]['charge']) for k in ['cells:16-16','cells:17-17','cells:18-18','cells:19-146','cells:147-274','cells:275-402','cells:403-530','boundary'])
    worst_dominant=max(abs(rel(f'cells:{i}-{i}')) for i in range(9,16))
    final=dict(classification='Counterexample candidate',passed=True,diagnostic_only=True,physical_steps=0,
        verdict='DOMINANT_ERROR_IS_PRIMARY_INTERIOR_RADIAL_RESOLUTION',
        endpoint_free_charge_fine=total,endpoint_free_charge_coarse=rows64['all']['charge'],band_closure_relative=[closure,closure64],
        shares={k:share(k) for k in parts},charges_fine={k:rows[k]['charge'] for k in parts},charges_coarse={k:rows64[k]['charge'] for k in parts},
        coarse_vs_fine_relative={k:rel(k) for k in parts},negative_sum=neg,positive_sum=pos,positive_over_negative=pos/abs(neg),
        deepest_invisible_cell_max=deep,shallow_cells_atmosphere_boundary_abs=shallow,worst_dominant_cell_time_relative=worst_dominant,
        visible_depth_at_T_km=hist['visible_depth_at_T_km'],e_folding_depth_km_range=hist['e_folding_depth_km_range'],
        interior_cell_width_km=hist['interior_cell_widths_km'][0],setup_plus_propagate_seconds_per_band=sum(r['setup_seconds']+r['propagate_seconds'] for r in fine['rows'])/len(fine['rows']),
        first_probe_harness_error_preserved=True,refined_primary_rerun_done=False,final_charge_conclusion='unadjudicated',full_goal_complete=False,snapshot_KST=now)
    note.write_text(f'''# 수락된 primary 전하의 깊이 분해와 지배 오차의 식별

분류: Counterexample candidate. **단계265에서 유지된 조건부 음의 전하는 68.75km 폭의 내부 셀 몇 개가 결정한다.** 그 high 성분(단계248 primary GR장의 끝점 전하)을 원천의 깊이 대역으로 선형 분해했다. 대역 합은 전체와 같았다(fine 상대 {closure:.1e}, coarse {closure64:.1e}). 가리지 않은 원천은 저장된 끝점 자유장 값 {total:.16e}를 비트 단위로 재현했다. 새 물리 단계는 없으며 저장 원천과 production 전파기만 사용했다. 게시={now}.

분류: Counterexample candidate. fine 시계의 대역별 비중은 다음과 같다.

| 셀 | 깊이 [km] | 전하 기여 | 비중 |
|---|---|---|---|
''' + ''.join(f"| {k.replace('cells:','')} | {rows[k]['depth_km'][0]:.1f}–{rows[k]['depth_km'][1]:.1f} | {rows[k]['charge']:.4e} | {share(k)*100:.2f}% |\n" for k in parts if k!='boundary')
      + f'''| 경계 누적 에너지 | – | {rows['boundary']['charge']:.4e} | {share('boundary')*100:.2e}% |

분류: Counterexample candidate. 전하는 셀 11(깊이 276–345km, {share('cells:11-11')*100:.1f}%)과 셀 12(207–276km, {share('cells:12-12')*100:.1f}%)가 대부분을 차지하고, 셀 10과 13이 그 다음이다. 끝점에서 가장 깊이 보이는 셀 9·8은 **양**의 기여로 합 {(share('cells:9-9')+share('cells:8-8'))*100:.1f}%다. 551km보다 깊은 셀의 기여는 셀당 최대 {deep:.1e}로 전하의 {deep/abs(total):.0e}배이며, 끝점의 가시 깊이 cT/2={hist['visible_depth_at_T_km']:.1f}km와 맞는다. 얇은 셀·대기·경계 항의 합은 {shallow/abs(total)*100:.2f}%다. 전하를 지배하는 셀 9–15의 64/128 시계 차이는 최대 {worst_dominant:.1e}이고, 양의 전선 셀 8에서는 {rel('cells:8-8'):.1e}다. 그러므로 이 분해는 시간 격자에 둔감하다. 저장 이력의 |q|는 후반부에 로그 기울기 {hist['last_half_log_slope_per_s'][0]:.0f}–{hist['last_half_log_slope_per_s'][1]:.0f}/s로 커지며, 가시 깊이로 환산한 e-배 길이는 {hist['e_folding_depth_km_range'][0]:.0f}–{hist['e_folding_depth_km_range'][1]:.0f}km로 셀 폭 {hist['interior_cell_widths_km'][0]:.2f}km보다 작다.

분류: Conjectural. 이 분해와 기존 기록으로 남은 오차원을 비교했다. 부호를 바꿀 수 없는 항목은 다음과 같다. 1회 GR 반환(low/high 3.2e−22), 외부 광자 경계(3.9e−8)와 그 편향·측도 항(약 2e−13), 반환 원천 시계 대조(3.9e−8 항에만 들어감), 판독 구적(1.5e−12 이하), 시간 격자(끝점 4.2e−5)다. 외부 스칼라 연산자의 배경 변분은 단계159의 조건부 상계 0.714% 안이다. 반면 원천의 반응은 약 3개 셀에 걸쳐 있다. 이 셀들은 구동 계산에서 한 번도 세분되지 않았다. 단계119는 마지막 내부 셀만 4개로 나눴고, 그 대조 대상도 무구동 직접 전하였다. 따라서 현재 지배 오차는 구동된 primary 이력의 내부 반경 해상도로 판단한다. 크기 오차는 수십 퍼센트일 수 있다(추정이며 계산하지 않았다).

분류: Conjectural. 음의 기여 합은 {neg:.4e}, 양의 기여 합은 {pos:.4e}로 비는 {pos/abs(neg):.3f}다. 부호가 바뀌려면 세분 뒤 이 균형이 약 {abs(neg)/pos:.0f}배 바뀌어야 한다. 이는 부호가 견고할 가능성을 보이는 수치 논거일 뿐 세분 계산이 아니다.

분류: Conjectural. 결정적 시험은 원천 셀 8–15(깊이 5–551km)를 세분해 무구동·구동 primary 이력을 같은 최종 방정식으로 다시 진화하고, 같은 전파기로 끝점 전하를 읽는 것이다. 판정 규칙은 다음과 같다. 세분한 해의 끝점 전하가 음이면 조건부 음의 전하 결론이 이 지배 오차를 고쳐도 유지된 것이다. 크기 변화는 2% 기준에 대한 해상도 수렴으로 따로 보고한다. 1회 GR 반환과 외부 광자 항은 위 크기로 부호에 영향이 없으므로 재계산 여부를 별도로 정한다. 이 경로는 단계181–244의 생산 체인을 새 격자에서 다시 만드는 일이다. 생산 경로 조사와 대표 구간 실측으로 예산을 정한 뒤에 시작한다.

분류: Counterexample candidate. 첫 분해 스크립트는 setup에 넘긴 배열을 가렸다. 그러나 production setup은 단계244 원천 파일을 직접 읽으므로 모든 대역이 전체 값을 돌려줬다. 이 하니스 오류 기록(depth-*-ineffective)을 보존했다. 가린 원천 파일을 임시 입력 디렉터리에 두어 production setup을 그대로 쓰도록 고쳤다. 대역당 비용은 약 {final['setup_plus_propagate_seconds_per_band']:.0f}초, 최대 RSS는 약 5.4GB였다. 퍼텐셜 반환 항은 저장 끝점의 2.5e−11이라 분해에서 뺐다.
''',encoding='utf-8',newline='\n')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        p=root/f'docs/{name}.md';prefixes[p.relative_to(root).as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
        text=(f'\n\n## 단계266 — 수락된 primary 전하의 깊이 분해와 지배 오차\n\n'
            f'분류: Counterexample candidate. 단계265에서 유지된 조건부 음의 전하의 high 성분을 원천 깊이 대역으로 선형 분해했다. 대역 합은 전체와 같았고 저장값을 비트 단위로 재현했다. '
            f'셀 11(276–345km)이 {share("cells:11-11")*100:.0f}%, 셀 12가 {share("cells:12-12")*100:.0f}%를 차지하며, 가장 깊이 보이는 셀 9·8은 양의 기여 {(share("cells:9-9")+share("cells:8-8"))*100:.0f}%다. '
            f'대역별 시간 격자 차이는 지배 셀에서 1e−4 이하다. 반응이 약 3개의 미세분 68.75km 셀에 몰려 있으므로, 현재 지배 오차는 구동된 primary 이력의 내부 반경 해상도로 판단한다(Conjectural). '
            f'음·양 기여의 비는 약 {abs(neg)/pos:.0f}:1이다. 세분 격자에서 primary 이력을 다시 진화하는 결정적 시험이 다음이며, 최종 전하는 미판정이다. [근거](../notes/{note.name}).\n')
        with p.open('ab') as h:h.write(text.encode())
    write(out/'result.json',final)
    write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes))
    files=[note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes,sha256={p.relative_to(root).as_posix():sha(p) for p in files})
    write(manifest,final);m=read(master);m['sha256'].update(final['sha256'])
    m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_depth_decomposition']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


def git(*args):return subprocess.check_output(['git',*args],cwd=root)


def check(mode):
    m=read(manifest);a=read(master)
    for p,h in m['sha256'].items():assert sha(root/p)==h and a['sha256'][p]==h,p
    for p,h in read(out/'publication.json')['previous_nondoc'].items():assert sha(root/p)==h,p
    for p,v in m['document_prefixes'].items():assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest()==v['sha256'],p
    assert a['sha256'][manifest.relative_to(root).as_posix()]==sha(manifest)
    paths=list(m['sha256'])+[manifest.relative_to(root).as_posix(),'paper/revision-manifest.json']
    (here.parent/'phase266-paths.txt').write_bytes(b'\0'.join(p.encode() for p in paths)+b'\0')
    if mode in ['staged','head']:
        if mode=='staged':
            staged=set(filter(None,git('-c','core.quotepath=off','diff','--cached','--name-only','-z').decode().split('\0')))
            assert staged==set(paths),sorted(staged^set(paths))[:10]
        ref=':' if mode=='staged' else 'HEAD:'
        for p in paths:assert hashlib.sha256(git('show',ref+p)).hexdigest()==sha(root/p),p
    print(json.dumps(dict(mode=mode,bound_files=len(m['sha256']),paths=len(paths),verdict=m['verdict'],final_charge_conclusion=m['final_charge_conclusion'])))


if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1])
