"""Archive the one-return photon-boundary pair (second long-double attempt), its 575-time input, its attempts and its same-solution charge.

Adapted from the main checkout's .phase265-publish.py (written for the H-gate run, which never completed).
This session commits in its own worktree, so root is that worktree; the main checkout is read only (its
publication chain for the helpers and its Windows-side diagnostic logs). Scripts come from the runtime,
where they ran.
"""
from pathlib import Path
from decimal import Decimal
import datetime,hashlib,importlib.util,json,os,subprocess,sys
legacy=Path('E:/lab/self-mass-unobservability')
here=Path(__file__).resolve();os.chdir(legacy)  # the chain loads its predecessors by relative path
spec=importlib.util.spec_from_file_location('prior',legacy/'.phase264-publish.py');b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root=Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
runtime=b.runtime;master=root/'paper/revision-manifest.json';read,write,sha,copy=b.read,b.write,b.sha,b.copy
previous_manifest=root/'outputs/direct-eos-gr33/native-primitive-photon-full-manifest.json'
out=root/'outputs/direct-eos-gr33/native-photon-boundary'
manifest=out.parent/'native-photon-boundary-manifest.json'
note=root/'notes/REQUEST265_ONE_RETURN_PHOTON_BOUNDARY_KO.md'
MODULES=['propagate_geometric_boundary_clock.py','apply_photon_geometric_boundary.py','apply_photon_geometric_boundary_newton.py',
         'apply_photon_geometric_boundary_hgate.py','apply_photon_geometric_boundary_ld.py','solve_long_double_fgmres.py','read_photon_boundary_exterior.py']
SCRIPTS=['.phase265-geometric-controller.py','.phase265-geometric-resume.py','.phase265-geometric-queue.py','.phase265-pack.py',
         '.phase265-followthrough.py','.phase265-chain.py','.phase265-newton-followthrough.py','.phase265-newton-chain.py',
         '.phase265-newtonb-followthrough.py','.phase265-newtonb-chain.py','.phase265-hgate-followthrough.py','.phase265-hgate-chain.py',
         '.phase265-floor-inspect.py','.phase265-hstate-inspect.py','.phase265-hgate-dryrun.py','.phase265-complete-test.py',
         '.phase265-fail-inspect2.py','.phase265-hgate-coarse-inspect.py','.phase265-late-cost-inspect.py','.phase265-linear-fail-inspect.py',
         '.phase265-hgate-status.py','.phase265-ld-probe.py','.phase265-ld-chain.py','.phase265-ld-followthrough.py','.phase265-ld-status.py',
         '.phase265-ld-refine-test-first.py','.phase265-ld-refine-test.py','.phase265-ld-fine231-test.py',
         '.phase265-ld2-chain.py','.phase265-ld2-followthrough.py','.phase251-reader-launch.py']
DIAGNOSTICS=[legacy/n for n in ['.phase265-floor-inspect.log','.phase265-hstate-inspect.log','.phase265-hgate-dryrun.stdout.log','.phase265-complete-test.stdout.log',
             '.phase265-ld-chain.stderr.log','.phase265-hgate-chain.stderr.log']]
DIAGNOSTICS+=[runtime/n for n in ['.phase265-ld-probe.json','.phase265-ld-refine-test-first.json','.phase265-ld-refine-test.json','.phase265-ld-fine231-test.json']]
ATTEMPTS=dict(first=('native-photon-boundary265-work','native-photon-boundary265-chain.json'),
              name_error=('native-photon-boundary-newton265-work','native-photon-boundary-newton265-chain.json'),
              newton12=('native-photon-boundary-newton265b-work','native-photon-boundary-newton265b-chain.json'),
              hgate=('native-photon-boundary-hgate265-work','native-photon-boundary-hgate265-chain.json'),
              long_double_log=('native-photon-boundary-ld265-work','native-photon-boundary-ld265-chain.json'))
SKIP=('sweep-0','initialization','__pycache__')


def files_of(folder,suffixes,skip=SKIP):
    return [p for p in sorted(folder.rglob('*')) if p.is_file() and p.suffix in suffixes and not any(k in p.parts for k in skip)]


def package():
    geom=runtime/'native-geometric-clock265-work';actual=runtime/'native-photon-boundary-ld2-265-work'
    charge=runtime/'native-photon-boundary-ld2-charge265-work/full';exterior=runtime/'native-photon-boundary-ld2-exterior265-work'
    chain=read(runtime/'native-photon-boundary-ld2-265-chain.json');status=read(actual/'pipeline-status.json')
    g=read(geom/'result.json');check=read(actual/'boundary-check.json');physical=read(actual/'result.json');plan=read(actual/'plan.json')
    compact=read(charge/'charge/result.json');ext=read(exterior/'full/result.json');audit=read(exterior/'full/audit.json')
    complete=read(exterior/'complete/result.json')
    assert chain['state']=='completed' and status['state']=='completed' and all(x['returncode']==0 for x in status['completed'])
    assert g['passed'] and all(g['knot_values_exact']) and g['pilot_repeat_exact'] and check['passed']
    assert physical['passed'] and not physical['scientific_gates_changed'] and physical['photon_geometric_boundary_applied']
    assert physical['internal_acceptance_changed_for_H'] and physical['internal_material_gate']=={'Etilde':1e-13,'H':2e-13,'B':1e-13,'S':1e-13}
    assert [v['actual_completed_steps'] for v in physical['rows']]==[119,231]
    assert compact['charge_comparison_admitted'] and compact['same_actual_returned_solution_read']
    assert ext['passed'] and audit['passed'] and complete['passed']
    events=physical['long_double_fallback_events'];assert all(e['accepted'] for v in events.values() for e in v)
    assert physical['hgate_prefix_reused'] and physical['first_long_double_attempt_preserved'] and physical['restart_interval']==15
    assert all(read(actual/f'restart-regression-{n}.json')['passed'] for n in [64,128])
    tests=dict(refine=read(runtime/'.phase265-ld-refine-test.json'),refine_first=read(runtime/'.phase265-ld-refine-test-first.json'),
               fine231=read(runtime/'.phase265-ld-fine231-test.json'),probe=read(runtime/'.phase265-ld-probe.json'))
    assert tests['refine']['passed'] and tests['refine']['json_written'] and tests['refine_first']['error']=="KeyError('linear')"
    assert 'finished' in tests['fine231']
    plans=[geom/'plan.json',geom/'execution-plan.json',actual/'plan.json',actual/'controller-start.json',charge/'endpoint/plan.json',
           charge/'dense/plan.json',charge/'charge/plan.json',exterior/'full/plan.json',exterior/'full/audit.json',exterior/'complete/plan.json']
    bindings={}
    for p_ in plans:
        for p,h in read(p_)['bindings'].items():
            assert sha(runtime/p)==h,(str(p_),p);assert p not in bindings or bindings[p]==h,p;bindings[p]=h
    for p,h in chain.get('bindings',{}).items():assert sha(runtime/p)==h,p
    for folder in [actual,charge,exterior]:
        for p in folder.rglob('*receipt.json'):
            if 'initialization' not in p.parts:assert read(p).get('error') is None,p
    for p in geom.glob('queue-q*-receipt.json'):assert read(p)['error'] is None,p
    assert read(charge/'charge/polynomial-audit.json')['passed']
    dense={str(n):read(charge/f'dense/source-{n}-check.json') for n in [64,128]}
    assert all(x['passed'] and x['dense_stage_max']<1e-12 for x in dense.values())
    for m in MODULES:assert sha(root/'verification'/m)==sha(runtime/'verification'/m),m
    for s in SCRIPTS:assert (runtime/s).is_file(),s
    for p in DIAGNOSTICS:assert p.is_file(),p
    preserved={p:h for p,h in read(previous_manifest)['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    assert not out.exists();out.mkdir()
    # 575-time input: summaries, stop and queue records, merged values; per-item values are packed once.
    for p in files_of(geom,['.json','.npz']):
        if (p.name.startswith('queue-') and p.name[6:7].isdigit()) or p.name.startswith('resume-'):continue
        copy(p,out/'geometric'/p.relative_to(geom))
    # Long-double arrays are read only in WSL: .phase265-pack.py packed the per-item queue values and
    # verified them bitwise against the merged boundary; here that pack is checked and copied as bytes.
    pack=read(geom/'queue-items-check.json');assert pack['passed'] and pack['items']==len(list(geom.glob('queue-[0-9]*.json')))
    assert pack['queue_items_npz_sha256']==sha(geom/'queue-items.npz') and pack['boundary_sha256']==sha(geom/'boundary-575.npz')
    copy(runtime/'native-geometric-clock265-status.json',out/'geometric/status.json')
    for p in sorted(runtime.glob('.phase265-geometric-*.log')):
        if p.stat().st_size:copy(p,out/'geometric/logs'/p.name.lstrip('.'))
    # Preserved attempts: every record; their large state arrays stay in the runtime and are hashed.
    attempt_sha={}
    for label,(name,chain_name) in ATTEMPTS.items():
        folder=runtime/name
        for p in files_of(folder,['.json','.log','.py']):copy(p,out/'attempts'/label/p.relative_to(folder))
        for p in files_of(folder,['.npz'],skip=SKIP+('gr',)):attempt_sha[f'{name}/{p.relative_to(folder).as_posix()}']=sha(p)
        copy(runtime/chain_name,out/'attempts'/label/'chain-status.json')
    for pattern in ['.phase265-chain-*.log','.phase265-newton*-chain-*.log','.phase265-hgate-chain-*.log','.phase265-ld-chain-*.log']:
        for p in sorted(runtime.glob(pattern)):
            if p.stat().st_size:copy(p,out/'attempts/logs'/p.name.lstrip('.'))
    # Final coupled pair: every record, the accepted states and captures. Metric arrays are reproducible
    # from the259metric and boundary-575; their hashes are recorded instead of copying them.
    for p in files_of(actual,['.json','.log','.py']):copy(p,out/'actual265'/p.relative_to(actual))
    for n in [64,128]:
        for name in [f'recovered-{n}.npz',f'sweep-1/photons/return-{n}.npz']:copy(actual/name,out/'actual265'/name)
        for name in [f'endpoint/gr/endpoint-{n}.npz',f'dense/gr/source-{n}.npz']:copy(charge/name,out/'charge265'/name)
    metric_sha={f'metric-{n}-g{q}.npz':sha(actual/f'metric/metric-{n}-g{q}.npz') for n,q in [(64,8),(128,4),(128,8)]}
    # The reused prefix copies and the last saved pairs are reproducible from the H-gate run; hashed only.
    actual_sha={p.relative_to(actual).as_posix():sha(p) for p in sorted(actual.glob('*.npz')) if not p.name.startswith('recovered-')}
    for p in files_of(charge,['.json','.log','.py']):copy(p,out/'charge265'/p.relative_to(charge))
    copy(charge/'dense/geometry.npz',out/'charge265/dense/geometry.npz')
    for n,q in [(64,8),(128,4),(128,8)]:copy(charge/f'charge/gr/fields-{n}-g{q}.npz',out/f'charge265/charge/gr/fields-{n}-g{q}.npz')
    for part in ['full','complete']:
        for p in files_of(exterior/part,['.json','.npz','.log']):copy(p,out/'exterior265'/part/p.relative_to(exterior/part))
    for p in files_of(exterior,['.json','.log'],skip=SKIP+('full','complete')):copy(p,out/'exterior265'/p.relative_to(exterior))
    copy(runtime/'native-photon-boundary-ld2-265-chain.json',out/'chain-status.json')
    for p in sorted(runtime.glob('.phase265-ld2-chain-*.log')):
        if p.stat().st_size:copy(p,out/'logs'/p.name.lstrip('.'))
    for s in SCRIPTS:copy(runtime/s,out/'scripts'/s.lstrip('.'))
    for p in DIAGNOSTICS:copy(p,out/'diagnostics'/p.name.lstrip('.'))
    copy(here,out/'publication-producer.py')
    # Values reported below.
    now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    c=next(v for v in compact['components'] if v['clock']==128)['values']['endpoint_compact_with_metric']
    old_c=read(root/'outputs/direct-eos-gr33/native-corrected-charge/charge259/charge/result.json')
    old_c=next(v for v in old_c['components'] if v['clock']==128)['values']['endpoint_compact_with_metric']
    e=next(v for v in ext['components'] if v['clock']==128);f=next(v for v in complete['components'] if v['clock']==128)
    old_e=read(root/'outputs/direct-eos-gr33/native-corrected-charge/result.json')['conditional_fine_exterior']
    low_change=float((Decimal(c['low'])-Decimal(old_c['low']))/Decimal(old_c['low']))
    high_same=Decimal(c['high'])==Decimal(old_c['high'])
    total_change=float((Decimal(f['total'])-Decimal(old_e['total']))/Decimal(old_e['total']))
    xc=complete['cross_checks'];rows={(r['component'],r['clock'],r['angular'],r['radial']):r for r in complete['rows']}
    hr=rows['high',128,8,8]
    fro={(r['component'],r['clock'],r['angular'],r['radial']):r for r in ext['rows']}['high',128,8,8]
    # G/c^4 from the frozen reader's own definitions: kappa=homogeneous/M, epsilon=G/c^4*arrived/M.
    fac=Decimal(fro['epsilon'])*Decimal(fro['homogeneous_mass_cm'])/(Decimal(fro['kappa'])*Decimal(fro['arrived_energy_erg']))
    last264=read(runtime/'native-primitive-photon263-work/fine/result.json')['rows'][-1]
    low_measure=float(Decimal(last264['photon_J_source_cm'])/fac/Decimal(last264['physical_launch_energy_increment_erg']))
    receipts={a:read(actual/f'{a}-receipt.json')['seconds'] for a in ['coarse','fine','audit']}
    hrows=physical['H_gate_rows'];term=read(geom/'termination.json');sched=g['scheduling']['items_by_process']
    tries=plan['attempt_H_material_relative'];hg=plan['hgate_attempt'];ld1=plan['first_long_double_attempt'];probe=plan['long_double_fallback']['probe']
    def event_text(tag,e):
        r=e['final'];gas=r['material_relative']
        return (f"{tag} {e['actual_step']}단계 {e['corrections']}회 보정({e['seconds']:.0f}초, 벡터 {r['relative']:.2e}, 물리 모멘트 최대 {max(r['moments']):.2e}, "
                f"H {max(gas[1],gas[5]):.2e}, 다른 물질 성분 최대 {max(x for k,x in enumerate(gas) if k%4!=1):.2e})")
    ev_text='; '.join(event_text(t,e) for t,v in [('coarse',events['coarse']),('fine',events['fine'])] for e in v) or '없음'
    rt=tests['refine'];ft=tests['fine231'];fev=ft.get('long_double_events',[])
    if ft.get('stage_passed'):
        f231=(f"재생 stage는 원 기준을 통과했다({ft['stage_seconds']:.0f}초). 첫 풀이 호출은 기록된 호출 {ft['first_call_matches_recorded']}번과 "
              f"비트 일치했고, 그 뒤 {ft.get('compared_calls')}개 중 {ft.get('bitwise_calls_from_match')}개가 비트 일치했다. "
              f"재생의 long-double 사건은 "+('; '.join(event_text('fine',dict(e,final=e['rows'][-1],actual_step=231)) for e in fev) or '없음')+"이다"
              f"(재생 모델은 새로 만들어 그 안의 단계 번호가 1로 기록된다).")
    else:
        f231=f"재생 stage는 통과하지 못했다({ft.get('stage_error') or ft.get('error')})."
    verdict='NEGATIVE_CHARGE_RETAINED_WITH_ONE_RETURN_PHOTON_BOUNDARY' if f['sign_negative'] else 'CHARGE_SIGN_CHANGED_WITH_ONE_RETURN_PHOTON_BOUNDARY'
    final=dict(classification='Counterexample candidate',passed=True,verdict=verdict,
        actual_same_joint_solution_GR_computed=True,actual_steps=[119,231],full_declared_period=True,
        same_horizon_seconds=physical['same_horizon_seconds'],photon_geometric_boundary_applied=True,
        geometric_boundary_times=g['output_times'],geometric_boundary_knots_bitwise=True,
        maximum_geometric_lapse_over_applied_outer_lapse=g['maximum_geometric_lapse_over_applied_outer_lapse'],
        internal_material_gate=physical['internal_material_gate'],internal_acceptance_changed_for_H=True,
        H_gate_user_approval=plan['gate_change']['user_approval'],H_gate_rows=hrows,attempt_H_material_relative=tries,
        attempts={k:v[0] for k,v in ATTEMPTS.items()},attempt_npz_sha256=attempt_sha,scientific_gates_changed=False,
        hgate_attempt=hg,first_long_double_attempt=ld1,long_double_fallback=plan['long_double_fallback'],long_double_fallback_events=events,
        long_double_tests=dict(refine_first=tests['refine_first'],refine=rt,fine231_replay={k:v for k,v in ft.items() if k!='stage_audit'}),
        reused_prefix_intervals=15,evolved_interval=16,reused_actual_steps=plan['reused_actual_steps'],actual_npz_sha256=actual_sha,
        paired_time_relative=max(physical['time_relative']),dense_stage_max={k:v['dense_stage_max'] for k,v in dense.items()},
        compact_high=c['high'],compact_low=c['low'],compact_high_unchanged_from_phase260=bool(high_same),
        compact_low_relative_change_from_phase260=low_change,frozen_exterior_total=e['total'],complete_total=f['total'],
        complete_total_relative_change_from_phase260=total_change,final_sign_negative=f['sign_negative'],
        cross_checks=dict(xc,phase264_low_photon_mass_over_launch_energy=low_measure),metric_sha256=metric_sha,
        coarse_seconds=receipts['coarse'],fine_seconds=receipts['fine'],scheduling=dict(items_by_process=sched,termination=term),
        exterior_scalar_operator_variation_complete=False,self_GR_return_closed=False,
        physical_final_charge_solved=False,final_charge_conclusion='unadjudicated',full_goal_complete=False,snapshot_KST=now)
    sign='음수로 유지됐다' if f['sign_negative'] else '양수로 바뀌었다'
    first_h=', '.join(f'{v:.4e}' for v in tries['first']);retry_h=f"{min(tries['newton12']):.5e}~{max(tries['newton12']):.5e}"
    h64,h128=hrows['64'],hrows['128']
    note.write_text(f'''# 1회 반환의 외부 광자 경계를 적용한 같은 결합 해의 전하

분류: Counterexample candidate. **1회 반환·고정 외부 범위에서 외부 광자 경계를 닫은 같은 결합 해의 전하는 {sign}.** 단계261이 계산했지만 적용 계량에 빠져 있던 입사 계량의 배경 광자 기하 lapse를 575개 적용 시각 전부에서 계산해 반환 계량 outer lapse에 넣었다. 원 119/231 결합 쌍을 새 계량으로 t=0부터 다시 진화했다. 마지막 구간은 H 기준 실행이 수락한 15/16 canonical 구간을 정확한 재시작·native 속도 재생 검사 뒤 재사용하고, 두 double 선형 풀이가 모두 실패할 때만 쓰는 long-double 대체 풀이를 더해 다시 진화했다. 그 해 자체의 compact·고정 외부·질량 정규화에 배경 광자의 발사 에너지 증분까지 더해 읽었다. 미세 시계의 최종 합은 {float(f['total']):.12e}로, 단계260의 {float(old_e['total']):.12e}에서 상대 {total_change:.3e} 바뀌었다. 사용자 승인에 따라 단계257 내부 물질 기준의 중성 H 성분만 1e−13에서 2e−13으로 바꿨고, 원래의 과학적 수락 기준은 모두 그대로다. 전체 물리 전하는 아래 남은 항목 때문에 아직 판정 불가다. 게시={now}.

분류: Proven. 적용 계량의 outer lapse는 photon+scalar−residual·K(r0)이고, photon 항의 광선 핵은 K(r0)−(1+μ²)K(r_p)다. 따라서 기존 경계는 −(물질 질량)K(r0)−G/c⁴Σe(1+μ²)K(r_p)+scalar이며, 단계261 특수해 −G/c⁴ΣE₀K(r_p)[(1+μ²)(δlnH−δν−δλ−δr/(rb))+2μδμ]는 같은 식을 배경 패킷에 대해 변분한 것이다. 새 경계는 기존 경계에 이 특수해를 더한 값이다. 질량 잔차에 광자 질량을 다시 더하면 이중 계산이므로 더하지 않았다. Lapse.run은 경계에서 안쪽으로 적분하므로 δν·그 면값·δlnN·δln(속도)에 같은 시간 함수가 더해지고 δλ·δν′·δu는 그대로다. 저장 배열에서 이 관계가 비트 단위로 성립함을 확인했다.

분류: Counterexample candidate. 575시각 경계는 16개 기존 매듭값을 비트 단위로 재현했고, 대표 측정값의 재실행도 비트 일치했다. 패킷 에너지 항등식 최대 {g['maximum_energy_identity_relative']:.3e}, 각 불변량 최대 {g['maximum_angular_invariant_relative']:.3e}로 원 기준 안이다. 새 lapse 항은 적용 outer lapse의 최대 {g['maximum_geometric_lapse_over_applied_outer_lapse']:.4e}배이고 종단 광자 질량은 {g['terminal_photon_geometric_mass_cm']:.6e} cm다. 계산 중 네 작업자가 같은 입력에서 다른 작업자보다 약 3배 느려 등록 상한을 넘을 것으로 보였다. 이들만 멈추려고 원 컨트롤러를 종료하자 WSL이 같은 세션의 14개 작업자를 모두 종료했고 수신증은 남지 않았다. 이미 저장된 {sum(term['kept_items'].values())}개 항목은 보존했고 나머지 {term['remaining_items']}개는 공유 대기열로 끝냈다. 같은 모듈 함수·패킷·기준을 쓰고 스케줄만 바꿨으며, 중단 기록과 첫 재개 시도의 실패를 보존했다.

분류: Counterexample candidate. 새 계량은 시간 대조 최대 {max(check['metric_controls']['time'].values()):.4e}, 구적 대조 최대 {max(check['metric_controls']['quadrature'].values()):.4e}로 원 2%·0.2% 기준을 통과했고, stage driver가 새 파일을 읽는 것도 확인했다. 첫 재진화는 coarse 3단계(t=5.366e−5 s)에서 단계257 성분별 내부 물질 기준(1e−13)에 걸렸다. 벡터 잔차 3.8e−16과 물리 모멘트는 통과했지만 H가 세 Newton 제안에서 {first_h}였다. 같은 단계의 단계259 값은 7.46e−14였고, 259 전체의 H 최댓값은 coarse 7.46e−14·fine 8.79e−14였다. Newton 제안을 12회로 늘린 첫 재시도는 반복 횟수를 이름으로 넣은 구현 오류(NameError)로 물리 단계 전에 멈췄다. 숫자로 고친 재시도에서 벡터 잔차는 5.4e−17까지 줄었지만 H는 {retry_h}에 머물렀다. 결함의 약 77%가 H가 있는 안쪽 6개 셀에 몰려 있었고, 이 단계에 들어가는 low H 상태는 259와 L1 0.36% 달랐다. 그래서 수렴 부족이 아니라 stiff H 속도의 long-double 반올림 바닥으로 판단했다(Conjectural). 사용자가 H 내부 기준만 2e−13으로 바꾸는 안을 승인했다. 이 값은 선형 제안과 비선형 수락에 모두 적용했다. Etilde·B·S 1e−13과 원 벡터·물리 모멘트·구성·면 수지·port·짝 시간·dense 1e−12 기준은 바꾸지 않았다. 세 실패 시도의 기록과 거부 단계, 진단 스크립트·출력은 보존했다.

분류: Counterexample candidate. H 기준 실행({hg['directory']})은 coarse {hg['coarse_accepted_steps']}/119단계와 fine {hg['fine_accepted_steps']}/231단계까지 수락했다. 이전에 막힌 3단계를 포함한 모든 단계가 새 H 기준과 나머지 원 기준을 통과했다. coarse 119단계에서는 짧은 오른쪽 전처리 제안 12회와 원 double 풀이의 보정 12회 뒤에도 long-double 참 잔차가 벡터 1.25e−12, 광자 수 모멘트 5.6e−9로 원 기준(1e−14·1e−13)에 못 미쳤다. double 보정마다 참 잔차 비가 {hg['double_refinement_ratio_range'][0]}~{hg['double_refinement_ratio_range'][1]}였으므로, GMRES 안의 double 연산자가 long-double 수락 연산자를 더는 충분히 대표하지 못한다고 판단했다. H 기준과는 무관한 실패다. 저장된 마지막 수락 쌍으로 그 선형계를 재구성해 기록된 첫 GMRES 호출을 비트 단위로 재현했다. 같은 double 전처리기에 Arnoldi·Hessenberg 최소제곱·갱신을 long double로 하는 flexible GMRES를 쓰자 {probe['corrections']}회 보정으로 벡터 {probe['final']['relative']:.2e}·모멘트 최대 {max(probe['final']['moments']):.2e}에 도달했고, 그 단계는 원 비선형 기준도 통과했다. 이 풀이는 두 double 풀이가 모두 실패할 때만 쓰는 세 번째 대체 경로다. double 풀이가 수락하는 선형계는 이전과 똑같이 계산되므로, H 기준 실행의 15/16 canonical 구간을 재사용했다.

분류: Counterexample candidate. 첫 long-double 시도({ld1['directory']})는 재시작 검사를 통과한 뒤 coarse {ld1['coarse_accepted_steps']}단계와 fine {ld1['fine_accepted_steps']}단계까지 수락했다. fine 231단계에서 두 double 풀이가 실패했고(마지막 double 반복 벡터 1.48e−13, 모멘트 최대 4.64e−9), long-double 보정이 실행됐다. 그러나 그 로그를 쓰는 중 numpy bool이 JSON에 들어가 구현 오류로 멈췄다. 이 bool은 벡터 잔차가 이미 1e−14 아래이고 물리 모멘트가 1e−13 이상인 보정 행에서만 생기므로, 적어도 한 번의 보정이 실행됐지만 결과는 기록되지 않았다. 수락 판정을 bool()로 고쳤다. 고친 보정을 재구성한 coarse 119단계 선형계에 돌린 scratch 시험의 첫 실행은 보정이 반환된 뒤 시험 하니스의 마지막 기록 줄에서 KeyError로 멈췄다. 그 한 줄만 고쳐 다시 실행하자 {rt['corrections']}회 보정으로 벡터 {rt['final']['relative']:.2e}에 도달했고 로그도 기록됐다. 계산이 결정적이므로 두 번째 시도는 fine 231단계에서 첫 시도의 미기록 보정을 그대로 반복한다. 그래서 첫 시도의 저장된 230단계 쌍으로 그 stage 전체를 scratch에서 production 경로로 재생하는 진단을 병행했다. {f231} 세 시험·진단과 두 실패 시도의 기록은 보존했다.

분류: Counterexample candidate. 원 119/231단계를 모두 완료했고 짝 시간 대조 최대는 {max(physical['time_relative']):.6e}다. 두 번째 long-double 시도의 대체 풀이 사건은 {ev_text}이며, 모두 원 선형 기준을 통과했다. H가 1e−13 이상인 채 수락된 단계는 coarse {h64['steps_with_H_above_1e_13']}개(최대 {h64['maximum_H']:.4e})와 fine {h128['steps_with_H_above_1e_13']}개(최대 {h128['maximum_H']:.4e})다. 다른 물질 성분의 최댓값은 coarse {h64['maximum_other_material']:.3e}, fine {h128['maximum_other_material']:.3e}다. Newton 제안 수 분포는 coarse {h64['proposals_used']}, fine {h128['proposals_used']}다. 같은 해의 dense 원천 단계 잔차는 coarse {dense['64']['dense_stage_max']:.4e}, fine {dense['128']['dense_stage_max']:.4e}로 원 1e−12를 통과했다. 마지막 구간만 다시 진화한 두 경로는 병행했으며 coarse {receipts['coarse']:.0f}초, fine {receipts['fine']:.0f}초가 걸렸다. 격자·시계·기간·허용오차는 바꾸지 않았다.

분류: Counterexample candidate. compact high는 primary 이력의 값이므로 {'단계260과 같다' if high_same else '단계260과 다르다'}({float(c['high']):.12e}). 새 경계로 다시 진화한 low는 {float(c['low']):.12e}로, 단계260 대비 상대 {low_change:.3e} 바뀌었다. 고정 외부·질량 정규화만 적용한 합은 {float(e['total']):.12e}다. 배경 광자의 발사 에너지 증분을 같은 고정 외부 연산자의 방출로 더하면 high의 외부 스칼라가 {hr['geometric_exterior']:.4e}, 방출 에너지가 {hr['geometric_exited_energy_erg']:.6e} erg, 도착 에너지가 {hr['geometric_arrived_energy_erg']:.6e} erg만큼 더해진다. 그 결과 합은 {float(f['total']):.12e}가 되며, 고정 외부만의 합 대비 상대 {f['relative_change_from_frozen_only']:.3e}다. 시간·각도·반경 대조와 두 표현 격자의 대조는 원 기준을 통과했다.

분류: Conjectural. 질량 장부는 다음과 같이 두었다. 배경 광자의 발사 에너지 증분 중 외부로 나간 몫은 외부 광자 질량으로 κ에, 무한대에 도착한 몫은 ε에 넣었다. 이것이 완전한 ADM·Bondi 질량 장부라는 증명은 없다. 발사 에너지 모형은 단계261 전 전파의 발사 에너지와 16개 매듭에서 {xc['high_launch_energy_vs_phase261_knots']:.3e}, 575개 시각에서 {xc['high_launch_energy_vs_575_time_propagation']:.3e} 이내로 같다. 전 전파의 광자 질량과 발사 에너지의 차이는 high에서 {xc['full_photon_mass_minus_launch_energy_relative']:.3e}로, 측도·편향·계량 일 항이 작다. low의 발사 에너지는 단계264(이전 계량) 값과 {xc['low_launch_energy_vs_phase264_old_metric']:.3e} 다르다. 단계264 종단에서 low 광자 질량은 발사 에너지의 {low_measure:.2f}배로 측도 항이 크지만, 이 항은 low κ의 약 1e−10 수준이라 합에는 영향이 없다. 이 차이는 모형 한계로 보고하며 흡수하지 않았다.

분류: Conjectural. 남은 항목은 외부 스칼라 연산자의 배경 변분, 기하 광자 스칼라 원천의 편향·측도 항, 두 번째 반환과 자기GR, 반환 원천의 별도 궤적 시계 대조(단계264 전 전파), EOS·공간·균일 오차, 비선형, 정적 EFT·관측 비교다. 이번 결과는 1회 반환·고정 외부 범위에서 외부 광자 경계의 누락을 고친 같은 해의 판정이며, 사용자 승인에 따른 H 내부 기준 조정에 조건부다. 마지막 단계들의 선형계는 원 기준을 long-double 대체 풀이로 만족했다. 이는 수락 기준의 변경이 아니지만, double 연산자가 수락 연산자를 대표하지 못하는 구간이 생긴다는 수치적 경계로 기록한다. 이를 전역적·관측적 확정으로 확대하지 않는다.
''',encoding='utf-8',newline='\n')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        p=root/f'docs/{name}.md';prefixes[p.relative_to(root).as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
        text=(f'\n\n## 단계265 — 1회 반환 외부 광자 경계를 적용한 같은 해의 전하\n\n'
            f'분류: Counterexample candidate. 입사 계량의 배경 광자 기하 lapse를 575개 적용 시각에서 계산해 반환 계량 outer lapse에 넣고, 원 119/231 결합 쌍을 t=0부터 다시 진화했다. '
            f'16개 매듭은 단계261과 비트 일치했고 계량·결합·판독의 원 기준을 통과했다. 첫 재진화는 단계257 내부 물질 기준의 H 성분(1e−13)이 반올림 바닥(1.00e−13)에 걸려 멈췄고, 사용자 승인으로 H만 2e−13으로 바꿔 다시 진화했다. '
            f'그 실행은 coarse 마지막 단계에서 double 선형 풀이가 long-double 수락 연산자를 대표하지 못해 멈췄다. 두 double 풀이가 모두 실패할 때만 쓰는 long-double flexible GMRES를 더하고, 15/16 구간을 정확한 재시작 검사 뒤 재사용해 마지막 구간을 다시 진화했다. '
            f'같은 해의 compact·고정 외부·질량 정규화에 배경 광자의 발사 에너지까지 더한 미세 시계 합은 {float(f["total"]):.6e}로 {sign}(단계260 대비 상대 {total_change:.2e}). '
            f'다섯 실패 시도(H 바닥, NameError, Newton 12회, double 선형 풀이, 대체 풀이 로그의 JSON 오류)와 작업자 스케줄 실패(컨트롤러 종료 시 WSL 세션 작업자 전체 종료)는 보존했다. 외부 스칼라 연산자 변분·자기GR·EOS/공간·관측 폐쇄가 남아 전체 물리 전하는 미판정이다. [근거](../notes/{note.name}).\n')
        with p.open('ab') as h:h.write(text.encode())
    write(out/'result.json',final)
    write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes,verified_runtime_bindings=bindings))
    files=[root/'verification'/m for m in MODULES]+[note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes,sha256={p.relative_to(root).as_posix():sha(p) for p in files})
    write(manifest,final);m=read(master);m['sha256'].update(final['sha256'])
    m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_photon_boundary']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


def git(*args):return subprocess.check_output(['git',*args],cwd=root)


def check(mode):
    """The phase162/186 check, bound to this worktree: manifest, master, previous files, document prefixes; staged or HEAD bytes."""
    m=read(manifest);a=read(master)
    for p,h in m['sha256'].items():assert sha(root/p)==h and a['sha256'][p]==h,p
    for p,h in read(out/'publication.json')['previous_nondoc'].items():assert sha(root/p)==h,p
    for p,v in m['document_prefixes'].items():assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest()==v['sha256'],p
    assert a['sha256'][manifest.relative_to(root).as_posix()]==sha(manifest)
    paths=list(m['sha256'])+[manifest.relative_to(root).as_posix(),'paper/revision-manifest.json']
    (here.parent/'phase265-ld2-paths.txt').write_bytes(b'\0'.join(p.encode() for p in paths)+b'\0')
    if mode in ['staged','head']:
        if mode=='staged':
            staged=git('-c','core.quotepath=off','diff','--cached','--name-only','-z').decode().split('\0')
            assert set(filter(None,staged))==set(paths),sorted(set(filter(None,staged))^set(paths))[:10]
        ref=':' if mode=='staged' else 'HEAD:'
        for p in paths:assert hashlib.sha256(git('show',ref+p)).hexdigest()==sha(root/p),p
    print(json.dumps(dict(mode=mode,bound_files=len(m['sha256']),paths=len(paths),prefixes_preserved=len(m['document_prefixes']),
        final_sign_negative=m['final_sign_negative'],final_charge_conclusion=m['final_charge_conclusion'])))


if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1])
