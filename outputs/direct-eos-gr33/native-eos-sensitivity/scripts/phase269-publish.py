"""Archive phase 269: EOS physics sensitivity of the compact charge (PL/MHD vs PL-off native EOS and H level libraries).

Counterexample candidate. Copies the small artifacts (stage-1 same-state comparison, library build receipts and
identity/physics checks, preparation logs including the preserved first attempt, segment rows, fallbacks, registered
last-segment exceptions, existence audits, readouts, source/state comparisons) into
outputs/direct-eos-gr33/native-eos-sensitivity, records SHA-256 of the libraries and large runtime arrays without copying
them, binds the reused phase-267/268 scripts by their published SHA, appends the result to the phase-269 note and the six
dynamic-chi documents, and binds everything in the manifest and paper/revision-manifest.json.
Usage: python phase269-publish.py package|check|staged|head
"""
from pathlib import Path
import datetime, hashlib, json, math, subprocess, sys
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
WSL = Path('//wsl.localhost/Ubuntu-22.04/home/lpaiu/work')
O, E, A1, LIBS = WSL/'native-retained-tail-runtime', WSL/'native-eos269-runtime', WSL/'native-eos269-runtime-attempt1', WSL/'direct-eos-gr33'
S = Path(__file__).resolve().parent
master = root/'paper/revision-manifest.json'; previous_manifest = root/'outputs/direct-eos-gr33/native-quad-refined-primary-manifest.json'
out = root/'outputs/direct-eos-gr33/native-eos-sensitivity'; manifest = out.parent/'native-eos-sensitivity-manifest.json'
note = root/'notes/REQUEST269_EOS_PHYSICS_SENSITIVITY_KO.md'
P267, P268 = root/'outputs/direct-eos-gr33/native-refined-primary/scripts', root/'outputs/direct-eos-gr33/native-quad-refined-primary/scripts'
DOCS = ['model-definition', 'observable-targets', 'adiabatic-limit', 'nonadiabatic-regime', 'failure-ledger-dynamic-chi', 'dynamic-charge-completion']
read = lambda p: json.loads(Path(p).read_text(encoding='utf-8'))
def write(p, v): Path(p).parent.mkdir(parents=True, exist_ok=True); Path(p).write_text(json.dumps(v, indent=1, ensure_ascii=False) + '\n', encoding='utf-8', newline='\n')
def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for block in iter(lambda: f.read(1 << 22), b''): h.update(block)
    return h.hexdigest()
def copy(src, dst):
    data = Path(src).read_bytes(); assert len(data) < 4 << 20, (src, len(data)); Path(dst).parent.mkdir(parents=True, exist_ok=True); Path(dst).write_bytes(data)
    assert hashlib.sha256(Path(dst).read_bytes()).hexdigest() == hashlib.sha256(data).hexdigest()
SCRIPTS = ['phase269-eos.py', 'phase269-build.py', 'run269.sh', 'run269-build.sh', 'phase269-checks.py', 'phase269-bank.py', 'run-prepare269.sh',
           'run-eos64.sh', 'run-readouts269.sh', 'compare269.py', 'seg-table.py', 'state-diff269.py', 'source-diff269.py', 'run-source-diff269.sh',
           'wait-any.sh', 'wait-run.sh', 'eos-option.sh', 'eos-options.sh', 'eos-options2.sh', 'eos-lineage.sh', 'eos-libs.sh', 'eos-build-inputs.sh',
           'layers.sh', 'list269.sh', 'phase269-publish.py']
REUSED = {**{n: P267 for n in ['phase267-template.py', 'phase267-undriven.py', 'phase267-p150-setup.py', 'phase267-p150.py', 'phase267-born.py',
                                'phase267-placeholders.py', 'phase267-readout.py', 'phase267-depth.py', 'phase267-audited.py']},
          **{n: P268 for n in ['phase268-driver.py', 'phase268-thermal.py', 'phase268-initial.py', 'phase268-initial-audit.py', 'launch268.sh']}}
def runtime_files():
    files = {}
    for name in ['.phase269-shared.npz', '.phase269-mhd.npz', '.phase269-shared-v1.npz', '.phase269-mhd-v1.npz', '.phase269-compare.log']:
        files['stage1/' + name.lstrip('.')] = O/name
    for p in sorted(O.glob('.phase269-build-*.log')) + sorted(O.glob('.phase269-s2-*')) + [O/'.phase269-stage2a.log']: files['stage2a/' + p.name.lstrip('.')] = p
    for d in ['native-cold-population-identity', 'native-cold-population-ploff', 'photon-eos-levels-identity', 'photon-eos-levels-ploff']:
        files[f'libraries/{d}-phase269-build.json'] = LIBS/d/'phase269-build.json'
    for name in ['.phase269-prepare.log', '.phase269-bank.log', '.phase269-thermal.log']: files['prepare-attempt1/' + name.lstrip('.')] = A1/name
    for name in ['.phase269-prepare.log', '.phase269-copied-inputs.json', '.phase269-bank.log', '.phase269-thermal.log', '.phase269-initial.log',
                 '.phase269-initial-audit.log', '.phase269-undriven.log', '.phase269-p150-audit.json', '.phase269-p150.log', '.phase269-born.log',
                 '.phase267-placeholders.json', '.phase269-exists-smoke.json', '.phase269-smoke.stdout.log', '.phase269-smoke.stderr.log']:
        files['prepare/' + name.lstrip('.')] = E/name
    W = E/'primary269-eos64-work'
    for pattern in ['seg-*-driver.json', 'seg-*-fallbacks.json', 'seg-*-nonlinear-exceptions.json', 'integer-64.json', 'polish.json']:
        for p in sorted(W.glob(pattern)): files['run/' + p.name] = p
    files['run/status.txt'] = E/'.phase269-eos64-logs/status.txt'
    for p in sorted((E/'.phase269-eos64-logs').glob('exists-seg-*.json')): files['exists/' + p.name] = p
    for n in [61, 62, 63, 64]:
        for p in sorted((E/f'readout269-eos{n}-work').glob('*.json')): files[f'readout-{n}/{p.name}'] = p
        for p in sorted(E.glob(f'.phase269-readout-eos{n}-*')): files[f'readout-{n}/{p.name.lstrip(".")}'] = p
        files[f'comparison/source-diff-{n}.json'] = E/f'.phase269-source-diff-{n}.json'
    files['comparison/state-diff.txt'] = E/'.phase269-state-diff.txt'; files['readouts.log'] = E/'.phase269-readouts.log'; files['compare.json'] = E/'readout269-compare.json'
    return files
LARGE = {**{f'direct-eos-gr33/{p}': LIBS/p for p in ['native-cold-population-ploff/libfree_eos_native_cold_stable.so', 'native-cold-population-ploff/stable/gas.so',
                                                   'native-cold-population-ploff/gas-stable-verbose.so', 'photon-eos-levels-ploff/levels.so',
                                                   'native-cold-population-identity/libfree_eos_native_cold_stable.so', 'photon-eos-levels-identity/levels.so']},
         **{f'native-eos269-runtime/{p}': E/p for p in ['primary269-eos64-work/sweep-1/photons/seg-61.npz', 'primary269-eos64-work/sweep-1/photons/seg-62.npz',
            'primary269-eos64-work/sweep-1/photons/seg-63.npz', 'primary269-eos64-work/sweep-1/photons/seg-64.npz', 'readout269-eos64-work/gr/source-64.npz',
            'readout269-eos64-work/gr/field-source-64.npz', 'readout269-eos64-work/recovered-64.npz',
            'outputs/direct-eos-gr33/def-native-boundary-layer/geometry.npz', 'outputs/direct-eos-gr33/def-native-boundary-layer/bank.npz',
            'outputs/direct-eos-gr33/def-native-conservative-rates/thermal-refined/bank.npz',
            'outputs/direct-eos-gr33/def-native-initial-constraints/finite-volume/balanced-initial-state.npz',
            'outputs/direct-eos-gr33/native-retained-completion/evolution/coupled-128.npz', 'outputs/direct-eos-gr33/native-retained-completion/evolution/source-128.npz',
            'native-incident-drive155-work/fields/born-g8.npz']}}
pct = lambda x: f'{100*x:+.2e}%'


def texts(f, now):
    r, T = f['rows'], f['rows']['64']; d = f['depth_T']; run = f['run']; le, ne = run['linear_exceptions'], run['nonlinear_exceptions']
    comp = f['source_T']['components']; big = comp['metric_stress_erg']
    thermal = {k: comp[k] for k in ['gas_nonrest_energy_erg', 'nonrest_trace_erg', 'photon_energy_erg', 'photon_radial_pressure_erg', 'pressure_volume_erg']}
    tmax = max(v['max_abs_plmhd'] for v in thermal.values()); tdiff = (min(v['max_diff_over_max'] for v in thermal.values()), max(v['max_diff_over_max'] for v in thermal.values()))
    kept = f['endpoint_negative_ploff']
    tbl = ''.join(f"| {lab} | {r[n]['q_plmhd_2x']:.10e} | {r[n]['q_ploff_2x']:.10e} | {r[n]['relative_change_vs_plmhd_2x']:+.2e} | {gate} |\n"
                  for n, lab, gate in [('61', 't₆₁=61/64·T', '모든 원 기준'), ('62', 't₆₂', '마지막 구간 규칙'), ('63', 't₆₃', '마지막 구간 규칙'), ('64', 'T(끝점)', '마지막 구간 규칙')])
    dtab = ''.join(f"| {c} | {d[str(c)][0]:.6e} | {d[str(c)][1]:.6e} | {(d[str(c)][1] - d[str(c)][0])/abs(d[str(c)][0]):+.2e} |\n" for c in range(9, 14))
    section = (f"\n\n## 최종 결과 — PL-off EOS·광학에서의 끝점 전하 ({now[:16].replace('T', ' ')} KST)\n\n"
        f"분류: Counterexample candidate. **끝점 compact 전하는 PL-off EOS·광학에서 {T['q_ploff_2x']:.15e}로 {'음이다' if kept else '음이 아니다'}. "
        f"{'사전 등록 규칙에 따라 조건부 음의 전하 결론은 이 EOS·광학 처리 선택에도 유지된다' if kept else '조건부 음의 전하 결론은 유지되지 않는다'}.** "
        f"같은 2배 격자·같은 방정식의 PL/MHD 해와의 상대 차이는 T에서 {T['relative_change_vs_plmhd_2x']:+.2e}이고, 네 시각 모두 {max(abs(v['relative_change_vs_plmhd_2x']) for v in r.values()):.1e} 이하다. "
        f"이 EOS·광학 선택의 영향은 현재 해상도 수렴 수준(2×→4× −1.17%)보다 {int(math.log10(0.011739559728022599/max(abs(v['relative_change_vs_plmhd_2x']) for v in r.values())))}자릿수 이상 작으므로 지배 오차가 아니다.\n\n"
        "| 시각 | PL/MHD 2배 | PL-off 2배 | 상대 차이 | 수락 |\n|---|---|---|---|---|\n" + tbl +
        f"\n분류: Counterexample candidate. 이유는 전하 원천의 구성에 있다. T 판독 원천에서 가장 큰 성분인 metric stress(최대 {big['max_abs_plmhd']:.2e})는 PL-off에서 {big['max_diff_over_max']:.1e}만 다르다. "
        f"EOS·광학에 민감한 열·광자 성분(기체 비정지 에너지, 비정지 트레이스, 광자 에너지·복사압, 압력·부피 일)은 {tdiff[0]*100:.0f}–{tdiff[1]*100:.0f}% 바뀌지만, 최대 크기가 {tmax:.1e}로 지배 성분의 약 {tmax/big['max_abs_plmhd']:.0e}배다. "
        "구동 해의 광자 모멘트와 충돌률도 수십 % 이상 달랐다(상태 비교 기록). 따라서 이 모형의 compact 전하는 스칼라 힘에 대한 정지질량(바리온) 응답이 결정하며, 열·광자 성분의 비중은 이 비교에서 약 1e−9 이하다. "
        "전하 층의 P/(ρc²)가 약 2.3e−9(셀 11)라는 점과 맞는다.\n\n"
        f"분류: Counterexample candidate. 끝점 깊이 분해(원 셀별 2배 두 칸 합)는 다음과 같다(PL-off 대역 합과 전체의 상대 차 {abs(d['closure_ploff']):.1e}).\n\n"
        "| 원 셀 | PL/MHD | PL-off | 상대 차이 |\n|---|---|---|---|\n" + dtab +
        f"\n분류: Counterexample candidate. 수락 기록: 구간 {len(run['segments'])}개, 실제 {run['segments'][-1]['actual_steps']}단계, {run['total_seconds']:.0f}초, 대체 풀이 {run['total_fallbacks']}회. "
        f"macro 0–61은 모든 원 기준을 지켰다. 마지막 구간 규칙으로 수락한 선형 해는 {le['count']}개(벡터 최대 {le['max_vector']:.1e}, 물리 {le['max_physical']:.1e}, 물질 {le['max_material']:.1e}), 비선형 반복은 {ne['count']}개다. "
        "PL-off 층에서 풀이가 느려져 macro 52 이후 구간이 PL/MHD보다 2–4배 오래 걸렸다.\n\n"
        "분류: Conjectural. 범위와 남은 오차: 이 시험은 내부 셀 0–26의 EOS·광학 처리만 바꿨고, 배경 구조(밀도·온도 분포)와 대기 표는 PL/MHD 그대로다. "
        "EOS와 불투명도가 정역학·복사 평형을 거쳐 배경 밀도 분포를 바꾸는 경로는 시험하지 않았다. 전하가 정지질량 응답으로 정해지므로, 남은 물리 오차 후보는 배경 밀도 구조, 스칼라–물질 결합 모형, 입사 펄스다. "
        "최종 물리 전하(자기GR, 비선형, 정적 EFT·관측 폐쇄)는 미판정이다.\n")
    head = (f"분류: Counterexample candidate. **결과 요약(2026-09-27): PL-off EOS·광학에서도 끝점 전하 {T['q_ploff_2x']:.4e}로 {'음의 부호 유지' if kept else '음의 부호 미유지'}. "
            f"PL/MHD 2배 해와의 상대 차이 {T['relative_change_vs_plmhd_2x']:+.1e}(네 시각 최대 {max(abs(v['relative_change_vs_plmhd_2x']) for v in r.values()):.1e}). "
            "전하는 정지질량 응답이 결정하며 EOS·광학 선택은 지배 오차가 아니다. 상세는 맨 아래 최종 결과 절.**\n\n")
    doc = (f"\n\n## 단계269 — 전하 지배 층의 EOS 물리 민감도\n\n"
        "분류: Counterexample candidate. 계산 줄기의 FreeEOS(option 11: Planck–Larkin + MDH)와 수소 준위 라이브러리를 PL만 끈 변형으로 다시 빌드했다(동일성 빌드는 원본과 비트 단위로 일치). "
        "그 라이브러리로 내부 셀 0–26의 EOS·광학 배열을 다시 만들어 2배 격자 primary를 끝점까지 진화했다. 전하 층의 중성 분율은 최대 3배, 일부 광자 계수는 수백 배 바뀌었다. "
        f"그런데도 끝점 compact 전하는 {T['q_ploff_2x']:.4e}로 PL/MHD 해와 상대 {T['relative_change_vs_plmhd_2x']:+.1e}만 다르다. "
        f"{'사전 등록 규칙에 따라 조건부 음의 전하 결론은 유지된다' if kept else '조건부 음의 전하 결론은 유지되지 않는다'}. 전하 원천은 정지질량(metric stress) 성분이 지배하고, EOS·광학에 민감한 열·광자 성분은 그 1e−7–1e−8 규모다. "
        "배경 밀도 구조의 EOS·불투명도 의존성은 시험하지 않았다. 최종 물리 전하는 미판정이다. [근거](../notes/REQUEST269_EOS_PHYSICS_SENSITIVITY_KO.md).\n")
    return head, section, doc


def package():
    assert not out.exists() and not manifest.exists()
    files = runtime_files()
    for dst, src in files.items(): assert src.is_file(), src
    for name in SCRIPTS: assert (S/name).is_file(), name
    reused = {}
    for name, where in REUSED.items(): assert sha(S/name) == sha(where/name), name; reused[name] = dict(sha256=sha(where/name), published=str((where/name).relative_to(root).as_posix()))
    preserved = {p: h for p, h in read(previous_manifest)['sha256'].items() if not p.startswith('docs/')}
    for p, h in preserved.items(): assert sha(root/p) == h, p
    cmp = read(E/'readout269-compare.json'); rows = cmp['rows']; assert cmp['decided'] and set(rows) == {'61', '62', '63', '64'}
    assert all(r['same_time_as_2x'] for r in rows.values())
    for n in rows:
        assert read(E/f'readout269-eos{n}-work/readout-field.json')['endpoint_compact_charge'] == rows[n]['q_ploff_2x']
        assert read(E/f'readout269-eos{n}-work/source-64-check.json')['passed'] and read(E/f'readout269-eos{n}-work/readout-endpoints.json')['assembly']['passed']
    assert abs(cmp['depth_T']['all'][1] - rows['64']['q_ploff_2x']) <= 1e-15*abs(rows['64']['q_ploff_2x'])
    assert all(s['passed'] for s in cmp['run']['segments']) and len(cmp['run']['segments']) == 18
    identity = [l for l in (O/'.phase269-stage2a.log').read_text().splitlines() if l.startswith('{"identical"')]; assert identity and json.loads(identity[0])['identical']
    large = {name: sha(p) for name, p in LARGE.items()}
    captures = sorted((E/'primary269-eos64-work/captures').glob('captured-64-*.npz')); steps = cmp['run']['segments'][-1]['actual_steps']; assert len(captures) == 2*steps
    large[f'native-eos269-runtime/primary269-eos64-work/captures ({len(captures)}, sha256 of the ordered list of file hashes)'] = hashlib.sha256('\n'.join(sha(p) for p in captures).encode()).hexdigest()
    out.mkdir(parents=True)
    for dst, src in files.items(): copy(src, out/dst)
    for name in SCRIPTS: copy(S/name, out/'scripts'/name)
    now = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    kept = cmp['endpoint_negative_ploff']; small = max(abs(v['relative_change_vs_plmhd_2x']) for v in rows.values())
    verdict = ('CONDITIONAL_NEGATIVE_CHARGE_KEPT_EOS_OPTICS_NOT_DOMINANT' if kept and small < 1e-6 else 'CONDITIONAL_NEGATIVE_CHARGE_KEPT_EOS_OPTICS_SENSITIVE' if kept
               else 'CONDITIONAL_NEGATIVE_CHARGE_NOT_KEPT_UNDER_PL_OFF')
    source_T = read(E/'.phase269-source-diff-64.json')
    final = dict(classification='Counterexample candidate', passed=True, verdict=verdict, rule=cmp['rule'], endpoint_negative_ploff=kept, rows=rows,
        max_abs_relative_change=small, eos_variant=dict(call='free_eos_modern(ifoption=3, ifmodified=11, ifion=-2): EOS1 without radiation pressure, ifcoulomb=5, ifpi=3 (PL+MDH)',
            ploff='ifpi=-3 (MDH only) and PL arrays zeroed before use; H level library qstar_calc PL flag off', identity_rebuild_bitwise=True,
            scope='EOS-derived arrays of interior cells 0-26 (2x grid) from the PL-off libraries; background geometry and atmosphere tables PL/MHD'),
        grid=dict(cells=539, interior=27, refined_cells='8-15 split 2x (phase 267 grid)'),
        run=dict(segments=len(cmp['run']['segments']), actual_steps=steps, seconds=cmp['run']['total_seconds'], fallbacks=cmp['run']['total_fallbacks'],
                 linear_exceptions=cmp['run']['linear_exceptions'], nonlinear_exceptions={k: v for k, v in cmp['run']['nonlinear_exceptions'].items() if k != 'rows'},
                 long_double_corrections=24, first_preparation_attempt='stopped at the thermal step: kept-cell EOS arrays were PL/MHD while the native y0 check requires the running EOS (preserved)'),
        depth_T=cmp['depth_T'], source_components_T=source_T['components'], reused_scripts_sha256=reused, large_arrays_sha256=large,
        not_tested=['background structure (density/temperature profile) EOS and opacity dependence', 'atmosphere tables'],
        final_charge_conclusion=('conditional negative charge kept under the PL-off EOS/optics (2x, 64 clock); ' if kept else 'conditional negative charge not kept under PL-off; ') + 'final physical charge unadjudicated',
        full_goal_complete=False, snapshot_KST=now)
    head, section, doc_text = texts(dict(final, source_T=source_T, run=cmp['run']), now)
    text = note.read_text(encoding='utf-8')
    title, rest = text.split('\n\n', 1); note.write_text(title + '\n\n' + head + rest.rstrip('\n') + section, encoding='utf-8', newline='\n')
    prefixes = {}
    for name in DOCS:
        p = root/f'docs/{name}.md'; prefixes[p.relative_to(root).as_posix()] = dict(bytes=p.stat().st_size, sha256=sha(p))
        with p.open('ab') as h: h.write(doc_text.encode())
    write(out/'result.json', final); write(out/'publication.json', dict(previous_nondoc=preserved, document_prefixes=prefixes))
    bound = [note] + [root/p for p in prefixes] + [p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes, sha256={p.relative_to(root).as_posix(): sha(p) for p in bound})
    write(manifest, final)
    raw = master.read_text(encoding='utf-8'); m = json.loads(raw); indent = 2 if raw.startswith('{\n  "') else 1
    m['sha256'].update(final['sha256']); m['sha256'][manifest.relative_to(root).as_posix()] = sha(manifest)
    m['native_eos_sensitivity'] = {k: v for k, v in final.items() if k != 'sha256'}
    master.write_text(json.dumps(m, indent=indent, ensure_ascii=False) + '\n', encoding='utf-8', newline='\n')


def git(*args): return subprocess.check_output(['git', *args], cwd=root)


def check(mode):
    m = read(manifest); a = read(master)
    for p, h in m['sha256'].items(): assert sha(root/p) == h and a['sha256'][p] == h, p
    for p, h in read(out/'publication.json')['previous_nondoc'].items(): assert sha(root/p) == h, p
    for p, v in m['document_prefixes'].items(): assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest() == v['sha256'], p
    for name, v in m['reused_scripts_sha256'].items(): assert sha(root/v['published']) == v['sha256'], name
    assert a['sha256'][manifest.relative_to(root).as_posix()] == sha(manifest)
    paths = list(m['sha256']) + [manifest.relative_to(root).as_posix(), 'paper/revision-manifest.json']
    (S/'phase269-paths.txt').write_bytes(b'\0'.join(p.encode() for p in paths) + b'\0')
    if mode in ['staged', 'head']:
        if mode == 'staged':
            staged = set(filter(None, git('-c', 'core.quotepath=off', 'diff', '--cached', '--name-only', '-z').decode().split('\0')))
            assert staged == set(paths), sorted(staged ^ set(paths))[:10]
        ref = ':' if mode == 'staged' else 'HEAD:'
        for p in paths: assert hashlib.sha256(git('show', ref + p)).hexdigest() == sha(root/p), p
    print(json.dumps(dict(mode=mode, bound_files=len(m['sha256']), paths=len(paths), verdict=m['verdict'], endpoint=m['rows']['64']['q_ploff_2x'],
                          relative_change={n: r['relative_change_vs_plmhd_2x'] for n, r in m['rows'].items()}, final_charge_conclusion=m['final_charge_conclusion'])))


if __name__ == '__main__':
    if sys.argv[1] == 'package': package()
    check(sys.argv[1] if sys.argv[1] != 'package' else 'check')
