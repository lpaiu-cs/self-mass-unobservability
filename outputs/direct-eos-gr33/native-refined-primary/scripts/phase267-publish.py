"""Archive phase 267: the refined-grid primary (interior cells 8-15 split 2x) evolved from t=0 with the final equation.

Counterexample candidate. Copies the small artifacts (drivers, audits, validations, segment rows, fallbacks, approved
exceptions, readouts, depth bands) into outputs/direct-eos-gr33/native-refined-primary, records SHA-256 of the large
runtime arrays without copying them, appends the result to the phase-267 note and the six dynamic-chi documents, and
binds everything in the manifest and paper/revision-manifest.json. Usage: python phase267-publish.py package|check|staged|head
"""
from pathlib import Path
import datetime, hashlib, json, subprocess, sys
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
G = Path('//wsl.localhost/Ubuntu-22.04/home/lpaiu/work/native-retained-tail-runtime')
R = Path('//wsl.localhost/Ubuntu-22.04/home/lpaiu/work/native-refined267-runtime')
S = Path(__file__).resolve().parent
master = root/'paper/revision-manifest.json'; previous_manifest = root/'outputs/direct-eos-gr33/native-depth-decomposition-manifest.json'
out = root/'outputs/direct-eos-gr33/native-refined-primary'; manifest = out.parent/'native-refined-primary-manifest.json'
note = root/'notes/REQUEST267_REFINED_PRIMARY_KO.md'
DOCS = ['model-definition', 'observable-targets', 'adiabatic-limit', 'nonadiabatic-regime', 'failure-ledger-dynamic-chi', 'dynamic-charge-completion']
read = lambda p: json.loads(Path(p).read_text(encoding='utf-8'))
def write(p, v): Path(p).parent.mkdir(parents=True, exist_ok=True); Path(p).write_text(json.dumps(v, indent=1, ensure_ascii=False) + '\n', encoding='utf-8', newline='\n')
def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for block in iter(lambda: f.read(1 << 22), b''): h.update(block)
    return h.hexdigest()
def copy(src, dst):
    data = Path(src).read_bytes(); Path(dst).parent.mkdir(parents=True, exist_ok=True); Path(dst).write_bytes(data)
    assert hashlib.sha256(Path(dst).read_bytes()).hexdigest() == hashlib.sha256(data).hexdigest()
SCRIPTS = ['phase267-audit.py', 'phase267-bank.py', 'phase267-thermal.py', 'phase267-initial.py', 'phase267-initial-audit.py', 'phase267-template.py',
           'phase267-undriven.py', 'phase267-p150-setup.py', 'phase267-p150.py', 'phase267-trace.py', 'phase267-identity.py', 'identity-compare.py',
           'phase267-born.py', 'phase267-placeholders.py', 'phase267-driver-evolve-v1.py', 'phase267-driver-v2.py', 'phase267-driver-v3.py',
           'phase267-driver-v3b.py', 'phase267-driver-v3c-probe.py', 'phase267-driver-v4.py', 'phase267-driver.py', 'phase267-readout.py',
           'phase267-depth.py', 'phase267-audited.py', 'validate-compare.py', 'capture-compare.py', 'exists-compare.py', 'exists-segments.py',
           't61-bands.py', 'T-bands.py', 'inverse-diag.py', 'run-summary.py', 'run-refined64.sh', 'run-t61.sh', 'run-readout-62b.sh', 'run-final-readouts.sh',
           'stop-attempt1.sh', 'stop-attempt7.sh', 'phase267-publish.py']
def runtime_files():
    files = {}
    for name in ['.phase267-trace.json', '.phase267-exists-original.json', '.phase267-born-identity.stdout.log',
                 'primary267-validate-work/validate-04-driver.json', 'primary267-validate-work/failure-64.json', 'primary267-validate2-work/validate-04-driver.json']:
        files['original/' + name.lstrip('.')] = G/name
    for v in ['base', 'unused', 'metric']: files[f'original/identity-{v}.json'] = G/f'.phase267-identity-results/{v}/identity.json'
    for p in sorted((G/'readout267-identity2-work').glob('*.json')): files['original/readout/' + p.name] = p
    for name in ['.phase267-born-refined.stdout.log', '.phase267-placeholders.json', 'primary267-smoke3-work/smoke-01-driver.json',
                 '.phase267-refined64-logs/status.txt', 'primary267-probe-work/probe-64-fallbacks.json',
                 '.phase267-readout-refined62-failed-exact-inverse-endpoints.stderr.log']:
        files['refined/' + name.lstrip('.')] = R/name
    W = R/'primary267-refined64-work'
    for pattern in ['seg-*-driver.json', 'seg-*-fallbacks.json', 'seg-*-nonlinear-exceptions.json', 'run-summary.json']:
        for p in sorted(W.glob(pattern)): files['refined/run/' + p.name] = p
    for p in sorted(W.glob('attempt*-seg-64/*.json')): files[f'refined/run/{p.parent.name}/{p.name}'] = p
    for p in sorted((R/'.phase267-refined64-logs').glob('exists-seg-*.json')): files['refined/exists/' + p.name] = p
    for n in [61, 62, 63, 64]:
        for p in sorted((R/f'readout267-refined{n}-work').glob('*.json')): files[f'refined/readout-{n}/{p.name}'] = p
    return files
LARGE = ['primary267-refined64-work/sweep-1/photons/seg-61.npz', 'primary267-refined64-work/sweep-1/photons/seg-62.npz',
         'primary267-refined64-work/sweep-1/photons/seg-63.npz', 'primary267-refined64-work/sweep-1/photons/seg-64.npz',
         'readout267-refined64-work/gr/source-64.npz', 'readout267-refined64-work/gr/field-source-64.npz', 'readout267-refined64-work/recovered-64.npz',
         'outputs/direct-eos-gr33/def-native-boundary-layer/bank.npz', 'outputs/direct-eos-gr33/def-native-conservative-rates/thermal-refined/bank.npz',
         'outputs/direct-eos-gr33/def-native-initial-constraints/finite-volume/balanced-initial-state.npz',
         'outputs/direct-eos-gr33/native-retained-completion/evolution/coupled-128.npz', 'outputs/direct-eos-gr33/native-retained-completion/evolution/source-128.npz',
         'native-incident-drive155-work/fields/born-g8.npz']


def package():
    assert not out.exists() and not manifest.exists()
    files = runtime_files()
    for dst, src in files.items(): assert src.is_file(), src
    for name in SCRIPTS: assert (S/name).is_file(), name
    preserved = {p: h for p, h in read(previous_manifest)['sha256'].items() if not p.startswith('docs/')}
    for p, h in preserved.items(): assert sha(root/p) == h, p
    field = lambda p: read(p)['endpoint_compact_charge']
    times = {61: 0.003273442159275794, 62: 0.00332710514549343, 63: 0.0033807681317110662}
    original = {n: field(G/f'readout267-identity2-work/readout-field-{repr(t)}.json') for n, t in times.items()}
    original[64] = field(G/'readout267-identity2-work/readout-field.json')
    refined = {n: field(R/f'readout267-refined{n}-work/readout-field.json') for n in [61, 62, 63, 64]}
    change = {n: (refined[n] - original[n])/abs(original[n]) for n in refined}
    magnitude = {n: (abs(refined[n]) - abs(original[n]))/abs(original[n]) for n in refined}
    depth_r = {r['band']: r['charge'] for r in read(R/'readout267-refined64-work/depth-bands.json')['rows']}
    depth_o = {r['band']: r['charge'] for r in read(G/'.phase266-depth-b64.json')['rows']}
    assert abs(depth_r['all'] - refined[64]) <= 1e-15*abs(refined[64]), (depth_r['all'], refined[64])
    assert abs(depth_o['all'] - original[64]) <= 1e-12*abs(original[64]), (depth_o['all'], original[64])
    cells = {c: (depth_o[f'cells:{c}-{c}'], depth_r[f'cells:{8+2*(c-8)}-{8+2*(c-8)}'] + depth_r[f'cells:{9+2*(c-8)}-{9+2*(c-8)}']) for c in range(8, 16)}
    summary = read(R/'primary267-refined64-work/run-summary.json')
    checks = {n: read(R/f'readout267-refined{n}-work/source-64-check.json') for n in [61, 62, 63, 64]}
    assert all(c['passed'] for c in checks.values()) and all(read(R/f'readout267-refined{n}-work/readout-endpoints.json')['assembly']['passed'] for n in [61, 62, 63, 64])
    inexact = max(read(R/f'readout267-refined{n}-work/readout-{s}.json').get('inexact_inverse_max') or 0. for n in [62, 63, 64] for s in ['endpoints', 'source'])
    assert refined[64] < 0 and all(s['passed'] for s in summary['segments'])
    large = {name: sha(R/name) for name in LARGE}
    captures = sorted((R/'primary267-refined64-work/captures').glob('captured-64-*.npz')); assert len(captures) == 256
    large['primary267-refined64-work/captures (256, sha256 of the ordered list of file hashes)'] = hashlib.sha256('\n'.join(sha(p) for p in captures).encode()).hexdigest()
    out.mkdir(parents=True)
    for dst, src in files.items(): copy(src, out/dst)
    for name in SCRIPTS: copy(S/name, out/'scripts'/name)
    now = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    le, ne = summary['linear_exceptions'], summary['nonlinear_exceptions']
    final = dict(classification='Counterexample candidate', passed=True, verdict='CONDITIONAL_NEGATIVE_CHARGE_KEPT_RESOLUTION_NOT_CONVERGED',
        refined_endpoint_compact_charge=refined[64], original_endpoint_compact_charge=original[64], endpoint_relative_change=change[64], endpoint_magnitude_change=magnitude[64], magnitude_change=magnitude,
        charges_original=original, charges_refined=refined, relative_change=change, t61_all_original_gates=True,
        convergence_criterion=0.02, convergence_met=abs(change[64]) <= 0.02, grid=dict(cells=539, interior=27, refined_cells='8-15 split 2x (34.375/32.375 km)'),
        run=dict(segments=len(summary['segments']), actual_steps=128, seconds=summary['total_seconds'], fallbacks=summary['total_fallbacks']),
        approved_exceptions=dict(scope='macro 61-64 of the refined 64-clock run', user_decisions=['vector linear <=1e-10', 'nonlinear vector defect <=1e-10', 'physical/material gates only, vectors recorded'],
            linear=le, nonlinear=dict(count=len(ne), max_defect=max(e['vector_defect'] for e in ne), max_physical=max(e['physical_moment_max'] for e in ne),
                                      max_material=max(e['material_component_max'] for e in ne))),
        readout_inexact_coordinate_inverse_max=inexact, depth_cells_original_refined=cells,
        not_recomputed=['128 clock', 'one GR return (low, 3.2e-22 of high)', 'exterior photon boundary term (3.9e-8)'],
        large_arrays_sha256=large, final_charge_conclusion='conditional negative charge kept (high, 64 clock); final physical charge unadjudicated',
        full_goal_complete=False, snapshot_KST=now)
    tbl = ''.join(f"| {label} | {original[n]:.4e} | {refined[n]:.4e} | {100*magnitude[n]:+.2f}% | {gate} |\n" for n, label, gate in
                  [(61, 't₆₁=61/64·T', '모든 원 기준'), (62, 't₆₂', '단계 수 246 예외 포함'), (63, 't₆₃', '예외 3단계'), (64, 'T(끝점)', '예외 5단계')])
    dtab = ''.join(f"| {c} | {cells[c][0]:.4e} | {cells[c][1]:.4e} | {100*(abs(cells[c][1])-abs(cells[c][0]))/abs(cells[c][0]):+.1f}% |\n" for c in range(9, 14))
    section = (f"\n\n## 최종 결과 — 세분 64 시계의 끝점 전하 ({now[:16].replace('T', ' ')} KST)\n\n"
        f"분류: Counterexample candidate. **끝점 compact 전하는 세분 격자에서도 음이다({refined[64]:.15e}). 사전 등록한 판정 규칙에 따라, 지배 오차(내부 반경 해상도)를 고친 해에서 조건부 음의 전하 결론은 유지된다.** "
        f"원 64 시계 값 {original[64]:.15e}보다 크기가 {100*abs(change[64]):.2f}% 작다. 2% 기준을 넘으므로 이 수준의 해상도 수렴은 미달이며, 추가 세분은 새 결정 없이 하지 않는다. "
        f"이 끝점은 마지막 구간(macro 61–64)의 벡터 선형·비선형 기준을 사용자 승인에 따라 물리 모멘트·물질 성분 기준(1e−13)으로 대체한 결과다. 모든 원 기준을 지킨 t₆₁ 비교도 같은 방향(−{100*abs(change[61]):.1f}%)이다.\n\n"
        "| 시각 | 원 격자 | 세분 격자 | 크기 변화 | 수락 |\n|---|---|---|---|---|\n" + tbl +
        f"\n분류: Counterexample candidate. 끝점의 깊이 분해(대역 합은 전체와 상대 {abs(read(R/'readout267-refined64-work/depth-bands.json')['band_sum_relative']):.1e}로 일치)에서 변화는 지배 셀에 몰려 있다.\n\n"
        "| 원 셀 | 원 격자 | 세분 격자(두 칸 합) | 크기 변화 |\n|---|---|---|---|\n" + dtab +
        f"\n분류: Counterexample candidate. 수락 기록은 다음과 같다. 구간 {len(summary['segments'])}개, 실제 128단계, {summary['total_seconds']:.0f}초가 걸렸다(실패 시도 제외). macro 0–61 전부와 단계 수 244는 원 기준을 지켰다. "
        f"승인된 예외는 선형 {le['count']}회(벡터 최대 {le['max_vector']:.1e}, 물리 {le['max_physical']:.1e}, 물질 {le['max_material']:.1e})와 비선형 {len(ne)}회다. 비선형 결함은 {min(e['vector_defect'] for e in ne):.1e}–{max(e['vector_defect'] for e in ne):.1e}, 물리 모멘트는 {max(e['physical_moment_max'] for e in ne):.1e} 이하, 물질 성분은 {max(e['material_component_max'] for e in ne):.1e} 이하였다. "
        "에너지 수지 4.7e−16, 종 수지 2.9e−11, 속도 jet 6.2e−9, native 구성·면 수지 등 나머지 기준은 모두 통과했다.\n\n"
        f"분류: Counterexample candidate. 판독 과정은 다음과 같다. 네 시각의 단계244 원천 검사(밀집 단계, 다항식, 끝점, 압력 mapping)와 캡처 조립이 통과했다. 생산 polish가 원 기준으로 수락한 단계 수 244–245의 저장 좌표에는 long double 기체 변수로 정확히 되돌릴 원상이 없었다. 이 판독 첫 시도의 실패는 보존했다. "
        f"판독에서만 가장 가까운 원상을 쓰고 불일치(최대 {inexact:.1e}, 약 1 ulp)를 기록했다. 해 자체는 바꾸지 않았다.\n\n"
        "분류: Conjectural. 남은 오차와 범위는 다음과 같다. 128 시계의 세분 재진화, 1회 GR 반환(low)과 외부 광자 경계 항의 세분 재계산(원 격자 크기 3.2e−22와 3.9e−8), 셀 8–15의 추가 세분이 남았다. 지금 셀 폭 34km는 반응의 e-배 길이 35–57km와 비슷하다. "
        "마지막 단계들에서 벡터 바닥이 커진 원인(약 5e−12에서 1e−8)이 수치 모드인지는 확인하지 않았다. 최종 물리 전하(자기GR, EOS, 비선형, 정적·관측 폐쇄)는 미판정이다.\n")
    text = note.read_text(encoding='utf-8')
    head = f"분류: Counterexample candidate. **결과 요약(2026-09-27): 세분 격자 끝점 전하 {refined[64]:.4e}(원 격자 {original[64]:.4e}), 음의 부호 유지, 크기 {100*magnitude[64]:+.1f}%(2% 수렴 미달). 마지막 구간은 사용자 승인 예외 아래에서 수락했다. 상세는 맨 아래 최종 결과 절.**\n\n"
    title, rest = text.split('\n\n', 1); note.write_text(title + '\n\n' + head + rest.rstrip('\n') + section, encoding='utf-8', newline='\n')
    prefixes = {}
    doc_text = (f"\n\n## 단계267 — 내부 셀을 세분한 primary 재진화와 끝점 전하\n\n"
        f"분류: Counterexample candidate. 단계266이 지배 오차로 지목한 내부 반경 해상도를 고쳤다. 셀 8–15를 2배 세분한 격자(539셀)에서 구동 primary(64 시계)를 같은 최종 방정식으로 t=0부터 끝점까지 진화했다. "
        f"끝점 compact 전하는 {refined[64]:.4e}로 음이고, 원 격자 {original[64]:.4e}보다 크기가 {100*abs(change[64]):.1f}% 작다. 사전 등록 규칙에 따라 조건부 음의 전하 결론은 유지된다. 2% 기준을 넘으므로 이 수준의 해상도 수렴은 미달이다. "
        f"마지막 2.5 macro 단계(실제 5단계)는 벡터 선형·비선형 잔차가 표현 바닥(증폭 약 5×10⁸)에 걸렸다. 그래서 사용자 승인에 따라 물리 모멘트·물질 성분 1e−13 기준으로 수락했다. 모든 원 기준을 지킨 t=61/64·T 비교는 −{100*abs(change[61]):.1f}%다. "
        f"변화는 셀 11(−{100*abs((cells[11][1]-cells[11][0])/cells[11][0]):.1f}%)과 셀 12(−{100*abs((cells[12][1]-cells[12][0])/cells[12][0]):.1f}%)가 주도한다. "
        "재생성 범위의 교훈: 구동 primary는 단계149–157 보정을 쓰지 않았다(입력 2배 교란에도 첫 단계 비트 동일). 존재 검사가 산술 경로를 바꾸는 경우도 감사해야 한다. "
        "최종 물리 전하(128 시계, 1회 GR 반환·외부 광자 항의 세분 재계산, 자기GR·EOS·비선형·관측 폐쇄)는 미판정이다. [근거](../notes/REQUEST267_REFINED_PRIMARY_KO.md).\n")
    for name in DOCS:
        p = root/f'docs/{name}.md'; prefixes[p.relative_to(root).as_posix()] = dict(bytes=p.stat().st_size, sha256=sha(p))
        with p.open('ab') as h: h.write(doc_text.encode())
    write(out/'result.json', final); write(out/'publication.json', dict(previous_nondoc=preserved, document_prefixes=prefixes))
    bound = [note] + [root/p for p in prefixes] + [p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes, sha256={p.relative_to(root).as_posix(): sha(p) for p in bound})
    write(manifest, final); m = read(master); m['sha256'].update(final['sha256'])
    m['sha256'][manifest.relative_to(root).as_posix()] = sha(manifest)
    m['native_refined_primary'] = {k: v for k, v in final.items() if k != 'sha256'}; write(master, m)


def git(*args): return subprocess.check_output(['git', *args], cwd=root)


def check(mode):
    m = read(manifest); a = read(master)
    for p, h in m['sha256'].items(): assert sha(root/p) == h and a['sha256'][p] == h, p
    for p, h in read(out/'publication.json')['previous_nondoc'].items(): assert sha(root/p) == h, p
    for p, v in m['document_prefixes'].items(): assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest() == v['sha256'], p
    assert a['sha256'][manifest.relative_to(root).as_posix()] == sha(manifest)
    paths = list(m['sha256']) + [manifest.relative_to(root).as_posix(), 'paper/revision-manifest.json']
    (S/'phase267-paths.txt').write_bytes(b'\0'.join(p.encode() for p in paths) + b'\0')
    if mode in ['staged', 'head']:
        if mode == 'staged':
            staged = set(filter(None, git('-c', 'core.quotepath=off', 'diff', '--cached', '--name-only', '-z').decode().split('\0')))
            assert staged == set(paths), sorted(staged ^ set(paths))[:10]
        ref = ':' if mode == 'staged' else 'HEAD:'
        for p in paths: assert hashlib.sha256(git('show', ref + p)).hexdigest() == sha(root/p), p
    print(json.dumps(dict(mode=mode, bound_files=len(m['sha256']), paths=len(paths), verdict=m['verdict'], endpoint=m['refined_endpoint_compact_charge'],
                          magnitude_change=m['endpoint_magnitude_change'], final_charge_conclusion=m['final_charge_conclusion'])))


if __name__ == '__main__':
    if sys.argv[1] == 'package': package()
    check(sys.argv[1] if sys.argv[1] != 'package' else 'check')
