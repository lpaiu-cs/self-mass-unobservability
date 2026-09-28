"""Archive phase 268: the 4x interior-refined primary (cells 8-15 split in four) and the three-resolution comparison.

Counterexample candidate. Copies the small artifacts (preparation logs and checks, segment rows, fallbacks, approved
last-segment exceptions, existence audits, readouts, depth bands, comparison) into
outputs/direct-eos-gr33/native-quad-refined-primary, records SHA-256 of the large runtime arrays without copying them,
binds the phase-267 scripts reused unchanged by their published SHA, appends the result to the phase-268 note and the six
dynamic-chi documents, and binds everything in the manifest and paper/revision-manifest.json.
Usage: python phase268-publish.py package|check|staged|head
"""
from pathlib import Path
import datetime, hashlib, json, subprocess, sys
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
R4 = Path('//wsl.localhost/Ubuntu-22.04/home/lpaiu/work/native-refined268-runtime')
S = Path(__file__).resolve().parent
master = root/'paper/revision-manifest.json'; previous_manifest = root/'outputs/direct-eos-gr33/native-refined-primary-manifest.json'
out = root/'outputs/direct-eos-gr33/native-quad-refined-primary'; manifest = out.parent/'native-quad-refined-primary-manifest.json'
note = root/'notes/REQUEST268_QUADRUPLE_REFINED_PRIMARY_KO.md'; published267 = root/'outputs/direct-eos-gr33/native-refined-primary/scripts'
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
SCRIPTS = ['phase268-bank.py', 'phase268-thermal.py', 'phase268-initial.py', 'phase268-initial-audit.py', 'phase268-checks.py', 'run-prepare268.sh',
           'run-prepare268-resume.sh', 'run-prepare268-resume2.sh', 'run-finish-prep268.sh', 'launch268.sh', 'wait268.sh', 'wait-run268.sh', 'sha-runtime268.sh',
           'run-quad64.sh', 'run-readouts268.sh', 'compare268.py', 'seg268.py', 'phase268-driver.py', 'run-quad64-retry.sh', 'phase268-publish.py']
REUSED = ['phase267-template.py', 'phase267-undriven.py', 'phase267-p150-setup.py', 'phase267-p150.py', 'phase267-born.py', 'phase267-placeholders.py',
          'phase267-driver.py', 'phase267-readout.py', 'phase267-depth.py', 'phase267-audited.py']
def runtime_files():
    files = {}
    for name in ['.phase268-prepare.log', '.phase268-copied-inputs.json', '.phase268-bank.log', '.phase268-thermal-attempt1.log', '.phase268-thermal.log',
                 '.phase268-initial.log', '.phase268-initial-audit.log', '.phase268-p150-audit.json', '.phase268-p150.log', '.phase268-born-attempt1.log',
                 '.phase268-born.log', '.phase267-placeholders.json', '.phase268-exists-smoke.json', 'primary268-smoke-work/smoke-01-driver.json']:
        files['prepare/' + name.split('/')[-1].lstrip('.')] = R4/name
    W = R4/'primary268-quad64-work'
    for pattern in ['seg-*-driver.json', 'seg-*-fallbacks.json', 'seg-*-nonlinear-exceptions.json']:
        for p in sorted(W.glob(pattern)): files['run/' + p.name] = p
    files['run/status.txt'] = R4/'.phase268-quad64-logs/status.txt'
    for p in sorted((W/'attempt1-seg-64').glob('*.json')): files['run/attempt1-seg-64/' + p.name] = p
    for p in sorted((R4/'.phase268-quad64-logs').glob('attempt1-*')): files['run/attempt1-seg-64/' + p.name] = p
    for p in sorted((R4/'.phase268-quad64-logs').glob('exists-seg-*.json')): files['exists/' + p.name] = p
    for n in [61, 62, 63, 64]:
        for p in sorted((R4/f'readout268-quad{n}-work').glob('*.json')): files[f'readout-{n}/{p.name}'] = p
        for p in sorted(R4.glob(f'.phase268-readout-quad{n}-*')): files[f'readout-{n}/{p.name.lstrip(".")}'] = p
    files['readouts.log'] = R4/'.phase268-readouts.log'; files['compare.json'] = R4/'readout268-compare.json'
    return files
LARGE = ['primary268-quad64-work/sweep-1/photons/seg-61.npz', 'primary268-quad64-work/sweep-1/photons/seg-62.npz',
         'primary268-quad64-work/sweep-1/photons/seg-63.npz', 'primary268-quad64-work/sweep-1/photons/seg-64.npz',
         'readout268-quad64-work/gr/source-64.npz', 'readout268-quad64-work/gr/field-source-64.npz', 'readout268-quad64-work/recovered-64.npz',
         'outputs/direct-eos-gr33/def-native-boundary-layer/bank.npz', 'outputs/direct-eos-gr33/def-native-conservative-rates/thermal-refined/bank.npz',
         'outputs/direct-eos-gr33/def-native-initial-constraints/finite-volume/balanced-initial-state.npz',
         'outputs/direct-eos-gr33/native-retained-completion/evolution/coupled-128.npz', 'outputs/direct-eos-gr33/native-retained-completion/evolution/source-128.npz',
         'native-incident-drive155-work/fields/born-g8.npz']
pct = lambda x: f'{100*x:+.2f}%'


def texts(f, now):
    """Result text for the note and the documents, from the final record f."""
    r, conv, neg = f['rows'], f['converged_2x_at_2pct'], f['endpoint_negative_4x']
    e = r['64']; t61 = r['61']; d = f['depth_T']; run = f['run']; le, ne = run['linear_exceptions'], run['nonlinear_exceptions']
    sign = '음이다' if neg else '음이 아니다'
    keep = '사전 등록 규칙에 따라 조건부 음의 전하 결론은 4배 격자에서도 유지된다' if neg else '조건부 음의 전하 결론은 4배 격자에서 유지되지 않는다'
    verdict = ('두 시각 모두 2% 이하이므로, 2배 결과를 2% 수준의 해상도 수렴으로 판정한다' if conv else
               '2% 기준을 넘는 시각이 있으므로, 2배 결과는 2% 수준의 해상도 수렴 미달로 판정한다')
    mono = all(r[n]['monotone'] for n in r)
    order = ', '.join(f"{'T' if n == '64' else 't' + n} {r[n]['observed_order']:.2f}" for n in ['61', '62', '63', '64'])
    rich = (f"관측 차수로 외삽한 Richardson 값은 T에서 {e['richardson']:.4e}, t61에서 {t61['richardson']:.4e}다" if mono and e['richardson'] is not None and t61['richardson'] is not None
            else '차이의 부호가 해상도마다 바뀌어 Richardson 외삽은 의미가 없다')
    tbl = ''.join(f"| {lab} | {r[n]['q1']:.4e} | {r[n]['q2']:.4e} | {r[n]['q4']:.4e} | {pct(r[n]['change_2x_to_4x'])} | {pct(r[n]['change_1x_to_4x'])} | {gate} |\n"
                  for n, lab, gate in [('61', 't₆₁=61/64·T', '모든 원 기준'), ('62', 't₆₂', '마지막 구간 규칙'), ('63', 't₆₃', '마지막 구간 규칙'), ('64', 'T(끝점)', '마지막 구간 규칙')])
    dtab = ''.join(f"| {c} | {d[str(c)][0]:.4e} | {d[str(c)][1]:.4e} | {d[str(c)][2]:.4e} | {pct((abs(d[str(c)][2]) - abs(d[str(c)][1]))/abs(d[str(c)][1]))} |\n" for c in range(9, 14))
    sub10 = d['10 4x sub-cells']; node = min(sub10) < 0 < max(sub10)
    shares = {c: (abs(d[str(c)][2]) - abs(d[str(c)][1]))/abs(d[str(c)][1]) for c in (9, 10, 11)}
    depth_note = (f"지배 셀 11의 2×→4× 변화는 {pct(shares[11])}다. 서로 부호가 반대인 셀 9·10은 변화가 더 크지만({pct(shares[9])}, {pct(shares[10])}) 대부분 상쇄된다. "
                  + ("셀 10의 4배 하위 셀 사이에서 기여의 부호가 바뀐다. " if node else ''))
    inexact = f['readout_inexact_coordinate_inverse_max']
    inexact_text = f"부정확 역변환 최대값은 {inexact:.1e}다" if inexact else '판독 좌표의 부정확 역변환은 없었다'
    exc = (f"마지막 구간 규칙으로 수락한 선형 해는 {le['count']}회(벡터 최대 {le['max_vector']:.1e}, 물리 {le['max_physical']:.1e}, 물질 {le['max_material']:.1e}), 비선형 반복은 {ne['count']}회다"
           if le['count'] or ne['count'] else '마지막 구간 규칙은 한 번도 쓰이지 않았고 모든 단계가 원 기준을 지켰다')
    section = (f"\n\n## 최종 결과 — 4배 세분 64 시계의 끝점 전하 ({now[:16].replace('T', ' ')} KST)\n\n"
        f"분류: Counterexample candidate. **끝점 compact 전하는 4배 격자에서 {e['q4']:.15e}로 {sign}. {keep}.** "
        f"2×→4× 크기 변화는 t₆₁에서 {pct(t61['change_2x_to_4x'])}, T에서 {pct(e['change_2x_to_4x'])}다. {verdict}. 4배를 넘는 세분은 계획대로 하지 않는다.\n\n"
        "| 시각 | 1배 | 2배 | 4배 | 2×→4× 크기 | 1×→4× 크기 | 수락 |\n|---|---|---|---|---|---|---|\n" + tbl +
        f"\n분류: Conjectural. 세 해상도의 관측 차수 p=log₂(|q₂−q₁|/|q₄−q₂|)는 {order}다. {rich}. 이 외삽은 세 격자가 점근 영역에 있다는 가정에 기대며, 연속 극한의 인증이 아니다.\n\n"
        f"분류: Counterexample candidate. 끝점 깊이 분해에서 원 셀별 합은 다음과 같다(4배 대역 합과 전체의 상대 차 {abs(d['closure_4x']):.1e}).\n\n"
        "| 원 셀 | 1배 | 2배(두 칸 합) | 4배(네 칸 합) | 2×→4× 크기 |\n|---|---|---|---|---|\n" + dtab + "\n" + depth_note.rstrip() + "\n" +
        f"\n분류: Counterexample candidate. 수락 기록: 구간 {len(run['segments'])}개, 실제 {run['segments'][-1]['actual_steps']}단계, {run['total_seconds']:.0f}초, 대체 풀이 {run['total_fallbacks']}회. "
        f"macro 0–61은 모든 원 기준을 지켰다. {exc}. seg-64 시도 1은 long double 보정 12회 상한에서 멈췄다(물질 성분 1.1e−12, 보정당 약 10배 감소 중). "
        f"2026-09-24 자원 정책에 따라 보정 상한만 24회로 늘린 드라이버 v6로 다시 실행했고 수락 기준은 바꾸지 않았다(시도 1 보존). 판독의 단계244 원천 검사와 캡처 조립은 네 시각 모두 통과했고, {inexact_text}.\n\n"
        "분류: Conjectural. 남은 오차와 범위: 128 시계의 세분 재진화, 1회 GR 반환(low)과 외부 광자 경계 항의 세분 재계산, 마지막 단계 벡터 바닥 증가의 원인 확인이 남았다. "
        "최종 물리 전하(자기GR, EOS, 비선형, 정적·관측 폐쇄)는 미판정이다.\n")
    head = (f"분류: Counterexample candidate. **결과 요약(2026-09-27): 4배 격자 끝점 전하 {e['q4']:.4e}(2배 {e['q2']:.4e}, 1배 {e['q1']:.4e}), "
            f"{'음의 부호 유지' if neg else '음의 부호 미유지'}, 2×→4× 크기 변화 t₆₁ {pct(t61['change_2x_to_4x'])}·T {pct(e['change_2x_to_4x'])} → "
            f"{'2배 결과 2% 수렴' if conv else '2배 결과 2% 수렴 미달'}. 상세는 맨 아래 최종 결과 절.**\n\n")
    doc = (f"\n\n## 단계268 — 내부 셀 4배 세분의 해상도 수렴 시험\n\n"
        f"분류: Counterexample candidate. 셀 8–15를 4배 세분한 격자(555셀)에서 구동 primary(64 시계)를 같은 최종 방정식으로 t=0부터 끝점까지 진화했다. "
        f"끝점 compact 전하는 {e['q4']:.4e}로 {sign}(2배 {e['q2']:.4e}, 1배 {e['q1']:.4e}). {keep}. "
        f"2×→4× 크기 변화는 t₆₁ {pct(t61['change_2x_to_4x'])}, T {pct(e['change_2x_to_4x'])}이며, {verdict}. "
        f"관측 차수는 {order}다. 최종 물리 전하(128 시계, 1회 GR 반환·외부 광자 항의 세분 재계산, 자기GR·EOS·비선형·관측 폐쇄)는 미판정이다. "
        "[근거](../notes/REQUEST268_QUADRUPLE_REFINED_PRIMARY_KO.md).\n")
    return head, section, doc


def package():
    assert not out.exists() and not manifest.exists()
    files = runtime_files()
    for dst, src in files.items(): assert src.is_file(), src
    for name in SCRIPTS: assert (S/name).is_file(), name
    reused = {}
    for name in REUSED: assert sha(S/name) == sha(published267/name), name; reused[name] = sha(published267/name)
    preserved = {p: h for p, h in read(previous_manifest)['sha256'].items() if not p.startswith('docs/')}
    for p, h in preserved.items(): assert sha(root/p) == h, p
    cmp = read(R4/'readout268-compare.json'); rows = cmp['rows']; assert cmp['decided'] and set(rows) == {'61', '62', '63', '64'}
    assert all(r['same_time_as_2x'] for r in rows.values())
    for n in rows:
        assert read(R4/f'readout268-quad{n}-work/readout-field.json')['endpoint_compact_charge'] == rows[n]['q4']
        assert read(R4/f'readout268-quad{n}-work/source-64-check.json')['passed'] and read(R4/f'readout268-quad{n}-work/readout-endpoints.json')['assembly']['passed']
    assert abs(cmp['depth_T']['all'][2] - rows['64']['q4']) <= 1e-15*abs(rows['64']['q4'])
    assert all(s['passed'] for s in cmp['run']['segments']) and len(cmp['run']['segments']) == 18
    inexact = max(read(R4/f'readout268-quad{n}-work/readout-{s}.json').get('inexact_inverse_max') or 0. for n in rows for s in ['endpoints', 'source'])
    large = {name: sha(R4/name) for name in LARGE}
    captures = sorted((R4/'primary268-quad64-work/captures').glob('captured-64-*.npz')); steps = cmp['run']['segments'][-1]['actual_steps']
    assert len(captures) == 2*steps, (len(captures), steps)
    large[f'primary268-quad64-work/captures ({len(captures)}, sha256 of the ordered list of file hashes)'] = hashlib.sha256('\n'.join(sha(p) for p in captures).encode()).hexdigest()
    out.mkdir(parents=True)
    for dst, src in files.items(): copy(src, out/dst)
    for name in SCRIPTS: copy(S/name, out/'scripts'/name)
    now = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    neg, conv = cmp['endpoint_negative_4x'], cmp['converged_2x_at_2pct']
    verdict = ('CONDITIONAL_NEGATIVE_CHARGE_KEPT_2X_CONVERGED_AT_2PCT' if neg and conv else 'CONDITIONAL_NEGATIVE_CHARGE_KEPT_2X_NOT_CONVERGED_AT_2PCT' if neg
               else 'CONDITIONAL_NEGATIVE_CHARGE_NOT_KEPT_AT_4X')
    final = dict(classification='Counterexample candidate', passed=True, verdict=verdict, rule=cmp['rule'], converged_2x_at_2pct=conv, endpoint_negative_4x=neg,
        rows=rows, grid=dict(cells=555, interior=43, refined_cells='8-15 split 4x (17.19/16.19 km), 2x edges contained bitwise'),
        run=dict(segments=len(cmp['run']['segments']), actual_steps=steps, seconds=cmp['run']['total_seconds'], fallbacks=cmp['run']['total_fallbacks'],
                 linear_exceptions=cmp['run']['linear_exceptions'], nonlinear_exceptions={k: v for k, v in cmp['run']['nonlinear_exceptions'].items() if k != 'rows'},
                 last_segment_rule='registered before the run (REQUEST268): physical moments and material components 1e-13, vectors recorded',
                 seg64_attempt1='stopped at stage 252 after 12 long-double corrections (physical 7.5e-15, material 1.1e-12 > 1e-13, falling about tenfold per correction); '
                                'rerun with driver v6 (24 corrections, resource policy 2026-09-24), every gate unchanged'),
        depth_T=cmp['depth_T'], readout_inexact_coordinate_inverse_max=inexact, reused_phase267_scripts_sha256=reused,
        not_recomputed=['128 clock', 'one GR return (low)', 'exterior photon boundary term'], large_arrays_sha256=large,
        final_charge_conclusion=('conditional negative charge kept at 4x (high, 64 clock); ' if neg else 'conditional negative charge not kept at 4x; ') + 'final physical charge unadjudicated',
        full_goal_complete=False, snapshot_KST=now)
    head, section, doc_text = texts(dict(final, rows=rows, depth_T=cmp['depth_T'], run=cmp['run']), now)
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
    m['native_quad_refined_primary'] = {k: v for k, v in final.items() if k != 'sha256'}
    master.write_text(json.dumps(m, indent=indent, ensure_ascii=False) + '\n', encoding='utf-8', newline='\n')


def git(*args): return subprocess.check_output(['git', *args], cwd=root)


def check(mode):
    m = read(manifest); a = read(master)
    for p, h in m['sha256'].items(): assert sha(root/p) == h and a['sha256'][p] == h, p
    for p, h in read(out/'publication.json')['previous_nondoc'].items(): assert sha(root/p) == h, p
    for p, v in m['document_prefixes'].items(): assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest() == v['sha256'], p
    for name, h in m['reused_phase267_scripts_sha256'].items(): assert sha(published267/name) == h, name
    assert a['sha256'][manifest.relative_to(root).as_posix()] == sha(manifest)
    paths = list(m['sha256']) + [manifest.relative_to(root).as_posix(), 'paper/revision-manifest.json']
    (S/'phase268-paths.txt').write_bytes(b'\0'.join(p.encode() for p in paths) + b'\0')
    if mode in ['staged', 'head']:
        if mode == 'staged':
            staged = set(filter(None, git('-c', 'core.quotepath=off', 'diff', '--cached', '--name-only', '-z').decode().split('\0')))
            assert staged == set(paths), sorted(staged ^ set(paths))[:10]
        ref = ':' if mode == 'staged' else 'HEAD:'
        for p in paths: assert hashlib.sha256(git('show', ref + p)).hexdigest() == sha(root/p), p
    print(json.dumps(dict(mode=mode, bound_files=len(m['sha256']), paths=len(paths), verdict=m['verdict'], endpoint=m['rows']['64']['q4'],
                          change_2x_to_4x={n: r['change_2x_to_4x'] for n, r in m['rows'].items()}, final_charge_conclusion=m['final_charge_conclusion'])))


if __name__ == '__main__':
    if sys.argv[1] == 'package': package()
    check(sys.argv[1] if sys.argv[1] != 'package' else 'check')
