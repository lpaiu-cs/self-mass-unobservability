"""Archive phases 270-271: physical composition of the endpoint charge and its free-fall reproduction.

Counterexample candidate. Copies the diagnostic records (component/part decomposition of the 4x endpoint charge,
baryon conservation and frozen-memory diagnostics, free-fall comparisons and charges on the 1x/2x/4x grids, the
rejected delay-arrival test) into outputs/direct-eos-gr33/native-charge-mechanism, records SHA-256 of the readout
sources used, appends the result to the phase-271 note and the six dynamic-chi documents, and binds everything in the
manifest and paper/revision-manifest.json. The phase-270 note is bound as committed.
Usage: python phase271-publish.py package|check|staged|head
"""
from pathlib import Path
import datetime, hashlib, json, math, subprocess, sys
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
WSL = Path('//wsl.localhost/Ubuntu-22.04/home/lpaiu/work')
R1, R2, R4 = WSL/'native-retained-tail-runtime', WSL/'native-refined267-runtime', WSL/'native-refined268-runtime'
S = Path(__file__).resolve().parent
master = root/'paper/revision-manifest.json'; previous_manifest = root/'outputs/direct-eos-gr33/native-eos-sensitivity-manifest.json'
out = root/'outputs/direct-eos-gr33/native-charge-mechanism'; manifest = out.parent/'native-charge-mechanism-manifest.json'
note270, note = root/'notes/REQUEST270_CHARGE_COMPOSITION_KO.md', root/'notes/REQUEST271_FREEFALL_REPRODUCTION_KO.md'
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
SCRIPTS = ['phase270-components.py', 'run-components270.sh', 'phase270-baryon.py', 'phase270-memory.py', 'run-diag270.sh', 'find-setup.sh', 'find-setup2.sh',
           'find-setup3.sh', 'inspect-ff.sh', 'phase271-freefall.py', 'run-freefall271.sh', 'run-freefall271-grids.sh', 'run-freefall271-all.sh', 'phase271-publish.py']
GRIDS = {'1x': (R1, 'readout267-identity2-work', 'original'), '2x': (R2, 'readout267-refined64-work', 'refined64'), '4x': (R4, 'readout268-quad64-work', 'quad64')}
def runtime_files():
    files = {'phase270/components-T.json': R4/'readout268-quad64-work/components-T.json', 'phase270/components.stdout.log': R4/'.phase270-components.stdout.log',
             'phase270/baryon-conservation.txt': R4/'.phase270-baryon.txt', 'phase270/frozen-memory.txt': R4/'.phase270-memory.txt',
             'phase271/4x-delay-arrival-compare.json': R4/'readout268-quad64-work/freefall-compare-quad64-delay.json'}
    for g, (rt, folder, label) in GRIDS.items():
        for m in ['compare', 'charge']: files[f'phase271/{g}-{m}.json'] = rt/folder/f'freefall-{m}-{label}.json'
    return files


def package():
    assert not out.exists() and not manifest.exists()
    files = runtime_files()
    for dst, src in files.items(): assert src.is_file(), src
    for name in SCRIPTS: assert (S/name).is_file(), name
    preserved = {p: h for p, h in read(previous_manifest)['sha256'].items() if not p.startswith('docs/')}
    for p, h in preserved.items(): assert sha(root/p) == h, p
    comp = read(files['phase270/components-T.json']); rows = {r['spec']: r['charge'] for r in comp['rows']}
    assert abs(comp['state_plus_geometry_relative']) < 1e-12
    charges = {g: read(files[f'phase271/{g}-charge.json']) for g in GRIDS}
    q_ff = [charges[g]['freefall_interior'] for g in GRIDS]; q_cb = [charges[g]['chain_baryon_only'] for g in GRIDS]
    q_ch = [read(R1/'readout267-identity2-work/readout-field.json')['endpoint_compact_charge'], read(R2/'readout267-refined64-work/readout-field.json')['endpoint_compact_charge'],
            read(R4/'readout268-quad64-work/readout-field.json')['endpoint_compact_charge']]
    def rich(q):
        d1, d2 = q[1] - q[0], q[2] - q[1]; p = math.log2(abs(d1)/abs(d2)); return dict(order=p, limit=q[2] + d2/(2**p - 1), monotone=d1*d2 > 0)
    rf, rc = rich(q_ff), rich(q_ch); assert rf['monotone'] and rc['monotone']
    compare = {g: read(files[f'phase271/{g}-compare.json']) for g in GRIDS}
    large = {f'{g}: {folder}/gr/{n}': sha(rt/folder/'gr'/n) for g, (rt, folder, _) in GRIDS.items() for n in ['source-64.npz', 'field-source-64.npz']}
    out.mkdir(parents=True)
    for dst, src in files.items(): copy(src, out/dst)
    for name in SCRIPTS: copy(S/name, out/'scripts'/name)
    now = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    final = dict(classification='Counterexample candidate', passed=True, verdict='ENDPOINT_CHARGE_IS_FREE_FALL_BARYON_MEMORY_REPRODUCED',
        composition_4x_T=dict(total=rows['all'], state=rows['state'], geometry=rows['geometry'], baryon_state=rows['state:baryon_g'], metric_stress_state=rows['state:metric_stress_erg'],
                              nonrest_trace_state=rows['state:nonrest_trace_erg'], state_share=comp['state_share'], geometry_share=comp['geometry_share'],
                              closure=comp['state_plus_geometry_relative']),
        freefall=dict(model='d2xi/dt2 = -c^2 (a^2/B^2) d/dr[alpha(phi0) dphi], alpha=-4 phi, dphi = eta (R_d/r) f((t+x/c)/D), f=(4p(1-p))^4; Phi=4 pi r^2 B rho0 xi',
                      charges={g: dict(chain_total=q_ch[i], chain_baryon_only=q_cb[i], freefall=q_ff[i], relative=q_ff[i]/q_cb[i] - 1) for i, g in enumerate(GRIDS)},
                      richardson_freefall=rf, richardson_chain_total=rc, limits_relative=rf['limit']/rc['limit'] - 1,
                      dM_T_relative_ranges={g: [min(abs(r['relative_T']) for r in compare[g]['rows'] if r['cell'] >= {'1x': 10, '2x': 12, '4x': 12}[g]),
                                                max(abs(r['relative_T']) for r in compare[g]['rows'] if r['cell'] >= {'1x': 10, '2x': 12, '4x': 12}[g])] for g in GRIDS},
                      rejected=['pulse arrival from the delay array (constant +4.0e-6 s): larger deviations']),
        sources_sha256=large, not_modelled=['atmosphere cells (chain values kept; charge share 6e-4)'],
        final_charge_conclusion='conditional negative charge; continuum limit about -2.335e-51 from two independent routes; final physical charge unadjudicated',
        full_goal_complete=False, snapshot_KST=now)
    pct = lambda x: f'{100*x:+.2f}%'
    tbl = ''.join(f"| {g} | {q_ch[i]:.6e} | {q_ff[i]:.6e} | {pct(q_ff[i]/q_cb[i] - 1)} |\n" for i, g in enumerate(GRIDS))
    section = (f"\n\n## 게시 요약 ({now[:16].replace('T', ' ')} KST)\n\n"
        f"분류: Counterexample candidate. 단계270 분해: 끝점 전하의 상태 부분 비중 {comp['state_share']:.9f}, 바리온 상태 성분 {rows['state:baryon_g']/rows['all']:.6f}, 펄스×배경 기하 {comp['geometry_share']:.1e}. "
        f"단계271 자유낙하: 연속 극한 {rf['limit']:.6e}(p={rf['order']:.2f}), 계산 줄기 {rc['limit']:.6e}(p={rc['order']:.2f}), 상대 {final['freefall']['limits_relative']:+.1e}.\n\n"
        "| 격자 | 계산 줄기 전체 | 자유낙하 | 자유낙하/줄기 바리온 성분 |\n|---|---|---|---|\n" + tbl +
        "\n기록은 `outputs/direct-eos-gr33/native-charge-mechanism/`에 있다.\n")
    note.write_text(note.read_text(encoding='utf-8').rstrip('\n') + section, encoding='utf-8', newline='\n')
    doc = (f"\n\n## 단계270–271 — 끝점 전하의 물리적 구성과 자유낙하 재현\n\n"
        f"분류: Counterexample candidate. 4배 해의 끝점 compact 전하를 원천 성분·부분별로 정확히 분해했다(닫힘 {abs(comp['state_plus_geometry_relative']):.0e}). "
        f"전하는 바리온 질량 섭동의 상태 응답이 {rows['state:baryon_g']/rows['all']:.6f}를 차지하고, 입사장×배경의 직접 결합은 {comp['geometry_share']:.0e}이며, 열·광자 성분은 3e−9 이하다. "
        "바리온 섭동은 질량을 보존하는 재배치이며, 끝점 전하의 99.9%는 펄스가 이미 떠난 층의 동결 변위(기억)에서 온다. "
        "선언 이론(A=exp(−2φ²))에서 유도한 무압력 자유낙하 모형은 지배 셀의 끝점 δM을 4배에서 6e−6 이내로 재현했다. "
        f"그 전하의 연속 극한 {rf['limit']:.4e}는 계산 줄기의 극한 {rc['limit']:.4e}와 상대 {final['freefall']['limits_relative']:+.0e}로 일치한다. "
        "전하는 배경 ρ₀·φ₀와 입사 펄스의 명시적 범함수이며, 최종 물리 전하는 미판정이다. "
        "[근거 270](../notes/REQUEST270_CHARGE_COMPOSITION_KO.md), [근거 271](../notes/REQUEST271_FREEFALL_REPRODUCTION_KO.md).\n")
    prefixes = {}
    for name in DOCS:
        p = root/f'docs/{name}.md'; prefixes[p.relative_to(root).as_posix()] = dict(bytes=p.stat().st_size, sha256=sha(p))
        with p.open('ab') as h: h.write(doc.encode())
    write(out/'result.json', final); write(out/'publication.json', dict(previous_nondoc=preserved, document_prefixes=prefixes))
    bound = [note270, note] + [root/p for p in prefixes] + [p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes, sha256={p.relative_to(root).as_posix(): sha(p) for p in bound})
    write(manifest, final)
    raw = master.read_text(encoding='utf-8'); m = json.loads(raw); indent = 2 if raw.startswith('{\n  "') else 1
    m['sha256'].update(final['sha256']); m['sha256'][manifest.relative_to(root).as_posix()] = sha(manifest)
    m['native_charge_mechanism'] = {k: v for k, v in final.items() if k != 'sha256'}
    master.write_text(json.dumps(m, indent=indent, ensure_ascii=False) + '\n', encoding='utf-8', newline='\n')


def git(*args): return subprocess.check_output(['git', *args], cwd=root)


def check(mode):
    m = read(manifest); a = read(master)
    for p, h in m['sha256'].items(): assert sha(root/p) == h and a['sha256'][p] == h, p
    for p, h in read(out/'publication.json')['previous_nondoc'].items(): assert sha(root/p) == h, p
    for p, v in m['document_prefixes'].items(): assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest() == v['sha256'], p
    assert a['sha256'][manifest.relative_to(root).as_posix()] == sha(manifest)
    paths = list(m['sha256']) + [manifest.relative_to(root).as_posix(), 'paper/revision-manifest.json']
    (S/'phase271-paths.txt').write_bytes(b'\0'.join(p.encode() for p in paths) + b'\0')
    if mode in ['staged', 'head']:
        if mode == 'staged':
            staged = set(filter(None, git('-c', 'core.quotepath=off', 'diff', '--cached', '--name-only', '-z').decode().split('\0')))
            assert staged <= set(paths) and set(paths) - staged <= {'notes/REQUEST270_CHARGE_COMPOSITION_KO.md'}, sorted(staged ^ set(paths))[:10]
        ref = ':' if mode == 'staged' else 'HEAD:'
        for p in paths: assert hashlib.sha256(git('show', ref + p)).hexdigest() == sha(root/p), p
    f = m['freefall']
    print(json.dumps(dict(mode=mode, bound_files=len(m['sha256']), paths=len(paths), verdict=m['verdict'], limit_freefall=f['richardson_freefall']['limit'],
                          limit_chain=f['richardson_chain_total']['limit'], limits_relative=f['limits_relative'])))


if __name__ == '__main__':
    if sys.argv[1] == 'package': package()
    check(sys.argv[1] if sys.argv[1] != 'package' else 'check')
