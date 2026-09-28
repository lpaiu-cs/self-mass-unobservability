"""Archive phases 272-273: background-structure density kernel of the free-fall charge, and the static-EFT comparison with
the adiabatic collapse boundary.

Counterexample candidate. Copies the kernel records, the charge-history decomposition, the restoring-timescale and
long-wavelength force computation into outputs/direct-eos-gr33/native-structure-eft-boundary, binds the two notes, appends
one section to the six dynamic-chi documents, and binds everything in the manifest and paper/revision-manifest.json.
Usage: python phase273-publish.py package|check|staged|head
"""
from pathlib import Path
import datetime, hashlib, json, subprocess, sys
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
R4 = Path('//wsl.localhost/Ubuntu-22.04/home/lpaiu/work/native-refined268-runtime'); O = R4/'readout268-quad64-work'
S = Path(__file__).resolve().parent
master = root/'paper/revision-manifest.json'; previous_manifest = root/'outputs/direct-eos-gr33/native-charge-mechanism-manifest.json'
out = root/'outputs/direct-eos-gr33/native-structure-eft-boundary'; manifest = out.parent/'native-structure-eft-boundary-manifest.json'
notes = [root/'notes/REQUEST272_BACKGROUND_STRUCTURE_KO.md', root/'notes/REQUEST273_STATIC_EFT_BOUNDARY_KO.md']
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
SCRIPTS = ['phase272-kernel.py', 'run-kernel272.sh', 'phase273-history.py', 'run-history273.sh', 'phase273-history-summary.py', 'phase273-modes.py',
           'run-modes273.sh', 'inspect-bg.sh', 'phase273-publish.py']
def runtime_files():
    files = {'phase272/kernel.json': O/'kernel-quad64.json', 'phase272/kernel-full.json': O/'kernel-quad64-full.json', 'phase272/analyze.log': R4/'.phase272-analyze.log'}
    for p in sorted(O.glob('kernel-quad64-*-*.json')): files['phase272/' + p.name] = p
    for p in sorted(R4.glob('.phase272-*-quad64.log')): files['phase272/' + p.name.lstrip('.')] = p
    for m in ['all', 'state', 'geometry']: files[f'phase273/history-{m}.json'] = O/f'history-{m}.json'
    files['phase273/history-summary.json'] = O/'history-summary.json'; files['phase273/modes.json'] = R4/'.phase273-modes.json'; files['phase273/modes.log'] = R4/'.phase273-modes.log'
    return files


def package():
    assert not out.exists() and not manifest.exists()
    files = runtime_files()
    for dst, src in files.items(): assert src.is_file(), src
    for name in SCRIPTS: assert (S/name).is_file(), name
    preserved = {p: h for p, h in read(previous_manifest)['sha256'].items() if not p.startswith('docs/')}
    for p, h in preserved.items(): assert sha(root/p) == h, p
    k = read(files['phase272/kernel.json']); hs = read(files['phase273/history-summary.json']); md = read(files['phase273/modes.json'])
    assert abs(k['additivity_relative']) < 1e-12 and all(v['sign_kept'] for v in k['families'].values()) and hs['closure'] < 1e-12
    out.mkdir(parents=True)
    for dst, src in files.items(): copy(src, out/dst)
    for name in SCRIPTS: copy(S/name, out/'scripts'/name)
    now = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    fam = {n: v['relative'] for n, v in k['families'].items()}
    final = dict(classification='Counterexample candidate', passed=True, verdict='SIGN_STRUCTURALLY_ROBUST_MAGNITUDE_DENSITY_SET__MEMORY_COLLAPSES_AT_ORBITAL_TIMESCALES',
        background_structure=dict(interior_freefall_charge=k['interior_freefall_charge'], additivity=k['additivity_relative'], sign_margin_Linf=k['sign_margin_Linf'],
                                  opposite_sign_share=k['positive_share'], centroid_depth_km=k['centroid_depth_km'], families_relative=fam),
        static_eft=dict(history_closure=hs['closure'], state_share_min=min(r['share'] for r in hs['state_share_at'] if r['share'] is not None),
                        radial_mode_periods_s=md['radial_modes']['period_s'], omega0=md['radial_modes']['omega'][0], dynamical_time_s=md['radial_modes']['dynamical_time_s'],
                        charge_layers=md['charge_layers'], long_wavelength_force_ratio=md['long_wavelength_force_ratio']['ratio'], orbital_drives=md['orbital_drives']),
        not_covered=['non-radial (tidal, l>=2) displacement channels', 'dissipation (quadrature)', 'background metric/scalar response to density changes (<1e-5 in the envelope)'],
        final_charge_conclusion='conditional negative charge: sign structurally robust to envelope structure; magnitude set by the envelope density at 240-465 km; the free-fall memory collapses to static coefficients at orbital timescales',
        full_goal_complete=False, snapshot_KST=now)
    doc = (f"\n\n## 단계272–273 — 배경 구조 민감도와 정적 EFT 붕괴 경계\n\n"
        f"분류: Counterexample candidate. 자유낙하 전하는 면 밀도에 정확히 선형이다. 4배 격자의 면 핵(합의 닫힘 {abs(k['additivity_relative']):.0e})은 거의 모두 전하와 같은 부호이며(반대 부호 몫 {abs(k['positive_share']):.1e}), 깊이 240–465 km에 모인다. "
        f"부호를 뒤집으려면 면 밀도의 최대노름 상대 변화가 {k['sign_margin_Linf']:.4f}여야 하므로, 밀도가 양수인 한 어떤 봉투 분포도 부호를 뒤집지 못한다. "
        f"크기는 이 깊이의 밀도에 비례하여, 분포가 10 km 어긋나면 약 {100*fam['shift +10 km']:.0f}%(바깥쪽)/{100*fam['shift -10 km']:.0f}%(안쪽) 바뀐다. "
        "전하 이력(128시각)에서 동적 변위 부분의 비중은 모든 시각에 1−1e−8이며, 전하는 펄스가 표면을 떠난 뒤에도 더 깊은 층의 변위로 약 600배 자란다. "
        f"배경 항성의 가장 낮은 반경 단열 모드 주기는 {md['radial_modes']['period_s'][0]:.0f} s, 전하 층의 음향 차단 주기는 28–47 s다. "
        f"긴 파장(λ≫R) 단극 힘은 펄스 힘의 {md['long_wavelength_force_ratio']['ratio'][0]:.1e}배이고, J0337 내측 궤도에서 (ω/ω₀)²={md['orbital_drives']['1.629 d (J0337 inner)']['over_omega0_squared']:.1e}이다. "
        "따라서 이 자유낙하 기억은 궤도 시간척도에서 정적 계수로 붕괴한다(no-go 경계, 조석 채널·소산은 미포함). "
        "[근거 272](../notes/REQUEST272_BACKGROUND_STRUCTURE_KO.md), [근거 273](../notes/REQUEST273_STATIC_EFT_BOUNDARY_KO.md).\n")
    prefixes = {}
    for name in DOCS:
        p = root/f'docs/{name}.md'; prefixes[p.relative_to(root).as_posix()] = dict(bytes=p.stat().st_size, sha256=sha(p))
        with p.open('ab') as h: h.write(doc.encode())
    write(out/'result.json', final); write(out/'publication.json', dict(previous_nondoc=preserved, document_prefixes=prefixes))
    bound = notes + [root/p for p in prefixes] + [p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes, sha256={p.relative_to(root).as_posix(): sha(p) for p in bound})
    write(manifest, final)
    raw = master.read_text(encoding='utf-8'); m = json.loads(raw); indent = 2 if raw.startswith('{\n  "') else 1
    m['sha256'].update(final['sha256']); m['sha256'][manifest.relative_to(root).as_posix()] = sha(manifest)
    m['native_structure_eft_boundary'] = {k_: v for k_, v in final.items() if k_ != 'sha256'}
    master.write_text(json.dumps(m, indent=indent, ensure_ascii=False) + '\n', encoding='utf-8', newline='\n')


def git(*args): return subprocess.check_output(['git', *args], cwd=root)


def check(mode):
    m = read(manifest); a = read(master)
    for p, h in m['sha256'].items(): assert sha(root/p) == h and a['sha256'][p] == h, p
    for p, h in read(out/'publication.json')['previous_nondoc'].items(): assert sha(root/p) == h, p
    for p, v in m['document_prefixes'].items(): assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest() == v['sha256'], p
    assert a['sha256'][manifest.relative_to(root).as_posix()] == sha(manifest)
    paths = list(m['sha256']) + [manifest.relative_to(root).as_posix(), 'paper/revision-manifest.json']
    (S/'phase273-paths.txt').write_bytes(b'\0'.join(p.encode() for p in paths) + b'\0')
    if mode in ['staged', 'head']:
        if mode == 'staged':
            staged = set(filter(None, git('-c', 'core.quotepath=off', 'diff', '--cached', '--name-only', '-z').decode().split('\0')))
            assert staged == set(paths), sorted(staged ^ set(paths))[:10]
        ref = ':' if mode == 'staged' else 'HEAD:'
        for p in paths: assert hashlib.sha256(git('show', ref + p)).hexdigest() == sha(root/p), p
    print(json.dumps(dict(mode=mode, bound_files=len(m['sha256']), paths=len(paths), verdict=m['verdict'], sign_margin=m['background_structure']['sign_margin_Linf'],
                          omega0=m['static_eft']['omega0'], force_ratio=m['static_eft']['long_wavelength_force_ratio'])))


if __name__ == '__main__':
    if sys.argv[1] == 'package': package()
    check(sys.argv[1] if sys.argv[1] != 'package' else 'check')
