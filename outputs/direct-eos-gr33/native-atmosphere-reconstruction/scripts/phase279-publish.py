"""Archive phase 279: direct reconstruction of the observed photosphere (spherical gray envelope with the declared EOS and
Rosseland tables, validated on the declared background) and the non-gray LTE correction of the charge magnitude.

Counterexample candidate (history) / Proven (long-time limit of a stable linear system). Copies the bound computation and the
transit/relaxation record into outputs/direct-eos-gr33/native-closure-transit, binds the two notes, appends one section to the six
dynamic-chi documents, and binds everything in the manifest and paper/revision-manifest.json.
Usage: python phase278-publish.py package|check|staged|head
"""
from pathlib import Path
import datetime, hashlib, json, subprocess, sys
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
R4 = Path('//wsl.localhost/Ubuntu-22.04/home/lpaiu/work/native-refined268-runtime')
S = Path(__file__).resolve().parent
master = root/'paper/revision-manifest.json'; previous_manifest = root/'outputs/direct-eos-gr33/native-closure-transit-manifest.json'
out = root/'outputs/direct-eos-gr33/native-atmosphere-reconstruction'; manifest = out.parent/'native-atmosphere-reconstruction-manifest.json'
notes = [root/'notes/REQUEST279_ATMOSPHERE_RECONSTRUCTION_KO.md']
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
SCRIPTS = ['phase279-gray-pp-cx.py', 'phase279-gray-sph-attempt1.py', 'phase279-gray-sph.py', 'run279a.sh', 'run279s.sh', 'envprobe279.sh', 'diag279.sh', 'diag279b.py',
           'rundiag279b.sh', 'diag279c.sh', 'phase279-nongray.py', 'test279ng.py', 'test279ng2.py', 'test279gray.py', 'phase279-combine.py', 'kappa279.py', 'phase279-publish.py']
FILES = {'phase279/gray-pp-cx.json': R4/'.phase279-gray.json', 'phase279/gray-pp-cx.log': R4/'.phase279-gray.log', 'phase279/gray-pp-cx.npz': R4/'.phase279-gray.npz',
         'phase279/gray-sph-attempt1.log': R4/'.phase279-gray-sph-attempt1.log', 'phase279/gray-sph.json': R4/'.phase279-gray-sph.json',
         'phase279/gray-sph.log': R4/'.phase279-gray-sph.log', 'phase279/gray-sph.npz': R4/'.phase279-gray-sph.npz', 'phase279/diag-matched-pressure-pp-cx.log': R4/'.diag279b.log',
         'phase279/nongray.json': S/'phase279-nongray.json', 'phase279/nongray.log': S/'phase279-nongray.log', 'phase279/combined.json': S/'phase279-combined.json',
         'phase279/kappa.json': S/'phase279-kappa.json'}
FILES.update({f'phase279/{p.name}': p for p in sorted(S.glob('phase279-nongray-*.npz'))})


def package():
    assert not out.exists() and not manifest.exists()
    for src in FILES.values(): assert src.is_file(), src
    for name in SCRIPTS: assert (S/name).is_file(), name
    preserved = {p: h for p, h in read(previous_manifest)['sha256'].items() if not p.startswith('docs/')}
    for p, h in preserved.items(): assert sha(root/p) == h, p
    g = read(FILES['phase279/gray-sph.json']); cb = read(FILES['phase279/combined.json']); kp = read(FILES['phase279/kappa.json'])
    assert g['validation']['kernel_weighted'] < 0.02 and cb['sign_kept_everywhere'] and all(v['nongray_flux_error_tau_le_100'] < 1e-3 for v in cb['cases'].values())
    out.mkdir(parents=True)
    for dst, src in FILES.items(): copy(src, out/dst)
    for name in SCRIPTS: copy(S/name, out/'scripts'/name)
    now = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    kc = cb['cases']['Kaplan central']
    final = dict(classification='Counterexample candidate', passed=True,
        verdict='OBSERVED_PHOTOSPHERE_RECONSTRUCTED__GRAY_TABLE_OPACITY_x3.2__NONGRAY_LTE_x6.9__SIGN_KEPT',
        gray=dict(model='spherical gray Eddington envelope, declared EOS and Rosseland tables, gravitating density cx rho, fixed mass, depth origin at P_gas = 1 dyn/cm^2',
                  validation_max_charge_layers=g['validation']['max_abs_rel_charge_layers'], validation_kernel_weighted=g['validation']['kernel_weighted'],
                  failed_attempts=['pp, origin tau=0: 87%', 'pp, origin P_gas=1: 8.6%', 'pp with cx: 3.6%', 'spherical with depth measured from tau=1e-9: 5.0%'],
                  cases={n: dict(ratio=v['ratio'], tau_at_363km=v['tau_at_363km'], depth_tau1_km=v['depth_tau1_km']) for n, v in g['cases'].items()}),
        nongray=dict(model='plane-parallel LTE H+He continuum, Eddington closure, damped Unsold-Lucy; ratio to the gray model with the same opacity',
                     flux_error_tau_le_100=max(v['nongray_flux_error_tau_le_100'] for v in cb['cases'].values()),
                     factors={n: v['nongray_factor'] for n, v in cb['cases'].items()}, T0_over_Teff=kc['T0_over_Teff_nongray'],
                     simple_over_table_rosseland_charge_layers=kp['charge_layer_ratio_range'], failed_attempt='undamped Lucy with 5% clip oscillated (400 iterations, flux 5e-3)'),
        kaplan_central=dict(gray=kc['gray_ratio'], nongray=kc['nongray_ratio']), one_sigma=cb['one_sigma'], two_sigma=cb['two_sigma'],
        sign_kept_everywhere=cb['sign_kept_everywhere'],
        remaining=['non-LTE (upper-layer heating) and line/metal opacity in the non-gray correction (opposite directions, size not determined)',
                   'plane-parallel geometry of the non-gray ratio'],
        full_goal_complete=False, snapshot_KST=now)
    doc = ("\n\n## 단계279 — 관측 광구의 직접 재구성과 비회색 LTE 보정\n\n"
        f"분류: Counterexample candidate. 선언 EOS·Rosseland 표·중력 밀도 cx·ρ로 구면 회색 Eddington 봉투를 적분했다. 깊이 원점은 선언 봉투와 같은 P_gas=1 dyn/cm² 자름점이다. 이 봉투는 선언 배경의 밀도를 전하 층에서 핵 가중 {100*g['validation']['kernel_weighted']:.2f}%로 재현했다(2% 관문 통과, 앞선 네 번의 실패 87%·8.6%·3.6%·5.0%는 보존). "
        f"질량을 고정한 관측 대기에서 끝점 전하는 선언 모형의 {kc['gray_ratio']:.2f}배(Kaplan 중심), 1σ {cb['one_sigma']['gray'][0]:.2f}–{cb['one_sigma']['gray'][1]:.2f}배, 2σ {cb['two_sigma']['gray'][0]:.2f}–{cb['two_sigma']['gray'][1]:.1f}배다. 이는 단계275의 상사 변환 극한과 맞으며, 두 불투명도 극한 가정을 대체한다. "
        f"같은 H·He 연속 불투명도의 LTE 비회색 복사평형(흐름 일정성 5.5e−5, 표면 T₀/T_eff=0.64)은 전하를 {min(cb['cases'][n]['nongray_factor'] for n in cb['cases']):.1f}–{max(cb['cases'][n]['nongray_factor'] for n in cb['cases']):.1f}배 키운다. 그래서 관측 대기의 전하는 {kc['nongray_ratio']:.1f}배(1σ {cb['one_sigma']['nongray'][0]:.1f}–{cb['one_sigma']['nongray'][1]:.1f}배)다. "
        "부호는 모든 경우에 유지된다. 비회색 보정은 LTE, 선·금속 생략(단순 불투명도의 Rosseland 평균은 표의 0.55–0.74배), Eddington 닫힘에 조건부다. 따라서 크기는 회색 값과 LTE 비회색 값을 함께 적는다. "
        "[근거 279](../notes/REQUEST279_ATMOSPHERE_RECONSTRUCTION_KO.md).\n")
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
    m['native_atmosphere_reconstruction'] = {k_: v for k_, v in final.items() if k_ != 'sha256'}
    master.write_text(json.dumps(m, indent=indent, ensure_ascii=False) + '\n', encoding='utf-8', newline='\n')


def git(*args): return subprocess.check_output(['git', *args], cwd=root)


def check(mode):
    m = read(manifest); a = read(master)
    for p, h in m['sha256'].items(): assert sha(root/p) == h and a['sha256'][p] == h, p
    for p, h in read(out/'publication.json')['previous_nondoc'].items(): assert sha(root/p) == h, p
    for p, v in m['document_prefixes'].items(): assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest() == v['sha256'], p
    assert a['sha256'][manifest.relative_to(root).as_posix()] == sha(manifest)
    paths = list(m['sha256']) + [manifest.relative_to(root).as_posix(), 'paper/revision-manifest.json']
    (S/'phase279-paths.txt').write_bytes(b'\0'.join(p.encode() for p in paths) + b'\0')
    if mode in ['staged', 'head']:
        if mode == 'staged':
            staged = set(filter(None, git('-c', 'core.quotepath=off', 'diff', '--cached', '--name-only', '-z').decode().split('\0')))
            assert staged <= set(paths), sorted(staged - set(paths))[:10]  # bound paths already committed unchanged need not be staged
        ref = ':' if mode == 'staged' else 'HEAD:'
        for p in paths: assert hashlib.sha256(git('show', ref + p)).hexdigest() == sha(root/p), p
    print(json.dumps(dict(mode=mode, bound_files=len(m['sha256']), paths=len(paths), verdict=m['verdict'], validation=m['gray']['validation_kernel_weighted'],
                          central=m['kaplan_central'], one_sigma=m['one_sigma'])))


if __name__ == '__main__':
    if sys.argv[1] == 'package': package()
    check(sys.argv[1] if sys.argv[1] != 'package' else 'check')
