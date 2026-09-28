"""Archive phases 277-278: closure of the self-consistent GR / ADM / infinity normalization and nonlinearity bounds, and the
whole-star transit and long-time relaxation history of the compact charge.

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
master = root/'paper/revision-manifest.json'; previous_manifest = root/'outputs/direct-eos-gr33/native-final-closure-manifest.json'
out = root/'outputs/direct-eos-gr33/native-closure-transit'; manifest = out.parent/'native-closure-transit-manifest.json'
notes = [root/'notes/REQUEST277_GR_NONLINEAR_CLOSURE_KO.md', root/'notes/REQUEST278_TRANSIT_LONG_TIME_KO.md']
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
SCRIPTS = ['phase277-bounds.py', 'run277.sh', 'phase278-transit.py', 'run278b.sh', 'phase278-publish.py']
FILES = {'phase277/bounds.json': R4/'.phase277-bounds.json', 'phase277/bounds.log': R4/'.phase277-bounds.log',
         'phase278/transit.json': R4/'.phase278-transit.json', 'phase278/transit.log': R4/'.phase278-transit.log'}


def package():
    assert not out.exists() and not manifest.exists()
    for src in FILES.values(): assert src.is_file(), src
    for name in SCRIPTS: assert (S/name).is_file(), name
    preserved = {p: h for p, h in read(previous_manifest)['sha256'].items() if not p.startswith('docs/')}
    for p, h in preserved.items(): assert sha(root/p) == h, p
    b = read(FILES['phase277/bounds.json']); t = read(FILES['phase278/transit.json'])
    assert b['largest_correction'] < 1e-4 and abs(t['validation']['relative']) < 0.01 and t['validation']['freefall_displacement_rel_Linf_charge_layers'] < 1e-4
    out.mkdir(parents=True)
    for dst, src in FILES.items(): copy(src, out/dst)
    for name in SCRIPTS: copy(S/name, out/'scripts'/name)
    now = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    tr, rl = t['transit'], t['relaxation']
    final = dict(classification='Counterexample candidate (history); Proven (long-time limit of a stable linear system)', passed=True,
        verdict='GR_NONLINEAR_CLOSED__ENDPOINT_CHARGE_IS_EARLY_RETARDED_SCATTERING__PERMANENT_CHARGE_ZERO',
        closure=dict(fixed_point_residual=b['fixed_point']['relative_residual_bound'], adm_cross=b['adm']['exterior_cross_energy_relative'],
                     infinity_tail_bound=b['infinity']['schwarzschild_tail_bound'], infinity_measured=b['infinity']['phase252_exterior_scalar_relative'],
                     nonlinear=b['nonlinear']['dphi_over_phi0'], mass_normalization=b['adm']['mass_normalization_measured'], continuum_mass_normalized=b['continuum']['mass_normalized'],
                     phi_inf='q/eta even in phi_inf (parity); state part proportional to phi_inf^2 (source weight alpha(phi0)/r); direct part (1e-8) exponent open'),
        history=dict(validation=t['validation'], center_arrival_s=t['center_arrival_s'], exit_s=t['exit_s'], transit_max_abs=tr['q_max_abs'], transit_t_max=tr['t_max_abs'],
                     transit_sign_changes=tr['sign_changes'], relaxation_rms=rl['q_rms'], relaxation_max=rl['q_max_abs'], relaxation_mean_2s_1e4s=rl['mean_first_to_1e4'],
                     dominant_periods_s=[m_['period_s'] for m_ in rl['dominant_modes'][:4]], lowest_periods_s=rl['period_lowest'][:3], relaxation_check=t.get('relaxation_check')),
        final_charge_conclusion=('endpoint (T=2D) charge negative, and negative for every endpoint up to 0.35 s; it is the earliest photospheric part of the '
                                 'matter-mediated monopole scattering of the pulse, which oscillates in sign on the outgoing pass (max 3.5e-36, 1.5e-11 of the incident '
                                 'amplitude) and then rings in radial p modes with zero mean; the permanent charge is zero at first order in eta for positive mode damping'),
        remaining=['radial-mode damping rates not computed (stability assumed)', 'interior Born scattering neglected (core compactness ~1e-4)',
                   'direct-part phi_inf exponent', 'atmosphere: gray, two opacity limits (phase 279)'],
        full_goal_complete=False, snapshot_KST=now)
    doc = ("\n\n## 단계277–278 — 자기 일관 GR·비선형의 폐쇄와 전하의 전체 시간 이력\n\n"
        f"분류: Counterexample candidate. 자기 일관 GR 고정점 잔차는 {b['fixed_point']['relative_residual_bound']:.1e} 이하(부등식 Proven), ADM 교차 에너지는 {b['adm']['exterior_cross_energy_relative']:.1e}, 무한대 꼬리는 {b['infinity']['schwarzschild_tail_bound']:.1e} 이하(측정 {b['infinity']['phase252_exterior_scalar_relative']:.1e}), 비선형은 {b['nonlinear']['dphi_over_phi0']:.1e}(α=βφ 정확, Proven)다. "
        f"질량 정규화를 적용한 연속 극한은 {b['continuum']['mass_normalized']:.5e}다. 판독 원천의 셀 가중치가 α(φ₀)/r이므로 상태 부분은 φ_∞²에 비례하고, 짝함수 성질로 q/η는 φ_∞의 짝함수다. "
        f"선형 단열 반경 모형(전체 별, 중심 반사 포함)은 끝점 전하를 {100*t['validation']['relative']:+.2f}%로 재현했다. 끝점 값은 정적 창 값의 {abs(t['validation']['q_T']/t['validation']['static_window_T']):.0f}배인 지연 단극장이다. "
        f"들어오는 구간(0–0.35 s)에서 전하는 음수로 약 12자릿수 자라고(중심 도달 −7.5e−41), 나가는 통과에서 부호가 진동하며 출사 시각 {tr['t_max_abs']:.3f} s에 최대 {tr['q_max_abs']:.1e}(입사 진폭의 1.5e−11)다. "
        f"펄스가 떠난 뒤 반경 p모드(주요 주기 44–59 s)로 진동하며(rms {rl['q_rms']:.1e}), 시간 평균은 0이다. 양의 모드 감쇠에서 영구 전하는 η의 1차에서 0이다(Proven). "
        "따라서 끝점 T=2D에 묶인 조건을 닫는다. 음의 끝점 결론은 0.35 s까지의 모든 끝점으로 넓어지며, 영구 전하는 0이다. "
        "[근거 277](../notes/REQUEST277_GR_NONLINEAR_CLOSURE_KO.md), [근거 278](../notes/REQUEST278_TRANSIT_LONG_TIME_KO.md).\n")
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
    m['native_closure_transit'] = {k_: v for k_, v in final.items() if k_ != 'sha256'}
    master.write_text(json.dumps(m, indent=indent, ensure_ascii=False) + '\n', encoding='utf-8', newline='\n')


def git(*args): return subprocess.check_output(['git', *args], cwd=root)


def check(mode):
    m = read(manifest); a = read(master)
    for p, h in m['sha256'].items(): assert sha(root/p) == h and a['sha256'][p] == h, p
    for p, h in read(out/'publication.json')['previous_nondoc'].items(): assert sha(root/p) == h, p
    for p, v in m['document_prefixes'].items(): assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest() == v['sha256'], p
    assert a['sha256'][manifest.relative_to(root).as_posix()] == sha(manifest)
    paths = list(m['sha256']) + [manifest.relative_to(root).as_posix(), 'paper/revision-manifest.json']
    (S/'phase278-paths.txt').write_bytes(b'\0'.join(p.encode() for p in paths) + b'\0')
    if mode in ['staged', 'head']:
        if mode == 'staged':
            staged = set(filter(None, git('-c', 'core.quotepath=off', 'diff', '--cached', '--name-only', '-z').decode().split('\0')))
            assert staged <= set(paths), sorted(staged - set(paths))[:10]  # bound paths already committed unchanged need not be staged
        ref = ':' if mode == 'staged' else 'HEAD:'
        for p in paths: assert hashlib.sha256(git('show', ref + p)).hexdigest() == sha(root/p), p
    print(json.dumps(dict(mode=mode, bound_files=len(m['sha256']), paths=len(paths), verdict=m['verdict'], validation=m['history']['validation']['relative'],
                          transit_max=m['history']['transit_max_abs'], relax_rms=m['history']['relaxation_rms'])))


if __name__ == '__main__':
    if sys.argv[1] == 'package': package()
    check(sys.argv[1] if sys.argv[1] != 'package' else 'check')
