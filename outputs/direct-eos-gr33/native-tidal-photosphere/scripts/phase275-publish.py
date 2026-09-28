"""Archive phases 274-275: tidal (l>=2) and dissipative channels of the free-fall charge at orbital timescales, and the
photospheric-density uncertainty of the charge magnitude.

Counterexample candidate. Copies the tidal/dissipation computation, the atmosphere mapping of the charge layers and the
photosphere variants into outputs/direct-eos-gr33/native-tidal-photosphere, binds the two notes, appends one section to the six
dynamic-chi documents, and binds everything in the manifest and paper/revision-manifest.json.
Usage: python phase275-publish.py package|check|staged|head
"""
from pathlib import Path
import datetime, hashlib, json, subprocess, sys
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
R4 = Path('//wsl.localhost/Ubuntu-22.04/home/lpaiu/work/native-refined268-runtime')
S = Path(__file__).resolve().parent
master = root/'paper/revision-manifest.json'; previous_manifest = root/'outputs/direct-eos-gr33/native-structure-eft-boundary-manifest.json'
out = root/'outputs/direct-eos-gr33/native-tidal-photosphere'; manifest = out.parent/'native-tidal-photosphere-manifest.json'
notes = [root/'notes/REQUEST274_TIDAL_DISSIPATION_KO.md', root/'notes/REQUEST275_PHOTOSPHERE_DENSITY_UNCERTAINTY_KO.md']
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
SCRIPTS = ['bgkeys274.sh', 'phase274-tides.py', 'run274.sh', 'atmos274.py', 'runatm274.sh', 'phase275-photosphere.py', 'run275.sh', 'phase275-publish.py']
FILES = {'phase274/tides.json': R4/'.phase274-tides.json', 'phase274/tides.log': R4/'.phase274-tides.log', 'phase274/atmos.json': R4/'.phase274-atmos.json',
         'phase274/atmos.log': R4/'.phase274-atmos.log', 'phase275/photosphere.json': R4/'.phase275-photosphere.json', 'phase275/photosphere.log': R4/'.phase275-photosphere.log'}


def package():
    assert not out.exists() and not manifest.exists()
    for src in FILES.values(): assert src.is_file(), src
    for name in SCRIPTS: assert (S/name).is_file(), name
    preserved = {p: h for p, h in read(previous_manifest)['sha256'].items() if not p.startswith('docs/')}
    for p, h in preserved.items(): assert sha(root/p) == h, p
    t = read(FILES['phase274/tides.json']); at = read(FILES['phase274/atmos.json']); ph = read(FILES['phase275/photosphere.json'])
    assert abs(t['love']['n1_polytrope_check'] - 0.2599088773) < 1e-6 and t['structural']['uniform_gamma_check_mass_weighted'] < 1e-3
    assert ph['sign_kept_everywhere'] and t['j0337']['margin'] > 1e6 and t['drives']['2 n_in (non-rotating tide)']['tau_wave_thick_envelope'] > 1e3
    out.mkdir(parents=True)
    for dst, src in FILES.items(): copy(src, out/dst)
    for name in SCRIPTS: copy(S/name, out/'scripts'/name)
    now = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    tau_at = {int(round(row['depth_km'])): row['gray_tau'] for row in at['rows']}
    kap = ph['cases']['kramers (a=1/2, b=1.25)']; con = ph['cases']['constant opacity (a=1, b=-1)']
    final = dict(classification='Counterexample candidate', passed=True,
        verdict='NO_GO_EXTENDED_TO_TIDAL_AND_DISSIPATIVE_CHANNELS__CHARGE_MAGNITUDE_SET_BY_PHOTOSPHERIC_GRAVITY',
        tidal=dict(selection_rule='Proven: the monopole compact charge responds at linear order only to the l=0 part of the drive (static spherical background)',
                   k2_apsidal=t['love']['k2_apsidal'], love_number=t['love']['love_number'], central_to_mean_density=t['love']['central_to_mean_density'],
                   DeltaPi_l1_s=t['buoyancy']['DeltaPi_l1_s'], DeltaPi_l2_s=t['buoyancy']['DeltaPi_l2_s'],
                   drives={k: dict(gmode_order_l2=v['gmode_order_l2'], tau_wave_thick_envelope=v['tau_wave_thick_envelope']) for k, v in t['drives'].items()},
                   nonstatic_inphase_tide=t['j0337']['nonstatic_inphase_tide'], j0337_tides=t['j0337']['tides']),
        dissipation=dict(charge_layer_gray_tau={'240 km': tau_at[240], '362 km': tau_at[362], '466 km': tau_at[466]}, tau_th_charge_layers_s=t['thermal']['tau_th_charge_layers_s'],
                         tau_th_envelope_base_s=t['thermal']['tau_th_base_s'], dq_per_eps=t['structural']['dq_per_eps'], kappa_struct_max=t['structural']['kappa_struct_max'],
                         kappa_over_beta=t['structural']['over_direct_beta'], lagged_monopole_delta_bound=t['j0337']['lagged_monopole_delta_bound'],
                         paper_b_lag_limit=t['j0337']['paper_b_lag_limit'], margin=t['j0337']['margin'], second_order_monopole_periodic=t['j0337']['second_order_monopole_periodic']),
        photosphere=dict(declared=ph['declared'], observed_Kaplan2014=ph['observed'], kaplan_central_ratio=[kap['Kaplan central']['ratio'], con['Kaplan central']['ratio']],
                         one_sigma_ratio_range=ph['one_sigma_ratio_range'], two_sigma_ratio_range=[min(kap['log g -2 sigma']['ratio'], con['log g -2 sigma']['ratio']), max(kap['log g +2 sigma']['ratio'], con['log g +2 sigma']['ratio'])],
                         sign_kept_everywhere=ph['sign_kept_everywhere'], face_background_kernel_weighted_deviation=ph['face_background_kernel_weighted_deviation']),
        assumptions=['linear response', 'white-dwarf spin period longer than about 45 min', 'WKB and diffusion approximations, ideal-gas heat capacity (tau_w is a lower bound)',
                     'single-relaxation quadrature bound', 'gray atmosphere, two opacity limits, fixed pulse kernel, radiation pressure and mu changes neglected'],
        final_charge_conclusion='conditional negative charge: sign robust to envelope structure and to the observed photosphere within 2 sigma; magnitude set by the photospheric density (x3.4-3.7 for the observed log g, 1 sigma x1.6-7.7); the free-fall memory collapses to static coefficients at orbital timescales in the radial, tidal and dissipative channels',
        full_goal_complete=False, snapshot_KST=now)
    doc = (f"\n\n## 단계274–275 — 조석·소산 경계와 광구 밀도 불확실성\n\n"
        "분류: Counterexample candidate. 정적 구대칭 배경에서 compact 전하는 선형 차수로 구동의 l=0 성분에만 반응한다(Proven 선택 규칙). "
        f"선언 배경은 중심 밀도가 평균의 {t['love']['central_to_mean_density']:.0f}배라 조석 근점 상수가 k₂={t['love']['k2_apsidal']:.2e}이고, J0337 내측 궤도에서 정적 조석의 상대 힘은 {t['j0337']['tides']['fractional_tidal_force']:.1e}다. "
        f"J0337 조석 구동은 l=2 g모드 차수 약 {t['drives']['2 n_in (non-rotating tide)']['gmode_order_l2']:.0f}에 해당하며, 광학적으로 두꺼운 봉투만으로 잰 감쇠 깊이가 {t['drives']['2 n_in (non-rotating tide)']['tau_wave_thick_envelope']:.1e}이라 이산 공명 없는 진행파 영역이다. "
        f"자유낙하 기억이 궤도 시간척도에서 붕괴해 들어가는 단극 정적 구조 계수는 |κ_struct|≤{t['structural']['kappa_struct_max']:.1e}(직접 계수 β의 {t['structural']['over_direct_beta']:.1e}, φ_∞=1e−3)이고, "
        f"이를 지연 상한으로 써도 J0337 SEP 진동은 {t['j0337']['lagged_monopole_delta_bound']:.1e}로 Paper B 한계 1.7e−9보다 {t['j0337']['margin']:.1e}배 작다. 따라서 no-go 경계를 조석·소산 채널로 확장한다. "
        f"전하 층은 광구다(회색 광학깊이 240 km {tau_at[240]:.3f}, 362 km {tau_at[362]:.2f}, 466 km {tau_at[466]:.1f}; 열 시간 1초 이하). "
        f"Kaplan et al. 2014의 log g 5.82±0.05, T_eff 15,800±100 K로 회색 대기를 정역학 상사 변환하면 전하 크기는 선언 대기의 {min(kap['Kaplan central']['ratio'], con['Kaplan central']['ratio']):.1f}–{max(kap['Kaplan central']['ratio'], con['Kaplan central']['ratio']):.1f}배(1σ {ph['one_sigma_ratio_range'][0]:.1f}–{ph['one_sigma_ratio_range'][1]:.1f}배)이고 부호는 2σ 범위에서 유지된다. "
        "전하 크기의 지배 불확실성은 표면중력이 정하는 광구 밀도다. "
        "[근거 274](../notes/REQUEST274_TIDAL_DISSIPATION_KO.md), [근거 275](../notes/REQUEST275_PHOTOSPHERE_DENSITY_UNCERTAINTY_KO.md).\n")
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
    m['native_tidal_photosphere'] = {k_: v for k_, v in final.items() if k_ != 'sha256'}
    master.write_text(json.dumps(m, indent=indent, ensure_ascii=False) + '\n', encoding='utf-8', newline='\n')


def git(*args): return subprocess.check_output(['git', *args], cwd=root)


def check(mode):
    m = read(manifest); a = read(master)
    for p, h in m['sha256'].items(): assert sha(root/p) == h and a['sha256'][p] == h, p
    for p, h in read(out/'publication.json')['previous_nondoc'].items(): assert sha(root/p) == h, p
    for p, v in m['document_prefixes'].items(): assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest() == v['sha256'], p
    assert a['sha256'][manifest.relative_to(root).as_posix()] == sha(manifest)
    paths = list(m['sha256']) + [manifest.relative_to(root).as_posix(), 'paper/revision-manifest.json']
    (S/'phase275-paths.txt').write_bytes(b'\0'.join(p.encode() for p in paths) + b'\0')
    if mode in ['staged', 'head']:
        if mode == 'staged':
            staged = set(filter(None, git('-c', 'core.quotepath=off', 'diff', '--cached', '--name-only', '-z').decode().split('\0')))
            assert staged == set(paths), sorted(staged ^ set(paths))[:10]
        ref = ':' if mode == 'staged' else 'HEAD:'
        for p in paths: assert hashlib.sha256(git('show', ref + p)).hexdigest() == sha(root/p), p
    print(json.dumps(dict(mode=mode, bound_files=len(m['sha256']), paths=len(paths), verdict=m['verdict'], k2=m['tidal']['k2_apsidal'],
                          kappa=m['dissipation']['kappa_struct_max'], margin=m['dissipation']['margin'], one_sigma=m['photosphere']['one_sigma_ratio_range'])))


if __name__ == '__main__':
    if sys.argv[1] == 'package': package()
    check(sys.argv[1] if sys.argv[1] != 'package' else 'check')
