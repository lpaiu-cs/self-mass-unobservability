"""Archive phase 287: thermal-relaxation strength of the white dwarf's structural monopole response at orbital frequencies
(per-layer relaxation bracket with isothermal relaxed exponents), plus analytic notes on the non-common timing channels and
the Section 4.3 stiffness condition. The manuscript is not changed.
Usage: python phase287-publish.py package|check|staged|head
"""
from pathlib import Path
import datetime, hashlib, json, subprocess, sys
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
S = Path(__file__).resolve().parent
master = root/'paper/revision-manifest.json'; previous_manifest = root/'outputs/direct-eos-gr33/native-manuscript-reintegration-manifest.json'
out = root/'outputs/direct-eos-gr33/native-thermal-relaxation'; manifest = out.parent/'native-thermal-relaxation-manifest.json'
notes = [root/'notes/REQUEST287_THERMAL_RELAXATION_STRENGTH_KO.md']
DOCS = ['model-definition', 'observable-targets', 'adiabatic-limit', 'nonadiabatic-regime', 'failure-ledger-dynamic-chi', 'dynamic-charge-completion']
FILES = {'phase287/relax.json': '.phase287-relax.json', 'phase287/relax.log': '.phase287-relax.log',
         'scripts/phase287-relax.py': 'phase287-relax.py', 'scripts/run287.sh': 'run287.sh', 'scripts/phase287-publish.py': 'phase287-publish.py'}
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
def git(*args): return subprocess.check_output(['git', '-c', 'gc.auto=0', '-c', 'maintenance.auto=false', *args], cwd=root)


def package():
    assert not out.exists() and not manifest.exists()
    for name in FILES.values(): assert (S/name).is_file(), name
    assert (S/'.phase287-relax.log').read_text(encoding='utf-8').rstrip().endswith('exit 0')
    preserved = {p: h for p, h in read(previous_manifest)['sha256'].items() if not p.startswith('docs/')}
    for p, h in preserved.items(): assert sha(root/p) == h, p
    for dst, name in FILES.items(): copy(S/name, out/dst)
    r = read(out/'phase287/relax.json')
    assert abs(r['dq_eps_adiabatic']/-1.840953484516629e-07 - 1) < 1e-9 and not r['stopped_by_positivity']
    lag = {k: dict(debye_quadrature_over_S_struct=r[k]['debye_quadrature_rel'], window_TV_half=r[k]['shells_within_1e4_of_1_over_omega_TV_half'],
                   depth_km_tau_1_over_omega=r[k]['depth_km_at_tau_1_over_omega'], mass_fraction_above=r[k]['mass_fraction_above'],
                   deeper_than_scan_suppression=r[k]['deeper_than_scan_suppression']) for k in ('n_in', '2 n_in', 'n_out')}
    assert max(v['debye_quadrature_over_S_struct'] for v in lag.values()) < 1e-2
    now = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    final = dict(classification='Imported from prior work (computation); Conjectural (per-layer relaxation model)', passed=True,
        verdict='THERMAL_RELAXATION_LAG_COMPUTED__4E-9_INNER_3E-7_OUTER_OF_S_STRUCT__ACCEPTANCE_1PCT_MET',
        acceptance='Debye-weighted quadrature <= 1e-2 |S_struct| replaces the strength assumption within the per-layer model; >= 1 would be loophole progress',
        model='overlying thermal time tau_th = int c_p T dm / L (ideal monatomic, no ionization energy); relaxed limit bracketed by isothermal Gamma_T = P_gas/P; 161 cuts 1e-3..1e13 s; positive-definite operator throughout',
        scanned=dict(last_cut_s=r['last_cut_s'], depth_km=r['last_cut_depth_km'], mass_fraction=r['last_cut_mass_fraction'], dS_rel=r['dS_rel_at_last_cut'], total_variation_rel=r['total_variation_rel']),
        lag=lag, lagged_pair_factor_bound='<= 3.3e-7 x 2.13e-18 ~ 7e-25',
        analytic_notes=['non-common channels enter timing through the same inner-coupling and SEP-like differential-acceleration paths as the common template; O(1)-O(10) weights expected; not computed (Nutimo work excluded by AGENTS.md)',
                        'Section 4.3 stiffness shift with weak-field kappa ~ 1/(|beta_p| m_p): |beta_s| m_wd m_p |beta_p| / a_in^2 ~ 1.1e-13 |beta_p|'],
        repository_classification='theorem progress: computed quantitative boundary for this state within the per-layer relaxation model; no observational exclusion; not a general A4',
        manuscript_changed=False, full_goal_complete=False, snapshot_KST=now)
    doc = ("\n\n## 단계287 — 깊은 층 열 완화 세기의 판정\n\n"
        "분류: Imported from prior work. 사용자 승인(2026-09-28, 부분 진행)으로 열 완화 세기를 싼 계산으로 판정했다. 모형은 층별 열 완화다: 각 층의 완화 시간은 위쪽 층의 열 시간이고, 완화 극한은 등온 Γ_T=P_gas/P로 괄호를 쳤다. "
        "단계274 정적 풀이기를 1e−3–1e13 s의 절단 161개로 다시 풀었다. 연산자는 모두 양정치였고, 단열 dq/ε는 그대로 재현됐다. "
        "궤도 진동수의 Debye 가중 지연은 |𝒮_struct|의 4.0e−9(내측 궤도)와 3.3e−7(외측 궤도)이고, 창 상한으로도 ≤8.9e−5다. 수락 기준(1%)을 충족한다. "
        "τ_th=1/ω인 층은 깊이 2,600–6,600 km, 위쪽 질량 몫 1e−9–1e−7로 가볍다. 비공통 채널 타이밍과 §4.3 강성 조건은 해석 추정만 남겼다(강성 비 약 1e−13|β_p|). 원고는 고치지 않았다.\n\n"
        "분류: Conjectural. 이 상태의 궤도 지연 결합은 층별 완화 모형에서 단열 구조 척도의 ≲3.3e−7이고, 쌍 인자로는 ≲7e−25다. "
        "분류는 이 상태에 대한 계산된 정량 경계(theorem progress)다. 관측 배제나 일반 A4는 아니다. 남은 최소 계산은 두 가지다: 완전한 비단열 반경 확산 응답, 비공통 채널의 타이밍 응답. "
        "[근거 287](../notes/REQUEST287_THERMAL_RELAXATION_STRENGTH_KO.md).\n")
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
    m['native_thermal_relaxation'] = {k: v for k, v in final.items() if k not in ('sha256', 'document_prefixes')}
    master.write_text(json.dumps(m, indent=indent, ensure_ascii=False) + '\n', encoding='utf-8', newline='\n')


def check(mode):
    m = read(manifest); a = read(master)
    for p, h in m['sha256'].items(): assert sha(root/p) == h and a['sha256'][p] == h, p
    for p, h in read(out/'publication.json')['previous_nondoc'].items(): assert sha(root/p) == h, p
    for p, v in m['document_prefixes'].items(): assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest() == v['sha256'], p
    assert a['sha256'][manifest.relative_to(root).as_posix()] == sha(manifest)
    paths = list(m['sha256']) + [manifest.relative_to(root).as_posix(), 'paper/revision-manifest.json']
    (S/'phase287-paths.txt').write_bytes(b'\0'.join(p.encode() for p in paths) + b'\0')
    if mode in ['staged', 'head']:
        if mode == 'staged':
            staged = set(filter(None, git('-c', 'core.quotepath=off', 'diff', '--cached', '--name-only', '-z').decode().split('\0')))
            assert staged <= set(paths), sorted(staged - set(paths))[:10]
        ref = ':' if mode == 'staged' else 'HEAD:'
        for p in paths: assert hashlib.sha256(git('show', ref + p)).hexdigest() == sha(root/p), p
    print(json.dumps(dict(mode=mode, bound_files=len(m['sha256']), paths=len(paths), verdict=m['verdict'])))


if __name__ == '__main__':
    if sys.argv[1] == 'package': package()
    check(sys.argv[1] if sys.argv[1] != 'package' else 'check')
