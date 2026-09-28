"""Archive phase 276: final-charge synthesis and J0337 observational closure of the free-fall compact charge.

Counterexample candidate / no-go boundary. Copies the closure record and scripts into outputs/direct-eos-gr33/native-final-closure,
binds the synthesis note and the English section draft, appends one section (with the failure-ledger entry) to the six
dynamic-chi documents, and binds everything in the manifest and paper/revision-manifest.json.
Usage: python phase276-publish.py package|check|staged|head
"""
from pathlib import Path
import datetime, hashlib, json, subprocess, sys
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
S = Path(__file__).resolve().parent
master = root/'paper/revision-manifest.json'; previous_manifest = root/'outputs/direct-eos-gr33/native-tidal-photosphere-manifest.json'
out = root/'outputs/direct-eos-gr33/native-final-closure'; manifest = out.parent/'native-final-closure-manifest.json'
notes = [root/'notes/REQUEST276_FINAL_CHARGE_CLOSURE_KO.md', root/'docs/white-dwarf-free-fall-charge-section.md']
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
SCRIPTS = ['phase276-closure.py', 'phase276-publish.py']


def package():
    assert not out.exists() and not manifest.exists()
    for name in SCRIPTS + ['phase276-closure.json']: assert (S/name).is_file(), name
    preserved = {p: h for p, h in read(previous_manifest)['sha256'].items() if not p.startswith('docs/')}
    for p, h in preserved.items(): assert sha(root/p) == h, p
    c = read(S/'phase276-closure.json'); cl = c['closure']; acc = c['accepted']; ob = c['observed_gravity_charge']
    assert c['symbolic_freefall_residual'] == '0' and cl['freefall_delta_bound_cassini'] < 1e-15 and cl['full_beta_relaxing_delta_bound'] < cl['paper_b_lag_limit']
    out.mkdir(parents=True)
    copy(S/'phase276-closure.json', out/'closure.json')
    for name in SCRIPTS: copy(S/name, out/'scripts'/name)
    now = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    final = dict(classification='Counterexample candidate at pulse timescales; no-go boundary (theorem progress) at orbital timescales', passed=True,
        verdict='FINAL_CHARGE_CONCLUSION_KEPT__NEGATIVE_FREE_FALL_MEMORY__NO_J0337_OBSERVABLE__A4_HOLDS_FOR_THIS_STATE',
        final_charge=dict(continuum=acc['richardson'], per_eta_phi_inf_squared=c['normalization']['per_eta_phi2'], observed_gravity_kaplan_central=ob['kaplan_central'],
                          observed_gravity_one_sigma=ob['one_sigma'], observed_gravity_two_sigma=ob['two_sigma'], sign_kept_all_accepted_variations=True,
                          normalization_scaling='eta * phi_inf^2 (derived: force and source coupling both proportional to alpha(phi0); not rerun at another phi_inf)'),
        closure=cl, accepted=acc,
        failure_ledger=dict(failing_step='omega << omega0 = 0.035 1/s, omega << omega_ac, g-mode traveling-wave regime: the free-fall memory becomes the static structural response |kappa_struct| <= 8.8e-9 (2.2e-9 of beta)',
                            minimal_missing_assumption='an internal state with relaxation time ~ orbital period coupled to the inner white-dwarf monopole charge with |kappa_lag| >~ 1.6e3/|alpha_p|; a weak-field white dwarf offers |beta| = 4 (instantaneous) plus a ~1e-8 structural term'),
        remaining=['one GR return; no ADM or infinity charge normalization', 'no fully nonlinear run and no rerun at another phi_inf',
                   'endpoint protocol T = 2D; long-time relaxation over 28-178 s not computed', 'gray atmosphere, two opacity limits, fixed pulse kernel',
                   'spin period >~ 45 min, WKB and diffusion approximations', 'manuscript integration pending Pandoc/TeX'],
        full_goal_complete=False, snapshot_KST=now)
    doc = ("\n\n## 단계276 — 최종 전하 판정과 관측 폐쇄\n\n"
        f"분류: Counterexample candidate. 최종 전하 결론은 유지된다. 지배 수치 오차를 해결한 같은 결합 해의 끝점 compact 전하는 연속 극한 {acc['richardson']:.4e}(φ_∞=1e−3, η=1e−30; q/(ηφ_∞²)={c['normalization']['per_eta_phi2']:.3e})로 음수다. "
        f"부호는 EOS·광학, 봉투 밀도 분포, 관측 광구(2σ)에 견고하다. 크기는 관측 표면중력에서 {ob['kaplan_central'][0]:.2e}–{ob['kaplan_central'][1]:.2e}(1σ {ob['one_sigma'][0]:.1e}–{ob['one_sigma'][1]:.1e})로 수정된다. 자유낙하 닫힌 식의 기호 잔차는 0이다. "
        f"분류: Conjectural. 관측 폐쇄: Cassini 2σ로 |α₀|≤{cl['cassini_alpha0_max_2sigma']:.2e}이면, 내측 백색왜성 단극 전하의 선형·수동 지연이 J0337 지연 한계 1.7e−9에 닿으려면 |κ_lag α_p|≥{cl['required_kappa_lag_times_alpha_p']:.2e}여야 한다. "
        f"정적 감수율 전체(|β|=4)가 완화되어도 {cl['full_beta_relaxing_delta_bound']:.1e}|α_p|, 자유낙하 기작은 {cl['freefall_delta_bound_cassini']:.1e}|α_p|다. 따라서 이 감도의 지연 신호는 백색왜성이 아니라 중성자별 전하 쪽이어야 한다. "
        "분류: Conjectural. 미션 분류는 theorem progress(no-go 경계)이며, 이 최소 상태에 대해 A4가 유지된다. "
        "실패 원장: 정확한 붕괴 단계는 ω≪ω₀, ω≪ω_ac, g모드 진행파 영역에서 자유낙하 기억이 정적 구조 응답(|κ_struct|≤8.8e−9)으로 바뀌는 단계다. 최소 누락 가정은 궤도 주기 완화 시간과 |κ_lag|≳1.6e3/|α_p|의 결합을 가진 내부 상태이며, 약한장 백색왜성에서는 성립하지 않는다. "
        "남은 조건: GR 반환 1회(ADM·무한대 정규화 없음), φ_∞² 스케일은 유도, 끝점 프로토콜 T=2D, 회색 대기, 느린 자전, 원고 통합은 Pandoc·TeX 부재로 보류. "
        "[근거 276](../notes/REQUEST276_FINAL_CHARGE_CLOSURE_KO.md), [영문 절 초안](white-dwarf-free-fall-charge-section.md).\n")
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
    m['native_final_closure'] = {k_: v for k_, v in final.items() if k_ != 'sha256'}
    master.write_text(json.dumps(m, indent=indent, ensure_ascii=False) + '\n', encoding='utf-8', newline='\n')


def git(*args): return subprocess.check_output(['git', *args], cwd=root)


def check(mode):
    m = read(manifest); a = read(master)
    for p, h in m['sha256'].items(): assert sha(root/p) == h and a['sha256'][p] == h, p
    for p, h in read(out/'publication.json')['previous_nondoc'].items(): assert sha(root/p) == h, p
    for p, v in m['document_prefixes'].items(): assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest() == v['sha256'], p
    assert a['sha256'][manifest.relative_to(root).as_posix()] == sha(manifest)
    paths = list(m['sha256']) + [manifest.relative_to(root).as_posix(), 'paper/revision-manifest.json']
    (S/'phase276-paths.txt').write_bytes(b'\0'.join(p.encode() for p in paths) + b'\0')
    if mode in ['staged', 'head']:
        if mode == 'staged':
            staged = set(filter(None, git('-c', 'core.quotepath=off', 'diff', '--cached', '--name-only', '-z').decode().split('\0')))
            assert staged == set(paths), sorted(staged ^ set(paths))[:10]
        ref = ':' if mode == 'staged' else 'HEAD:'
        for p in paths: assert hashlib.sha256(git('show', ref + p)).hexdigest() == sha(root/p), p
    print(json.dumps(dict(mode=mode, bound_files=len(m['sha256']), paths=len(paths), verdict=m['verdict'], continuum=m['final_charge']['continuum'],
                          required=m['closure']['required_kappa_lag_times_alpha_p'], freefall_bound=m['closure']['freefall_delta_bound_cassini'])))


if __name__ == '__main__':
    if sys.argv[1] == 'package': package()
    check(sys.argv[1] if sys.argv[1] != 'package' else 'check')
