"""Archive phase 280: closure of the open conditions of the phase-276 synthesis and the revised final charge statement.

Binds the synthesis note and the revised English section draft, appends one section to the six dynamic-chi documents, and binds
everything in the manifest and paper/revision-manifest.json. No new computation: values are read from the phase 276-279 manifests.
Usage: python phase280-publish.py package|check|staged|head
"""
from pathlib import Path
import datetime, hashlib, json, subprocess, sys
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
S = Path(__file__).resolve().parent
master = root/'paper/revision-manifest.json'; previous_manifest = root/'outputs/direct-eos-gr33/native-atmosphere-reconstruction-manifest.json'
out = root/'outputs/direct-eos-gr33/native-conditions-closed'; manifest = out.parent/'native-conditions-closed-manifest.json'
notes = [root/'notes/REQUEST280_OPEN_CONDITIONS_CLOSED_KO.md', root/'docs/white-dwarf-free-fall-charge-section.md']
DOCS = ['model-definition', 'observable-targets', 'adiabatic-limit', 'nonadiabatic-regime', 'failure-ledger-dynamic-chi', 'dynamic-charge-completion']
SCRIPTS = ['phase280-publish.py']
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


def package():
    assert not out.exists() and not manifest.exists()
    for name in SCRIPTS: assert (S/name).is_file(), name
    preserved = {p: h for p, h in read(previous_manifest)['sha256'].items() if not p.startswith('docs/')}
    for p, h in preserved.items(): assert sha(root/p) == h, p
    O = root/'outputs/direct-eos-gr33'
    b = read(O/'native-closure-transit-manifest.json'); a = read(O/'native-atmosphere-reconstruction-manifest.json'); f = read(O/'native-final-closure-manifest.json')
    out.mkdir(parents=True)
    for name in SCRIPTS: copy(S/name, out/'scripts'/name)
    now = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    qc = b['closure']['continuum_mass_normalized']; kg, kn = a['kaplan_central']['gray'], a['kaplan_central']['nongray']
    final = dict(classification='Counterexample candidate (endpoint and history); Proven (permanent charge zero for positive damping); theorem progress (no-go)', passed=True,
        verdict='OPEN_CONDITIONS_CLOSED__ENDPOINT_CHARGE_NEGATIVE_TRANSIENT__PERMANENT_ZERO__A4_HOLDS',
        final_charge=dict(endpoint_continuum_mass_normalized=qc, per_eta_phi_inf2=qc/1e-36, observed_gray=qc*kg, observed_nongray_lte=qc*kn,
                          observed_one_sigma_gray=[qc*x for x in a['one_sigma']['gray']], observed_one_sigma_nongray=[qc*x for x in a['one_sigma']['nongray']],
                          history=dict(transit_max_abs=b['history']['transit_max_abs'], relaxation_rms=b['history']['relaxation_rms'],
                                       permanent='zero at first order in eta for positive mode damping')),
        closed=dict(gr_fixed_point=b['closure']['fixed_point_residual'], adm=b['closure']['adm_cross'], infinity=b['closure']['infinity_tail_bound'],
                    nonlinear=b['closure']['nonlinear'], phi_inf=b['closure']['phi_inf'], endpoint_and_long_time='phase 278',
                    atmosphere=dict(validation=a['gray']['validation_kernel_weighted'], nongray_flux=a['nongray']['flux_error_tau_le_100'])),
        j0337_closure=f['closure'],
        remaining=['positive radial-mode damping assumed', 'interior Born scattering neglected', 'non-LTE and line/metal opacity in the non-gray correction',
                   'phi_inf exponent of the 1e-8 direct part', 'slow white-dwarf rotation', 'manuscript integration pending Pandoc (and TeX for the PDF)'],
        full_goal_complete=False, snapshot_KST=now)
    doc = ("\n\n## 단계280 — 남은 조건의 폐쇄와 최종 전하 진술의 개정\n\n"
        f"분류: Counterexample candidate. 끝점 T=2D의 compact 전하는 음수다. 질량 정규화 연속 극한은 {qc:.5e}이고, 관측 광구에서는 회색 {qc*kg:.2e}, LTE 비회색 {qc*kn:.2e}다. "
        f"다만 이 값은 입사 펄스의 단극 산란 신호 가운데 가장 이른 광구 부분이다. 신호는 들어오는 구간에서 음수로 자라고, 나가는 통과에서 최대 {b['history']['transit_max_abs']:.1e}로 진동하며, 펄스가 떠난 뒤 영평균 반경 모드 진동이 된다. 영구 전하는 η의 1차에서 0이다(양의 감쇠, Proven). "
        f"단계276에 남은 조건을 모두 닫았다: GR 고정점 {b['closure']['fixed_point_residual']:.1e}, ADM {b['closure']['adm_cross']:.1e}, 무한대 {b['closure']['infinity_tail_bound']:.1e}, 비선형 {b['closure']['nonlinear']:.1e}, φ_∞ 짝함수·φ_∞² 가중(단계277), 전 시간 이력(단계278), 관측 광구 재구성과 비회색 보정(단계279). "
        "남는 가정은 양의 모드 감쇠, 내부 Born 산란 무시, 비회색 보정의 LTE 연속 불투명도, 직접 결합 부분의 φ_∞ 지수, 느린 자전이다. 원고 통합은 Pandoc이 없어 보류했다. no-go 경계와 A4 유지 판정은 바뀌지 않는다. "
        "[근거 280](../notes/REQUEST280_OPEN_CONDITIONS_CLOSED_KO.md), [개정 영문 절 초안](white-dwarf-free-fall-charge-section.md).\n")
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
    m['native_conditions_closed'] = {k_: v for k_, v in final.items() if k_ != 'sha256'}
    master.write_text(json.dumps(m, indent=indent, ensure_ascii=False) + '\n', encoding='utf-8', newline='\n')


def git(*args): return subprocess.check_output(['git', *args], cwd=root)


def check(mode):
    m = read(manifest); a = read(master)
    for p, h in m['sha256'].items(): assert sha(root/p) == h and a['sha256'][p] == h, p
    for p, h in read(out/'publication.json')['previous_nondoc'].items(): assert sha(root/p) == h, p
    for p, v in m['document_prefixes'].items(): assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest() == v['sha256'], p
    assert a['sha256'][manifest.relative_to(root).as_posix()] == sha(manifest)
    paths = list(m['sha256']) + [manifest.relative_to(root).as_posix(), 'paper/revision-manifest.json']
    (S/'phase280-paths.txt').write_bytes(b'\0'.join(p.encode() for p in paths) + b'\0')
    if mode in ['staged', 'head']:
        if mode == 'staged':
            staged = set(filter(None, git('-c', 'core.quotepath=off', 'diff', '--cached', '--name-only', '-z').decode().split('\0')))
            assert staged <= set(paths), sorted(staged - set(paths))[:10]
        ref = ':' if mode == 'staged' else 'HEAD:'
        for p in paths: assert hashlib.sha256(git('show', ref + p)).hexdigest() == sha(root/p), p
    print(json.dumps(dict(mode=mode, bound_files=len(m['sha256']), paths=len(paths), verdict=m['verdict'], endpoint=m['final_charge']['endpoint_continuum_mass_normalized'],
                          observed=[m['final_charge']['observed_gray'], m['final_charge']['observed_nongray_lte']])))


if __name__ == '__main__':
    if sys.argv[1] == 'package': package()
    check(sys.argv[1] if sys.argv[1] != 'package' else 'check')
