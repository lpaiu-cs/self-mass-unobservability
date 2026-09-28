"""Archive phase 292: final independent review of the focused paper (Claude Opus 5.5, Claude Fable 5.1, GPT-6-Astra), the periastron
convention correction (ELL1: eta = e sin varpi) with regenerated outputs, and the manuscript, SM and submission revisions.
Usage: python phase292-publish.py package|check|staged|head
"""
from pathlib import Path
import csv, datetime, hashlib, json, math, re, subprocess, sys
import pypdfium2 as pdfium
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
S = Path(__file__).resolve().parent
master = root/'paper/revision-manifest.json'; previous_manifest = root/'outputs/direct-eos-gr33/native-focused-manuscript-manifest.json'
out = root/'outputs/direct-eos-gr33/native-focused-final-review'; manifest = out.parent/'native-focused-final-review-manifest.json'
RC = 'outputs/research-completion/'
REGENERATED = [RC + n for n in ('physical-matching.json', 'corrected-physical-drive.json', 'comparator-audit.json', 'simultaneous-validation.json', 'runtime12-analysis.json')] + \
              ['paper/figures/comparator-phase-validation.pdf', 'paper/figures/remaining-levers-validation.pdf']
CHANGED = [root/p for p in REGENERATED + [
    'symbolic/physical_matching.py', 'verification/physical_drive_completion.py', 'verification/verify_unified_paper.py',
    'paper/manuscript.md', 'paper/main.tex', 'paper/supplement.md', 'paper/supplement.tex', 'paper/references.bib', 'paper/README.md',
    'docs/white-dwarf-free-fall-charge-section.md', 'output/pdf/free-fall-identifiability.pdf', 'output/pdf/free-fall-identifiability-supplement.pdf',
    'output/submission/free-fall-identifiability-source.zip', 'output/submission/cover-letter-prd.tex', 'output/submission/cover-letter-prd.pdf',
    'output/submission/abstract-plain.txt', 'output/submission/submission-checklist-prd.md']]
WITHDRAWN = root/RC/'withdrawn-periastron-convention'
NEW = [root/'notes/REQUEST292_FINAL_INDEPENDENT_REVIEW_KO.md'] + sorted(p for p in WITHDRAWN.iterdir() if p.is_file())
DOCS = ['model-definition', 'observable-targets', 'adiabatic-limit', 'nonadiabatic-regime', 'failure-ledger-dynamic-chi', 'dynamic-charge-completion',
        'remaining-levers-2026-09-09']
SCRIPTS = ['compare292.py', 'phase292-validation-data.py', 'phase292-supplement.py', 'phase292-draft.py', 'find292.py', 'ransom-table.py',
           'repro292.py', 'phase292-publish.py']
REVIEWS = {'final292-prompt.md': 'review-prompt.md', 'final292-astra-prompt.md': 'review-gpt-6-astra-prompt.md',
           'final292-agent-prompt.md': 'review-claude-prompt.md', 'final292-astra.md': 'review-gpt-6-astra.md',
           'final292-astra.log': 'review-gpt-6-astra-cli.log', 'final292-opus.md': 'review-opus-5.5.md', 'final292-fable.md': 'review-fable-5.1.md',
           'compare292.json': 'compare292.json'}
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
pages = lambda p: len(pdfium.PdfDocument(str(p)))


def provenance(ms):
    """Every number the paper adds must match its source at the printed rounding."""
    m = read(root/RC/'physical-matching.json'); c = read(root/RC/'comparator-audit.json'); v = read(root/RC/'simultaneous-validation.json')
    live = read(root/RC/'runtime12-analysis.json')
    want = [f"{math.degrees(m['pericenter_radians'][k]):.2f}" for k in ('in', 'out')] + [f"{m['physical_closure_radians'] % (2*math.pi):.8f}"]
    for k, x in m['pericenter_radians'].items(): assert f'{x:.8f}' in ms, k
    rows = list(csv.DictReader(open(root/RC/'nuisance-intervals.csv', encoding='utf-8')))
    full = {float(r['tau']): float(r['U']) for r in rows if float(r['cut']) == 0 and r['rank'] == '90' and float(r['K']) == 1}
    trunc = {float(r['tau']): float(r['U']) for r in rows if r['cut'] == '0.001' and r['rank'] == '71' and float(r['K']) == 1}
    for tau in (2., 5., 18., 52., 200.):
        line = next(l for l in ms.splitlines() if l.startswith(f'{tau:g} & '))
        assert f"{full[tau]:.3e}".split('e')[0] in line, tau
    k1 = [full[t]/trunc[t] for t in full]; want += [f'{min(k1):.2f}', f'{max(k1):.1f}']
    curve = list(csv.DictReader(open(root/'request10_external/sep_dynamic/sep_limit_curve_10_8e.tsv', encoding='utf-8'), delimiter='\t'))
    k10 = [float(r['u95pm_K10_fullrank'])/float(r['u95pm_K10']) for r in curve]; assert len(k10) == 65; want += [f'{min(k10):.2f}', f'{max(k10):.2f}']
    mins = {}
    for r in c['rows']: mins[r['comparator']] = min(mins.get(r['comparator'], 1.), r['relative_information'])
    want += [f"{mins['P1']:.6f}"] + [f"{mins[k]:.2e}".split('e')[0] for k in ('P2', 'P3', 'P4')] + [f"{round(mins['P4']**-.5)}", f"{mins['even_P4']:.4f}"]
    p5 = [r['relative_information'] for r in c['rows'] if r['comparator'] == 'P5']
    want += [f"{min(p5):.1e}".replace('e-', '\\times10^{-') + '}', f"{max(p5):.1e}".replace('e-', '\\times10^{-') + '}']
    sections = v['data']['physical_lag_sections']; assert all(s['empty'] and s['minimum_statistic'] > v['threshold'] for s in sections)
    other = [s['minimum_statistic'] for s in sections if s['tau'] != 2.]
    want += [f"{sections[0]['minimum_statistic']:.2f}", f'{min(other):.2f}', f'{max(other):.2f}']
    h = v['data']['omnibus_null_statistic']/2; want.append(f'{math.exp(-h)*(1+h+h*h/2):.3f}')
    ratios = [f['sigma_beta_ratio'] for t in live['transient'] for f in t['fits']]; want += [f'{min(ratios):.6f}', f'{max(ratios):.5f}']
    proj = [r['new_sigma_over_archived'] for r in live['physical_diagonal_projection_comparison']]; want += [f'{min(proj):.4f}', f'{max(proj):.4f}']
    missing = [w for w in want if w not in ms]; assert not missing, missing
    return want


def package():
    assert not out.exists() and not manifest.exists()
    for name in SCRIPTS + list(REVIEWS): assert (S/name).is_file(), name
    for p in NEW: assert p.is_file(), p
    preserved = {p: h for p, h in read(previous_manifest)['sha256'].items() if not p.startswith('docs/')}
    for p, h in preserved.items():
        if root/p not in CHANGED: assert sha(root/p) == h, p
    for p in CHANGED: assert hashlib.sha256(git('show', 'HEAD:' + p.relative_to(root).as_posix())).hexdigest() != sha(p), p
    # The withdrawn copies are the committed pre-correction outputs.
    for rel in REGENERATED:
        assert hashlib.sha256(git('show', 'HEAD:' + rel)).hexdigest() == sha(WITHDRAWN/Path(rel).name), rel
    ms = (root/'paper/manuscript.md').read_bytes().decode('utf-8')
    assert '\n' not in ms.replace('\r\n', '') and '**Repository:**' not in ms and 'deprojected' not in ms
    checked = provenance(ms)
    assert (pages(root/'output/pdf/free-fall-identifiability.pdf'), pages(root/'output/pdf/free-fall-identifiability-supplement.pdf'),
            pages(root/'output/submission/cover-letter-prd.pdf')) == (13, 30, 1)
    for name in ('main', 'supplement'):
        log = (root/f'paper/build/{name}.log').read_text(encoding='utf-8', errors='replace')
        assert 'Overfull' not in log and 'undefined' not in log.lower(), name
    for name in SCRIPTS: copy(S/name, out/'scripts'/name)
    for src, dst in REVIEWS.items(): copy(S/src, out/dst)
    for name in ('main', 'supplement'): copy(root/f'paper/build/{name}.log', out/f'build/{name}.log')
    now = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    final = dict(classification='Imported from prior work (review, correction and regenerated outputs); Proven (convention check)', passed=True,
        verdict='FINAL_REVIEW_DONE__PERIASTRON_CONVENTION_CORRECTED__J0337_PHYSICAL_PLANE_REJECTED__FINAL_CHARGE_CONCLUSION_UNCHANGED',
        user_instruction='final independent review by Claude Opus 5.5, Claude Fable 5.1 and GPT-6-Astra (2026-09-28)',
        reviews={'GPT-6-Astra': 'major revision; no mathematical error or numerical mismatch', 'Claude Fable 5.1': 'major revision; mathematics verified with sympy, numbers matched',
                 'Claude Opus 5.5': 'major revision; blocking error B1 (eta/kappa convention reversed relative to the released code)'},
        blocking_error=dict(description='physical-drive pericenters used varpi = atan2(kappa, eta); the released Nutimo code uses eta = e sin(varpi), kappa = e cos(varpi) (ELL1)',
                            evidence=['Nutimo src/Parameters.cpp: omp = inversetrigo(kappap/ei, etap/ei); src/Utilities.cpp: inversetrigo(cosv, sinv)',
                                      'frozen parameter indices 2/4/8/10 = eta_p/kappa_p/eta_b/kappa_b = Nutimo etap/kappap/etaB/kappaB',
                                      'Ransom et al. 2014 Table 1: (e sin w)_I = 6.8567e-4, (e cos w)_I = -9.171e-5, w_I = 97.6182 deg, w_O = 95.619493 deg',
                                      'time origin: Nutimo shifts treference to zero and builds the initial state there; DopplerF = 1'],
                            cause='the parameter ordering "e cos w, e sin w" in the text of Voisin et al. 2025 was read as (eta, kappa)',
                            failing_step='symbolic/physical_matching.py and verification/physical_drive_completion.py pericenter angle',
                            missing_assumption='check of the runtime code convention'),
        corrected=dict(pericenter_radians=[1.69406454, 1.67084424], closure=3.16481295, physical_minus_archived=[0.12326821, 0.10004792, 'pi'],
                       physical_sections='all six empty; minimum statistics 12.84 (2 d) and 14.39-14.59 above threshold 12.8242',
                       comparator_minima=[0.005389, 6.29e-5, 1.17e-5, 4.63e-6], even_only=0.0311, transient_sigma_ratio=[0.999998, 1.00742],
                       halfstep_sigma_factor=[0.7858, 0.8722]),
        withdrawn=['every physical lag section includes beta = 0', 'endpoints 4.10217e-10 and 8.31672e-9', 'minima 0.08289/8.36e-6/3.19e-6/2.88e-6 (589x), even 0.003035',
                   'transient ratio 1.00096-1.00523', 'half-step factors 0.8680-0.8785'],
        unchanged=['archived dictionary closure violation and physical-beta withdrawal', 'omitted-input RMS 3.37427 percent', 'Tables 1-2 envelopes and coverage',
                   'registered six-coefficient inclusion 95.17-95.43 percent', 'omnibus 16.3525', 'positive rate-gap witnesses', 'white-dwarf final-charge conclusion'],
        reproducibility='rerunning simultaneous_inference validate reproduces the registered simulations only within Monte Carlo error when the thread count differs (eigenbasis of repeated eigenvalues); a fixed configuration reproduces itself; registered rows kept, data section replaced',
        responses='notes/REQUEST292_FINAL_INDEPENDENT_REVIEW_KO.md', new_number_checks=checked,
        pages=dict(paper=13, supplement=30, cover_letter=1),
        remaining=['confirmation review of the correction', 'affiliation, e-mail, ORCID', 'public snapshot update (author instruction)', 'cause of the six-coefficient excess'],
        full_goal_complete=False, snapshot_KST=now)
    doc = ("\n\n## 단계292 — 최종 독립 검토와 근점 규약 정정\n\n"
        "분류: Imported from prior work. Claude Opus 5.5, Claude Fable 5.1, GPT-6-Astra의 독립 심사가 모두 주요 수정을 권했다. "
        "Opus가 J0337 물리 구동의 η·κ 규약 오류를 찾았다. 코드는 η=e sin ϖ를 쓰는데 분석은 η=e cos ϖ로 읽었다. "
        "이를 바로잡고 영향받는 계산을 다시 돌렸다(각 수 초). 이전 출력은 `outputs/research-completion/withdrawn-periastron-convention/`에 보존했다.\n\n"
        "실패 기록: 실패한 단계는 `symbolic/physical_matching.py`와 `verification/physical_drive_completion.py`의 근점 각 계산이다. "
        "빠진 최소 가정은 런타임 코드의 매개변수 규약 확인이다. 이전 결론 '모든 물리 지연 단면이 β=0을 포함한다'를 철회한다.\n\n"
        "분류: Imported from prior work. 수정 뒤 여섯 단면이 모두 비어 있다. 기록된 반송파 초과(옴니버스 16.3525, 명목 p≈0.012)를 "
        "물리 구동의 순간+완화 응답이 재현하지 못한다. relaxation 검출은 아니며 원인은 분리하지 못했다. 백색왜성의 최종 전하 결론은 유지된다. "
        "문헌 위치, 두 표본 조건, 한정어, 정의 등 나머지 지적도 반영했다(본문 13쪽). "
        "[근거 292](../notes/REQUEST292_FINAL_INDEPENDENT_REVIEW_KO.md).\n")
    prefixes = {}
    for name in DOCS:
        p = root/f'docs/{name}.md'; prefixes[p.relative_to(root).as_posix()] = dict(bytes=p.stat().st_size, sha256=sha(p))
        with p.open('ab') as h: h.write(doc.encode())
    write(out/'result.json', final); write(out/'publication.json', dict(previous_nondoc=preserved, document_prefixes=prefixes))
    bound = NEW + CHANGED + [root/p for p in prefixes] + [p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes, sha256={p.relative_to(root).as_posix(): sha(p) for p in bound})
    write(manifest, final)
    raw = master.read_text(encoding='utf-8'); m = json.loads(raw); indent = 2 if raw.startswith('{\n  "') else 1
    m['sha256'].update(final['sha256']); m['sha256'][manifest.relative_to(root).as_posix()] = sha(manifest)
    m['native_focused_final_review'] = {k: v for k, v in final.items() if k not in ('sha256', 'document_prefixes')}
    master.write_text(json.dumps(m, indent=indent, ensure_ascii=False) + '\n', encoding='utf-8', newline='\n')


def check(mode):
    m = read(manifest); a = read(master)
    for p, h in m['sha256'].items(): assert sha(root/p) == h and a['sha256'][p] == h, p
    for p, h in read(out/'publication.json')['previous_nondoc'].items():
        if root/p not in CHANGED: assert sha(root/p) == h, p
    for p, v in m['document_prefixes'].items(): assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest() == v['sha256'], p
    assert a['sha256'][manifest.relative_to(root).as_posix()] == sha(manifest)
    provenance((root/'paper/manuscript.md').read_bytes().decode('utf-8'))
    paths = list(m['sha256']) + [manifest.relative_to(root).as_posix(), 'paper/revision-manifest.json']
    (S/'phase292-paths.txt').write_bytes(b'\0'.join(p.encode() for p in paths) + b'\0')
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
