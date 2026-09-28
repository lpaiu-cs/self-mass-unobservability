"""Archive phase 293: confirmation reviews of the phase-292 revision (Claude Opus 5.5, Claude Fable 5.1, GPT-6-Astra; all minor revision
then accept) and the resulting minor fixes to the paper, SM and submission files.
Usage: python phase293-publish.py package|check|staged|head
"""
from pathlib import Path
import datetime, hashlib, importlib.util, json, math, subprocess, sys
import pypdfium2 as pdfium
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
S = Path(__file__).resolve().parent
master = root/'paper/revision-manifest.json'; previous_manifest = root/'outputs/direct-eos-gr33/native-focused-final-review-manifest.json'
out = root/'outputs/direct-eos-gr33/native-confirmation-review'; manifest = out.parent/'native-confirmation-review-manifest.json'
RC = root/'outputs/research-completion'
CHANGED = [root/p for p in ('paper/manuscript.md', 'paper/main.tex', 'paper/supplement.md', 'paper/supplement.tex', 'paper/references.bib',
    'output/pdf/free-fall-identifiability.pdf', 'output/pdf/free-fall-identifiability-supplement.pdf', 'output/submission/free-fall-identifiability-source.zip',
    'output/submission/cover-letter-prd.tex', 'output/submission/cover-letter-prd.pdf', 'output/submission/abstract-plain.txt',
    'output/submission/submission-checklist-prd.md')]
NEW = [root/'notes/REQUEST293_CONFIRMATION_REVIEW_KO.md']
DOCS = ['model-definition', 'observable-targets', 'adiabatic-limit', 'nonadiabatic-regime', 'failure-ledger-dynamic-chi', 'dynamic-charge-completion',
        'remaining-levers-2026-09-09']
SCRIPTS = ['phase293-edits.py', 'phase293-bib.py', 'phase293-publish.py']
REVIEWS = {'confirm292-prompt.md': 'confirmation-prompt.md', 'confirm292-astra-prompt.md': 'confirmation-gpt-6-astra-prompt.md',
           'confirm292-opus-prompt.md': 'confirmation-opus-5.5-prompt.md', 'confirm292-fable-prompt.md': 'confirmation-fable-5.1-prompt.md',
           'confirm292-astra.md': 'confirmation-gpt-6-astra.md', 'confirm292-astra.log': 'confirmation-gpt-6-astra-cli.log',
           'confirm292-opus.md': 'confirmation-opus-5.5.md', 'confirm292-fable.md': 'confirmation-fable-5.1.md'}
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
def tex_exp(x):
    mantissa, exponent = f'{x:.1e}'.split('e')
    return f'{mantissa}\\times10^{{{int(exponent)}}}'


def provenance(ms):
    """Phase-292 number checks (archived publisher) plus the numbers this phase adds."""
    spec = importlib.util.spec_from_file_location('p292', root/'outputs/direct-eos-gr33/native-focused-final-review/scripts/phase292-publish.py')
    p292 = importlib.util.module_from_spec(spec); spec.loader.exec_module(p292)
    checked = p292.provenance(ms)
    cal = read(RC/'simultaneous-calibration.json'); v = read(RC/'simultaneous-validation.json'); d = read(RC/'corrected-physical-drive.json')
    stats = [r['order_statistic_95'] for r in cal['rows']]; q = cal['nominal_chi6']
    assert max(stats) == cal['threshold']
    density = q*q*math.exp(-q/2)/16; se = math.sqrt(.95*.05/8192)/density
    sections = v['data']['physical_lag_sections']
    betas = {s['tau']: s['beta'] for s in sections}
    assert all(b < 0 for t, b in betas.items() if t != 18.) and betas[18.] > 0
    assert max(s['minimum_statistic'] for s in sections) == next(s['minimum_statistic'] for s in sections if s['tau'] == 18.)
    sig = [f['sigma'] for f in d['fits']]; u = d['drive']['Ustar']
    want = [f'{min(stats):.2f}', f'{max(stats):.2f}', f'{se:.2f}', f'{q:.2f}', tex_exp(betas[18.]), tex_exp(min(sig)), tex_exp(max(sig)),
            f'about {round(min(sig)/u)} to {round(max(sig)/u)}']
    missing = [w for w in want if w not in ms]; assert not missing, missing
    for gone in ('Every relaxing response to this drive lies in the rejected plane', 'empirical promotion', 'factors of 2.1--27.5'):
        assert gone not in ms, gone
    return checked + want


def package():
    assert not out.exists() and not manifest.exists()
    for name in SCRIPTS + list(REVIEWS): assert (S/name).is_file(), name
    for p in NEW: assert p.is_file(), p
    preserved = {p: h for p, h in read(previous_manifest)['sha256'].items() if not p.startswith('docs/')}
    for p, h in preserved.items():
        if root/p not in CHANGED: assert sha(root/p) == h, p
    for p in CHANGED: assert hashlib.sha256(git('show', 'HEAD:' + p.relative_to(root).as_posix())).hexdigest() != sha(p), p
    ms = (root/'paper/manuscript.md').read_bytes().decode('utf-8'); assert '\n' not in ms.replace('\r\n', '')
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
    final = dict(classification='Imported from prior work (confirmation review and wording fixes); Proven (sign-constrained section bound)', passed=True,
        verdict='CONFIRMATION_REVIEW_DONE__ALL_MINOR_REVISION_THEN_ACCEPT__FIXES_APPLIED__FINAL_CHARGE_CONCLUSION_UNCHANGED',
        reviews={'GPT-6-Astra': 'minor revision then accept; all earlier points resolved; B1 confirmed', 'Claude Fable 5.1': 'minor revision then accept; B1 confirmed',
                 'Claude Opus 5.5': 'minor revision then accept; B1 resolved; scope and marginal 2-day rejection flagged'},
        fixes=['conclusion limited to six evaluated lags and an instantaneous plus single-relaxation response', 'marginal 2-day rejection quantified (order statistics 12.27-12.82, SE 0.13)',
               'sign-constrained bound: with beta >= 0 every section minimum >= 14.59 (convex quadratic, beta = 0 line in every section)',
               'archived-derivative condition', 'relation of p = 0.26 and p = 0.012', 'ratio range wording', 'diagonal weighting of both K=1 columns',
               'unit amplitudes in comparator phase sets', 'REML and timing-noise citations', 'sensitivity scale sigma_B about 2-10', 'internal wording',
               'SM alignment (Repository line, Theorem 3, Section 3.5, Section 5.3, Section 5.9, AI section, process wording)'],
        declined=['eta/kappa as symbols (kappa already denotes the charge stiffness)', 'section minima on all 65 lags (new computation)', 'Cholesky draws (would change the registered procedure)',
                  'six-coefficient estimate table (not in the stored record)'],
        new_number_checks=checked, pages=dict(paper=13, supplement=30, cover_letter=1),
        remaining=['affiliation, e-mail, ORCID', 'public snapshot update (author instruction)', 'cause of the six-coefficient excess'],
        full_goal_complete=False, snapshot_KST=now)
    doc = ("\n\n## 단계293 — 확인 심사와 경미 수정\n\n"
        "분류: Imported from prior work. 단계292 개정판을 세 심사자(Opus 5.5, Fable 5.1, GPT-6-Astra)가 확인 심사했다. 모두 경미 수정 후 수락을 권했고, "
        "세 명 모두 근점 규약 정정을 독립적으로 확인했다. J0337 결론의 범위를 평가한 여섯 지연과 순간항+단일 완화 모형으로 한정했다. "
        "2일 단면의 근소한 기각을 정량화했고(보정 순서통계량 12.27–12.82, 표준오차 약 0.13), 보관된 도함수 조건을 명시했다.\n\n"
        "분류: Proven. 통계량은 볼록 이차식이고 β=0 직선이 모든 단면에 들어 있다. 제약 없는 β̂은 18일을 빼면 음수다. "
        "따라서 등전하 실현의 β≥0 아래에서는 여섯 단면 모두 최솟값이 14.59 이상이다. 백색왜성의 최종 전하 결론은 유지된다. "
        "[근거 293](../notes/REQUEST293_CONFIRMATION_REVIEW_KO.md).\n")
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
    m['native_confirmation_review'] = {k: v for k, v in final.items() if k not in ('sha256', 'document_prefixes')}
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
    (S/'phase293-paths.txt').write_bytes(b'\0'.join(p.encode() for p in paths) + b'\0')
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
