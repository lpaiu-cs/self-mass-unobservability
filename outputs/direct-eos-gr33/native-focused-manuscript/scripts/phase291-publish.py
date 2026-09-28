"""Archive phase 291: focused paper (paper/manuscript.md, 11 pages) and Supplemental Material (paper/supplement.md, the previous full
manuscript with a new title and abstract, 30 pages); builder, README, submission files; no new computation.
Usage: python phase291-publish.py package|check|staged|head
"""
from pathlib import Path
import datetime, hashlib, json, re, subprocess, sys
import pypdfium2 as pdfium
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
S = Path(__file__).resolve().parent
PREVIOUS_COMMIT = '39962a753'
master = root/'paper/revision-manifest.json'; previous_manifest = root/'outputs/direct-eos-gr33/native-final-review-manifest.json'
out = root/'outputs/direct-eos-gr33/native-focused-manuscript'; manifest = out.parent/'native-focused-manuscript-manifest.json'
CHANGED = [root/p for p in ('paper/manuscript.md', 'paper/main.tex', 'paper/build_manuscript.py', 'paper/README.md',
                            'output/pdf/free-fall-identifiability.pdf', 'output/submission/free-fall-identifiability-source.zip',
                            'output/submission/cover-letter-prd.tex', 'output/submission/cover-letter-prd.pdf',
                            'output/submission/abstract-plain.txt', 'output/submission/submission-checklist-prd.md')]
NEW = [root/p for p in ('notes/REQUEST291_FOCUSED_MANUSCRIPT_KO.md', 'paper/supplement.md', 'paper/supplement.tex',
                        'output/pdf/free-fall-identifiability-supplement.pdf')]
DOCS = ['model-definition', 'observable-targets', 'adiabatic-limit', 'nonadiabatic-regime', 'failure-ledger-dynamic-chi', 'dynamic-charge-completion']
SCRIPTS = ['phase291-split.py', 'check291-numbers.py', 'check291-tex.py', 'phase291-package.py', 'phase291-publish.py']
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


def package():
    assert not out.exists() and not manifest.exists()
    for name in SCRIPTS: assert (S/name).is_file(), name
    for p in NEW: assert p.is_file(), p
    preserved = {p: h for p, h in read(previous_manifest)['sha256'].items() if not p.startswith('docs/')}
    for p, h in preserved.items():
        if root/p not in CHANGED: assert sha(root/p) == h, p
    for p in CHANGED: assert hashlib.sha256(git('show', 'HEAD:' + p.relative_to(root).as_posix())).hexdigest() != sha(p), p
    old = git('show', f'{PREVIOUS_COMMIT}:paper/manuscript.md').decode('utf-8')
    ms = (root/'paper/manuscript.md').read_bytes().decode('utf-8'); sm = (root/'paper/supplement.md').read_bytes().decode('utf-8')
    # The SM is the previous full text below its abstract; the paper is CRLF throughout and adds no number absent from the previous text.
    assert sm.startswith('# Supplemental Material for ``Identifying a relaxing') and sm.split('---\r\n', 1)[1] == old.split('---\r\n', 1)[1]
    assert ms.startswith('# Identifying a relaxing internal state') and '\n' not in ms.replace('\r\n', '')
    tok = lambda s: re.findall(r'(?<![\w.])\d+(?:[.,]\d+)*(?![\w])', s)
    new_numbers = sorted(set(tok(ms)) - set(tok(old))); assert not new_numbers, new_numbers
    n_main, n_sm = pages(root/'output/pdf/free-fall-identifiability.pdf'), pages(root/'output/pdf/free-fall-identifiability-supplement.pdf')
    assert (n_main, n_sm, pages(root/'output/submission/cover-letter-prd.pdf')) == (11, 30, 1)
    for name in ('main', 'supplement'):
        log = (root/f'paper/build/{name}.log').read_text(encoding='utf-8', errors='replace')
        assert 'Overfull' not in log and 'undefined' not in log.lower(), name
    for name in SCRIPTS: copy(S/name, out/'scripts'/name)
    for name in ('main', 'supplement'): copy(root/f'paper/build/{name}.log', out/f'build/{name}.log')
    now = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    final = dict(classification='Imported from prior work (presentation only; no new computation or numerical claim)', passed=True,
        verdict='FOCUSED_PAPER_AND_SUPPLEMENT_BUILT__NO_NEW_CLAIMS__FINAL_CHARGE_CONCLUSION_UNCHANGED',
        user_instruction='narrow the focus; do not repeat known facts (2026-09-28)',
        title='Identifying a relaxing internal state in free-fall timing: finite-frequency boundaries and an application to PSR J0337+1715',
        pages=dict(paper=n_main, supplement=n_sm, previous_paper=30, cover_letter=1),
        paper_keeps=['Theorem 1 (previous Theorem 3)', 'rate-gap comparator inequality', 'two-sample moment-variance identification and close-pole limit',
                     'transient insufficiency of a finite periodic record', 'damped scalar-charge realization and phase closure',
                     'J0337 application: Tables 1-2, Figure 1, nuisance-span dependence, coverage, comparator information, corrected drive, record limits'],
        moved_to_supplement=['static operator catalog (SM Section 2, Appendices A-C)', 'nuisance-projected rank and pair-coupling benchmark (SM 4.1-4.2)',
                             'white-dwarf calculation (SM 4.6)', 'full audits and provenance (SM 5.2, 5.5-5.10, Appendices D-E)', 'Figures 2-4'],
        checks=[f'all numbers of the paper occur in {PREVIOUS_COMMIT}:paper/manuscript.md', 'SM body below the abstract is byte-identical to the previous manuscript',
                'no undefined references, no overfull boxes, no dropped inline math', 'builder regenerates the previous main.tex byte for byte',
                'source archive compiles standalone (11 pages)', 'verify_unified_paper table checks pass on the new paper'],
        final_charge_conclusion='unchanged: no new computation; the Discussion summarizes SM Section 4.6 with the same numbers',
        remaining=['affiliation, e-mail, ORCID (author)', 'public snapshot update: the paper points to paper/supplement.md (needs author instruction to push)', 'submission'],
        full_goal_complete=False, snapshot_KST=now)
    doc = ("\n\n## 단계291 — 원고 초점 축소\n\n"
        "분류: Imported from prior work. 사용자 지시로 PRD 제출 원고의 초점을 좁혔다. "
        "본문(11쪽, 새 제목 *Identifying a relaxing internal state in free-fall timing: finite-frequency boundaries and an application to PSR J0337+1715*)은 "
        "식별성 경계(유한 반송파 보간 경계, 속도 간극 부등식, 두 표본 모멘트-분산 식별과 가까운 극의 한계, 유한 주기 기록의 과도 응답 불충분성), "
        "감쇠 스칼라 전하 실현과 위상 폐합, J0337 조건부 적용만 남긴다. 알려진 사실은 한 줄이나 인용으로 줄였다. "
        "이전 전체 원고는 제목·초록만 바꿔 보충 자료(`paper/supplement.md`, 30쪽)로 옮겼다. "
        "새 계산·새 수치는 없고(본문 수치가 모두 이전 원고에 있다), 최종 전하 결론은 유지된다. "
        "[근거 291](../notes/REQUEST291_FOCUSED_MANUSCRIPT_KO.md).\n")
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
    m['native_focused_manuscript'] = {k: v for k, v in final.items() if k not in ('sha256', 'document_prefixes')}
    master.write_text(json.dumps(m, indent=indent, ensure_ascii=False) + '\n', encoding='utf-8', newline='\n')


def check(mode):
    m = read(manifest); a = read(master)
    for p, h in m['sha256'].items(): assert sha(root/p) == h and a['sha256'][p] == h, p
    for p, h in read(out/'publication.json')['previous_nondoc'].items():
        if root/p not in CHANGED: assert sha(root/p) == h, p
    for p, v in m['document_prefixes'].items(): assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest() == v['sha256'], p
    assert a['sha256'][manifest.relative_to(root).as_posix()] == sha(manifest)
    paths = list(m['sha256']) + [manifest.relative_to(root).as_posix(), 'paper/revision-manifest.json']
    (S/'phase291-paths.txt').write_bytes(b'\0'.join(p.encode() for p in paths) + b'\0')
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
