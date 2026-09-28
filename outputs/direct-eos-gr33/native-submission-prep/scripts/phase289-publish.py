"""Archive phase 289: Physical Review D submission preparation (PDF built with Tectonic 0.17.0, math fix, AI-use section, status
line removed, journal source archive, cover letter, plain abstract, checklist, README update).
Usage: python phase289-publish.py package|check|staged|head
"""
from pathlib import Path
import datetime, hashlib, json, subprocess, sys
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
S = Path(__file__).resolve().parent
master = root/'paper/revision-manifest.json'; previous_manifest = root/'outputs/direct-eos-gr33/native-section46-relaxation-manifest.json'
out = root/'outputs/direct-eos-gr33/native-submission-prep'; manifest = out.parent/'native-submission-prep-manifest.json'
CHANGED = [root/p for p in ('paper/manuscript.md', 'paper/main.tex', 'paper/README.md', 'docs/white-dwarf-free-fall-charge-section.md',
                            'output/pdf/free-fall-identifiability.pdf', 'output/submission/free-fall-identifiability-source.zip')]
NEW = [root/p for p in ('notes/REQUEST289_PRD_SUBMISSION_PREP_KO.md', 'output/submission/cover-letter-prd.tex', 'output/submission/cover-letter-prd.pdf',
                        'output/submission/abstract-plain.txt', 'output/submission/submission-checklist-prd.md')]
DOCS = ['model-definition', 'observable-targets', 'adiabatic-limit', 'nonadiabatic-regime', 'failure-ledger-dynamic-chi', 'dynamic-charge-completion']
SCRIPTS = ['phase289-fixmath.py', 'phase289-manuscript.py', 'phase289-package.py', 'find289-math.py', 'render-pages.py', 'phase289-publish.py']
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
    for name in SCRIPTS: assert (S/name).is_file(), name
    for p in NEW: assert p.is_file(), p
    preserved = {p: h for p, h in read(previous_manifest)['sha256'].items() if not p.startswith('docs/')}
    for p, h in preserved.items():
        if root/p not in CHANGED: assert sha(root/p) == h, p
    for p in CHANGED: assert hashlib.sha256(git('show', 'HEAD:' + p.relative_to(root).as_posix())).hexdigest() != sha(p), p
    ms = (root/'paper/manuscript.md').read_bytes()
    assert b'Use of AI tools' in ms and b'Unified revised manuscript' not in ms and b'(K\\(\\approx\\)934)' not in ms
    assert len([l for l in ms.split(b'\n')[:-1] if not l.endswith(b'\r')]) == 1
    for name in SCRIPTS: copy(S/name, out/'scripts'/name)
    copy(root/'paper/build/main.log', out/'build/main.log'); copy(S/'cover-build.log', out/'build/cover-letter.log')
    now = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    final = dict(classification='Imported from prior work', passed=True, verdict='PRD_SUBMISSION_PREPARED__AUTHOR_ACTIONS_AND_REPOSITORY_PUSH_PENDING',
        venue=dict(primary='Physical Review D, Regular Article', alternative='Classical and Quantum Gravity, Research Paper',
                   reasons=['scope: gravitation, compact objects, scalar-tensor theories', 'APS AI policy (June 2026) allows substantive AI use with in-paper disclosure',
                            'PDF suffices for initial review; REVTeX recommended, not required; no length limit; no mandatory charges']),
        tex='Tectonic 0.17.0 (official GitHub release, windows-msvc zip, SHA-256 digest matched), user-level install',
        fixes=['(K\\(\\approx\\)934) -> (\\(K\\approx934\\)): Pandoc does not read math whose closing $ is followed by a digit',
               'Status metadata line removed (title block read "Unified revised manuscript")', 'Use of AI tools section added before Data and code availability'],
        pdf=dict(pages=30, warnings='one underfull hbox; no undefined citations or references', visual_check='title page, Section 4.6 pages 9-15, AI section page 25 (pypdfium2 render)'),
        source_archive='main.tex, main.bbl, references.bib, four figures, README.txt; compiles on its own to a PDF of identical size',
        author_actions=['confirm the AI-use section and add earlier model versions', 'affiliation, e-mail, ORCID', 'push the repository (local main 1137 commits ahead of origin, last push 2026-07-12) or deposit a DOI snapshot',
                        'confirm no prior submission and no competing interests', 'optional suggested referees', 'submit at authors.aps.org'],
        package_revision_hazard='paper/package_revision.py rewrites revision-manifest.json from its 9 September list; not used, README warns',
        full_goal_complete=False, snapshot_KST=now)
    doc = ("\n\n## 단계289 — Physical Review D 제출 준비\n\n"
        "분류: Imported from prior work. 사용자 지시(2026-09-28)로 PDF를 만들고 투고처를 조사해 PRD(Regular Article)를 1순위로 정했다. "
        "근거는 범위와 APS의 2026년 6월 AI 정책이다. 이 정책은 실질적 AI 사용을 논문 안에 공개하는 조건으로 허용한다. 대안은 CQG다. "
        "Tectonic 0.17.0을 설치해 30쪽 PDF를 빌드했다. 빌드 중 Pandoc이 수식을 놓친 §4.6의 `K≈934` 한 곳을 고쳤다. "
        "원고에는 'Use of AI tools' 절을 넣고 Status 줄을 지웠다. 저널용 소스 zip(자체 컴파일 확인), 1쪽 커버레터, 평문 초록, 체크리스트를 `output/submission/`에 두었다. "
        "남은 일은 저자 몫이다: AI 절 확인, 소속·ORCID, 공개 저장소 push(로컬이 1,137커밋 앞섬), 제출. "
        "`paper/package_revision.py`는 manifest를 덮어쓰므로 쓰지 않았다. [근거 289](../notes/REQUEST289_PRD_SUBMISSION_PREP_KO.md).\n")
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
    m['native_submission_prep'] = {k: v for k, v in final.items() if k not in ('sha256', 'document_prefixes')}
    master.write_text(json.dumps(m, indent=indent, ensure_ascii=False) + '\n', encoding='utf-8', newline='\n')


def check(mode):
    m = read(manifest); a = read(master)
    for p, h in m['sha256'].items(): assert sha(root/p) == h and a['sha256'][p] == h, p
    for p, h in read(out/'publication.json')['previous_nondoc'].items():
        if root/p not in CHANGED: assert sha(root/p) == h, p
    for p, v in m['document_prefixes'].items(): assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest() == v['sha256'], p
    assert a['sha256'][manifest.relative_to(root).as_posix()] == sha(manifest)
    paths = list(m['sha256']) + [manifest.relative_to(root).as_posix(), 'paper/revision-manifest.json']
    (S/'phase289-paths.txt').write_bytes(b'\0'.join(p.encode() for p in paths) + b'\0')
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
