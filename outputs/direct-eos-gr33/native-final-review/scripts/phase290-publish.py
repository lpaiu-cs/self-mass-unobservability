"""Archive phase 290: final pre-submission review (fixes, AI model versions, bibliography note removed and DOIs added, public-snapshot
data-availability paragraph, PUBLIC_SNAPSHOT.md, repackaged PDF and source archive, acceptance assessment).
Usage: python phase290-publish.py package|check|staged|head
"""
from pathlib import Path
import datetime, hashlib, json, subprocess, sys
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
S = Path(__file__).resolve().parent
master = root/'paper/revision-manifest.json'; previous_manifest = root/'outputs/direct-eos-gr33/native-submission-prep-manifest.json'
out = root/'outputs/direct-eos-gr33/native-final-review'; manifest = out.parent/'native-final-review-manifest.json'
CHANGED = [root/p for p in ('paper/manuscript.md', 'paper/main.tex', 'paper/references.bib', 'docs/white-dwarf-free-fall-charge-section.md',
                            'output/pdf/free-fall-identifiability.pdf', 'output/submission/free-fall-identifiability-source.zip',
                            'output/submission/submission-checklist-prd.md')]
NEW = [root/p for p in ('notes/REQUEST290_FINAL_REVIEW_KO.md', 'PUBLIC_SNAPSHOT.md')]
DOCS = ['model-definition', 'observable-targets', 'adiabatic-limit', 'nonadiabatic-regime', 'failure-ledger-dynamic-chi', 'dynamic-charge-completion']
SCRIPTS = ['phase290-ai.py', 'phase290-fixes.py', 'final290-checks.py', 'phase290-snapshot.py', 'phase290-publish.py']
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
    assert b'GPT-6-Astra (OpenAI' in ms and b'public snapshot' in ms and b'centre' not in ms
    assert b'Source of the' not in (root/'paper/references.bib').read_bytes() and b'Source of the' not in (root/'paper/build/main.bbl').read_bytes()
    assert len([l for l in ms.split(b'\n')[:-1] if not l.endswith(b'\r')]) == 1
    for name in SCRIPTS: copy(S/name, out/'scripts'/name)
    copy(root/'paper/build/main.log', out/'build/main.log')
    now = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    final = dict(classification='Imported from prior work (review); Conjectural (acceptance assessment)', passed=True,
        verdict='FINAL_PRE_SUBMISSION_REVIEW_DONE__NO_BLOCKING_ERRORS__AFFILIATION_ORCID_PENDING',
        author_confirmations=['AI-use section factual; model versions Claude Opus 5.5 and GPT-6-Astra', 'no prior submission, no competing interests',
                              'affiliation, e-mail, ORCID to be filled later', 'public snapshot instead of full-history push (51.6 GB exceeds GitHub limits)'],
        checks=['full reading of Sections 1-6, AI use, data availability, appendices A-E', 'bibliography: 16 cited entries all resolve (arXiv, Crossref, DataCite)',
                'no undefined references; overfull 0, underfull 1', 'rendered pages: title, 9-15, 25-26, 30'],
        fixes=['printed internal bibliography note removed', 'DOIs added: voisin2025planet, damour1992tensor', 'AI model versions',
               'American spelling (center x3, meters)', 'inline math c_chi and pi x4', 'Nançay', 'en dashes x3', 'internal ledger reference removed',
               'Figure 1 reference', 'data availability for the public snapshot'],
        remaining_recommendations=['affiliation and ORCID (required)', 'dense provenance-heavy presentation, 30 pages', 'nonstandard claim labels',
                                   'state novelty against interpolation/moment-problem literature', 'non-detection and conditional results', 'cited commits not publicly verifiable'],
        acceptance_assessment=dict(base_rate='PRD acceptance not published by APS; third-party estimates 50-65%', desk_rejection='about 25-40%',
                                   accept_if_reviewed='about 35-55% after major revision', overall='about 15-35% (central about 25%)',
                                   levers=['focus and shorten to 15-20 pages', 'explicit novelty statement', 'separate the white-dwarf calculation', 'human expert pre-read', 'REVTeX']),
        full_goal_complete=False, snapshot_KST=now)
    doc = ("\n\n## 단계290 — 제출 전 최종 검토\n\n"
        "분류: Imported from prior work. 원고 전문을 정독하고 자동 점검(참고문헌 실재·교차 참조·조판·표기)과 쪽 렌더링을 했다. 제출을 막는 오류는 없었다. "
        "고친 것은 다음과 같다: 참고문헌에 인쇄되던 내부 메모 삭제, DOI 2건 추가, AI 모델 버전(Claude Opus 5.5, GPT-6-Astra), 미국식 철자·수식·en dash 정리, 저장소 내부 표현 삭제, 데이터 가용성의 공개 스냅숏 문구. "
        "전체 이력(51.6 GB)은 GitHub 한도를 넘어, 저자 결정에 따라 대형 배열을 뺀 공개 스냅숏으로 올린다(`PUBLIC_SNAPSHOT.md`). 남은 필수 항목은 소속·ORCID뿐이다.\n\n"
        "분류: Conjectural. PRD 게재 가능성은 약 15–35%(중심 약 25%)로 본다. 초기 반려 약 25–40%, 심사로 가면 약 35–55%다. "
        "주된 약점은 분량·문체, 핵심 정리의 제한된 새로움, 비검출·조건부 결과다. 초점을 좁히고 새로움의 위치를 명시하면 가능성이 오른다. "
        "[근거 290](../notes/REQUEST290_FINAL_REVIEW_KO.md).\n")
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
    m['native_final_review'] = {k: v for k, v in final.items() if k not in ('sha256', 'document_prefixes')}
    master.write_text(json.dumps(m, indent=indent, ensure_ascii=False) + '\n', encoding='utf-8', newline='\n')


def check(mode):
    m = read(manifest); a = read(master)
    for p, h in m['sha256'].items(): assert sha(root/p) == h and a['sha256'][p] == h, p
    for p, h in read(out/'publication.json')['previous_nondoc'].items():
        if root/p not in CHANGED: assert sha(root/p) == h, p
    for p, v in m['document_prefixes'].items(): assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest() == v['sha256'], p
    assert a['sha256'][manifest.relative_to(root).as_posix()] == sha(manifest)
    paths = list(m['sha256']) + [manifest.relative_to(root).as_posix(), 'paper/revision-manifest.json']
    (S/'phase290-paths.txt').write_bytes(b'\0'.join(p.encode() for p in paths) + b'\0')
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
