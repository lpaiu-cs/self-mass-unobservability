"""Archive phase 286: confirmation review of revision 5, revision 6 of the white-dwarf section, and its reintegration as Section 4.6
of the unified manuscript (manuscript, main.tex regenerated with Pandoc 3.11, references, README).

Binds the paper files, the draft, the response note and the confirmation-review records; appends one section to the six
dynamic-chi documents; binds everything in the manifest and paper/revision-manifest.json. No new computation.
Usage: python phase286-publish.py package|check|staged|head
"""
from pathlib import Path
import datetime, hashlib, json, subprocess, sys
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
S = Path(__file__).resolve().parent
OLD = S.parent.parent/'c8fbf92c-e1e4-431b-9665-aff6bfed2aa3'/'scratchpad'
master = root/'paper/revision-manifest.json'; previous_manifest = root/'outputs/direct-eos-gr33/native-manuscript-rereview2-manifest.json'
out = root/'outputs/direct-eos-gr33/native-manuscript-reintegration'; manifest = out.parent/'native-manuscript-reintegration-manifest.json'
PAPER = [root/f'paper/{n}' for n in ('manuscript.md', 'main.tex', 'references.bib', 'README.md')]
draft = root/'docs/white-dwarf-free-fall-charge-section.md'
notes = [root/'notes/REQUEST286_MANUSCRIPT_REINTEGRATION_KO.md']
DOCS = ['model-definition', 'observable-targets', 'adiabatic-limit', 'nonadiabatic-regime', 'failure-ledger-dynamic-chi', 'dynamic-charge-completion']
SCRIPTS = {'phase286-integrate.py': S/'phase286-integrate.py', 'phase286-publish.py': S/'phase286-publish.py'}
REVIEWS = {'confirm-prompt.md': 'confirm285-prompt.md', 'confirm-gpt-6-astra-prompt.md': 'confirm285-astra-prompt.md',
           'confirm-gpt-6-astra.md': 'confirm285-astra.md', 'confirm-gpt-6-astra-cli.log': 'confirm285-astra.log',
           'confirm-opus-5.5.md': 'confirm285-opus.md', 'confirm-fable-5.1.md': 'confirm285-fable.md'}
PRE_EDIT_MAIN_TEX = '29bb4eaf94c500132a6a178b576866fa79c34871b5b3726bd9b9b2022c100d95'
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
    for p in SCRIPTS.values(): assert p.is_file(), p
    for name in REVIEWS.values(): assert (OLD/name).is_file(), name
    preserved = {p: h for p, h in read(previous_manifest)['sha256'].items() if not p.startswith('docs/')}
    for p, h in preserved.items(): assert sha(root/p) == h, p
    head_main = hashlib.sha256(git('show', 'HEAD:paper/main.tex')).hexdigest(); assert head_main == PRE_EDIT_MAIN_TEX
    for p in PAPER + [draft]: assert hashlib.sha256(git('show', 'HEAD:' + p.relative_to(root).as_posix())).hexdigest() != sha(p), p
    out.mkdir(parents=True)
    for name, p in SCRIPTS.items(): copy(p, out/'scripts'/name)
    for dst, name in REVIEWS.items(): copy(OLD/name, out/dst)
    now = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    final = dict(classification='Imported from prior work', passed=True, verdict='CONFIRMATION_REVIEW_PASSED__REVISION_6_INTEGRATED_AS_SECTION_4_6',
        confirmation_review_of_revision_5={'fable5.1 (subagent)': 'accept', 'gpt-6-astra (Codex CLI, read-only)': 'accept after minor revision',
                                           'opus5.5 (subagent)': 'accept after minor revision'},
        revision_6=['instantaneous response: full bound 1.24e-9 a_p^2 + 8.2e-10 |a_p a_o|; reaches only the smallest stored scales for |a_p| >~ 0.5; static SEP limits (<= 2.6e-6) give |a_p| <~ 4e-3 for Cassini-scale a_o, or need |a_o| <~ 5e-6',
                    'outer white dwarf has the analogous zero-lag response (leading-order estimate 1.86e-9 a_p^2)',
                    'structural parts each <= 2.13e-18 (4.3e-18 together); about eight orders',
                    'summary: deeper layers with orbital-scale local thermal times expected; pole structure and coupling not computed',
                    'abstract: permanent charge zero only as t -> infinity with slow thermal modes damped',
                    'data availability: verdict, classification and derived-bound fields superseded; rerun limit narrowed; Born source'],
        integration=dict(manuscript_lines=[697, 850], pandoc='3.11', main_tex_reproduced_before_edit_sha256=PRE_EDIT_MAIN_TEX,
                         edits=['date 28 September 2026', 'abstract: two Conjectural sentences', 'Section 4.6 after Section 4.5',
                                'Section 6: one paragraph after the physical-target paragraph', 'data availability: paragraph, list, paragraph',
                                'references: three entries restored from a322c6462', 'README: PDF and submission archive predate Section 4.6'],
                         pdf_and_submission_archive='predate Section 4.6 (no TeX on this host)'),
        errata=['16 manifest numeric and classification fields', '17 REQUEST276 line 47 last sentence', '18 failure ledger lines 2578 and 2588',
                '19 REQUEST285 details (tracked arrays, run commands, 1.8 ms, line numbers)'],
        repository_classification='conditional theorem progress on the lag amplitude within the chosen response family; no observational exclusion; A4 not established',
        full_goal_complete=False, snapshot_KST=now)
    doc = ("\n\n## 단계286 — 반영 확인 재심 통과와 원고 재통합\n\n"
        "분류: Imported from prior work. 5판의 반영 확인 재심에서 fable5.1은 수락, gpt-6-astra와 opus5.5는 경미 수정 후 수락이었다. "
        "경미 지적을 반영한 6판을 통합 원고 §4.6으로 다시 넣었다. 6판은 순간 응답의 전체 상한과 정적 SEP 양립 조건(Cassini 수준 a_o에서 |a_p|≲4e−3), 외측 백색왜성의 같은 순간 응답, "
        "구조 변조의 합계(각 2.13e−18, 합 4.3e−18), 초록의 t→∞ 한정을 담는다. 초록·§6·데이터 가용성 문장을 함께 넣었고, 참고문헌 세 항목을 복원했다. "
        "main.tex는 Pandoc 3.11로 다시 만들었다(고치기 전 원고에서 바이트 재현을 먼저 확인). PDF와 제출 zip은 §4.6 이전 판이다(TeX 없음). "
        "저장소 분류는 선택한 응답족의 지연 진폭 상한에 관한 조건부 theorem progress다. 관측 배제도 A4 유지도 아니다. "
        "[근거 286](../notes/REQUEST286_MANUSCRIPT_REINTEGRATION_KO.md).\n")
    prefixes = {}
    for name in DOCS:
        p = root/f'docs/{name}.md'; prefixes[p.relative_to(root).as_posix()] = dict(bytes=p.stat().st_size, sha256=sha(p))
        with p.open('ab') as h: h.write(doc.encode())
    write(out/'result.json', final); write(out/'publication.json', dict(previous_nondoc=preserved, document_prefixes=prefixes))
    bound = notes + PAPER + [draft] + [root/p for p in prefixes] + [p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes, sha256={p.relative_to(root).as_posix(): sha(p) for p in bound})
    write(manifest, final)
    raw = master.read_text(encoding='utf-8'); m = json.loads(raw); indent = 2 if raw.startswith('{\n  "') else 1
    m['sha256'].update(final['sha256']); m['sha256'][manifest.relative_to(root).as_posix()] = sha(manifest)
    m['native_manuscript_reintegration'] = {k: v for k, v in final.items() if k not in ('sha256', 'document_prefixes')}
    master.write_text(json.dumps(m, indent=indent, ensure_ascii=False) + '\n', encoding='utf-8', newline='\n')


def check(mode):
    m = read(manifest); a = read(master)
    for p, h in m['sha256'].items(): assert sha(root/p) == h and a['sha256'][p] == h, p
    for p, h in read(out/'publication.json')['previous_nondoc'].items(): assert sha(root/p) == h, p
    for p, v in m['document_prefixes'].items(): assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest() == v['sha256'], p
    assert a['sha256'][manifest.relative_to(root).as_posix()] == sha(manifest)
    paths = list(m['sha256']) + [manifest.relative_to(root).as_posix(), 'paper/revision-manifest.json']
    (S/'phase286-paths.txt').write_bytes(b'\0'.join(p.encode() for p in paths) + b'\0')
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
