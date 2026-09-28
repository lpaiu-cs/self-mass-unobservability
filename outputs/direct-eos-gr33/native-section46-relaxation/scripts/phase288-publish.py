"""Archive phase 288: the layer-by-layer thermal-relaxation estimate of phase 287 added to Section 4.6 after a single self-review
(user decision); draft revision 7, manuscript and main.tex regenerated with Pandoc 3.11.
Usage: python phase288-publish.py package|check|staged|head
"""
from pathlib import Path
import datetime, hashlib, json, subprocess, sys
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
S = Path(__file__).resolve().parent
master = root/'paper/revision-manifest.json'; previous_manifest = root/'outputs/direct-eos-gr33/native-thermal-relaxation-manifest.json'
out = root/'outputs/direct-eos-gr33/native-section46-relaxation'; manifest = out.parent/'native-section46-relaxation-manifest.json'
PAPER = [root/'paper/manuscript.md', root/'paper/main.tex']
draft = root/'docs/white-dwarf-free-fall-charge-section.md'
notes = [root/'notes/REQUEST288_RELAXATION_ESTIMATE_SECTION46_KO.md']
DOCS = ['model-definition', 'observable-targets', 'adiabatic-limit', 'nonadiabatic-regime', 'failure-ledger-dynamic-chi', 'dynamic-charge-completion']
SCRIPTS = ['phase288-apply.py', 'phase288-refine.py', 'phase288-publish.py']
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
    preserved = {p: h for p, h in read(previous_manifest)['sha256'].items() if not p.startswith('docs/')}
    for p, h in preserved.items(): assert sha(root/p) == h, p
    for p in PAPER + [draft]: assert hashlib.sha256(git('show', 'HEAD:' + p.relative_to(root).as_posix())).hexdigest() != sha(p), p
    lines = (root/'paper/manuscript.md').read_bytes().split(b'\n'); assert len([l for l in lines[:-1] if not l.endswith(b'\r')]) == 1
    for name in SCRIPTS: copy(S/name, out/'scripts'/name)
    now = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    final = dict(classification='Imported from prior work', passed=True, verdict='SECTION_4_6_LAYER_BY_LAYER_RELAXATION_ESTIMATE_ADDED__SINGLE_SELF_REVIEW',
        review='single self-review by user decision (three reviewers judged excessive for a sentence-level change)',
        edits=['Deeper layers bullet: layer-by-layer estimate (depths ~2,600 and 6,600 km, mass fractions ~2e-9 and 1.4e-7, Debye-weighted lag 4.0e-9 and 3.3e-7 of |S_struct|); Imported, with Conjectural model limits',
               'summary: bounded under the two assumptions and estimated layer by layer; no non-adiabatic calculation',
               'assumptions list: for the lag, either the two assumptions or the layer-by-layer estimate',
               'data availability: notes through REQUEST288; native-thermal-relaxation manifest listed'],
        unchanged='abstract, Section 6 paragraph, Lag bound and Structural bullets, Cassini paragraph (still true)',
        checks=['numbers asserted against relax.json by the apply script', 'Pandoc: 10 lists, 4 equations, all labels referenced, no missing citation keys',
                'manuscript CRLF kept (single original LF line kept)', 'main.tex regenerated with Pandoc 3.11'],
        pdf_and_submission_archive='predate Section 4.6 (no TeX on this host)', full_goal_complete=False, snapshot_KST=now)
    doc = ("\n\n## 단계288 — §4.6에 층별 열 완화 추정 반영(단독 검토)\n\n"
        "분류: Imported from prior work. 사용자 결정(2026-09-28)으로 문장 수준 변경은 외부 심사 없이 단독 검토로 반영했다. "
        "통합 원고 §4.6의 \"Deeper layers\" 항목, 요약, 가정 목록, 데이터 가용성에 단계287의 층별 열 완화 추정을 넣었다. "
        "추정 지연은 |𝒮_struct|의 4.0e−9(내측)와 3.3e−7(외측)이고, 계산 결과는 Imported, 모형 한계는 Conjectural로 표지했다. "
        "두 가정 상한, 초록, §6 문장은 여전히 참이라 그대로 두었다. main.tex는 Pandoc 3.11로 다시 만들었고 논문 검증을 통과한다. "
        "[근거 288](../notes/REQUEST288_RELAXATION_ESTIMATE_SECTION46_KO.md).\n")
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
    m['native_section46_relaxation'] = {k: v for k, v in final.items() if k not in ('sha256', 'document_prefixes')}
    master.write_text(json.dumps(m, indent=indent, ensure_ascii=False) + '\n', encoding='utf-8', newline='\n')


def check(mode):
    m = read(manifest); a = read(master)
    for p, h in m['sha256'].items(): assert sha(root/p) == h and a['sha256'][p] == h, p
    for p, h in read(out/'publication.json')['previous_nondoc'].items(): assert sha(root/p) == h, p
    for p, v in m['document_prefixes'].items(): assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest() == v['sha256'], p
    assert a['sha256'][manifest.relative_to(root).as_posix()] == sha(manifest)
    paths = list(m['sha256']) + [manifest.relative_to(root).as_posix(), 'paper/revision-manifest.json']
    (S/'phase288-paths.txt').write_bytes(b'\0'.join(p.encode() for p in paths) + b'\0')
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
