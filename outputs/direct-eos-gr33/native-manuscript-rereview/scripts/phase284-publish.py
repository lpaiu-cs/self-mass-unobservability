"""Archive phase 284: the re-review of revision 3 of the white-dwarf section, revision 4 of the draft, and errata to earlier notes.

Binds the response note (with errata), the previously unbound integration-plan note REQUEST281, the revision-4 draft and the
re-review records; appends one section to the six dynamic-chi documents; binds everything in the manifest and
paper/revision-manifest.json. No new computation: values are read from stored manifests and checked by short arithmetic.
Usage: python phase284-publish.py package|check|staged|head
"""
from pathlib import Path
import datetime, hashlib, json, subprocess, sys
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
S = Path(__file__).resolve().parent
master = root/'paper/revision-manifest.json'; previous_manifest = root/'outputs/direct-eos-gr33/native-manuscript-review-manifest.json'
out = root/'outputs/direct-eos-gr33/native-manuscript-rereview'; manifest = out.parent/'native-manuscript-rereview-manifest.json'
draft = root/'docs/white-dwarf-free-fall-charge-section.md'
notes = [root/'notes/REQUEST281_MANUSCRIPT_INTEGRATION_KO.md', root/'notes/REQUEST283_REVISION_RESPONSE_KO.md', root/'notes/REQUEST284_REVISION4_RESPONSE_KO.md']
DOCS = ['model-definition', 'observable-targets', 'adiabatic-limit', 'nonadiabatic-regime', 'failure-ledger-dynamic-chi', 'dynamic-charge-completion']
SCRIPTS = ['phase284-publish.py', 'check284-pandoc.py']
REVIEWS = {'rereview-prompt.md': 'rereview283-prompt.md', 'rereview-gpt-6-astra-prompt.md': 'rereview283-astra-prompt.md',
           'rereview-gpt-6-astra.md': 'rereview283-astra.md', 'rereview-gpt-6-astra-cli.log': 'rereview283-astra.log',
           'rereview-opus-5.5.md': 'rereview283-opus.md', 'rereview-fable-5.1.md': 'rereview283-fable.md'}
REV3_SHA = 'cf5460324d19230beaee42f64863d43fc30a216e800df889f2425917fa0b5be4'
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
    for name in REVIEWS.values(): assert (S/name).is_file(), name
    preserved = {p: h for p, h in read(previous_manifest)['sha256'].items() if not p.startswith('docs/')}
    for p, h in preserved.items(): assert sha(root/p) == h, p
    assert read(master)['sha256'][draft.relative_to(root).as_posix()] == REV3_SHA and sha(draft) != REV3_SHA
    out.mkdir(parents=True)
    for name in SCRIPTS: copy(S/name, out/'scripts'/name)
    for dst, name in REVIEWS.items(): copy(S/name, out/dst)
    now = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    final = dict(classification='Imported from prior work', passed=True, verdict='REVISION_3_REREVIEWED__REVISION_4_DRAFTED__ERRATA_RECORDED__REREVIEW_PENDING',
        rereview_of_revision_3={'gpt-6-astra (Codex CLI, read-only)': 'major revision; blocking item kept (Section 5 sensitivity conclusion)',
                                'opus5.5 (subagent)': 'accept after minor revision; three sentence-level majors (instantaneous response, relaxation strength, sign sentence)',
                                'fable5.1 (subagent)': 'accept after minor revision'},
        revision_4=['sensitivity-exclusion conclusions removed; scale comparison only; timing response of the two non-common channels not computed',
                    'relaxation strength <= S_struct and single relaxation stated as uncomputed assumptions; A4 inference removed',
                    'instantaneous response beta_s*dphi_mod separated (no lag; Section 4.3 feedback-stiffness condition)',
                    'sign sentence corrected: no face-density change below 99.87% (max norm) reverses the sign',
                    'sources and norms: 3.5/4.0 ms and static window in the whole-star paragraph; Born comparison for orientation only; peak from the 0.2 ms grid',
                    'numbers: 2.13e-18, 7.53e-21, |S_struct|<=8.84e-9, non-gray 1.9-2.4 (1 sigma), -8.0e-41 at 0.23 s',
                    'notation, labels, blank lines before lists (Pandoc: 10 lists, 4 equations), companion edits incl. bib restore and PDF note'],
        errata=['REQUEST272 line 19 sign statement false', 'REQUEST278 history table values and time labels', 'REQUEST282 peak grid; core conclusion conditional',
                'REQUEST273/274/276/280 orbital-timescale no-go and A4 statements unproven: conditional theorem progress only',
                'REQUEST274 |kappa_struct|<=8.8e-9 rounded down (8.8366e-9)', 'REQUEST283 verifier table-row reason inaccurate'],
        repository_classification='conditional theorem progress (no-go only under the two stated relaxation assumptions; A4 remains an assumption)',
        previous_draft_sha256=REV3_SHA, full_goal_complete=False, snapshot_KST=now)
    doc = ("\n\n## 단계284 — 재심 응답(4판)과 결론 정정\n\n"
        "분류: Imported from prior work. 개정 3판의 재심에서 gpt-6-astra는 주요 수정을 권고했다(차단: §5 감도 결론). opus5.5와 fable5.1은 경미 수정 후 수락이었다. "
        "4판(`docs/white-dwarf-free-fall-charge-section.md`)은 다음을 고쳤다. 감도 배제 결론은 척도 비교로 낮췄다. 열 완화 세기 ≤ 𝒮_struct와 단일 완화는 계산하지 않은 가정으로 명시했다. "
        "순간 응답 β_sδφ_mod(지연 없음, §4.3 조건)를 분리했다. 부호 문장은 '최대노름 99.87% 미만의 변화는 부호를 바꾸지 못한다'로 고쳤다. 새 결합 계산은 없다.\n\n"
        "분류: Conjectural. 결론 정정: 궤도 시간척도에서 이 자유낙하 상태가 관측량을 만들지 않는다는 이전 진술과 A4가 유지된다는 진술(단계273·274·276·280)은 증명되지 않았다. "
        "성립하는 것은 두 가정 아래의 척도 비교다. 그 아래에서 구조 변조의 척도가 §5 저장 척도보다 8자릿수 이상 작다. 두 비공통 채널의 타이밍 감도는 계산하지 않았다. "
        "분류는 조건부 theorem progress다. 붕괴를 피하는 최소 추가 조건은 궤도 주기와 비슷한 완화 시간을 가진 상태가 있고, 그 완화 세기가 𝒮_struct를 넘거나 단일 완화가 아닌 것이다. "
        "이전 노트의 정오표(부호 문장, 단계278 이력 값, 반올림)는 [근거 284](../notes/REQUEST284_REVISION4_RESPONSE_KO.md)에 있다. 같은 세 심사자로 다시 재심한다.\n")
    prefixes = {}
    for name in DOCS:
        p = root/f'docs/{name}.md'; prefixes[p.relative_to(root).as_posix()] = dict(bytes=p.stat().st_size, sha256=sha(p))
        with p.open('ab') as h: h.write(doc.encode())
    write(out/'result.json', final); write(out/'publication.json', dict(previous_nondoc=preserved, document_prefixes=prefixes))
    bound = notes + [draft] + [root/p for p in prefixes] + [p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes, sha256={p.relative_to(root).as_posix(): sha(p) for p in bound})
    write(manifest, final)
    raw = master.read_text(encoding='utf-8'); m = json.loads(raw); indent = 2 if raw.startswith('{\n  "') else 1
    m['sha256'].update(final['sha256']); m['sha256'][manifest.relative_to(root).as_posix()] = sha(manifest)
    m['native_section_revision_4'] = dict(date='2026-09-27', draft='docs/white-dwarf-free-fall-charge-section.md (revision 4)',
        response='notes/REQUEST284_REVISION4_RESPONSE_KO.md', status='pending re-review by gpt-6-astra, opus5.5, fable5.1',
        rereview_of_revision_3=final['rereview_of_revision_3'], previous_draft_sha256=REV3_SHA)
    master.write_text(json.dumps(m, indent=indent, ensure_ascii=False) + '\n', encoding='utf-8', newline='\n')


def git(*args): return subprocess.check_output(['git', '-c', 'gc.auto=0', '-c', 'maintenance.auto=false', *args], cwd=root)


def check(mode):
    m = read(manifest); a = read(master)
    for p, h in m['sha256'].items(): assert sha(root/p) == h and a['sha256'][p] == h, p
    for p, h in read(out/'publication.json')['previous_nondoc'].items(): assert sha(root/p) == h, p
    for p, v in m['document_prefixes'].items(): assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest() == v['sha256'], p
    assert a['sha256'][manifest.relative_to(root).as_posix()] == sha(manifest)
    paths = list(m['sha256']) + [manifest.relative_to(root).as_posix(), 'paper/revision-manifest.json']
    (S/'phase284-paths.txt').write_bytes(b'\0'.join(p.encode() for p in paths) + b'\0')
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
