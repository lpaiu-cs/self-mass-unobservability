"""Archive phase 282: the three independent reviews of the white-dwarf charge section, the revert of its manuscript integration,
and the revision plan.

Binds the synthesis note and the revised English section draft, appends one section to the six dynamic-chi documents, and binds
everything in the manifest and paper/revision-manifest.json. No new computation: values are read from the phase 276-279 manifests.
Usage: python phase280-publish.py package|check|staged|head
"""
from pathlib import Path
import datetime, hashlib, json, subprocess, sys
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
S = Path(__file__).resolve().parent
master = root/'paper/revision-manifest.json'; previous_manifest = root/'outputs/direct-eos-gr33/native-conditions-closed-manifest.json'
out = root/'outputs/direct-eos-gr33/native-manuscript-review'; manifest = out.parent/'native-manuscript-review-manifest.json'
notes = [root/'notes/REQUEST282_INDEPENDENT_REVIEW_REVISION_KO.md']
DOCS = ['model-definition', 'observable-targets', 'adiabatic-limit', 'nonadiabatic-regime', 'failure-ledger-dynamic-chi', 'dynamic-charge-completion']
SCRIPTS = ['phase282-publish.py']
REVIEWS = {'review-prompt.md': 'review281-prompt.md', 'review-gpt-6-astra.md': 'review281-astra.md', 'review-opus-5.5.md': 'review281-opus.md', 'review-fable-5.1.md': 'review281-fable.md', 'review-gpt-6-astra-cli.log': 'review281-astra.log'}
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
    for name in REVIEWS.values(): assert (S/name).is_file(), name
    out.mkdir(parents=True)
    for name in SCRIPTS: copy(S/name, out/'scripts'/name)
    for dst, name in REVIEWS.items(): copy(S/name, out/dst)
    now = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    final = dict(classification='Imported from prior work', passed=True, verdict='THREE_INDEPENDENT_REVIEWS_MAJOR_REVISION__INTEGRATION_REVERTED__REVISION_PLANNED',
        reviews={'gpt-6-astra (Codex CLI 0.154.0, read-only)': 'major revision; one blocking item (Section 5 coefficient interval used as white-dwarf observational sensitivity)',
                 'opus5.5 (subagent)': 'major revision then re-review; no blocking item', 'fable5.1 (subagent)': 'major revision; no blocking item'},
        integration_commit='a322c6462', revert_commit='b40864a22', verify_after_revert='PASS',
        consensus=['J0337 comparison: coefficient envelope vs amplitude, two pair channels, scale range, Cassini rescaling, drop NS-sector conclusion',
                   'claim labels per manuscript convention; define A4; per-item labels; soften abstract and discussion',
                   'history in readout time; static-window sign; kernel-weighted 0.25% (max 0.38%)',
                   'tidal: selection rule and equilibrium-tide bound only; dynamical tide dissipative; damping-depth formula outside validity',
                   'symbol collisions with Sections 3-5; define local notation'],
        full_goal_complete=False, snapshot_KST=now)
    doc = ("\n\n## 단계282 — 독립 심사와 통합 되돌림\n\n"
        "분류: Imported from prior work. 통합 원고 §4.6(커밋 a322c6462)을 gpt-6-astra(Codex, 읽기 전용), opus5.5, fable5.1이 서로 모른 채 심사했다. 셋 모두 주요 수정을 권고했다(astra는 차단 1건). "
        "수치는 기록과 일치했고, 핵심 결론을 뒤집는 결함은 없었다. 합의된 지적은 다음과 같다: J0337 비교(계수 구간 대 진폭, 두 쌍 채널, 척도 범위, Cassini 재척도), 주장 표지, 판독 시각 기준 이력, 조석 서술, 기호 충돌. "
        "사용자 결정에 따라 통합을 커밋 b40864a22로 되돌렸고 논문 검증은 통과한다. 수정 뒤 같은 세 심사자로 재심한다. "
        "[근거 282](../notes/REQUEST282_INDEPENDENT_REVIEW_REVISION_KO.md).\n")
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
    m['native_manuscript_review'] = {k_: v for k_, v in final.items() if k_ != 'sha256'}
    master.write_text(json.dumps(m, indent=indent, ensure_ascii=False) + '\n', encoding='utf-8', newline='\n')


def git(*args): return subprocess.check_output(['git', *args], cwd=root)


def check(mode):
    m = read(manifest); a = read(master)
    for p, h in m['sha256'].items(): assert sha(root/p) == h and a['sha256'][p] == h, p
    for p, h in read(out/'publication.json')['previous_nondoc'].items(): assert sha(root/p) == h, p
    for p, v in m['document_prefixes'].items(): assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest() == v['sha256'], p
    assert a['sha256'][manifest.relative_to(root).as_posix()] == sha(manifest)
    paths = list(m['sha256']) + [manifest.relative_to(root).as_posix(), 'paper/revision-manifest.json']
    (S/'phase282-paths.txt').write_bytes(b'\0'.join(p.encode() for p in paths) + b'\0')
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
