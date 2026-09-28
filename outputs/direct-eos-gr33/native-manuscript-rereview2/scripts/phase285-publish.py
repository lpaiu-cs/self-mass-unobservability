"""Archive phase 285: the re-review of revision 4 of the white-dwarf section, revision 5 of the draft, further errata and the
reproduction map.

Binds the response note, the revision-5 draft and the re-review records; appends one section to the six dynamic-chi documents;
binds everything in the manifest and paper/revision-manifest.json. No new computation.
Usage: python phase285-publish.py package|check|staged|head
"""
from pathlib import Path
import datetime, hashlib, json, subprocess, sys
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
S = Path(__file__).resolve().parent
master = root/'paper/revision-manifest.json'; previous_manifest = root/'outputs/direct-eos-gr33/native-manuscript-rereview-manifest.json'
out = root/'outputs/direct-eos-gr33/native-manuscript-rereview2'; manifest = out.parent/'native-manuscript-rereview2-manifest.json'
draft = root/'docs/white-dwarf-free-fall-charge-section.md'
notes = [root/'notes/REQUEST285_REVISION5_RESPONSE_KO.md']
DOCS = ['model-definition', 'observable-targets', 'adiabatic-limit', 'nonadiabatic-regime', 'failure-ledger-dynamic-chi', 'dynamic-charge-completion']
SCRIPTS = ['phase285-publish.py', 'scan285-inputs.py']
REVIEWS = {'rereview-prompt.md': 'rereview284-prompt.md', 'rereview-gpt-6-astra-prompt.md': 'rereview284-astra-prompt.md',
           'rereview-gpt-6-astra.md': 'rereview284-astra.md', 'rereview-gpt-6-astra-cli.log': 'rereview284-astra.log',
           'rereview-opus-5.5.md': 'rereview284-opus.md', 'rereview-fable-5.1.md': 'rereview284-fable.md'}
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
    for name in REVIEWS.values(): assert (S/name).is_file(), name
    preserved = {p: h for p, h in read(previous_manifest)['sha256'].items() if not p.startswith('docs/')}
    for p, h in preserved.items(): assert sha(root/p) == h, p
    rev4 = hashlib.sha256(git('show', 'HEAD:docs/white-dwarf-free-fall-charge-section.md')).hexdigest()
    assert read(master)['sha256'][draft.relative_to(root).as_posix()] == rev4 and sha(draft) != rev4
    out.mkdir(parents=True)
    for name in SCRIPTS: copy(S/name, out/'scripts'/name)
    for dst, name in REVIEWS.items(): copy(S/name, out/dst)
    now = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    final = dict(classification='Imported from prior work', passed=True, verdict='REVISION_4_REREVIEWED__REVISION_5_DRAFTED__ERRATA_EXTENDED__CONFIRMATION_REVIEW_PENDING',
        rereview_of_revision_4={'gpt-6-astra (Codex CLI, read-only)': 'major revision; blocking item lifted; errata no-go and necessary condition wrong (algebraic counterexample)',
                                'opus5.5 (subagent)': 'accept after minor revision; sentence-level majors: instantaneous-response attribution, strength assumption in summary and abstract, errata scope',
                                'fable5.1 (subagent)': 'accept after minor revision; one sentence-level major: instantaneous-response attribution and scale'},
        revision_5=['instantaneous response: zero-lag coefficient omitted by the fixed-companion reduction; up to about 1.24e-9 a_p^2, reaching Section 5 scales for |a_p| >~ 0.5; timing response and stiffness condition not evaluated',
                    'both relaxation assumptions in abstract and summary; lag bound limits size, does not remove the lag',
                    'Q_c = -Psi_out/m_wd vs Q = Q_c - a_i0 dm/m; window formula and whole-star history are Q_c',
                    'fixed-kernel condition on the 99.87% sign statement; Born wording; a_p = Q_p/m_p; modulation order; permanent zero only as t -> infinity',
                    'data availability: manifest file names, superseded verdict fields, reproduction map in REQUEST285'],
        errata=['REQUEST284 errata 4 and docs phase-284 sections: no-go under two assumptions wrong',
                'REQUEST284 and docs: instantaneous response attributed to the Section 4.3 stiffness condition',
                'REQUEST272 lines 8 and 19', 'REQUEST274 lines 48, 54, 58 (4.5e-18 with the outer term), 60', 'REQUEST276 lines 37, 39, 43, 47',
                'REQUEST277 and REQUEST280 closure labels', 'REQUEST278 cause: preceding stored samples', 'manifest verdict fields superseded',
                'REQUEST284 self-coupling sum 1.3e-7 may double count'],
        repository_classification='conditional theorem progress on the lag amplitude within the chosen response family; no observational exclusion; A4 not established; deep thermal relaxation remains an uncomputed loophole candidate',
        previous_draft_sha256=rev4, full_goal_complete=False, snapshot_KST=now)
    doc = ("\n\n## 단계285 — 4판 재심 응답(5판)과 분류 정정\n\n"
        "분류: Imported from prior work. 개정 4판의 재심에서 gpt-6-astra는 차단을 해제했지만 주요 수정을 권고했다. opus5.5와 fable5.1은 경미 수정 후 수락이었다. "
        "5판(`docs/white-dwarf-free-fall-charge-section.md`)은 다음을 고쳤다. 순간 응답 β_sδφ_mod는 Section 3의 지연 없는 계수로 적었다. 이 항은 §4.3의 고정 동반성 축약이 빠뜨리며, 펄서–내측 쌍에서 약 1.24e−9 a_p² 이하다. "
        "초록과 요약에는 두 완화 가정을 모두 적었다. 판독은 compact 부분 𝒬_c와 질량 정규화 𝒬로 나눴다. 새 계산은 없다.\n\n"
        "분류: Conjectural. 분류 정정: 단계284 항목의 '두 가정 아래 no-go'는 틀렸다. 두 가정은 지연의 크기(≤|𝒮_struct|/2)를 묶을 뿐 지연을 없애지 않는다. "
        "반례는 H=s/2+(s/2)/(1+iωτ)로, 두 가정을 만족하면서 궤도 진동수의 직교 성분이 s/4다. "
        "분류는 '선택한 응답족의 지연 진폭 상한에 관한 조건부 theorem progress'이며, 관측 배제나 A4 유지는 성립하지 않는다. "
        "no-go가 무너지는 정확한 단계는 열 완화 세기다. 단열 정적 계산은 이를 묶지 못한다. "
        "빠진 최소 계산은 두 가지다: 깊은 층의 비단열 열 응답(세기와 극점 구조), 두 비공통 채널의 타이밍 응답. 깊은 층의 열 완화 상태는 세기를 계산하지 않은 loophole 후보로 남는다. "
        "단계272–284 항목의 부호·no-go·A4 진술은 [근거 285](../notes/REQUEST285_REVISION5_RESPONSE_KO.md)의 정오표로 대체된다. 같은 세 심사자가 반영 여부를 확인한다.\n")
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
    m['native_section_revision_5'] = dict(date='2026-09-27', draft='docs/white-dwarf-free-fall-charge-section.md (revision 5)',
        response='notes/REQUEST285_REVISION5_RESPONSE_KO.md', status='pending confirmation review by gpt-6-astra, opus5.5, fable5.1',
        rereview_of_revision_4=final['rereview_of_revision_4'], superseded_verdict_fields='see REQUEST285 errata 14', previous_draft_sha256=rev4)
    master.write_text(json.dumps(m, indent=indent, ensure_ascii=False) + '\n', encoding='utf-8', newline='\n')


def check(mode):
    m = read(manifest); a = read(master)
    for p, h in m['sha256'].items(): assert sha(root/p) == h and a['sha256'][p] == h, p
    for p, h in read(out/'publication.json')['previous_nondoc'].items(): assert sha(root/p) == h, p
    for p, v in m['document_prefixes'].items(): assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest() == v['sha256'], p
    assert a['sha256'][manifest.relative_to(root).as_posix()] == sha(manifest)
    paths = list(m['sha256']) + [manifest.relative_to(root).as_posix(), 'paper/revision-manifest.json']
    (S/'phase285-paths.txt').write_bytes(b'\0'.join(p.encode() for p in paths) + b'\0')
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
