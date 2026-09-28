"""Phase 290: final pre-submission fixes (manuscript CRLF, draft LF, references CRLF with a few LF lines kept).

Manuscript: American spelling (center, meters), inline math for c_chi and pi, Nançay, en dashes in ranges, drop an internal ledger
reference, 'the figure' -> 'Figure 1', and the data-availability paragraph for the public snapshot (author decision 2026-09-28).
References: remove the internal note that printed in the bibliography; add verified DOIs (Crossref) to two entries.
"""
from pathlib import Path
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')

def edit(rel, pairs):
    p = root/rel; t = p.read_bytes().decode('utf-8')
    for old, new, n in pairs:
        assert t.count(old) == n, (rel, old[:60], t.count(old)); t = t.replace(old, new)
    p.write_bytes(t.encode('utf-8')); print(rel, 'ok')

centre = [('through the centre, and', 'through the center, and', 1), ('reaches the centre at', 'reaches the center at', 1), ('from the centre arrival', 'from the center arrival', 1)]
edit('paper/manuscript.md', centre + [
    (') metres. In this convention', ') meters. In this convention', 1),
    ('can remain when c_chi is nonzero.', r'can remain when \(c_\chi\) is nonzero.', 1),
    ('adds a phase of pi.', r'adds a phase of \(\pi\).', 1),
    ('1.67084424 and pi radians.', r'1.67084424 and \(\pi\) radians.', 1),
    ('cosine and plus-pi/2 drives at', r'cosine and plus-\(\pi/2\) drives at', 1),
    ('The plus-pi/2 column', r'The plus-\(\pi/2\) column', 1),
    ('public Nancay pulse times', 'public Nançay pulse times', 1),
    ('days in 2013-2021, a published', 'days in 2013--2021, a published', 1),
    ('by about 0.52-0.60 percent.', 'by about 0.52--0.60 percent.', 1),
    ('the factor 2.1-17.4 full/truncated', 'the factor 2.1--17.4 full/truncated', 1),
    (' The potential premise is A9 in the earlier repository ledger.', '', 1),
    ('Table 1 and the figure read', 'Table 1 and Figure 1 read', 1),
    (r"Manuscript source, symbolic checks and stored artifacts are maintained at \url{https://github.com/lpaiu-cs/self-mass-unobservability}. "
     r"The accompanying \path{paper/revision-manifest.json} identifies this revision's inputs by SHA-256. "
     r"The unified input state is commit `4897038`; the ordered follow-through designs and results are recorded in the revision manifest. "
     r"Remote availability must be checked at submission; this local revision does not claim a newly deposited DOI archive.",
     r"The manuscript source, code, notes, manifests and numerical results are available in a public snapshot at \url{https://github.com/lpaiu-cs/self-mass-unobservability}. "
     r"The accompanying \path{paper/revision-manifest.json} identifies this revision's inputs by SHA-256. "
     r"Binary runtime arrays, about 50 GB in total, exceed the hosting limits and are not deposited; the manifests identify them by SHA-256, and they are available from the author on request. "
     r"Commit identifiers cited in this paper, such as the unified input state `4897038`, refer to the author's full repository history, of which the public repository holds a snapshot; "
     r"\path{PUBLIC_SNAPSHOT.md} there states what the snapshot omits. The ordered follow-through designs and results are recorded in the revision manifest.", 1),
])
edit('docs/white-dwarf-free-fall-charge-section.md', centre)
edit('paper/references.bib', [
    ("  archivePrefix = {arXiv},\r\n  note    = {Source of the $|\\Delta| < 1.5\\times10^{-6}$ (planet hypothesis) and $|\\Delta| < 2.3\\times10^{-6}$ (red-noise hypothesis) 95\\% CL inputs}\r\n}",
     "  archivePrefix = {arXiv},\r\n  doi     = {10.1051/0004-6361/202452100}\r\n}", 1),
    ("  volume  = {9},\r\n  pages   = {2093},\r\n  year    = {1992}\r\n}",
     "  volume  = {9},\r\n  pages   = {2093--2176},\r\n  year    = {1992},\r\n  doi     = {10.1088/0264-9381/9/9/015}\r\n}", 1),
])
