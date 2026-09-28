"""Phase 291: narrow the paper's focus.

paper/manuscript.md becomes the focused paper (scratchpad manuscript-focused.md, CRLF). The previous full manuscript
(39962a753:paper/manuscript.md) moves unchanged into paper/supplement.md as the Supplemental Material; only its title and abstract
are replaced by an SM title and a description of what it contains.
"""
import subprocess
from pathlib import Path
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
S = Path(__file__).resolve().parent
old = subprocess.check_output(['git', 'show', '39962a753:paper/manuscript.md'], cwd=root).decode('utf-8')
new_title = 'Identifying a relaxing internal state in free-fall timing: finite-frequency boundaries and an application to PSR J0337+1715'

# Supplement: previous full text, new title and abstract.
title_line = old.split('\r\n', 1)[0]
assert title_line == '# Static response and dynamical identifiability in free-fall tests: finite-order boundaries and a pulsar-triple application'
head, rest = old.split('## Abstract\r\n\r\n', 1)
abstract_old, tail = rest.split('\r\n\r\n---\r\n', 1)
assert '\r\n' not in abstract_old and abstract_old.startswith('**Proven.** We separate')
sm_abstract = ("This Supplemental Material (SM) is the complete technical account of the paper named in the title. "
               "It is the full earlier version of the work; apart from this title and abstract, its text is unchanged. "
               "It contains the finite static operator catalog (Section 2 and Appendices A--C), the nuisance-projected rank and the pair-coupling benchmark (Sections 4.1--4.2), "
               "the complete scalar-charge matching (Sections 4.3--4.5), the worked calculation for the inner white dwarf of PSR J0337+1715 (Section 4.6), "
               "and all audits of the stored timing analysis with their provenance and reproduction commands (Sections 5.1--5.10 and Appendices D--E). "
               "Its section, equation, table and figure numbers are its own; Theorem 3 here is Theorem 1 of the main text. "
               "The status labels are defined in Section 1.")
sm = (head.replace(title_line, f"# Supplemental Material for ``{new_title}''", 1)
      + '## Abstract\r\n\r\n' + sm_abstract + '\r\n\r\n---\r\n' + tail)
assert sm.count('\r\n') == old.count('\r\n') and sm.split('---\r\n', 1)[1] == old.split('---\r\n', 1)[1]
(root/'paper/supplement.md').write_bytes(sm.encode('utf-8'))

# Focused paper, CRLF.
focused = (S/'manuscript-focused.md').read_text(encoding='utf-8').replace('\r\n', '\n')
assert focused.startswith('# ' + new_title + '\n')
(root/'paper/manuscript.md').write_bytes(focused.replace('\n', '\r\n').encode('utf-8'))
print('supplement bytes', len(sm.encode()), 'manuscript bytes', (root/'paper/manuscript.md').stat().st_size)
