"""Phase 286: integrate the revision-6 white-dwarf section (Section 4.6) and its companion edits into the unified manuscript.

Byte-level edit of paper/manuscript.md that keeps its CRLF line endings (and its single bare-LF line untouched), restores the
three bibliography entries exactly as in commit a322c6462, and adds the PDF-predates note to paper/README.md (LF).
main.tex is regenerated afterwards with paper/build_manuscript.py (Pandoc 3.11). Run once on a clean tree.
"""
from pathlib import Path
import subprocess
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
draft = (root/'docs/white-dwarf-free-fall-charge-section.md').read_text(encoding='utf-8')
ms_path, bib_path, readme_path = root/'paper/manuscript.md', root/'paper/references.bib', root/'paper/README.md'


def between(text, start, end):
    i = text.index(start) + len(start); return text[i:text.index(end, i)]


section = draft[draft.index('### 4.6'):draft.index('\n---\n\n## Proposed companion edits')].strip('\n')
abstract_add = between(draft, '**Abstract, after the damped scalar-charge sentence:**\n', '\n\n').strip()
discussion_add = between(draft, '**Section 6, after the paragraph on the physical target:**\n', '\n\n').strip()
availability_add = between(draft, '**Data and code availability, new sentences:**\n', '\n\n**References:**').strip('\n')
for s in (abstract_add, discussion_add): assert s.startswith('**Conjectural.**') and '\n' not in s
assert availability_add.startswith('The white-dwarf calculation') and '\n- endpoint readout' in availability_add


def crlf(block): return [l.encode('utf-8') + b'\r' for l in block.split('\n')]


raw = ms_path.read_bytes(); lines = raw.split(b'\n'); n0 = len(lines)
bare = [i for i, l in enumerate(lines[:-1]) if not l.endswith(b'\r')]
assert len(bare) == 1 and not any(b'4.6 A worked white-dwarf' in l for l in lines)
find = lambda prefix: [i for i, l in enumerate(lines) if l.startswith(prefix.encode('utf-8'))]

# date
(i,) = find('**Date:** 9 September 2026'); lines[i] = b'**Date:** 28 September 2026\r'
# abstract: after the damped scalar-charge sentence
(i,) = find('**Proven.** We separate finite-order'); anchor = 'controlled one-pole reduction.'
s = lines[i].decode('utf-8'); assert s.count(anchor) == 1
lines[i] = s.replace(anchor, anchor + ' ' + abstract_add).encode('utf-8')
# Section 4.6 before Section 5
(i,) = find('## 5. Conditional application to PSR J0337+1715'); assert lines[i - 1] == b'\r'
lines[i:i] = crlf(section) + [b'\r']
# Section 6: after the physical-target paragraph
(i,) = find('**Counterexample candidate.** The physical target is a shared transfer relation'); assert lines[i + 1] == b'\r'
lines[i + 2:i + 2] = crlf(discussion_add) + [b'\r']
# data availability: before the reproduction-commands paragraph
(i,) = find('Reproduction commands and their scope are in'); assert lines[i - 1] == b'\r'
lines[i:i] = crlf(availability_add) + [b'\r']
out = b'\n'.join(lines)
assert [j for j, l in enumerate(out.split(b'\n')[:-1]) if not l.endswith(b'\r')].__len__() == 1
ms_path.write_bytes(out)

# bibliography: exactly the entries restored in a322c6462
old_bib = subprocess.check_output(['git', 'show', 'a322c6462^:paper/references.bib'], cwd=root)
new_bib = subprocess.check_output(['git', 'show', 'a322c6462:paper/references.bib'], cwd=root)
assert bib_path.read_bytes() == old_bib and new_bib.startswith(old_bib)
bib_path.write_bytes(new_bib)

# README (LF): PDF and archive predate Section 4.6
readme = readme_path.read_text(encoding='utf-8'); anchor = '- `../output/pdf/free-fall-identifiability.pdf` is the reviewed rendered output.\n'
assert readme.count(anchor) == 1 and 'Section 4.6' not in readme
note = ('- Section 4.6 (the worked white-dwarf calculation) was integrated on 28 September 2026 after independent review, with Pandoc 3.11, '
        'which reproduces the previous `main.tex` byte for byte from the previous source. The rendered PDF and the submission archive below '
        'predate Section 4.6; rebuilding them requires a TeX distribution, which is not installed on the current host.\n')
readme_path.write_text(readme.replace(anchor, anchor + note), encoding='utf-8', newline='\n')
print(dict(manuscript_lines=(n0, len(out.split(b'\n'))), section_lines=section.count('\n') + 1, bib_added=len(new_bib) - len(old_bib)))
