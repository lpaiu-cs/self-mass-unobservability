"""Convert the revision-4 section draft with the manuscript's Pandoc path and count lists, equations, citations and references."""
import os, re, sys
from pathlib import Path
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
S = Path(__file__).resolve().parent
sys.path.insert(0, str(root/'paper'))
os.environ.setdefault('PANDOC', 'C:/Users/lpaiu/AppData/Local/Pandoc/pandoc.exe')
import build_manuscript as b

t = (root/'docs/white-dwarf-free-fall-charge-section.md').read_text(encoding='utf-8')
sec = t[t.index('### 4.6'):t.index('\n---\n\n## Proposed')]
comp = t[t.index('## Proposed companion edits'):]
bib = (root/'paper/references.bib').read_text(encoding='utf-8')
for name, s in [('section', sec), ('companion', comp)]:
    tex = b.pandoc_fragment(b.preprocess_math(s))
    (S/f'rev4-{name}.tex').write_text(tex, encoding='utf-8')
    keys = sorted(set(k.strip() for c in re.findall(r'\\cite\{([^}]+)\}', tex) for k in c.split(',')))
    print(name, 'itemize', tex.count(r'\begin{itemize}'), 'items', tex.count(r'\item'), 'equations', tex.count(r'\begin{equation}'),
          'cite keys', keys, 'missing in bib', [k for k in keys if '{' + k + ',' not in bib])
labels = set(re.findall(r'\\label\{([^}]+)\}', sec)); refs = set(re.findall(r'\\ref\{([^}]+)\}', sec))
print('labels', sorted(labels), 'unreferenced', sorted(labels - refs), 'dangling', sorted(refs - labels))
