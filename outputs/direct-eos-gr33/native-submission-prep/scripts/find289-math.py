"""List escaped dollars in main.tex and inline-math closings followed by a digit in manuscript.md (Pandoc drops such math)."""
import re
from pathlib import Path
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
tex = (root/'paper/main.tex').read_text(encoding='utf-8').split('\n')
for i, l in enumerate(tex, 1):
    if re.search(r'(?<!\\)\\\$', l): print('main.tex', i, l[:150])
md = (root/'paper/manuscript.md').read_text(encoding='utf-8').split('\n')
for i, l in enumerate(md, 1):
    for m in re.finditer(r'\\\)[0-9]', l): print('manuscript.md', i, l[max(0, m.start() - 60):m.end() + 20])
