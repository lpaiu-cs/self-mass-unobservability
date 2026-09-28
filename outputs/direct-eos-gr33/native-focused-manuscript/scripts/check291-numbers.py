"""Phase 291: every number in the focused manuscript must occur in the previous full manuscript (no new quantitative claims)."""
import re, subprocess, sys
from pathlib import Path
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
old = subprocess.check_output(['git', 'show', '39962a753:paper/manuscript.md'], cwd=root).decode('utf-8')
new = Path(sys.argv[1]).read_text(encoding='utf-8')
tok = lambda s: re.findall(r'(?<![\w.])\d+(?:[.,]\d+)*(?![\w])', s)
oldset = set(tok(old))
miss = sorted({t for t in tok(new) if t not in oldset}, key=lambda t: new.index(t))
print('numbers', len(set(tok(new))), 'not in previous manuscript:', miss)
