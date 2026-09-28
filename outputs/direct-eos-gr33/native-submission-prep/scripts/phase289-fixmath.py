"""Phase 289: move the digit inside the inline math so Pandoc keeps it as math (closing $ followed by a digit is not math)."""
from pathlib import Path
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
old, new = r"(K\(\approx\)934)", r"(\(K\approx934\))"
for p in ('docs/white-dwarf-free-fall-charge-section.md', 'paper/manuscript.md'):
    t = (root/p).read_bytes().decode('utf-8'); assert t.count(old) == 1, p
    (root/p).write_bytes(t.replace(old, new).encode('utf-8')); print(p, 'ok')
