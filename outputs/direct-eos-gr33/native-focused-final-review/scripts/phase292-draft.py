"""Phase 292: drop the withdrawn physical-drive endpoint scale from the white-dwarf section draft (bytes kept otherwise)."""
from pathlib import Path
p = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672/docs/white-dwarf-free-fall-charge-section.md')
b = p.read_bytes()
old = rb'- \(4.1\times10^{-10}\) (physical-drive endpoint)' + b'\n'
assert b.count(old) == 1
p.write_bytes(b.replace(old, b''))
print('removed', len(old), 'bytes')
