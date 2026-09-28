"""Phase 292: locate phase-dependent numbers and statements in the SM, main text, drafts and docs; summarize nuisance-intervals.csv."""
import csv, re
from pathlib import Path
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
needles = ['3.11837236', '1.69406454', '1.67084424', '0.12326821', '0.10004792', '0.08289', '589', '0.003035', '1.00096', '1.00523',
           '0.8680', '0.8785', '4.10217', '8.31672', r'4.1\times10^{-10}', 'eta as', r'includes \(\beta=0\)', 'physical-drive endpoint',
           '8.36', '3.19', '2.88', '-0.12326821', 'e\\cos']
for rel in ['paper/supplement.md', 'paper/manuscript.md', 'docs/white-dwarf-free-fall-charge-section.md', 'paper/README.md']:
    lines = (root/rel).read_text(encoding='utf-8').splitlines()
    for i, line in enumerate(lines, 1):
        hits = [n for n in needles if n in line]
        if hits: print(f'{rel}:{i}', hits, '|', line[:110])
rows = list(csv.DictReader(open(root/'outputs/research-completion/nuisance-intervals.csv', encoding='utf-8')))
print('csv rows', len(rows), 'cuts', sorted({r['cut'] for r in rows}), 'K', sorted({r['K'] for r in rows}), 'ranks', sorted({r['rank'] for r in rows}))
for r in rows:
    if r['K'] in ('1.0', '1') and r['rank'] in ('90', '71'): print(r['cut'], r['rank'], r['tau'], r['K'], r['U'])
