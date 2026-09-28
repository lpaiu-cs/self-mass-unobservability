"""Phase 291: citations and dangling references in the generated TeX files."""
import re
root = 'E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672/'
for f in ('paper/main.tex', 'paper/supplement.tex'):
    m = open(root + f, encoding='utf-8').read()
    cites = sorted({k.strip() for c in re.findall(r'\\cite[pt]?\{([^}]*)\}', m) for k in c.split(',')})
    dangling = sorted(set(re.findall(r'\\ref\{([^}]*)\}', m)) - set(re.findall(r'\\label\{([^}]*)\}', m)))
    print(f, 'cites', len(cites), cites, 'dangling refs', dangling, 'escaped dollars', m.count('\\$'))
