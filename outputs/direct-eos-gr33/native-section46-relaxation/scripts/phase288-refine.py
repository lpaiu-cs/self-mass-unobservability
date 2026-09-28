"""Phase 288 refinement: define the shell strength in the Deeper-layers bullet (draft and manuscript alike)."""
from pathlib import Path
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
old = (r"Each layer relaxes on the thermal time of the layers above it, its fully relaxed limit is bracketed by the isothermal exponent \(P_{\rm gas}/P\), and the static structural response is solved again with the faster layers relaxed. "
       r"The layers whose thermal time equals \(1/\omega\) lie about 2,600 km deep for the inner orbit and 6,600 km for the outer one, above about \(2\times10^{-9}\) and \(1.4\times10^{-7}\) of the mass. "
       r"Summing the shell strengths, of either sign, with the Debye weight")
new = (r"Each layer relaxes on the thermal time of the layers above it, and its fully relaxed limit is bracketed by the isothermal exponent \(P_{\rm gas}/P\). "
       r"Solving the static structural response again with all layers faster than a given time relaxed gives the relaxation strength of each deeper shell. "
       r"The layers whose thermal time equals \(1/\omega\) lie about 2,600 km deep for the inner orbit and 6,600 km for the outer one, above about \(2\times10^{-9}\) and \(1.4\times10^{-7}\) of the mass. "
       r"Summing the absolute shell strengths with the Debye weight")
texts = {p: (root/p).read_bytes().decode('utf-8') for p in ('docs/white-dwarf-free-fall-charge-section.md', 'paper/manuscript.md')}
for p, t in texts.items(): assert t.count(old) == 1, p
for p, t in texts.items(): (root/p).write_bytes(t.replace(old, new).encode('utf-8')); print(p, 'ok')
