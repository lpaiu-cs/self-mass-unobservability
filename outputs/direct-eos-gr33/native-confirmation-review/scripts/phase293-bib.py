"""Phase 293: add the two REML/timing-noise references (verified on Crossref); CRLF kept."""
from pathlib import Path
p = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672/paper/references.bib')
t = p.read_bytes().decode('utf-8')
assert 'patterson1971reml' not in t and 'vanhaasteren2013noise' not in t
new = ("\r\n@article{patterson1971reml,\r\n"
       "  author  = {Patterson, H. D. and Thompson, R.},\r\n"
       "  title   = {Recovery of inter-block information when block sizes are unequal},\r\n"
       "  journal = {Biometrika},\r\n  volume  = {58},\r\n  pages   = {545--554},\r\n  year    = {1971},\r\n"
       "  doi     = {10.1093/biomet/58.3.545}\r\n}\r\n"
       "\r\n@article{vanhaasteren2013noise,\r\n"
       "  author  = {van Haasteren, Rutger and Levin, Yuri},\r\n"
       "  title   = {Understanding and analysing time-correlated stochastic signals in pulsar timing},\r\n"
       "  journal = {Monthly Notices of the Royal Astronomical Society},\r\n  volume  = {428},\r\n  pages   = {1147--1159},\r\n  year    = {2013},\r\n"
       "  doi     = {10.1093/mnras/sts097}\r\n}\r\n")
p.write_bytes((t.rstrip('\r\n') + '\r\n' + new).encode('utf-8'))
b = p.read_bytes(); print('entries', b.count(b'\n@'), 'bare LF', b.count(b'\n') - b.count(b'\r\n'))
