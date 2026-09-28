"""Reuse the same failed system; admit B columns through the existing ULP test."""
from pathlib import Path
source=Path('.phase215-polish-probe.py').read_text()
changes=[
    ("OUT=Path('native-broad-polish215-work')","OUT=Path('native-broad-polish215-work/basis')"),
    (".replace('if norm>goal*10000:break','if norm>goal*1000000:break')", ".replace('if norm>goal*10000:break','if norm>goal*1000000:break').replace('if value and col%4!=2','if value')"),
    ("sol=polish(m,op,rhs,z['solution'])", "sol=polish(m,op,rhs,np.load(OUT.parent/'proposal.npz')['solution'])"),
    ("Keep24columns,4passes and all full-vector/physical improvement requirements.", "Keep24columns,4passes and all full-vector/physical improvement requirements. Remove the blanket exclusion of B columns: the existing coefficient-times-state-ULP gate still filters every candidate. Reuse the prior polished proposal; do not fit a physical output."),
    ("Test actual improvement at the unchanged1e-14linear and1e-12nonlinear gates before any new evolution.","The broader trigger reduced the defect but left3.301e-13linear and3.788e-7actual failures. Test the quantum-filtered correction basis at unchanged1e-14linear and1e-12nonlinear gates before any new evolution.")]
for a,b in changes:assert source.count(a)==1,(a,source.count(a));source=source.replace(a,b)
exec(compile(source,__file__,'exec'),globals())
