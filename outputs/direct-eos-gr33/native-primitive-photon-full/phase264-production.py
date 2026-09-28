"""Run the unchanged phase-263 production after the completed phase-264 reference.

Only three things differ from phase 263: the admitted reference result, the
execution plan that binds this continuation, and the per-path wall cap (the
production cannot resume from partial snapshots, so a stop forces a full rerun).
"""
from pathlib import Path
source=Path('.phase263-production.py').read_text()
changes={
    "out/'direct-reference/result.json'":"out/'direct-reference-264/result.json'",
    "'deadline(14400)':'deadline(28800)'":"'deadline(14400)':'deadline(43200)'",
    "changes={\n":"changes={\n    \"b.read(out/'execution-plan.json')\":\"b.read(out/'execution-plan-264.json')\",\n",
}
for old,new in changes.items():
    assert source.count(old)==1,(old,source.count(old));source=source.replace(old,new)
exec(compile(source,__file__,'exec'))
