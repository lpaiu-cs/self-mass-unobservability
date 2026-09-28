"""Finish only the22 missing density columns after a preserved call-cap stop."""
from pathlib import Path
import inspect
import json
import time
import numpy as np
import def_native_energy_columns as prior


def main():
    out=prior.OUT;d=np.load(out/'temperature-columns-progress.npz');completed=len(d['offsets'])-1
    assert completed==125 and int(d['native_calls'])==5999
    assert not (out/'column-completion-plan.json').exists()
    prior.old.write(out/'column-budget-stop.json',dict(classification='Counterexample candidate',native_call_cap=6000,native_calls_spent=6000,
        completed_columns=125,total_columns=147,saved_native_calls=5999,unsaved_native_calls=1,
        exception='Native call budget',finished_columns_reused=True,original_failure_preserved=True))
    prior.old.write(out/'column-completion-plan.json',dict(classification='Counterexample candidate',
        reason='The bounded refinement completed125 of147 columns and saved7249 states; only22 specified dense columns remain. Reuse every finished column instead of restarting or increasing resolution.',
        additional_native_calls=3000,additional_seconds=60,total_column_native_cap=9000,
        forecast='Up to101 states per saved column suggests about2000-2500 further native calls; dense-column curvature is an assumption, capped at3000 calls and60seconds.',
        unchanged=['density grid','source temperature span','midpoint threshold','depth limit','independent controls','fluid grid','horizon'],
        stop='If this completion budget is insufficient, retain the incomplete bank and reassess. No automatic further expansion.',
        bindings={str(p.relative_to(prior.old.old.ROOT)):prior.old.old.photons.digest(p) for p in [Path(__file__),Path(prior.__file__),out/'temperature-columns-progress.npz']}))
    source=inspect.getsource(prior.main)
    replacements={
        'signal.alarm(160)':'signal.alarm(60)',
        "OUT/'column-reassessment.json'":"OUT/'column-completion-execution.json'",
        'additional_native_calls=6000,seconds=160':'additional_native_calls=3000,seconds=60',
        'call_cap=6000':'call_cap=3000',
        'columns=[];temps=[];offsets=[0];records=[]':"saved=np.load(OUT/'temperature-columns-progress.npz');columns=list(saved['raw']);temps=list(saved['logT']);offsets=list(saved['offsets']);records=[];finished=len(offsets)-1",
        "for j,xx in enumerate(d['x']):":"for j,xx in enumerate(d['x']):\n        if j<finished:continue",
        'native_calls=fan.calls,seconds=time.monotonic()-start,states=len(temps)':'native_calls=fan.calls,total_native_calls=6000+fan.calls,reused_finished_columns=125,seconds=time.monotonic()-start,states=len(temps)'}
    for before,after in replacements.items():assert source.count(before)==1,before;source=source.replace(before,after)
    (out/'reused-column-completion.py').write_text(source)
    namespace=dict(vars(prior));exec(compile(source,__file__,'exec'),namespace);namespace['main']()


if __name__=='__main__':main()
