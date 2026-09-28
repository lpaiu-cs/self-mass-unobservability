"""Execute the frozen radial correction with the copied function's defaults.

The first launcher failed before computing any path: FunctionType does not
copy __defaults__. Preserve that failure and bind this execution-only repair.
"""
from pathlib import Path
import json
import def_scalar_regular_radius as model


def run():
    parent=model.OUT
    plan=model.bindings()
    failure=json.loads((parent/'failure.json').read_text())
    assert 'duration' in failure['error']
    model.evolve.__defaults__=model.prior.evolve.__defaults__
    model.OUT=parent/'execution-fixed'
    assert not model.OUT.exists();model.OUT.mkdir()
    plan['bindings'].update({p.relative_to(model.ROOT).as_posix():model.digest(p)
        for p in [Path(__file__),parent/'plan.json',parent/'failure.json']})
    plan['execution_repair']='Restore the existing duration default on the copied function. No scientific equations, initial data, grids, duration, thresholds or budget change.'
    model.save('plan.json',plan);model.save('pilot.json',dict(passed=True,reused_timing=2.82))
    model.run.__globals__['OUT']=model.OUT
    model.run()


if __name__=='__main__':run()
