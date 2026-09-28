"""Use the saved metric to distinguish the fluid buffer from scalar reach."""
import json
import signal
import time
import numpy as np
from scipy.optimize import brentq
import def_native_reactive_interface as task
from def_native_release_charge import Geometry


def main():
    out=task.PHYSICAL;assert not (out/'causal-scope.json').exists();start=time.monotonic();signal.alarm(15)
    f=task.PhysicalFlow(448);geom=Geometry(f.base);end=float(f.base.hist_t[-1])
    def travel(x):return -float(geom(np.array([x]))[1][0])
    depth=-brentq(lambda x:travel(x)-end,-2*task.C*end,0.,xtol=1e-5,rtol=1e-12)
    roundtrip=abs(travel(-depth)/end-1);assert roundtrip<1e-10
    result=dict(classification='Counterexample candidate',passed=True,coordinate_interval_seconds=end,
        scalar_light_cone_depth_cm=depth,scalar_light_cone_depth_km=depth/1e5,
        hydro_buffer_width_cm=20000.,old_interface_depth_cm=20000.,actual_inner_domain_depth_cm=40000.,
        scalar_travel_to_inner_edge_seconds=travel(-40000.),root_relative=roundtrip,seconds=time.monotonic()-start,
        inference='The matched reactive200m buffer isolates the hydrodynamic old interface over this interval. It does not enclose the scalar past light cone. Inner thermochemical/radiative source changes cannot be replaced by a mass-port-only acoustic profile without further evidence.',
        next_decisive_requirement='Use the actual reactive patch stress in the scalar source, and determine or bound the deeper source through the same physical radiation/background. Do not expand a fluid buffer repeatedly as a substitute for that source calculation.',
        full_stellar_interior=False,full_goal_complete=False,
        bindings={str(p):task.ex.old.cold.sha(p) for p in [task.Path(__file__),task.Path(task.__file__),out/'result.json']})
    task.write(out/'causal-scope.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
