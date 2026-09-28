"""Separate one-off worker setup from measured native per-cell throughput."""
from pathlib import Path
import inspect
import def_reactive_direct_gate as direct

ns=dict(direct.ns,OUT=direct.ns['OUT'].parent/'warm-budget',__file__=__file__,direct=direct)
text=inspect.getsource(direct.paired.old.run)
text=text.replace('files=[Path(__file__),','files=[Path(__file__),Path(direct.__file__),Path(direct.paired.__file__),Path(direct.paired.old.__file__),Path(direct.paired.old.thermal.g.s.__file__),')
old="""        forecast=seconds*5735/len(indices)"""
new="""        warm_start=time.monotonic();warm=[]
        for part in pool.map(block,np.array_split(indices,8)):warm.extend(part)
        warm_seconds=time.monotonic()-warm_start
        # Reuse the accepted warm records in the production result.
        rows=warm
        forecast=seconds+1.5*warm_seconds*(5735-len(indices))/len(indices)"""
assert text.count(old)==1;text=text.replace(old,new)
text=text.replace('seconds=seconds,projection_seconds=forecast','seconds=seconds,warm_seconds=warm_seconds,projection_seconds=forecast')
text=text.replace("claim='Compute", "prior_budget_failure='Cold-pool pilot extrapolation included worker initialization once per small batch; preserve both rejected forecasts. Measure a second 32-cell batch on the same live pool and count setup once. Keep 200 seconds and all science gates unchanged.',\n        claim='Compute")
exec(compile(text,__file__,'exec'),ns)

if __name__=='__main__':ns['run']()
