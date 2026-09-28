"""Apply updated actual photon transfers to free material and compact GR."""
from pathlib import Path
from types import FunctionType
import json,resource,sys,time
import numpy as np
import verify_native_updated_gr_return as prior
import def_native_material_branch_response as branch
import def_native_feedback_return as previous
import def_native_parallel_return as parallel

OUT=prior.OUT/'material-return';GR=OUT/'gr';write=prior.write;sha=prior.sha;AMP=branch.AMP
PATHS=[(64,128),(128,128),(128,64)]


def configure():
    prior.configure();branch.base.photons.OUT=prior.run.RESPONSE


def photon_path(steps,reference):return prior.run.RESPONSE/f'steps-{steps}-reference-{reference}.npz'


class Material(branch.Material):
    def __init__(self,reference,steps=128):
        super().__init__(reference);self.steps=steps
        p=np.load(photon_path(steps,reference));ids=[int(np.argmin(abs(p['t']-t))) for t in self.t]
        assert np.max(abs(p['t'][ids]-self.t))<1e-18
        collision=p['collision_transfer'][ids]
        self.transfer=np.stack([np.zeros_like(collision[:,:,0]),p['moments'][ids,3]/self.a,collision[:,:,0],collision[:,:,1]],axis=1)/AMP

    def point(self,k):
        cached=k in self.cache;row=super().point(k)
        # The moving deep gas has integrated radial kinetic stress2*T.
        # It was previously only bounded as omitted metric work.
        if not cached:row['Pr'][:self.nb]+=2*self.model.kinetic()
        return row

    def rhs(self,t,z,probe=1.):
        return super().rhs(t,z,probe*(8 if self.reference==64 else 1))

    run=FunctionType(branch.base.Material.run.__code__,dict(vars(branch.base),OUT=OUT),argdefs=branch.base.Material.run.__defaults__)


dispatch=FunctionType(parallel.dispatch.__code__,dict(vars(parallel),OUT=OUT,__file__=__file__))


def worker(steps,reference,label,limit,restart):
    cap=2*1024**3;resource.setrlimit(resource.RLIMIT_AS,(cap,cap));configure()
    started=time.monotonic();cpu=time.process_time();m=Material(reference,steps)
    row=m.run(steps,label,limit,restart);row.update(physical_branch_ratio=m.physical_branch_ratio,
        maximum_owner_error=max(p['owner_error'] for p in m.cache.values()),
        worker_wall_seconds=time.monotonic()-started,worker_CPU_seconds=time.process_time()-cpu,
        peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
    if row['completed']==steps and reference==64:
        d=np.load(OUT/f'{label}.npz');z=d['delta_scaled'];t=float(d['time']);a=m.rhs(t,z);b=m.rhs(t,z,2.)
        row['reverse_probe_8_16_relative']=(np.sum(abs(a[0]-b[0]),axis=1)/np.maximum(np.sum(abs(b[0]),axis=1),1.)).tolist()
        row['passed']=row['passed'] and max(row['reverse_probe_8_16_relative'])<.002
    row['passed']=row['passed'] and row['maximum_owner_error']<1e-8 and row['physical_branch_ratio']<.01
    write(OUT/f'{label}.json',row);print(json.dumps(row),flush=True);assert row['passed'],row


def prepare():
    assert not OUT.exists();OUT.mkdir();GR.mkdir()
    assert json.loads((prior.OUT/'audit.json').read_text())['passed']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',before_checkpoint='2469970b1',
        claim='Apply updated actual photon collision energy/H/momentum to free shared material mass/momentum/energy/H, then apply its pressure/trace and photons to compact GR and quantify the charge correction.',
        reuse='Completed Phase132 full531-cell spectra, thermal/H histories, corrected EOS, GR/lapse and banks. Original64/128 clocks and3.434ms; no native calls, physical background replay or extra paths.',
        method='Reuse branch-aware actual material fluxes, SSP/CFL and direct conservative primitive map. Keep1x probes on reference128 and the previously established8x arithmetic repair on reference64, with4/8 and8/16 checks. Restore the known deep radial kinetic stress2*T to metric work.',
        limits='This returns actual transfers to material and compact GR. New material motion must still return to photons; a small finite correction is not a contraction bound, complete exterior treatment or final charge.',
        budgets=dict(pilot_wall_seconds=45,production_wall_seconds=350,CPU_processes=3,threads_each=1,total_virtual_memory_GiB=6,source_seconds=75,GR_seconds=120),
        forecast='Use the previous completed material paths raw-call counts including late CFL subcycles, times each concurrent pilot seconds/raw-call plus15s; require2x maximum below350s. Measure real contention. Later changed-background subcycles remain uncertain.',
        gates=dict(owner=1e-8,conservation=1e-8,directional=.002,branch=.01,small_state=1e-6,time=.02,background=.02,pressure=.002,GR_quadrature=.002,independent_GR=1e-9),
        stop='Stop all workers on failure or350s. No automatic extra path, clock, mesh, horizon, native support or relaxed gate. Preserve all outputs and original failures.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(prior.__file__),Path(branch.__file__),Path(branch.base.__file__),Path(previous.__file__),prior.OUT/'audit.json',previous.OUT/'material-audit.json']}))


def pilot():
    assert not (OUT/'pilot.json').exists()
    rows,seconds=dispatch([(n,r,f'pilot-{n}-{r}',2,None) for n,r in PATHS],45)
    old=json.loads((previous.OUT/'material-audit.json').read_text())['paths'];forecasts=[]
    for row in rows:
        ref=next(r for r in old if (r['steps'],r['reference'])==(row['steps'],row['reference']))
        forecasts.append(ref['raw_owner_calls']*row['seconds']/row['raw_owner_calls']+15)
    p=dict(classification='Counterexample candidate',rows=rows,forecast_each_seconds=forecasts,upper_seconds=2*max(forecasts),
        eligible=all(r['passed'] for r in rows) and 2*max(forecasts)<350,seconds=seconds)
    write(OUT/'pilot.json',p);print(json.dumps(p),flush=True)
    if p['eligible']:write(OUT/'execution-plan.json',dict(classification='Counterexample candidate',eligible=True,hard_seconds=350,
        bindings={str(f):sha(f) for f in [Path(__file__),Path(prior.__file__),Path(branch.__file__),Path(branch.base.__file__),OUT/'plan.json',OUT/'pilot.json']}))


def production():
    assert not (OUT/'production.json').exists();p=json.loads((OUT/'execution-plan.json').read_text());assert p['eligible']
    for f,h in p['bindings'].items():assert sha(f)==h,f
    rows,seconds=dispatch([(n,r,f'steps-{n}-reference-{r}',None,f'pilot-{n}-{r}') for n,r in PATHS],350)
    write(OUT/'production.json',dict(classification='Counterexample candidate',passed=all(r['passed'] for r in rows),paths=rows,seconds=seconds,
        CPU_seconds=sum(r['worker_CPU_seconds'] for r in rows),sum_peak_RSS_bytes=sum(r['peak_RSS_bytes'] for r in rows),additional_material_motion_evolved=True,final_charge_solved=False))
    print((OUT/'production.json').read_text(),flush=True)


if __name__=='__main__':
    if sys.argv[1]=='worker':worker(int(sys.argv[2]),int(sys.argv[3]),sys.argv[4],None if sys.argv[5]=='None' else int(sys.argv[5]),None if sys.argv[6]=='None' else sys.argv[6])
    else:globals()[sys.argv[1]]()
