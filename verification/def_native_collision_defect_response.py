"""Actual native collision defect -> simultaneous photon/thermal/H response."""
from pathlib import Path
from types import FunctionType
import inspect,json,signal,sys,time
import numpy as np
import def_native_history_collision_forcing as forcing

updated=forcing.history.updated;old=updated.old
OUT=forcing.OUT/'response';write=forcing.write;sha=forcing.sha;AMP=old.AMPLITUDE


def configure():updated.configure();old.OUT=OUT


class Coefficients(old.Response):
    def restore(self,k):
        # velocity() reads face momentum j, not the cell Pi field.
        self.model.j=self.data['snapshot_j'][self.keep[k]].copy()
        super().restore(k)
        assert np.array_equal(self.bulk_beta,self.model.velocity())


class Response(updated.Response):
    def __init__(self,reference=128):
        assert reference==128;super().__init__(reference)
        for key in self.g:
            if key.startswith('delta_'):self.g[key]=np.zeros_like(self.g[key])
        self.defects={};self.zero=np.zeros_like(self.I[0])

    def defect(self,k):
        if k not in self.defects:
            self.defects[k]=dict(np.load(forcing.OUT/f'point-{k}.npz'))
            for previous in list(self.defects):
                if previous not in [k-1,k,k+1]:del self.defects[previous]
        return self.defects[k]

    def local(self,t):
        c=super().local(t);j=max(0,min(np.searchsorted(self.t,t,side='right')-1,15));f=(t-self.t[j])/(self.t[j+1]-self.t[j])
        a,b=self.defect(j),self.defect(j+1)
        for key,source in [('q','photon'),('qb','bound'),('qe','escape')]:
            c[key]=((1-f)*a[source]+f*b[source])/AMP
            if key!='qe':c[key]/=self.scale
        return c

    def source(self,t):return self.zero,np.zeros(5),0.
    run=FunctionType(updated.Response.run.__code__,dict(updated.Response.run.__globals__,OUT=OUT),argdefs=updated.Response.run.__defaults__)


def prepare():
    assert json.loads((forcing.OUT/'production.json').read_text())['passed'];assert not OUT.exists();OUT.mkdir()
    paths=[Path(__file__),Path(forcing.__file__),Path(updated.__file__),Path(old.__file__),
           forcing.OUT/'production.json',forcing.history.OUT/'production-samples.json']
    paths+=list(forcing.OUT.glob('point-*.npz'))
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        claim='Evolve the actual native-minus-table collision defect in the same simultaneous moving photon/thermal/H equations, then deliver the paired histories and impulse to free material and GR.',
        operator_repair='Rebuild the17 fine-background coefficient jets with saved face momentum j restored before velocity(). The previous response-bank restore also omitted j. Same EOS interpolation and original sampled derivative gate1e-4; no native table or physical replay.',
        model='Zero additional GR forcing for this isolated constitutive defect. Fixed saved density/motion in this radiation/thermal/H sweep; paired force is retained for subsequent free material. This is a linear defect response, not nonlinear native-EOS re-evolution or a certified error bound.',
        source='Actual photon occupation derivative, bound-free number derivative and frequency-exit moments at17 stored knots. Same existing weights form exactly opposite gas energy/H and momentum. No escaped-photon energy is reassigned as gas heat.',
        paths=[[64,128],[128,128]],reuse='Existing531 cells,8 angles,152 frequencies,3.434ms,stage-time SDIRK and1e-12 full linear residual; resume4/8 pilot prefixes. Native forcing is available only on background128, so no native background64 certificate is claimed.',
        budgets=dict(bank_seconds=90,pilot_seconds=90,production_seconds=900,CPU_threads=1,virtual_GiB=3,new_native_states=0),
        forecast='Actual Phase132 two matching complete radiation paths took about314s. Measure4/8 prefixes and all-point setup; use2x of the larger prefix extrapolation and completed-path late-cost scaling, plus20s. No dispatch above900s.',
        gates=dict(derivative=1e-4,forcing_number=1e-12,linear=1e-12,energy=1e-8,species=1e-8,time=.02),
        stop='No automatic finer clocks, support expansion or budget increase. Preserve failed paths. The coupled objective remains open until free-material/GR return and its required feedback/error controls are evaluated.',
        bindings={str(p):sha(p) for p in paths}))


def bank():
    configure();source=inspect.getsource(old.bank).replace('for reference in [128,64]:','for reference in [128]:')
    scope=dict(vars(old),OUT=OUT,Response=Coefficients);exec(compile(source,__file__,'exec'),scope)
    (OUT/'expanded-bank.py').write_text(source);scope['bank']()


def forcing_check(m):
    rows=[]
    for k in range(17):
        d=m.defect(k);q=d['photon'];b=d['bound'];e=d['escape'];N=m.Nweight/m.scale
        defect=np.sum((q-b)*N,axis=(1,2))+e[0]
        scale=np.maximum(np.sum((abs(q)+abs(b))*N,axis=(1,2))+abs(e[0]),1.)
        relative=float(np.max(abs(defect)/scale));assert relative<1e-12,('Native source number pairing',k,relative)
        row=np.load(OUT/f'bank-128/point-{k}.npz');m.restore(k)
        # Coefficient restore is checked independently from the stored j.
        m.model.j=m.data['snapshot_j'][m.keep[k]].copy();beta=m.model.velocity()
        assert np.array_equal(row['beta'][:m.nb],beta)
        rows.append(dict(k=k,number_relative=relative,deep_velocity_max=float(np.max(abs(beta)))))
    write(OUT/'forcing-check.json',dict(classification='Counterexample candidate',passed=True,rows=rows,
        same_native_forcing_weights=True,actual_saved_j_in_operator=True,uniform_derivative_bound=False))


def pilot():
    assert json.loads((OUT/'bank-result.json').read_text())['passed'];assert not (OUT/'pilot.json').exists();configure()
    start=time.monotonic();signal.alarm(90);m=Response();forcing_check(m);rows=[]
    for n,k in [(64,4),(128,8)]:rows.append((m if n==64 else Response()).run(n,f'pilot-{n}',k))
    previous=json.loads((updated.run.RESPONSE/'result.json').read_text())['paths']
    estimates=[]
    for r in rows:
        before=next(v for v in previous if v['steps']==r['steps'] and v['reference']==128)
        point=r['operator_point_seconds']/max(r['operator_points'],1)
        step=r['stepping_seconds']/r['new_steps']
        # Include the accepted completed path's late mean step cost as a floor.
        step=max(step,before['stepping_seconds']/before['new_steps'])
        estimates.append(17*point+(r['steps']-r['completed_steps'])*step+20)
    a=np.load(OUT/'pilot-64.npz')['moments'][-1,[0,1,2,3,5,6]];b=np.load(OUT/'pilot-128.npz')['moments'][-1,[0,1,2,3,5,6]]
    errors=(np.sum(abs(a-b),axis=1)/np.maximum(np.sum(abs(b),axis=1),1.)).tolist();upper=2*sum(estimates)
    result=dict(classification='Counterexample candidate',rows=rows,equal_horizon=errors,forecast_each_seconds=estimates,
        upper_seconds=upper,eligible=all(r['passed'] for r in rows) and max(errors)<.02 and upper<900,seconds=time.monotonic()-start)
    write(OUT/'pilot.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


def production():
    p=json.loads((OUT/'pilot.json').read_text());assert p['eligible'];assert not (OUT/'result.json').exists();configure()
    start=time.monotonic();signal.alarm(900);rows=[]
    for n in [64,128]:rows.append(Response().run(n,f'steps-{n}-reference-128',restart=f'pilot-{n}'))
    a=np.load(OUT/'steps-64-reference-128.npz')['moments'][:,[0,1,2,3,5,6]]
    b=np.load(OUT/'steps-128-reference-128.npz')['moments'][:,[0,1,2,3,5,6]]
    errors=(np.max(np.sum(abs(a-b),axis=2),axis=0)/np.maximum(np.max(np.sum(abs(b),axis=2),axis=0),1.)).tolist()
    result=dict(classification='Counterexample candidate',passed=all(r['passed'] for r in rows) and max(errors)<.02,
        paths=rows,time_comparison=errors,order=['photon_energy','material_energy','neutral_H','momentum_impulse','photon_radial_pressure','material_pressure'],
        seconds=time.monotonic()-start,actual_native_defect_in_finite_response=True,free_material_return=False,GR_readout=False,
        native_background_time_audit=False,uniform_EOS_derivative_bound=False,coupled_fixed_point_verified=False,final_charge_solved=False)
    write(OUT/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True);assert result['passed']


if __name__=='__main__':
    cap=3*1024**3;forcing.resource.setrlimit(forcing.resource.RLIMIT_AS,(cap,cap));signal.signal(signal.SIGALRM,forcing.history.flow.old.optical.timeout)
    action=sys.argv[1]
    try:globals()[action]()
    except Exception as exc:write(OUT/f'{action}-failure.json',dict(error=repr(exc)));raise
