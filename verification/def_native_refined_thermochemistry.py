"""Apply improved native thermal interpolation to actual coupled evolution.

Counterexample candidate: retain physical support, chemical planes, mesh,
clocks and initial conserved inventories. Preserve the rejected old bank.
"""
from pathlib import Path
from types import FunctionType
import inspect,json,signal,sys,time
import numpy as np
import def_native_conservative_rates as failed
import def_native_stage_energy_history as capture

flow=failed.flow;chem=failed.chem;C=flow.C;write=flow.write;sha=flow.sha
OUT=failed.OUT/'thermal-refined';EV=OUT/'evolution';GR=EV/'gr'


def prepare():
    assert not OUT.exists();OUT.mkdir();EV.mkdir();GR.mkdir()
    old=json.loads((failed.OUT/'result.json').read_text());assert not old['passed']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='095fd9cf6',
        failure='Positive native rates differ by less than9e-8, but net neutral exchange differs6.875percent on the actual conserved trajectory. Large-channel agreement does not control their cancellation.',
        evidence='Eight withheld native evaluations separate chemical-plane interpolation at four dominant cells: at most3.34e-10 relative. Temperature interpolation is the actionable larger approximation.',
        repair='Reuse all190 native nodes; add the four original-interval midpoints at both existing H planes,152states. Rebuild the SAME thermal EOS and Boltzmann-factored spectral cubic. Retain physical T/rho/y support, density/inventory/moving spectral approximations, geometry and all gates.',
        decision='Require original0.2percent positive-rate and2percent net-exchange controls on all323 withheld conserved states, then install the new constitutive bank in actual material/photon64/128 paths and read the resulting GR charge.',
        initial='Preserve original baryon,energy,H,photon and metric data. Recover initial temperature with the new thermal EOS; retain initial pressure mismatch against the old mechanical/scalar-balanced baseline. No equilibrium subtraction.',
        arithmetic='Evaluate collision differences through the linear collision owner on signed coefficient differences. Subtracting two large collision arrays failed the original1e-10 relative-defect check; preserve it and use a difference-first identity with the same threshold.',
        budget=dict(bank_seconds=60,bank_native_calls=1000,controls_seconds=40,pilot_seconds=40,production_seconds=650,charge_seconds=60,CPU_threads=1,memory_GB=3),
        forecast='Native dimension test45calls4.05s including setup;152 new bank states expected10-45s. Full old identical-step capture cost401s. Use new/old pilot cost ratio and1.5x the measured full old cost plus20s; require below650s before dispatch. State-dependent late cost remains unmeasured.',
        gates=dict(positive_rates=.002,net_exchange=.02,owner=1e-10,energy=1e-8,time_charge=.02,quadrature=.002),
        stop='No extra physical cells, clocks, horizon, support or repeated table expansion. Stop on failed native controls, budget or coupled gate; no time-only fitted source.',
        limits='Updated thermal/rate representation of the same restricted H/Thomson model. No complete opacity, exact advected other inventories, atmospheric native inversion, nonlinear GR or final-charge certificate.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(failed.__file__),Path(capture.__file__),flow.previous.OUT/'thermal-support/bank.npz',failed.OUT/'dimension-result.json',failed.OUT/'result.json']}))


def bank():
    assert not (OUT/'bank.json').exists();started=time.monotonic()
    signal.signal(signal.SIGALRM,flow.old.optical.timeout);signal.alarm(60)
    d=np.load(flow.previous.OUT/'geometry.npz');old=dict(np.load(flow.previous.OUT/'thermal-support/bank.npz'))
    offsets=(old['offsets'][:-1]+old['offsets'][1:])/2;n=chem.old.Native(cap=1000)
    raw=np.zeros((19,2,4,21));rates=np.zeros((19,2,4,2,len(d['Einf'])));done=np.zeros((19,2,4),bool)
    try:
        for j in range(19):
            chem.setup(n,d,j)
            for iy,ratio in enumerate(old['ratios']):
                for it,dt in enumerate(offsets):
                    s=n.state(0.,float(np.log(d['T'][j])+dt),float(n.y0*ratio))
                    chi,em=chem.prior.coefficients(n,s,d['Einf']/d['a'][j]);raw[j,iy,it]=s['raw']
                    rates[j,iy,it]=[(chi+em)/s['y'],em/(1-s['y'])];done[j,iy,it]=True
        order=np.argsort(np.r_[old['offsets'],offsets])
        np.savez_compressed(OUT/'bank.npz',raw=np.concatenate([old['raw'],raw],axis=2)[:,:,order],
            rates=np.concatenate([old['rates'],rates],axis=2)[:,:,order],offsets=np.r_[old['offsets'],offsets][order],
            ratios=old['ratios'],native_calls=n.ion.calls)
        write(OUT/'bank.json',dict(classification='Counterexample candidate',passed=True,new_states=int(done.sum()),
            reused_states=190,native_calls=n.ion.calls,seconds=time.monotonic()-started))
    except Exception as exc:
        np.savez_compressed(OUT/'partial-bank.npz',raw=raw,rates=rates,done=done,offsets=offsets)
        write(OUT/'bank-failure.json',dict(error=repr(exc),completed=int(done.sum()),native_calls=n.ion.calls,seconds=time.monotonic()-started));raise
    finally:signal.alarm(0)
    print((OUT/'bank.json').read_text(),flush=True)


# Reuse the actual thermal constructor, with the added bank's zero-T anchor
# located by value instead of the old five-node grid's hardcoded index1.
init=flow.previous.thermal_init.replace("OUT/'thermal-support/bank.npz'","BANK")
init=init.replace('rho=d[',"zero=int(np.argmin(abs(v['offsets'])));assert v['offsets'][zero]==0\n    rho=d[",1)
init=init.replace('self.raw[:,:,1,1]','self.raw[:,:,zero,1]').replace('self.raw[:,:,1,2]','self.raw[:,:,zero,2]').replace('self.rates[:,:,1]','self.rates[:,:,zero]')
ns=dict(flow.previous.tns,BANK=OUT/'bank.npz');exec(compile(init,__file__,'exec'),ns)
class Thermal(flow.previous.Thermal):__init__=ns['__init__']


class Coupled(capture.Capture):
    def __init__(self):
        super().__init__();self.bulk.eos.base=Thermal()
        # b.u0 and f0 remain the original conserved Cauchy reference; gas
        # recovery (including time zero) supplies the new consistent T.


def controls():
    target=OUT/'controls';assert not target.exists();target.mkdir()
    src=inspect.getsource(failed.audit).replace('def audit():','def run():')
    src=src.replace("for p,h in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert sha(p)==h,p",'')
    src=src.replace('m=flow.Coupled();','m=Coupled();')
    src=src.replace("theta=z['snapshot_theta'][k];eta=z['snapshot_eta'][k];", "eta=z['snapshot_eta'][k];theta=e.recover(z['snapshot_u'][k],eta,z['snapshot_theta'][k]);")
    src=src.replace('defect=after[0]-before[0];estimate=raw+after[1]-before[1]+ds', '''subtractive=after[0]-before[0]
        def delta_gas(t,y):
            v=list(gas(t,y));v[6]=ne-oldgas[6];return v
        def delta_rad(t,y):
            v=list(rad(t,y));v[0]=nr[:,0]-oldrad[0];v[1]=nr[:,1]-oldrad[1];return v
        e.gas=delta_gas;e.radiation=delta_rad;scatter_error=m.scatter_number_error
        try:delta=m.collision(I,theta,eta,beta)
        finally:e.gas=gas;e.radiation=rad;m.scatter_number_error=scatter_error
        defect=delta[0];estimate=raw+delta[1]+ds''')
    src=src.replace("(oldrad[:2] if v is before else (nr[:,0],nr[:,1]))", "(oldrad[:2] if v is before else ((nr[:,0]-oldrad[0],nr[:,1]-oldrad[1]) if v is delta else (nr[:,0],nr[:,1])))")
    src=src.replace('original=integrated(before);actual=integrated(after);moments.append([original,actual,actual-original])',
        'original=integrated(before);change=integrated(delta);actual=original+change;moments.append([original,actual,change])')
    src=src.replace('np.trapz(','np.trapezoid(')
    scope=dict(vars(failed),OUT=target,Coupled=Coupled);exec(compile(src,__file__,'exec'),scope)
    (target/'expanded-check.py').write_text(src);scope['run']()


runner=inspect.getsource(capture.Capture.run_capture)
import textwrap
runner=textwrap.dedent(runner)
runner=runner.replace('theta=np.zeros(b.n);eta=np.zeros(b.n);','eta=np.zeros(b.n);theta=b.eos.recover(u,eta,np.zeros(b.n));')
runner=runner.replace('replay_state_relative=state_errors,legacy_source_relative=max(matched,default=0.)',
    'prior_state_change=state_errors,prior_source_change=max(matched,default=0.)')
runner=runner.replace('passed=balance<1e-8 and max(state_errors.values())<1e-8 and max(matched,default=0.)<1e-8',
    'passed=bool(balance<1e-8 and abs(ports[0]-ports[1]-owner_port)/scale<1e-10)')
run_ns=dict(vars(capture),OUT=EV);exec(compile(runner,__file__,'exec'),run_ns);Coupled.run_capture=run_ns['run_capture']


def pilot():
    assert json.loads((OUT/'controls/result.json').read_text())['passed'];assert not (OUT/'pilot.json').exists()
    (EV/'expanded-run.py').write_text(runner);started=time.monotonic();signal.signal(signal.SIGALRM,flow.old.optical.timeout);signal.alarm(40)
    rows=[]
    for n in [64,128]:
        checkpoint=EV/f'checkpoint-{n}.npz';begin=int(np.load(checkpoint)['completed']) if checkpoint.exists() else 0
        row=Coupled().run_capture(n,begin+2,resume=begin>0);row['new_pilot_steps']=2;rows.append(row)
    old=json.loads((capture.OUT/'pilot.json').read_text())['paths'];full=json.loads((capture.OUT/'result.json').read_text())['seconds']
    ratio=sum(r['seconds'] for r in rows)/sum(r['seconds'] for r in old);forecast=full*ratio
    upper=1.5*forecast+20;eligible=all(r['passed'] for r in rows) and upper<650
    write(OUT/'pilot.json',dict(classification='Counterexample candidate',paths=rows,ratio=ratio,forecast_seconds=forecast,upper_seconds=upper,eligible=eligible,seconds=time.monotonic()-started))
    if eligible:write(OUT/'execution-plan.json',dict(classification='Counterexample candidate',eligible=True,cap_seconds=650,forecast_seconds=forecast,upper_seconds=upper,
        bindings={str(p):sha(p) for p in [Path(__file__),Path(capture.__file__),OUT/'bank.npz',OUT/'controls/result.json',OUT/'pilot.json']}))
    signal.alarm(0);print((OUT/'pilot.json').read_text(),flush=True)


def production():
    assert not (EV/'result.json').exists();p=json.loads((OUT/'execution-plan.json').read_text());assert p['eligible']
    for path,h in p['bindings'].items():assert sha(path)==h,path
    start=time.monotonic();signal.signal(signal.SIGALRM,flow.old.optical.timeout);signal.alarm(p['cap_seconds']);rows=[]
    try:
        for steps in [64,128]:
            row=Coupled().run_capture(steps,resume=True);rows.append(row)
            if not row['passed']:break
        write(EV/'result.json',dict(classification='Counterexample candidate',passed=len(rows)==2 and all(r['passed'] for r in rows),paths=rows,seconds=time.monotonic()-start,
            new_native_thermal_bank_in_actual_evolution=True,full_GR_feedback=False,final_charge_solved=False))
    except Exception as exc:write(EV/'failure.json',dict(error=repr(exc),completed=rows,seconds=time.monotonic()-start));raise
    finally:signal.alarm(0)


if __name__=='__main__':globals()[sys.argv[1]]()
