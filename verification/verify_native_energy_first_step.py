"""Native primitive recovery of the first conservative release step."""
import json
import time
import runpy
import numpy as np
import def_native_energy_flow as flow


def main():
    out=flow.OUT;assert not (out/'first-step-native.json').exists();start=time.monotonic()
    flow.write(out/'first-step-plan.json',dict(classification='Counterexample candidate',
        prior_full_entropy_bank_failed=True,claim='Measure the actual first conservative energy state with direct native EOS recovery before choosing any table repair. Use only the previously controlled initial entropy range for the first flux.',
        native_calls=50,seconds=10,new_evolution_steps=1,continue_after_step=False))
    m=runpy.run_path(str(out/'flow-first-step-source.py'))['Flow'](224)
    old=flow.bank.prior.EOS();old.sunit=float(old.d['sunit']);m.eos=old
    rho,v,sigma=m.base.primitive(m.base.initial);m.initial=m.conserved(rho,v,sigma,m.base.a)[0]
    k,ledger,dt=m.rhs(m.initial,0.);trial=m.initial+dt*k;indices=np.where(trial[0]>=old.floor)[0][-5:]
    fan=flow.bank.task.prior.Fan(call_cap=50,reuse=True);rows=[]
    for i in indices:
        D,S,K=trial[:,i];tau=(K-(m.base.a[i]-m.base.a0)*old.cx*D)/m.base.a[i];p=0.;lt=np.log(fan.T)
        for iteration in range(16):
            vv=S/(old.cx*D+tau+p);root=np.sqrt(1-vv*vv);rr=D*root;W=1/root;wm=vv*vv/(root*(1+root))
            a=fan.call(np.log(rr*old.rho0),lt);p=a[1]/(old.rho0*flow.C**2);u=a[2]/flow.C**2
            nr=old.cx*D*wm+(rr*u+p)*W*W-p;error=nr-tau
            if abs(error)/abs(tau)<2e-11:break
            lt-=np.clip(error/(rr*W*W*a[10]/flow.C**2),-.3,.3)
        else:raise AssertionError('Native energy root')
        ss=(a[3]-float(old.d['s0']))/old.sunit
        rows.append(dict(cell=int(i),x_cm=float(m.base.x[i]),density_ratio=float(rr),temperature_K=float(np.exp(lt)),sigma=float(ss),energy_relative=float(abs(error)/abs(tau)),native_iterations=iteration+1))
    np.savez_compressed(out/'first-step-native.npz',trial=trial,initial=m.initial,dt=dt,rate=k,ledger_rate=ledger)
    result=dict(classification='Counterexample candidate',passed=True,dt_seconds=dt,rows=rows,native_calls=fan.calls,seconds=time.monotonic()-start,
        actual_required_maximum_sigma=max(r['sigma'] for r in rows),full_horizon_completed=False)
    flow.write(out/'first-step-native.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
