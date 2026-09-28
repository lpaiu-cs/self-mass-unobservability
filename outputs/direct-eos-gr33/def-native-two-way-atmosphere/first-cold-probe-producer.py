import sys,json,signal,time
import numpy as np
sys.path.insert(0,'verification')
import def_native_two_way_atmosphere as task
p=task.OUT;f=task.Flow(896);m=f.base;d=np.load(p/'root-failure-state.npz');U=d['U'];i=416;D=U[0,i];tau=(U[2,i]-(m.a[i]-m.a0)*f.eos.cx*D)/m.a[i];y=U[3,i]/D
v=U[1,i]/(f.eos.cx*D+tau);root=np.sqrt(1-v*v);rho=D*root;W=1/root;wm=v*v/(root*(1+root));f.eos.y=np.array([y]);tab=[]
for T in [240,310,700]:
    pp,u,*_=f.eos(np.array([rho]),np.array([np.log(T)]));res=(f.eos.cx*D*wm+(rho*u[0]+pp[0])*W*W-pp[0]-tau)/tau;tab.append([T,float(res)])
print(json.dumps(dict(cell=i,rho_relative=rho,y=y,table_residuals=tab,linear_T_root=240-tab[0][1]*(310-240)/(tab[1][1]-tab[0][1]))),flush=True)
plan=dict(classification='Counterexample candidate',claim='Determine whether the actual failing dilute cell has a positive-temperature native EOS root below the retained240K table.',decision='Native positive cold root supports a minimal constitutive extension; absence of such a root requires fixing conservative time/spatial admissibility. No production or support extension authorized by this probe alone.',temperatures_K=[240,160,80],native_call_cap=16,wall_seconds=20,reuse='The exact failing conserved cell; no fluid rerun.',physical_support_not_changed=True)
task.write(p/'cold-root-probe-plan.json',plan);signal.signal(signal.SIGALRM,task.optical.timeout);signal.alarm(20);start=time.monotonic();native=task.optical.ex.Native(cap=16);rows=[]
for T in plan['temperatures_K']:
    state=native.state(float(np.log(rho)),float(np.log(T)),float(y));pp=state['raw'][1]/(f.eos.rho0*task.C**2);u=state['raw'][2]/task.C**2;res=(f.eos.cx*D*wm+(rho*u+pp)*W*W-pp-tau)/tau
    rows.append(dict(T=T,native_relative_energy_residual=float(res),native_population_error=state['population_error']));print(json.dumps(rows[-1]),flush=True)
task.write(p/'cold-root-probe.json',dict(classification='Counterexample candidate',table_residuals=tab,native=rows,native_calls=native.ion.calls,seconds=time.monotonic()-start,finite_native_positive_root_bracketed=rows[0]['native_relative_energy_residual']*rows[-1]['native_relative_energy_residual']<0));signal.alarm(0)
