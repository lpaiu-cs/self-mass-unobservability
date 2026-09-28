"""Input-defined scalar cascade at the actual primary front; no evolution replay."""
import time
start=time.monotonic()
import numpy as np
import compare_full_incident_fluid_time as run

run.prior.joint.previous.original.inf.incident.native.deadline(8)
spent=run.read(run.OUT/'audit-receipt.json')['seconds']+run.read(run.OUT/'time-localization.json')['seconds']
assert spent+8<15
d=run.prior.drive.Driver(8);T=d.T;end=T/32;h0=T/64
a,b,c=run.prior.joint.A,run.prior.joint.B,run.prior.joint.C
clock=np.arange(3)*h0;arrivals=-d.xc/run.prior.C
old=np.load(run.OUT/'sweep-1/photons/pilot-64.npz')
rows=[]
co=np.r_[np.zeros(4),256*np.array([1,-4,6,-4,1])];integ=np.polynomial.polynomial.polyint(co)
for cell in [15,16,17,18,143,144,145,146,147,148,149]:
    arrival=arrivals[cell];values=[]
    for n in [2,4,8]:
        h=end/n;e=q=0.
        for j in range(n):
            now=(j+c)*h;s=(now-arrival)/d.D
            forcing=run.prior.drive.pulse(s,1)/d.D
            es=e+h*(a@forcing);qs=q+h*(a@es);e=es[-1];q=qs[-1]
        values.append([e,q])
    s=max(0,min(1,(end-arrival)/d.D));exact=np.array([run.prior.drive.pulse(s),d.D*np.polynomial.polynomial.polyval(s,integ)])
    v=np.array(values);relative=abs(v-exact)/np.maximum(abs(exact),1e-290)
    pairs=[(abs(v[i]-v[i+1])/np.maximum(abs(v[i+1]),1e-290)).tolist() for i in [0,1]]
    rows.append(dict(cell=cell,arrival_seconds=float(arrival),arrival_inside_prefix=bool(0<arrival<end),
        pair_relative=pairs,exact_relative=relative.tolist(),exact=exact.tolist()))
new_flags={}
for n in [64,128]:
    h=T/n;starts=np.r_[arrivals,arrivals+d.D]
    new_flags[str(n)]=[bool(any(t<(k+1)*h and t+2*h0>k*h and 0<=t<T for t in starts)) for k in range(n//32)]
result=dict(classification='Conjectural',no_trajectory_replay=True,rows=rows,
    old_split_macros=old['split_macro_steps'][:2].tolist(),candidate_material_front_split_macros=new_flags,
    toy='Eprime=primary_pulse_prime, Bprime=E with zero state and measured unmodified input arrival. Exact E=primary_pulse, B=its polynomial primitive. No fit to solved amplitudes.',
    interpretation='A material response cascade can be underresolved at a front excluded by the photon-slow-cell split rule. Scalar agreement motivates a material-front rule, not a full-system error theorem or a passing physical time verdict.',
    seconds=time.monotonic()-start,source_sha256=run.sha(__file__))
run.write(run.OUT/'front-cascade.json',result);print(result)
