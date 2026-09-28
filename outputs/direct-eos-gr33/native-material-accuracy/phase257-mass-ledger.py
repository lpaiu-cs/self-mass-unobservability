"""Locate the failed mass clock comparison in the same completed solution."""
from pathlib import Path
from decimal import Decimal,localcontext
import json,resource,time
import numpy as np
import read_returned_exterior as p

OUT=Path('native-material-exterior257-work');SOURCE=Path('native-material-charge257-work/full/charge/gr')
assert not p.read(OUT/'full/audit.json')['passed']
resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,)*2)
p.common.prior.prior.joint.previous.original.inf.incident.native.deadline(180)
start=time.monotonic()
p.write(OUT/'mass-ledger-plan.json',dict(classification='Conjectural',
    decision='Distinguish endpoint source/port cancellation and floating arithmetic from the original failed2percent clock comparison before choosing any additional physical solve.',
    budget_seconds=180,virtual_GiB=8,physical_steps=0,new_clock_paths=0,
    original_mass_time_gate=.02,original_verdict_preserved=True,
    scope='No replacement of the original audit, no continuum extrapolation and no post-result sign certificate.'))
p.bind(p.charge.base.endpoint.initialize,OUT=OUT/'full/initialization')()
m=p.exterior.Exterior();con=p.common.prior.geometry.constraints;rows=[]
terms={};history={}
for n in [64,128]:
    d=dict(np.load(SOURCE/f'source-{n}.npz'));q=con.flow.initial.Quadrature(d['edges'],8)
    _,_,a,B,_=m.response.geo(q.r.ravel()-m.response.model.m.RJ)
    measure=q.r*q.r*B.reshape(q.r.shape);mean=q.h*((measure*a.reshape(q.r.shape))@q.w)/(q.h*(measure@q.w))
    rest=np.asarray(d['baryon_g'],p.LD)*p.LD(d['cx'])*p.LD(p.C)**2
    labels=['rest','gas_nonrest','photons','inner_port','outer_port']
    values=[np.sum(v*mean,axis=1,dtype=p.LD) for v in [rest,d['gas_nonrest_energy_erg'],d['photon_energy_erg']]]
    values += [-d['inner_cumulative_energy_erg'],d['outer_cumulative_energy_erg']]
    v=np.array(values,p.LD).T;total=v.sum(1,dtype=p.LD);terms[n]=v[-1];history[n]=total
    end=next(r for r in p.read(OUT/'full/result.json')['rows'] if (r['component'],r['clock'],r['angular'],r['radial'])==('low',n,8,8))
    mass=p.G/p.C**4*total[-1]
    reproduction=float(abs(mass-p.LD(end['homogeneous_mass_cm']))/abs(mass));assert reproduction<1e-12,reproduction
    def D(x):
        a,b=p.LD(x).as_integer_ratio();return Decimal(a)/Decimal(b)
    with localcontext() as ctx:
        ctx.prec=80
        exact=sum((D(w)*D(v) for source in [rest,d['gas_nonrest_energy_erg'],d['photon_energy_erg']] for w,v in zip(mean,source[-1])),Decimal(0))
        exact-=D(d['inner_cumulative_energy_erg'][-1]);exact+=D(d['outer_cumulative_energy_erg'][-1])
        rounding=float(abs(exact-D(total[-1]))/abs(exact))
    rows.append(dict(clock=n,energy_terms_erg={k:float(x) for k,x in zip(labels,v[-1])},total_erg=float(total[-1]),
        homogeneous_cm=float(mass),stored_mass_reproduction_relative=reproduction,
        exact80_input_sum_erg=str(exact),summation_relative=rounding,
        cancellation_condition=float(np.sum(abs(v[-1]))/abs(total[-1]))))
    np.savez_compressed(OUT/f'mass-ledger-{n}.npz',t=d['t'],terms=v,total_erg=total)
diff=terms[64]-terms[128]
result=dict(classification='Counterexample candidate',rows=rows,
    clock_difference_terms_erg={k:float(x) for k,x in zip(labels,diff)},
    clock_difference_erg=float(sum(diff)),relative_clock_difference=float(abs(sum(diff))/abs(history[128][-1])),
    original_gate=.02,original_audit_passed=False,original_verdict_preserved=True,
    physical_final_charge_solved=False,final_charge_conclusion='unadjudicated',full_goal_complete=False,
    seconds=time.monotonic()-start,peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    bindings={str(x):p.sha(x) for x in [Path(__file__),OUT/'mass-ledger-plan.json',OUT/'full/result.json',OUT/'full/audit.json',SOURCE/'source-64.npz',SOURCE/'source-128.npz']})
p.write(OUT/'mass-ledger.json',result);print(json.dumps({k:v for k,v in result.items() if k!='bindings'}))
