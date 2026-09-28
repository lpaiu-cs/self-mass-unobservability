"""Independent original time-jet trajectories for two preregistered packets."""
from pathlib import Path
import json,resource,time
import numpy as np
import propagate_primitive_characteristics as current
b=current.b;p=current.p;root=current.OUT;out=root/'direct-reference'
read,write,sha=current.read,current.write,current.sha
start=time.monotonic();error=None
resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,)*2);p.incident.native.deadline(3600)
try:
    assert not out.exists();out.mkdir()
    write(out/'plan.json',dict(classification='Conjectural',packets=[0,31],emission_cell=0,
        claim='Compare unchanged direct time-jet propagation with the primitive method on the same earliest cohort, at near-grazing and near-radial physical quadrature directions.',
        gate=.002,rtol=2e-8,atol=2e-11,maximum_seconds=3600,
        scope='Independent actual-trajectory control for the exact variable transformation. Two representative rays are not a uniform continuum certificate.',
        bindings={str(q):sha(q) for q in [Path(__file__),Path(current.__file__),Path(b.__file__),Path(p.__file__),root/'plan.json',root/'check.json']}))
    p.OUT=out; p.Metric=b.Metric; b.Moments=current.base.previous.Moments
    p.initialize();m=p.Photons(8,8);original=m.cohorts;rows=[]
    for index in [0,31]:
        def cohorts(now,order,cells=None):
            values=original(now,order,[0]);return tuple(v[index:index+1] for v in values)
        m.cohorts=cohorts
        z,row=m.propagate(m.d.T,8,[0]);row['packet_index']=index;row['work_identity_independent']=True
        np.savez_compressed(out/f'packet-{index}.npz',**z);write(out/f'packet-{index}.json',row);rows.append(row)
        write(out/'progress.json',dict(rows=rows))
    fine=np.load(root/'pilot-0.npz');errors={}
    for index in [0,31]:
        direct=np.load(out/f'packet-{index}.npz')
        for key in ['emission_t','owner','background_packet_energy_erg','radius_cm','direction']:
            assert np.array_equal(direct[key],fine[key][index:index+1]),(index,key)
        for key in ['delta_radius_cm','delta_direction','delta_log_H','delta_arrival_seconds','integrated_log_H_work']:
            a=direct[key];v=fine[key][index:index+1]
            errors[f'{index}:{key}']=float(np.max(abs(a-v))/max(np.max(abs(a)),1e-290))
    passed=max(errors.values())<.002
    write(out/'result.json',dict(classification='Counterexample candidate',passed=passed,relative=errors,rows=rows,
        original_energy_integral_checked=True,physical_final_charge_solved=False));assert passed,errors
except BaseException as exc:error=repr(exc);raise
finally:
    if out.exists():write(out/'receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
