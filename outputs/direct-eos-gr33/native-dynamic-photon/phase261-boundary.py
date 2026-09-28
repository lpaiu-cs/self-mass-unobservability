"""Project the accepted physical photon stress onto the same vacuum GR constraint.

This exports the missing photon particular source. It is not the reciprocal
scalar source, a solved total boundary, or a post-hoc charge correction.
"""
from pathlib import Path
import json,resource,time
import numpy as np
import sympy as sp
import propagate_exterior_vacuum as adapter

p=adapter.base;root=p.OUT;out=root/'photon-boundary';start=time.monotonic();error=None
resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,)*2);p.incident.native.deadline(300)
try:
    result=p.read(root/'result.json');assert result['passed'];assert not out.exists();out.mkdir()
    metric_path=p.ACTUAL/'metric/metric-128-g8.npz'
    files=[Path(__file__),Path(p.__file__),Path(adapter.__file__),root/'result.json',metric_path]
    files += [q for k in p.SETTINGS for q in sorted((root/k).glob('snapshot-*.npz'))]
    p.write(out/'plan.json',dict(classification='Conjectural',
        claim='Apply the new actual photon stress distributions to the mass and lapse constraint particular source at the existing outer face, preserving the already-computed matter mass and scalar boundary owners.',
        method='In the same scalar vacuum K(r)=1/(r*N*sqrt(b)) is the exact unit-mass lapse tail. Differentiate the packet mass/pressure measure including its moving radius and angle; do not vary only its energy.',
        decision='Export the missing signed photon mass/lapse source for the next complete mixed-constraint and feedback solve. Do not add this to a charge or relabel the total physical boundary as solved.',
        gates=dict(quadrature=.002,mass_kernel=1e-10,symbolic=True),budget_seconds=300,
        reused_same_history=True,new_physical_steps=0,
        boundaries_not_yet_closed=['Reciprocal scalar stress and variation of the exterior background operator.','Returned-metric geometric photon response.','Mass normalization and actual application of the complete new boundary to the coupled solution.'],
        bindings={str(q):p.sha(q) for q in files}))
    r,a,b,m,dr,mu,dmu,dH,nu,lam,s=sp.symbols('r a b m dr mu dmu dH nu lam s',nonzero=True)
    aprime=2*m*a/(r*r*b)
    assert sp.simplify((-1/(r*r*a)-aprime/(r*a*a)+1/(r*r*a*b)).subs(b,1-2*m/r))==0
    local=dH-nu-lam;K=1/(r*a);Kprime=-K/(r*b)
    expr=sp.exp(s*local)*(1+(mu+s*dmu)**2)*(K+s*dr*Kprime)
    expected=K*((1+mu*mu)*(local-dr/(r*b))+2*mu*dmu)
    assert sp.simplify(sp.diff(expr,s).subs(s,0)-expected)==0
    p.write(out/'symbolic.json',dict(classification='Proven',passed=True,
        vacuum='K(r)=int_r^infinity ds/(N*b^(3/2)*s^2)=1/(r*N*sqrt(b)). The scalar-vacuum equations give a_prime/a=2*m/(r^2*b).',
        packet='delta_nu_photon_particular(r0)=-G/c^4*sum E0*K(r)*[(1+mu^2)*(delta_lnH-delta_nu-delta_lambda-delta_r/(r*b))+2*mu*delta_mu].',
        mass='delta_J_photon(infinity)=G/c^4*sum E0*(delta_lnH-delta_nu-delta_lambda). A photon source is not a complete ADM balance.',
        scope='Particular response of the original vacuum linear constraint to the new physical photon stress only. Scalar reciprocal/cross-operator terms and the boundary initial value are not set to zero.'))
    p.initialize();model=adapter.Metric(8);d=model.d;paths={};rows=[]
    for label in p.SETTINGS:
        values=[]
        for j in range(1,17):
            z=np.load(root/label/f'snapshot-{j}.npz');r=z['radius_cm'];mu=z['direction'];dr=z['delta_radius_cm'];dm=z['delta_direction']
            _,N,b,a,_=d.bg.metric(r/d.model.m.R);K=1/(r*a)
            local=z['delta_log_H']-z['metric_nu']-z['metric_lambda'];energy=z['background_packet_energy_erg']
            mass=p.LD(p.G)/p.LD(p.C)**4*np.sum(energy*local,dtype=p.LD)
            pieces=np.array([-p.LD(p.G)/p.LD(p.C)**4*np.sum(energy*K*v,dtype=p.LD) for v in
                [(1+mu*mu)*local,-(1+mu*mu)*dr/(r*b),2*mu*dm]])
            values.append(np.r_[mass,pieces.sum(dtype=p.LD),pieces])
        paths[label]=np.asarray(values,p.LD)
        rows.append(dict(path=label,terminal_photon_J_source_cm=float(paths[label][-1,0]),
            terminal_photon_lapse_source=float(paths[label][-1,1]),terminal_lapse_energy_radius_angle_parts=np.asarray(paths[label][-1,2:],float).tolist()))
    norm=np.maximum(np.max(abs(paths['fine'][:,:2]),axis=0),p.LD('1e-290'))
    controls={k:np.asarray(np.max(abs(v[:,:2]-paths['fine'][:,:2]),axis=0)/norm,float).tolist() for k,v in paths.items() if k!='fine'}
    met=np.load(metric_path);t=np.linspace(0,d.T,17);ids=np.array([np.argmin(abs(met['t']-v)) for v in t])
    assert np.max(abs(met['t'][ids]-t))<1e-18
    actual=np.asarray(met['delta_nu_faces'][ids,-1],p.LD)
    residual=met['asymptotic_mass_residual_cm'][ids];mask=abs(residual)>0
    inferred=-met['outer_ADM_residual_lapse'][ids][mask]/residual[mask]
    K0=1/(d.r0*d.z0['lapse'][0]*np.sqrt(d.z0['b'][0]))
    kernel=float(np.max(abs(inferred/K0-1)));assert kernel<1e-10,kernel
    passed=all(max(v)<.002 for v in controls.values())
    np.savez_compressed(out/'source.npz',t=t,photon_J_source_cm=np.r_[p.LD(0),paths['fine'][:,0]],
        photon_lapse_source=np.r_[p.LD(0),paths['fine'][:,1]],lapse_energy_radius_angle_parts=paths['fine'][:,2:],
        existing_applied_total_lapse=actual)
    verdict=dict(classification='Counterexample candidate',passed=passed,controls=controls,rows=rows,
        independent_saved_mass_kernel_relative=kernel,
        maximum_particular_lapse_over_existing_applied_lapse=float(np.max(abs(paths['fine'][:,1]))/np.max(abs(actual))),
        actual_photon_stress_inserted_in_GR_constraint=True,complete_physical_boundary=False,
        corrected_boundary_applied_to_matter=False,physical_final_charge_solved=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    p.write(out/'result.json',verdict);print(json.dumps(verdict),flush=True);assert passed,controls
except BaseException as exc:error=repr(exc);raise
finally:
    if out.exists():p.write(out/'receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=p.sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
