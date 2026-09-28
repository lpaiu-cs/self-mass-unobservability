"""Audit saved material responses and export nonduplicated GR sources."""
from pathlib import Path
import json,signal,time
import numpy as np
import sympy as sp
import def_native_material_branch_response as run
import def_native_anisotropic_gr as gr

OUT=run.OUT;write=run.write;sha=run.sha;AMP=run.AMP;C=run.base.C;LD=np.longdouble


def pressure(m,k,z,field,eps):
    raw=m.raw(k,z,field,eps);model=m.model;b=model.bulk;f=model.flow;nb=m.nb
    pp=b.eos.gas(raw[3]['theta'],raw[3]['eta'])[0]*b.volume
    rho,v,lt,y=raw[3]['primitive'];f.eos.y=y;p,u,*_=f.eos(rho,lt)
    volume=4*np.pi*model.m.RJ**2*model.m.vol;pg=p*f.eos.rho0*C*C*volume
    pr=pg+volume*(rho*(f.eos.cx+u)+p)*f.eos.rho0*C*C*v*v/(1-v*v)
    return np.array([np.r_[pp,pg],np.r_[pp+2*model.kinetic(),pr]])


def main():
    assert not (OUT/'result.json').exists();start=time.monotonic();production=json.loads((OUT/'production.json').read_text())
    assert production['passed'];write(OUT/'audit-plan.json',dict(classification='Counterexample candidate',
        claim='Compare the saved three material histories, independently recompute their four balances, and recover pressure/anisotropic stress including transported rest mass for the next GR source.',
        controls='Same canonical17 times. Componentwise maximum-in-time spatial L1 differences, divided by the fine maximum-in-time L1. Also check endpoint differences. Native pressure differences use a second independent half-probe.',
        source='Convert integrated reference energy using delta(VE)=(deltaEref+a_surface*cx*c^2*deltaB)/a_ref. For fixed-radius density sources subtract the evolved background density times deltaV. Remove the Phase122 initial canonical gas response using its exact Eg,Pg,Kg coefficients; retain the difference between initial and evolved backgrounds explicitly.',
        limits='One material sweep driven by fixed Phase123 metric and Phase124 photon transfers. New motion is not returned into photons or GR. Sampled directional derivatives are not a uniform certificate.',
        budget_seconds=40,gates=dict(time=.02,background=.02,conservation=1e-8,pressure_probe=.002),
        bindings={str(p):sha(p) for p in [Path(__file__),Path(run.__file__),Path(run.base.__file__),OUT/'production.json']}))
    signal.signal(signal.SIGALRM,run.base.flow.old.optical.timeout);signal.alarm(40)
    g=gr.Response()
    histories=[];source_rows=[];balance=0.;pressure_checks=[];branch_ratio=0.;deep_missing_work=0.
    for steps,ref in [[64,128],[128,128],[128,64]]:
        m=run.Material(ref);coeff=g.coeff(m.rE);d=np.load(OUT/f'steps-{steps}-reference-{ref}.npz');ids=[int(np.argmin(abs(d['t']-t))) for t in m.t];assert np.max(abs(d['t'][ids]-m.t))<1e-18
        hist=d['history_scaled'][ids];histories.append(hist)
        defect=np.sum(d['history_scaled'],axis=2,dtype=LD)+d['discards_scaled']-d['ledgers_scaled']
        balance=max(balance,float(np.max(abs(defect)/np.maximum(d['norms_scaled'],1.))))
        Pg0=coeff['Pg']*C**4/gr.G*m.V;Eg0=coeff['Eg']*C**4/gr.G*m.V;Kg0=coeff['Kg']*C**4/gr.G*m.V
        output=[];probes=[];proper=[]
        for k,t in enumerate(m.t):
            z=hist[k];point=m.point(k);field=m.fields(t)[2];m.min_probe=np.inf;m.rhs(t,z);eps=m.min_probe
            p0=pressure(m,k,np.zeros_like(z),np.zeros_like(field),0.)
            pp=[pressure(m,k,z,field,eps/factor) for factor in [1,2,4]]
            derivative=[(p.astype(LD)-p0.astype(LD))/(eps/factor) for p,factor in zip(pp,[1,2,4])]
            first=np.asarray(2*derivative[1]-derivative[0],float);second=np.asarray(2*derivative[2]-derivative[1],float);probes.append(first-second)
            volume=3*field[0]+field[2];q=point['Q'];backgroundE=(q[2]+m.rest*q[0])/m.a
            deltaE=(z[2]+m.rest*z[0])/m.a
            eforcing=deltaE+(Eg0+Pg0-backgroundE)*volume
            pforcing=second+(Kg0[None]-p0)*volume
            trace=eforcing-pforcing[1]-2*pforcing[0]
            output.append(np.array([eforcing,pforcing[1],trace,pforcing[0]]));proper.append(np.array([deltaE,second[1],second[0]]))
            rates=m.fields(t)[3];m.raw(k,np.zeros_like(z),np.zeros_like(field),0.)
            deep_missing_work=max(deep_missing_work,float(np.sum(abs(2*m.model.kinetic()*(rates[0,:m.nb]+rates[1,:m.nb])*m.a[:m.nb]))*AMP*m.t[-1]))
        output=np.array(output)*AMP;proper=np.array(proper)*AMP
        numerator=np.max(np.sum(abs(np.array(probes)),axis=2),axis=0)
        denominator=np.maximum(np.max(np.sum(abs(proper[:,[2,1]]/AMP),axis=2),axis=0),1.)
        pressure_checks.append((numerator/denominator).tolist());branch_ratio=max(branch_ratio,m.physical_branch_ratio)
        source_rows.append(output)
        np.savez_compressed(OUT/f'material-source-{steps}-reference-{ref}.npz',t=m.t,radius_E=m.rE,reference_lapse=m.a,volume=m.V,
            additional_gas_energy_erg=output[:,0],additional_gas_radial_pressure_erg=output[:,1],additional_gas_trace_erg=output[:,2],additional_gas_tangential_pressure_erg=output[:,3],
            total_integrated_gas_energy_erg=proper[:,0],total_integrated_gas_radial_pressure_erg=proper[:,1],total_integrated_gas_tangential_pressure_erg=proper[:,2],delta_baryon_g=hist[:,0]*AMP,
            delta_radial_momentum_c_erg=hist[:,1]*AMP,delta_reference_material_energy_erg=hist[:,2]*AMP,delta_neutral_number=hist[:,3]*AMP)
    def compare(a,b):return (np.max(np.sum(abs(a-b),axis=2),axis=0)/np.maximum(np.max(np.sum(abs(b),axis=2),axis=0),1e-300)).tolist()
    comparisons=dict(time=compare(histories[0],histories[1]),background=compare(histories[2],histories[1]),stress_time=compare(source_rows[0],source_rows[1]),stress_background=compare(source_rows[2],source_rows[1]),
        endpoint_time=compare(histories[0][-1:],histories[1][-1:]),endpoint_background=compare(histories[2][-1:],histories[1][-1:]))
    errors=[v for row in comparisons.values() for v in row];pressure_error=max(v for row in pressure_checks for v in row)
    # Integrated orthonormal momentum: conserved covariant momentum is h*V*Jhat.
    h,P,hdot,Pdot=sp.symbols('h P hdot Pdot');assert sp.simplify((hdot*P+h*Pdot).subs(Pdot,-hdot/h*P))==0
    write(OUT/'momentum-symbolic.json',dict(classification='Proven',passed=True,scope='Homogeneous zero-shift cell: conservation of h*V*Jhat gives d(V*Jhat)/dt=-(h_t/h)*(V*Jhat). This is not a full nonlinear numerical GR theorem.'))
    result=dict(classification='Counterexample candidate',passed=bool(max(errors)<.02 and balance<1e-8 and pressure_error<.002),comparisons=comparisons,conserved_order=['baryon','momentum_c','reference_energy','neutral_H'],stress_order=['energy','radial_pressure','trace','tangential_pressure'],
        independent_balance_relative=balance,pressure_probe_relative=pressure_checks,maximum_physical_branch_ratio=branch_ratio,
        omitted_deep_kinetic_metric_work_sampled_envelope_erg=deep_missing_work,
        endpoint_additional_source_sum=np.sum(source_rows[1][-1],axis=1).tolist(),maximum_additional_source_L1=np.max(np.sum(abs(source_rows[1]),axis=2),axis=0).tolist(),
        additional_material_motion_evolved=True,additional_GR_source_exported=True,additional_GR_source_applied=False,full_photon_material_feedback=False,final_charge_solved=False,full_goal_complete=False,seconds=time.monotonic()-start)
    write(OUT/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
