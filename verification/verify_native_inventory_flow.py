"""Audit conservative histories and the missing-tail trace identity."""
import json
import time
import numpy as np
import sympy as s
import def_native_inventory_flow as task


def main():
    out=task.OUT;assert not (out/'audit.json').exists();begin=time.monotonic()
    rho,p,c,u0,du,W=s.symbols('rho p c u0 du W',positive=True)
    u=u0+s.Rational(3,2)*p/rho+du;D=rho*W
    energy=(rho*(c+u)+p)*W**2-p-c*D
    trace=rho*(c+u)-3*p-c*D
    defect=(c+u0)*rho*(W-1)**2+s.Rational(5,2)*p*(W**2-1)+rho*du*(1+W**2)
    assert s.simplify(trace+energy-2*u0*D-defect)==0
    task.write(out/'symbolic.json',dict(classification='Proven',passed=True,
        definitions='D=rho*W; E_nr=(rho*(c_rest+u)+p)*W^2-p-c_rest*D; T_nr=rho*(c_rest+u)-3*p-c_rest*D; u=u0+3p/(2rho)+du.',
        identity='T_nr=-E_nr+2*u0*D+(c_rest+u0)*rho*(W-1)^2+(5/2)*p*(W^2-1)+rho*du*(1+W^2).',
        limit='For fixed u0, monatomic ideal gas and nonrelativistic flow, the volume-integrated trace is determined by conserved energy and mass up to order v^4 rest and v^2 thermal corrections. Deleted dilute inventory must participate in that identity.',
        scope='Local algebra in one orthonormal frame. This is not a retarded scalar-charge identity; lapse, geometric weights, boundary histories and tail locations matter.'))
    eos=task.EOS();rows=[];native=task.Native(cap=120);checks=[]
    for n in [224,448,896]:
        label=('pilot-' if n==224 else 'cells-')+str(n);d=np.load(out/(label+'.npz'));flow=task.Flow(n);m=flow.base
        h=d['history'];discard=d['conserved_discard'];vol=m.vol;scale=4*np.pi*m.RJ**2*eos.rho0
        mass=np.sum((d['U'][0]-d['initial'][0])*vol)+discard[0]-h[-1,3]
        en=np.sum((d['U'][2]-d['initial'][2])*vol)+discard[2]-h[-1,4]
        r=d['rho']/eos.rho0;v=d['velocity_cm_s']/task.C;theta=d['logT'];p,u,_,T,_,_,_=eos.evaluate(r,theta)
        active=r>=eos.floor;root=np.sqrt(1-v*v);W=1/root;wm=v*v/(root*(1+root))
        specific_defect=u-eos.u0/task.C**2-1.5*p/np.maximum(r,eos.floor)
        remainder=(eos.cx+eos.u0/task.C**2)*r*wm**2+2.5*p*(W*W-1)+r*specific_defect*(1+W*W)
        Enr=(d['U'][2]-(m.a-m.a0)*eos.cx*d['U'][0])/m.a
        Tnr=r*(eos.cx+u)-3*p-eos.cx*d['U'][0]
        identity_error=float(np.max(abs(Tnr+Enr-2*eos.u0/task.C**2*d['U'][0]-remainder)))
        # Evaluate the trace without subtracting two O(rest energy) numbers.
        stable_trace=-eos.cx*d['U'][0]*v*v/(1+root)+r*u-3*p
        stable_identity=float(np.max(abs(stable_trace+Enr-2*eos.u0/task.C**2*d['U'][0]-remainder)))
        assert stable_identity<1e-17
        result=json.loads((out/(label+'.json')).read_text());dm=discard[0]*scale;dk=discard[2]*scale*task.C**2
        # Diagnostic completion under a common lapse and zero defect only.
        # Do not turn this into the registered failed spatial verdict.
        correction=2*eos.u0*dm-dk/m.a0
        rows.append(dict(cells=n,baryon_residual=float(mass),energy_residual=float(en),
            point_identity_absolute=identity_error,stable_identity_absolute=stable_identity,
            resolved_trace_erg=result['integrated_trace_energy_erg'],discarded_baryon_g=float(dm),
            common_lapse_ideal_tail_trace_erg=float(correction),
            diagnostic_completed_trace_erg=float(result['integrated_trace_energy_erg']+correction),
            discarded_to_outside_baryon=float(dm/result['gas_outside_original_radius_g']),
            maximum_saved_velocity_over_c=float(max(abs(v))),maximum_saved_specific_nonideal_defect_erg_g=float(max(abs(specific_defect[active]))*task.C**2)))
        assert result['passed'] and abs(mass)/np.sum(d['initial'][0]*vol)<1e-10
        if n==896:
            # Independent actual final states, not the bank's preselected controls.
            ids=np.flatnonzero(active);chosen=ids[np.linspace(0,len(ids)-1,5).astype(int)]
            for i in chosen:
                a=native(float(np.log(r[i])),float(theta[i]));errors=[abs(p[i]*eos.rho0*task.C**2/a[1]-1),abs(u[i]*task.C**2/a[2]-1)]
                checks.append(dict(cell=int(i),rho=float(d['rho'][i]),T=float(T[i]),relative=list(map(float,errors))))
    assert max(max(c['relative']) for c in checks)<.002
    contrast=abs(rows[-2]['diagnostic_completed_trace_erg']/rows[-1]['diagnostic_completed_trace_erg']-1)
    task.write(out/'trace-tail-diagnosis.json',dict(classification='Counterexample candidate',rows=rows,
        common_lapse_ideal_tail_relative=float(contrast),u0_erg_g=eos.u0,
        interpretation='The stored deleted inventory carries trace comparable to the resolved signal. The exact local identity explains why resolved-only trace convergence fails. Common lapse and ideal-tail completion are diagnostic assumptions, not a physical tail solution or a post-result rescue of the frozen gate.',
        original_spatial_verdict_preserved=False if json.loads((out/'result.json').read_text())['passed'] else True,
        physical_tail_certified=False,retarded_charge_computed=False))
    task.write(out/'actual-state-audit.json',dict(classification='Counterexample candidate',passed=True,checks=checks,
        EOS_calls=native.ion.calls,source_sha256=task.cold.sha(__file__)))
    # Save exact native calls for independent state reuse.
    np.savez_compressed(out/'actual-state-native.npz',**{k:np.array([v[k] for v in native.ion.states]) for k in native.ion.states[0]})
    result=dict(classification='Counterexample candidate',passed=True,conserved_histories_checked=True,
        sampled_actual_native_states_passed=True,registered_spatial_trace_passed=False,seconds=time.monotonic()-begin,
        full_goal_complete=False)
    task.write(out/'audit.json',result);print(json.dumps(dict(audit=result,rows=rows,controls=checks)),flush=True)


def dilute():
    out=task.OUT/'dilute';assert not (out/'audit.json').exists();start=time.monotonic();eos=task.DiluteEOS();native=task.Native(cap=100);checks=[];rows=[]
    for n in [448,896]:
        d=np.load(out/f'cells-{n}.npz');flow=task.DiluteFlow(n);m=flow.base;h=d['history'];discard=d['conserved_discard'];vol=m.vol
        mass=np.sum((d['U'][0]-d['initial'][0])*vol)+discard[0]-h[-1,3]
        energy=np.sum((d['U'][2]-d['initial'][2])*vol)+discard[2]-h[-1,4]
        relative=abs(mass)/np.sum(d['initial'][0]*vol)
        result=json.loads((out/f'cells-{n}.json').read_text());assert result['passed'] and relative<1e-10
        rho=d['rho']/eos.rho0;p,u,*_=eos(rho,d['logT']);active=rho>=eos.floor
        rows.append(dict(cells=n,baryon_relative=float(relative),energy_residual=float(energy),discarded_to_outside_baryon=result['dilute_baryon_g']/result['gas_outside_original_radius_g']))
        if n==896:
            ids=np.flatnonzero(active);chosen=ids[np.linspace(0,len(ids)-1,6).astype(int)]
            for i in chosen:
                a=native(float(np.log(rho[i])),float(d['logT'][i]));errors=[abs(p[i]*eos.rho0*task.C**2/a[1]-1),abs(u[i]*task.C**2/a[2]-1)]
                checks.append(dict(cell=int(i),rho=float(d['rho'][i]),T=float(d['T'][i]),relative=list(map(float,errors))))
    assert max(max(c['relative']) for c in checks)<.002
    old=json.loads((task.previous.OUT/'cells-896.json').read_text());new=json.loads((out/'cells-896.json').read_text())
    comparison={k:dict(frozen=new[k],LTE=old[k],ratio=new[k]/old[k]) for k in ['gas_outside_original_radius_g','integrated_trace_energy_erg','scattering_work_into_gas_erg']}
    result=dict(classification='Counterexample candidate',passed=True,rows=rows,actual_native_checks=checks,
        comparison_same_896_grid=comparison,EOS_calls=native.ion.calls,seconds=time.monotonic()-start,
        earlier_failed_cutoff_verdict_preserved=True,physical_chemistry_closed=False,final_charge_solved=False,full_goal_complete=False,
        source_sha256=task.cold.sha(__file__))
    task.write(out/'audit.json',result)
    np.savez_compressed(out/'audit-native.npz',**{k:np.array([v[k] for v in native.ion.states]) for k in native.ion.states[0]})
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    import sys
    globals()[sys.argv[1] if len(sys.argv)>1 else 'main']()
