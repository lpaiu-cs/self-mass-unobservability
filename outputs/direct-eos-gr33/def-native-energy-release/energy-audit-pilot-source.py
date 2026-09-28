"""Independent native-state and global-energy checks on completed paths."""
import json
import time
import numpy as np
import sympy as s
import def_native_energy_flow as task


def main():
    out=task.OUT;assert not (out/'energy-audit.json').exists();start=time.monotonic()
    paths=[p for p in [out/'pilot-224.npz',out/'cells-896.npz',out/'cells-1792.npz'] if p.exists()]
    assert paths,'No completed energy path'
    task.write(out/'energy-audit-plan.json',dict(classification='Counterexample candidate',seconds=15,native_calls=18,
        paths=[p.name for p in paths],gates=dict(native_state=.002,energy_ledger=1e-8),
        claim='Check accepted energy/baryon/tail/inner/photon histories and direct native EOS at actual final states. This does not certify chemical equilibrium or full metric feedback.'))
    aa,a0,cx,D,E,p,v=s.symbols('a a0 cx D E p v');K=aa*(E-cx*D)+(aa-a0)*cx*D
    assert s.expand((K+aa*p)*v-(aa*(E+p)*v-a0*cx*D*v))==0
    fan=task.bank.task.prior.Fan(call_cap=18,reuse=True);rows=[]
    for path in paths:
        d=np.load(path);n=d['U'].shape[1];m=task.bank.prior.Flow(n);history=d['history'];scale=4*np.pi*m.RJ**2*m.eos.rho0
        # Independent long-double summation of saved conserved variables;
        # never reconstruct total rest energy and subtract it afterward.
        delta=d['U'].astype(np.longdouble)-d['initial'].astype(np.longdouble)
        eb=float(np.sum(delta[2]*d['volume'],dtype=np.longdouble)+d['conserved_discard'][2]-history[-1,4])
        bb=float(np.sum(delta[0]*d['volume'],dtype=np.longdouble)+d['conserved_discard'][0]-history[-1,3])
        active=np.flatnonzero(d['rho']>=m.eos.rho0*m.eos.floor)
        logs=np.log(d['rho'][active]/m.eos.rho0);ids=np.unique([active[np.argmin(abs(logs-z))] for z in [0,-2,-5,-10,-16.5]]+[active[np.argmax(abs(d['velocity_cm_s'][active]))]])
        controls=[]
        # Read the final cached temperature representation without extending it.
        eos=task.temperature.parent.EOS(out/'runtime-columns.npz')
        for i in ids:
            raw=fan.call(float(np.log(d['rho'][i])),float(np.log(d['T'][i])))
            pp,uu,gg,tt,kk=eos(np.array([d['rho'][i]/m.eos.rho0]),np.array([np.log(d['T'][i])]))
            errors=[float(abs(pp[0]*m.eos.rho0*task.C**2/raw[1]-1)),float(abs(uu[0]*task.C**2/raw[2]-1)),float(abs(gg[0]/raw[4]-1))]
            controls.append(dict(cell=int(i),log_density=float(np.log(d['rho'][i]/m.eos.rho0)),T_K=float(d['T'][i]),relative=errors))
        old=json.loads(path.with_suffix('.json').read_text());response=max(abs(float(history[-1,2])),float(np.sum(d['volume']*(d['rho']/m.eos.rho0*(d['velocity_cm_s']/task.C)**2+d['pressure']/(m.eos.rho0*task.C**2)))))
        rows.append(dict(cells=n,completed=True,independent_energy_residual_erg=eb*scale*task.C**2,
            independent_energy_relative=abs(eb)/response,independent_baryon_residual_g=bb*scale,
            energy_into_bulk_erg=-old['inner_Killing_nonrest_energy_into_layer_erg'],photon_energy_change_erg=-old['scattering_work_into_gas_erg'],
            native_controls=controls,passed=abs(eb)/response<1e-8 and max(max(c['relative']) for c in controls)<.002))
    result=dict(classification='Counterexample candidate',passed=all(r['passed'] for r in rows),rows=rows,native_calls=fan.calls,seconds=time.monotonic()-start,
        fixed_metric_energy_identity='Proven',energy_defect_added_as_heat=False,full_GR_scalar_feedback=False,full_goal_complete=False)
    task.write(out/'energy-audit.json',result);print(json.dumps(result),flush=True);assert result['passed']


if __name__=='__main__':main()
