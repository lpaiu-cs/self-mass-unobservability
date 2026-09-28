"""Resolve the failed trace comparison into stored conservative budgets."""
import json
import numpy as np
import def_native_reactive_flow as flow


def main():
    out=flow.exchange.OUT;rows=[]
    for n in [448,896]:
        f=flow.Flow(n);e=f.eos;m=f.base;d=np.load(flow.OUT/f'cells-{n}.npz');scale=4*np.pi*m.RJ**2*e.rho0*flow.C**2
        terms=[]
        for U in [d['initial'],d['U']]:
            rho,v,t,y=f.primitive(U);p,u,*_=e(rho,t);root=np.sqrt(1-v*v);W=1/root;wm=v*v/(root*(1+root))
            u0=(e.u0[0]+(e.u0[1]-e.u0[0])*(y-e.ys[0])/np.diff(e.ys)[0])/flow.C**2
            du=u-u0-1.5*p/np.maximum(rho,e.floor)
            remainder=(e.cx+u0)*rho*wm**2+2.5*p*(W*W-1)+rho*du*(1+W*W)
            nr=(U[2]-(m.a-m.a0)*e.cx*U[0])/m.a
            trace=-e.cx*U[0]*v*v/(1+root)+rho*u-3*p
            assert max(abs(trace+nr-2*u0*U[0]-remainder))<1e-17
            terms.append([float(np.sum(z*m.vol)*scale) for z in [-nr,2*u0*U[0],remainder,trace]])
        delta=np.subtract(terms[1],terms[0]);h=d['history'][-1];discard=d['conserved_discard']
        photon=-h[5]*scale/m.a0;port=-(h[4]-h[5])*scale/m.a0;cut=discard[2]*scale/m.a0
        lapse=delta[0]-(photon+port+cut)
        pieces=dict(photon=photon,inner_boundary=port,discarded_Killing_energy=cut,lapse_conversion=lapse,chemical_inventory=delta[1],nonideal_relativistic_remainder=delta[2])
        reconstructed=sum(pieces.values());reported=json.loads((flow.OUT/f'cells-{n}.json').read_text())['integrated_trace_energy_erg']
        error=abs(reconstructed-reported)/abs(reported);assert error<1e-8
        rows.append(dict(cells=n,terms_erg=pieces,reconstructed_trace_erg=reconstructed,reported_trace_erg=reported,relative=error))
    differences={key:rows[1]['terms_erg'][key]-rows[0]['terms_erg'][key] for key in rows[0]['terms_erg']}
    total=rows[1]['reported_trace_erg']-rows[0]['reported_trace_erg']
    result=dict(classification='Counterexample candidate',passed=True,rows=rows,fine_minus_coarse_erg=differences,total_difference_erg=total,
        signed_fraction_of_grid_difference={key:value/total for key,value in differences.items()},
        interpretation='An exact saved-budget decomposition, not a corrected trace or a new accepted verdict. The largest differing term identifies the next interface/energy target; it does not prove a particular boundary prescription is physically correct.',
        registered_flow_verdict_preserved=False,final_charge_solved=False)
    flow.write(out/'trace-budget-decomposition.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
