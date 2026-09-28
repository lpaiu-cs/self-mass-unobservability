"""Saved-path audit using direct quadrature sources and adjoint reciprocity."""
from pathlib import Path
import json
import signal
import time
import numpy as np
import sympy as sp
from scipy.sparse import diags,coo_matrix
import def_orbital_conductive_feedback as run


def main():
    out=run.OUT;assert not (out/'audit.json').exists();signal.alarm(90);start=time.monotonic()
    plan=json.loads((out/'plan.json').read_text())
    for p,h in plan['bindings'].items():assert run.go.task.digest(Path(p))==h,p
    # Algebraic reaction: A qa=-kL, A dq=B E, A=A^T.
    A,b,f,k,E=sp.symbols('A b f k E',nonzero=True)
    qa=-k/A;dq=b*E/A
    assert sp.simplify(k*dq-f*E+(qa*b+f)*E)==0
    q,u,H,Tq,TE=sp.symbols('q u H Tq TE')
    assert sp.expand((Tq*u+(TE-Tq*H)*E).subs(u,q+H*E)-(Tq*q+TE*E))==0
    p=run.Problem();m=p.model;active=p.active;data=m.data
    src=run.go.task.fem.source_points(m.heat,m.points)
    B=data[:,4:8].reshape(-1,2,2);gs=data[:,16:22].reshape(-1,2,3);hs=data[:,22:28].reshape(-1,2,3)
    g=[sum(diags(gs[:,i,j])@src[j] for j in range(3)) for i in range(2)]
    h=[sum(diags(hs[:,i,j])@src[j] for j in range(3)) for i in range(2)]
    # Reassemble F directly. Do not recover it by subtracting K H from load.
    direct=sum(m.cov[i].astype(np.clongdouble).T@diags(m.weights.astype(np.longdouble)*B[:,i,j])@g[j].astype(np.clongdouble) for i in range(2) for j in range(2))
    direct-=sum(m.V[i].astype(np.clongdouble).T@diags(m.weights.astype(np.longdouble))@h[i].astype(np.clongdouble) for i in range(2))
    direct=direct[:,active].tocsc()
    bank=np.load(run.go.task.BANK/'fine-bank.npz');faces=bank['faces'];nn=len(p.theta)
    diff=coo_matrix((np.tile([-1.,1.],len(faces)),(np.repeat(np.arange(len(faces)),2),np.c_[nn-faces,nn-1-faces].ravel())),shape=(len(faces),nn)).tocsr()
    GE=(diff@diags(p.theta)@p.TE)[:,active]
    rows=[]
    for n in [1,2,3]:
        saved=np.load(out/f'p4-{n}.npz');qa=saved['qa'];dq=saved['dq'];E=saved['E'];z=np.clongdouble(-1j*n*p.omega)
        op=direct-z*z*(p.Mx@p.H)
        residual=p.Kx@dq+z*z*(p.Mx@dq)-op@E
        error=float(np.max(abs(residual)/(p.absK@abs(dq)+abs(z*z)*(p.absM@abs(dq))+abs(op)@abs(E)+1e-100)))
        direct_reaction=(p.kl+z*z*p.ml)@dq-p.fL@E
        adjoint_reaction=-(op.T@qa+p.fL)@E
        relative=float(abs(direct_reaction-adjoint_reaction)/abs(adjoint_reaction))
        gain=np.sum(p.conductance*m.heat.geometry.tc*p.lam/(z*(z+p.lam)),axis=1)
        defect=E-gain*(p.Gq@(qa+dq)+p.gL+GE@E)
        heat_error=float(abs(defect).max()/abs(E).max())
        cancellation=float((abs((p.kl+z*z*p.ml)@dq)+abs(p.fL@E))/abs(direct_reaction))
        rows.append(dict(harmonic=n,direct_quadrature_GR_residual=error,independent_heat_map_residual=heat_error,
            adjoint_reaction_relative=relative,reaction_cancellation_ratio=cancellation))
    result=json.loads((out/'result.json').read_text())
    passed=all(r['direct_quadrature_GR_residual']<1e-9 and r['independent_heat_map_residual']<1e-9 and r['adjoint_reaction_relative']<1e-5 for r in rows)
    row=dict(classification='Counterexample candidate',passed=passed,rows=rows,symbolic_identities_passed=True,
        sum_three_harmonic_charge_contributions=sum(r['delta_alpha_radiative_over_phi0'] for r in result['cases']['p4']),
        seconds=time.monotonic()-start,total_accounted_compute_seconds=time.monotonic()-start+result['total_compute_seconds'],
        full_objective_complete=False,scope='Direct original heat quadrature and reciprocity audit of the stored coupled solutions; not an EOS/continuum certificate or full thermal-background solution.')
    run.write(out/'audit.json',row);print(json.dumps(row),flush=True)
    assert passed


def project():
    import def_orbital_charge_audit as check
    check.mp.mp.dps=70;out=run.OUT
    thermal=json.loads((out/'result.json').read_text())['cases']['p4']
    old=json.loads((run.orbit.OUT.parent/'def-full-orbital-exterior/comparator-result.json').read_text())
    baseline=old['cases']['p4-tight'];combined=[baseline[0]];heat=[dict(charge=[0,0])]
    for a,b in zip(baseline[1:],thermal):
        total=[check.mp.mpf(str(x))+check.mp.mpf(str(y)) for x,y in zip(a['charge'],b['radiative_charge_gain'])]
        combined.append(dict(charge=[check.mp.nstr(v,60) for v in total]));heat.append(dict(charge=b['radiative_charge_gain']))
    amplitudes=[check.mp.mpf('0')]+[check.mp.mpf(str(r['actual_drive_amplitude'])) for r in thermal]
    phi0=check.mp.mpf('.001')
    norm=check.mp.sqrt(sum(abs(check.mp.mpc(*[str(v) for v in r['charge']])*amplitudes[n]/phi0)**2 for n,r in enumerate(heat)))
    fits=[check.projection(combined,d,amplitudes,phi0) for d in [2,3,4,5]]
    contribution=check.projection(heat,4,amplitudes,phi0)
    assert contribution['residual_norm']<=float(norm)*(1+1e-12)
    comparison=old['degree4_control_difference']
    assert fits[2]['residual_norm']<comparison
    value=dict(classification='Counterexample candidate',combined_fits=fits,
        unprojected_thermal_charge_norm=float(norm),degree4_thermal_projection=contribution,
        previous_degree4_control_difference=comparison,
        thermal_norm_over_previous_control_difference=float(norm)/comparison,
        whole_response_degree4_nonabsorption_resolved=False,full_objective_complete=False,
        identity='Proven: an orthogonal least-squares nuisance projection is nonexpansive. This known conductive addition cannot increase its six-coefficient residual norm by more than the unprojected addition.',
        scope='Same frozen three-harmonic charge coefficient norm and comparator; no timing likelihood, missing outer transport or thermal-background bound.')
    run.write(out/'combined-projection.json',value);print(json.dumps(value),flush=True)


if __name__=='__main__':
    import sys
    project() if len(sys.argv)>1 and sys.argv[1]=='project' else main()
