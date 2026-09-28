"""Read-only checks of saved local heat paths and the native initial local lift."""
from pathlib import Path
import json
import numpy as np
import def_heat_coupling as model


def finalize_saved():
    plan=model.bindings();failure=json.loads((model.OUT/'failure-result.json').read_text())
    assert 'not JSON serializable' in failure['error']
    read=lambda name:json.loads((model.OUT/name).read_text())
    rows=[read(f'steps-{k}-damping-1-zeroheat-0.json') for k in [8,16,32]]
    zero=read('steps-8-damping-1-zeroheat-1.json');reversible=read('steps-8-damping-0-zeroheat-0.json')
    endpoints=np.array([r['endpoint'] for r in rows]);diff=np.max(abs(np.diff(endpoints,axis=0)),axis=1)
    order=float(np.log2(diff[0]/diff[1]))
    gates=dict(time_order=order>=plan['gates']['refinement_order'],
        entropy_positive=all(r['entropy_change_scaled']>0 for r in rows),
        entropy_budget=rows[-1]['entropy_budget_relative']<=plan['gates']['entropy_budget_relative'],
        scalar_turning=all(r['momentum_crossings']>0 for r in rows),
        zero_heat_invariant=bool(max(abs(np.array(zero['endpoint'])[1:4]))<1e-12))
    result=dict(classification='Counterexample candidate',passed=all(gates.values()),gates=gates,rows=rows,
        time_order=order,endpoint_differences=diff.tolist(),zero_heat=zero,reversible=reversible,
        native_calls=failure['native_calls'],seconds=failure['seconds'],full_GR=False,
        finite_step_entropy_exact=False,physical_EOS_certified=False,observational_closure=False,
        execution_repair='All five paths were saved before a numpy.bool_ serialization failure. Reconstruct the unchanged original gates and cast that one flag to bool. No trajectory rerun or scientific threshold change.',
        finalizer_sha256=model.digest(Path(__file__)))
    model.save('result.json',result)
    model.save('completion-manifest.json',dict(sha256={p.relative_to(model.ROOT).as_posix():model.digest(p)
        for p in model.OUT.iterdir() if p.is_file() and p.name!='completion-manifest.json'}))
    assert result['passed'],result


def native_lift():
    saved=np.load(model.INITIAL);base=saved['base'];aux=saved['aux']
    rho,T=np.exp(base[:,:2]).T;w=saved['qscale'];ct=rho*aux[:,10]
    q=np.column_stack([base[:,3]-base[:,4],base[:,4]])*w[:,None]
    K=16*model.ld('5.670400e-5')*T[:,None]**3/(3*rho[:,None]*aux[:,[27,24]])
    tau=np.column_stack([np.full(len(w),model.two.e.TAU),1/(model.prior.C*rho*aux[:,24])])
    C=K*T[:,None]/(model.prior.C**2*tau)
    br=np.column_stack([1+aux[:,28],np.zeros(len(w))])
    bt=np.column_stack([-5+aux[:,29],np.full(len(w),-5)])
    Q=q.sum(axis=1);cr=rho*aux[:,9]-aux[:,1]
    # v=0, imposed d ln a/dt=1, reversible heat part, phi=0.
    # Unknowns: d ln T/dt, dv/dt, dqcond/(w dt), dqrad/(w dt).
    M=np.zeros((len(w),4,4),dtype=model.ld)
    M[:,0,0]=1;M[:,0,1]=2*Q/ct
    M[:,1,1:]=1
    for j in range(2):
        M[:,2+j,0]=q[:,j]*bt[:,j]/(2*w)
        M[:,2+j,1]=C[:,j]/w;M[:,2+j,2+j]=1
    rhs=np.column_stack([cr/ct,-2*Q/w,q*(br-1)/(2*w[:,None])])
    rates=np.linalg.solve(M.astype(float),rhs.astype(float)[:,:,None])[:,:,0].astype(model.ld)
    residual=float(np.max(abs(np.einsum('nij,nj->ni',M,rates)-rhs)))
    margin=(w-C.sum(axis=1)+Q*np.sum(q*bt,axis=1)/ct)/w
    # Energy identity, evaluated with physical material derivatives.
    edot=-(w+cr)+ct*rates[:,0]+2*Q*rates[:,1]
    energy_residual=float(np.max(abs(edot+w)/w))
    assert np.all(base[:,2]==0) and np.min(margin)>0 and residual<1e-12 and energy_residual<1e-12
    return dict(classification='Counterexample candidate',passed=True,cells=len(w),
        minimum_scaled_rest_determinant=float(np.min(margin)),maximum_local_linear_residual=residual,
        maximum_relative_energy_identity_residual=energy_residual,
        scope='All saved native initial states, local reversible homogeneous metric deformation only. No spatial evolution, new EOS call or global well-posedness claim.')


def replay():
    plan=model.bindings();result=json.loads((model.OUT/'result.json').read_text());assert result['passed']
    for rel,sha in json.loads((model.OUT/'completion-manifest.json').read_text())['sha256'].items():
        assert model.digest(model.ROOT/rel)==sha,rel
    cell=model.Cell();rows=result['rows']+[result['zero_heat'],result['reversible']]
    max_work=0.;max_entropy_difference=0.;paths=0
    for row in rows:
        history=np.load(model.OUT/row['history_file'])['history']
        y0=cell.y0.copy()
        if row['zero_heat']:y0[2:4]=0
        first=cell.state(y0,0)
        for point in history:
            z=cell.state(point[1:6].astype(model.ld),model.ld(point[0]))
            p=model.ld(point[6])*cell.phiscale
            scalar=cell.I/2*(p*p+(z['phi']+first['phi'])*(z['phi']-first['phi']))
            work_defect=abs((cell.energy_increment(first,z)+scalar)/cell.thermal-point[7])
            entropy=(z['entropy_cell']-first['entropy_cell'])*first['T']/cell.thermal
            max_work=max(max_work,float(work_defect));max_entropy_difference=max(max_entropy_difference,float(abs(entropy-point[8])))
        paths+=1
    # Histories were stored in binary64; fresh endpoint reconstruction is not bitwise replay.
    assert max_work<1e-10 and max_entropy_difference<1e-10,(max_work,max_entropy_difference)
    points=np.array([np.r_[np.load(model.OUT/row['history_file'])['endpoint'],
        np.load(model.OUT/row['history_file'])['history'][-1,6]] for row in result['rows']])
    differences=np.max(abs(np.diff(points,axis=0)),axis=1)
    order=float(np.log2(differences[0]/differences[1]))
    assert abs(order-result['time_order'])<1e-12
    return dict(classification='Counterexample candidate',passed=True,paths=paths,
        maximum_reconstructed_work_defect=float(max_work),maximum_entropy_record_difference=float(max_entropy_difference),
        observed_time_order=order,fresh_native_calls=cell.provider.calls,
        scope='Stored binary64 path states reconstructed with fresh native EOS/opacity; energy, entropy records and endpoint refinement checked. Does not independently recompute every nonlinear stage residual.')


if __name__=='__main__':
    if not (model.OUT/'result.json').exists():finalize_saved()
    target=model.OUT/'replay.json';assert not target.exists()
    value=dict(symbolic=model.symbolic(),native_initial_lift=native_lift(),saved_paths=replay(),
        source_sha256=model.digest(Path(__file__)))
    target.write_text(json.dumps(value,indent=2)+'\n');print(json.dumps(value,indent=2))
