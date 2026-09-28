"""Nonuniform fixed-metric cell-moment inverse with native EOS at every node."""
import json,sys
from functools import lru_cache
from types import FunctionType
import numpy as np
import gr_subcell_reference as reference
import gr_heat_primitive_newton as point

g=reference.g;OUT=g.OUT/'gr-subcell-inverse';traces=[]


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();point.verify()
    plan=json.loads((point.OUT/'plan.json').read_text())
    plan.update(checkpoint='6338252',cells=[0,1,2,2972],nodes=16,tolerance_index=1,
        change='Retain nonuniform rho0,T0 and positive coordinate/proper-volume weights. Apply constant cell log-density and log-temperature shifts and a common cell velocity; fixed Q equals the saved midpoint heat density on all nodes. Native EOS is evaluated on every node for each residual/Jacobian. Invert B, U-C*kappa*B, Pi under the same decreasing Newton and original error tolerances.',
        scale='Use the unperturbed node heat-capacity and enthalpy integrals for residual scales; no truth-derived solver scale.',
        startup_gate='The full native node export must complete; all four 16-node fine-tolerance baryon, mass and internal-energy comparisons must pass their original gates.',
        scope='Three-parameter nonuniform reconstruction in four fixed-metric cells. Constructed moment targets vary baryon inventory in a fixed geometry for inverse testing; this is not a material time path. No uniform-cell substitution, full-star evolution, arbitrary subcell closure or continuous root/EOS error certificate.')
    plan['bindings'].update({p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in [
        g.ROOT/'verification/gr_subcell_inverse.py',reference.OUT/'plan.json',point.OUT/'manifest.json',
        g.OUT/'gr-heat-entropy-closure/initial-rates.npz']})
    save('plan.json',plan)


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    reference.verify()
    quadrature=json.loads((reference.OUT/'subcell-quadrature.json').read_text())['rows']
    accepted=[r for r in quadrature if r['nodes']==plan['nodes'] and r['tolerance_index']==plan['tolerance_index']]
    assert len(accepted)==4 and all(r['passed'] for r in accepted)
    save('startup-binding.json',dict(classification='Proven',reference_manifest_sha256=g.c.sha(reference.OUT/'manifest.json')))
    state=dict(np.load(g.OUT/'initial-state-17-4.npz'));Qvalues=np.load(g.OUT/'gr-heat-entropy-closure/initial-rates.npz')['Q']
    eos=g.EOS();ld=np.longdouble;rows=[];evaluations=0
    solve=FunctionType(point.root.__code__,dict(point.root.__globals__,OUT=OUT,save=save,traces=traces))
    for i in plan['cells']:
        data=dict(np.load(reference.OUT/f"cell-{i}-tolerance-{plan['tolerance_index']}-nodes-{plan['nodes']}.npz"))
        rho0=data['eos'][:,0].astype(ld);T0=data['lnT'];weights=data['coordinate_weights_cm3'].astype(ld)
        proper=data['proper_weights_cm3'].astype(ld);B0=np.sum(rho0*proper)
        C=ld(data['C_X'])*(ld(g.c.gr.C)*100)**2;Q=ld(Qvalues[i]);X=state['X'][i]
        capacity=np.sum(rho0*data['eos'][:,10]*weights)
        enthalpy=np.sum((rho0*(C+data['eos'][:,2])+data['eos'][:,1])*proper)
        scale=np.array([capacity,enthalpy]);assert np.all(scale>0)
        def values(rho,theta,v):
            nonlocal evaluations
            a=np.array([eos(2,float(np.log(r)),float(t+theta),X) for r,t in zip(rho,T0)])
            evaluations+=len(rho);assert np.all(a[:,10]>0)
            moments=np.array([point.original.forward(r,b,v,Q,C) for r,b in zip(rho,a)])
            return np.array([np.sum(moments[:,0]*proper),np.sum(moments[:,1]*weights),np.sum(moments[:,2]*proper)]),a
        for case,(eta,theta,velocity) in enumerate(plan['probes']):
            target,truth=values(rho0*np.exp(ld(eta)),theta,velocity);Dprofile=rho0*(target[0]/B0)
            @lru_cache(maxsize=1)
            def evaluate(theta_new,v):
                W=1/np.sqrt(1-ld(v)**2);rho=Dprofile/W;value,a=values(rho,theta_new,v)
                residual=np.asarray((value[1:]-target[1:])/scale,float)
                P=a[:,1];w=rho*(C+a[:,2])+P;er=rho*(C+a[:,2]+a[:,9]);et=rho*a[:,10]
                pr=P*a[:,5];pt=P*a[:,6]
                j00=np.sum(W*W*(et+pt*v*v)*weights)
                j01=np.sum(W**4*(2*(v*w+Q*(1+v*v))-v*(er+pr*v*v))*weights)
                j10=np.sum(v*W*W*(et+pt)*proper)
                j11=np.sum(W**4*((1+v*v)*w+4*v*Q-v*v*(er+pr))*proper)
                jac=np.asarray(np.array([[j00,j01],[j10,j11]])/scale[:,None],float)
                return residual,jac,value,a
            solution=solve(lambda x:evaluate(*map(float,x))[0],[theta+.003,0.],
                jac=lambda x:evaluate(*map(float,x))[1],method='coupled Newton',options={})
            residual,jac,recovered,a=evaluate(*map(float,solution.x));t,v=solution.x
            recovered_eta=float(np.log(target[0]/B0)+ld(.5)*np.log1p(-ld(v)**2))
            errors=np.abs(np.array([recovered_eta-eta,t-theta,v-velocity]))
            row=dict(cell=i,case=case,log_density_shift_error=float(errors[0]),log_temperature_shift_error=float(errors[1]),
                velocity_error=float(errors[2]),baryon_relative_error=float(abs(recovered[0]/target[0]-1)),
                scaled_moment_residual=float(abs(residual).max()),solver_success=bool(solution.success))
            row['passed']=bool(errors[0]<=plan['root_log_density_tolerance'] and errors[1]<=plan['root_logT_tolerance'] and
                errors[2]<=plan['root_velocity_tolerance'] and row['baryon_relative_error']<=plan['root_log_density_tolerance'] and
                row['scaled_moment_residual']<=plan['scaled_residual_tolerance'])
            rows.append(row);np.savez_compressed(OUT/f'cell-{i}-case-{case}.npz',target=target,recovered_moments=recovered,
                true_parameters=np.array([eta,theta,velocity]),recovered_parameters=np.r_[recovered_eta,solution.x],
                EOS=a,truth_EOS=truth,residual=residual,jacobian=jac)
            save('progress.json',dict(classification='Counterexample candidate',rows=rows));print('SUBCELL INVERSE',row,flush=True)
    save('result.json',dict(classification='Counterexample candidate',completed=True,rows=rows,all_passed=all(r['passed'] for r in rows),
        EOS_evaluations=evaluations,nonuniform_reference_retained=True,metric_fixed=True,physical_EOS_certified=False,
        arbitrary_cell_moment_closure=False,continuous_root_error_certified=False,full_GR_evolution=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify()


def verify():
    for rel,digest in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert g.c.sha(reference.OUT/'manifest.json')==json.loads((OUT/'startup-binding.json').read_text())['reference_manifest_sha256']
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'result.json').read_text())['completed']
    for row in json.loads((OUT/'Newton-traces.json').read_text())['rows']:
        assert all(b['merit']<a['merit'] for a,b in zip(row['history'],row['history'][1:]))
    print('PASS nonuniform subcell inverse bindings; prescribed reconstruction family only',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
