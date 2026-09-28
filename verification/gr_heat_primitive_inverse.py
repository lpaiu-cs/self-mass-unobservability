"""Pointwise nonlinear primitive recovery with a comoving radial heat flux."""
import json,sys
from fractions import Fraction as F
from functools import lru_cache
import numpy as np
import sympy as s
from scipy.optimize import root
import gr_heat_entropy_closure as heat

g=heat.g;OUT=g.OUT/'gr-heat-primitive-inverse'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();heat.verify()
    paths=[g.ROOT/'verification/gr_heat_primitive_inverse.py',g.OUT/'initial-state-17-4.npz',
        g.OUT/'gr-microphysics/auxiliaries.npz',heat.OUT/'initial-rates.npz',
        heat.old.char.OUT/'coefficients.npz',heat.OUT/'manifest.json']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='a2d50b8',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        cells=[0,1,2,64,1024,2048,2972,4096,5734],
        probes=[[-1e-4,2e-4,-1e-4],[0.,0.,0.],[1e-4,-2e-4,1e-4]],
        probe_order='log density shift, log temperature shift, v/c; Q retains the saved comoving heat value.',
        frozen_velocity_interval=[-.5,.5],root_logT_tolerance=1e-9,root_velocity_tolerance=1e-12,
        root_log_density_tolerance=1e-10,scaled_residual_tolerance=1e-10,
        root_guess='lnT_true+0.003, v=0, density inferred from fixed normal baryon D=rho*W.',
        scope='Nonlinear point-state Lorentz/EOS inverse with explicit rest-energy subtraction. Manufactured/perturbed native states, not a GR trajectory, observational fit or closure for unresolved cell moments. Frozen-coefficient determinant bounds do not enclose native EOS changes over a neighbourhood.'))
    symbolic()


def symbolic():
    rho,w,b,r,d,e,j,v=s.symbols('rho w b r d e j v',real=True);W=1/s.sqrt(1-v*v)
    matrix=s.Matrix([[rho*W,0,rho*v*W**3],
        [w*W**2*(1-e+r*v*v),w*W**2*(b+d*v*v),2*w*W**4*(v+j*(1+v*v))],
        [w*v*W**2*(1-e+r),w*v*W**2*(b+d),w*W**4*(1+v*v+4*v*j)]])
    polynomial=b-v*v*(b*r+d*e)+2*j*v*(b-d)
    assert s.simplify(matrix.det()-rho*w*w*W**5*polynomial)==0
    assert s.simplify(matrix.det().subs(v,0)-rho*w*w*b)==0
    eps,P,Q=s.symbols('eps P Q',real=True)
    E=(eps+P*v*v+2*v*Q)/(1-v*v);J=((eps+P)*v+Q*(1+v*v))/(1-v*v)
    assert s.simplify(J-Q-v*(E+P))==0
    assert s.simplify(eps-E+v*(J+Q))==0
    assert s.simplify(W-1-v*v*W*W/(W+1))==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        primitive_relations='At fixed comoving Q and composition, rho=D/W, J-Q=v*(E+P), epsilon=E-v*(J+Q). A complete heat evolution must also supply Q; these three conserved values alone do not determine an extra unconstrained heat variable.',
        Jacobian='For (ln rho,ln T,v) -> (D,E,J), fixed Q, let w=epsilon+P, b=rho*cv*T/w, r=P_lnrho/w, d=P_lnT/w, e=(P-rho*u_lnrho)/w, j=Q/w. The determinant is rho*w^2*W^5*[b-v^2*(b*r+d*e)+2*j*v*(b-d)]. At v=0 it is rho^2*cv*T*w>0 under positive rho,cv,T,w.',
        conditional_local_inverse='A differentiable EOS and nonzero determinant imply a local pointwise inverse by the inverse function theorem. For frozen g=r+d*e/b>=0 and |v|<=V<1, a sufficient bound is b*(1-g*V^2)-2*abs(j)*V*abs(b-d)>0. This is not a proof of global uniqueness or native-neighbourhood invertibility.',
        stable_energy='Use E-D*C = D*C*(W-1)+rho*u*W^2+P*v^2*W^2+2*v*Q*W^2, C=C_X*c^2. Evaluate W-1=v^2*W^2/(W+1) directly. This avoids subtracting two rest-energy-sized numbers; it does not create information absent from rounded conserved input.',
        cell_caution='A pointwise inverse does not identify an unresolved nonuniform cell with a uniform EOS state while preserving every subcell moment or original entropy. This test supplies a necessary nonlinear local component, not that reconstruction.'))


def forward(rho,a,v,Q,C):
    ld=np.longdouble;rho,u,P,v,Q,C=map(ld,[rho,a[2],a[1],v,Q,C]);W=1/np.sqrt(1-v*v)
    D=rho*W;energy=D*C*(v*v*W*W/(W+1))+rho*u*W*W+P*v*v*W*W+2*v*Q*W*W
    J=((rho*(C+u)+P)*v+Q*(1+v*v))*W*W
    return np.array([D,energy,J],dtype=ld)


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    fields=np.load(heat.old.char.OUT/'coefficients.npz')['coefficients'];heats=dict(np.load(heat.OUT/'initial-rates.npz'))
    lower=[]
    for values,j in zip(fields,heats['j']):
        b,r,d,e,h=[F(float(x)) for x in values];j=F(float(j));V=F(1,2);sound=r+d*e/b
        assert b>0 and sound>=0
        lower.append((b*(1-sound*V*V)-2*abs(j)*V*abs(b-d))/b)
    frozen=dict(classification='Proven',cells=len(lower),all_positive=all(x>0 for x in lower),
        minimum_scaled_lower_bound=str(min(lower)),display=float(min(lower)),
        scope='Exact rational polynomial inequalities for the frozen stored coefficients and |v|<=1/2 only.')
    save('frozen-determinants.json',frozen)
    control=forward(1.,[1.,0.,1e-25],0.,0.,1.)
    assert control[1]==np.longdouble(1e-25) and np.longdouble(1)+np.longdouble(1e-25)-1==0
    state=dict(np.load(g.OUT/'initial-state-17-4.npz'));eos=g.EOS();rows=[];evaluations=0;c=g.c.gr.C*100
    for i in plan['cells']:
        C=state['CX'][i]*c*c;Q=heats['Q'][i]
        for case,(dr,dT,velocity) in enumerate(plan['probes']):
            lr=float(state['lnd'][i]+dr);lt=float(state['lnT'][i]+dT);truth=eos(2,lr,lt,state['X'][i]);evaluations+=1
            rho=np.exp(lr);target=forward(rho,truth,velocity,Q,C);D=target[0]
            scale=np.array([rho*truth[10],rho*(C+truth[2])+truth[1]],dtype=np.longdouble)
            @lru_cache(maxsize=1)
            def evaluate(t,v):
                nonlocal evaluations
                assert abs(v)<.5,(i,case,'primitive speed outside bracket',v)
                W=1/np.sqrt(1-np.longdouble(v)**2);rh=D/W
                a=eos(2,float(np.log(rh)),float(t),state['X'][i]);evaluations+=1
                value=forward(rh,a,v,Q,C)[1:];residual=np.asarray((value-target[1:])/scale,float)
                P=a[1];w=rh*(C+a[2])+P;er=rh*(C+a[2]+a[9]);et=rh*a[10]
                pr=P*a[5];pt=P*a[6]
                jac=np.array([[W*W*(et+pt*v*v),W**4*(2*(v*w+Q*(1+v*v))-v*(er+pr*v*v))],
                    [v*W*W*(et+pt),W**4*((1+v*v)*w+4*v*Q-v*v*(er+pr))]],dtype=np.longdouble)
                return residual,np.asarray(jac/scale[:,None],float),float(np.log(rh)),a
            solution=root(lambda x:evaluate(*map(float,x))[0],[lt+.003,0.],
                jac=lambda x:evaluate(*map(float,x))[1],method='hybr',options={'xtol':1e-10})
            residual,jac,recovered_lr,a=evaluate(*map(float,solution.x))
            row=dict(cell=i,case=case,solver_success=bool(solution.success),message=solution.message,
                logT_error=float(abs(solution.x[0]-lt)),velocity_error=float(abs(solution.x[1]-velocity)),
                log_density_error=abs(recovered_lr-lr),scaled_residual=float(abs(residual).max()))
            row['passed']=bool(row['logT_error']<=plan['root_logT_tolerance'] and
                row['velocity_error']<=plan['root_velocity_tolerance'] and row['log_density_error']<=plan['root_log_density_tolerance'] and
                row['scaled_residual']<=plan['scaled_residual_tolerance'])
            rows.append(row);np.savez_compressed(OUT/f'cell-{i}-case-{case}.npz',conserved=target,
                true_primitive=np.array([lr,lt,velocity]),recovered_primitive=np.r_[recovered_lr,solution.x],residual=residual,jacobian=jac)
            save('progress.json',dict(classification='Counterexample candidate',rows=rows))
            print('NONLINEAR HEAT PRIMITIVE',row,flush=True)
    save('result.json',dict(classification='Counterexample candidate',completed=True,rows=rows,all_passed=all(r['passed'] for r in rows),
        EOS_evaluations=evaluations,rest_subtraction_positive_control=True,physical_EOS_certified=False,
        continuous_root_error_certified=False,conserved_cell_average_closure=False,full_GR_evolution=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify()


def verify():
    for rel,digest in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'symbolic.json').read_text())['passed']
    assert json.loads((OUT/'result.json').read_text())['completed']
    print('PASS nonlinear point-primitive result bindings; full inverse/GR error remains open',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
