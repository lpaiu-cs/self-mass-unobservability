"""Native-EOS cell inversion with the radial metric recomputed from energy.

This is a conservative reconstruction component, not a stellar time path.
The represented coordinate-volume Jacobian and heat profile are prescribed;
the metric and its complete energy chain rule are not frozen.
"""
from concurrent.futures import ProcessPoolExecutor,as_completed
from functools import lru_cache
from pathlib import Path
import json,re,subprocess,sys,traceback
import numpy as np
import sympy as sp
import gr_subcell_precision_resume as coverage
import gr_heat_primitive_newton as point

g=coverage.g;ROOT=g.ROOT;OUT=g.OUT/'gr-metric-coupled-subcell';sha=g.c.sha;ld=np.longdouble


def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2)+'\n')


def runtime():
    library=g.d.CACHE/'direct_ion_integral_full_bridge.so'
    paths=[library]
    for line in subprocess.check_output(['ldd',str(library)],text=True).splitlines():
        match=re.search(r'(?:=>\s+)?(/\S+)\s+\(',line)
        if match:paths.append(Path(match.group(1)).resolve())
    return {str(p):sha(p) for p in paths}


def symbolic():
    eps,P,Q,v,r,m,dm,a0,z=sp.symbols('eps P Q v r m dm a0 z',real=True)
    E=(eps+P*v*v+2*Q*v)/(1-v*v);S=((eps+P)*v+Q*(1+v*v))/(1-v*v)
    assert sp.simplify(sp.diff(E,v)-2*(v*(eps+P)+Q*(1+v*v))/(1-v*v)**2)==0
    assert sp.simplify(sp.diff(S,v)-((1+v*v)*(eps+P)+4*Q*v)/(1-v*v)**2)==0
    a=(1-2*m/r)**-sp.Rational(1,2)
    assert sp.simplify(sp.diff(a,m)-a**3/r)==0
    assert sp.simplify((1/sp.sqrt(1-z)-1)-z/(sp.sqrt(1-z)*(1+sp.sqrt(1-z))))==0
    E0,de=sp.symbols('E0 de',real=True)
    assert sp.simplify((E-E0).subs(eps,E0+de)-(de+(P+E0)*v*v+2*Q*v)/(1-v*v))==0
    assert sp.simplify(S-Q-((eps+P)*v+2*Q*v*v)/(1-v*v))==0
    x=sp.symbols('x');checks=0
    for n in range(1,17):
        assert sp.simplify(sp.diff(sp.legendre(n+1,x)-sp.legendre(n-1,x),x)-(2*n+1)*sp.legendre(n,x))==0;checks+=1
    save('symbolic.json',dict(classification='Proven',passed=True,Legendre_antiderivatives=checks,
        geometry='At fixed radius nodes and a prescribed coordinate-volume Jacobian V(s), project f(s)=E(s)V(s) on the Gauss interpolation polynomial. m_j=m_inner+(G/c^4)*sum_k A_jk E_k V_k; a_j=(1-2m_j/r_j)^(-1/2). A_jk is the integrated cardinal polynomial from the inner face to node j.',
        chain_rule='For each primitive parameter x, m_x=(G/c^4)*A*(V*E_x), a_x=a^3*m_x/r. Therefore B_x=sum wV*(a*D_x+D*a_x), M_x=sum wV*E_x, Pi_x=sum wV*(a*S_x+S*a_x). No metric derivative is discarded.',
        increments='Store the baseline and changes separately. Delta a=a0*z/[sqrt(1-z)*(1+sqrt(1-z))], z=2*a0^2*Delta m/r. Delta E=[Delta epsilon+(P+epsilon0)*v^2+2Qv]/(1-v^2). Delta S=[(epsilon+P)*v+2Qv^2]/(1-v^2).',
        centre='In the central cell use v(r)=v_amplitude*r/r_outer and a luminosity profile linear in enclosed baryon fraction, so regular finite-density data have v=O(r),Q=O(r). Other cells use a common velocity amplitude.',
        boundary='These are exact identities for a declared differentiable EOS and the chosen collocation representation. Its native EOS, interpolation, centre seed, continuum Hamiltonian error and physical transport are not certified. Q is prescribed; its own evolution is still required.'))


def integration_matrix(number):
    x,w=np.polynomial.legendre.leggauss(number);s=(x+1)/2;w=w/2
    P=np.polynomial.legendre.legvander(x,number)
    A=np.tile(s[:,None],(1,number))
    for n in range(1,number):A+=.5*(P[:,n+1]-P[:,n-1])[:,None]*P[:,n][None,:]
    A*=w[None,:]
    errors=[float(np.max(abs(A@(s**k)-s**(k+1)/(k+1)))) for k in range(number)]
    assert max(errors)<1e-13
    return s.astype(ld),w.astype(ld),A.astype(ld),errors


def prepare():
    assert not OUT.exists();coverage.verify();assert json.loads((coverage.OUT/'result.json').read_text())['all_passed']
    original=json.loads((point.OUT/'plan.json').read_text());sources=json.loads((coverage.OUT/'result.json').read_text())['sources']
    cell_sources={};files=[Path(__file__),ROOT/'verification/direct_eos_gr.py',ROOT/'verification/direct_ion_eos.py',
        Path(g.c.__file__),Path(g.c.gr.__file__),coverage.OUT/'manifest.json',point.OUT/'plan.json',g.OUT/'initial-state-17-4.npz',
        g.OUT/'gr-increment-structure/path-4.npz',g.OUT/'gr-transport/diagnostics.npz',
        g.OUT/'gr-microphysics/auxiliaries.npz',
        g.d.OUT/'model-data.json',g.d.OUT/'full-integral-build.json']
    for source in sources:
        path=ROOT/source['path'];files.append(path)
        for n in [8,16]:files.append(path.with_name(path.stem+f'-nodes-{n}.npz'))
        for i in source['cells']:cell_sources[str(i)]=path.relative_to(ROOT).as_posix()
    assert sorted(map(int,cell_sources))==list(range(5735))
    OUT.mkdir();plan=dict(original,classification='Counterexample candidate',checkpoint='b6048be8',
        bindings={p.relative_to(ROOT).as_posix():sha(p) for p in files},runtime_sha256=runtime(),
        cells=list(range(5735)),pilot_cells=[0,1,2,2688,2972,5734],nodes=[8,16],workers=2,block_size=16,
        cell_sources=cell_sources,inner_log_mass_shifts=[-1e-4,0.,1e-4],logrho_step_cap=.1,
        finite_Jacobian_steps=[1e-5,5e-6],finite_Jacobian_relative_gate=1e-4,
        reference_metric_relative_gate=1e-8,finite_quadrature_gate=1e-5,
        root_guess='Fixed [0,0.003,0] in log-density shift, log-temperature shift and velocity amplitude, without truth-derived scales.',
        representation='Keep every original positive coordinate-volume weight and nonuniform native reference node. Reorient nodes inner-to-outer. Integrate the Legendre projection of E*dV/ds to obtain each metric value from a supplied inner mass. The inner mass has a separately prescribed perturbation for each manufactured probe. Do not rescale quadrature weights to manufacture mass conservation.',
        heat='One shared original luminosity per material face, affine interpolation in enclosed baryon fraction, original isentropic enthalpy lapse reconstruction. Proper Q=L/(4*pi*r^2*N^2*c) is prescribed during the inverse. In the central cell velocity is proportional to radius. Neither a constant central heat density nor a constant nonzero central velocity is used.',
        conserved='Invert the changes of B=sum aD*dV, M=sum E*dV, Pi=sum aS*dV. Recompute a from the full nodal energy, including rest energy and moving heat. Use a fixed invertible energy-row subtraction C*kappa0*B for conditioning; kappa0=sum rho0*dV/B0. This does not assume the current coordinate/proper baryon ratio is constant.',
        gates='Retain original primitive and scaled-residual gates. Require metric agreement, full analytic-Jacobian finite controls and finite 8/16 moment controls before a full inventory run. At least one nonzero preflight target must fail the original residual gate when the metric is deliberately frozen. Preserve failed cases and do not skip them.',
        scope='Nonlinear metric-coupled conservative reconstruction tests on the actual native-EOS profile. Targets are manufactured, not an evolved star. No arbitrary global closure, continuous native/metric error, physical heat calibration, atmosphere, time trajectory or observation inference.')
    save('plan.json',plan);symbolic()
    save('quadrature-controls.json',dict(classification='Counterexample candidate',passed=True,
        rows=[dict(nodes=n,polynomial_integral_errors=integration_matrix(n)[3]) for n in plan['nodes']]))


def bindings():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert sha(ROOT/rel)==digest,rel
    for path,digest in plan['runtime_sha256'].items():assert sha(Path(path))==digest,path
    return plan


def initialize():
    global PLAN,EOS,STATE,GRID,LUMINOSITY,AUX
    PLAN=bindings();EOS=g.EOS();STATE=dict(np.load(g.OUT/'initial-state-17-4.npz'))
    GRID=dict(np.load(g.OUT/'gr-increment-structure/path-4.npz'))
    AUX=np.load(g.OUT/'gr-microphysics/auxiliaries.npz')['eos']
    LUMINOSITY=np.r_[0.,np.load(g.OUT/'gr-transport/diagnostics.npz')['interior_Linf'],0.].astype(ld)


@lru_cache(maxsize=2)
def node_file(path):return dict(np.load(path))


class Cell:
    def __init__(self,index,number):
        path=ROOT/PLAN['cell_sources'][str(index)];data=node_file(str(path.with_name(path.stem+f'-nodes-{number}.npz')))
        j=list(data['cells']).index(index);r=data['radius_cm'][j].astype(ld)
        order=np.arange(number) if np.all(np.diff(r)>0) else np.arange(number-1,-1,-1)
        self.r=r[order];assert np.all(np.diff(self.r)>0)
        self.weights=data['coordinate_weights_cm3'][j,order].astype(ld);assert np.all(self.weights>0)
        self.rho0=data['eos'][j,order,0].astype(ld);self.logrho0=np.log(self.rho0).astype(float)
        self.logT0=data['lnT'][j,order].astype(float);self.X=STATE['X'][index]
        self.c=ld(g.c.gr.C)*100;self.gravity=ld(g.c.gr.G)*1000/self.c**4
        self.C=ld(data['C_X'][j])*self.c*self.c
        self.base_eos=self.native(0.,0.);self.eps0=self.rho0*(self.C+self.base_eos[:,2])
        s,w,A,_=integration_matrix(number);self.integral=A*(self.weights/w)[None,:]
        self.inner_mass=ld(GRID['mass_geom_m'][index+1])*100
        self.m0=self.inner_mass+self.gravity*(self.integral@self.eps0)
        self.a0=1/np.sqrt(1-2*self.m0/self.r);assert np.all(np.isfinite(self.a0))
        shell=self.gravity*np.sum(self.eps0*self.weights)
        assert np.all(self.m0>=self.inner_mass) and np.all(self.m0<=self.inner_mass+shell)
        self.metric_reference_difference=float(np.max(abs(self.a0/data['metric_a'][j,order]-1)))
        assert self.metric_reference_difference<=PLAN['reference_metric_relative_gate'],(index,number,self.metric_reference_difference)
        outside=index<len(GRID['outer'])-1;branch=GRID['outer'] if outside else GRID['inner']
        k=index if outside else len(STATE['dm'])-1-index
        low=ld(0) if k==0 else ld(branch[k,0]);high=ld(branch[k+1,0]);knots,_=np.polynomial.legendre.leggauss(number)
        if outside:fraction=(1-knots.astype(ld))/2
        else:
            left=np.cbrt(low);right=np.cbrt(high);q=((left+right)/2+(right-left)*knots.astype(ld)/2)**3
            fraction=(q-low)/(high-low)
        fraction=fraction[order];assert np.all((fraction>0)&(fraction<1))
        eos_mid=AUX[index]
        rho_mid=np.exp(ld(STATE['lnd'][index]));Hmid=ld(eos_mid[2])+ld(eos_mid[1])/rho_mid
        H=self.base_eos[:,2]+self.base_eos[:,1]/self.rho0
        lapse=np.exp(ld(STATE['nu'][index])+np.log1p((Hmid-H)/(self.C+H)))
        luminosity=LUMINOSITY[index+1]+fraction*(LUMINOSITY[index]-LUMINOSITY[index+1])
        self.Q=luminosity/(4*np.pi*self.r*self.r*lapse*lapse*self.c)
        self.velocity_shape=self.r/(ld(GRID['radius_m'][index])*100) if index==5734 else np.ones(number,dtype=ld)
        self.B0=np.sum(self.rho0*self.a0*self.weights)
        self.kappa0=np.sum(self.rho0*self.weights)/self.B0
        self.capacity=np.sum(self.rho0*self.base_eos[:,10]*self.weights)
        self.enthalpy=np.sum((self.eps0+self.base_eos[:,1])*self.a0*self.weights)
        self.scale=np.array([self.B0,self.capacity,self.enthalpy],dtype=ld);assert np.all(self.scale>0)

    def native(self,eta,theta):
        values=np.array([EOS(2,float(r+eta),float(t+theta),self.X) for r,t in zip(self.logrho0,self.logT0)],dtype=ld)
        assert np.all(values[:,1]>0) and np.all(values[:,10]>0);return values

    def transform(self,value):
        result=value.copy();result[1]-=self.C*self.kappa0*result[0]
        return result/self.scale if result.ndim==1 else result/self.scale[:,None]

    def evaluate(self,x,inner_shift,frozen=False):
        eta,theta,amplitude=map(ld,x);v=amplitude*self.velocity_shape;W=1/np.sqrt(1-v*v)
        aEOS=self.native(float(eta),float(theta));rho=self.rho0*np.exp(eta);P=aEOS[:,1];u=aEOS[:,2]
        eps=rho*(self.C+u);enthalpy=eps+P;D=rho*W
        deps=self.rho0*((self.C+self.base_eos[:,2])*np.expm1(eta)+np.exp(eta)*(u-self.base_eos[:,2]))
        dE=(deps+(P+self.eps0)*v*v+2*self.Q*v)*W*W
        dS=(enthalpy*v+2*self.Q*v*v)*W*W;S=self.Q+dS
        dm=self.inner_mass*np.expm1(ld(inner_shift))+self.gravity*(self.integral@dE)
        z=2*self.a0*self.a0*dm/self.r;assert np.all(1-z>0)
        da=self.a0*z/(np.sqrt(1-z)*(1+np.sqrt(1-z)))
        if frozen:da=np.zeros_like(da)
        a=self.a0+da
        dD=self.rho0*(W*np.expm1(eta)+v*v/(np.sqrt(1-v*v)*(1+np.sqrt(1-v*v))))
        changes=np.array([np.sum((a*dD+self.rho0*da)*self.weights),np.sum(dE*self.weights),
            np.sum((a*dS+self.Q*da)*self.weights)],dtype=ld)
        er=rho*(self.C+u+aEOS[:,9]);et=rho*aEOS[:,10];pr=P*aEOS[:,5];pt=P*aEOS[:,6]
        Ex=np.array([W*W*(er+pr*v*v),W*W*(et+pt*v*v),
            2*W**4*(v*enthalpy+self.Q*(1+v*v))*self.velocity_shape]).T
        Sx=np.array([W*W*v*(er+pr),W*W*v*(et+pt),
            W**4*((1+v*v)*enthalpy+4*self.Q*v)*self.velocity_shape]).T
        Dx=np.array([D,np.zeros_like(D),rho*W**3*v*self.velocity_shape]).T
        ax=a[:,None]**3/self.r[:,None]*self.gravity*(self.integral@Ex)
        if frozen:ax=np.zeros_like(ax)
        J=np.array([np.sum(self.weights[:,None]*(a[:,None]*Dx+D[:,None]*ax),axis=0),
            np.sum(self.weights[:,None]*Ex,axis=0),np.sum(self.weights[:,None]*(a[:,None]*Sx+S[:,None]*ax),axis=0)])
        details=dict(maximum_relative_metric_change=float(np.max(abs(da/self.a0))),
            metric_chain_B_eta_relative=float(np.sum(self.weights*D*ax[:,0])/J[0,0]),
            minimum_metric_factor=float(np.min(a)))
        return self.transform(changes),self.transform(J),changes,details


def newton(cell,target,inner_shift):
    x=np.array([0.,.003,0.]);history=[];success=False
    for iteration in range(PLAN['max_iterations']):
        value,J,_,_=cell.evaluate(x,inner_shift);residual=np.asarray(value-target,float);J=np.asarray(J,float)
        merit=float(np.max(abs(residual)));row=dict(iteration=iteration,x=x.tolist(),merit=merit,
            condition=float(np.linalg.cond(J)),trials=[]);history.append(row)
        if merit<=PLAN['Newton_scaled_residual']:success=True;break
        step=np.linalg.solve(J,-residual);step*=min(1.,.1/max(abs(step[0]),1e-300),.1/max(abs(step[1]),1e-300),.01/max(abs(step[2]),1e-300))
        for backtrack in range(PLAN['max_backtracks']):
            proposed=x+step*2.**(-backtrack)
            if abs(proposed[2])>=.5:continue
            new_merit=float(np.max(abs(cell.evaluate(proposed,inner_shift)[0]-target)))
            row['trials'].append(dict(backtrack=backtrack,merit=new_merit))
            if new_merit<merit:x=proposed;row['accepted_backtrack']=backtrack;break
        else:break
    assert all(b['merit']<a['merit'] for a,b in zip(history,history[1:]))
    return x,success,history


def one_cell(index,controls=False):
    rows=[];targets={};matrices=[]
    for number in PLAN['nodes']:
        cell=Cell(index,number);targets[number]=[]
        for case,(truth,inner_shift) in enumerate(zip(PLAN['probes'],PLAN['inner_log_mass_shifts'])):
            target,J,moments,detail=cell.evaluate(truth,inner_shift);targets[number].append(target)
            recovered,success,history=newton(cell,target,inner_shift)
            residual,_,_,_=cell.evaluate(recovered,inner_shift);error=abs(recovered-np.array(truth))
            frozen=cell.evaluate(truth,inner_shift,True)[0]
            row=dict(cell=index,nodes=number,case=case,solver_success=success,errors=error.tolist(),
                scaled_residual=float(np.max(abs(residual-target))),maximum_frozen_metric_scaled_defect=float(np.max(abs(frozen-target))),
                reference_metric_relative_difference=cell.metric_reference_difference,geometry=detail,history=history,
                target_changes_exact=[[str(n),str(d)] for n,d in (value.as_integer_ratio() for value in moments)])
            row['passed']=bool(success and error[0]<=PLAN['root_log_density_tolerance'] and error[1]<=PLAN['root_logT_tolerance']
                and error[2]<=PLAN['root_velocity_tolerance'] and row['scaled_residual']<=PLAN['scaled_residual_tolerance'])
            rows.append(row)
            if controls and number==16 and case==2:
                for h in PLAN['finite_Jacobian_steps']:
                    fd=[]
                    for k in range(3):
                        delta=np.zeros(3);delta[k]=h
                        fd.append((cell.evaluate(np.array(truth)+delta,inner_shift)[0]-cell.evaluate(np.array(truth)-delta,inner_shift)[0])/(2*ld(h)))
                    fd=np.array(fd).T;den=np.maximum(np.max(abs(J),axis=1),ld('1e-30'))
                    score=float(np.max(abs(fd-J)/den[:,None]));matrices.append(dict(cell=index,step=h,score=score,passed=score<=PLAN['finite_Jacobian_relative_gate']))
        assert len(targets[number])==3
    gaps=[]
    for a,b in zip(targets[8],targets[16]):gaps.append(float(np.max(abs(a-b)/np.maximum(abs(b),ld('1e-10')))))
    return dict(classification='Counterexample candidate',cell=index,rows=rows,Jacobian_controls=matrices,
        finite_8_16_target_differences=gaps,finite_quadrature_passed=max(gaps)<=PLAN['finite_quadrature_gate'],
        all_inverse_passed=all(r['passed'] for r in rows),all_Jacobian_controls_passed=all(r['passed'] for r in matrices),
        metric_fixed=False,physical_EOS_certified=False,full_GR_evolution=False)


def preflight():
    initialize();assert not (OUT/'preflight.json').exists();rows=[];failures=[]
    for i in PLAN['pilot_cells']:
        try:rows.append(one_cell(i,True));print('METRIC CELL PREFLIGHT',i,rows[-1]['all_inverse_passed'],rows[-1]['finite_quadrature_passed'],flush=True)
        except Exception:failures.append(dict(cell=i,traceback=traceback.format_exc()));print('METRIC CELL PREFLIGHT ERROR',i,failures[-1]['traceback'],flush=True)
        save('preflight-progress.json',dict(rows=rows,failures=failures))
    negative=any(r['maximum_frozen_metric_scaled_defect']>PLAN['scaled_residual_tolerance'] for cell in rows for r in cell['rows'])
    passed=not failures and negative and all(r['all_inverse_passed'] and r['all_Jacobian_controls_passed'] and r['finite_quadrature_passed'] for r in rows)
    save('preflight.json',dict(classification='Counterexample candidate',passed=passed,rows=rows,failures=failures,frozen_metric_negative_control_passed=negative))
    files=[OUT/'plan.json',OUT/'symbolic.json',OUT/'quadrature-controls.json',OUT/'preflight.json']
    save('preflight-manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in files}))
    assert passed,'Original preflight failure retained; full-grid start is gated';verify_preflight()


def verify_preflight():
    plan=bindings()
    for rel,digest in json.loads((OUT/'preflight-manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    r=json.loads((OUT/'preflight.json').read_text());assert r['passed'] and [a['cell'] for a in r['rows']]==plan['pilot_cells']
    print('PASS native-EOS, nonuniform, metric-coupled primitive preflight',flush=True)


def block(cells):
    folder=OUT/f'block-{cells[0]:04d}';folder.mkdir();rows=[];failures=[]
    for i in cells:
        try:rows.append(one_cell(i))
        except Exception:failures.append(dict(cell=i,traceback=traceback.format_exc()))
    result=dict(classification='Counterexample candidate',cells=cells,rows=rows,failures=failures,
        passed=not failures and all(r['all_inverse_passed'] and r['finite_quadrature_passed'] for r in rows))
    (folder/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    return dict(start=cells[0],cells=len(cells),passed=result['passed'],sha256=sha(folder/'result.json'))


def run():
    verify_preflight();plan=bindings();assert not any(OUT.glob('block-*'));records=[]
    jobs=[list(range(i,min(i+plan['block_size'],len(plan['cells'])))) for i in range(0,len(plan['cells']),plan['block_size'])]
    with ProcessPoolExecutor(max_workers=plan['workers'],initializer=initialize) as pool:
        futures=[pool.submit(block,cells) for cells in jobs]
        for f in as_completed(futures):
            records.append(f.result());save('progress.json',dict(classification='Counterexample candidate',blocks=records))
            print('METRIC CELL BLOCK',records[-1],flush=True)
    save('result.json',dict(classification='Counterexample candidate',completed=True,cells=sum(r['cells'] for r in records),
        passed=all(r['passed'] for r in records),blocks=records,physical_EOS_certified=False,full_GR_evolution=False))
    save('manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in OUT.rglob('*') if p.is_file()}));verify()


def verify():
    verify_preflight()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());assert r['completed'] and r['cells']==5735
    print('PASS complete metric-coupled inverse inventory bindings; inspect passed flag and retained failures',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
