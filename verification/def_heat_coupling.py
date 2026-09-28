"""Regular two-heat local matter/scalar coupling. No spherical GR evolution.

Proven: homogeneous metric-work/entropy identities within the declared closure.
Counterexample candidate: native EOS and opacity source-step controls.
"""
from pathlib import Path
import argparse
import hashlib
import json
import time
import numpy as np
import sympy as sp
from scipy.optimize import root
import def_native_coupling as prior
import gr_two_carrier_evolution as two

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/direct-eos-gr33/def-heat-coupling'
INITIAL=ROOT/'outputs/direct-eos-gr33/gr-molecular-shell-initial/initial.npz'
ld=np.longdouble


def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def symbolic():
    v,T,dv,dt,ell=sp.symbols('v T dv dt ell',real=True)
    q1,q2,dq1,dq2,b1,b2,db1,db2,tau1,tau2,W=sp.symbols('q1 q2 dq1 dq2 b1 b2 db1 db2 tau1 tau2 W',positive=True)
    q=q1+q2;dq=dq1+dq2;W2=1/(1-v*v)
    # Gibbs + energy/momentum and baryon conservation; scalar work cancels.
    rho_ds=-(v*dq+2*q*W2*dv+2*q*v*ell)/T
    # d[a W (rho*s+v*Q/T-sum(beta*q^2)/2)]/(a W).
    entropy=rho_ds+(v*dq+q*dv-v*q*dt+v*q*(ell+W2*v*dv))/T
    entropy-=b1*q1*dq1+b2*q2*dq2+(b1*q1*q1*(db1+ell+W2*v*dv)+b2*q2*q2*(db2+ell+W2*v*dv))/2
    laws=[dq1+(W2*dv+v*(ell+dt))/(b1*T)+q1*(db1+ell+W2*v*dv)/2,
          dq2+(W2*dv+v*(ell+dt))/(b2*T)+q2*(db2+ell+W2*v*dv)/2]
    assert sp.simplify(entropy+b1*q1*laws[0]+b2*q2*laws[1])==0
    # With reversible law=0 entropy is constant; law=-q/(tau*W) gives production.
    production=b1*q1*q1/(tau1*W)+b2*q2*q2/(tau2*W)
    assert sp.simplify((-b1*q1*(-q1/(tau1*W))-b2*q2*(-q2/(tau2*W)))-production)==0
    ct,w,C1,C2,qt,br1,br2,bt1,bt2=sp.symbols('ct w C1 C2 qt br1 br2 bt1 bt2')
    M=sp.Matrix([[ct,2*(q1+q2),0,0],[0,w,1,1],
        [q1*bt1/2,C1,1,0],[q2*bt2/2,C2,0,1]])
    determinant=sp.factor(M.det())
    assert sp.simplify(determinant-(ct*(w-C1-C2)+(q1+q2)*(q1*bt1+q2*bt2)))==0
    a,ar0,ar1,g0,g1,dp,phi0,phi1,p0,p1,h,I=sp.symbols('a ar0 ar1 g0 g1 dp phi0 phi1 p0 p1 h I')
    dphi=h*(p0+p1)/2
    scalar_increment=I*(p1-p0)*(p1+p0)/2+I*(phi1-phi0)*(phi1+phi0)/2
    assert sp.simplify(scalar_increment.subs(p1-p0,h*(-(phi1+phi0)/2+(g0+g1)/(2*I))).subs(phi1-phi0,dphi)-(g0+g1)*dphi/2)==0
    return dict(classification='Proven',passed=True,
        entropy='Sigma=a*W*(rho*s+v*Q/T-sum(beta_j*q_j^2)/2); Sigma_dot=a*sum(beta_j*q_j^2/tau_j)>=0',
        rest_matrix_determinant=str(determinant),
        scope='Homogeneous radial metric deformation, fixed composition, Gibbs relation and full diagonal quadratic entropy closure. A dynamical heat/temperature lift, not a fixed-entropy canonical heat coordinate or a full GR theorem.')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    files=[Path(__file__),Path(prior.__file__),Path(two.__file__),INITIAL,
           ROOT/'outputs/direct-eos-gr33/gr-regular-scalar-milestone-manifest.json',Path(prior.initial.molecular.model.LIB)]
    save('plan.json',dict(classification='Counterexample candidate',bindings={str(p.resolve()):digest(p) for p in files},
        cell=2762,steps=[8,16,32],duration_scaled=6.283185307179586,
        initial='Saved native shell rho,T,26-species composition. v=0, scalar p=0; phi=.01*sqrt(thermal derivative/enthalpy), each heat loading q_j=.01*sqrt(C_j*thermal derivative). Each channel therefore starts with the same small quadratic entropy deficit. Loading is a declared stress control, not the stellar heat flux.',
        metric='a(t)=exp(.001*sin(t)); external homogeneous radial strain, not an Einstein solution. Unit time is initial radiation relaxation time.',
        method='Conserve B=a*rho*W and J=a^2*S, solve finite aE energy work, both full midpoint heat laws and reciprocal scalar midpoint simultaneously. Scalar inertia=16*initial enthalpy. No inverse scalar momentum or inverse field increment.',
        gates=dict(stage_residual=2e-8,work_defect_over_initial_thermal=1e-8,
            refinement_order=1.7,entropy_budget_relative=0.15,entropy_increase_positive=True,
            maximum_native_calls=6000),
        budget=dict(cpu=1,blas_threads=1,pilot_steps=8,pilot_timeout_seconds=90,
            full_timeout_seconds=240,launch_only_if_projected_seconds_below=180,maximum_full_runs=1,
            no_automatic_grid_or_duration_extension=True),
        native_EOS=True,native_opacity=True,full_GR=False,physical_transport_calibrated=False,observational_closure=False))
    save('symbolic.json',symbolic());print('PREPARED heat coupling',flush=True)


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for name,sha in p['bindings'].items():assert digest(name)==sha,name
    return p


class Cell:
    def __init__(self):
        self.plan=bindings();data=np.load(INITIAL);i=self.plan['cell']
        self.base=data['base'][i];self.X=self.base[5:];self.provider=prior.Matter()
        self.opacity=two.radiative.tables.Opacity();self.cache={}
        zero=self.provider.state(self.base[0],self.base[1],0,self.X)
        self.thermal=zero['denergy_dlogT'];self.w=zero['epsilon']+zero['P']
        self.phiscale=ld('.01')*np.sqrt(self.thermal/self.w)
        T=np.exp(self.base[1]);rho=np.exp(self.base[0]);aux=data['aux'][i]
        K=16*ld('5.670400e-5')*T**3/(3*rho*np.array([aux[27],aux[24]],dtype=ld))
        tau=np.array([two.e.TAU,1/(prior.C*rho*aux[24])],dtype=ld)
        C=K*T/(prior.C**2*tau)
        self.qscale=ld('.01')*np.sqrt(C*self.thermal)
        self.total_qscale=self.qscale.sum()
        self.vscale=self.total_qscale/self.w
        self.B=rho;self.I=16*self.w
        self.tunit=1/(prior.C*rho*aux[24]);self.Tscale=ld('.001')
        self.y0=np.array([0,0,1,1,1],dtype=ld)
        self.z0=self.state(self.y0,0);self.J=self.z0['a']**2*self.z0['S']

    def state(self,y,t):
        key=tuple(map(float,y))+(float(t),)
        if key in self.cache:return self.cache[key]
        a=np.exp(ld('.001')*np.sin(ld(t)));v=y[1]*self.vscale
        assert abs(v)<ld('.1')
        W=1/np.sqrt(1-v*v);rho=self.B/(a*W)
        logT=self.base[1]+self.Tscale*y[0];phi=self.phiscale*y[4]
        z=self.provider.state(np.log(rho),logT,phi,self.X)
        assert self.provider.calls<=self.plan['gates']['maximum_native_calls']
        parts=two.opacity_parts(self.opacity,(np.log(rho)-3*z['logA'],logT-z['logA'],self.X))
        TJ=np.exp(logT-z['logA']);rhoJ=z['rhoJ'];factor=16*ld('5.670400e-5')*TJ**3/(3*rhoJ)
        K=z['A']**2*factor/np.asarray([parts[3],parts[0]],ld)
        tau=np.array([two.e.TAU,1/(prior.C*rhoJ*parts[0])],dtype=ld)/z['A']
        T=np.exp(logT);C=K*T/(prior.C**2*tau);beta=1/(C*T)
        q=self.qscale*y[2:4];Q=q.sum();w=z['epsilon']+z['P']
        E=(z['epsilon']+z['P']*v*v+2*Q*v)*W*W
        S=(w*v+Q*(1+v*v))*W*W
        R=(z['epsilon']*v*v+z['P']+2*Q*v)*W*W
        entropy=a*W*(rho*z['entropy']+v*Q/T-np.sum(beta*q*q)/2)
        z.update(a=a,v=v,W=W,T=T,q=q,Q=Q,C=C,beta=beta,tau=tau/self.tunit,E=E,S=S,R=R,
            entropy_cell=entropy,g=a*(-4*phi)*z['trace'])
        self.cache[key]=z;return z

    def energy_increment(self,z0,z1):
        # B is fixed. Subtract rest energies through expm1, not two total energies.
        dlog=z1['logA']-z0['logA']+np.log(z1['W']/z0['W'])
        rest=self.B*z0['A']*z0['W']*(np.expm1(dlog)*(z0['restJ']+z0['uJ'])+np.exp(dlog)*(z1['uJ']-z0['uJ']))
        kinetic=lambda z:z['a']*(z['P']*z['v']**2+2*z['Q']*z['v'])*z['W']**2
        return rest+kinetic(z1)-kinetic(z0)

    def step(self,y,p,t,h,damping=1):
        old=self.state(y,t)
        def residual(next_y):
            z=self.state(np.asarray(next_y,ld),t+h)
            dphi=z['phi']-old['phi'];p1=2*dphi/h-p
            g=(z['g']+old['g'])/2
            dell=np.log(z['a']/old['a']);dlnT=z['logT']-old['logT']
            wm=(z['W']+old['W'])/2;vm=(z['v']+old['v'])/2
            cm=(z['C']+old['C'])/2;qm=(z['q']+old['q'])/2;taum=(z['tau']+old['tau'])/2
            heat=z['q']-old['q']+cm*(wm*wm*(z['v']-old['v'])+vm*(dell+dlnT))
            heat+=qm/2*(np.log(z['beta']/old['beta'])+np.log(z['W']/old['W'])+dell)+h*damping*qm/(taum*wm)
            work=-(z['a']*z['R']+old['a']*old['R'])/2*dell
            return np.array([(self.energy_increment(old,z)+g*dphi-work)/(self.thermal*self.Tscale),
                (z['a']**2*z['S']-self.J)/self.total_qscale,*list(heat/self.qscale),
                (dphi-h*p+h*h/2*((z['phi']+old['phi'])/2-g/self.I))/self.phiscale],float)
        solved=root(residual,np.asarray(y,float),method='hybr',options=dict(xtol=1e-9,eps=1e-8,maxfev=100))
        error=float(np.max(abs(residual(solved.x))))
        assert error<self.plan['gates']['stage_residual'],(t,error,solved.message)
        next_y=np.asarray(solved.x,ld);z=self.state(next_y,t+h);p1=2*(z['phi']-old['phi'])/h-p
        work=-(z['a']*z['R']+old['a']*old['R'])/2*np.log(z['a']/old['a'])
        production=lambda s:s['a']*np.sum(s['beta']*s['q']**2/s['tau'])
        return next_y,p1,work,h*damping*(production(old)+production(z))/2,error


def path(cell,steps,damping=1,zero_heat=False):
    began=time.monotonic();h=ld(str(cell.plan['duration_scaled']))/steps
    y=cell.y0.copy();p=ld(0);work=ld(0);entropy_integral=ld(0);maxres=0.;crossings=0
    if zero_heat:y[2:4]=0
    initial=cell.state(y,0);oldJ=cell.J;cell.J=initial['a']**2*initial['S']
    rows=[]
    for k in range(steps):
        y,p1,dwork,ds,res=cell.step(y,p,ld(k)*h,h,damping)
        crossings+=int(p*p1<0);p=p1;work+=dwork;entropy_integral+=ds;maxres=max(maxres,res)
        z=cell.state(y,ld(k+1)*h)
        scalar=cell.I/2*(p*p+(z['phi']+initial['phi'])*(z['phi']-initial['phi']))
        defect=float(abs(cell.energy_increment(initial,z)+scalar-work)/cell.thermal)
        assert defect<cell.plan['gates']['work_defect_over_initial_thermal'],defect
        rows.append([float((k+1)*h),*map(float,y),float(p/cell.phiscale),float(work/cell.thermal),
            float((z['entropy_cell']-initial['entropy_cell'])*initial['T']/cell.thermal),
            float(entropy_integral*initial['T']/cell.thermal),defect])
    cell.J=oldJ;history=np.asarray(rows)
    stem=f'steps-{steps}-damping-{damping}-zeroheat-{int(zero_heat)}'
    np.savez_compressed(OUT/(stem+'.npz'),history=history,endpoint=y)
    production=float(entropy_integral*initial['T']/cell.thermal)
    entropy=float((z['entropy_cell']-initial['entropy_cell'])*initial['T']/cell.thermal)
    result=dict(classification='Counterexample candidate',steps=steps,damping=damping,zero_heat=zero_heat,
        seconds=time.monotonic()-began,maximum_stage_residual=maxres,maximum_work_defect=float(history[:,-1].max()),
        entropy_change_scaled=entropy,entropy_production_integral_scaled=production,
        entropy_budget_relative=abs(entropy-production)/max(abs(production),1e-30),momentum_crossings=crossings,
        endpoint=np.r_[y,p/cell.phiscale].astype(float).tolist(),history_file=stem+'.npz')
    save(stem+'.json',result);return result


def run(pilot=False):
    plan=bindings();began=time.monotonic();cell=Cell();target='pilot.json' if pilot else 'result.json'
    assert not (OUT/target).exists()
    try:
        if pilot:
            row=path(cell,8);estimate=row['seconds']*9
            result=dict(classification='Counterexample candidate',passed=estimate<180,row=row,
                projected_remaining_seconds=estimate,native_calls=cell.provider.calls,seconds=time.monotonic()-began,
                scales=dict(phi=float(cell.phiscale),q_over_w=(cell.qscale/cell.w).astype(float).tolist(),thermal_over_w=float(cell.thermal/cell.w),unit_seconds=float(cell.tunit)))
        else:
            preview=json.loads((OUT/'pilot.json').read_text());assert preview['passed']
            rows=[preview['row'],path(cell,16),path(cell,32)]
            endpoints=np.array([r['endpoint'] for r in rows]);diff=np.max(abs(np.diff(endpoints,axis=0)),axis=1)
            order=float(np.log2(diff[0]/diff[1]))
            zero=path(cell,8,zero_heat=True)
            reversible=path(cell,8,damping=0)
            gates=dict(time_order=order>=plan['gates']['refinement_order'],
                entropy_positive=all(r['entropy_change_scaled']>0 for r in rows),
                entropy_budget=rows[-1]['entropy_budget_relative']<=plan['gates']['entropy_budget_relative'],
                scalar_turning=all(r['momentum_crossings']>0 for r in rows),
                zero_heat_invariant=max(abs(np.array(zero['endpoint'])[1:4]))<1e-12)
            result=dict(classification='Counterexample candidate',passed=all(gates.values()),gates=gates,rows=rows,
                time_order=order,endpoint_differences=diff.tolist(),zero_heat=zero,reversible=reversible,
                native_calls=cell.provider.calls,seconds=time.monotonic()-began,full_GR=False,
                finite_step_entropy_exact=False,physical_EOS_certified=False,observational_closure=False)
        save(target,result);print(json.dumps(result),flush=True)
        assert result['passed'],result
    except Exception as error:
        save('failure-'+target,dict(error=repr(error),native_calls=cell.provider.calls,seconds=time.monotonic()-began));raise
    save('manifest.json',dict(sha256={str(p.relative_to(ROOT)):digest(p) for p in OUT.iterdir() if p.is_file() and p.name!='manifest.json'}))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['symbolic','prepare','pilot','run'])
    action=parser.parse_args().action
    if action in ['pilot','run']:run(action=='pilot')
    else:print(globals()[action]())
