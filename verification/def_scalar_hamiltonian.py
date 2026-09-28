"""Regular discrete-gradient evolution of the spherical vacuum scalar sector.

Proven (conditional algebra): the polarized mass constraint supplies a smooth,
symmetric discrete gradient of the outer mass. No division by field momentum,
field increments, or an energy residual is used. G=c=1, outer Dirichlet phi=0.
Counterexample candidate: finite numerical controls below; no native matter.
"""
from pathlib import Path
import argparse
import hashlib
import json
import time
import numpy as np
import sympy as sp
from scipy.linalg import solve_banded

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/direct-eos-gr33/def-scalar-hamiltonian'
ld=np.longdouble
PI=ld(str(np.pi))


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n')


class Scalar:
    def __init__(self,n,gravity=True):
        self.n,self.gravity=n,gravity
        self.rf=np.linspace(ld(0),ld(1),n+1,dtype=ld)
        self.r=(self.rf[:-1]+self.rf[1:])/2
        self.volume=4*PI/3*np.diff(self.rf**3)
        self.w=self.volume/(4*PI)
        self.area=4*PI*self.rf**2
        self.distance=np.r_[ld(1),np.diff(self.r),self.rf[-1]-self.r[-1]]
        self.face_weight=self.area*self.distance
        self.fraction=(self.r**3-self.rf[:-1]**3)/np.diff(self.rf**3)
        self.inverse_faces=np.r_[ld(0),1/self.rf[1:]]

    def gradient(self,phi):
        return np.r_[ld(0),np.diff(phi)/np.diff(self.r),-phi[-1]/self.distance[-1]]

    def divergence(self,flux):return np.diff(self.area*flux)/self.volume

    def coefficients(self,kin,left,right):
        if not self.gravity:return np.ones(self.n,dtype=ld),np.ones(self.n,dtype=ld)
        C=1-2*kin*(1-self.fraction)/self.r-2*left*self.inverse_faces[:-1]
        D=1+2*kin*self.fraction/self.r+2*right*self.inverse_faces[1:]
        assert np.min(C)>0 and np.min(D)>0,'Outside the declared positive radial-adjoint branch'
        return C,D

    def state(self,y):
        phi,p=y
        Phi=self.gradient(phi)
        kin=self.w*p*p/2
        face=self.face_weight*Phi*Phi/(8*PI)
        left,right=face[:-1]/2,face[1:]/2
        right[-1]*=2
        C,D=self.coefficients(kin,left,right)
        mf=np.zeros(self.n+1,dtype=ld)
        for i in range(self.n):mf[i+1]=(C[i]*mf[i]+kin[i]+left[i]+right[i])/D[i]
        mass=mf[:-1]+self.fraction*np.diff(mf)
        b=1-2*mass/self.r if self.gravity else np.ones(self.n,dtype=ld)
        bf=1-2*mf*self.inverse_faces if self.gravity else np.ones(self.n+1,dtype=ld)
        assert np.min(b)>0 and np.min(bf)>0,'No trapped surface permitted in this control'
        return dict(phi=phi,p=p,Phi=Phi,kin=kin,left=left,right=right,C=C,D=D,mf=mf,m=mass,b=b,bf=bf)

    def pair(self,old,new):
        C,D=self.coefficients((old['kin']+new['kin'])/2,(old['left']+new['left'])/2,(old['right']+new['right'])/2)
        H=np.empty(self.n,dtype=ld);adjoint=ld(1)
        for i in range(self.n-1,-1,-1):
            H[i]=adjoint/D[i]
            adjoint*=C[i]/D[i]
        b,bf=(old['b']+new['b'])/2,(old['bf']+new['bf'])/2
        Hf=np.r_[H[0],(H[:-1]+H[1:])/2,H[-1]]
        return dict(C=C,D=D,H=H,k=H*b,kf=Hf*bf,b=b,bf=bf)

    def discrete_gradient(self,old,new):
        a=self.pair(old,new)
        p=(old['p']+new['p'])/2
        Phi=(old['Phi']+new['Phi'])/2
        return np.array([-self.w*self.divergence(a['kf']*Phi),self.w*a['k']*p])

    def mass_increment(self,old,new):
        a=self.pair(old,new)
        dkin=self.w*(new['p']+old['p'])*(new['p']-old['p'])/2
        dface=self.face_weight*(new['Phi']+old['Phi'])*(new['Phi']-old['Phi'])/(8*PI)
        dl,dr=dface[:-1]/2,dface[1:]/2;dr[-1]*=2
        dm=ld(0)
        for i in range(self.n):
            dm=(a['C'][i]*dm+a['b'][i]*dkin[i]+a['bf'][i]*dl[i]+a['bf'][i+1]*dr[i])/a['D'][i]
        return dm

    def step(self,y,h,amplitude):
        old=self.state(y);guess=y.copy()
        for iteration in range(16):
            new=self.state(guess);a=self.pair(old,new)
            conductance=self.area*a['kf']/self.distance;conductance[0]=0
            f=h*h*a['k']/(4*self.volume)
            lower=-f[1:]*conductance[1:-1];upper=-f[:-1]*conductance[1:-1]
            diagonal=1+f*(conductance[:-1]+conductance[1:])
            rhs=y[0]+h*a['k']*y[1]/2
            band=np.zeros((3,self.n));band[1]=diagonal
            band[0,1:],band[2,:-1]=upper,lower
            mid=solve_banded((1,1),band,np.asarray(rhs,float)).astype(ld)
            for _ in range(2):
                error=rhs-diagonal*mid
                error[1:]-=lower*mid[:-1];error[:-1]-=upper*mid[1:]
                mid+=solve_banded((1,1),band,np.asarray(error,float)).astype(ld)
            p=y[1]+h*self.divergence(a['kf']*self.gradient(mid))
            phi=y[0]+h*a['k']*(y[1]+p)/2
            candidate=np.array([phi,p])
            difference=float(np.max(abs(candidate-guess))/amplitude)
            guess=candidate
            if difference<2e-17:break
        else:raise RuntimeError(('Nonlinear iteration cap',difference))
        new=self.state(guess);a=self.pair(old,new)
        defect=np.array([guess[0]-y[0]-h*a['k']*(y[1]+guess[1])/2,
            guess[1]-y[1]-h*self.divergence(a['kf']*(old['Phi']+new['Phi'])/2)])
        residual=float(np.max(abs(defect))/amplitude)
        assert residual<2e-15,residual
        return guess,dict(iterations=iteration+1,residual=residual,
            stable_mass_increment=float(self.mass_increment(old,new)),
            mass=float(new['mf'][-1]),minimum_b=float(np.min(new['b'])))


def initial(model,amplitude):
    return np.array([ld(amplitude)*np.sinc(model.r),np.zeros(model.n,dtype=ld)])


def symbolic():
    a0,a1,L0,L1,R0,R1,x0,x1,y0,y1,r,l,u,f=sp.symbols('a0 a1 L0 L1 R0 R1 x0 x1 y0 y1 r l u f',nonzero=True)
    constraint=lambda a,L,R,x,y:y-x-a*(1-2*((1-f)*x+f*y)/r)-L*(1-2*x/l)-R*(1-2*y/u)
    a,L,R=(a0+a1)/2,(L0+L1)/2,(R0+R1)/2
    C=1-2*a*(1-f)/r-2*L/l;D=1+2*a*f/r+2*R/u
    b=1-((1-f)*(x0+x1)+f*(y0+y1))/r
    polarized=D*(y1-y0)-C*(x1-x0)-b*(a1-a0)-(1-(x0+x1)/l)*(L1-L0)-(1-(y0+y1)/u)*(R1-R0)
    assert sp.expand(constraint(a1,L1,R1,x1,y1)-constraint(a0,L0,R0,x0,y0)-polarized)==0
    q0,q1=sp.symbols('q0 q1');assert sp.expand(q1*q1-q0*q0-(q1+q0)*(q1-q0))==0
    return dict(classification='Proven',passed=True,
        scope='Exact polarized radial mass-constraint chain rule. Radial adjoint telescoping and shared-face summation by parts yield the discrete gradient; symmetric smooth coefficients give second-order time consistency on the positive regular branch. No physical matter/heat closure or global PDE error bound.')


def check():
    symbolic();s=Scalar(12);r=s.r
    y0=np.array([ld('.04')*np.sinc(r),ld('.01')*np.sin(2*PI*r)])
    y1=np.array([ld('.031')*np.sinc(r)+ld('.003')*np.sin(PI*r),-ld('.008')*np.sin(PI*r)])
    old,new=s.state(y0),s.state(y1)
    delta=s.mass_increment(old,new);gradient=s.discrete_gradient(old,new)
    assert abs(delta-np.sum(gradient*(y1-y0)))<ld('2e-18')*old['mf'][-1]
    assert abs(delta-(new['mf'][-1]-old['mf'][-1]))<ld('2e-17')*old['mf'][-1]
    assert np.array_equal(gradient,s.discrete_gradient(new,old))
    stationary=np.zeros((2,s.n),dtype=ld)
    same,info=s.step(stationary,ld('.01'),ld('.04'));assert np.all(same==0)
    turning=initial(s,'.04');forward,_=s.step(turning,ld('.01'),ld('.04'))
    restored,_=s.step(forward,ld('-.01'),ld('.04'))
    reversal=float(np.max(abs(restored-turning))/ld('.04'));assert reversal<2e-14
    assert np.all(np.isfinite(s.discrete_gradient(s.state(turning),s.state(turning))))
    return dict(classification='Counterexample candidate',passed=True,
        secant_relative_defect=float(abs(delta-np.sum(gradient*(y1-y0)))/old['mf'][-1]),
        time_reversal_relative_defect=reversal,zero_state_exact=True,all_momenta_zero_start_passed=True)


def evolve(n,steps,gravity,amplitude,duration=ld('1.7')):
    start=time.monotonic();s=Scalar(n,gravity);y=initial(s,amplitude);M0=s.state(y)['mf'][-1]
    h=duration/steps;rows=[];crossings=0
    for step in range(steps):
        next_y,row=s.step(y,h,ld(amplitude))
        crossings+=int(np.sum(y[1]*next_y[1]<0))
        y=next_y;rows.append(row)
    maximum_energy=max(abs(ld(row['stable_mass_increment']))/M0 for row in rows)
    global_energy=abs(s.state(y)['mf'][-1]-M0)/M0
    assert maximum_energy<5e-14 and global_energy<2e-12,(maximum_energy,global_energy)
    return y,dict(classification='Counterexample candidate',cells=n,steps=steps,gravity=gravity,
        amplitude=float(amplitude),duration=float(duration),seconds=time.monotonic()-start,
        maximum_equation_residual=max(row['residual'] for row in rows),
        maximum_stable_relative_mass_increment=float(maximum_energy),relative_final_mass_drift=float(global_energy),
        maximum_iterations=max(row['iterations'] for row in rows),momentum_sign_crossings=crossings,
        minimum_b=min(row['minimum_b'] for row in rows))


def prepare():
    assert not OUT.exists();OUT.mkdir()
    save('plan.json',dict(classification='Counterexample candidate',source_sha256=digest(Path(__file__)),
        predecessor_milestone_sha256=digest(ROOT/'outputs/direct-eos-gr33/gr-spherical-coupling-milestone-manifest.json'),
        decision='Replace the inverse-momentum energy projection with a smooth polarized mass-constraint discrete gradient; verify zero/turning states, exact conservation, nonlinear GR evolution and time/space consistency before native stellar reuse.',
        symbolic=symbolic(),operator_check=check(),
        controls=dict(pilot=[32,32,.4],flat_space_cells=[32,64,128],flat_steps_per_cell=4,
            nonlinear_time_cells=64,nonlinear_time_steps=[128,256,512],
            nonlinear_space_cells=[32,64,128],nonlinear_steps_per_cell=4,
            amplitude=.05,duration=1.7,initial='phi=A*sin(pi*r)/(pi*r), Pi=0',boundary='regular center, fixed outer phi=0'),
        gates=dict(equation_residual=2e-15,stable_step_relative_mass_increment=5e-14,
            final_relative_mass_drift=2e-12,time_order_minimum=1.8,
            flat_exact_L2_order_minimum=1.7,nonlinear_space_L2_order_minimum=1.6,
            turning_crossings_positive=True,minimum_nontrivial_compactness=1e-5),
        budget=dict(cpu=1,blas_threads=1,gpu=False,pilot_timeout_seconds=30,
            full_timeout_seconds=180,maximum_runs=1,launch_only_if_pilot_scaled_seconds_below=120,
            estimate='Measure the 32-cell 32-step nonlinear pilot; extrapolate by total cell-steps with a factor two safety margin. Reuse identical n=64,steps=256 path. No new native EOS calls, automatic grids, amplitudes or duration changes.'),
        native_stellar_matter_included=False,physical_EOS_certified=False,observational_closure=False))
    print('PREPARED regular scalar Hamiltonian controls',flush=True)


def bindings():
    p=json.loads((OUT/'plan.json').read_text());assert digest(Path(__file__))==p['source_sha256'];return p


def pilot():
    bindings();assert not (OUT/'pilot.json').exists()
    _,row=evolve(32,32,True,'.05',ld('.4'))
    cellsteps=sum(n*4*n for n in [32,64,128])*2+64*(128+512)
    estimated=2*row['seconds']*cellsteps/(32*32)
    result=dict(classification='Counterexample candidate',passed=estimated<120,control=row,estimated_full_seconds=estimated)
    save('pilot.json',result);assert result['passed'],result
    print('PILOT',json.dumps(result),flush=True)


def weighted_norm(model,delta):return np.sqrt(np.sum(delta*delta*model.volume,axis=1)/np.sum(model.volume))


def run():
    plan=bindings();assert json.loads((OUT/'pilot.json').read_text())['passed']
    assert not (OUT/'result.json').exists();start=time.monotonic();paths={};metrics=[]
    jobs=[(False,n,4*n) for n in [32,64,128]]+[(True,64,k) for k in [128,256,512]]+[(True,32,128),(True,128,512)]
    try:
        for gravity,n,steps in jobs:
            y,row=evolve(n,steps,gravity,'.05');paths[gravity,n,steps]=y;metrics.append(row)
            stem=f'{"gr" if gravity else "flat"}-{n}-{steps}'
            np.savez_compressed(OUT/f'{stem}.npz',state=y)
            save(stem+'.json',row);print('PATH',stem,row['seconds'],row['relative_final_mass_drift'],flush=True)
        exact_errors=[]
        for n in [32,64,128]:
            s=Scalar(n,False);profile=ld('.05')*np.sinc(s.r)
            exact=np.array([profile*np.cos(PI*ld('1.7')),-PI*profile*np.sin(PI*ld('1.7'))])
            exact_errors.append(weighted_norm(s,paths[False,n,4*n]-exact))
        exact_orders=np.log2(np.array(exact_errors[:-1])/np.array(exact_errors[1:]))
        s=Scalar(64);dt1=weighted_norm(s,paths[True,64,128]-paths[True,64,256])
        dt2=weighted_norm(s,paths[True,64,256]-paths[True,64,512]);time_orders=np.log2(dt1/dt2)
        space_differences=[]
        for n in [32,64]:
            fine=paths[True,2*n,8*n];restriction=(fine[:,::2]+fine[:,1::2])/2
            space_differences.append(weighted_norm(Scalar(n),paths[True,n,4*n]-restriction))
        space_orders=np.log2(space_differences[0]/space_differences[1])
        gates=dict(flat_exact_orders=bool(np.min(exact_orders)>=1.7),time_orders=bool(np.min(time_orders)>=1.8),
            nonlinear_space_orders=bool(np.min(space_orders)>=1.6),
            turning_crossings=all(row['momentum_sign_crossings']>0 for row in metrics),
            nonlinear_metric=all(1-row['minimum_b']>1e-5 for row in metrics if row['gravity']))
        result=dict(classification='Counterexample candidate',passed=all(gates.values()),gates=gates,
            flat_exact_errors=np.array(exact_errors).astype(float).tolist(),flat_exact_orders=exact_orders.astype(float).tolist(),
            nonlinear_time_orders=time_orders.astype(float).tolist(),nonlinear_space_orders=space_orders.astype(float).tolist(),
            nonlinear_time_differences=[dt1.astype(float).tolist(),dt2.astype(float).tolist()],
            nonlinear_space_differences=np.array(space_differences).astype(float).tolist(),
            maximum_relative_mass_drift=max(row['relative_final_mass_drift'] for row in metrics),
            maximum_stable_step_relative_mass_increment=max(row['maximum_stable_relative_mass_increment'] for row in metrics),
            maximum_equation_residual=max(row['maximum_equation_residual'] for row in metrics),
            minimum_momentum_sign_crossings=min(row['momentum_sign_crossings'] for row in metrics),
            seconds=time.monotonic()-start,native_stellar_matter_included=False,
            uniform_PDE_error_bound=False,physical_EOS_certified=False,observational_closure=False)
        save('result.json',result);print('VERDICT',json.dumps(result),flush=True)
    except Exception as error:
        save('failure.json',dict(classification='Counterexample candidate',error=repr(error),seconds=time.monotonic()-start));raise
    save('manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():digest(p) for p in OUT.iterdir() if p.is_file() and p.name!='manifest.json'}))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['check','prepare','pilot','run'])
    print(globals()[parser.parse_args().action]())
