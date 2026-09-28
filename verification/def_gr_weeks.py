"""Original GR resolvent in a Laguerre time basis (Weeks inversion).

No projected spectrum, mode removal, changed inertia or smoothed input.
Counterexample candidate until every original output and contrast passes.
"""
from pathlib import Path
import argparse
import json
import resource
import signal
import time
import numpy as np
import sympy as sp
from scipy.sparse.linalg import splu
import def_gr_matrix_propagation as setup

prior=setup.prior
OUT=setup.OUT/'weeks';write=setup.write
COUNT=4096;BETA=512.;DEGREES=[512,1024,2048];TIMES=prior.TIMES


def laguerre(sigma):
    L=np.zeros((len(TIMES),DEGREES[-1]));L[:,0]=np.exp(-BETA*TIMES)
    L[:,1]=(1-2*BETA*TIMES)*L[:,0]
    for n in range(1,L.shape[1]-1):L[:,n+1]=((2*n+1-2*BETA*TIMES)*L[:,n]-n*L[:,n-1])/(n+1)
    return np.exp(sigma*TIMES[:,None])*L


def coefficients(half):
    # Uniform circle, including exact G(1)=0. Endpoint Nyquist is real.
    full=np.concatenate([half,half[-2:0:-1].conj()],axis=0)
    return (np.fft.fft(full,axis=0)/len(full)).real


def control():
    theta=2*np.pi*np.arange(1,COUNT//2+1)/COUNT;z=np.exp(1j*theta)
    s=12+BETA*(1+z)/(1-z);errors={};L=laguerre(12.)
    for omega in [3.,50.,300.,600.,800.,2000.]:
        G=np.r_[0.,2*BETA/(1-z)/(s*(s*s+omega*omega))]
        G[-1]=G[-1].real;a=coefficients(G[:,None])[:,0]
        exact=(1-np.cos(omega*TIMES))/omega**2
        errors[str(omega)]=[float(max(abs((L[:,:n]@a[:n])[1:]-exact[1:]))*omega**2) for n in DEGREES]
    assert max(errors[str(w)][-1] for w in [3.,50.,300.,600.,800.])<1e-6
    t,sig,b,s=sp.symbols('t sig b s',positive=True)
    for n in range(4):
        expr=sp.exp((sig-b)*t)*sp.laguerre(n,2*b*t)
        assert sp.simplify(sp.laplace_transform(expr,t,s,noconds=True)-(s-sig-b)**n/(s-sig+b)**(n+1))==0
    return dict(classification='Proven',symbolic_identity_passed=True,
        identity='Laplace[exp((sigma-beta)t) L_n(2 beta t)] = (s-sigma-beta)^n/(s-sigma+beta)^(n+1)',
        numerical_classification='Counterexample candidate',oscillator_errors=errors,
        unresolved_high_frequency_control='omega2000 remains unresolved; never claim a uniform all-mode guarantee.')


class Problem:
    def __init__(self,bank,outer=2):
        self.p=prior.Problem(bank,outer);self.L,self.source_check=setup.source_map(self.p)
        self.error=0.;self.K=self.p.K;self.D=self.p.D
    def transform(self,s):
        p=self.p;heat=p.heat;lam=heat.rates*p.radiation.geometry.tc
        energy=np.zeros(len(heat.edges),complex)
        energy[heat.face_ids]=np.sum(heat.amplitude*p.radiation.geometry.tc*lam/(s*s*(s+lam)),axis=1)
        rhs=self.L@energy;matrix=(self.K-s*s*self.D).tocsc();scale=np.asarray(abs(matrix).sum(1)).ravel()
        lu=splu(matrix.multiply((1/scale)[:,None]).tocsc());answer=lu.solve(rhs/scale)
        extended=matrix.astype(np.clongdouble);r=rhs.astype(np.clongdouble)
        for _ in range(2):answer+=lu.solve(np.asarray((r-extended@answer)/scale,complex))
        residual=float(np.max(abs(r-extended@answer)/(abs(extended)@abs(answer)+abs(r)+1e-100)))
        self.error=max(self.error,residual);assert self.error<1e-9
        return answer


def contour(k,sigma):
    z=np.exp(2j*np.pi*k/COUNT)
    return sigma+BETA*(1+z)/(1-z),2*BETA/(1-z)


def prepare():
    assert not OUT.exists();OUT.mkdir()
    paths=[Path(__file__),Path(setup.__file__),Path(prior.__file__),prior.prior.OUT/'fine-bank.npz',
           prior.prior.OUT/'coarse-bank.npz',prior.OUT/'result.json']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        rationale='Unaccelerated Fourier cutoff and finite-precision Q-D fail the local wave; projected full/finite Krylov matrices introduce fast growing modes. Use a linear Laguerre representation of the unchanged resolvent, validated on known frequencies through800 before production.',
        method='sigma12,beta512;4096 uniform mapped-contour nodes (2048 finite solves plus conjugacy and G(1)=0). Coefficient counts512/1024/2048; same65 physical times. This is a time-basis truncation, not an assertion of passing old16/32/64 RK steps.',
        retained='Same K,D, exact microscopic heat poles, face-energy debit, momentum lift, initial condition, radial grid, duration and all four outputs. No eigenvalue clipping, modal filter or window.',
        decision='Run only one fine-bank contour after measured pilot. If any component or nested contour comparison fails, stop. Only if it passes, execute the three fixed coefficient/outer/sigma14 contrasts, within the remaining measured budget. No automatic degree/contour/grid expansion.',
        gates=dict(relative=.02,order=1.5,coefficient=.02,outer=.002,contour=.0002,abscissa=.0002,linear_residual=1e-9,heat_balance=2e-13),
        budget=dict(pilot_resolvents=16,maximum_total_resolvents=8192,main_hard_seconds=180,total_hard_seconds=600,CPU_threads=1,memory_GB=3,new_EOS_calls=0),
        cost_boundary='A fixed replacement propagator, not additional physical trajectories. Reuse pilot points and nested coarse contour. Main forecast must fit180s; all four contours must fit600s before production.',
        bindings={str(p):prior.prior.digest(p) for p in paths}))
    write(OUT/'control.json',control())
    signal.alarm(60);start=time.monotonic();p=Problem(prior.prior.OUT/'fine-bank.npz');setup_seconds=time.monotonic()-start
    ids=np.linspace(1,COUNT//2,16,dtype=int);start=time.monotonic()
    answers=np.array([p.transform(contour(k,12.)[0]) for k in ids]);seconds=time.monotonic()-start
    np.savez_compressed(OUT/'pilot.npz',ids=ids,answers=answers)
    per=seconds/len(ids);main=setup_seconds+per*2048*1.5+20;total=4*main
    write(OUT/'pilot.json',dict(classification='Counterexample candidate',setup_seconds=setup_seconds,solve_seconds=seconds,
        main_forecast_seconds=main,total_forecast_seconds=total,linear_residual=p.error,
        assumption='Representative logarithmically mapped frequencies,1.5 times measured solve mean plus20s inversion per contour. Other banks and enlarged exterior unmeasured.',source_check=p.source_check))
    signal.alarm(0);print('PILOT',main,total,flush=True)


def solve(bank,outer,sigma,label,reuse=False):
    start=time.monotonic();problem=Problem(bank,outer);p=problem.p;native=p.native;N=len(native)
    half=np.zeros((COUNT//2+1,2*N),complex)
    L=laguerre(sigma);end=np.zeros((3,len(p.bg.grid)*4),np.longdouble);endv=end.copy()
    kernels=[]
    for n in DEGREES:
        ell=np.zeros(COUNT);ell[:n]=L[-1,:n];kernels.append(np.fft.fft(ell)[:COUNT//2+1]/COUNT)
    kernels=np.array(kernels)
    pilot=dict(np.load(OUT/'pilot.npz')) if reuse else None
    saved={int(k):v for k,v in zip(pilot['ids'],pilot['answers'])} if reuse else {}
    for k in range(1,COUNT//2+1):
        s,factor=contour(k,sigma);answer=saved[k] if k in saved else problem.transform(s)
        q=answer.reshape(-1,4)
        half[k,:N]=factor*s*np.interp(native,p.bg.grid,q[:,0]*p.speed)
        half[k,N:]=factor*np.interp(native,p.bg.grid,q[:,2]-p.bg.grid*p.bg.nodes['v']*q[:,0])
        mult=(1 if k==COUNT//2 else 2)*kernels[:,k]*factor
        for j in range(3):
            end[j]+=np.asarray((mult[j]*answer).real,np.longdouble)
            endv[j]+=np.asarray((mult[j]*s*answer).real,np.longdouble)
    half[-1]=half[-1].real
    histories=np.zeros((3,65,2*N));coarse=np.zeros((65,2*N))
    # Column blocks keep the full source resolvent and FFT memory bounded.
    for first in range(0,2*N,256):
        last=min(first+256,2*N);a=coefficients(half[:,first:last]);c=coefficients(half[::2,first:last])
        for j,n in enumerate(DEGREES):histories[j,:,first:last]=L[:,:n]@a[:n]
        coarse[:,first:last]=L@c[:DEGREES[-1]]
    del half
    lift=np.array([np.interp(native,p.bg.grid,p.speed*p.heat.lift(t,p.bg.nodes,True).reshape(-1,4)[:,0]) for t in TIMES])
    flux,energy=p.heat.faces(1.);balance=float(abs(np.sum(-np.diff(energy),dtype=np.longdouble))/max(abs(energy).max(),1e-100))
    cases={}
    def readouts(state):
        velocity=state[:,:N]-lift;scalar=state[:,N:];velocity[0]=0;scalar[0]=0
        return np.column_stack([np.sqrt((velocity*velocity)@p.weights),np.sqrt((scalar*scalar)@p.weights)]+
            [np.sqrt((velocity[:,m]**2)@p.weights[m]/p.weights[m].sum()) for m in p.masks]),velocity,scalar
    for j,n in enumerate(DEGREES):
        values,velocity,scalar=readouts(histories[j]);cases[n]=values
        q=np.asarray(end[j],float).reshape(-1,4)-p.heat.lift(1.,p.bg.nodes).reshape(-1,4)
        v=np.asarray(endv[j],float).reshape(-1,4)-p.heat.lift(1.,p.bg.nodes,True).reshape(-1,4)
        np.savez_compressed(OUT/f'{label}-{n}.npz',times=TIMES,readouts=values,velocity_native=velocity,scalar_native=scalar,
            grid=p.bg.grid,response=q,velocity=v,heat_energy=energy,heat_luminosity=flux)
    coarse_values,_,_=readouts(coarse);comparison={}
    for j,field in enumerate(prior.FIELDS):
        a,b,c=[cases[n][:,j] for n in DEGREES];norm=max(abs(c).max(),1e-100)
        d1=float(max(abs(a-b))/norm);d2=float(max(abs(b-c))/norm)
        comparison[field]=dict(previous=d1,last=d2,order=float(np.log2(d1/d2)),contour=float(max(abs(c-coarse_values[:,j]))/norm))
    row=dict(classification='Counterexample candidate',comparisons=comparison,seconds=time.monotonic()-start,
        source_check=problem.source_check,linear_residual=problem.error,heat_balance=balance,sigma=sigma,outer=outer,
        endpoint=cases[2048][-1].tolist(),original_failure_resolved=False)
    write(OUT/(label+'.json'),row);print('PATH',label,json.dumps(row),flush=True)
    return cases[2048],row


def run():
    assert not (OUT/'result.json').exists()
    plan=json.loads((OUT/'plan.json').read_text());pilot=json.loads((OUT/'pilot.json').read_text())
    for p,h in plan['bindings'].items():assert prior.prior.digest(Path(p))==h,p
    assert pilot['main_forecast_seconds']<180 and pilot['total_forecast_seconds']<600
    signal.alarm(180);resource.setrlimit(resource.RLIMIT_AS,(int(3e9),int(3e9)));start=time.monotonic()
    fine,row=solve(prior.prior.OUT/'fine-bank.npz',2,12.,'fine',True);rows={'fine':row}
    passed=all(v['last']<.02 and v['order']>1.5 and v['contour']<.0002 for v in row['comparisons'].values())
    contrasts={}
    if passed:
        signal.alarm(max(1,int(600-(time.monotonic()-start))))
        for key,bank,outer,sigma in [('coefficient',prior.prior.OUT/'coarse-bank.npz',2,12.),('outer',prior.prior.OUT/'fine-bank.npz',3,12.),('abscissa',prior.prior.OUT/'fine-bank.npz',2,14.)]:
            values,r=solve(bank,outer,sigma,key);rows[key]=r
            contrasts[key]={field:float(max(abs(fine[:,j]-values[:,j]))/max(abs(fine[:,j]).max(),1e-100)) for j,field in enumerate(prior.FIELDS)}
        passed=all(max(v.values())<{'coefficient':.02,'outer':.002,'abscissa':.0002}[k] for k,v in contrasts.items())
    passed=passed and max(r['linear_residual'] for r in rows.values())<1e-9 and max(r['heat_balance'] for r in rows.values())<2e-13
    result=dict(classification='Counterexample candidate',passed=passed,comparisons=row['comparisons'],contrasts=contrasts,
        seconds=time.monotonic()-start,memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
        actual_GR_descriptor_propagated=True,original_same_input_failure_resolved=passed,
        full_temperature_feedback=False,photon_transport=False,whole_star_heat_closed=False,full_dynamic_charge_solved=False)
    write(OUT/'result.json',result);signal.alarm(0);print('RESULT',json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
