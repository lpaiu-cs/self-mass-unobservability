"""Unreduced hierarchical GR with extended original resolvent and Weeks."""
from pathlib import Path
import argparse
import json
import signal
import time
import resource
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import splu
from scipy.fft import fft
import def_gr_weeks as legacy
import def_gr_response_basis as prior

task=prior.task;space=prior.space;transfer=prior.transfer
OUT=task.OUT.parent/'def-gr-full-weeks';write=task.write
COUNT=4096;BETA=512;DEGREES=[512,1024,2048];TIMES=np.arange(65,dtype=np.longdouble)/64


def contour(k,sigma):
    z=np.exp(2j*np.longdouble(np.pi)*k/COUNT)
    return np.clongdouble(sigma)+BETA*(1+z)/(1-z),2*BETA/(1-z)


def laguerre(sigma):
    L=np.empty((65,2048),dtype=np.longdouble);L[:,0]=np.exp(-BETA*TIMES)
    L[:,1]=(1-2*BETA*TIMES)*L[:,0]
    for n in range(1,2047):L[:,n+1]=((2*n+1-2*BETA*TIMES)*L[:,n]-n*L[:,n-1])/(n+1)
    return np.exp(sigma*TIMES[:,None])*L


def coefficients(half):
    values=np.concatenate([half,half[-2:0:-1].conj()])
    return (fft(values,axis=0)/len(values)).real


class Problem:
    def __init__(self,degree=4,bank=None,outer=2,quadrature=6):
        self.model=m=space.Model(degree,bank or task.BANK/'fine-bank.npz',outer,quadrature)
        self.K=(m.K+m.K.T)/2;self.M=(m.M+m.M.T)/2
        self.Kx=self.K.astype(np.clongdouble);self.Mx=self.M.astype(np.clongdouble)
        self.load=(m.load+(self.K-m.K)@m.H).astype(np.clongdouble)
        self.error=0.;r=m.original.native;N,a=m.original.radiation.geometry.metric(r)
        self.speed=a/N*r*task.fem.base.task.h.gr.C

    def transform(self,z):
        m=self.model;heat=m.heat;tc=np.longdouble(m.original.radiation.geometry.tc)
        lam=heat.rates.astype(np.longdouble)*tc;E=np.zeros(len(heat.edges),np.clongdouble)
        E[heat.face_ids]=np.sum(heat.amplitude.astype(np.longdouble)*tc*lam/(z*z*(z+lam)),axis=1)
        rhs=self.load@E
        matrix=(self.K+complex(z*z)*self.M).tocsc();scale=np.sqrt(abs(matrix.diagonal()));D=diags(1/scale)
        lu=splu((D@matrix@D).tocsc());u=(lu.solve(np.asarray(rhs/scale,complex))/scale).astype(np.clongdouble)
        # Keep K and z^2 M separate in the residual; small inertia must not be
        # lost while adding it to a large stored stiffness diagonal.
        for _ in range(3):
            defect=rhs-self.Kx@u-z*z*(self.Mx@u)
            u+=(lu.solve(np.asarray(defect/scale,complex))/scale).astype(np.clongdouble)
        defect=rhs-self.Kx@u-z*z*(self.Mx@u)
        error=float(np.max(abs(defect)/(abs(self.Kx)@abs(u)+abs(z*z)*(abs(self.Mx)@abs(u))+abs(rhs)+1e-100)))
        self.error=max(self.error,error);assert self.error<1e-9
        return u


def prepare():
    assert not OUT.exists();OUT.mkdir();signal.alarm(90)
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        claim='Remove the measured projection error by inverting the unreduced same-source GR resolvent in the existing Laguerre time basis.',
        method='Reuse original Weeks4096 mapped-contour points,sigma12,beta512 and512/1024/2048 coefficient comparison. Same degree4 matrix,source,65 times and horizon. Symmetric diagonal LU preconditions extended original K*u+s^2*M*u residual; accepted resolvents,FFT and Laguerre recurrence retain extended precision. No modes are clipped or projected.',
        gates=dict(propagation_relative=.02,propagation_order=1.5,contour=.0002,abscissa=.0002,spatial=.02,spatial_decrease=True,coefficient=.02,outer=.002,quadrature=.002,linear_residual=1e-9,heat_balance=2e-13),
        decision='Measure16 representative resolvents,cache and reuse them. Run one degree4 contour only if forecast fits450s. Stop on any of four propagation or contour failures. Fixed degree2/1 and original contrasts are conditional and require cost re-evaluation within1200s total; never expand contour,degree,period or coefficient count after failure.',
        budget=dict(pilot_seconds=90,first_case_seconds=450,total_seconds=1200,CPU_threads=1,memory_GB=4,new_EOS_calls=0),
        bindings={str(p):task.digest(p) for p in [Path(__file__),Path(legacy.__file__),Path(space.__file__),Path(task.__file__),prior.OUT/'stage-result.json',task.BANK/'fine-bank.npz',task.BANK/'coarse-bank.npz']}))
    write(OUT/'control.json',legacy.control());start=time.monotonic();p=Problem();setup=time.monotonic()-start
    ids=np.linspace(1,COUNT//2,16,dtype=int);start=time.monotonic()
    answers=np.array([p.transform(contour(k,12)[0]) for k in ids]);elapsed=time.monotonic()-start
    np.savez_compressed(OUT/'pilot.npz',ids=ids,answers=answers)
    # Measure an actual long-precision native-block inversion, rather than
    # assuming BLAS double speed for long-double matrix multiplication.
    L=laguerre(12);sample=np.ones((COUNT//2+1,128),dtype=np.clongdouble)
    start=time.monotonic();a=coefficients(sample)
    for n in DEGREES:test=L[:,:n]@a[:n]
    coarse=L@coefficients(sample[::2])[:2048];block=time.monotonic()-start
    nr=len(p.model.original.native);inversion=block*(2*nr/128)
    forecast=1.4*(setup+elapsed/16*2048+inversion+20)
    write(OUT/'pilot-budget.json',dict(classification='Counterexample candidate',setup_seconds=setup,resolvent16_seconds=elapsed,
        inversion128_columns_seconds=block,inversion_forecast_seconds=inversion,first_case_forecast_seconds=forecast,linear_residual=p.error,
        assumptions='Representative16 mapped frequencies scaled to2048,actual128-column long-precision inversion scaled to all native fields,20s endpoint/readout/storage allowance and40 percent margin. Other degrees/outer domain unmeasured.'))
    signal.alarm(0);print('FULL WEEKS FORECAST',forecast,flush=True)


def solve(degree=4,bank=None,outer=2,quadrature=6,sigma=12,label='p4',reuse=False):
    started=time.monotonic();p=Problem(degree,bank,outer,quadrature);m=p.model;nr=len(m.original.native)
    half=np.zeros((COUNT//2+1,2*nr),dtype=np.clongdouble);L=laguerre(sigma);kernels=[]
    for n in DEGREES:
        ell=np.zeros(COUNT,dtype=np.longdouble);ell[:n]=L[-1,:n];kernels.append(fft(ell)[:COUNT//2+1]/COUNT)
    kernels=np.array(kernels);end=np.zeros((3,m.size),dtype=np.longdouble);endv=end.copy()
    pilot=np.load(OUT/'pilot.npz') if reuse else None
    saved={int(i):u for i,u in zip(pilot['ids'],pilot['answers'])} if reuse else {}
    for k in range(1,COUNT//2+1):
        z,factor=contour(k,sigma);u=saved[k] if k in saved else p.transform(z)
        half[k,:nr]=factor*z*p.speed*(m.nativeV[0]@u);half[k,nr:]=factor*(m.nativeV[1]@u)
        mult=(1 if k==COUNT//2 else 2)*kernels[:,k]*factor
        for j in range(3):
            end[j]+=(mult[j]*u).real;endv[j]+=(mult[j]*z*u).real
        if k%256==0:print('FULL WEEKS',label,k,round(time.monotonic()-started,2),flush=True)
    half[-1]=half[-1].real
    native=np.empty((3,65,2*nr),dtype=np.longdouble);coarse=np.empty((65,2*nr),dtype=np.longdouble)
    for lo in range(0,2*nr,128):
        hi=min(lo+128,2*nr);a=coefficients(half[:,lo:hi]);b=coefficients(half[::2,lo:hi])
        for j,n in enumerate(DEGREES):native[j,:,lo:hi]=L[:,:n]@a[:n]
        coarse[:,lo:hi]=L@b[:2048]
    del half
    lift=np.array([p.speed*(m.nativeV[0]@(m.H@(m.heat.faces(t)[0]*m.original.radiation.geometry.tc))) for t in TIMES])
    flux,energy=m.heat.faces(1.);balance=float(abs(np.sum(-np.diff(energy),dtype=np.longdouble))/max(abs(energy).max(),1e-100))
    w=m.original.weights;masks=m.original.masks
    def observe(values):
        cv=values[:,:nr]-lift;cf=values[:,nr:];cv[0]=0;cf[0]=0
        fields=np.column_stack([np.sqrt(np.sum(w*cv*cv,axis=1)),np.sqrt(np.sum(w*cf*cf,axis=1)),
            *[np.sqrt(np.sum(w[mask]*cv[:,mask]**2,axis=1)/w[mask].sum()) for mask in masks]])
        return fields,cv,cf
    rows=[]
    for j,n in enumerate(DEGREES):
        values,cv,cf=observe(native[j]);q=end[j]-m.H@energy;qt=endv[j]-m.H@(flux*m.original.radiation.geometry.tc)
        history=[dict(tau=float(t),**dict(zip(task.FIELDS,values[i].astype(float)))) for i,t in enumerate(TIMES)]
        np.savez_compressed(OUT/f'{label}-{n}.npz',grid=m.grid,cells=m.cells,indices=m.indices,q=q,qt=qt,native_radius=m.original.native,
            native_velocity=cv,native_scalar=cf,weights=w,masks=masks,heat_energy=energy,heat_flux=flux)
        row=dict(classification='Counterexample candidate',degree=degree,coefficients=n,history=history,heat_telescoping=balance,linear_residual=p.error)
        write(OUT/f'{label}-{n}.json',row);rows.append(row)
    cmp=common.comparison(rows);coarse_values,_,_=observe(coarse)
    for i,key in enumerate(task.FIELDS):
        fine=np.array([r[key] for r in rows[-1]['history']]);cmp[key]['contour']=float(max(abs(fine-coarse_values[:,i]))/max(abs(fine).max(),1e-100))
    passed=all(v['last']<.02 and v['order']>1.5 and v['contour']<.0002 for v in cmp.values())
    result=dict(classification='Counterexample candidate',actual_same_input_evolved=True,propagation_passed=passed,
        comparisons=cmp,seconds=time.monotonic()-started,linear_residual=p.error,heat_balance=balance,
        memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,original_failure_resolved=False,full_dynamic_charge_solved=False)
    write(OUT/f'{label}-result.json',result);print('FULL WEEKS RESULT',json.dumps(result),flush=True);return result


def run():
    assert not (OUT/'p4-result.json').exists();plan=json.loads((OUT/'plan.json').read_text())
    for p,h in plan['bindings'].items():assert task.digest(Path(p))==h,p
    assert json.loads((OUT/'pilot-budget.json').read_text())['first_case_forecast_seconds']<450
    signal.alarm(450);resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)));solve(reuse=True);signal.alarm(0)


common=prior.common
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
