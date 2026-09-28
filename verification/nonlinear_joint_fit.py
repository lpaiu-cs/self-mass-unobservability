"""Live, constrained 28-parameter timing + profiled Gaussian red/white noise fit.

Every accepted timing step uses the full engine. Broyden proposals do not certify
stationarity: a fresh full Jacobian check is required separately. No EOS signal
constraint is inferred by this preliminary timing/noise calculation.
"""
import argparse
import json
import os
from pathlib import Path
import sys

sys.path.insert(0,str(Path.home()/'work/nutimo_pilot/request13_deps'))
import numpy as np
from scipy.optimize import minimize
from runtime_remediation import ROOT,OUT,live


class Noise:
    def __init__(self,t,err):
        self.sw=1/err; self.n=len(t)
        phase=2*np.pi*(t-t.min())/np.ptp(t)
        f=np.column_stack([v for k in range(1,31) for v in [np.cos(k*phase),np.sin(k*phase)]])
        self.q,self.r=np.linalg.qr(np.column_stack([self.sw[:,None]*f,self.sw]),mode='reduced')
        self.rf=self.r[:,:60]
        self.offset=self.q.T@self.sw
        self.k=np.repeat(np.arange(1,31),2)
        self.last=np.array([0.,4.])

    def profile(self,res):
        y=self.sw*res; z=self.q.T@y
        rest=max(0.,float(y@y-z@z))
        def evaluate(x,details=False):
            ratio=np.exp(x[0]); gamma=x[1]
            eig,rot=np.linalg.eigh((self.rf*self.k**(-gamma))@self.rf.T)
            eig=np.maximum(eig,0)
            inv=1/(1+ratio*ratio*eig)
            v=rot.T@z; u=rot.T@self.offset
            offset=float(np.dot(u*inv,v)/np.dot(u*inv,u))
            v=v-offset*u
            rss=rest+float(v@(inv*v))
            sigma2=rss/self.n
            assert sigma2>0
            logdet=float(np.log1p(ratio*ratio*eig).sum())
            dev=self.n*(1+np.log(sigma2))+logdet
            if details:
                return dict(deviance=dev,sigma=np.sqrt(sigma2),red_amplitude_us=ratio*np.sqrt(sigma2),
                    gamma=gamma,ratio=ratio,offset=offset,root=np.sqrt(inv),rot=rot)
            return dev
        # Several starts for the initial hyperparameter fit; warm start thereafter.
        starts=[self.last,np.array([-5.,1.]),np.array([3.,6.])]
        fits=[minimize(evaluate,start,method='L-BFGS-B',bounds=[(-16.,16.),(0.,10.)],
                       options={'ftol':1e-11,'maxiter':100}) for start in starts]
        best=min(fits,key=lambda r:r.fun)
        self.last=best.x.copy()
        result=evaluate(best.x,True)
        result['noise_optimizer_success']=bool(best.success)
        white=evaluate(np.array([-100.,0.]),True)
        if white['deviance']<=result['deviance']+1e-7:
            result=white
            result.update(red_amplitude_us=0.,ratio=0.,noise_optimizer_success=True)
        return result

    def whiten(self,x,m):
        qx=self.q.T@x
        correction=m['rot']@((m['root']-1)[:,None]*(m['rot'].T@qx)) if x.ndim==2 else m['rot']@((m['root']-1)*(m['rot'].T@qx))
        return (x+self.q@correction)/m['sigma']


def feasible(p):
    # Positive periods/scales/mass ratio; all three eccentricity disks.
    e=[p[6]**2+p[8]**2,p[12]**2+p[14]**2,p[21]**2+p[22]**2]
    return np.array([*(1-1e-6-x for x in e),p[7],p[11],p[13],p[17],p[18],p[23],p[26],
                     np.pi/2-p[1],np.pi/2+p[1]])


def run(args):
    fit,base,r0=live('nonlinear_'+args.label)
    pars=base['params'][base['fmap']].astype(float)
    scale=base['scales'][base['fmap']].astype(float)
    h=np.asarray(json.loads((ROOT/'request10_external/finite_jacobian_v2_meta.json').read_text())['abs_steps'])
    j=np.load(ROOT/'request10_external/finite_jacobian_v2.npy')*h
    t=base['toas'].astype(float)
    ts=np.sort(t); cut=(ts[np.argmax(np.diff(ts))]+ts[np.argmax(np.diff(ts))+1])/2
    counts=args.slip*(t>cut).astype(float)
    noise=Noise(t,base['errs'])
    z=np.zeros(28)
    if args.resume:
        saved=np.load(OUT/('nonlinear-'+args.resume+'.npz'))
        z=saved['scaled_parameters'].copy()
        j=saved['jacobian_proposal'].copy()
        if args.fresh:
            cols=[]
            for k in range(28):
                col=np.load(OUT/f'gradient-{args.resume}-{k:02d}.npz')
                assert np.array_equal(col['scaled_parameters'],z)
                cols.append(col['dcol'])
            j=np.column_stack(cols)
    if args.start == 'displaced' and not args.resume:
        z[21]=4.; z[22]=-3.; z[27]=2.
    assert np.min(feasible(pars+h*z))>0
    ncall=0
    def residual(zz):
        nonlocal ncall
        p=pars+h*zz
        assert np.min(feasible(p))>0, 'Never pass invalid orbit to engine'
        fit.Set_fitted_parameter_relativeshifts(h*zz/scale)
        fit.Compute_lnposterior(0) # Preserve the fixed turn numbers, no phase wrapping.
        r=fit.Get_time_residuals().copy()
        assert np.isfinite(r).all()
        ncall+=1
        return r+counts*(86400e6/p[4])
    r=residual(z); metric=noise.profile(r)
    history=[]; trust=10.; stagnations=0
    path=OUT/('nonlinear-'+args.label+'.json')
    for iteration in range(args.iterations):
        wj=noise.whiten(noise.sw[:,None]*j,metric)
        wr=noise.whiten(noise.sw*r,metric)
        off=noise.whiten(noise.sw,metric); off/=np.linalg.norm(off)
        wj-=off[:,None]*(off@wj); wr-=off*(off@wr)
        gram=wj.T@wj; score=wj.T@wr
        ridge=1e-7*max(np.max(np.diag(gram)),1.)
        def objective(d):
            return .5*d@gram@d+score@d+.5*ridge*(d@d)
        def gradient(d):
            return gram@d+score+ridge*d
        constraints={'type':'ineq','fun':lambda d:feasible(pars+h*(z+d))}
        proposal=minimize(objective,np.zeros(28),jac=gradient,method='SLSQP',
            bounds=[(-trust,trust)]*28,constraints=constraints,
            options={'maxiter':200,'ftol':1e-9})
        d=proposal.x
        accepted=False
        previous=metric['deviance']
        for factor in [1.,.5,.25,.125]:
            dz=factor*d
            if np.min(feasible(pars+h*(z+dz)))<=0: continue
            rr=residual(z+dz); mm=noise.profile(rr)
            if mm['deviance'] < previous-1e-6:
                denom=float(dz@dz)
                if denom>1e-20:
                    # ponytail: Broyden saves full live Jacobians; fresh derivatives
                    # and a stationarity check are required before promotion.
                    j+=np.outer(rr-r-j@dz,dz)/denom
                z+=dz; r=rr; metric=mm; accepted=True
                trust=min(100.,trust*1.5)
                break
        if not accepted:
            trust*=.25; stagnations+=1
        else: stagnations=0
        row=dict(iteration=iteration,evaluations=ncall,accepted=accepted,trust=trust,
            deviance=metric['deviance'],sigma=metric['sigma'],red_amplitude_us=metric['red_amplitude_us'],
            gamma=metric['gamma'],noise_success=metric['noise_optimizer_success'],
            proposal_success=bool(proposal.success),eccentricity_extra=float(np.hypot(*(pars+h*z)[21:23])),
            maximum_step=float(np.max(abs(d))),scaled_parameters=z.tolist())
        history.append(row)
        data=dict(status='Counterexample candidate',label=args.label,pulse_step=args.slip,start=args.start,resumed_from=args.resume,
            covariance='sigma^2 diag(TOA_error^2) + a^2 F diag(k^-gamma) F^T; 30 Fourier pairs',
            linear_nuisance='profiled constant offset',timing='all 28 live nonlinear parameters',
            stationarity_certified=False,physical_signal_matched=False,history=history)
        path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
        np.savez(OUT/('nonlinear-'+args.label+'.npz'),res=r,scaled_parameters=z,parameters=pars+h*z,
                 pulse_counts=counts,jacobian_proposal=j)
        print('ITER',row,flush=True)
        if stagnations>=3 or trust<1e-3: break
    print('DONE',args.label,'candidate fit; stationarity not certified',flush=True)


def check():
    rng=np.random.default_rng(1313); t=np.linspace(0,100,128); err=np.linspace(.8,1.2,128)
    n=Noise(t,err); y=rng.normal(size=128)+3*np.cos(2*np.pi*t/100)
    m=n.profile(y)
    assert m['red_amplitude_us']>0, 'Exercise the correlated-noise path'
    w=n.whiten(n.sw*y,m); off=n.whiten(n.sw,m)
    w-=off*(off@w)/(off@off)
    assert abs(w@w-128)<1e-7
    # Dense covariance independently checks Woodbury whitening and determinant.
    phase=2*np.pi*(t-t.min())/np.ptp(t)
    f=np.column_stack([v for k in range(1,31) for v in [np.cos(k*phase),np.sin(k*phase)]])
    cov=m['sigma']**2*np.diag(err**2)+(f*n.k**(-m['gamma']))@f.T*m['red_amplitude_us']**2
    rr=y-m['offset']; direct=rr@np.linalg.solve(cov,rr)
    assert np.isclose(direct,w@w,rtol=1e-9)
    assert np.isclose(np.linalg.slogdet(cov)[1]-2*np.log(err).sum()+direct,m['deviance'],rtol=1e-9)
    print('PASS: dense covariance/profile/likelihood crosscheck')


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--check',action='store_true')
    p.add_argument('--label',default='zero'); p.add_argument('--slip',type=int,default=0)
    p.add_argument('--start',choices=['baseline','displaced'],default='baseline')
    p.add_argument('--iterations',type=int,default=16)
    p.add_argument('--resume')
    p.add_argument('--fresh',action='store_true')
    a=p.parse_args(); check() if a.check else run(a)
