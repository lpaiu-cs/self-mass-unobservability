"""Reuse saved resolvents: de Hoog acceleration, no new GR solves.

Vectorized Q-D/continued-fraction algorithm from installed mpmath 1.3.0,
de Hoog, Knight & Stokes (1982), DOI 10.1137/0903022. Binary64 arithmetic
is independently checked; it is not arbitrary-precision mpmath.
"""
import json
import signal
import time
from pathlib import Path
import numpy as np
import def_gr_laplace as base

OUT=base.OUT


def invert(fp,times=base.TIMES,sigma=6.,period=4.):
    M=(len(fp)-1)//2
    scale=np.max(abs(fp),axis=0);assert np.all(scale>0)
    fp=fp/scale
    e=np.zeros_like(fp);q=fp[1:]/fp[:-1];q[0]*=2
    d=[fp[0]/2]
    for r in range(1,M+1):
        count=2*(M-r)+1
        ee=q[1:count+1]-q[:count]+e[1:count+1]
        d.extend([-q[0].copy(),-ee[0].copy()])
        if r<M:q=q[1:count]*ee[1:count]/ee[:count-1]
        e=ee
    d=np.array(d);assert np.all(np.isfinite(d))
    z=np.exp(2j*np.pi*times[:,None]/period)
    a0=np.zeros((len(times),fp.shape[1]),complex);a1=np.broadcast_to(d[0],a0.shape).copy()
    b0=np.ones_like(a0);b1=b0.copy()
    for j in range(1,2*M):
        aa=a1+d[j]*z*a0;bb=b1+d[j]*z*b0
        normalization=np.maximum.reduce([abs(aa),abs(bb),abs(a1),abs(b1)])
        assert np.all(normalization>0)
        a0=a1/normalization;a1=aa/normalization;b0=b1/normalization;b1=bb/normalization
    brem=(1+(d[-2]-d[-1])*z)/2
    # Algebraically equivalent to brem*(sqrt(1+d[-1]*z/brem)-1).
    rem=d[-1]*z/(np.sqrt(1+d[-1]*z/brem)+1)
    value=(a1+rem*a0)/(b1+rem*b0)
    out=2/period*np.exp(sigma*times[:,None])*value.real*scale
    assert np.all(np.isfinite(out))
    return out


def control():
    errors={}
    for omega in [3.,50.,300.,600.]:
        z=6+2j*np.pi*np.arange(513)/4
        fp=1/(z*(z*z+omega*omega))
        exact=(1-np.cos(omega*base.TIMES))/omega**2
        row=[]
        for n in [128,256,512]:
            actual=invert(fp[:n+1,None])[:,0]
            row.append(float(max(abs(actual[1:]-exact[1:]))*omega**2))
        errors[str(omega)]=row
    # Known slow oscillator must pass; unresolved high-frequency controls are
    # retained, not promoted to an automatic uniform oscillatory guarantee.
    assert max(errors['3.0'])<1e-7
    return dict(classification='Counterexample candidate',errors=errors)


def main():
    assert not (OUT/'dehoog-result.json').exists()
    base.write(OUT/'dehoog-plan.json',dict(classification='Counterexample candidate',
        claim='Separate unaccelerated Fourier-tail error from actual old-interface motion using the same513 saved GR resolvents.',
        method='de Hoog Q-D continued fraction, cutoffs128/256/512; all65 same times, same four outputs and thresholds. No added poles, filters or changed initial state.',
        budget=dict(new_GR_solves=0,new_EOS_calls=0,hard_seconds=180),
        source=base.prior.digest(Path(__file__)),input=base.prior.digest(OUT/'fine-transform.npy'),
        acceptance='Same four component relative2% and order1.5. Failure stops; coefficient/outer/alias contrasts remain unrun.'))
    signal.alarm(180);start=time.monotonic();base.write(OUT/'dehoog-control.json',control())
    p=base.Problem(base.prior.OUT/'fine-bank.npz');z=np.load(OUT/'fine-transform.npy').reshape(513,-1,4)
    s=6+2j*np.pi*np.arange(513)/4
    cv=np.array([np.interp(p.native,p.bg.grid,row[:,0]*p.speed) for row in z])*s[:,None]
    cf=np.array([np.interp(p.native,p.bg.grid,row[:,2]) for row in z])
    spectral=np.c_[cv,cf];N=len(p.native)
    lift=np.array([np.interp(p.native,p.bg.grid,p.heat.lift(t,p.bg.nodes,True).reshape(-1,4)[:,0]*p.speed) for t in base.TIMES])
    cases={}
    for n in [128,256,512]:
        state=invert(spectral[:n+1]);state[0]=0
        velocity=state[:,:N]-lift;scalar=state[:,N:]
        readouts=np.column_stack([np.sqrt((velocity*velocity)@p.weights),np.sqrt((scalar*scalar)@p.weights)]+
            [np.sqrt((velocity[:,m]**2)@p.weights[m]/p.weights[m].sum()) for m in p.masks])
        cases[n]=readouts
        np.savez_compressed(OUT/f'dehoog-{n}.npz',times=base.TIMES,velocity_m_s=velocity,Eulerian_scalar=scalar,readouts=readouts)
        print('DEHOOG',n,'SECONDS',time.monotonic()-start,'END',readouts[-1].tolist(),flush=True)
    comparisons={}
    for j,field in enumerate(base.FIELDS):
        a,b,c=[cases[n][:,j] for n in [128,256,512]];norm=max(abs(c).max(),1e-100)
        e1=float(max(abs(a-b))/norm);e2=float(max(abs(b-c))/norm)
        comparisons[field]=dict(previous=e1,last=e2,order=float(np.log2(e1/e2)))
    passed=all(v['last']<.02 and v['order']>1.5 for v in comparisons.values())
    result=dict(classification='Counterexample candidate',cutoff_gates_passed=passed,comparisons=comparisons,
        seconds=time.monotonic()-start,new_GR_solves=0,original_failure_resolved=False,
        full_dynamic_charge_solved=False,remaining='Coefficient, outer, alias and causal-tail checks still required before adopting this propagator.')
    base.write(OUT/'dehoog-result.json',result);signal.alarm(0);print('RESULT',json.dumps(result),flush=True)


if __name__=='__main__':main()
