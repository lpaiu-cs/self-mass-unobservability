"""Bounded native temperature-column repair for conservative energy recovery."""
from pathlib import Path
import json
import signal
import time
import numpy as np
from scipy.interpolate import CubicHermiteSpline
import def_native_energy_release as old

OUT=old.OUT


def main():
    assert not (OUT/'temperature-columns.npz').exists();start=time.monotonic();signal.alarm(160)
    old.write(OUT/'column-reassessment.json',dict(classification='Counterexample candidate',
        reason='The actual first conservative step needs sigma0.799, inside the failed wide entropy interpolant. Temperature recovery avoids that ill-shaped coordinate, but neighboring density columns need real native support and resolved ionization curvature.',
        decision='Reuse1470 native states; pad each temperature column by0.08 in lnT and refine only native midpoint disagreements. No fluid grid or time enlargement. The original wide-entropy audit stays failed.',
        additional_native_calls=6000,seconds=160,maximum_depth=7,column_midpoint_gate=.0005,final_independent_gate=.002,
        stop='Persist each finished column and stop on budget, depth, native domain or final independent control failure. No extrapolation, tolerance change or automatic further run.',
        bindings={str(p.relative_to(old.old.ROOT)):old.old.photons.digest(p) for p in [Path(__file__),OUT/'eos.npz',OUT/'eos.json',OUT/'first-step-native.json']}))
    d=np.load(OUT/'eos.npz');fan=old.task.prior.Fan(call_cap=6000,reuse=True);columns=[];temps=[];offsets=[0];records=[]
    for j,xx in enumerate(d['x']):
        cache={float(np.log(t)):a.copy() for t,a in zip(d['T'][:,j],d['raw'][:,j])}
        def call(t):
            if t not in cache:cache[t]=fan.call(np.log(fan.rho)+xx,t)
            return cache[t]
        knots=sorted(cache);left=knots[0]-.08;right=knots[-1]+.08;call(left);call(right);knots=[left]+knots+[right]
        worst=0.;maxdepth=0
        def refine(l,r,depth):
            nonlocal worst,maxdepth
            a,b=call(l),call(r);mid=(l+r)/2;native=call(mid);rho=native[0]
            lp=CubicHermiteSpline([l,r],np.log([a[1],b[1]]),[a[6],b[6]])
            u=CubicHermiteSpline([l,r],[a[2],b[2]],[a[10],b[10]])
            pp=float(np.exp(lp(mid)));uu=float(u(mid));cv=float(u(mid,1));chi=float(lp(mid,1));cr=(a[5]+b[5])/2
            gm=cr+pp/rho*chi*chi/cv
            err=max(abs(pp/native[1]-1),abs(uu/native[2]-1),abs(cv/native[10]-1),abs(gm/native[4]-1),abs(cr/native[5]-1))
            worst=max(worst,err);maxdepth=max(maxdepth,depth)
            if err>.0005:
                assert depth<7,('Column refinement depth',j,l,r,err)
                refine(l,mid,depth+1);refine(mid,r,depth+1)
        for l,r in zip(knots[:-1],knots[1:]):refine(l,r,0)
        ts=sorted(cache);temps.extend(ts);columns.extend([cache[t] for t in ts]);offsets.append(len(temps))
        records.append(dict(density_index=j,states=len(ts),original_midpoint_maximum=worst,depth=maxdepth))
        np.savez_compressed(OUT/'temperature-columns-progress.npz',x=d['x'],logT=temps,raw=columns,offsets=offsets,native_calls=fan.calls)
        if j==3:
            elapsed=time.monotonic()-start;forecast=elapsed*len(d['x'])/4*1.5
            old.write(OUT/'column-budget.json',dict(first4_columns_seconds=elapsed,forecast_seconds=forecast,assumption='First low-density columns times147/4 with50percent margin; denser ionization regions may be more expensive. Both time and native-call hard caps remain active.'))
    np.savez_compressed(OUT/'temperature-columns.npz',x=d['x'],logT=temps,raw=columns,offsets=offsets,rho0=d['rho0'],cx=d['cx'],s0=d['s0'],sunit=d['sunit'])
    result=dict(classification='Counterexample candidate',native_calls=fan.calls,seconds=time.monotonic()-start,states=len(temps),reused_states=1470,columns=records,independent_controls_pending=True)
    old.write(OUT/'column-result.json',result);signal.alarm(0);print(json.dumps({k:v for k,v in result.items() if k!='columns'}),flush=True)


if __name__=='__main__':main()
