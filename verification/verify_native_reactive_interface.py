"""Check saved physical initialization and both sides of the reactive interface."""
from pathlib import Path
import json
import signal
import time
import numpy as np
import sympy as s
import def_native_reactive_interface as task


def main():
    out=task.PHYSICAL;assert not (out/'audit.json').exists();start=time.monotonic();signal.alarm(20)
    # Telescoping one common face flux makes the transfer cancel before
    # evaluating either EOS. Sources retain their separate photon ledger.
    F0,Fj,FN,Sl,Sr=s.symbols('F0 Fj FN Sl Sr')
    left=F0-Fj+Sl;right=Fj-FN+Sr
    assert s.expand(left+right-(F0-FN+Sl+Sr))==0
    task.write(out/'symbolic.json',dict(classification='Proven',passed=True,
        scope='For each conservative variable, the shared interior face flux cancels between adjacent subdomains. Separate material/photon energy sources remain paired. This does not certify a whole-star physical boundary.'))
    rows=[];checks=[];native=task.ex.Native(cap=74)
    for n,label,boundary in [(448,'cells-448','copy'),(896,'cells-896','copy'),(448,'boundary-448','saved')]:
        f=task.PhysicalFlow(n,boundary);m=f.base;j=f.join;d=np.load(out/f'full-{label}.npz');local=np.load(out/f'{label}.npz')
        assert np.array_equal(d['U'][:,j:],local['U']) and np.array_equal(d['initial'][:,j:],local['initial'])
        assert np.array_equal(d['snapshots'][:,:,j:],local['snapshots']) and np.array_equal(d['volume'][j:],local['volume'])
        actual_t=np.interp(np.minimum(m.R+m.x/m.As,m.R),m.env['r'],np.log(m.env['T']));active=d['initial'][0]>0
        assert np.array_equal(d['initial_logT'][active],actual_t[active])
        assert np.array_equal(d['initial'],f.initial)
        h=d['history'];lh=local['history'];discard=d['conserved_discard']-local['conserved_discard'];volume=d['volume'][:j]
        difference=np.sum((d['U'][:,:j]-d['initial'][:,:j])*volume,axis=1)+discard
        difference[0]-=h[-1,3]-lh[-1,3]
        difference[2]-=h[-1,4]-lh[-1,4]
        difference[3]-=(h[-1,6]+h[-1,7])-(lh[-1,6]+lh[-1,7])
        denom=np.sum(d['initial'][0,:j]*volume)
        normalized=[abs(difference[0])/denom,abs(difference[2])/max(abs(h[-1,2]-lh[-1,2]),1e-100),abs(difference[3])/(denom*f.eos.y0)]
        assert normalized[0]<1e-10 and normalized[1]<1e-8 and normalized[2]<1e-9
        result=json.loads((out/f'{label}.json').read_text())
        rows.append(dict(label=label,saved_physical_logT_exact=True,interface_buffer_balance=list(map(float,normalized)),hydrodynamic_characteristic_distance_cm=result['inner_characteristic_distance_cm'],buffer_width_cm=result['buffer_width_cm']))
        if n==896:
            e=f.eos;rho=d['rho']/e.rho0;lt=d['logT'];y=np.divide(d['U'][3],d['U'][0],out=np.full(f.n,e.y0),where=d['U'][0]>=e.floor);e.y=y
            p,u,g,T,kap,cv,_=e.evaluate(rho,lt);ids=np.flatnonzero(rho>=e.floor)
            chosen=np.unique(np.r_[j//2,j,j+32,ids[np.linspace(0,len(ids)-1,4).astype(int)]])
            for i in chosen:
                x=float(np.log(rho[i]));t=float(lt[i]);yy=float(y[i]);a=native.state(x,t,yy);r=native.rates(a,float(e.d['Trad']),32)
                e.y=np.array([yy]);rr=e.reactions(np.array([rho[i]]),np.array([t]))[0]
                errors=[abs(p[i]*e.rho0*task.C**2/a['raw'][1]-1),abs(u[i]*task.C**2/a['raw'][2]-1)]
                if i in [j//2,j]:
                    step=1e-4;ar=(native.state(x+step,t,yy)['raw']-native.state(x-step,t,yy)['raw'])/(2*step);at=(native.state(x,t+step,yy)['raw']-native.state(x,t-step,yy)['raw'])/(2*step)
                    gamma=ar[1]/a['raw'][1]+at[1]/a['raw'][1]*(a['raw'][1]/a['raw'][0]-ar[2])/at[2]
                    errors.extend([abs(g[i]/gamma-1),abs(cv[i]*task.C**2/at[2]-1)])
                checks.append(dict(cell=int(i),x_cm=float(m.x[i]),T=float(T[i]),y=yy,constitutive=list(map(float,errors)),rates=float(np.max(abs(rr/r-1)))))
    passed=bool(max(max(z['constitutive']) for z in checks)<.002 and max(z['rates'] for z in checks)<.002)
    np.savez_compressed(out/'audit-native-states.npz',**{k:np.array([z[k] for z in native.ion.states]) for k in native.ion.states[0]})
    result=dict(classification='Counterexample candidate',passed=passed,rows=rows,actual_native_controls=checks,EOS_calls=native.ion.calls,seconds=time.monotonic()-start,source_sha256=task.ex.old.cold.sha(__file__),full_stellar_interior=False,full_goal_complete=False)
    task.write(out/'audit.json',result);signal.alarm(0);print(json.dumps(result),flush=True);assert passed


if __name__=='__main__':main()
