"""Counterexample candidate: coupled Radau stages on the unchanged source.

Reuse the actual photon/thermal/H and physical ledger owners. This bounded
prefix tests the time bottleneck; no full path or charge is authorized here.
"""
from pathlib import Path
from types import SimpleNamespace
import gc,json,resource,shutil,sys,time
import numpy as np
import sympy as sp
from scipy import sparse
from scipy.sparse.linalg import splu,LinearOperator
import apply_native_conservative_redshift as prior

OUT=Path('native-radau-transfer168-work');OLD=prior.OUT
read,write,sha=prior.read,prior.write,prior.sha
AMP=prior.reuse.AMP;LD=np.longdouble
RK_A=np.array([[5/12,-1/12],[3/4,1/4]])
RK_B=np.array([3/4,1/4]);RK_C=np.array([1/3,1.])
CAPS=dict(prepare=30,pilot=240,finish=30);TOTAL=300
CAPS['profile']=20


def stages(m,t,h,x,g,lus):
    """Solve both collocation stages as one linear physical system."""
    cs=[];sources=[];ls=[];errors=[]
    for fraction in RK_C:
        now=t+fraction*h;c=m.local(now);s,l,e=m.source(now)
        cs.append(c);sources.append(s/(m.scale*AMP));ls.append(l/AMP);errors.append(e)
    # Canonical coefficient knots are step boundaries. The closing stage
    # takes the mechanical derivative from this same interval, as before.
    mechanical=cs[0]['mechanical'].copy()
    for c in cs:c['mechanical']=mechanical
    def stream(xx):return (m.A@xx.reshape(m.n*m.q,m.nf)).reshape(xx.shape)
    def L(j,v):
        xx,gg=m.unpack(v);p,q,*_=m.collision(cs[j],xx,gg)
        return m.pack(stream(xx)+p,q)
    v=m.pack(x,g);dim=len(v)
    q=np.array([m.pack(s+c['q'],m.gas(c['q'],c['qb'],c['qe'])+mechanical) for s,c in zip(sources,cs)])
    rhs=np.tile(v,(2,1))+h*(RK_A@q)
    def mat(value):
        u=value.reshape(2,dim);rates=np.array([L(j,u[j]) for j in range(2)])
        return (u-h*(RK_A@rates)).ravel()
    inverses=[m.inverse(c,h*RK_A[j,j]) for j,c in enumerate(cs)]
    def pre(value):
        result=[]
        for j,vv in enumerate(value.reshape(2,dim)):
            xx,gg=m.unpack(vv);xx=lus[j].solve(xx.reshape(m.n*m.q,m.nf)).reshape(xx.shape)
            result.append(inverses[j](m.pack(xx,gg)))
        return np.array(result).ravel()
    def unpack(value):
        pairs=[m.unpack(vv) for vv in value.reshape(2,dim)]
        return tuple(np.concatenate([p[j] for p in pairs]) for j in range(2))
    moment_owner=SimpleNamespace(unpack=unpack,Nweight=np.tile(m.Nweight,(2,1,1)),
        Eweight=np.tile(m.Eweight,(2,1,1)),eu=np.tile(m.eu,2),nu=np.tile(m.nu,2))
    op=LinearOperator((2*dim,)*2,mat,dtype=float);P=LinearOperator(op.shape,pre,dtype=float)
    iterations=[]
    sol,info=prior.owner.moment_gmres(op,rhs.ravel(),owner=moment_owner,x0=np.tile(v,2),M=P,
        rtol=1e-14,atol=0.,restart=20,maxiter=5,callback=iterations.append,callback_type='pr_norm')
    residual=float(np.linalg.norm(mat(sol)-rhs.ravel())/max(np.linalg.norm(rhs),1e-290))
    m.max_residual=max(m.max_residual,residual);m.max_iterations=max(m.max_iterations,len(iterations))
    assert info==0 and residual<1e-12,('Coupled Radau stage',info,residual,len(iterations))
    result=[]
    for j,vv in enumerate(sol.reshape(2,dim)):
        xx,gg=m.unpack(vv);p,q,e,_=m.collision(cs[j],xx,gg,True)
        result.append((xx,gg,p,q,e,ls[j],errors[j]))
    return result,mechanical


def initialize():
    prior.OUT=OUT;prior.METRIC=OUT/'metric';prior.install();prior.initialize(1)
    model=prior.c.Response;runner=(OUT/'sweep-1/expanded-moment-run.py').read_text()
    def change(a,b):
        nonlocal runner
        assert runner.count(a)==1,(a,runner.count(a));runner=runner.replace(a,b)
    change('gamma=old.prior.GAMMA','gamma=1/4')
    change("lu=splu(sparse.eye(self.n*self.q,format='csc')-gamma*h*self.A)",
           "lus=[splu(sparse.eye(self.n*self.q,format='csc')-v*h*self.A) for v in [5/12,1/4]]")
    a=runner.index('        t=k*h;c=');b=runner.index('        transfer+=',a)
    runner=runner[:a]+'''        t=k*h
        pair,mechanical=stages(self,t,h,x,g,lus)
        (y,gy,p1,g1,es1,l1,e1),(z,gz,p2,g2,es2,l2,e2)=pair
'''+runner[b:]
    change('self.boundary_ports(t+gamma*h,y)','self.boundary_ports(t+h/3,y)')
    change('accepted_angular_luminosity=self.angular)\n    row=',
           "accepted_angular_luminosity=self.angular,accepted_angular_quadrature_weights=h*np.tile(RK_B,count),time_integrator='RadauIIA2',RK_A=RK_A,RK_b=RK_B,RK_c=RK_C)\n    row=")
    namespace=dict(model.run.__globals__,stages=stages,RK_A=RK_A,RK_B=RK_B,RK_C=RK_C)
    exec(compile(runner,__file__,'exec'),namespace);model.run=namespace['run']
    (OUT/'sweep-1/expanded-radau-run.py').write_text(runner)


def prepare():
    assert not OUT.exists();assert not read(OLD/'result.json')['passed'];OUT.mkdir()
    files=[Path(__file__),Path(prior.__file__),OLD/'result.json',OLD/'primitive-check.json',OLD/'plan.json']
    for s in [0,1]:
        for folder in ['photons','material']:(OUT/f'sweep-{s}/{folder}').mkdir(parents=True)
    for p in (OLD/'sweep-0').rglob('*.npz'):
        shutil.copyfile(p,OUT/p.relative_to(OLD));files.append(p)
    for name in ['normalization.json','photon-conservation-plan.json','lift-plan.json','primitive-plan.json']:
        shutil.copyfile(OLD/name,OUT/name);files.append(OLD/name)
    for p,h in read(OLD/'plan.json')['bindings'].items():
        target=OLD/'registered-producer.py' if Path(p).name==Path(prior.__file__).name else Path(p)
        assert sha(target)==h,target
    # Collocation coefficients, order conditions, and stiff accuracy.
    z=sp.symbols('z');cc=sp.Matrix([sp.Rational(1,3),1]);bb=sp.Matrix([[sp.Rational(3,4),sp.Rational(1,4)]])
    aa=sp.Matrix([[sp.Rational(5,12),-sp.Rational(1,12)],[sp.Rational(3,4),sp.Rational(1,4)]])
    ell=[(z-1)/(cc[0]-1),(z-cc[0])/(1-cc[0])]
    assert aa==sp.Matrix([[sp.integrate(p,(z,0,c)) for p in ell] for c in cc])
    assert aa*sp.ones(2,1)==cc and aa[1,:]==bb
    assert (bb*sp.ones(2,1))[0]==1 and (bb*cc)[0]==sp.Rational(1,2)
    assert (bb*cc.applyfunc(lambda v:v*v))[0]==sp.Rational(1,3) and (bb*aa*cc)[0]==sp.Rational(1,6)
    R=sp.factor(1+z*(bb*(sp.eye(2)-z*aa).inv()*sp.ones(2,1))[0])
    assert sp.limit(R,z,sp.oo)==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,stability_function=str(R),
        scope='Two-stage collocation coefficients, classical order3 conditions and stiff accuracy. These do not certify this stiff physical problem or its time error.'))
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='a06c4fa04',
        claim='Test whether simultaneous two-stage RadauIIA order3 removes the failed photon-to-material time controls at the same corrected source and same64/128 clocks.',
        decision='Only a passed4/8-step equal-horizon pair can justify adapting full-path consumers and a fresh measured production forecast. No full path or GR charge is dispatched here.',
        reuse='Same531cells8angles152frequencies,17-knot background and exact167source primitive. Existing collision, physical moment refinement, lifted-state output, gas/impulse/frequency work and signed angular port owners. No new EOS or background trajectory.',
        method='Both physical stages solved together with actual local coefficients at1/3 and1. Every ledger/transfer/port uses weights3/4 and1/4. Saved angular quadrature weights explicitly accompany samples; legacy SDIRK-weight readers must not consume these paths unmodified.',
        gates=dict(time=.02,energy_H=1e-8,stage=1e-12,physical_stage_moment=1e-13,port=1e-12,velocity=1e-4,mapping=1e-10),
        forecast='The167pair used66.5action seconds including41.3run seconds. A coupled two-stage solve can use more Krylov iterations and live coefficient memory;240s hard prefix cap is an assumption of less than4x action cost, not a measured forecast. Release each clock object before constructing the next. No production admission before actual timing and memory.',
        budget=dict(actions=CAPS,total_action_seconds=TOTAL,CPU_threads=1,virtual_GiB=3),
        stop='Stop on error, any gate or cap. Preserve167failures. No finer clocks, new knots, changed source, relaxed gate or automatic extra scheme.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)}))


def pilot():
    initialize();folder=OUT/'sweep-1/photons';rows=[];then=time.monotonic()
    for n in [64,128]:
        m=prior.c.Response(n);row=m.run(n,f'pilot-{n}',n//16)
        with np.load(folder/f'pilot-{n}.npz') as p:
            tt=p['accepted_angular_times'];w=p['accepted_angular_quadrature_weights']
            h=m.t[-1]/n;expected=(h*(np.arange(n//16)[:,None]+RK_C)).ravel()
            assert tt.shape==w.shape==expected.shape and np.max(abs(tt-expected))<1e-18
            flux=p['accepted_angular_luminosity']@(np.arange(1,8,2)/32)
            error=float(abs(w@flux-p['radial_ports'][-1,1,1])/max(w@abs(flux),1e-290))
            assert error<1e-12,error
        row['actual_angular_quadrature_relative']=error;row['time_integrator']='RadauIIA2'
        write(folder/f'pilot-{n}.json',row);rows.append(row);assert row['passed'],row
        del m;gc.collect()
    data=[np.load(folder/f'pilot-{n}.npz')['moments'][:,[0,1,2,3,5,6]] for n in [64,128]]
    errors=np.max(np.sum(abs(data[0]-data[1]),axis=2),axis=0)/np.maximum(np.max(np.sum(abs(data[1]),axis=2),axis=0),1e-290)
    result=dict(classification='Counterexample candidate',passed=bool(max(errors)<.02),rows=rows,time_comparison=errors.astype(float).tolist(),
        seconds=time.monotonic()-then,full_horizon_completed=False,physical_final_charge_solved=False,full_goal_complete=False,
        remaining='A passing prefix is not a full-path result. Adapt all angular/time consumers and measure a justified production budget before continuing these exact prefixes.')
    write(folder/'pilot.json',result);print(json.dumps(result),flush=True);assert result['passed'],result


def finish():
    folder=OUT/'sweep-1/photons';rows=[read(folder/f'pilot-{n}.json') for n in [64,128]]
    assert all(r['passed'] and r['actual_angular_quadrature_relative']<1e-12 for r in rows)
    data=[np.load(folder/f'pilot-{n}.npz')['moments'][:,[0,1,2,3,5,6]] for n in [64,128]]
    errors=np.max(np.sum(abs(data[0]-data[1]),axis=2),axis=0)/np.maximum(np.max(np.sum(abs(data[1]),axis=2),axis=0),1e-290)
    p=dict(classification='Counterexample candidate',passed=bool(max(errors)<.02),rows=rows,time_comparison=errors.astype(float).tolist(),
        original_time_gate=.02,full_horizon_completed=False,physical_final_charge_solved=False,full_goal_complete=False,
        export_repair='Read both completed physical prefixes without repeating any stage. Convert longdouble comparison scalars to binary64 for JSON only; preserve the original failed export and producer.',
        remaining='A prefix alone is not a full-path result. Full-path angular/time consumers and a measured budget must be checked separately; the original physical gates remain unchanged.')
    write(OUT/'result.json',p);print(json.dumps(p),flush=True)


def profile():
    """Locate the saved error before selecting another physical calculation."""
    folder=OUT/'sweep-1/photons';a=np.load(folder/'pilot-64.npz');b=np.load(folder/'pilot-128.npz')
    earlier=np.load(OLD/'sweep-1/photons/pilot-128.npz');rows=[]
    assert np.array_equal(a['t'],b['t']) and np.array_equal(b['t'],earlier['t'])
    for j,name in zip([0,1,2,3,5,6],['photon_energy','material_energy','neutral_H','impulse','photon_pressure','material_pressure']):
        lo=a['moments'][-1,j];hi=b['moments'][-1,j];diff=abs(lo-hi);ids=np.argsort(diff)[-8:][::-1]
        norm=max(np.sum(abs(hi),dtype=LD),LD('1e-290'));den=max(np.sum(diff,dtype=LD),LD('1e-290'))
        rows.append(dict(component=name,time_difference=float(np.sum(diff,dtype=LD)/norm),
            fine_method_difference=float(np.sum(abs(hi-earlier['moments'][-1,j]),dtype=LD)/norm),
            top8_fraction=float(np.sum(diff[ids],dtype=LD)/den),
            cells=[dict(index=int(k),radius_cm=float(b['radius_E'][k]),error_fraction=float(diff[k]/den),
                        coarse=float(lo[k]),fine=float(hi[k])) for k in ids]))
    write(OUT/'error-profile.json',dict(classification='Counterexample candidate',rows=rows,
        scope='Read-only localization of saved short-prefix differences. Neither a continuum error bound nor a demonstrated unique cause. No new stages.',
        inputs={str(p):sha(p) for p in [folder/'pilot-64.npz',folder/'pilot-128.npz',OLD/'sweep-1/photons/pilot-128.npz']}))
    print(json.dumps([{k:v for k,v in r.items() if k!='cells'} for r in rows]),flush=True)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS
    receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));prior.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():
                target=OUT/'registered-producer.py' if p==str(Path(__file__)) and (OUT/'registered-producer.py').exists() else p
                assert sha(target)==h,p
            assert sum(read(p)['seconds'] for p in OUT.glob('*-receipt.json'))+CAPS[action]<=TOTAL
        globals()[action]()
    except Exception as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(action=action,seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
