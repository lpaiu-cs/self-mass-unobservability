"""Keep conserved recovery, spline evaluation and packet differences extended."""
from pathlib import Path
from types import FunctionType,MethodType
import inspect,json,signal,sys,textwrap,time
import numpy as np
import def_native_finite_collision_remainder as old

run=old.run;LD=np.longdouble;OUT=old.OUT/'extended';write=old.write;sha=old.sha


class Spline:
    def __init__(self,s):
        self.original=s;self.tx,self.ty=[np.asarray(v,LD) for v in s.get_knots()]
        self.c=np.asarray(s.get_coeffs(),LD).reshape(len(self.tx)-4,len(self.ty)-4)

    def ev(self,x,y,dx=0,dy=0):
        x,y=np.broadcast_arrays(np.asarray(x,LD),np.asarray(y,LD));shape=x.shape;x=x.ravel();y=y.ravel()
        tx,ty=self.tx,self.ty;c=self.c;kx=ky=3
        for _ in range(dx):
            c=kx*np.diff(c,axis=0)/(tx[kx+1:-1]-tx[1:-kx-1])[:,None];tx=tx[1:-1];kx-=1
        for _ in range(dy):
            c=ky*np.diff(c,axis=1)/(ty[ky+1:-1]-ty[1:-ky-1])[None,:];ty=ty[1:-1];ky-=1
        ix=np.clip(np.searchsorted(tx,x,side='right')-1,kx,len(tx)-kx-2)
        iy=np.clip(np.searchsorted(ty,y,side='right')-1,ky,len(ty)-ky-2)
        a=c[(ix[:,None,None]-kx+np.arange(kx+1)[None,:,None]),(iy[:,None,None]-ky+np.arange(ky+1)[None,None,:])].copy()
        for r in range(1,kx+1):
            for j in range(kx,r-1,-1):
                lo=ix-kx+j;alpha=(x-tx[lo])/(tx[ix+j-r+1]-tx[lo]);a[:,j]=(1-alpha[:,None])*a[:,j-1]+alpha[:,None]*a[:,j]
        a=a[:,kx]
        for r in range(1,ky+1):
            for j in range(ky,r-1,-1):
                lo=iy-ky+j;alpha=(y-ty[lo])/(ty[iy+j-r+1]-ty[lo]);a[:,j]=(1-alpha)*a[:,j-1]+alpha*a[:,j]
        return a[:,ky].reshape(shape)


def promote(method,edits):
    src=textwrap.dedent(inspect.getsource(method.__func__))
    for a,b in edits:
        assert a in src,a;src=src.replace(a,b)
    ns=dict(method.__func__.__globals__,LD=LD);exec(compile(src,__file__,'exec'),ns)
    return MethodType(ns[method.__name__],method.__self__)


class State(old.State):
    def __init__(self):
        super().__init__();m=self.m;e=m.model.flow.eos;s=m.model.spectrum
        self.splines=[];self.owners=[]
        for model,attributes in [(e.warm,['f']),(e.cold,['f']),(s,['f','rev']),(s.cold,['f','rev'])]:
            for key in attributes:
                groups=[]
                for group in getattr(model,key):
                    row=[Spline(v) for v in group];self.splines.extend(row);groups.append(row)
                self.owners.append((model,key,getattr(model,key),groups));setattr(model,key,groups)
        originals=[(e,'evaluate',e.evaluate),(s,'levels',s.levels),(s,'cross',s.cross),(m,'coefficients',m.coefficients),(m,'scattering_matrix',m.scattering_matrix)]
        e.evaluate=promote(e.evaluate,[('np.zeros((7,len(rho)))','np.zeros((7,len(rho)),dtype=LD)')])
        s.levels=promote(s.levels,[('np.zeros((len(rho),10))','np.zeros((len(rho),10),dtype=LD)')])
        cross=s.cross
        # The existing native cross-section routine has a binary64 ABI.
        s.cross=lambda energy,level:cross(np.asarray(energy,float),level).astype(LD)
        m.coefficients=promote(m.coefficients,[('np.zeros((n,self.q,self.nf))','np.zeros((n,self.q,self.nf),dtype=LD)'),('np.zeros(n)','np.zeros(n,dtype=LD)')])
        m.scattering_matrix=promote(m.scattering_matrix,[('np.zeros((3,n,q,nf))',"np.zeros((3,n,q,nf),dtype=c['beta'].dtype)")])
        self.owners.extend((obj,key,fn,getattr(obj,key)) for obj,key,fn in originals);self.precision(False)

    def precision(self,extended):
        for obj,key,before,after in self.owners:setattr(obj,key,after if extended else before)

    def coefficients(self,k,z,factor):
        m=self.m;mat=m.material;model=m.model;b=model.bulk;f=model.flow;nb=m.nb;row=mat.point(k)
        eps=LD(old.AMP)*LD(factor);Q=row['Q'].astype(LD)+eps*z.astype(LD)*row['active'][None]
        h=row['h'].astype(LD);h[1:]-=eps*np.cumsum(z[0,:nb],dtype=LD)
        u,m.theta,m.eta=model.recover_material([h,Q[1,:nb]/LD(run.C),Q[2,:nb],Q[3,:nb]],mat.t[k],row['theta'].astype(LD))
        U=Q[:,nb:]/(mat.V[nb:].astype(LD)*LD(f.eos.rho0)*np.array([1,LD(run.C)**2,LD(run.C)**2,f.eos.nH],LD)[:,None])
        f.seed=row['seed'].astype(LD);m.rho,m.beta,m.lt,m.y=f.primitive(U)
        m.bulk_beta=model.velocity();m.bulk_x=b.eos.x.copy();m.active=m.rho>=f.eos.floor
        return m.coefficients()


source=inspect.getsource(old.State.point)
source=textwrap.dedent(source).replace('c=run.base.Response.local(m,t);',
    "self.precision(False);m.model.flow.seed=np.asarray(m.model.flow.seed,float);target=m.point(k);point=m.point;m.point=lambda index:target\n    try:c=run.base.Response.local(m,t)\n    finally:m.point=point\n    ")
source=source.replace('zero=np.zeros_like(z);','self.precision(True);zero=np.zeros_like(z);')
source=source.replace("rounding=16*np.finfo(float).eps", "rounding=16*np.finfo(LD).eps")
source=source.replace("remainder=np.asarray(dp-factor*linear,float);br=np.asarray(db-factor*bound,float);er=de-factor*escape",
    "remainder=dp-factor*linear;br=db-factor*bound;er=de-factor*escape")
source=source.replace("number=np.sum((remainder-br)*N,axis=(1,2))+er[0]",
    "before=remainder.copy();elo,ehi=LD(m.E[0]),LD(m.E[-1])\n        for _ in range(2):\n            number=np.sum((remainder-br)*N,axis=(1,2))+er[0]\n            remainder[:,0,0]-=number*ehi/(ehi-elo)/N[:,0,0]\n            remainder[:,0,-1]+=number*elo/(ehi-elo)/N[:,0,-1]\n        projection=float(np.sum(abs(remainder-before)*weight)/norm)\n        number=np.sum((remainder-br)*N,axis=(1,2))+er[0]")
source=source.replace("number_relative=float(np.max(abs(number)/ns)),", "number_relative=float(np.max(abs(number)/ns)),projection_over_linear=projection,")
source=source.replace("and max(r['number_relative'] for r in rows)<1e-10", "and max(r['number_relative'] for r in rows)<1e-10 and max(r['projection_over_linear'] for r in rows)<1e-8")
scope=dict(vars(old),OUT=OUT);exec(compile(source,__file__,'exec'),scope);State.point=scope['point']


def prepare():
    assert not OUT.exists();OUT.mkdir();p=json.loads((old.OUT/'pilot.json').read_text());assert not p['eligible']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        failure='The finite collision remainder was below the binary64 gross-rate arithmetic diagnostic and its scattering number cancellation failed. Preserve that pilot; it is not a measured physical remainder.',
        repair='Reuse the same spline knots/coefficients with longdouble de Boor evaluation, carry conserved recovery and collision coefficients in extended precision, and restore the exact known scattering-number invariant at zero energy/momentum with a bounded two-frequency projection. Evaluate only the one needed operator at a canonical knot.',
        native_ABI='Native photoionization cross sections still use their existing binary64 ABI; no microscopic EOS or uniform cross-section derivative certificate is claimed.',
        budgets=dict(repair_pilot_seconds=60,production_seconds=150,prior_pilot_charged_seconds=p['seconds'],original_total_seconds=245),
        reallocation='Original pilot+production allowance65+180=245s. Charge the failed pilot and use60s repair pilot plus150s production; do not enlarge the aggregate allowance.',
        gates=dict(base_owner=1e-9,number=1e-10,projection_over_linear=1e-8,rounding_over_linear=.002),
        stop='One bounded repair pilot before dispatch; no repeated trajectory, new native states or threshold change.',
        bindings={str(f):sha(f) for f in [Path(__file__),Path(old.__file__),old.OUT/'pilot.json',old.OUT/'point-8.npz']}))
    (OUT/'expanded-point.py').write_text(source)


def pilot():
    write(OUT/'tangent-plan.json',dict(classification='Counterexample candidate',
        correction='Freeze the accepted binary64 tangent owner while promoting only the finite nonlinear reference. Restore original spline/method objects before each tangent read, then switch to extended objects for finite evaluation. No prior tangent or trajectory is changed.',
        budgets_unchanged=True,bindings={str(p):sha(p) for p in [Path(__file__),OUT/'plan.json',OUT/'first-producer.py']}))
    start=time.monotonic();signal.alarm(60);s=State();errors=[]
    for spl in s.splines:
        x=np.array([(spl.tx[3]+spl.tx[-4])/2]);y=np.array([(spl.ty[3]+spl.ty[-4])/2])
        for dx,dy in [(0,0),(1,0),(0,1)]:
            before=spl.original.ev(np.asarray(x,float),np.asarray(y,float),dx=dx,dy=dy);after=spl.ev(x,y,dx,dy)
            errors.append(float(np.max(abs(before-after))/max(np.max(abs(before)),1.)))
    assert max(errors)<1e-10
    rows=[]
    for k in [8,16]:
        rows.append(s.point(k))
        if not rows[-1]['passed']:break
    forecast=2*max(r['seconds'] for r in rows)*14+10
    p=dict(classification='Counterexample candidate',rows=rows,spline_owner_relative=max(errors),forecast_upper_seconds=forecast,
        eligible=len(rows)==2 and all(r['passed'] for r in rows) and forecast<150,seconds=time.monotonic()-start)
    write(OUT/'pilot.json',p);signal.alarm(0);print(json.dumps(p),flush=True)


def complete_pilot():
    before=json.loads((OUT/'pilot-failure.json').read_text());spent=json.loads((old.OUT/'pilot.json').read_text())['seconds']+before['seconds']
    assert spent+45<245
    write(OUT/'completion-plan.json',dict(classification='Counterexample candidate',
        failure=before,repair='Restore the binary64 primitive seed before returning to the original tangent and its native cross-section ABI. The successful midpoint is reused; complete only the missing endpoint.',
        pilot_seconds=45,prior_charged_seconds=spent,remaining_original_total_seconds=245-spent,
        bindings={str(p):sha(p) for p in [Path(__file__),OUT/'second-producer.py',OUT/'pilot-failure.json',OUT/'point-8.json',OUT/'point-8.npz']}))
    started=time.monotonic();signal.alarm(45);s=State();row=s.point(16)
    rows=[json.loads((OUT/'point-8.json').read_text()),row];forecast=2*max(r['seconds'] for r in rows)*14+10
    p=dict(classification='Counterexample candidate',rows=rows,forecast_upper_seconds=forecast,
        eligible=all(r['passed'] for r in rows) and forecast<245-spent-(time.monotonic()-started),
        original_remaining_seconds=245-spent-(time.monotonic()-started),seconds=time.monotonic()-started,midpoint_reused=True)
    write(OUT/'completed-pilot.json',p);signal.alarm(0);print(json.dumps(p),flush=True)


if __name__=='__main__':
    signal.signal(signal.SIGALRM,run.native.forcing.history.flow.old.optical.timeout)
    action=sys.argv[1];start=time.monotonic()
    try:globals()[action]()
    except Exception as exc:write(OUT/f'{action}-failure.json',dict(error=repr(exc),seconds=time.monotonic()-start));raise
