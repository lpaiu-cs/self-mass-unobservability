"""Counterexample candidate: compensated GR return about saved joint stages.

The lower component solves the difference of the SAME stage equations. Native
branches are selected from both components; this is not a separate zero-state
nonlinear response subsequently added to a charge. Frozen producers stay intact.
"""
from pathlib import Path
from types import FunctionType
import inspect,json,os,resource,sys,time
import numpy as np
import sympy as sp
import return_joint_gr_geometry as prior

OUT=Path('native-return-evolution188-work')
owner=prior.prior.run.owner;joint=owner.joint
LD,AMP=prior.LD,prior.AMP
read,write,sha=prior.read,prior.write,prior.sha
CAPS=dict(prepare=20,coarse=150,fine=250,audit=25)


def choose(kind,hi,lo=None):
    """Select from a high/low pair without first rounding the combined state."""
    if lo is None:lo=[np.zeros_like(v) for v in hi]
    a,da=hi[0],lo[0]
    if kind=='positive':return a+da>0
    if kind=='nonnegative':return a+da>=0
    if kind in ['min','max']:
        d=(a-hi[1])+(da-lo[1]);return d<=0 if kind=='min' else d>=0
    if kind=='same':return np.sign(a+da)*np.sign(hi[1]+lo[1])>0
    if kind=='absmin':
        sa=np.sign(a+da);sb=np.sign(hi[1]+lo[1])
        return (sa*a-sb*hi[1])+(sa*da-sb*lo[1])<=0
    assert kind in ['argmin','argmax']
    mask=hi[1];assert lo[1].shape==mask.shape
    index=np.argmax(mask,axis=0)
    for i in range(len(a)):
        high=np.take_along_axis(a,index[None],axis=0)[0]
        low=np.take_along_axis(da,index[None],axis=0)[0]
        d=(a[i]-high)+(da[i]-low)
        better=(d<0 if kind=='argmin' else d>0)&mask[i]
        index=np.where(better,i,index)
    return index


def paired_tangent():
    """Reuse the live directional owners, recording actual branch operands."""
    engine=owner.engine;events=[];decisions=[];desired=[];position=0;mode='high'
    def select(kind,*args):
        nonlocal position
        args=[np.array(v,copy=True) for v in args]
        if mode=='high':
            events.append((kind,args))
            if position==len(decisions):decisions.append(choose(kind,args))
        elif mode=='low':
            old,hi=events[position];assert old==kind
            desired.append(choose(kind,hi,args))
        answer=decisions[position];position+=1;return answer
    def minimum(a,b):return np.where(select('min',a,b),a,b)
    def maximum(a,b):return np.where(select('max',a,b),a,b)
    def minmod(a,b):return np.where(select('same',a,b),np.where(select('absmin',a,b),a,b),0.)
    direction,_=owner.clone(engine.face.minmod_direction,[
        ('np.minimum(da,db)','minimum(da,db)'),('np.maximum(da,db)','maximum(da,db)'),
        ('(da*b>0)',"select('positive',da*b)"),('(db*a>0)',"select('positive',db*a)")],
        dict(minimum=minimum,maximum=maximum,minmod=minmod,select=select))
    reconstruction,_=owner.clone(engine.face.reconstruction,[('np.maximum(d[0],0)','maximum(d[0],0)')],
        dict(minmod_direction=direction,maximum=maximum))
    def extremum(values,derivatives,minimum):
        values=np.array(values);derivatives=np.array(derivatives)
        best=np.min(values,axis=0) if minimum else np.max(values,axis=0)
        index=select('argmin' if minimum else 'argmax',derivatives,values==best)
        return best,np.take_along_axis(derivatives,index[None],axis=0)[0]
    flux=FunctionType(engine.flux_direction.__code__,dict(engine.flux_direction.__globals__,extremum=extremum))
    deep,_=owner.clone(joint.previous.original.reuse.deep.deep_tangent,
        [('dmass>=0',"select('nonnegative',dmass)")],dict(select=select))
    source=engine.source.replace('flux[0]>=0',"select('nonnegative',flux[0])").replace('m.deep_tangent(k,z,field)','deep(m,k,z,field)')
    ns=dict(engine.tangent.__globals__,reconstruction=reconstruction,flux_direction=flux,deep=deep,select=select)
    exec(compile(source,__file__,'exec'),ns)
    def reset(next_mode):
        nonlocal mode,position
        assert position==len(decisions),(position,len(decisions))
        mode=next_mode;position=0
        if mode=='high':events.clear()
        if mode=='low':desired.clear()
    def settle():
        assert len(desired)==len(decisions)
        count=sum(np.count_nonzero(a!=b) for a,b in zip(desired,decisions))
        # A crossing away from a high-component tie requires a high-component
        # residual correction. Reject it instead of subtracting two large RHSs.
        for new,old,(kind,hi) in zip(desired,decisions,events):
            changed=new!=old
            if not np.any(changed):continue
            a=hi[0]
            if kind in ['positive','nonnegative']:tie=a==0
            elif kind in ['min','max']:tie=a==hi[1]
            elif kind=='same':tie=(a==0)|(hi[1]==0)
            elif kind=='absmin':tie=abs(a)==abs(hi[1])
            else:tie=np.take_along_axis(a,new[None],axis=0)[0]==np.take_along_axis(a,old[None],axis=0)[0]
            assert np.all(tie|~changed),('Nonzero high-branch crossing',kind,int(np.count_nonzero(~tie&changed)))
        if count:decisions[:]=[v.copy() for v in desired]
        return count
    return ns['tangent'],reset,settle


def prepare():
    assert not OUT.exists();OUT.mkdir();assert read(prior.OUT/'compensated-result.json')['passed']
    files=[];reused={}
    for s in [0,1]:
        for part in ['photons','material']:(OUT/f'sweep-{s}/{part}').mkdir(parents=True)
    for src in list((prior.OUT/'sweep-0').rglob('*.npz'))+[prior.OUT/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
        dst=OUT/src.relative_to(prior.OUT);os.link(src,dst);reused[str(dst)]=sha(src);files.append(src)
    files += [prior.prior.saved(n) for n in [64,128]]+[prior.OUT/'metric/metric-128-g8.npz',prior.OUT/'compensated-result.json']
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='25244662e',
        claim='Evolve the actual GR return as a compensated increment of the same saved photon/B/S/Etilde/H Radau equations, including branch ties, collisions, floor and boundary ledgers.',
        decision='Original stage, balance, branch and paired-time gates decide whether this return representation may be extended. Failure stops before full-horizon return or charge reading.',
        scope='Original64/128 clocks over the same firstT/64 interval:2/4 existing front-split steps. Fixed187metric from the same185base is one feedback iterate, not a converged self-GR solution.',
        base='Decode the actual stored Eref/B/S/H stage state; require its native RHS to reproduce the saved native rates. Base solver residual and state-storage error remain separate from the new lower-component residual. Tiny feedback is not thereby certified against total discretization error.',
        arithmetic='Evaluate branch operands as high/low pairs. Permit switches at exact high ties, reject crossings away from high ties. Reuse actual native directional owners and current odd photon source. Never round high+low first or evolve an independent zero-state native response.',
        gates=dict(stage=1e-12,physical_stage=1e-13,base_native=1e-12,constitutive=.002,conservation=1e-8,port=1e-12,time=.02,max_Newton=3,max_branch_passes=8),
        budgets=CAPS,CPU_threads=1,virtual_GiB=4,
        forecast='Measured183four/eight-step costs197/246s included setup. Here two/four steps plus paired branch evaluation are expected60..150s /100..250s; this added cost is unmeasured. Hard caps150/250s; no automatic expansion or retries.',
        stop='Any old gate, nonzero branch crossing, branch stabilization, source binding or time budget failure. Keep rejected states. Frozen185sources and plan are untouched.',
        final_charge_conclusion='unadjudicated',full_goal_complete=False,
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},reused=reused))
    x,d,a,b=sp.symbols('x d a b');assert sp.expand((a*(x+d)+b)-(a*x+b)-a*d)==0
    tiny=LD('1e-30');h=np.array([1.,1.],LD);l=np.array([tiny,0.],LD)
    assert np.all(h+l==h) and not choose('min',[h[:1],h[1:]],[l[:1],l[1:]])[0]
    assert choose('argmin',[h[:,None],np.ones((2,1),bool)],[l[:,None],np.ones((2,1),bool)])[0]==1
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        scope='Affine difference identity within a shared branch and high/low exact-tie selection only; no global native linearity or full-GR theorem.'))


def initialize():
    global Model
    prior.prior.run.prior.OUT=OUT;prior.prior.run.prior.initialize();Parent=owner.Model
    class Increment(Parent):
        def __init__(self,n):
            super().__init__(n);self.primary=self.driver;self.returned=prior.ReturnOnly(self.primary)
            self.select(self.returned);self.anchor=dict(np.load(prior.prior.saved(n)));self.selected={}
            self.anchor_checks=[];self.branch_checks=[];self.checked_times=set()
        def select(self,driver):self.driver=driver;self.redshift_driver=driver;self.material.driver=driver
        def anchor_state(self,t):
            k=int(np.argmin(abs(self.anchor['joint_stage_times']-t)))
            assert abs(self.anchor['joint_stage_times'][k]-t)<1e-18,('Not an actual saved stage',t)
            q=self.anchor['joint_stage_conserved_scaled'][k]
            return k,np.column_stack([(q[2]-self.kappa*q[0])/self.eu,q[3]/self.nu,q[0]/self.bu,q[1]/self.su])
        def native(self,t,g,probe=1.,details=False,tangent=None,metric=True):
            if tangent is not None:
                self.select(self.returned)
                return Parent.native(self,t,g,probe,details,tangent,metric)
            assert metric
            k,base=self.anchor_state(t);fn,reset,settle=paired_tangent();changes=[]
            try:
                for iteration in range(8):
                    if iteration:reset('high')
                    self.select(self.primary);high=Parent.native(self,t,base,probe,tangent=fn)
                    reset('low');self.select(self.returned);low=Parent.native(self,t,g,probe,True,fn)
                    changed=settle();changes.append(changed)
                    if not changed:break
                else:raise AssertionError(('Unsettled paired branches',changes))
                if t not in self.checked_times:
                    saved=self.anchor['joint_native_rates_scaled'][k]/self.units
                    error=(np.sum(abs(high-saved)*self.units,axis=0)/np.maximum(np.sum(abs(saved)*self.units,axis=0),LD('1e-290'))).astype(float)
                    self.anchor_checks.append(dict(time=float(t),native_relative=error.tolist()))
                    assert max(error)<1e-12,('Saved-stage anchor precision',error.tolist())
                    self.checked_times.add(t)
                self.selected[t]=(fn,lambda capture:reset('replay'))
                self.branch_checks.append(dict(time=float(t),iterations=len(changes),tie_switches=changes))
                return low if details else low[0]
            finally:self.select(self.returned)
    source=(OUT/'expanded-affine-jacobian.py').read_text()
    old="base=self.native(t,g);tangent,reset=selected_tangent()\n    selected=self.native(t,g,tangent=tangent)\n    assert np.array_equal(selected,base),'Selected branches must reproduce the original native direction'"
    assert source.count(old)==1
    source=source.replace(old,'base=self.native(t,g);tangent,reset=self.selected[t]')
    ns=dict(Parent.jacobian.__globals__);exec(compile(source,__file__,'exec'),ns);Increment.jacobian=ns['jacobian']
    (OUT/'expanded-increment-jacobian.py').write_text(source);Model=Increment


def evolve(n):
    initialize();m=Model(n);row=m.run(n,f'return-{n}',n//64)
    p=OUT/f'sweep-1/photons/return-{n}.npz'
    audit,_,_,_=prior.prior.run.verify(p,p)
    row.update(audit=audit,anchor_checks=m.anchor_checks,branch_checks=m.branch_checks,
        maximum_true_stage=max(v[-1]['relative'] for v in m.newton_iterations),
        maximum_true_physical_stage=max(max(v[-1]['moments']) for v in m.newton_iterations),
        same_saved_stage_equation=True,actual_return_time_evolved=True,self_GR_return_closed=False,
        final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/f'run-{n}.json',row);assert row['passed']
    write(OUT/f'checks-{n}.json',dict(classification='Counterexample candidate',newton=m.newton_iterations,stages=m.stage_log))


def audit():
    rows=[read(OUT/f'run-{n}.json') for n in [64,128]]
    values=[];clocks=[]
    for n in [64,128]:
        p=OUT/f'sweep-1/photons/return-{n}.npz'
        _,ph,gas,clock=prior.prior.run.verify(p,p);values.append((ph,gas));clocks.append(clock)
    assert np.array_equal(*clocks)
    relative=joint.previous.run.c.relative
    errors=relative(values[0][0],values[1][0])+relative(values[0][1],values[1][1])
    result=dict(classification='Counterexample candidate',passed=max(errors)<.02,
        same_horizon_seconds=float(clocks[0][-1]),time_relative=errors,rows=rows,
        actual_return_time_evolved=True,self_GR_return_closed=False,
        final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'result.json',result);print(json.dumps(result),flush=True);assert result['passed'],errors


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        if action!='prepare':
            for p,h in dict(read(OUT/'plan.json')['bindings'],**read(OUT/'plan.json')['reused']).items():assert sha(p)==h,p
        if action in ['coarse','fine']:evolve(64 if action=='coarse' else 128)
        else:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(action=action,seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
