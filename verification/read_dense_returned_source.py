"""Keep actual returned metric jumps in a general two-field source readout.

Counterexample candidate. Hermite interpolation is a declared finite time
representation, checked against held-out actual stages, not a uniform bound.
"""
from pathlib import Path
from types import SimpleNamespace
from math import factorial
import gc,json,os,resource,sys,time
import numpy as np
from scipy.interpolate import PPoly
import read_returned_joint_source as prior

CHECK=prior.CHECK;ROOT=Path('native-dense-returned246-work');OUT=ROOT/('check' if CHECK else 'full')
INPUT=prior.OUT;PRIMARY=Path('native-complete-radau224-check-work' if CHECK else 'native-complete-radau224-work')
base=prior.prior;read,write,sha,bind,LD,AMP=prior.read,prior.write,prior.sha,prior.bind,prior.LD,prior.AMP
KEYS=base.KEYS;StageDriver=prior.returned.StageDriver
CAPS=dict(prepare=600,geometry=1800,source=5400)


def prepare():
    assert read(INPUT/'result.json')['representation_controls_passed'];assert not OUT.exists();OUT.mkdir(parents=True)
    files=[]
    for part in ['sweep-0','sweep-1/photons','sweep-1/material','gr','metric']:(OUT/part).mkdir(parents=True)
    srcs=list((INPUT/'sweep-0').rglob('*.npz'))
    srcs += [INPUT/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json','metric/metric-128-g8.npz','endpoint-sources.json']]
    srcs += [INPUT/f'gr/endpoint-{n}.npz' for n in [64,128]]
    for p in srcs:
        dst=OUT/p.relative_to(INPUT);dst.parent.mkdir(parents=True,exist_ok=True);os.link(p,dst);files += [p,dst]
    files += [prior.saved(n) for n in [64,128]]+[prior.INPUT/f'recovered-{n}.npz' for n in [64,128]]
    files += [INPUT/n for n in ['plan.json','result.json','source-receipt.json']]
    files += [PRIMARY/'gr/source-128.npz',base.OUT/'expanded-source.py']
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Carry the actual GR-return low component into a continuous source using its accepted Radau gas/photon states and its own metric, while preserving mass-constraint floor jumps.',
        method='Decompose the same stable readout into cubic zero-geometry state plus two independently computed affine coefficient maps for u and lambda. Hermite metric pieces use actual u/lambda values and propagated u_t/analytic lambda rates. Derive left/right lambda jumps and derivative changes from the SAME primary-source mass constraint. Preserve actual post-floor values at all applied nodes.',
        controls='Exact source times and original endpoint reproduction; source polynomial/dense-stage/mapping1e-12, pressure0.002. Hold out actual fine-stage metric samples from a coarser subset retaining all known floor knots; require spatialL1 metric error below0.002. This sampled interpolation check is not a uniform temporal error certificate. Keep that limitation in every charge verdict.',
        decision='Use the completed short actual227solve to check the adapter, then the admitted236/245same long solve. No independent response, added old charge, repeated material/photon/GR field, finer physical mesh or changed input.',
        budgets=CAPS,CPU_threads=1,CPU_affinity=8,virtual_GiB=16,
        forecast='24544endpoint readouts took56.86s. Existing224350-stage dense source took570s for326steps. Short dense map approximately2..5minutes; full10..25minutes if marginal cost holds. Geometry algebra adds several fixed source maps, not a field solve. Allow30minutes geometry and90minutes source.',
        stop='Any source binding, endpoint/jump, held-out interpolation, stage, pressure, mapping or polynomial gate fails; preserve rejected representation, no automatic refinement or gate relaxation.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},final_charge_conclusion='unadjudicated',full_goal_complete=False))
    import sympy as s
    x,h,a,b,c,d=s.symbols('x h a b c d');p=a+h*c*x+(3*(b-a)-h*(2*c+d))*x*x+(2*(a-b)+h*(c+d))*x**3
    assert s.expand(p.subs(x,0)-a)==s.expand(p.subs(x,1)-b)==0
    assert s.expand(s.diff(p,x).subs(x,0)-h*c)==s.expand(s.diff(p,x).subs(x,1)-h*d)==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        scope='One-sided Hermite value/derivative identities; the right endpoint point value may have a separate jump. No continuum or physical error theorem.'))


def hermite(t,y,right_rate,jump,left_rate,ids):
    h=np.diff(t[ids])[:,None];a=y[ids[:-1]];b=y[ids[1:]]-jump[ids[1:]]
    c=right_rate[ids[:-1]];d=left_rate[ids[1:]];slope=(b-a)/h
    return PPoly(np.array([(c+d-2*slope)/(h*h),(3*slope-2*c-d)/h,c,a]),t[ids])


def geometry():
    bind(base.endpoint.initialize,OUT=OUT)();m=base.base.gr.Response()
    g=dict(np.load(OUT/'metric/metric-128-g8.npz'));t=g['t'];raw=dict(np.load(PRIMARY/'gr/source-128.npz'))
    knots,co=base.coefficients(raw);poly={k:PPoly(np.asarray(v[::-1],float),knots) for k,v in co.items()}
    at=t.copy();near=np.array([np.argmin(abs(raw['t']-v)) for v in t]);is_edge=abs(raw['t'][near]-t)<1e-18
    at[is_edge]=raw['t'][near[is_edge]]
    right={k:p(at) for k,p in poly.items()};dr={k:p(at,nu=1) for k,p in poly.items()};left={};dl={}
    for k,p in poly.items():
        idx=np.clip(np.searchsorted(p.x,at,side='left')-1,0,len(p.x)-2);dt=(at-p.x[idx]).reshape((-1,)+(1,)*(p.c.ndim-2))
        v=np.zeros_like(right[k]);dv=np.zeros_like(v)
        for c in p.c:dv=dv*dt+v;v=v*dt+c[idx]
        left[k]=v;dl[k]=dv;right[k][is_edge]=raw[k][near[is_edge]]
    # The scalar field is continuous; only the instantaneous constraint source
    # contributes to a lambda jump. Its derivative has analogous one-sided data.
    zero=np.zeros((len(t),len(raw['radius'])+1));field=dict(delta_phi=zero,delta_Phi=zero)
    values=[]
    for sources in [right,left,dr,dl]:
        d=dict(raw,t=t,**sources)
        values.append(prior.returned.geometry.constraints.centers(m,d,field,8)['delta_lambda'][:,:-1])
    jump=values[0]-values[1];jump[~is_edge]=0.;jump[0]=0.
    lam_left_rate=g['actual_delta_lambda_rate']+values[3]-values[2]
    ids=np.arange(len(t));u=hermite(t,g['delta_u'],g['delta_u_t'],np.zeros_like(jump),g['delta_u_t'],ids)
    lam=hermite(t,g['delta_lambda'],g['actual_delta_lambda_rate'],jump,lam_left_rate,ids)
    p=np.load(prior.saved(64));coarse_times=np.r_[p['actual_step_edges'],p['joint_stage_times'],raw['t']]
    subset=np.unique([int(np.argmin(abs(t-v))) for v in coarse_times]);assert np.max([np.min(abs(t-v)) for v in coarse_times])<1e-18
    held=np.setdiff1d(ids,subset);assert len(held)>0
    errors={};endpoint={}
    for key,y,rate,jmp,left_rate,full in [('u',g['delta_u'],g['delta_u_t'],np.zeros_like(jump),g['delta_u_t'],u),
        ('lambda',g['delta_lambda'],g['actual_delta_lambda_rate'],jump,lam_left_rate,lam)]:
        coarse=hermite(t,y,rate,jmp,left_rate,subset)
        norm=max(np.max(np.sum(abs(y),axis=-1)),1e-290)
        errors[key]=float(np.max(np.sum(abs(coarse(t[held])-y[held]),axis=-1))/norm)
        endpoint[key]=float(np.max(np.sum(abs(full(t[:-1])-y[:-1]),axis=-1))/norm)
    row=dict(classification='Counterexample candidate',passed=max(errors.values())<.002 and max(endpoint.values())<1e-12,
        actual_metric_samples=len(t),held_out_actual_samples=len(held),held_out_spatialL1=errors,node_relative=endpoint,
        maximum_lambda_jump=float(np.max(abs(jump))),retained_floor_knots=int(np.count_nonzero(is_edge)),
        uniform_temporal_error_certificate=False,final_charge_conclusion='unadjudicated')
    np.savez_compressed(OUT/'geometry.npz',times=t,u_coeff=u.c,lambda_coeff=lam.c,lambda_jump=jump,
        lambda_left_rate=lam_left_rate,coarse_subset=subset,held_out=held)
    write(OUT/'geometry-result.json',row);assert row['passed'],row


class DenseDriver(StageDriver):
    def __init__(self,primary):
        super().__init__(primary);d=np.load(OUT/'geometry.npz')
        self.u=PPoly(d['u_coeff'],d['times']);self.lam=PPoly(d['lambda_coeff'],d['times'])
    def returned(self,t):
        j=int(np.argmin(abs(self.clock-t)))
        if abs(self.clock[j]-t)<1e-18:return StageDriver.returned(self,t)
        assert self.clock[0]<t<self.clock[-1];i=np.clip(np.searchsorted(self.clock,t)-1,0,len(self.clock)-2)
        w=LD((t-self.clock[i])/(self.clock[i+1]-self.clock[i]))
        row={k:(1-w)*v[i].astype(LD)+w*v[i+1] for k,v in self.g.items() if k.startswith('delta_') and v.shape==(len(self.clock),self.n)}
        # Only u and lambda enter this source map. Other interpolated entries
        # support its pressure check; this driver must never evolve new states.
        row.update(delta_u=self.u(t),delta_lambda=self.lam(t),delta_u_t=self.u(t,nu=1),delta_lambda_rate=self.lam(t,nu=1))
        row['delta_log_volume']=3*row['delta_u']+row['delta_lambda'];row['delta_log_areal_radius']=row['delta_u']
        return row


def initialize():
    prior.returned.StageDriver=DenseDriver;bind(prior.initialize,OUT=OUT)()


def coefficients(d):
    knots=np.unique(np.r_[d['t'],d['geometry_times'],d['metric_times']]);knots=knots[(knots>=0)&(knots<=d['t'][-1])]
    left=knots[:-1];gt=d['geometry_times'];idx=np.clip(np.searchsorted(gt,left,side='right')-1,0,len(gt)-2);result={}
    fields={k:PPoly(d[k+'_coeff'],d['metric_times']) for k in ['u','lambda']}
    for key in KEYS:
        p=base.prior.prior.poly(d['t'],d['state_coeff_'+key]);out=np.zeros((5,len(left))+d[key].shape[1:],LD)
        for j in range(4):out[j]=p(left,nu=j)/factorial(j)
        if key not in KEYS[-2:]:
            for field,q in fields.items():
                a,b=d['geometry_'+field+'_'+key];slope=b[idx];value=a[idx]+(left-gt[idx])[:,None]*slope
                for j in range(4):
                    v=q(left,nu=j)/(factorial(j)*AMP);out[j]+=value*v;out[j+1]+=slope*v
        result[key]=out
    return knots,result


def source():
    assert read(OUT/'geometry-result.json')['passed'];s=(base.OUT/'expanded-source.py').read_text();replace=base.replace
    old="unit=np.zeros((5,m.n),LD);unit[0]=m.driver.zc['alpha'];unit[2]=m.driver.centers*m.driver.zc['Phi']"
    s=replace(s,old,"unit=np.zeros((5,m.n),LD)")
    old="gt=np.unique(np.r_[m.t[m.t<T],T]);gg=[readout(v,zero,zero_m,zero_p,unit,False) for v in gt]\n        geometry={k:np.array([[a[k] for a in gg[:-1]],[(b[k]-a[k])/(v-u) for a,b,u,v in zip(gg[:-1],gg[1:],gt[:-1],gt[1:])]]) for k in KEYS}"
    new="gt=np.unique(np.r_[m.t[m.t<T],T]);geometry={}\n        for name,index in [('u',0),('lambda',2)]:\n            unit[:]=0;unit[index]=1;gg=[readout(v,zero,zero_m,zero_p,unit,False) for v in gt]\n            geometry[name]={k:np.array([[a[k] for a in gg[:-1]],[(b[k]-a[k])/(v-u) for a,b,u,v in zip(gg[:-1],gg[1:],gt[:-1],gt[1:])]]) for k in KEYS}"
    s=replace(s,old,new)
    s=replace(s,"phi=m.driver.wave(float(t+h*u),m.driver.xc)[0]/(m.driver.centers*AMP)","metric=m.driver.returned(float(t+h*u))")
    s=replace(s,"geometry_at(geometry,gt,np.array([t+h*u]),k)[0]*phi","sum(geometry_at(geometry[name],gt,np.array([t+h*u]),k)[0]*metric['delta_'+name]/AMP for name in ['u','lambda'])")
    s=replace(s,"d.update({'geometry_coeff_'+k:v for k,v in geometry.items()})","d.update({'geometry_'+name+'_'+k:v for name,g in geometry.items() for k,v in g.items()})")
    old="d.update(geometry_times=gt,drive_x=m.driver.xc,drive_duration=m.driver.D,drive_centers=m.driver.centers,drive_radius=m.driver.r0,drive_amplitude=incident.ETA,drive_times=m.driver.times,drive_born=np.array([np.interp(m.driver.xc,m.driver.tx,v) for v in m.driver.born['U']]))"
    s=replace(s,old,"d.update(geometry_times=gt,metric_times=m.driver.clock,u_coeff=m.driver.u.c,lambda_coeff=m.driver.lam.c)")
    def geometry_at(g,gt,t,key):
        ids=np.clip(np.searchsorted(gt,t,side='right')-1,0,len(gt)-2);a,b=g[key];return a[ids]+(t-gt[ids])[:,None]*b[ids]
    endpoint=SimpleNamespace(**dict(vars(base.endpoint),initialize=initialize))
    ns=dict(base.prior.prior.source.__globals__,OUT=OUT,INPUT=prior.INPUT,prior=endpoint,
        base=SimpleNamespace(**dict(vars(base.base),saved=prior.saved)),coefficients=coefficients,
        restored_gas=base.recovery.prior.restored_gas,geometry_at=geometry_at,background_contrasts=base.background_contrasts)
    exec(compile(s,__file__,'exec'),ns);(OUT/'expanded-source.py').write_text(s);ns['source']()
    r=read(OUT/'sources.json');r.update(actual_returned_metric=True,metric_floor_jumps_preserved=True,
        sampled_metric_interpolation=read(OUT/'geometry-result.json'),uniform_temporal_error_certificate=False,
        new_material_steps=0,new_photon_solves=0,new_GR_fields=0,full_goal_complete=False)
    write(OUT/'result.json',r)


if __name__=='__main__':
    action=sys.argv[1].removeprefix('check_');assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    os.sched_setaffinity(0,{8});resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));base.endpoint.evolution.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
