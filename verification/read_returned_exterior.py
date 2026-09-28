"""Same-solution Radau emission through the existing frozen exterior operator.

This reads separate high/low components. It does not certify an evolving
exterior, a background mass remap, or the final observable charge.
"""
from pathlib import Path
from decimal import Decimal, localcontext
import json, os, resource, sys, time
import numpy as np
import sympy as sp
import read_compensated_return_charge as charge
import return_complete_history_gr as common
import extend_retarded_history as full
import def_native_global_scalar_closure as exterior

ROOT=Path('native-returned-exterior252-work')
read,write,sha,bind=charge.read,charge.write,charge.sha,charge.bind
LD=np.longdouble; C=exterior.C; G=exterior.G


class Emission:
    def __init__(self,edges,lum):
        self.edges=np.asarray(edges);self.h=np.diff(self.edges).astype(LD)
        assert self.edges[0]==0 and np.all(self.h>0)
        lum=np.asarray(lum,LD).reshape(-1,2,4)
        assert len(lum)==len(self.h)
        self.a=LD('1.5')*lum[:,0]-LD('.5')*lum[:,1]
        self.b=LD('1.5')*(lum[:,1]-lum[:,0])/self.h[:,None]
        h=self.h[:,None]
        increments=h*(LD('.75')*lum[:,0]+LD('.25')*lum[:,1])
        self.q=np.vstack([np.zeros(4,LD),np.cumsum(increments,axis=0)])
        self.r=np.vstack([np.zeros(4,LD),np.cumsum(h*self.q[:-1]+h*h*self.a/2+h*h*h*self.b/6,axis=0)])
        self.absolute=np.sum(h*(LD('.75')*abs(lum[:,0])+LD('.25')*abs(lum[:,1])),axis=0)

    def primitives(self,at):
        x=np.clip(np.asarray(at),0,self.edges[-1])
        j=np.clip(np.searchsorted(self.edges,x,side='right')-1,0,len(self.h)-1)
        s=(x-self.edges[j]).astype(LD)[...,None]
        return (self.q[j]+s*self.a[j]+s*s*self.b[j]/2,
                self.r[j]+s*self.q[j]+s*s*self.a[j]/2+s*s*s*self.b[j]/6)

    def check(self):
        query=np.unique(np.r_[self.edges,(self.edges[:-1]+self.edges[1:])/2])
        # Independent integration of every interval, including the second
        # primitive's elapsed-time factor. No cumulative recurrence reused.
        dt=np.maximum(query[:,None].astype(LD)-self.edges[:-1].astype(LD),0)
        s=np.minimum(dt,self.h);integral=s[...,None]*self.a+s[...,None]**2*self.b/2
        second=(dt*s-s*s/2)[...,None]*self.a+(dt*s*s/2-s*s*s/3)[...,None]*self.b
        direct=[integral.sum(1,dtype=LD),second.sum(1,dtype=LD)]
        errors=[float(np.max(abs(a-b))/max(np.max(abs(b)),LD('1e-290'))) for a,b in zip(self.primitives(query),direct)]
        assert max(errors)<1e-12,errors
        return errors


def self_check():
    x,h,L0,L1=sp.symbols('x h L0 L1',positive=True)
    p=sp.Rational(3,2)*(1-x/h)*L0+(sp.Rational(3,2)*x/h-sp.Rational(1,2))*L1
    assert sp.simplify(sp.integrate(p,(x,0,h))-h*(3*L0+L1)/4)==0
    m=Emission([0.,.1,.37,1.],np.arange(24).reshape(6,4)-12)
    m.check();zero=Emission([0.,1.],np.zeros((2,4)))
    assert all(np.count_nonzero(v)==0 for v in zero.primitives([0.,.5,1.]))
    a,s,k,e=sp.symbols('a s k e')
    assert sp.factor((a+s)/(1+k-e)-a-(s+a*(e-k))/(1+k-e))==0
    return dict(classification='Proven',passed=True,scope='Radau dense emission primitives and rational readout identity only.')


def run(scope):
    assert scope in ['common','full'];out=ROOT/scope;assert not out.exists();out.mkdir(parents=True)
    actual=common.OUT if scope=='common' else Path('native-full-return249-work/full')
    high=common.OUT if scope=='common' else full.OUT
    low=charge.OUT if scope=='common' else Path('native-full-charge251-work/full/charge')
    primary_saved=common.saved if scope=='common' else full.source.saved
    accepted=read(low/'result.json');assert accepted['charge_comparison_admitted']
    assert read(actual/'result.json')['passed']
    files=[Path(__file__),Path(exterior.__file__),Path(common.prior.__file__),low/'result.json',actual/'result.json']
    data={};port_checks=[]
    for n in [64,128]:
        for label,folder,path in [('high',high,primary_saved(n)),('low',low,actual/f'sweep-1/photons/return-{n}.npz')]:
            p=dict(np.load(path));d=dict(np.load(folder/f'gr/source-{n}.npz'))
            edges=p['actual_step_edges'];count=2*(len(edges)-1)
            assert np.array_equal(p['accepted_angular_times'][:count],p['joint_stage_times'])
            emission=Emission(edges,p['accepted_angular_luminosity'][:count]);checks=emission.check()
            aw=np.arange(1,8,2,dtype=LD)/32
            q=emission.primitives(edges)[0]@aw
            reference=p['radial_ports'][:,1,1]
            ids=np.array([np.argmin(abs(edges-t)) for t in p['t']]);assert np.max(abs(edges[ids]-p['t']))<1e-18
            norm=max(emission.absolute@aw,LD('1e-290'))
            port=float(np.max(abs(q[ids]-reference))/norm)
            source_port=float(abs(q[-1]-d['outer_cumulative_energy_erg'][-1])/norm)
            assert max(port,source_port)<1e-12,(label,n,port,source_port)
            data[label,n]=(emission,d)
            port_checks.append(dict(component=label,clock=n,port=port,source_port=source_port,primitives=checks))
            files += [path,folder/f'gr/source-{n}.npz']
            files += [folder/f'gr/fields-{n}-g{q}.json' for q in ([4,8] if n==128 else [8])]
            if label=='high':files += [actual/f'metric/metric-{n}-g{q}.npz' for q in ([4,8] if n==128 else [8])]
    plan=dict(classification='Conjectural',claim='Read the same accepted high/low histories through the existing frozen null-infinity exterior kernel, retaining their actual dense Radau emission and signed mass-constraint ports.',
        decision='Measure exterior and mass-normalization changes of the same compact endpoint. Physical null-infinity charge remains unadjudicated until evolving exterior, physical energy conversion, background remap, source errors and selfGR are closed.',
        gates=dict(port=1e-12,primitive=1e-12,time=.02,angular=.002,radial=.002,mass_owner=1e-12),
        budget=dict(seconds=1800,virtual_GiB=16,CPU_threads=1,new_physical_steps=0),
        reuse='Completed coupled histories and compact fields; existing 4/8 exterior quadrature with17causal cuts. No new physical bins, time clocks, GR return iterate or tolerance relaxation.',
        forecast='Historical same exterior kernel is seconds per path; variable Radau primitives add interval lookup only. One production is allowed30minutes including initialization; keep successful outputs on failure.',
        limits='Reference-energy emission on frozen initial rays, not the complete physical coordinate-energy flux on an evolving exterior. Nonzero homogeneous mass ports are applied, not fitted away. Initial mass normalization here is not the previous evolving-background normalization or an observational EFT comparison.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},full_goal_complete=False)
    write(out/'plan.json',plan);write(out/'symbolic.json',self_check())
    init=out/'initialization'
    for p in list((low/'sweep-0').rglob('*.npz'))+[low/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
        dst=init/p.relative_to(low);dst.parent.mkdir(parents=True,exist_ok=True);os.link(p,dst)
    for part in ['sweep-1/photons','sweep-1/material','gr']:(init/part).mkdir(parents=True,exist_ok=True)
    bind(charge.base.endpoint.initialize,OUT=init)()
    m=exterior.Exterior();T=data['high',128][0].edges[-1];m.T=float(T);m.t=np.linspace(0,T,17)
    assert m.M==float(data['high',128][1]['M_cm']) and abs(m.K/float(data['high',128][1]['K_cm'])-1)<1e-12
    constraints=common.prior.geometry.constraints
    rows=[];paths={}
    for n,a,r in [(128,8,8),(128,4,8),(128,8,4),(64,8,8)]:
        begin=time.monotonic();kernel=m.kernel(a,r)
        for label in ['high','low']:
            emission,d=data[label,n];values=[]
            for now in m.t:
                h,hh=emission.primitives(np.maximum(now-kernel['delay'],0));idx=np.arange(len(h));bins=kernel['node_bins']
                mass=-G/C**3*np.sum(kernel['weights']*kernel['mass']*hh[idx,bins],dtype=LD)
                stress=-G/(2*C**4)*np.sum(kernel['weights']*kernel['stress']*h[idx,bins],dtype=LD)
                h,_=emission.primitives(np.maximum(now-kernel['infinity'],0));idx=np.arange(len(h))
                arrived=np.sum(kernel['mw']*kernel['mu']*h[idx,kernel['bins']],dtype=LD)
                values.append([-(mass+stress)/m.M,arrived])
            values=np.array(values,LD);dend=dict(d)
            for key in ['t',*charge.KEYS]:dend[key]=d[key][-2:]
            zero=np.zeros((2,len(d['edges'])))
            c=constraints.centers(m.response,dend,dict(delta_phi=zero,delta_Phi=zero),r)
            homogeneous=c['J'][-1,-1]*m.response.tz['lapse'][-1]/np.sqrt(m.response.tz['b'][-1])+G/C**4*d['outer_cumulative_energy_erg'][-1]
            owner_error=None
            if label=='high':
                metric=actual/f'metric/metric-{n}-g{r}.npz'
                owner=np.load(metric)['asymptotic_mass_residual_cm'][-1]
                owner_error=float(abs(homogeneous-owner)/max(abs(owner),1e-290));assert owner_error<1e-12,owner_error
            folder=high if label=='high' else low
            compact=read(folder/f'gr/fields-{n}-g{r}.json')['endpoint_compact_with_metric']
            with localcontext() as ctx:
                ctx.prec=100
                def D(v):
                    num,den=LD(v).as_integer_ratio();return Decimal(num)/Decimal(den)
                alpha=-D(m.K)/D(m.M);epsilon=D(G)/D(C)**4*D(values[-1,1])/D(m.M);kappa=D(homogeneous)/D(m.M)
                scalar=D(compact)+D(values[-1,0]);normalized=(scalar+alpha*(epsilon-kappa))/(1+kappa-epsilon)
            row=dict(component=label,clock=n,angular=a,radial=r,compact=compact,
                exterior=float(values[-1,0]),arrived_energy_erg=float(values[-1,1]),homogeneous_mass_cm=float(homogeneous),
                mass_owner_relative=owner_error,normalized_standalone=str(normalized),relative_change_from_compact=float((normalized-D(compact))/D(compact)))
            row.update(scalar_numerator=str(scalar),epsilon=str(epsilon),kappa=str(kappa))
            paths[label,n,a,r]=values;rows.append(row)
            np.savez_compressed(out/f'{label}-{n}-a{a}-r{r}.npz',t=m.t,exterior=values[:,0],arrived_energy_erg=values[:,1])
            write(out/f'{label}-{n}-a{a}-r{r}.json',row)
        print(json.dumps(dict(path=[n,a,r],seconds=time.monotonic()-begin)),flush=True)
    controls={}
    for label in ['high','low']:
        fine=paths[label,128,8,8];norm=np.maximum(np.max(abs(fine),axis=0),LD('1e-290'))
        controls[label]={name:[float(v) for v in np.max(abs(paths[(label,)+setting]-fine),axis=0)/norm]
            for name,setting in [('angular',(128,4,8)),('radial',(128,8,4)),('time',(64,8,8))]}
    passed=all(max(v)<plan['gates'][name] for group in controls.values() for name,v in group.items())
    components=[]
    with localcontext() as ctx:
        ctx.prec=100
        for n in [64,128]:
            highrow,lowrow=[next(v for v in rows if v['component']==label and (v['clock'],v['angular'],v['radial'])==(n,8,8)) for label in ['high','low']]
            qh=Decimal(highrow['normalized_standalone']);den=1+Decimal(highrow['kappa'])-Decimal(highrow['epsilon'])
            dm=Decimal(lowrow['kappa'])-Decimal(lowrow['epsilon'])
            increment=(Decimal(lowrow['scalar_numerator'])-(alpha+qh)*dm)/(den+dm)
            total=qh+increment
            components.append(dict(clock=n,high=str(qh),same_solution_low_increment=str(increment),total=str(total),
                nominal_compact_sign_unchanged=bool(total*Decimal.from_float(highrow['compact'])>0)))
    result=dict(classification='Counterexample candidate',passed=passed,controls=controls,rows=rows,port_checks=port_checks,
        components=components,
        same_actual_returned_solution_read=True,conditional_frozen_exterior_only=True,
        physical_coordinate_energy_conversion_complete=False,complete_physical_exterior=False,background_mass_remap_complete=False,self_GR_return_closed=False,
        physical_final_charge_solved=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(out/'result.json',result);assert passed,result


if __name__=='__main__':
    action=sys.argv[1];ROOT.mkdir(exist_ok=True);start=time.monotonic();error=None
    assert action in ['check','common','full']
    if action!='check':
        assert read(ROOT/'check.json')['passed'] and read(ROOT/'check-receipt.json')['source_sha256']==sha(__file__)
        assert not (ROOT/f'{action}-start.json').exists()
        write(ROOT/f'{action}-start.json',dict(pid=os.getpid(),process_start_ticks=Path('/proc/self/stat').read_text().split()[21],
            boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),started_unix=time.time(),source_sha256=sha(__file__)))
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,)*2)
    common.prior.prior.joint.previous.original.inf.incident.native.deadline(1800)
    try:
        if action=='check':write(ROOT/'check.json',self_check())
        else:run(action)
    except BaseException as exc:error=repr(exc);raise
    finally:write(ROOT/f'{action}-receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
