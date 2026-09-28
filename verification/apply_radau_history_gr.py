"""Counterexample candidate: carry actual Radau polynomials into retarded GR.

Uses the original188 horizon and both saved trajectories. No new fluid steps,
time paths, fitted moments or source-gate relaxation. Floor jumps stay one-sided.
"""
from pathlib import Path
from types import FunctionType
import gc,inspect,json,os,resource,sys,time
import numpy as np
from numpy.polynomial import legendre as leg
from scipy.interpolate import PPoly
import return_resolved_joint_history as prior
import complete_joint_photon_recovery as recovery

OUT=Path('native-radau-gr211-work');OLD=prior.OUT
base,run=prior.base,prior.run
LD,AMP,C=prior.LD,prior.AMP,prior.C
read,write,sha=prior.read,prior.write,prior.sha
CAPS=dict(prepare=180,source=1200,fields=1200)
KEYS=['baryon_g','gas_nonrest_energy_erg','nonrest_trace_erg','nonrest_stress_erg','pressure_volume_erg','photon_energy_erg','photon_radial_pressure_erg','metric_stress_erg','inner_cumulative_energy_erg','outer_cumulative_energy_erg']


def cubic(y):
    """Coefficients in theta from values at 0,1/3,2/3,1; no least-squares fit."""
    a,b,c,d=np.asarray(y,LD)
    return np.array([a,(-11*a+18*b-9*c+2*d)/2,(18*a-45*b+36*c-9*d)/2,(-9*a+27*b-27*c+9*d)/2])


def poly(edges,coefficients):
    h=np.diff(edges);a=np.asarray(coefficients,float)
    power=np.arange(len(a)).reshape((-1,)+(1,)*(a.ndim-1))
    scale=h.reshape((1,-1)+(1,)*(a.ndim-2))**power
    return PPoly((a/scale)[::-1],edges)


def evaluate(coefficients,theta):
    out=np.zeros_like(coefficients[0])
    for c in coefficients[::-1]:out=out*LD(theta)+c
    return out


def quadratic(initial,pair,theta):
    a,b=pair;u=LD(theta)
    return initial+u*((9*a-b-8*initial)/2+u*(3*b-9*a+6*initial)/2)


class Response(base.gr.Response):
    def setup(self,d,order):
        samples=dict(d)
        for k in KEYS:samples[k]=d['radau_'+k].reshape((-1,)+d[k].shape[1:])
        samples['t']=np.concatenate([a+(b-a)*np.arange(4)/3 for a,b in zip(d['t'][:-1],d['t'][1:])])
        super().setup(samples,order)
        source=cubic(self.source.reshape(-1,4,len(self.ids)).swapaxes(0,1))
        direct=cubic(self.direct.reshape(-1,4,len(self.ids)).swapaxes(0,1))
        super().setup(d,order)
        self.source=poly(self.t,source);self.direct=poly(self.t,direct)

    def propagate(self,source):
        p=source if isinstance(source,PPoly) else base.gr.base.flow.green.polynomial(self.t,source)
        H=p.antiderivative();hc=H.c.reshape(len(H.c),len(self.t)-1,-1,self.order)/self.dx.reshape(-1,self.order)
        coefficients=np.einsum('cij,ptcj->ptci',self.inverse,hc)
        gx,gw=leg.leggauss((self.order+len(p.c)+1)//2)
        U=[];Ut=[];Ux=[];cols=np.arange(len(self.ids))[None,:]
        for t in self.t:
            ret=t-self.distance;cut=np.clip(ret,0,self.t[-1]);idx=np.clip(np.searchsorted(self.t,cut,side='right')-1,0,len(self.t)-2);dt=cut-self.t[idx]
            hv=np.zeros_like(ret);sv=np.zeros_like(ret)
            for co in H.c:hv=hv*dt+co[idx,cols]
            for co in p.c:sv=sv*dt+co[idx,cols]
            hv[ret<=0]=0;sv[ret<=0]=0
            hv=hv.reshape(len(self.tx),-1,self.order).sum(2)
            dv=(sv*self.sign).reshape(len(self.tx),-1,self.order).sum(2);sv=sv.reshape(len(self.tx),-1,self.order).sum(2)
            owner,cell,lo,hi=self.segments(t)
            if len(cell):
                xx=(lo[:,None]+hi[:,None])/2+(hi-lo)[:,None]*gx/2
                rr=t-abs(self.tx[owner,None]-xx)/C
                jj=np.clip(np.searchsorted(self.t,np.clip(rr,0,self.t[-1]),side='right')-1,0,len(self.t)-2)
                tt=np.clip(rr,0,self.t[-1])-self.t[jj]
                vv=leg.legvander((xx-self.mid[cell,None])/self.half[cell,None],self.order-1)
                cc=np.array([np.sum(co[jj,cell[:,None]]*vv,axis=-1) for co in coefficients])
                hh=np.zeros_like(rr);ss=np.zeros_like(rr)
                for co in cc:hh=hh*tt+co
                for k,co in enumerate(cc[:-1]):ss=ss*tt+(len(cc)-1-k)*co
                hh[rr<=0]=0;ss[rr<=0]=0;ww=(hi-lo)[:,None]*gw/2
                ii,jj=np.unique(np.column_stack([owner,cell]),axis=0).T
                hv[ii,jj]=0;sv[ii,jj]=0;dv[ii,jj]=0
                np.add.at(hv,(owner,cell),np.sum(ww*hh,axis=1))
                np.add.at(sv,(owner,cell),np.sum(ww*ss,axis=1))
                np.add.at(dv,(owner,cell),np.sum(ww*ss*np.sign(self.tx[owner,None]-xx),axis=1))
            U.append(C/2*hv.sum(1,dtype=LD));Ut.append(C/2*sv.sum(1,dtype=LD));Ux.append(-.5*dv.sum(1,dtype=LD))
        return np.asarray(U,float),np.asarray(Ut,float),np.asarray(Ux,float)


def symbolic():
    import sympy as s
    x=s.symbols('x');a,b,c=s.symbols('a b c')
    q=a+x*((9*b-c-8*a)/2+x*(3*c-9*b+6*a)/2)
    assert s.simplify(q.subs(x,0)-a)==0 and s.simplify(q.subs(x,s.Rational(1,3))-b)==0 and s.simplify(q.subs(x,1)-c)==0
    errors=[]
    for order in [4,8]:
        m=Response.__new__(Response);m.order=order;m.t=np.array([0.,.3,.7,1.])/C
        m.xfaces=np.array([-1.,0.,1.]);m.mid=np.array([-.5,.5]);m.half=np.array([.5,.5]);gx,gw=leg.leggauss(order)
        m.x=(m.mid[:,None]+m.half[:,None]*gx).ravel();m.ids=np.repeat(np.arange(2),order)
        m.dx=(m.half[:,None]*np.broadcast_to(gw,(2,order))).ravel();m.tx=np.array([-.75,0.,.3,2.])
        m.inverse=np.linalg.inv(leg.legvander(np.broadcast_to(gx,(2,order)),order-1))
        m.distance=abs(m.tx[:,None]-m.x[None,:])/C;m.sign=np.sign(m.tx[:,None]-m.x[None,:])
        for degree in range(4):
            values=np.array([[(C*(left+(right-left)*theta))**degree*m.dx for left,right in zip(m.t[:-1],m.t[1:])] for theta in np.arange(4)/3])
            coefficients=cubic(values);u,ut,ux=m.propagate(poly(m.t,coefficients));exact=[]
            for t in m.t*C:
                primitive=lambda y,power:np.sign(y)*(t**power-max(t-abs(y),0.)**power)/power
                exact.append([[(primitive(1-x,degree+2)-primitive(-1-x,degree+2))/(2*(degree+1)),
                    (primitive(1-x,degree+1)-primitive(-1-x,degree+1))/2,
                    (max(t-abs(x+1),0.)**(degree+1)-max(t-abs(x-1),0.)**(degree+1))/(2*(degree+1))] for x in m.tx])
            errors.append(float(np.max(abs(np.stack([u,ut/C,ux],axis=-1)-exact))))
    assert max(errors)<3e-14,errors
    return dict(classification='Proven',passed=True,polynomial_degrees=[0,1,2,3],box_U_Ut_over_c_Ux_errors=errors,scope='Radau interpolation identities and characteristic polynomial controls only; no physical error certificate.')


def prepare():
    assert not OUT.exists();OUT.mkdir();(OUT/'gr').mkdir()
    assert recovery.read(recovery.OUT/'result.json')['original_endpoint_and_ledger_passed']
    assert not read(OLD/'original-source-norm.json')['passed']
    files=list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']];reused={}
    for p in files:
        dst=OUT/p.relative_to(OLD);dst.parent.mkdir(parents=True,exist_ok=True);os.link(p,dst);reused[str(dst)]=sha(p)
    for s in ['sweep-1/photons','sweep-1/material']:(OUT/s).mkdir(parents=True,exist_ok=True)
    files += [OLD/'original-source-norm.json',prior.evolution.OUT/'result.json',Path('native-early-gr207-work/gr/fields-128-g8.npz')]
    files += [OLD/'gr'/f'source-{n}.npz' for n in [64,128]]+[base.saved(n) for n in [64,128]]
    files += [prior.REC/'recovered-64.npz',recovery.OUT/'recovered-128.npz',recovery.OUT/'result.json']
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='a72a3a12d',
        claim='Replace the dominant straight-line time representation by the actual Radau stage polynomial in the same GR source and characteristic operator, preserving physical source-state failures.',
        scope='Original188T/64,2coarse/4fineacceptedsteps and531cells. Four algebraic evaluations per existing step construct cubic source coefficients; three additional interior probes validate the readout. No new physical evolution or time path.',
        method='Quadratic dense gas/photon moments through initial and actual two Radau stages. Preserve post-floor initial states and pre-floor closing states on each interval. Integrate the original Radau boundary rates with their dense weights. Apply the same pressure/background map, cubic source polynomial, characteristic cuts, and independently integrated Jordan-radius direct readout.',
        decision='Apply the time-representation repair to GR even if the original source gate still fails, explicitly as an unadmitted readout. Only unchanged original source AND field time gates would admit actual GR return. Never count repaired representation as repaired physical state or final charge.',
        limitations='Potential-feedback time interpolation remains piecewise linear on the original edges and is measured separately. Missing exterior scalar, EOS, uniform derivatives, spatial/nonlinear/observational closure remain. Actual floor jumps are not smoothed.',
        gates=dict(dense_stage=1e-12,polynomial=1e-12,pressure=.002,mapping=1e-12,source_time=.02,field_time=.02,quadrature=.002,independent_GR=1e-9),budgets=CAPS,CPU_threads=1,virtual_GiB=6,
        forecast='206source38s and207fourGRreads10s. Additional algebraic pressure readouts and cubic integration have not been timed; allow20minutes each with6GiB. The live210producer and all its dependencies remain unchanged.',
        stop='Input/reconstruction/polynomial/pressure/numerical-control failure stops the repair. Preserve source-time failure separately from operator verification; no metric or physical return is executed here.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},reused=reused,final_charge_conclusion='unadjudicated',full_goal_complete=False))
    write(OUT/'symbolic.json',symbolic())


def source():
    FunctionType(prior.initialize.__code__,dict(prior.initialize.__globals__,OUT=OUT))()
    pressure=FunctionType(prior.precision.pressure_function.__code__,dict(prior.precision.pressure_function.__globals__,OUT=OUT))()
    model=base.gr.Response();checks=[]
    for n in [64,128]:
        m=run.owner.Model(n);p=dict(np.load(base.saved(n)));count=n//32
        photon=dict(np.load((prior.REC if n==64 else recovery.OUT)/f'recovered-{n}.npz'))
        d=dict(np.load(OLD/'gr'/f'source-{n}.npz'));assert len(d['t'])==count+1
        material=m.material;cf=model.coeff(material.rE)
        Eg0,Pg0,Kg0,Er0,Pr0=[cf[key]*C**4/base.gr.base.G*material.V for key in ['Eg','Pg','Kg','Er','Pr']]
        b=model.model.bulk;weights=4*np.pi*np.r_[b.W,model.model.W][:,None,None]*b.w[None,:,None]*b.d['num']*b.d['Einf']
        banks={k:dict(np.load(base.feedback.old.OUT/f'bank-{m.reference}/point-{k}.npz')) for k in [0,1]}
        probes=[];mapping=[];references=[];values={k:[] for k in KEYS};dense_errors=[];poly_errors=[];endpoint_errors=[]
        def readout(t,g,mom,ports):
            z=np.array([g[:,2]*m.bu,g[:,3]*m.su,g[:,0]*m.eu,g[:,1]*m.nu],LD)
            field=m.geometry(float(t))[0];vol=(3*field[0]+field[2])*AMP
            k=int(np.clip(np.searchsorted(m.t,t,side='right')-1,0,15));w=(t-m.t[k])/(m.t[k+1]-m.t[k]);assert k==0
            pp=[];qs=[]
            for j,v in [(k,1-w),(k+1,w)]:
                pp.append(v*np.array(pressure(material,j,z,field,banks[j])));qs.append(v*material.point(j)['Q'])
            p0,delta,error=sum(pp);Q=sum(qs);background=(Q[2]+material.rest*Q[0])/material.a
            nonrest=z[2]*AMP/material.a+(Eg0+Pg0-background)*vol;pg,pr=delta*AMP+(Kg0[None]-p0)*vol;B=z[0]*AMP
            local=m.local(float(t));reference=(m.pressure(float(t),g)+local['pressure_source'])*m.volume*AMP
            probes.append(error*AMP);references.append(reference);mapping.append(delta[0]*AMP-reference)
            I=(1-w)*m.I[k]+w*m.I[k+1];Ebg=np.einsum('nqf,nqf->n',I,weights)/material.a;Pbg=np.einsum('nqf,nqf,q->n',I,weights,b.mu2)/material.a
            phi,lam=field[0]*AMP,field[2]*AMP;em,pm=mom[:2]
            pe=em/material.a-Ebg*vol+4*Er0*phi+(Er0+Pr0)*lam
            pp=pm/material.a-Pbg*vol+4*Pr0*phi+(3*Pr0-model.ratio4*Er0)*lam
            return dict(zip(KEYS,[B,nonrest,nonrest-pr-2*pg,nonrest-pr,pg,pe,pp,B*LD(m.model.cx)*LD(C)**2+nonrest-pr+pe-pp,ports[0,1],ports[1,1]]))
        gas=lambda q:np.column_stack([(q[2]-m.kappa*q[0])/m.eu,q[3]/m.nu,q[0]/m.bu,q[1]/m.su])
        cumulative=np.zeros((2,2),LD)
        for step in range(count):
            times=p['joint_stage_times'][2*step:2*step+2];h=LD(4)*p['joint_stage_weights'][2*step+1];t=LD(p['actual_step_edges'][step])
            initial=np.zeros((m.n,4),LD) if step==0 else gas(p['joint_stage_conserved_scaled'][2*step-1]);initial[~material.active(t)]=0
            pair=np.array([gas(q) for q in p['joint_stage_conserved_scaled'][2*step:2*step+2]])
            g0=initial.copy();m0=np.zeros((3,m.n),LD) if step==0 else photon['photon_moments'][2*step-1]
            mp=photon['photon_moments'][2*step:2*step+2];port=photon['radial_ports'][2*step:2*step+2]
            def at(theta):
                u=LD(theta);q=np.array([LD('1.5')*u-LD('.75')*u*u,LD('.75')*u*u-LD('.5')*u])
                return readout(t+h*u,quadratic(g0,pair,u),quadratic(m0,mp,u),cumulative+h*np.einsum('j,jab->ab',q,port))
            rates=(p['joint_native_rates_scaled']+p['joint_collision_rates_scaled'])[2*step:2*step+2]/m.units
            for j,u in enumerate([LD(1)/3,LD(1)]):
                q=np.array([LD('1.5')*u-LD('.75')*u*u,LD('.75')*u*u-LD('.5')*u])
                defect=pair[j]-g0-h*np.einsum('j,jnk->nk',q,rates)
                dense_errors.append((np.sum(abs(defect)*m.units,axis=0)/np.maximum(np.sum(abs(pair[j])*m.units,axis=0),LD('1e-290'))).astype(float).tolist())
            rows=[at(LD(j)/3) for j in range(4)];co={k:cubic([r[k] for r in rows]) for k in KEYS}
            for u in [LD(1)/6,LD(1)/2,LD(5)/6]:
                actual=at(u)
                poly_errors.append({k:float(np.sum(abs(evaluate(co[k],u)-actual[k]))/max(max(np.sum(abs(r[k])) for r in rows),LD('1e-290'))) for k in KEYS})
            for k in KEYS:values[k].append([r[k] for r in rows])
            cumulative+=np.sum(p['joint_stage_weights'][2*step:2*step+2,None,None]*port,axis=0,dtype=LD)
            closing=pair[-1].copy();closing[~material.active(p['actual_step_edges'][step+1])]=0
            end=readout(p['actual_step_edges'][step+1],closing,mp[-1],cumulative)
            endpoint_errors.append({k:float(np.sum(abs(end[k]-d[k][step+1]))/max(np.max(np.sum(abs(d[k]),axis=-1)) if d[k].ndim>1 else np.max(abs(d[k])),LD('1e-290'))) for k in KEYS})
        norm=max(np.max(np.sum(abs(np.array(references)),axis=-1)),LD('1e-290'))
        row=dict(clock=n,dense_stage_max=float(np.max(dense_errors)),polynomial_max=max(max(r.values()) for r in poly_errors),endpoint_max=max(max(r.values()) for r in endpoint_errors),
            pressure_probe=float(np.max(np.sum(abs(np.array(probes)),axis=-1))/norm),pressure_mapping=float(np.max(np.sum(abs(np.array(mapping)),axis=-1))/norm),dense_stage=dense_errors,polynomial=poly_errors,endpoint=endpoint_errors)
        row['passed']=row['dense_stage_max']<1e-12 and row['polynomial_max']<1e-12 and row['endpoint_max']<1e-12 and row['pressure_probe']<.002 and row['pressure_mapping']<1e-12
        d.update({'radau_'+k:np.array(v) for k,v in values.items()});np.savez_compressed(OUT/'gr'/f'source-{n}.npz',**d);write(OUT/f'source-{n}-check.json',row)
        assert row['passed'],row;checks.append(row);del m,p;gc.collect()
    old=read(OLD/'original-source-norm.json')
    write(OUT/'sources.json',dict(classification='Counterexample candidate',representation_controls_passed=True,source_time_passed=old['passed'],unchanged_source_time=old['source_time_spatial_L1'],rows=checks,original_failure_preserved=True,GR_return_admitted=False,final_charge_conclusion='unadjudicated'))


def fields():
    assert read(OUT/'sources.json')['representation_controls_passed']
    FunctionType(prior.initialize.__code__,dict(prior.initialize.__globals__,OUT=OUT))()
    m=Response();fn=FunctionType(base.gr.base.Response.run.__code__,dict(base.gr.base.Response.run.__globals__,OUT=OUT/'gr'))
    rows=[fn(m,n,q) for n,q in [(128,8),(64,8),(128,4)]]
    fine,coarse,low=[dict(np.load(OUT/'gr'/f'fields-{n}-g{q}.npz')) for n,q in [(128,8),(64,8),(128,4)]]
    d=dict(np.load(OUT/'gr/source-128.npz'))
    trace=d['radau_baryon_g']*LD(d['cx'])*LD(C)**2+d['radau_nonrest_trace_erg']
    direct_source=inspect.getsource(base.charge.independent.direct)
    old="H=flow.green.polynomial(d['t'],trace).antiderivative()";assert direct_source.count(old)==1
    direct_source=direct_source.replace(old,"H=trace_poly.antiderivative()")
    ns=dict(base.charge.independent.direct.__globals__,trace_poly=poly(d['t'],cubic(trace.swapaxes(0,1))))
    exec(compile(direct_source,__file__,'exec'),ns);direct,coordinate=ns['direct'](m,d,8)
    controls=dict(quadrature=prior.aligned(low,fine,'U'),independent_GR=abs(direct-rows[0]['endpoint_direct'])/max(abs(direct),1e-290))
    times={key:prior.aligned(coarse,fine,key) for key in ['U','U_t','U_x']}
    old=dict(np.load(Path('native-early-gr207-work')/'gr/fields-128-g8.npz'))
    change={key:float(np.max(abs(fine[key]-old[key]))/max(np.max(abs(fine[key])),1e-290)) for key in ['U','U_t','U_x']}
    row=dict(classification='Counterexample candidate',numerical_controls_passed=controls['quadrature']<.002 and controls['independent_GR']<1e-9,controls=controls,time=times,
        repaired_source_time_representation_applied=True,change_from_straight_line=change,potential_relative=max(float(np.max(abs(fine['potential_U']))/max(np.max(abs(fine['U'])),1e-290)),0),
        source_time_passed=read(OUT/'sources.json')['source_time_passed'],field_time_passed=times['U']<.02,rows=rows,physical_steps=0,GR_return_evolution_executed=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    row['GR_return_admitted']=row['numerical_controls_passed'] and row['source_time_passed'] and row['field_time_passed']
    write(OUT/'result.json',row);print(json.dumps(row),flush=True);assert row['numerical_controls_passed'],row


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));prior.evolution.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();error=None
    try:
        if action!='prepare':
            for p,h in dict(read(OUT/'plan.json')['bindings'],**read(OUT/'plan.json')['reused']).items():assert sha(p)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
