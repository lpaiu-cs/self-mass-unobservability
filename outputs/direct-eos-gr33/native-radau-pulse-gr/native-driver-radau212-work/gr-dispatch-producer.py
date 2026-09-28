"""Carry Radau state AND the actual compact degree8 driver into the same GR."""
from pathlib import Path
from types import FunctionType
from math import comb,factorial
import inspect,json,os,resource,sys,time,textwrap
import numpy as np
from scipy.interpolate import PPoly
import apply_radau_history_gr as prior
import def_native_incident_drive as incident

OUT=Path('native-driver-radau212-work');OLD=prior.OUT
base,read,write,sha,LD,C,AMP=prior.base,prior.read,prior.write,prior.sha,prior.LD,prior.C,prior.AMP
KEYS=prior.KEYS;CAPS=dict(prior.CAPS,source_retry=1200)


def coefficients(d):
    """Exact algebra for the declared finite time representation; no resampling."""
    T=d['t'][-1];x=d['drive_x'];D=LD(d['drive_duration'])
    knots=np.unique(np.r_[d['t'],d['drive_times'],-x/C,D-x/C]);knots=knots[(knots>=0)&(knots<=T)]
    left=knots[:-1];mid=(left+knots[1:])/2
    phase=(left[:,None]+x[None]/C)/D;active=((mid[:,None]+x[None]/C)/D>0)&((mid[:,None]+x[None]/C)/D<1)
    pulse=np.zeros((9,len(left),len(x)),LD)
    for k,weight in zip(range(4,9),[256,-1024,1536,-1024,256]):
        for j in range(k+1):pulse[j]+=weight*comb(k,j)*phase**(k-j)/D**j
    pulse*=active[None]
    pulse*=d['drive_amplitude']*d['drive_radius']/(d['drive_centers'][None,None]*AMP)
    times=d['drive_times'];ids=np.clip(np.searchsorted(times,left,side='right')-1,0,len(times)-2)
    h=times[ids+1]-times[ids];slope=(d['drive_born'][ids+1]-d['drive_born'][ids])/h[:,None]
    pulse[0]+=(d['drive_born'][ids]+(left-times[ids])[:,None]*slope)/(d['drive_centers'][None]*AMP)
    pulse[1]+=slope/(d['drive_centers'][None]*AMP)
    output={}
    for key in KEYS:
        p=prior.poly(d['t'],d['state_coeff_'+key]);shape=(10,len(left))+d[key].shape[1:]
        co=np.zeros(shape,LD)
        for j in range(4):co[j]=p(left,nu=j)/factorial(j)
        if key not in KEYS[-2:]:
            a,b=d['geometry_coeff_'+key];g=a[None]+left[:,None]*b
            co[:9]+=pulse*g[None];co[1:]+=pulse*b[None,None]
        output[key]=co
    return knots,output


class Response(prior.Response):
    def setup(self,d,order):
        knots,co=coefficients(d);shape=co['baryon_g'].shape[:2];samples=dict(d)
        for key in KEYS:samples[key]=co[key].reshape((-1,)+d[key].shape[1:])
        samples['t']=np.zeros(shape[0]*shape[1])
        base.gr.Response.setup(self,samples,order)
        source=self.source.reshape(*shape,-1).copy();direct=self.direct.reshape(*shape,-1).copy()
        base.gr.Response.setup(self,d,order)
        self.source=PPoly(source[::-1],knots);self.direct=PPoly(direct[::-1],knots)


# Keep the checked characteristic algorithm; its integration cuts must follow
# each input polynomial, independently of the requested output edge times.
s=inspect.getsource(prior.Response.propagate)
s=s.replace('        H=p.antiderivative();','        knots=p.x\n        H=p.antiderivative();')
s=s.replace("len(self.t)-1","len(knots)-1").replace("np.clip(ret,0,self.t[-1])","np.clip(ret,0,knots[-1])")
s=s.replace("np.searchsorted(self.t,cut","np.searchsorted(knots,cut").replace("len(self.t)-2","len(knots)-2").replace("dt=cut-self.t[idx]","dt=cut-knots[idx]")
s=s.replace("owner,cell,lo,hi=self.segments(t)","saved=self.t;self.t=knots\n            owner,cell,lo,hi=self.segments(t);self.t=saved")
s=s.replace("np.searchsorted(self.t,np.clip(rr,0,self.t[-1])","np.searchsorted(knots,np.clip(rr,0,knots[-1])").replace("tt=np.clip(rr,0,self.t[-1])-self.t[jj]","tt=np.clip(rr,0,knots[-1])-knots[jj]")
ns=dict(prior.Response.propagate.__globals__);exec(compile(textwrap.dedent(s),__file__,'exec'),ns);Response.propagate=ns['propagate']


def symbolic():
    source=inspect.getsource(prior.symbolic)
    source=source.replace('range(4):','range(10):')
    before="values=np.array([[(C*(left+(right-left)*theta))**degree*m.dx for left,right in zip(m.t[:-1],m.t[1:])] for theta in np.arange(4)/3])\n            coefficients=cubic(values);u,ut,ux=m.propagate(poly(m.t,coefficients));exact=[]"
    after="co=np.array([comb(degree,j)*(C*m.t[:-1,None])**(degree-j)*C**j*m.dx[None] for j in range(degree+1)])\n            u,ut,ux=m.propagate(PPoly(co[::-1],m.t));exact=[]"
    assert before in source;source=source.replace(before,after).replace('polynomial_degrees=[0,1,2,3]','polynomial_degrees=list(range(10))')
    scope=dict(prior.symbolic.__globals__,Response=Response,comb=comb,PPoly=PPoly);exec(compile(source,__file__,'exec'),scope)
    return scope['symbolic']()


def prepare():
    assert not OUT.exists();OUT.mkdir();(OUT/'gr').mkdir()
    assert read(OLD/'source-receipt.json')['error'] is not None
    assert read(OLD/'source-64-check.json')['polynomial_max']>.04
    files=list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']];reused={}
    for p in files:
        dst=OUT/p.relative_to(OLD);dst.parent.mkdir(parents=True,exist_ok=True);os.link(p,dst);reused[str(dst)]=sha(p)
    for p in ['sweep-1/photons','sweep-1/material']:(OUT/p).mkdir(parents=True,exist_ok=True)
    plan=read(OLD/'plan.json');bindings=dict(plan['bindings'])
    for p in [Path(__file__),OLD/'source-receipt.json',OLD/'source-64-check.json']:bindings[str(p)]=sha(p)
    plan.update(checkpoint='driver-aware-repair-after-e2927ee40',claim='Apply the actual Radau dense state with the unmodified analytic moving degree8 incident pulse and stored Born term to the same GR source; remove the false cubic-full-source assumption.',
        method='Split only algebraic GR source polynomials at existing step/Born knots and actual cell-center pulse arrival/departure times. Native material/photon state remains the original quadratic Radau trajectory. Its zero-geometry readout is cubic; the geometry map is affine in background time times the exact degree8 pulse, hence degree9 on each support segment. Preserve original physical edges, floor jumps, ports and all gates.',
        original_failure='211zero-state-degree assumption failed up to4.75percent despite exact endpoints and dense state. Preserve it; do not accept a low-order fit to the pulse or smooth its support.',
        decision='Verify independent intermediate full readout, actual original endpoints, characteristic monomials0..9 and independent Jordan integral. Then apply the repaired source to GR; original source-state time failure still blocks physical return.',
        forecast='211coarse source took21.69s. Two constructors and two background geometry maps plus the existing small Radau set should take under2minutes; retain20minute source/field caps. Actual pulse cut count is measured, not a new physical grid.',
        bindings=bindings,reused=reused)
    write(OUT/'plan.json',plan);write(OUT/'symbolic.json',symbolic())


def source():
    s=inspect.getsource(prior.source)
    changes=[
        ('def readout(t,g,mom,ports):','def readout(t,g,mom,ports,field_override=None,record=True):'),
        ("field=m.geometry(float(t))[0];vol=","field=m.geometry(float(t))[0] if field_override is None else field_override;vol="),
        ('probes.append(error*AMP);references.append(reference);mapping.append(delta[0]*AMP-reference)','if record:probes.append(error*AMP);references.append(reference);mapping.append(delta[0]*AMP-reference)'),
        ('cumulative=np.zeros((2,2),LD)',"cumulative=np.zeros((2,2),LD)\n        unit=np.zeros((5,m.n),LD);unit[0]=m.driver.zc['alpha'];unit[2]=m.driver.centers*m.driver.zc['Phi']\n        zero=np.zeros((m.n,4),LD);zero_m=np.zeros((3,m.n),LD);zero_p=np.zeros((2,2),LD);T=d['t'][-1]\n        first=readout(0.,zero,zero_m,zero_p,unit,False);last=readout(T,zero,zero_m,zero_p,unit,False)\n        geometry={k:np.array([first[k],(last[k]-first[k])/T]) for k in KEYS}"),
        ('def at(theta):','def at(theta,zero_field=False):'),
        ("cumulative+h*np.einsum('j,jab->ab',q,port))","cumulative+h*np.einsum('j,jab->ab',q,port),np.zeros((5,m.n),LD) if zero_field else None,not zero_field)"),
        ('rows=[at(LD(j)/3) for j in range(4)]','rows=[at(LD(j)/3,True) for j in range(4)]'),
        ('actual=at(u)\n                poly_errors.append',"actual=at(u)\n                phi=m.driver.wave(float(t+h*u),m.driver.xc)[0]/(m.driver.centers*AMP)\n                prediction={k:evaluate(co[k],u)+(geometry[k][0]+(t+h*u)*geometry[k][1])*phi if k not in KEYS[-2:] else evaluate(co[k],u) for k in KEYS}\n                poly_errors.append"),
        ("abs(evaluate(co[k],u)-actual[k])","abs(prediction[k]-actual[k])"),
        ("max(max(np.sum(abs(r[k])) for r in rows),LD('1e-290'))","max(np.sum(abs(actual[k])),max(np.sum(abs(r[k])) for r in rows),LD('1e-290'))"),
        ("d.update({'radau_'+k:np.array(v) for k,v in values.items()})", "d.update({'state_coeff_'+k:cubic(np.array(v).swapaxes(0,1)) for k,v in values.items()});d.update({'geometry_coeff_'+k:v for k,v in geometry.items()})\n        d.update(drive_x=m.driver.xc,drive_duration=m.driver.D,drive_centers=m.driver.centers,drive_radius=m.driver.r0,drive_amplitude=AMP,drive_times=m.driver.times,drive_born=np.array([np.interp(m.driver.xc,m.driver.tx,v) for v in m.driver.born['U']]))\n        knots,polys=coefficients(d);row['source_polynomial_intervals']=len(knots)-1\n        del knots,polys")]
    for a,b in changes:assert s.count(a)==1,(a,s.count(a));s=s.replace(a,b)
    s=s.replace('drive_amplitude=AMP','drive_amplitude=incident.ETA')
    ns=dict(prior.source.__globals__,OUT=OUT,coefficients=coefficients,incident=incident);exec(compile(s,__file__,'exec'),ns)
    (OUT/'expanded-source.py').write_text(s);ns['source']()


def fields():
    verified=read(OUT/'driver-polynomial-fixed-audit.json');assert verified['passed']
    for p,h in verified['bindings'].items():assert sha(p)==h,p
    s=inspect.getsource(prior.fields)
    a="trace=d['radau_baryon_g']*LD(d['cx'])*LD(C)**2+d['radau_nonrest_trace_erg']"
    b="knots,co=coefficients(d);trace=co['baryon_g']*LD(d['cx'])*LD(C)**2+co['nonrest_trace_erg']"
    assert s.count(a)==1;s=s.replace(a,b)
    a="trace_poly=poly(d['t'],cubic(trace.swapaxes(0,1)))";assert s.count(a)==1;s=s.replace(a,"trace_poly=PPoly(np.asarray(trace[::-1],float),knots)")
    s=s.replace("ns['direct'](m,d,8)","ns['direct'](m,dict(d,t=knots),8)")
    ns=dict(prior.fields.__globals__,OUT=OUT,Response=Response,coefficients=coefficients,PPoly=PPoly)
    exec(compile(s,__file__,'exec'),ns);(OUT/'expanded-fields.py').write_text(s);ns['fields']()


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));prior.prior.evolution.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();error=None
    try:
        if action!='prepare':
            for p,h in dict(read(OUT/'plan.json')['bindings'],**read(OUT/'plan.json')['reused']).items():
                actual=OUT/'prepared-producer.py' if Path(p).resolve()==Path(__file__).resolve() else p
                assert sha(actual)==h,p
        if action=='source_retry':
            assert 'not subscriptable' in read(OUT/'source-receipt.json')['error']
            write(OUT/'dispatch-amendment.json',dict(change='Accept the same drive duration as a Python float before serialization or a scalar ndarray after loading. No physical equation or gate changes.',prepared_sha256=sha(OUT/'prepared-producer.py'),executed_sha256=sha(__file__),retry_cap_seconds=CAPS[action]))
        globals()['source' if action=='source_retry' else action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
