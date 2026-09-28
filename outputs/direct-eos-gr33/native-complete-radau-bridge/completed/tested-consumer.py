"""Read the accepted common photon/material history with its dense GR source.

Counterexample candidate. No new physical trajectory, source fit, relaxed gate,
or final charge claim. The declared final interval and self-GR remain separate.
"""
from pathlib import Path
from types import FunctionType,SimpleNamespace
from decimal import Decimal,localcontext
import gc,inspect,json,os,resource,shutil,sys,time
import numpy as np
from scipy.interpolate import PPoly
import apply_driver_aware_radau_gr as prior
import continue_captured_photon_history as recovery

OUT=Path('native-complete-radau224-work');INPUT=recovery.OUT
read,write,sha,LD,AMP,C=prior.read,prior.write,prior.sha,prior.LD,prior.AMP,prior.C
base=prior.base;endpoint=prior.prior.prior;KEYS=prior.KEYS
CAPS=dict(check=300,prepare=180,endpoints=1800,source=3600,fields=3600)
saved=recovery.prior.prior.saved


def background_contrasts(m,E,P):
    """Evaluate the same affine energy contrast before cancellation, at 80 digits."""
    def dec(v):
        a,b=np.longdouble(v).as_integer_ratio();return Decimal(a)/Decimal(b)
    rows=[]
    with localcontext() as ctx:
        ctx.prec=80;rest=dec(m.rest)
        for j in range(17):
            q=m.point(j)['Q']
            rows.append([LD(str(dec(E[i])+dec(P[i])-(dec(q[2,i])+rest*dec(q[0,i]))/dec(m.a[i]))) for i in range(m.n)])
    return np.array(rows,LD)


def replace(source,old,new):
    assert source.count(old)==1,(old,source.count(old))
    return source.replace(old,new)


def coefficients(d):
    source=inspect.getsource(prior.coefficients)
    source=replace(source,"d['drive_times'],-x/C","d['drive_times'],d['geometry_times'],-x/C")
    source=replace(source,"a,b=d['geometry_coeff_'+key];g=a[None]+left[:,None]*b",
        "a,b=d['geometry_coeff_'+key];gt=d['geometry_times'];ids=np.clip(np.searchsorted(gt,left,side='right')-1,0,len(gt)-2)\n            b=b[ids];g=a[ids]+(left-gt[ids])[:,None]*b")
    source=replace(source,'pulse*b[None,None]','pulse*b[None]')
    ns=dict(prior.coefficients.__globals__);exec(compile(source,__file__,'exec'),ns)
    return ns['coefficients'](d)


class Response(prior.Response):
    setup=FunctionType(prior.Response.setup.__code__,dict(prior.Response.setup.__globals__,coefficients=coefficients))


def geometry(d,t,key):
    gt=d['geometry_times'];ids=np.clip(np.searchsorted(gt,t,side='right')-1,0,len(gt)-2)
    a,b=d['geometry_coeff_'+key]
    return a[ids]+(t-gt[ids])[:,None]*b[ids]


def check():
    # Regression against the accepted short-interval consumer, without a solve.
    errors={}
    for n in [64,128]:
        d=dict(np.load(prior.OUT/'gr'/f'source-{n}.npz'));old_knots,old=prior.coefficients(d)
        d['geometry_times']=np.array([0,d['t'][-1]])
        for key in KEYS:d['geometry_coeff_'+key]=d['geometry_coeff_'+key][:,None]
        knots,new=coefficients(d);assert np.array_equal(knots,old_knots)
        assert all(np.array_equal(old[k],new[k]) for k in KEYS)
        # Split the same affine background; its polynomial must remain equal.
        T=d['t'][-1];d['geometry_times']=np.array([0,T/2,T])
        for k in KEYS:
            a,b=d['geometry_coeff_'+k][:,0];d['geometry_coeff_'+k]=np.array([[a,a+T/2*b],[b,b]])
        knots,new=coefficients(d);times=(knots[:-1]+knots[1:])/2;values={}
        for k in KEYS:
            a=PPoly(np.asarray(old[k][::-1],float),old_knots)(times)
            b=PPoly(np.asarray(new[k][::-1],float),knots)(times)
            values[k]=float(np.max(abs(a-b))/max(np.max(abs(a)),1e-290))
        assert max(values.values())<1e-12,values;errors[n]=values
    import sympy as sp
    e,p,a,b,w=sp.symbols('e p a b w')
    assert sp.expand(e+p-((1-w)*a+w*b)-((1-w)*(e+p-a)+w*(e+p-b)))==0
    result=dict(classification='Counterexample candidate',passed=True,original_single_interval_coefficients_exact=True,
        affine_background_contrast_identity=True,
        split_affine_errors=errors,scope='Saved212consumer regression and an affine partition identity, not full-history acceptance.',source_sha256=sha(__file__))
    print(json.dumps(result),flush=True);return result


def prepare():
    result=read(INPUT/'result.json');assert result['same_solution_accepted_history_recovered']
    assert read(INPUT/'controller-status.json')['state']=='completed'
    assert not OUT.exists();OUT.mkdir();files=[]
    for folder in ['sweep-0','sweep-1/photons','sweep-1/material','gr']:(OUT/folder).mkdir(parents=True)
    for src in list((INPUT/'sweep-0').rglob('*.npz'))+[INPUT/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
        dst=OUT/src.relative_to(INPUT);dst.parent.mkdir(parents=True,exist_ok=True);os.link(src,dst);files += [src,dst]
    horizons=[]
    for n in [64,128]:
        p=np.load(saved(n));q=np.load(INPUT/f'recovered-{n}.npz')
        assert np.array_equal(p['joint_stage_times'],q['times']) and np.array_equal(p['joint_stage_weights'],q['weights'])
        assert len(q['photon_moments'])==2*(len(p['actual_step_edges'])-1);horizons.append(float(p['actual_step_edges'][-1]))
        files += [saved(n),INPUT/f'recovered-{n}.npz',INPUT/f'recovered-{n}.json']
    assert abs(horizons[0]-horizons[1])<1e-18
    files += [INPUT/n for n in ['result.json','controller-status.json','coarse-receipt.json','fine-receipt.json','audit-receipt.json']]
    files += [prior.OUT/n for n in ['expanded-source.py','expanded-fields.py','driver-polynomial-fixed-audit.json','source_retry-receipt.json']]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Apply the complete accepted common15/16photon/material history to its own dense GR source and retarded field. Judge the unchanged time gate before physical return.',
        method='Reuse206endpoint map and212actual degree8pulse plus Born term and Radau state. Retain the original piecewise-affine background on all canonical intervals instead of extrapolating interval0. Restore conserved gas exactly, use the same accepted floor, photon moments, energy and radial ports. No added diagnostic charge.',
        gates=dict(dense_stage=1e-12,polynomial=1e-12,pressure=.002,mapping=1e-12,ledger=1e-8,source_time=.02,field_time=.02,quadrature=.002,independent_GR=1e-9),
        original_return_horizon_seconds=horizons[0],budgets=CAPS,CPU_threads=1,virtual_GiB=12,
        forecast='212source6totalsteps45.81s including setup.326existingsteps and16background intervals have unmeasured marginal cost; allow30minutes endpoints and60minutes dense source. Characteristic source cuts and larger time arrays increase cost; allow60minutes fields.12GiB leaves room for the degree9source arrays. No physical reintegration.',
        stop='Any provenance, mapping, dense-source or independent numerical gate fails. Source/field time failures are preserved even when diagnostic field evaluation completes. No automatic extra grid, interval, physical path or relaxed gate.',
        full_declared_period=False,missing_final_interval=True,self_GR_closed=False,final_charge_conclusion='unadjudicated',full_goal_complete=False,
        bindings={str(p):sha(p) for p in dict.fromkeys(files)}))
    write(OUT/'symbolic.json',prior.symbolic());write(OUT/'regression.json',check())


def aligned(a,b,key):
    ids=np.array([np.argmin(abs(b['t']-t)) for t in a['t']]);assert np.max(abs(a['t']-b['t'][ids]))<1e-18
    # The original source criterion is spatial L1, maximized over common times.
    return float(np.max(np.sum(abs(a[key]-b[key][ids]),axis=-1))/max(np.max(np.sum(abs(b[key][ids]),axis=-1)),LD('1e-290')))


def recovered(n):
    p=np.load(INPUT/f'recovered-{n}.npz');return len(p['times'])//2,p['photon_moments'],p['radial_ports']


def endpoints():
    source=inspect.getsource(endpoint.source)
    marker="        b=model.model.bulk;weights="
    source=replace(source,marker,"        contrasts=background_contrasts(material,Eg0,Pg0)\n"+marker)
    source=replace(source,'(Eg0+Pg0-backgroundE)','((1-w)*contrasts[k]+w*contrasts[k+1])')
    source=replace(source,"g=np.column_stack([(q[2]-m.kappa*q[0])/m.eu,q[3]/m.nu,q[0]/m.bu,q[1]/m.su]);off=",'g=restored_gas(m,q);off=')
    source=replace(source,";assert result['passed'],result",'')
    initialize=FunctionType(endpoint.initialize.__code__,dict(endpoint.initialize.__globals__,OUT=OUT))
    ns=dict(endpoint.source.__globals__,OUT=OUT,initialize=initialize,recovered=recovered,aligned=aligned,
        restored_gas=recovery.prior.restored_gas,background_contrasts=background_contrasts,base=SimpleNamespace(**dict(vars(base),saved=saved)))
    exec(compile(source,__file__,'exec'),ns);(OUT/'expanded-endpoints.py').write_text(source);ns['source']()
    os.rename(OUT/'sources.json',OUT/'endpoint-sources.json')
    for n in [64,128]:
        os.rename(OUT/f'source-{n}-check.json',OUT/f'endpoint-{n}-check.json')
        shutil.copyfile(OUT/'gr'/f'source-{n}.npz',OUT/'gr'/f'endpoint-{n}.npz')


def source():
    s=(prior.OUT/'expanded-source.py').read_text()
    changes=[('count=n//32',"count=len(p['actual_step_edges'])-1"),
        ("        b=model.model.bulk;weights=","        contrasts=background_contrasts(material,Eg0,Pg0)\n        b=model.model.bulk;weights="),
        ('(Eg0+Pg0-background)','((1-w)*contrasts[k]+w*contrasts[k+1])'),
        ("(prior.REC if n==64 else recovery.OUT)",'INPUT'),
        ("OLD/'gr'/f'source-{n}.npz'","OUT/'gr'/f'endpoint-{n}.npz'"),
        ('for k in [0,1]}','for k in range(17)}'),(';assert k==0',''),
        ("gas=lambda q:np.column_stack([(q[2]-m.kappa*q[0])/m.eu,q[3]/m.nu,q[0]/m.bu,q[1]/m.su])",'gas=lambda q:restored_gas(m,q)'),
        ("first=readout(0.,zero,zero_m,zero_p,unit,False);last=readout(T,zero,zero_m,zero_p,unit,False)\n        geometry={k:np.array([first[k],(last[k]-first[k])/T]) for k in KEYS}",
         "gt=np.unique(np.r_[m.t[m.t<T],T]);gg=[readout(v,zero,zero_m,zero_p,unit,False) for v in gt]\n        geometry={k:np.array([[a[k] for a in gg[:-1]],[(b[k]-a[k])/(v-u) for a,b,u,v in zip(gg[:-1],gg[1:],gt[:-1],gt[1:])]]) for k in KEYS}"),
        ('(geometry[k][0]+(t+h*u)*geometry[k][1])',"geometry_at(geometry,gt,np.array([t+h*u]),k)[0]"),
        ("d.update(drive_x=","d.update(geometry_times=gt,drive_x="),
        ('drive_amplitude=AMP','drive_amplitude=prior.incident.ETA'),
        ("old=read(OLD/'original-source-norm.json')","old=read(OUT/'endpoint-sources.json')"),
        ("old['source_time_spatial_L1']","old['source_time']")]
    for a,b in changes:s=replace(s,a,b)
    def geometry_at(g,gt,t,key):
        ids=np.clip(np.searchsorted(gt,t,side='right')-1,0,len(gt)-2);a,b=g[key]
        return a[ids]+(t-gt[ids])[:,None]*b[ids]
    ns=dict(prior.prior.source.__globals__,OUT=OUT,INPUT=INPUT,prior=endpoint,
        base=SimpleNamespace(**dict(vars(base),saved=saved)),
        coefficients=coefficients,restored_gas=recovery.prior.restored_gas,geometry_at=geometry_at,background_contrasts=background_contrasts)
    # The generated source uses the live driver constant, not response AMP.
    ns['incident']=prior.incident;s=s.replace('drive_amplitude=prior.incident.ETA','drive_amplitude=incident.ETA')
    exec(compile(s,__file__,'exec'),ns);(OUT/'expanded-source.py').write_text(s);ns['source']()
    m=base.run.owner.Model(64);rows=[]
    for n in [64,128]:
        d=dict(np.load(OUT/'gr'/f'source-{n}.npz'));knots,co=coefficients(d);times=(knots[:-1]+knots[1:])/2
        phi=np.array([m.driver.wave(float(t),m.driver.xc)[0]/(m.driver.centers*AMP) for t in times]);errors={}
        for k in KEYS:
            actual=PPoly(np.asarray(co[k][::-1],float),knots)(times)
            expected=prior.prior.poly(d['t'],d['state_coeff_'+k])(times)
            if k not in KEYS[-2:]:expected+=geometry(d,times,k)*phi
            measure=lambda v:np.sum(abs(v),axis=-1) if v.ndim>1 else abs(v)
            errors[k]=float(np.max(measure(actual-expected))/max(np.max(measure(expected)),1e-290))
        rows.append(dict(clock=n,interval_midpoints=len(times),errors=errors));assert max(errors.values())<1e-12,rows[-1]
    write(OUT/'driver-polynomial-audit.json',dict(classification='Counterexample candidate',passed=True,rows=rows))


def fields():
    assert read(OUT/'sources.json')['representation_controls_passed'] and read(OUT/'driver-polynomial-audit.json')['passed']
    s=(prior.OUT/'expanded-fields.py').read_text()
    old="old=dict(np.load(Path('native-early-gr207-work')/'gr/fields-128-g8.npz'))\n    change={key:float(np.max(abs(fine[key]-old[key]))/max(np.max(abs(fine[key])),1e-290)) for key in ['U','U_t','U_x']}"
    s=replace(s,old,'change={} # Earlier short-horizon arrays are not this full-history comparison.')
    #212finish reused three precomputed fields. This new history computes its own.
    s=replace(s,"rows=[read(OUT/'gr'/f'fields-{n}-g{q}.json') for n,q in [(128,8),(64,8),(128,4)]]","rows=[fn(m,n,q) for n,q in [(128,8),(64,8),(128,4)]]")
    ns=dict(prior.prior.fields.__globals__,OUT=OUT,Response=Response,coefficients=coefficients,PPoly=PPoly)
    exec(compile(s,__file__,'exec'),ns);(OUT/'expanded-fields.py').write_text(s);ns['fields']()
    r=read(OUT/'result.json');r.update(accepted_common_history_applied=True,physical_horizon_seconds=read(OUT/'plan.json')['original_return_horizon_seconds'],
        full_declared_period=False,missing_final_interval=True,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'result.json',r)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS
    resource.setrlimit(resource.RLIMIT_AS,(12*1024**3,12*1024**3));endpoint.evolution.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();error=None
    try:
        if action not in ['check','prepare']:
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(OUT/f'{action}-receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
