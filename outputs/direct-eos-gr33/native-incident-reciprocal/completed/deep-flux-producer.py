import resource,time
import numpy as np
import solve_native_incident_reciprocal as solve
start=time.monotonic();cpu=time.process_time();error=None;out=solve.OUT;LD=np.longdouble
solve.base.drive.native.deadline(30);resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
try:
    solve.initialize(2);m=solve.Material(128,128);rows=[]
    for k in [0,8,16]:
        q=m.point(k)['Q'].astype(LD);shape=1+np.sin(np.arange(m.n)*.37)/3
        z=q*.001*shape;z[1]=q[0]*m.model.cx*solve.base.C**2*.0002*shape
        z[2]=((m.a-m.model.m.a0)*m.model.cx*solve.base.C**2*z[0]+abs(q[2]-(m.a-m.model.m.a0)*m.model.cx*solve.base.C**2*q[0])*.001*shape)
        field=np.array([shape*.002,shape*.001,shape*-.001,shape*.0004,m.a/m.R*shape*.001])
        exact=m.deep_tangent(k,z,field);zero=m.raw(k,np.zeros_like(z),np.zeros_like(field),0.)
        eps=LD('1e-5');a=m.raw(k,z,field,eps);b=m.raw(k,z,field,eps/2)
        derivative=[2*(v.astype(LD)-w.astype(LD))/(eps/2)-(u.astype(LD)-w.astype(LD))/eps for u,v,w in zip(a[:2],b[:2],zero[:2])]
        measured=[derivative[0][:,:m.nb],derivative[1][:m.nb]]
        error_flux=(np.sum(abs(measured[0]-exact[0]),axis=1)/np.maximum(np.sum(abs(exact[0]),axis=1),1.)).astype(float).tolist()
        error_gravity=float(np.sum(abs(measured[1]-exact[1]))/max(np.sum(abs(exact[1])),1.))
        row=dict(k=k,deep_flux_components=error_flux,deep_gravity=error_gravity);rows.append(row)
        assert max(error_flux+[error_gravity])<.002,row
    solve.write(out/'deep-flux-check.json',dict(classification='Counterexample candidate',passed=True,rows=rows,
        scope='Independent positive Richardson probes of the actual deep flux and volume force, normalized separately without dilution by atmospheric forcing. Combined mass/momentum/thermal/H and geometry direction at three saved backgrounds; not uniform physical-EOS certification.'))
    print(rows,flush=True)
except Exception as exc:error=repr(exc);raise
finally:solve.write(out/'deep-flux-check-receipt.json',dict(seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,error=error,source_sha256=solve.sha(__file__)))
