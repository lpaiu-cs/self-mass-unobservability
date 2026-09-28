"""Measure the actual finite collision defect of the accepted radiation state."""
from pathlib import Path
import json,resource,signal,sys,time
import numpy as np
import def_native_finite_motion_feedback as run

OUT=run.OUT/'finite-collision';write=run.write;sha=run.sha;AMP=run.AMP;LD=np.longdouble


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        claim='Replace the unmeasured radiation linearization remainder by the actual finite collision difference on the saved simultaneous state, preserving its prescribed baryon/momentum and solved energy/H; determine whether a nonlinear correction is needed before another charge claim.',
        scope='Same declared native-interpolation EOS and moving packet owner. Fixed native-minus-table forcing cancels from this comparison; its derivative error relative to the underlying microscopic EOS is separate and still open.',
        exact_state='For each saved time use the prescribed Phase141 B/S/inventory but the actual Phase142 simultaneous E/H and angular photons. Do not silently replace prescribed B/S by the later free-material output.',
        reuse='All saved17 states,531 cells,8 angles,152 frequencies. No trajectory, native inverse, EOS table or coefficient bank replay.',
        budgets=dict(pilot_seconds=65,production_seconds=180,CPU_threads=1,virtual_GiB=3),
        forecast='Measure0/8/16 point costs, reuse them, and require2x maximum cost times14 remaining plus10s under180s. Do not run a response until the exported error is resolved and a measured response plan is registered.',
        method='Compute coefficient differences before applying them to the background occupation, retain the new coefficient acting on the actual occupation increment, and subtract the stored linear tangent response. Evaluate half amplitude to distinguish a finite remainder from unresolved arithmetic.',
        gates=dict(base_owner=1e-9,material_support=True,photon_negativity_energy=1e-12,finite_difference_resolution=.002,number=1e-10),
        stop='Stop on positivity, coefficient-owner disagreement or unresolved difference; repair that owner before propagation. No clipping, amplitude reduction as a replacement result, tolerance relaxation, grid expansion or automatic full-response replay.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(run.__file__),Path(run.base.__file__),run.OUT/'audit.json',run.photon_path(128,128),run.matter.OUT/'steps-128-reference-128.npz']}))


class State:
    def __init__(self):
        run.configure();self.m=run.Response();self.p=dict(np.load(run.photon_path(128,128)))
        self.m.model=self.m.material.model

    def coefficients(self,k,z,factor):
        m=self.m;mat=m.material;raw=mat.raw(k,z,np.zeros((5,m.n)),AMP*factor)[3]
        m.theta=raw['theta'];m.eta=raw['eta'];m.rho,m.beta,m.lt,m.y=raw['primitive']
        m.bulk_beta=m.model.velocity();m.bulk_x=m.model.bulk.eos.x.copy();m.active=m.rho>=m.model.flow.eos.floor
        return m.coefficients()

    def point(self,k):
        m=self.m;p=self.p;start=time.monotonic();t=m.t[k];assert p['t'][k]==t
        # Evaluate the tangent at exactly the state used by the accepted solve.
        c=run.base.Response.local(m,t);x=p['photon_history_scaled_occupation'][k]/(m.scale*AMP);g=p['material_history'][k]/AMP
        linear,_,escape,bound=m.collision(c,x,g,True)
        linear=linear*m.scale*AMP;bound=bound*m.scale*AMP;escape=escape*AMP
        z=m.motion[k].copy();z[[2,3]]=p['moments'][k,[1,2]]/AMP
        zero=np.zeros_like(z);a=self.coefficients(k,zero,0.);sa,ea=m.scattering_matrix(a)
        bank=dict(np.load(run.prior.OUT/f'bank-128/point-{k}.npz'))
        owner=max(float(np.max(abs(a[q]-bank[q]))/max(np.max(abs(bank[q])),1e-300)) for q in ['emit','loss','sc','rho','beta','u','p'])
        I=m.I[k];delta=p['photon_history_scaled_occupation'][k];X=I/m.scale;dx=delta/m.scale
        rows=[];full=[];floors=[]
        for factor in [1.,.5]:
            b=self.coefficients(k,z,factor);sb,eb=m.scattering_matrix(b)
            db=(b['emit'].astype(LD)-a['emit'])-(b['loss'].astype(LD)-a['loss'])*I-b['loss'].astype(LD)*delta*factor
            scatter=((sb-sa)@X.ravel()).reshape(X.shape)*m.scale+(sb@(dx*factor).ravel()).reshape(X.shape)*m.scale
            de=np.einsum('knqf,nqf->kn',eb-ea,X)+np.einsum('knqf,nqf->kn',eb,dx*factor)
            dp=db+scatter
            remainder=np.asarray(dp-factor*linear,float);br=np.asarray(db-factor*bound,float);er=de-factor*escape
            weight=m.Eweight/m.scale;N=m.Nweight/m.scale
            norm=max(np.sum(abs(factor*linear)*weight),1.)
            absolute=float(np.sum(abs(remainder)*weight));relative=absolute/norm
            rounding=16*np.finfo(float).eps*np.sum((abs(a['emit'])+abs(b['emit'])+(abs(a['loss'])+abs(b['loss']))*abs(I))*weight)
            number=np.sum((remainder-br)*N,axis=(1,2))+er[0]
            ns=np.maximum(np.sum((abs(remainder)+abs(br))*N,axis=(1,2))+abs(er[0]),1.)
            negative=float(np.sum(np.maximum(-(I+factor*delta),0)*weight)/max(np.sum(abs(I)*weight),1.))
            rows.append(dict(factor=factor,remainder_over_linear=float(relative),energy_weighted_L1=absolute,rounding_over_linear=float(rounding/norm),
                number_relative=float(np.max(abs(number)/ns)),negative_photon_energy_relative=negative,negative_packet_count=int(np.sum(I+factor*delta<0))))
            full.append((remainder,br,er));floors.append(float(rounding))
        ratio=rows[1]['energy_weighted_L1']/max(rows[0]['energy_weighted_L1'],1.)
        row=dict(classification='Counterexample candidate',k=k,t=float(t),base_owner_relative=owner,rows=rows,half_over_full=ratio,seconds=time.monotonic()-start,
                 passed=bool(owner<1e-9 and max(r['negative_photon_energy_relative'] for r in rows)<1e-12 and max(r['rounding_over_linear'] for r in rows)<.002 and max(r['number_relative'] for r in rows)<1e-10))
        np.savez_compressed(OUT/f'point-{k}.npz',t=t,photon=full[0][0],bound=full[0][1],escape=full[0][2],half_photon=full[1][0],
                            linear=linear,linear_bound=bound,linear_escape=escape)
        write(OUT/f'point-{k}.json',row);print(json.dumps(row),flush=True);return row


def pilot():
    started=time.monotonic();signal.alarm(65);s=State();rows=[]
    for k in [0,8,16]:
        rows.append(s.point(k))
        if not rows[-1]['passed']:break
    forecast=2*max(r['seconds'] for r in rows)*14+10
    p=dict(classification='Counterexample candidate',rows=rows,forecast_upper_seconds=forecast,
        eligible=len(rows)==3 and all(r['passed'] for r in rows) and forecast<180,seconds=time.monotonic()-started)
    write(OUT/'pilot.json',p);signal.alarm(0);print(json.dumps(p),flush=True)


if __name__=='__main__':
    cap=3*1024**3;resource.setrlimit(resource.RLIMIT_AS,(cap,cap));signal.signal(signal.SIGALRM,run.native.forcing.history.flow.old.optical.timeout)
    action=sys.argv[1];start=time.monotonic()
    try:globals()[action]()
    except Exception as exc:write(OUT/f'{action}-failure.json',dict(error=repr(exc),seconds=time.monotonic()-start));raise
