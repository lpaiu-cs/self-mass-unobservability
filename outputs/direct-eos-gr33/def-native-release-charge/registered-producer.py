"""Direct retarded scalar response of the saved conservative gas release.

Counterexample candidate: prescribed nonlinear gas on the saved metric, an
incoming local acoustic bulk response, and the direct (potential-free) Green
term. This is not the full coupled Einstein/scalar/matter solution.
"""
from pathlib import Path
import argparse
import json
import signal
import time
import numpy as np
from scipy.interpolate import PPoly, CubicSpline
import def_native_radial_thermo as prior

old=prior.task.old
OUT=prior.OUT.parent/'def-native-release-charge'
C=old.C
G=old.G
write=old.write


def prepare():
    assert not OUT.exists();OUT.mkdir()
    paths=[Path(__file__),Path(prior.__file__),Path(prior.task.__file__),
           prior.OUT/'eos.npz',prior.OUT/'result.json',prior.OUT/'flow-audit.json',
           old.prior.OUT/'background.npz',old.prior.OUT/'final-envelope.npz']
    paths += [prior.OUT/f'cells-{n}.{ext}' for n in [896,1792] for ext in ['npz','json']]
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='30c4641ab',
        claim='Apply the saved radial gas histories and opposite inner baryon flux to the direct retarded scalar Green term. Resolve rest-mass cancellation before summation and normalize its outgoing coefficient by the saved ADM mass.',
        model='Fixed native Jordan metric, nonlinear saved gas, local incoming linear acoustic bulk response. Exact background source weight and radial light delay. No scalar/metric potential iteration or gas feedback is claimed.',
        time='Same saved3.434ms horizon; readout stops before any radius would require a future source value. No center-return wave is in this causal window.',
        controls=['896/1792 saved gas paths','piecewise-linear versus natural-cubic history','12/24 Gauss points per acoustic history interval','manufactured conservative retarded source'],
        gates=dict(wave_grid_relative=.02,wave_history_relative=.02,acoustic_quadrature_relative=1e-7,manufactured_relative=1e-8),
        budget=dict(seconds=120,CPU_threads=1,memory_GB=2,native_calls=4,new_fluid_steps=0,new_whole_star_steps=0),
        tail='Unresolved mass stays nonnegative and within the old surface to+1200m, moves at<=0.001c, and has0<=u<=max native cut u and0<=p/rho<=max native cut p/rho. Conditional premises, not an EOS certificate.',
        stop='Preserve failures; no new fluid resolution, time horizon, EOS domain or relaxed gates. Report numerical comparisons separately from missing physical equations.',
        bindings={str(p.relative_to(old.ROOT)):old.photons.digest(p) for p in paths}))


def polynomial(t,y,kind='linear'):
    y=np.asarray(y,float)
    if kind=='cubic':return CubicSpline(t,y,axis=0,bc_type='natural',extrapolate=False)
    slope=np.diff(y,axis=0)/np.diff(t).reshape((-1,)+(1,)*(y.ndim-1))
    return PPoly(np.stack([slope,y[:-1]]),t,extrapolate=False)


def paired(poly,t):
    """Evaluate each cell's polynomial only at that cell's query time."""
    t=np.asarray(t,float);clipped=np.clip(t,0,poly.x[-1]);ids=np.clip(np.searchsorted(poly.x,clipped,side='right')-1,0,len(poly.x)-2)
    dt=clipped-poly.x[ids];value=np.zeros_like(t)
    for coeff in poly.c:
        c=coeff[ids] if coeff.ndim==1 else coeff[ids,np.arange(coeff.shape[-1])]
        value=value*dt+c
    assert np.max(t)<=poly.x[-1]+2e-15,'Future source requested'
    return np.where(t<=0,0.,value)


class Geometry:
    def __init__(self,m):
        self.m=m;self.bg=m.bg;self.gx,self.gw=np.polynomial.legendre.leggauss(16)
        self.phis=float(self.bg.d['phi'][-1])

    def metric(self,x):
        x=np.asarray(x);re=(self.m.RJ+x)/(self.m.As*self.m.R)
        for _ in range(3):
            p=self.bg.sample(re);phi=p['phi'].copy();outside=re>1
            # sample() retains the surface phi outside. Integrate the actual
            # saved vacuum gradient instead of extending that constant.
            rr=1+(re[outside,None]-1)*(self.gx+1)/2
            if rr.size:phi[outside]=self.phis+(re[outside]-1)/2*(self.bg.metric(rr.ravel())[-1].reshape(rr.shape)@self.gw)
            A=np.exp(-2*phi*phi);alpha=-4*phi;den=1+re*alpha*p['v']
            re-=(A*self.m.R*re-(self.m.RJ+x))/(A*self.m.R*den)
        a=A*p['N'];b=1-2*p['m']/re;B=1/(np.sqrt(b)*den)
        w=alpha*a/(re*self.m.R)
        return w,a,B,re,phi

    def __call__(self,x):
        x=np.asarray(x);w,a,B,re,phi=self.metric(x)
        xx=x[:,None]*(self.gx+1)/2;_,aa,bb,_,_=self.metric(xx.ravel())
        delay=x/(2*C)*((bb/aa).reshape(xx.shape)@self.gw)
        return w,delay,a,B,re


def native_bulk(m):
    fan=prior.task.prior.Fan(call_cap=4,reuse=True);env=m.env
    x=-20000-np.array([0.,2000.,4000.,8000.]);r=m.R+x/m.As
    rho=np.exp(np.interp(r,env['r'],np.log(env['rho'])));T=np.exp(np.interp(r,env['r'],np.log(env['T'])))
    raw=np.array([fan.call(np.log(rr),np.log(tt)) for rr,tt in zip(rho,T)])
    np.savez_compressed(OUT/'native-bulk.npz',x_cm=x,rho=rho,T=T,raw=raw)
    return raw


def acoustic(port,geom,m,raw,utimes,order):
    rr,p,u,gamma=raw[0,[0,1,2,4]];cx=m.eos.cx
    cs=np.sqrt(gamma*p/(rr*(cx*C*C+u)+p));w0=geom(np.array([0.]))[0][0]
    wi,di,ai,Bi,re=geom(np.array([-20000.]));di=float(di[0]);ai=float(ai[0]);Bi=float(Bi[0])
    speed=ai/Bi*cs*C;kbulk=u+p/rr-3*gamma*p/rr
    gx,gw=np.polynomial.legendre.leggauss(order);H=port.antiderivative();rest=[];nr=[]
    for t in utimes:
        upper=max(0.,(t+di)/(1+cs));cuts=np.unique(np.r_[0,upper,(t+di-port.x)/(1+cs)])
        cuts=cuts[(cuts>=0)&(cuts<=upper)]
        if len(cuts)>1:
            tau=(cuts[:-1,None]+np.diff(cuts)[:,None]*(gx+1)/2).ravel();weight=(np.diff(cuts)[:,None]*gw/2).ravel()
            w=geom.metric(-20000-speed*tau)[0];B=-paired(port,t+di-(1+cs)*tau)
            varying=np.sum(weight*(w-w0)*B,dtype=np.longdouble)
            nrt=np.sum(weight*w*B,dtype=np.longdouble)*kbulk
        else:varying=nrt=0.
        # Integral of the local acoustic profile is exactly -port(t). The
        # expression below subtracts that common-time baryon history before
        # multiplying by the large rest energy.
        delayed=-paired(H,np.array(t+di))/(1+cs)+paired(H,np.array(t))
        rest.append(-G*cx/(2*C)*(varying+w0*delayed));nr.append(-G/(2*C**3)*nrt)
    hnr=ai*(u+p/rr)+(ai-m.a0)*cx*C*C
    return np.array(rest,float),np.array(nr,float),dict(cs_over_c=float(cs),coordinate_speed_cm_s=float(speed),
        trace_nonrest_erg_per_g=float(kbulk),Killing_nonrest_erg_per_g=float(hnr),bulk_depth_cm=float(speed*utimes[-1]))


def readout(n,kind,utimes,raw,order=12):
    m=prior.Flow(n);geom=Geometry(m);saved=np.load(prior.OUT/f'cells-{n}.npz');h=saved['history'];t=np.r_[0,h[:,0]]
    states=np.concatenate([saved['initial'][None],saved['snapshots']]);scale=4*np.pi*m.RJ**2*m.eos.rho0;w,d,a,B,re=geom(m.x);w0=geom(np.array([0.]))[0][0]
    baryon=(states[:,0].astype(np.longdouble)-states[0,0])*m.vol*scale
    nr=[]
    for U in states:
        rho,v,sigma=m.primitive(U);p,u,*_=m.eos(rho,sigma)
        # rho-D=-D*v^2/(1+sqrt(1-v^2)), evaluated without cancellation.
        trace=-m.eos.cx*U[0]*v*v/(1+np.sqrt(1-v*v))+rho*u-3*p
        nr.append(trace*m.vol*scale*C*C)
    nr=np.array(nr);nr-=nr[0]
    port_values=np.r_[0,h[:,3]]*scale;port=polynomial(t,port_values,kind)
    # Conservation infers only the cumulative missing mass, not its position.
    missing=port_values-np.sum(baryon,axis=1,dtype=np.longdouble)
    tail_roundoff=float(max(0,-np.min(missing)))
    assert tail_roundoff<.01,'Negative unresolved mass exceeds roundoff'
    HB=polynomial(t,baryon,kind).antiderivative();HN=polynomial(t,nr,kind).antiderivative()
    qr=[];qn=[]
    for tt in utimes:
        at=np.full(n,tt);now=paired(HB,at);later=paired(HB,at+d)
        centered=np.sum(w*(later-now)+(w-w0)*now,dtype=np.longdouble)
        qr.append(-G*m.eos.cx/(2*C)*centered)
        qn.append(-G/(2*C**3)*np.sum(w*paired(HN,at+d),dtype=np.longdouble))
    br,bn,bulk=acoustic(port,geom,m,raw,utimes,order)
    qr=np.array(qr,float);qn=np.array(qn,float);Q=qr+qn+br+bn
    mass=m.bg.M*m.R;alpha=-Q/mass
    energy_bulk=-np.r_[0,h[:,4]-h[:,5]]*scale*C*C
    acoustic_energy=-port_values*bulk['Killing_nonrest_erg_per_g'];energy_mismatch=energy_bulk-acoustic_energy
    filename=f'cells-{n}-{kind}-g{order}'
    np.savez_compressed(OUT/f'{filename}.npz',u_seconds=utimes,charge_cm=Q,normalized_charge=alpha,
        layer_rest_charge_cm=qr,layer_nonrest_charge_cm=qn,bulk_rest_charge_cm=br,bulk_nonrest_charge_cm=bn,
        source_times=t,unresolved_baryon_g=np.asarray(missing,float),bulk_baryon_g=-port_values,
        bulk_Killing_nonrest_energy_erg=energy_bulk,acoustic_Killing_nonrest_energy_erg=acoustic_energy,
        acoustic_energy_mismatch_erg=energy_mismatch)
    row=dict(classification='Counterexample candidate',cells=n,history=kind,acoustic_order=order,
        endpoint_direct_normalized_charge=float(alpha[-1]),peak_direct_normalized_charge=float(max(abs(alpha))),
        endpoint_components_cm=[float(z[-1]) for z in [qr,qn,br,bn]],
        minimum_unresolved_baryon_g=float(min(missing)),maximum_unresolved_baryon_g=float(max(missing)),
        bulk=bulk,endpoint_bulk_energy_mismatch_erg=float(energy_mismatch[-1]),
        scope='Direct scalar response with prescribed source and local acoustic bulk response; tail center contribution is excluded and conditionally enclosed separately.')
    write(OUT/f'{filename}.json',row)
    return alpha,row


def run():
    assert not (OUT/'result.json').exists();start=time.monotonic();signal.alarm(120)
    plan=json.loads((OUT/'plan.json').read_text())
    for p,h in plan['bindings'].items():assert old.photons.digest(old.ROOT/p)==h,p
    m=prior.Flow(1792);geom=Geometry(m);raw=native_bulk(m);end=float(m.hist_t[-1]);maxdelay=float(geom(np.array([120000.]))[1][0])
    utimes=np.linspace(0,end-maxdelay,129)
    fine,fr=readout(1792,'linear',utimes,raw);pilot=time.monotonic()-start
    write(OUT/'measured-budget.json',dict(measured_native_and_first_readout_seconds=pilot,forecast_total_seconds=5*pilot,
        assumption='Remaining three readouts reuse the same native states; cubic and doubled acoustic quadrature may cost up to twice the first readout.'))
    assert 5*pilot<120,'Readout forecast exceeds registered budget'
    coarse,cr=readout(896,'linear',utimes,raw);cubic,cu=readout(1792,'cubic',utimes,raw);quadrature,qu=readout(1792,'linear',utimes,raw,24)
    peak=max(abs(fine));errors=dict(wave_grid_relative=float(max(abs(coarse-fine))/peak),wave_history_relative=float(max(abs(cubic-fine))/peak),acoustic_quadrature_relative=float(max(abs(quadrature-fine))/peak))
    # A missing positive tail's centered rest source is bounded by its maximum
    # mass, source-time span, weight variation and light-delay Lipschitz bound.
    x=np.linspace(0,120000,65);w,d,_,_,_=geom(x);w0=geom(np.array([0.]))[0][0];mt=fr['maximum_unresolved_baryon_g'];M=m.bg.M*m.R
    tail_rest=G*m.eos.cx/(2*C*M)*mt*(max(abs(w-w0))*end+abs(w0)*max(abs(d)))
    cut=m.eos.d['raw'][:,0];umax=float(max(cut[:,2]));prmax=float(max(cut[:,1]/cut[:,0]));vmax=.001
    tail_specific=umax+3*prmax+m.eos.cx*C*C*(1-np.sqrt(1-vmax*vmax))
    tail_nr=G/(2*C**3*M)*max(abs(w))*mt*end*tail_specific
    # This bounds Bondi mass normalization only in the declared optically thin
    # elastic channel; it is not a bound on metric-driven scalar forcing.
    flow=json.loads((prior.OUT/'cells-1792.json').read_text());tau=flow['maximum_electron_scattering_optical_depth']
    energy_bound=2*float(m.env['Linfinity'])*end*tau*np.exp(tau);alpha0=abs(m.bg.K/m.bg.M)
    mass_bound=alpha0*G*energy_bound/(C**4*M)
    write(OUT/'conditional-tail-and-mass.json',dict(classification='Proven',conditional=True,premises=plan['tail'],
        tail_rest_normalized_bound=float(tail_rest),tail_nonrest_normalized_bound=float(tail_nr),
        tail_speed_bound_over_c=vmax,tail_specific_trace_bound_erg_per_g=tail_specific,
        elastic_photon_energy_redistribution_bound_erg=float(energy_bound),elastic_mass_normalization_bound=float(mass_bound),
        elastic_premise='Same initial photon luminosity; optical depth bounded by the saved maximum in both histories, elastic scattering, no absorption, and all counted photons emitted during this horizon. Not a bound on pre-existing exterior photon populations.',
        whole_charge_enclosed=False))
    variation={name:float(np.max(abs(raw[:,i]/raw[0,i]-1))) for name,i in [('rho',0),('pressure',1),('specific_energy',2),('Gamma1',4)]}
    passed=all(errors[k]<plan['gates'][k] for k in errors)
    result=dict(classification='Counterexample candidate',passed=passed,controls=errors,
        endpoint_direct_normalized_charge=fr['endpoint_direct_normalized_charge'],peak_direct_normalized_charge=fr['peak_direct_normalized_charge'],
        readout_end_seconds=float(utimes[-1]),source_end_seconds=end,maximum_light_delay_seconds=maxdelay,
        local_bulk_native_variation=variation,native_calls=4,new_fluid_steps=0,seconds=time.monotonic()-start,
        outgoing_direct_scalar_computed=True,inner_baryon_acoustic_response_applied=True,inner_energy_residual_retained=True,
        full_inner_boundary_feedback=False,full_GR_scalar_feedback=False,physical_tail_certified=False,
        final_charge_solved=False,full_goal_complete=False)
    write(OUT/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
