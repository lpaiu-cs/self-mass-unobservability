"""Independent trapezoidal integration of the actual saved gas history."""
import json
import time
import numpy as np
import def_native_release_charge as task


def main():
    start=time.monotonic();out=task.OUT;assert not (out/'history-audit.json').exists()
    task.write(out/'history-audit-plan.json',dict(classification='Counterexample candidate',seconds=30,native_calls=0,new_evolution_steps=0,
        control='Integrate the original piecewise-linear baryon histories by split trapezoids, without producer antiderivatives; compare the actual endpoint. Measure initial causal leakage and the body-only quadrature difference.',
        numerical_gate=1e-8,producer_sha256=task.old.photons.digest(out/'cells-1792-linear-g12.npz')))
    m=task.prior.Flow(1792);g=task.Geometry(m);d=np.load(task.prior.OUT/'cells-1792.npz');q=np.load(out/'cells-1792-linear-g12.npz');q2=np.load(out/'cells-1792-linear-g24.npz')
    t=np.r_[0.,d['history'][:,0]];states=np.concatenate([d['initial'][None],d['snapshots']]);scale=4*np.pi*m.RJ**2*m.eos.rho0
    baryon=(states[:,0].astype(np.longdouble)-states[0,0])*m.vol*scale
    w,delay,*_=g(m.x);w0=g(np.array([0.]))[0][0];at=float(q['u_seconds'][-1]);terms=[]
    def value(y,x):
        if x<=0:return np.longdouble(0)
        k=min(np.searchsorted(t,x,side='right')-1,len(t)-2)
        return y[k]+(y[k+1]-y[k])*np.longdouble((x-t[k])/(t[k+1]-t[k]))
    def integral(y,a,b):
        sign=1
        if b<a:a,b,sign=b,a,-1
        a=max(0.,a);b=max(0.,b);cuts=np.r_[a,t[(t>a)&(t<b)],b]
        v=np.array([value(y,x) for x in cuts],np.longdouble)
        return sign*np.sum(np.diff(cuts).astype(np.longdouble)*(v[1:]+v[:-1])/2,dtype=np.longdouble)
    for i in range(len(w)):
        y=baryon[:,i]
        terms.append(w[i]*integral(y,at,at+delay[i])+(w[i]-w0)*integral(y,0,at))
    rest=-task.G*m.eos.cx/(2*task.C)*np.sum(terms,dtype=np.longdouble)
    error=float(abs(rest/q['layer_rest_charge_cm'][-1]-1))
    body=abs(q['bulk_rest_charge_cm'])+abs(q['bulk_nonrest_charge_cm'])
    body_error=float(np.max(abs(q['bulk_rest_charge_cm']-q2['bulk_rest_charge_cm'])+abs(q['bulk_nonrest_charge_cm']-q2['bulk_nonrest_charge_cm']))/np.max(body))
    initial_ratio=float(abs(q['normalized_charge'][0])/np.max(abs(q['normalized_charge'])))
    # Evaluate the local acoustic field at the endpoint, rather than replacing
    # transferred baryons by an arbitrary thin mass sheet. Its integrated mass
    # is the saved opposite port and its propagation depth is c_s*a/B*t.
    raw=np.load(out/'native-bulk.npz')['raw'][0];r0,p,u,gamma=raw[[0,1,2,4]]
    saved=json.loads((out/'cells-1792-linear-g12.json').read_text());speed=saved['bulk']['coordinate_speed_cm_s'];cs=saved['bulk']['cs_over_c']
    port_values=np.r_[0,d['history'][:,3]]*scale
    tau=np.linspace(0,at,257);emitted=at-tau;k=np.clip(np.searchsorted(t,emitted,side='right')-1,0,len(t)-2)
    mass_rate=-np.diff(port_values)[k]/np.diff(t)[k];x=-20000-speed*tau
    _,a,B,re,phi=g.metric(x);volume_per_cm=4*np.pi*(m.RJ+x)**2*B
    delta_rho=mass_rate/speed/volume_per_cm
    delta_p=gamma*p/r0*delta_rho;delta_v=-cs*task.C*delta_rho/r0
    np.savez_compressed(out/'bulk-acoustic-endpoint.npz',x_cm=x,delta_rho_cgs=delta_rho,delta_pressure_cgs=delta_p,delta_velocity_cm_s=delta_v,
        physical_model='Local constant-coefficient incoming acoustic response; background area-volume conversion only. Not a new nonlinear bulk evolution.')
    result=dict(classification='Counterexample candidate',passed=error<1e-8 and body_error<1e-8,
        independent_actual_rest_integral_relative=error,body_only_quadrature_relative=body_error,
        first_readout_over_peak=initial_ratio,maximum_bulk_density_fraction=float(max(abs(delta_rho/r0))),
        maximum_bulk_velocity_cm_s=float(max(abs(delta_v))),seconds=time.monotonic()-start,
        temporal_scope='The tiny u=0 value is interpolation/cell support leakage, retained without subtraction. Sparse histories do not certify pointwise causal support.',
        full_inner_energy_closure=False,whole_physical_charge_certified=False)
    task.write(out/'history-audit.json',result);print(json.dumps(result),flush=True);assert result['passed']


if __name__=='__main__':main()
