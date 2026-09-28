"""Independent Jordan-radius readout and exact center mass constraints."""
from pathlib import Path
import json
import signal
import sys
import time
import numpy as np
from numpy.polynomial import legendre as leg
import def_native_characteristic_gr as new

base=new.base;flow=base.flow;OUT=base.OUT;write=base.write;sha=base.sha
C=base.C;G=base.G;LD=base.LD


def prepare():
    assert not (OUT/'audit-plan.json').exists()
    write(OUT/'audit-plan.json',dict(classification='Counterexample candidate',
        claim='Check the corrected characteristic readout with a separately evaluated Jordan-radius integral and reconstruct J directly at each center instead of interpolating its sparse quadrature samples.',
        original_verdicts='Both fields.json files remain failed. The first fails field quadrature; the characteristic version fails compatibility with the old unsplit direct quadrature. No prior verdict or threshold is replaced.',
        independent_decision='A separate accuracy audit passes only if an independent Jordan-radius characteristic-split integral matches the corrected optical-coordinate answer at the original1e-9 equivalence tolerance, and the original0.02 time/0.002 radial gates pass. This does not make old compatibility or full physical closure pass.',
        budget=dict(seconds=30,new_fluid_steps=0,new_native_calls=0,CPU_threads=1),
        measured_basis='Characteristic three-path evaluation18.88s; this audit is one model load and two one-dimensional endpoint integrals, then531 center algebraic constraints.',
        gates=dict(independent_equivalence=1e-9,time=.02,quadrature=.002,exact_box=3e-15),
        bindings={str(p):sha(p) for p in [Path(__file__),Path(new.__file__),Path(base.__file__),OUT/'fields.json',new.OUT/'fields.json',new.OUT/'fields-128-g8.npz',OUT/'source-128.npz']}))


def direct(model,d,order):
    """Independent integral in physical Jordan radius, not interpolated x density."""
    geo=model.geo;edges=d['edges'];T=d['t'][-1];gx,gw=leg.leggauss(order)
    xf=C*geo(edges-model.model.m.RJ)[1];xo=xf[-1]
    cuts=xo-C*(T-d['t']);cuts=cuts[(cuts>xf[0])&(cuts<xf[-1])]
    rr=np.interp(cuts,xf,edges)
    for _ in range(4):
        _,delay,a,B,_=geo(rr-model.model.m.RJ);rr-=(C*delay-cuts)*a/B
    rerr=float(np.max(abs(C*geo(rr-model.model.m.RJ)[1]-cuts),initial=0))
    split=np.unique(np.r_[edges,rr]);r=(split[:-1,None]+split[1:,None])/2+np.diff(split)[:,None]*gx/2
    r=r.ravel();width=np.repeat(np.diff(split)/2,order)*np.tile(gw,len(split)-1)
    ids=np.clip(np.searchsorted(edges,r,side='right')-1,0,len(edges)-2)
    q=flow.initial.Quadrature(edges,order);_,_,_,BB,_=geo(q.r.ravel()-model.model.m.RJ)
    norm=q.h*((q.r*q.r*BB.reshape(q.r.shape))@q.w)
    w,delay,a,B,re=geo(r-model.model.m.RJ);ret=T-(xo/C-delay)
    rest=np.asarray(d['baryon_g'],LD)*LD(d['cx'])*LD(C)**2
    trace=np.asarray(rest+d['nonrest_trace_erg'],float);H=flow.green.polynomial(d['t'],trace).antiderivative()
    cut=np.clip(ret,0,T);ii=np.clip(np.searchsorted(d['t'],cut,side='right')-1,0,len(d['t'])-2);dt=cut-d['t'][ii]
    value=np.zeros_like(r)
    for co in H.c:value=value*dt+co[ii,ids]
    value[ret<=0]=0.
    result=G/(2*C**3*float(d['M_cm']))*np.sum(width*r*r*B/norm[ids]*w*value,dtype=LD)
    return float(result),rerr


def centers(model,d,field,order=8):
    model.setup(d,order);z=model.tz;r=model.tr;target=model.targets;edges=d['edges'];n=len(edges)-1
    for key in ['Eg','Pg','Er','Pr','Kg','R4']:z[key][-1]=0.
    z['nu_prime'][-1]=z['mass'][-1]/(r[-1]**2*z['b'][-1])+r[-1]*z['Phi'][-1]**2/2
    q=flow.initial.Quadrature(edges,order);gx,gw=q.x,q.w
    w,delay,a,B,re=model.geo(q.r.ravel()-model.model.m.RJ)
    measure=(q.r*q.r*B.reshape(q.r.shape));norm=q.h*(measure@gw)
    mean_a=q.h*((measure*a.reshape(q.r.shape))@gw)/norm
    part=(edges[:-1,None]+target[:-1,None])/2+(target[:-1]-edges[:-1])[:,None]*gx/2
    _,_,pa,pB,_=model.geo(part.ravel()-model.model.m.RJ)
    fraction=(target[:-1]-edges[:-1])/2*((part*part*(pa*pB).reshape(part.shape))@gw)/norm
    rest=np.asarray(d['baryon_g'],LD)*LD(d['cx'])*LD(C)**2
    energy=rest+d['gas_nonrest_energy_erg']+d['photon_energy_erg'];killing=energy*mean_a
    prefix=np.cumsum(killing,axis=1,dtype=LD)-killing
    bracket=np.column_stack([prefix+energy*fraction,killing.sum(1,dtype=LD)])-d['inner_cumulative_energy_erg'][:,None]
    J=np.asarray(G/C**4*np.sqrt(z['b'])/z['lapse']*bracket,float)
    f=field['delta_phi'];fr=field['delta_Phi'];dm=r*r*z['b']*z['Phi']*f+J;dl=dm/(r*z['b']);dv=3*z['alpha']*f+dl
    dp=-z['Kg']*dv-4*z['alpha']*z['Pr']*f-(3*z['Pr']-z['R4'])*dl
    prescribed=np.asarray((energy-d['metric_stress_erg'])/d['volume']*LD(G)/LD(C)**4,float)
    # The outer value is the vacuum-side pressure, not a copied last cell.
    pF=np.column_stack([prescribed,np.zeros(len(d['t']))]);dp+=pF
    P=z['Pg']+z['Pr'];N=z['lapse'];b=z['b'];A=np.exp(-2*z['phi']**2)
    dnu=(1+8*np.pi*r*r*z['A4']*P)*dm/(r*r*b*b)+4*np.pi*r*z['A4']*(dp+4*z['alpha']*P*f)/b+r*z['Phi']*fr
    # Variation of the proper static-frame gravitational acceleration.
    dg=C*C*np.sqrt(b)/A*(dnu-4*z['Phi']*f+z['alpha']*fr-(dl+z['alpha']*f)*(z['nu_prime']+z['alpha']*z['Phi']))
    return dict(t=d['t'],radius_E=r,J=J,delta_mass_cm=dm,delta_lambda=dl,delta_log_proper_volume=dv,
        delta_nu_prime=dnu,delta_proper_gravitational_acceleration=dg,delta_gas_energy_geom=-(z['Eg']+z['Pg'])*dv,
        delta_gas_radial_pressure_geom=-z['Kg']*dv,delta_photon_energy_geom=-4*z['alpha']*z['Er']*f-(z['Er']+z['Pr'])*dl,
        delta_photon_radial_pressure_geom=-4*z['alpha']*z['Pr']*f-(3*z['Pr']-z['R4'])*dl)


def audit():
    assert not (OUT/'audit.json').exists()
    plan=json.loads((OUT/'audit-plan.json').read_text())
    for p,h in plan['bindings'].items():assert sha(p)==h,p
    signal.signal(signal.SIGALRM,flow.old.optical.timeout);signal.alarm(30);start=time.monotonic()
    model=new.Response();d=dict(np.load(OUT/'source-128.npz'));field=np.load(new.OUT/'fields-128-g8.npz')
    refs=[direct(model,d,order) for order in [4,8]];nominal=json.loads((new.OUT/'fields-128-g8.json').read_text())['endpoint_direct']
    reference=max(abs(q/nominal-1) for q,_ in refs)
    cs=centers(model,d,field);np.savez_compressed(OUT/'center-constraints.npz',**cs)
    old=json.loads((OUT/'fields.json').read_text());fixed=json.loads((new.OUT/'fields.json').read_text());box=json.loads((new.OUT/'check.json').read_text())
    passed=reference<1e-9 and fixed['controls']['time']<.02 and fixed['controls']['quadrature']<.002 and box['passed']
    result=dict(classification='Counterexample candidate',passed=bool(passed),
        independent_Jordan_direct_endpoints=[q for q,_ in refs],optical_direct_endpoint=nominal,
        independent_equivalence_relative=reference,coordinate_inverse_residual_cm=max(x for _,x in refs),
        legacy_compatibility_passed=fixed['passed'],original_unresolved_quadrature_passed=old['passed'],
        old_unsplit_direct_compatibility_error=fixed['controls']['direct_reproduction'],
        field_time_relative=fixed['controls']['time'],field_quadrature_relative=fixed['controls']['quadrature'],
        exact_center_J_vs_old_interpolation_relative=float(np.max(abs(cs['J']-field['J_center_interpolated']))/np.max(abs(cs['J']))),
        maximum_delta_mass_cm=float(np.max(abs(cs['delta_mass_cm']))),maximum_delta_lambda=float(np.max(abs(cs['delta_lambda']))),
        maximum_delta_nu_prime_per_cm=float(np.max(abs(cs['delta_nu_prime']))),
        maximum_delta_proper_acceleration_cm_s2=float(np.max(abs(cs['delta_proper_gravitational_acceleration']))),
        endpoint_compact_with_metric=fixed['paths'][0]['endpoint_compact_with_metric'],
        exterior_or_deep_source_enclosed=False,spatial_material_photon_feedback=False,nonlinear_GR=False,final_charge_solved=False,
        interpretation='Independent accuracy audit only. Original failed raw gates are preserved; no finite-resolution comparison is a uniform error certificate. Center forces include represented first variations and exact cell-integral J. Lapse boundary, full spatial photon/material feedback and exterior/deep causal forcing remain.',
        seconds=time.monotonic()-start)
    write(OUT/'audit.json',result);print(json.dumps(result),flush=True);signal.alarm(0)


if __name__=='__main__':globals()[sys.argv[1]]()
