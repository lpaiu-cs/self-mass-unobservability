"""Counterexample candidate: actual two-way stresses in linear GR constraints.

The conserved-volume tangent closure is explicit. Exterior photon metric
terms and all scalar-potential iterations are enclosed, not silently zeroed.
This is still a prescribed-source response, not a nonlinear metric re-evolution.
"""
from pathlib import Path
import json
import signal
import sys
import time
import numpy as np
import mpmath as mp
import def_native_interior_feedback as flow
import def_native_conserved_wave as wave

prior=flow.prior;green=prior.green;C=flow.C;G=prior.G;write=flow.write;sha=flow.sha
OUT=flow.OUT.parent/'def-native-feedback-gr'


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='e1862c3be',
        claim='Apply the actual two-way material and photon energy/radial-stress histories to the linear Einstein mass constraint and scalar equation. Decide whether these missing GR terms can erase the surviving direct charge.',
        closure='Conserved coordinate inventory, entropy and composition under the additional infinitesimal metric-volume response, with the saved static-background Gamma1. Matter and radiation changes already evolved in Phase116 are prescribed forcing. This is not a second evolution of photons/material in a dynamic metric.',
        equations=wave.symbolic(),
        mass_boundary='The actual inner photon Killing-energy flux supplies the negative mass change below the represented inner face. Its direct scalar emission cannot reach the observer during this horizon. Independently compare both energy ports with the saved represented energy change.',
        exterior='No incoming photons or matter beyond the outer face. Use the emitted positive energy to bound both exterior anisotropic photon stress and causal mass-debit scalar forcing, for any outward angular distribution. Do not equate trace-free radiation with zero metric forcing.',
        scalar='Compute the forced interior metric-stress and mass terms. Bound all potential iterations by the retarded operator norm using absolute source envelopes; do not use the observed cancellation as the source norm.',
        reuse='Only saved64/128 histories, fixed background and existing response operators. No new EOS bank, fluid steps, grid, frequency, angle or horizon.',
        gates=dict(port_history=.02,direct_reproduction=1e-10,time_complete=.02,history=.02,quadrature=.002,
            GR_enclosure_over_direct=.02,positive_conditional_lower=True),
        budget=dict(source_seconds=25,readout_seconds=45,bound_seconds=30,CPU_threads=1,memory_GB=3,new_fluid_steps=0),
        measured_basis='The previous actual two-path source export and direct/ray readout took9.50s. Reuse those small source arrays; only boundary photon snapshots require decompression. No large evolution is launched.',
        stop='Preserve failed gates and reassess; no automatic extra paths, weaker gates or transfer of old EOS bounds. Keep whole-star and physical charge flags false.',
        limitations='First-order conserved-volume metric closure, frozen source histories and coefficients. Initial full GR constraint matching, freely matched mechanics, cumulative neutral reaction ledger, physical radial/spectral/inner angular errors and full nonlinear backreaction remain unclosed.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(wave.__file__),Path(flow.__file__),Path(prior.__file__),
            flow.OUT/'source-64.npz',flow.OUT/'source-128.npz',flow.OUT/'coupled-64.npz',flow.OUT/'coupled-128.npz',flow.OUT/'result.json']}))


def sources():
    assert not (OUT/'sources.json').exists();start=time.monotonic();signal.signal(signal.SIGALRM,flow.old.optical.timeout);signal.alarm(25)
    model=flow.Coupled();b=model.bulk;rows=[]
    for steps in [64,128]:
        d=dict(np.load(flow.OUT/f'source-{steps}.npz'));z=np.load(flow.OUT/f'coupled-{steps}.npz');times=z['snapshot_t'];t=d['t']
        keep=np.array([np.argmin(abs(times-tt)) for tt in t]);assert np.max(abs(times[keep]-t))<1e-18
        IB=z['snapshot_bulk_I'][keep];I=z['snapshot_I'][keep].sum(1);inner=[];outer=[]
        for ib,ii in zip(IB,I):
            inward=np.where((b.mu>0)[:,None],b.incoming,ib[0]);outward=np.where((b.mu>0)[:,None],ii[-1],0.)
            inner.append(4*np.pi*C*b.area[0]*((b.w*b.mu@inward)*(b.d['num']*b.d['Einf'])).sum())
            outer.append(4*np.pi*C*model.area[-1]*((b.w*b.mu@outward)*(b.d['num']*b.d['Einf'])).sum())
        inner=np.array(inner);outer=np.array(outer);assert np.min(outer)>=0
        Hi=green.polynomial(t,inner).antiderivative()(t);Ho=green.polynomial(t,outer).antiderivative()(t)
        LD=np.longdouble;rest=np.asarray(d['baryon_g'],LD)*LD(d['cx'])*LD(C)**2
        energy=(rest+d['gas_nonrest_energy_erg']+d['photon_energy_erg'])*np.asarray(d['a'],LD)
        stress=rest+d['nonrest_stress_erg']+d['photon_energy_erg']-d['photon_radial_pressure_erg']
        photon_escape=float(z['ledger'][5]+z['scalar_deep_escape'])
        # No correction is fitted to force this independent port audit to pass.
        mismatch=np.asarray(energy.sum(1,dtype=LD),float)-(Hi-Ho)
        bound=float(np.max(abs(mismatch))+abs(photon_escape))
        scale=max(float(np.max(abs(Hi))),float(np.max(abs(Ho))),1.)
        exact_boundary=float(z['boundary'][0]);quadrature=abs(Hi[-1]-Ho[-1]-exact_boundary)/scale
        d.update(Killing_cell_energy_erg=energy,metric_stress_erg=np.asarray(stress,float),inner_luminosity=inner,outer_luminosity=outer,
            inner_cumulative_energy_erg=Hi,outer_cumulative_energy_erg=Ho,port_mismatch_erg=mismatch,port_error_bound_erg=bound,spectral_escape_erg=photon_escape)
        np.savez_compressed(OUT/f'source-{steps}.npz',**d)
        row=dict(steps=steps,port_history_relative=bound/scale,port_quadrature_vs_exact_integrator=quadrature,
            endpoint_inner_energy_erg=float(Hi[-1]),endpoint_outgoing_energy_erg=float(Ho[-1]),endpoint_represented_energy_erg=float(energy[-1].sum()),
            maximum_port_error_erg=bound,spectral_escape_erg=photon_escape)
        rows.append(row);print(json.dumps(row),flush=True)
    passed=max(r['port_history_relative'] for r in rows)<.02
    write(OUT/'sources.json',dict(classification='Counterexample candidate',passed=passed,paths=rows,seconds=time.monotonic()-start));signal.alarm(0)


def read(d,model,order=8,kind='linear'):
    edges=d['edges'];gx,gw=np.polynomial.legendre.leggauss(order);frac=(gx+1)/2
    rJ=(edges[:-1,None]+np.diff(edges)[:,None]*frac).ravel();weights=np.tile(gw/2,len(edges)-1);ids=np.repeat(np.arange(len(edges)-1),order)
    geom=green.Geometry(model.m);w,delay,a,B,re=geom(rJ-model.m.RJ);r=re*model.m.R
    z=model.m.bg.sample(re);Phi=z['v']/model.m.R;N=z['N'];V,Kj,bb,NN=wave.coeff(model.m.bg,r)
    times=np.load(flow.OUT/'wave-128.npz')['t'];t=d['t'];M=float(d['M_cm']);cx=float(d['cx'])
    direct=np.asarray(d['baryon_g'],float)*cx*C*C+d['nonrest_trace_erg'];stress=d['metric_stress_erg']
    energy=d['Killing_cell_energy_erg'];prefix=np.cumsum(energy,axis=1,dtype=np.longdouble)-energy
    Hi=green.polynomial(t,d['inner_luminosity'],kind).antiderivative()(t)
    J=np.asarray((G/C**4)*(np.sqrt(bb)/N)[None,:]*(-Hi[:,None]+prefix[:,ids]+energy[:,ids]*np.tile(frac,len(edges)-1)),float)
    HJ=green.polynomial(t,J,kind).antiderivative();HD=green.polynomial(t,direct,kind).antiderivative();HS=green.polynomial(t,stress,kind).antiderivative()
    dx=np.repeat(np.diff(edges),order)*weights*B/a;rows=[]
    def evaluate(H,at,ii):
        cut=np.clip(at,0,t[-1]);j=np.clip(np.searchsorted(H.x,cut,side='right')-1,0,len(H.x)-2);dt=cut-H.x[j];value=np.zeros_like(dt)
        for c in H.c:value=value*dt+c[j,ii]
        value[at<=0]=0;assert at.max()<=t[-1]+2e-15
        return value
    for u in times:
        at=u+delay
        direct_q=G/(2*C**3*M)*np.sum(w*weights*evaluate(HD,at,ids),dtype=np.longdouble)
        stress_q=G/(2*C**3*M)*np.sum(a*Phi*weights*evaluate(HS,at,ids),dtype=np.longdouble)
        mass_q=-C/(2*M)*np.sum(dx*Kj*evaluate(HJ,at,np.arange(len(ids))),dtype=np.longdouble)
        rows.append([direct_q,stress_q,mass_q])
    values=np.asarray(rows,float).T
    return times,values,dict(maximum_J_cm=float(np.max(abs(J))),minimum_potential=float(V.min()),maximum_potential=float(V.max()),
        exact_direct_inner_cut_delay_seconds=float(geom(np.array([edges[0]-model.m.RJ]))[1][0]),
        represented_mass_boundary='Negative actual integrated inner photon Killing-energy flux; no fitted constant or artificial point scalar source.')


def readout():
    assert not (OUT/'result.json').exists() and json.loads((OUT/'sources.json').read_text())['passed']
    start=time.monotonic();signal.signal(signal.SIGALRM,flow.old.optical.timeout);signal.alarm(45);model=flow.Coupled();saved={};rows={}
    for steps in [64,128]:
        d=dict(np.load(OUT/f'source-{steps}.npz'));t,p,row=read(d,model);saved[steps]=p;rows[steps]=row
        np.savez_compressed(OUT/f'wave-{steps}.npz',t=t,components=p,component_names=['direct_trace','metric_stress','metric_mass'],normalized_forced=p.sum(0))
    d=dict(np.load(OUT/'source-128.npz'));p=saved[128];q=p.sum(0);scale=max(abs(q));original=np.load(flow.OUT/'wave-128.npz')['direct']
    errors=dict(direct_reproduction=float(max(abs(p[0]-original))/max(abs(original))),time_complete=float(max(abs(q-saved[64].sum(0)))/scale))
    for key,kw in [('history',dict(kind='cubic')),('quadrature',dict(order=4))]:
        _,other,_=read(d,model,**kw);errors[key]=float(max(abs(q-other.sum(0)))/scale)
    gates=json.loads((OUT/'plan.json').read_text())['gates'];passed=all(v<gates[k] for k,v in errors.items())
    row=dict(classification='Counterexample candidate',passed=passed,controls=errors,paths=rows,endpoint_components=p[:,-1].tolist(),
        endpoint_forced=float(q[-1]),seconds=time.monotonic()-start,actual_photon_and_matter_metric_stresses_applied=True,
        actual_inner_photon_mass_debit_applied=True,potential_and_exterior_certified=False,full_GR_feedback=False,final_charge_solved=False,full_goal_complete=False)
    write(OUT/'result.json',row);signal.alarm(0);print(json.dumps(row),flush=True)


def bound():
    assert not (OUT/'bound.json').exists();start=time.monotonic();signal.signal(signal.SIGALRM,flow.old.optical.timeout);signal.alarm(30)
    model=flow.Coupled();bg=model.m.bg;d=dict(np.load(OUT/'source-128.npz'));z=np.load(OUT/'wave-128.npz');geom=green.Geometry(model.m)
    mp.iv.dps=40;iv=mp.iv;I=iv.mpf
    def up(x):return float(np.nextafter(float(x.b),np.inf))
    def low(x):return float(np.nextafter(float(x.a),-np.inf))
    T=I(float(d['t'][-1]));M=I(float(d['M_cm']));K=I(abs(float(d['K_cm'])));R=I(bg.R);cc=I(C);gg=I(G)
    qlo,qhi=geom(np.array([d['edges'][0],d['edges'][-1]])-model.m.RJ)[1]*C
    qmin=I(qlo)-cc*T/2;rmin=R+qmin;rmax=I(float(geom.metric(np.array([d['edges'][-1]-model.m.RJ]))[3][0]*model.m.R))
    idx=max(0,int(np.searchsorted(bg.d['radius_cm'],float(rmin.a)))-1)
    E=I(float(max(bg.d['energy_cgs'][idx:])))*gg/cc**4;P=I(float(max(bg.d['pressure_cgs'][idx:])))*gg/cc**4
    Phi=I(float(max(abs(bg.d['phi_prime_cm'][idx:]))));alpha=4*I(float(max(abs(bg.d['phi'][idx:]))));gamma=I(float(max(bg.d['gamma1'][idx:])))
    # Explicit frozen piecewise-linear coefficient envelope; add adjacent
    # node through idx and allow every coefficient extreme independently.
    Nmin=I(float(min(bg.d['lapse'][idx:])));bmin=1-2*M/rmin;H=E+P;gp=gamma*P
    Kj=2*Phi/(rmin*bmin)*(1+4*iv.pi*rmax*rmax*H)+8*iv.pi/bmin*alpha*(E+3*P)
    D=4*iv.pi*(alpha*(H+3*gp)+rmax*Phi*(H+gp))
    Kje=Kj+D/bmin
    V=2*M/rmin**3+4*iv.pi*H+Kj*rmax*Phi+16*iv.pi*alpha*rmax*Phi*H+4*iv.pi*(4+4*alpha**2)*(E+3*P)+D*(3*alpha+rmax*Phi)
    b0=1-2*M/rmax;N0=I(float(bg.metric(np.array([float(rmax.a)/bg.R]))[1][0]));C0=N0*iv.sqrt(b0)
    Ivac=(M/rmax**2+2*K*K/(3*C0*C0*rmax**3))/iv.sqrt(b0)
    integral=(I(qhi)-qmin)*V+Ivac;eta=cc*T/2*integral;assert up(eta)<1
    # Absolute source norms, independent of cancellation in the readout.
    rest=np.asarray(d['baryon_g'],float)*float(d['cx'])*C*C
    trace=rest+d['nonrest_trace_erg'];stress=d['metric_stress_erg']
    abs_trace=I(float(np.sum(np.max(abs(trace),axis=0),dtype=np.longdouble)))*I('1.000000001')
    abs_stress=I(float(np.sum(np.max(abs(stress),axis=0),dtype=np.longdouble)))*I('1.000000001')
    direct_norm=gg/(2*cc**3)*T*alpha/rmin*abs_trace
    stress_norm=gg/(2*cc**3)*T*Phi*abs_stress
    energy_abs=I(float(np.sum(np.max(abs(d['Killing_cell_energy_erg']),axis=0),dtype=np.longdouble)))*I('1.000000001')
    inner_abs=I(float(np.max(abs(d['inner_cumulative_energy_erg']))))*I('1.000000001')
    Jmax=gg/cc**4/Nmin*(energy_abs+inner_abs)
    mass_norm=cc*T/2*(I(qhi)-I(qlo))*Kje*Jmax
    emitted=I(float(d['outer_cumulative_energy_erg'][-1]))*I('1.000000001')
    assert float(d['K_cm'])<0 and up(gg*emitted/(cc**4*M))<1
    # For all outward rays, r*nu'<=kappa<1 gives
    # mu(r)^2>=(1-kappa)*(1-r0^2/r^2), even for grazing emission.
    kappa=M/(rmax*b0)+K*K/(2*C0*C0*rmax*rmax);assert up(kappa)<1
    ext_stress=gg/(2*cc**4)*emitted*K/(C0*C0)*iv.pi/(2*rmax*iv.sqrt(1-kappa))
    ext_mass=cc*T/2*gg/cc**4*emitted/N0*K/(b0*b0*rmax*rmax)
    port=I(float(d['port_error_bound_erg']))
    port_effect=cc*T/2*(I(qhi)-I(qlo))*Kje*gg/cc**4/Nmin*port
    M0=direct_norm+stress_norm+mass_norm+ext_stress+ext_mass+port_effect
    potential=eta/(1-eta)*M0
    unknown=(ext_stress+ext_mass+potential+port_effect)/M
    end=float(z['normalized_forced'][-1]);correction=np.sum(z['components'][1:],axis=0)
    total_gr=I(float(max(abs(correction))))+unknown
    lower=I(end)-unknown;upper=I(end)+unknown;direct_scale=float(max(abs(z['components'][0])))
    row=dict(classification='Counterexample candidate',conditional=True,passed=bool(up(total_gr)/direct_scale<.02 and low(lower)>0),
        scalar_potential_contraction=up(eta),absolute_free_source_norm_cm=up(M0),all_orders_potential_normalized_bound=up(potential/M),
        exterior_photon_stress_normalized_bound=up(ext_stress/M),exterior_photon_mass_normalized_bound=up(ext_mass/M),
        stored_port_history_normalized_bound=up(port_effect/M),uncomputed_GR_normalized_bound=up(unknown),
        all_GR_change_bound_over_direct=up(total_gr)/direct_scale,endpoint_forced_normalized=end,
        endpoint_direct_with_GR_interval=[low(lower),up(upper)],unknown_positive_photon_mass_does_not_erase_conditional_lower=low(lower)>0,
        arbitrary_nonnegative_outward_exterior_angle=True,inner_direct_boundary_outside_observer_cone=float(z['t'][-1])+qlo/C<0,
        coefficient_enclosure='Frozen piecewise-linear background extrema with adjacent node, m<=M, N<=1, A<=1, explicit lower r,b,N, and exact positive-mass vacuum bounds. Binary source extrema enter with a1e-9 guard; interval arithmetic is outward after those reductions.',
        exclusions='Not an enclosure of source radial/frequency/inner angular errors, native constitutive truncation, initial nonlinear Einstein constraints, free mechanical matching, unrepresented reaction channels or dynamically changed photon/matter trajectories.',
        seconds=time.monotonic()-start,full_GR_feedback=False,final_charge_solved=False,full_goal_complete=False)
    write(OUT/'bound.json',row);signal.alarm(0);print(json.dumps(row),flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
