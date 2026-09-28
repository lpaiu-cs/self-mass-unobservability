"""Counterexample candidate: saved two-way sources and independent closure checks."""
from pathlib import Path
import json
import signal
import sys
import time
import numpy as np
import sympy as sp
import def_native_interior_feedback as run

OUT=run.OUT;prior=run.prior;C=run.C;G=prior.G;write=run.write;sha=run.sha


def prepare():
    assert not (OUT/'readout-plan.json').exists()
    write(OUT/'readout-plan.json',dict(classification='Counterexample candidate',
        claim='Use actual two-way64/128 histories to decide whether the earlier direct-charge cancellation survives. Use actual redistributed interior baryons; remove the old conditional point debit.',
        method='Same radial retarded integral and continuous vacuum rays as Phase114. Retain17 common saved times, stable conserved rest/nonrest sources, and all actual emitted photons. Compare direct charge and direct-plus-photon mass separately.',
        gates=dict(time_direct=.02,time_total=.02,time_deep_trace=.02,history_direct=.02,cadence_direct=.02,quadrature_direct=.02,ray_quadrature=.002,native_anchor=.002),
        budget=dict(saved_readout_seconds=60,independent_audit_seconds=30,native_call_cap=400,new_fluid_steps=0,CPU_threads=1,memory_GB=3),
        stop='Keep any failed gate; no automatic time/radial/frequency refinement. A cancelled direct residual does not justify normalizing its error by the much larger photon term.',
        limits='Fixed metric,16 deep cells, first-order density/inventory closure, finite microphysics and imposed mechanical interface. Frozen initial outgoing photons are only a mass-term reference, not a fully evolved static comparator. This is not final physical charge.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(run.__file__),Path(run.old.__file__),Path(prior.__file__)]}))


def source(model,steps):
    z=np.load(OUT/f'coupled-{steps}.npz');assert int(z['completed_steps'])==steps
    keys=['U','I','bulk_I','u','theta','eta','h','j','t']
    data={k:z['snapshot_'+k] for k in keys};t=np.arange(17)*run.old.END/16
    keep=np.array([np.argmin(abs(data['t']-tt)) for tt in t]);assert np.max(abs(data['t'][keep]-t))<1e-18
    data={k:v[keep] for k,v in data.items()};f=model.flow;m=model.m;b=model.bulk
    volumes=np.r_[b.volume,4*np.pi*m.RJ**2*m.vol];edges=np.r_[b.d['edges'][:-1],m.rf];radius=np.r_[b.d['r'],m.r]
    weight,delay,a,B,re=prior.green.Geometry(m)(radius-m.RJ)
    f.eos.y=f.eos.y0;ip=f.eos(f.initial[0],f.initial_temperature)[0];tau0=(f.initial[2]-(m.a-m.a0)*f.eos.cx*f.initial[0])/m.a
    bp0=b.eos.base.gas(np.zeros(b.n),np.zeros(b.n))[0]
    trace=[];stress=[];baryon=[];gasE=[];press=[];phE=[];phP=[];velocity=[];lum=[];native_u=[]
    for k,state in enumerate(data['U']):
        model.h=data['h'][k];model.j=data['j'][k];dm=-np.diff(model.h);model.mass=model.mass0-model.h[1:]+model.h[:-1];model.set_material(t[k])
        rho,v,lt,y=f.primitive(state);p,*_=f.eos(rho,lt);tau=(state[2]-(m.a-m.a0)*f.eos.cx*state[0])/m.a
        pp,uu,*_=b.eos.gas(data['theta'][k],data['eta'][k]);u=data['u'][k];native_u.append(float(np.max(abs(uu-u)/np.maximum(abs(u),1.))))
        beta=model.velocity();kinetic_trace=-model.mass*(model.cx*C*C+u)*beta**2/(1+np.sqrt(1-beta**2))
        de=dm*b.u0+model.mass*(u-b.u0);deep_nonrest=de+kinetic_trace;dpb=(pp-bp0)*b.volume
        unit=f.eos.rho0*C*C*volumes[b.n:];dp=p-ip
        trace.append(np.r_[deep_nonrest-3*dpb,((tau-tau0)-state[1]*v-3*dp)*unit])
        stress.append(np.r_[deep_nonrest-dpb,((tau-tau0)-state[1]*v-dp)*unit])
        gasE.append(np.r_[de+model.kinetic(),(tau-tau0)*unit]);press.append(np.r_[dpb,dp*unit])
        baryon.append(np.r_[dm,(state[0].astype(np.longdouble)-f.initial[0])*f.eos.rho0*volumes[b.n:]])
        velocity.append(np.r_[beta,v]);II=data['I'][k].sum(0);IB=data['bulk_I'][k];en=b.d['num']*b.d['Einf']
        phE.append(np.r_[np.einsum('iqf,q,f->i',IB,b.w,en)/b.d['a']**4,np.einsum('iqf,q,f->i',II,b.w,en)/m.a**4]*volumes)
        phP.append(np.r_[np.einsum('iqf,q,f->i',IB,b.w*b.mu2,en)/b.d['a']**4,np.einsum('iqf,q,f->i',II,b.w*b.mu2,en)/m.a**4]*volumes)
        lum.append(2*np.pi*C*model.area[-1]*(II[-1,-1]@(model.number*model.E)))
    f.seed=f.initial_temperature.copy();rho,v,lt,y=f.primitive(data['U'][-1]);p,*_=f.eos(rho,lt)
    check=((tau-tau0)-data['U'][-1,1]*v-3*(p-ip))*unit
    seed=float(abs((check-np.array(trace)[-1,b.n:])@m.a)/max(np.max(abs(np.array(trace)[:,b.n:]@m.a)),1.))
    assert seed<1e-8 and max(native_u)<1e-11
    d=dict(t=t,radius=radius,edges=edges,volume=volumes,weight=weight,delay=delay,a=a,B=B,re=re,
        nonrest_trace_erg=np.array(trace),nonrest_stress_erg=np.array(stress),baryon_g=np.array(baryon),
        gas_nonrest_energy_erg=np.array(gasE),pressure_volume_erg=np.array(press),photon_energy_erg=np.array(phE)-phE[0],
        photon_radial_pressure_erg=np.array(phP)-phP[0],velocity=np.array(velocity),luminosity_per_mu=np.array(lum),
        deep_cells=b.n,cx=f.eos.cx,M_cm=m.bg.M*m.R,K_cm=m.bg.K*m.R,inner_material_face_cm=m.rf[0],RJ=m.RJ)
    np.savez_compressed(OUT/f'source-{steps}.npz',**d)
    return d,dict(steps=steps,primitive_seed_relative=seed,saved_u_EOS_relative=max(native_u),maximum_primitive_residual=f.max_recovery)


def direct(d,model,**kw):
    t,q,parts,_=prior.integrate(d,model,**kw)
    # Reuse the already verified energy-source integral for the equivalent
    # rest energy in the deep cells. Never use its conditional point debit.
    rest=np.zeros_like(d['nonrest_trace_erg']);n=int(d['deep_cells'])
    rest[:,:n]=np.asarray(d['baryon_g'][:,:n],float)*float(d['cx'])*C*C
    r=dict(d,nonrest_trace_erg=rest,baryon_g=np.zeros_like(rest))
    _,_,rp,_=prior.integrate(r,model,**kw);pieces=np.vstack([parts[:3],rp[0]])
    return t,pieces.sum(0),pieces


def readout():
    assert not (OUT/'result.json').exists() and json.loads((OUT/'production.json').read_text())['passed']
    plan=json.loads((OUT/'readout-plan.json').read_text())
    for p,h in plan['bindings'].items():assert sha(p)==h,p
    start=time.monotonic();signal.signal(signal.SIGALRM,run.old.optical.timeout);signal.alarm(60)
    model=run.Coupled();datasets={};waves={};rows=[]
    for steps in [64,128]:
        d,row=source(model,steps);datasets[steps]=d;rows.append(row);t,q,parts=direct(d,model)
        E,ray=prior.rays(model,d['t'],d['luminosity_per_mu'],t);E-=E[0]
        frozen,_=prior.rays(model,d['t'],np.full(17,d['luminosity_per_mu'][0]),t);frozen-=frozen[0]
        M=float(d['M_cm']);alpha=-float(d['K_cm'])/M;eps=G*E/(C**4*M);eps0=G*frozen/(C**4*M)
        total=(q-q[0]+alpha*eps)/(1-eps);excess=total-alpha*eps0/(1-eps0)
        waves[steps]=dict(t=t,direct=q,direct_relative=q-q[0],components=parts,direct_plus_photon_mass=total,excess_over_frozen_photon_mass=excess,arrived_energy_erg=E,frozen_initial_energy_erg=frozen)
        np.savez_compressed(OUT/f'wave-{steps}.npz',**waves[steps],component_names=['deep_nonrest','atmosphere_nonrest','atmosphere_rest','deep_rest'])
    d=datasets[128];fine=waves[128];coarse=waves[64];q=fine['direct'];scale=max(abs(q));total=fine['direct_plus_photon_mass'];n=model.bulk.n
    def relative(a,b):return float(np.max(abs(a-b))/max(float(np.max(abs(a))),1e-300))
    errors=dict(time_direct=relative(q,coarse['direct']),time_total=relative(total,coarse['direct_plus_photon_mass']),
        time_deep_trace=relative(d['nonrest_trace_erg'][:,:n]@d['a'][:n],datasets[64]['nonrest_trace_erg'][:,:n]@d['a'][:n]))
    for key,kw in [('history_direct',dict(kind='cubic')),('cadence_direct',dict(stride=2)),('quadrature_direct',dict(order=4))]:
        _,other,_=direct(d,model,**kw);errors[key]=relative(q,other)
    other,_=prior.rays(model,d['t'],d['luminosity_per_mu'],t,vacuum_order=16,angular_order=4);other-=other[0]
    errors['ray_quadrature']=relative(fine['arrived_energy_erg'],other)
    upper=(1-model.bulk.edges_mu[-2]**2)/2*prior.green.polynomial(d['t'],d['luminosity_per_mu']).antiderivative()(np.minimum(t+ray['radial_face_delay_seconds'],d['t'][-1]))
    angular_bound=alpha*G*upper/(C**4*M)
    passed=all(value<plan['gates'][key] for key,value in errors.items())
    row=dict(classification='Counterexample candidate',passed=passed,controls=errors,sources=rows,
        endpoint_direct=float(q[-1]),endpoint_direct_relative=float(fine['direct_relative'][-1]),endpoint_components=fine['components'][:,-1].tolist(),
        endpoint_direct_plus_photon_mass=float(total[-1]),endpoint_excess_over_frozen_photon_mass=float(fine['excess_over_frozen_photon_mass'][-1]),
        endpoint_arrived_energy_erg=float(fine['arrived_energy_erg'][-1]),positive_subbin_mass_normalization_upper=float(angular_bound[-1]),
        seconds=time.monotonic()-start,new_fluid_steps=0,actual_two_way_material_photon_feedback=True,conditional_point_debit=False,
        old_kernel_trajectories_recertified=False,first_order_interior_constitutive_closure=True,full_angular_error_certified=False,
        spatial_frequency_continuum_certified=False,full_GR_feedback=False,final_charge_solved=False,full_goal_complete=False)
    write(OUT/'result.json',row);signal.alarm(0);print(json.dumps(row),flush=True)


def audit():
    assert not (OUT/'audit.json').exists();start=time.monotonic();signal.signal(signal.SIGALRM,run.old.optical.timeout);signal.alarm(30)
    m=run.Coupled();b=m.bulk;z=np.load(OUT/'coupled-128.npz');m.h=z['h'];m.j=z['j'];m.mass=m.mass0-m.h[1:]+m.h[:-1];m.set_material(run.old.END)
    # Independent native EOS anchors at the actual final density,T,H state.
    # This excludes the declared advected-inventory term; it cannot certify
    # that first-order inventory model or its unsampled derivatives globally.
    native=run.chem.old.Native(cap=400);theta=z['theta'];eta=z['eta'];xi=b.eos.xi.copy();b.eos.xi[:]=0
    pp,uu,*_=b.eos.gas(theta,eta);err=[]
    for j in range(b.n):
        run.chem.setup(native,b.d,j);state=native.state(np.log1p(b.eos.x[j]),np.log(b.d['T'][j])+theta[j],b.d['y0'][j]*(1+eta[j]))
        err.append([abs(pp[j]/state['p']-1),abs(uu[j]-state['u'])/max(abs(state['u']),1.)])
    b.eos.xi=xi
    # The same bound-free number collision is the opposite neutral-H source.
    I=z['bulk_I'];beta=m.velocity();ab,em,*_=b.eos.radiation(theta,eta)
    _,_,_,extra=m.collision(I,theta,eta,beta)
    bound=b.factor[:,None,None]*(em[:,None,:]-(ab-em)[:,None,:]*I)+extra
    photons=np.einsum('iqf,q,f->i',bound,b.w,b.d['num'])/b.d['a']**3
    # Independently sum per-ray absorption/emission photon packets.
    neutral=np.sum(bound*(b.w[None,:,None]*b.d['num'][None,None,:]),axis=(1,2))/b.d['a']**3
    species=float(np.max(abs(photons-neutral))/max(float(np.max(abs(photons))),1.))
    D,u,p,v,cx,c=sp.symbols('D u p v cx c',positive=True)
    E=D*(cx*c*c+u+p/(D*sp.sqrt(1-v*v)))/sp.sqrt(1-v*v)-p
    S=(E+p)*v
    assert sp.simplify(E-S*v-3*p-(D*sp.sqrt(1-v*v)*(cx*c*c+u)-3*p))==0
    # Exact positive two-node packet interpolation preserves all first moments.
    x,l,r=sp.symbols('x l r');wl=(r-x)/(r-l);wr=(x-l)/(r-l)
    assert sp.simplify(wl+wr-1)==0 and sp.simplify(wl*l+wr*r-x)==0
    native_error=float(np.max(err));passed=native_error<.002 and species<1e-9
    row=dict(classification='Counterexample candidate',passed=bool(passed),native_endpoint_pressure_energy_relative=native_error,
        native_calls=native.ion.calls,instantaneous_boundfree_number_pair_relative=species,trace_identity_symbolic=True,packet_interpolation_identity_symbolic=True,
        full_trajectory_neutral_species_ledger=False,limitation='Endpoint fixed-inventory anchors and an instantaneous paired-number check do not replace advected composition evolution, an accumulated neutral reaction ledger or uniform EOS derivative bounds.',seconds=time.monotonic()-start)
    write(OUT/'audit.json',row);signal.alarm(0);print(json.dumps(row),flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
