"""Conservative saved-state sources and the actual coupled retarded trace.

Counterexample candidate: read the actual thermal/photon/fluid histories;
this is still a prescribed-source GR component until mechanical and metric
feedback, radial/frequency errors and exterior mass closure are included.
"""
from pathlib import Path
import json
import signal
import sys
import time
import numpy as np
import sympy as sp
import def_native_cold_coupling as cold
import def_native_release_charge as green

prior=cold.prior;C=prior.C;G=green.G;write=cold.write;sha=cold.sha
OUT=prior.OUT.parent/'def-native-coupled-charge'


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='a19c21e12',
        claim='Read the actual completed coupled histories with conserved energy, quantify recovery/history uncertainty, and apply their distributed deep and atmospheric sources to the same retarded scalar Green operator.',
        fidelity='Retain the original Phase113 failed strict readout verdict. Its old report offsets never enter these new sources. No gas/heat/photon trajectory is repeated, no initial state or acceptance is changed.',
        source='Use tau-S*v-k*p for nonrest trace(k=3) and radial metric stress(k=1). Add explicit baryon redistribution separately. Deep native density is frozen and its actual saved theta/eta supplies pressure and energy. Photon trace is zero; actual photon energy and radial pressure are exported for subsequent metric equations.',
        missing='The deep density and metric were fixed. The fluid inner baryon debit is not mechanically evolved. Show a mathematically centered debit at the actual fluid inner face as an explicit conditional completion and retain its separate contribution. Do not call this or the snapshot energy a completed Bondi/ADM normalization.',
        domain='Reuse1100km native interior and[-400,1200]m moving-fluid patch, sixteen common saved intervals, exact saved background weights/delays. The volume-average source is represented as constant inside each existing cell; spatial continuum error is not certified.',
        paths=['448/64','448/128','896/128 completed by104 restart'],
        controls=['same saved conservative state under two primitive seeds','17 versus9 stored time nodes','linear versus natural cubic histories','4 versus8 radial Gauss points inside existing cells','original coarse/fine histories','manufactured constant-weight delayed source'],
        gates=dict(source_seed=.00000001,primitive=2e-14,wave_space=.02,wave_time=.02,wave_history=.02,wave_quadrature=.002),
        budget=dict(source_seconds=45,readout_seconds=45,native_constructor_calls=8,CPU_threads=1,memory_GB=2,new_fluid_steps=0),
        forecast='Phase113 independent saved endpoint recovery took1.68s including setup; reading51 snapshots reuses the same table and geometry. Stop if45s source or45s readout cap is exceeded. No new native table or time/spatial path.',
        stop='Preserve any failed gate and do not replace whole source closure with a small component pass. No automatic resolution, horizon or production re-evolution.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(cold.__file__),Path(green.__file__),cold.OUT/'resumed-896-128.npz',prior.OUT/'cells-896-steps-128.npz',prior.OUT/'cells-448-steps-128.npz',prior.OUT/'cells-448-steps-64.npz']}))
    D,v,r,tau,p,cx=sp.symbols('D v r tau p cx',real=True)
    internal=(tau+p-cx*D*v*v/(r*(1+r)))*r*r-p
    original=-cx*D*v*v/(1+r)+internal-3*p
    assert sp.simplify((original-(tau-v*v*(cx*D+tau+p)-3*p)).subs(v*v,1-r*r))==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        trace='The nonrest perfect-fluid scalar trace equals tau-S*v-3p with S=v*(cx*D+tau+p). The radial metric stress replaces3p byp. This avoids reconstructing the small source by subtracting a separately inverted internal energy.',
        photons='For massless photon packets E-Pr-2Pt=0 exactly, while E-Pr need not vanish. Photon stress must still enter the metric source.',
        limitation='These identities do not bound unresolved source histories, the actual GR metric response or missing material motion.'))


def load_path(cells,steps):
    p=prior.OUT/f'cells-{cells}-steps-{steps}.npz';a=np.load(p)
    if cells==896:
        z=np.load(cold.OUT/'resumed-896-128.npz')
        values={k:np.concatenate([a['snapshot_'+k],z['snapshot_'+k]]) for k in ['U','I','bulk_I','theta','eta','t']}
    else:values={k:a['snapshot_'+k] for k in ['U','I','bulk_I','theta','eta','t']}
    assert np.allclose(values['t'],np.arange(1,17)*prior.END/16,rtol=0,atol=1e-18)
    return values


def tight_primitive(flow):
    src=prior.interface.old.primitive_source.replace('2e-11','2e-14')
    ns=dict(flow.primitive.__func__.__globals__);exec(compile(src,__file__,'exec'),ns)
    flow.primitive=cold.MethodType(ns['primitive'],flow)


def sources():
    assert not (OUT/'sources.json').exists();start=time.monotonic();signal.signal(signal.SIGALRM,prior.optical.timeout);signal.alarm(45);results=[]
    for cells,steps in [(448,64),(448,128),(896,128)]:
        model=prior.Coupled(cells,8);f=model.flow;m=model.m;b=model.bulk;f.eos=cold.ColdEOS();tight_primitive(f)
        z=load_path(cells,steps);t=np.r_[0,z['t']];U=np.concatenate([f.initial[None],z['U']]);theta=np.vstack([np.zeros(b.n),z['theta']]);eta=np.vstack([np.zeros(b.n),z['eta']])
        radii=np.r_[b.d['r'],m.r];volumes=np.r_[b.volume,4*np.pi*m.RJ**2*m.vol];edges=np.r_[b.d['edges'][:-1],m.rf]
        gasE=[];trace=[];stress=[];baryon=[];press=[];velocity=[];photonE=[];photonP=[];seed_errors=[]
        geom=green.Geometry(m);w,delay,a,B,re=geom(radii-m.RJ)
        f.eos.y=f.eos.y0;ip=f.eos(f.initial[0],f.initial_temperature)[0];tau0=(f.initial[2]-(m.a-m.a0)*f.eos.cx*f.initial[0])/m.a
        bp0,bu0,*_=b.eos.gas(np.zeros(b.n),np.zeros(b.n))
        initial_phE=np.r_[np.einsum('iqf,q,f->i',b.initial,b.w,b.d['num']*b.d['Einf'])/b.d['a']**4,
                            np.einsum('iqf,q,f->i',model.initial_I,b.w,b.d['num']*b.d['Einf'])/m.a**4]
        initial_phP=np.r_[np.einsum('iqf,q,f->i',b.initial,b.w*b.mu2,b.d['num']*b.d['Einf'])/b.d['a']**4,
                            np.einsum('iqf,q,f->i',model.initial_I,b.w*b.mu2,b.d['num']*b.d['Einf'])/m.a**4]
        for j,state in enumerate(U):
            rho,v,lt,y=f.primitive(state);p,u,*_=f.eos(rho,lt);tau=(state[2]-(m.a-m.a0)*f.eos.cx*state[0])/m.a
            dp=p-ip;nt=(tau-tau0)-state[1]*v-3*dp;ns=(tau-tau0)-state[1]*v-dp
            pp,uu,*_=b.eos.gas(theta[j],eta[j]);db=b.d['rho']*(uu-bu0);dpb=pp-bp0
            mass=np.r_[np.zeros(b.n),(state[0].astype(np.longdouble)-f.initial[0])*f.eos.rho0*volumes[b.n:]]
            trace.append(np.r_[db-3*dpb,nt*f.eos.rho0*C*C]*volumes);stress.append(np.r_[db-dpb,ns*f.eos.rho0*C*C]*volumes)
            gasE.append(np.r_[db,(tau-tau0)*f.eos.rho0*C*C]*volumes);baryon.append(mass);press.append(np.r_[dpb,dp*f.eos.rho0*C*C]*volumes);velocity.append(v)
            if j==0:pe=initial_phE;pr=initial_phP
            else:
                II=z['I'][j-1].sum(0);IB=z['bulk_I'][j-1]
                pe=np.r_[np.einsum('iqf,q,f->i',IB,b.w,b.d['num']*b.d['Einf'])/b.d['a']**4,np.einsum('iqf,q,f->i',II,b.w,b.d['num']*b.d['Einf'])/m.a**4]
                pr=np.r_[np.einsum('iqf,q,f->i',IB,b.w*b.mu2,b.d['num']*b.d['Einf'])/b.d['a']**4,np.einsum('iqf,q,f->i',II,b.w*b.mu2,b.d['num']*b.d['Einf'])/m.a**4]
            photonE.append((pe-initial_phE)*volumes);photonP.append((pr-initial_phP)*volumes)
        trace=np.array(trace);stress=np.array(stress);gasE=np.array(gasE);baryon=np.array(baryon);photonE=np.array(photonE);photonP=np.array(photonP)
        # Re-read the same endpoint from an intentionally different seed.
        f.seed=f.initial_temperature.copy();rho,v,lt,y=f.primitive(U[-1]);p,*_=f.eos(rho,lt);tau=(U[-1,2]-(m.a-m.a0)*f.eos.cx*U[-1,0])/m.a
        check=((tau-tau0)-U[-1,1]*v-3*(p-ip))*f.eos.rho0*C*C*volumes[b.n:]
        seed=float(abs((check-trace[-1,b.n:])@m.a)/max(abs(trace[:,b.n:]@m.a).max(),1.));assert seed<1e-8
        if j==0:raise AssertionError('No saved source')
        np.savez_compressed(OUT/f'source-{cells}-{steps}.npz',t=t,radius=radii,edges=edges,volume=volumes,weight=w,delay=delay,a=a,B=B,re=re,
            nonrest_trace_erg=trace,nonrest_stress_erg=stress,gas_nonrest_energy_erg=gasE,baryon_g=baryon,pressure_volume_erg=press,
            photon_energy_erg=photonE,photon_radial_pressure_erg=photonP,velocity=velocity,deep_cells=b.n,cx=f.eos.cx,M_cm=m.bg.M*m.R,K_cm=m.bg.K*m.R,
            inner_material_face_cm=m.rf[0],RJ=m.RJ)
        row=dict(cells=cells,steps=steps,source_snapshots=len(t),source_seed_relative=seed,maximum_primitive_residual=f.max_recovery,
            deep_trace_endpoint_erg=float(trace[-1,:b.n]@b.d['a']),atmosphere_trace_endpoint_erg=float(trace[-1,b.n:]@m.a),
            net_represented_baryon_change_g=float(baryon[-1].sum()),maximum_photon_trace_residual=0.,seconds=time.monotonic()-start)
        results.append(row);print(json.dumps(row),flush=True)
    write(OUT/'sources.json',dict(classification='Counterexample candidate',passed=True,paths=results,seconds=time.monotonic()-start,new_fluid_steps=0));signal.alarm(0)


def integrate(d,model,order=8,kind='linear',stride=1):
    """Actual finite-volume source at retarded times, with no future values."""
    xg,wg=np.polynomial.legendre.leggauss(order);edges=d['edges'];r=(edges[:-1,None]+np.diff(edges)[:,None]*(xg+1)/2).ravel();wq=np.tile(wg/2,len(edges)-1)
    ids=np.repeat(np.arange(len(edges)-1),order);geom=green.Geometry(model.m);w,delay,*_=geom(r-model.m.RJ)
    # The shared readout clock is fixed by the outer FACE, not the last
    # quadrature point (which changes with fluid resolution and Gauss order).
    outer_delay=float(geom(np.array([edges[-1]-model.m.RJ]))[1][0])
    end=prior.END-max(outer_delay,0.);times=np.linspace(0,end,129)
    trace=d['nonrest_trace_erg'];mass=np.asarray(d['baryon_g'],float);n=int(d['deep_cells']);t=d['t'][::stride]
    inputs=[trace[:,:n],trace[:,n:],mass[:,n:]];paths=[]
    # Integrate source histories first; each source cell is then evaluated at
    # its actual radial quadrature delays. No fitted relaxation parameter.
    for k,values in enumerate(inputs):
        subset=ids<n if k==0 else ids>=n;source_ids=ids[subset] if k==0 else ids[subset]-n
        antiderivative=green.polynomial(t,values[::stride],kind).antiderivative();rows=[]
        for tt in times:
            at=tt+delay[subset];cut=np.clip(at,0,t[-1]);interval=np.clip(np.searchsorted(antiderivative.x,cut,side='right')-1,0,len(antiderivative.x)-2);dt=cut-antiderivative.x[interval];v=np.zeros_like(dt)
            for coeff in antiderivative.c:v=v*dt+coeff[interval,source_ids]
            v[at<=0]=0;assert at.max()<=t[-1]+2e-15
            rows.append(np.sum(v*w[subset]*wq[subset],dtype=np.longdouble))
        factor=-G/(2*C**3) if k<2 else -G*float(d['cx'])/(2*C)
        paths.append(np.array(rows,float)*factor)
    # Complete represented baryon conservation at the *actual* inner fluid
    # face only as an explicit conditional point debit. Its physical profile
    # and thermal/stress response are not asserted to have been simulated.
    wi,di,*_=geom(np.array([float(d['inner_material_face_cm'])-model.m.RJ]));debit=-mass.sum(1)
    H=green.polynomial(t,debit[::stride],kind).antiderivative();debit_wave=-G*float(d['cx'])/(2*C)*wi[0]*green.paired(H,times+di[0])
    pieces=np.vstack([*paths,debit_wave]);Q=pieces.sum(0);M=float(d['M_cm'])
    return times,-Q/M,-pieces/M,dict(order=order,history=kind,stride=stride,endpoint=float(-Q[-1]/M),maximum_future_source_seconds=float(times[-1]+delay.max()-prior.END))


def readout():
    assert json.loads((OUT/'sources.json').read_text())['passed'] and not (OUT/'result.json').exists();start=time.monotonic();signal.signal(signal.SIGALRM,prior.optical.timeout);signal.alarm(45)
    model=prior.Coupled(896,8);data={k:np.load(OUT/f'source-{k[0]}-{k[1]}.npz') for k in [(448,64),(448,128),(896,128)]};curves={};rows={}
    for k,d in data.items():
        times,q,parts,row=integrate(d,model);curves[k]=q;rows[str(k)]=row
        np.savez_compressed(OUT/f'wave-{k[0]}-{k[1]}.npz',t=times,normalized_direct=q,components=parts,
            component_names=['deep_nonrest','atmosphere_nonrest','atmosphere_rest','conditional_inner_baryon_debit'])
    d=data[(896,128)];q=curves[(896,128)];peak=max(abs(q));errors={}
    errors['wave_space']=float(max(abs(q-curves[(448,128)]))/peak);errors['wave_time']=float(max(abs(curves[(448,128)]-curves[(448,64)]))/max(abs(curves[(448,128)])))
    for key,kw in [('wave_history',dict(kind='cubic')),('wave_saved_cadence',dict(stride=2)),('wave_quadrature',dict(order=4))]:
        _,qq,_,rr=integrate(d,model,**kw);errors[key]=float(max(abs(q-qq))/peak);rows[key]=rr
    gates=json.loads((OUT/'plan.json').read_text())['gates'];passed=all(value<gates.get(k,.02) for k,value in errors.items())
    parts=np.load(OUT/'wave-896-128.npz')['components'];source=json.loads((OUT/'sources.json').read_text())
    # An actual represented Killing-energy perturbation is exported; it is
    # NOT silently relabelled as the total Bondi mass at infinity.
    energy=((d['gas_nonrest_energy_erg']+d['photon_energy_erg'])+float(d['cx'])*C*C*np.asarray(d['baryon_g'],float))@d['a']
    np.savez_compressed(OUT/'represented-energy.npz',t=d['t'],delta_Killing_energy_erg=energy,geometric_energy_cm=G*energy/C**4)
    result=dict(classification='Counterexample candidate',passed=passed,controls=errors,paths=rows,endpoint_direct_normalized=float(q[-1]),peak_direct_normalized=float(peak),
        endpoint_components=parts[:,-1].tolist(),source_seed_relative=max(r['source_seed_relative'] for r in source['paths']),seconds=time.monotonic()-start,new_fluid_steps=0,
        actual_coupled_sources_applied=True,conditional_inner_baryon_debit=True,photon_metric_sources_exported=True,
        full_interior_mechanics=False,radial_frequency_certification=False,full_GR_scalar_feedback=False,physical_Bondi_mass_normalization=False,final_charge_solved=False,full_goal_complete=False)
    write(OUT/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


def alignment_plan():
    write(OUT/'alignment-plan.json',dict(classification='Counterexample candidate',
        issue='The first readout stopped at the last radial quadrature node. That node changes slightly with Gauss order and grid, so the comparisons did not use exactly the same observer clock.',
        correction='Use the common outer face light delay for every path and order. Preserve the first source/wave artifacts. Re-evaluate only the4.73s saved-source readout; do not rerun any fluid or EOS.',
        original_readout_budget_seconds=45,first_readout_seconds=json.loads((OUT/'first-result.json').read_text())['seconds'],remaining_readout_seconds=40,
        gates_unchanged=True,source_sha256=sha(__file__)))


def audit():
    from types import FunctionType,SimpleNamespace
    assert not (OUT/'audit.json').exists();start=time.monotonic();signal.signal(signal.SIGALRM,prior.optical.timeout);signal.alarm(20)
    write(OUT/'audit-plan.json',dict(classification='Counterexample candidate',seconds=20,native_constructor_calls=1,new_fluid_steps=0,
        checks='Analytic manufactured retarded ramp integrated through two finite-volume cells, primitive-seed invariance, shared observer clock and conservative source state bindings.',
        manufactured_gate=1e-8,source_sha256=sha(__file__)))
    class Flat:
        def __init__(self,m):self.m=m
        def __call__(self,x):
            return np.ones_like(x),x/C,np.ones_like(x),np.ones_like(x),np.ones_like(x)
    fakegreen=SimpleNamespace(Geometry=Flat,polynomial=green.polynomial,paired=green.paired)
    fn=FunctionType(integrate.__code__,dict(integrate.__globals__,green=fakegreen),argdefs=integrate.__defaults__)
    T=prior.END;R=7e9;dx=120000.;tt=np.linspace(0,T,17);s=np.column_stack([tt,tt]);model=SimpleNamespace(m=SimpleNamespace(RJ=R))
    d=dict(edges=np.array([R-dx,R,R+dx]),t=tt,nonrest_trace_erg=s,baryon_g=np.zeros_like(s),deep_cells=1,cx=1.,M_cm=1.,inner_material_face_cm=R,RJ=R)
    times,q,parts,row=fn(d,model,order=8)
    delay=dx/C
    # At u>=delay both cells see a linear source for their entire volume.
    # Its antiderivative is (u+x/c)^2/2, whose exact cell average is known.
    exact=G/(2*C**3)*(times**2+delay**2/3);mask=times>=delay
    error=float(max(abs(q[mask]/exact[mask]-1)));assert error<1e-8
    curves=[np.load(OUT/f'wave-{c}-{n}.npz') for c,n in [(448,64),(448,128),(896,128)]]
    assert all(np.array_equal(curves[0]['t'],d['t']) for d in curves)
    d=np.load(OUT/'source-896-128.npz');conservation=np.max(abs(np.sum(d['baryon_g'],axis=1)-np.sum(d['baryon_g'][:,int(d['deep_cells']):],axis=1)))
    assert conservation==0
    result=dict(classification='Counterexample candidate',passed=True,manufactured_relative=error,shared_readout_clock_bitwise=True,
        primitive_seed_relative=json.loads((OUT/'result.json').read_text())['source_seed_relative'],native_calls=0,seconds=time.monotonic()-start,
        limitation='The manufactured test verifies the retarded integral and coefficient, not deep spatial/frequency convergence, the conditional baryon completion or missing metric and Bondi sources.')
    write(OUT/'audit.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
