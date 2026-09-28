"""Apply the current interior metric to the actual incoming primary.

Counterexample candidate. This is an explicit missing scalar-operator return,
not the reciprocal mixed Einstein constraints or a complete stellar response.
Keep the exterior domain boundary paired with the interior boundary.
"""
from pathlib import Path
from types import FunctionType
import json, resource, sys, time
import numpy as np
import sympy as sp
from numpy.polynomial import legendre as leg
import def_retained_static_response as current
import def_retained_metric_fields as retained
import read_native_incident_infinity as infinity
import couple_native_exterior_incident as exterior

OUT=Path('native-interior-incident160-work');FIELDS=OUT/'background';METRIC=OUT/'metric'
read,write,sha=current.read,current.write,current.sha
wave=infinity.wave;C=infinity.C;G=infinity.G;LD=np.longdouble
CAPS=dict(prepare=30,background=180,coefficients=90,readout=120,audit=60)
SETTINGS=[(128,8), (128,4), ('128-17',8)]
PACKETS=[Path('retained-motion-return151-work/photons/steps-128-reference-128.npz'),
         Path('retained-metric-return152-work/photons/steps-128-reference-128.npz'),
         current.ACOUSTIC/'photons/steps-128-reference-128.npz']


def symbolic():
    r,b,m,A,N,Phi,alpha,Eg,Pg,Er,Pr,Kg,R4,f,J,eF,pF=sp.symbols(
        'r b m A N Phi alpha Eg Pg Er Pr Kg R4 f J eF pF',nonzero=True)
    E=Eg+Er;P=Pg+Pr;H=Eg+Pg;lam=r*Phi*f+J/(r*b)
    de=eF-H*(3*alpha*f+lam)-4*alpha*Er*f-(Er+Pr)*lam
    dp=pF-Kg*(3*alpha*f+lam)-4*alpha*Pr*f-(3*Pr-R4)*lam
    physical=2*lam/(r*b)+4*sp.pi*r*A/b*(dp-de+(4*alpha*f+2*lam)*(P-E))
    stable=(2/(r*b)+4*sp.pi*r*A/b*(3*Pg-Eg-Kg-Er+R4))*lam
    stable+=4*sp.pi*r*A/b*(alpha*(7*Pg-Eg-3*Kg)*f+pF-eF)
    assert sp.simplify(physical-stable)==0
    x,t,c=sp.symbols('x t c',positive=True);Z=sp.Function('Z')(t,x);g=sp.Function('g');U=g(t+x/c)
    derivative=2*Z*sp.diff(U,x,2)+sp.diff(Z,x)*sp.diff(U,x)+sp.diff(Z,t)*sp.diff(U,t)/c**2
    assert sp.simplify(derivative-sp.diff(Z*sp.diff(U,x),x)-sp.diff(Z*sp.diff(U,x),t)/c)==0
    return dict(classification='Proven',passed=True,
        constraint='zeta_r=[2/(r*b)+4*pi*r*A4/b*(3Pg-Eg-Kg-Er+R4)]*lambda +4*pi*r*A4/b*[alpha*(7Pg-Eg-3Kg)*f-(eF-pF)]. This is the first polar constraint on the declared momentarily balanced initial operator.',
        reduction='For incoming U=g(t+x/c), the derivative source is (partial_x+partial_t/c)(zeta*U_x). Interior volume density is (-2*zeta*Vgeom-a^2*zeta_r/r)*U in dx dt. Its outer boundary is +integral zeta(0,t)*U_x(0,t)dt; the exterior term is its negative.',
        scope='Exact identities for this selected operator. Reciprocal mixed scalar-metric sources, incident action on the generated background scalar, new matter feedback and full physical errors remain open.')


def prepare():
    assert not OUT.exists();FIELDS.mkdir(parents=True);METRIC.mkdir()
    files=[Path(__file__),Path(current.__file__),Path(retained.__file__),Path(wave.__file__),Path(wave.base.__file__),
           Path(infinity.__file__),Path(exterior.__file__),wave.base.flow.INPUT/'balanced-20.npz',
           current.BASE/'fields/source-128.npz',infinity.prior.EV/'accepted-ports-128.npz',
           infinity.BACKGROUND,infinity.OUT/'charge-parts.npz',exterior.OUT/'bound-result.json',
           exterior.OUT/'charge-128-a8-t8-r8.npz']+PACKETS
    files += [folder/'source-128-reference-128.npz' for folder in [current.BASE/'gr',current.ACOUSTIC/'gr-precision']]
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='e725cf704',
        previous_goal_turn='Progress: Phase159 applied the missing exterior interaction and bounded selected exterior terms, preserving the failed original point controls.',
        claim='Apply the current saved interior background metric to the actual incoming primary, pair the interface trace with the exterior domain, and measure the resulting null-infinity charge contribution.',
        decision='Does this previously omitted interior scalar-operator return change the stored matter charge sign or magnitude? A material contribution changes the next mixed-constraint/feedback priority; a failed control remains a failure, not a license to refine.',
        current_background='Use the Phase154 current-source composition: Phase152 full retained/native/motion sources plus its GR return and Phase153 native-acoustic return once. Build the newly combined scalar/metric field on the existing33 knots. Reuse actual original128 emission intervals and all three saved signed SDIRK packet histories.',
        operator='Initial polar constraint coefficients; current background f,lambda,zeta are separate increments. Apply -deltaL_background to the exact existing incident primary. Integrate the derivative divergence analytically and the remaining volume on its characteristic cuts; do not approximate the primary by saved waveform samples.',
        representation='Existing531 cells,33 source knots and4/8 radial reconstruction. The17-knot subset is an interpolation control, not a new evolution path. Five time Gauss points and degree-dependent spatial Gauss points integrate the declared polynomial times the degree8 pulse; no physical mesh or source degree is enlarged.',
        boundary='The full current background interface trace is used in both domain terms, which cancel exactly. Reuse the Phase159 failed exterior point only as an explicitly unaccepted displayed estimate; use its selected-term envelope when bounding that unknown correction.',
        exclusions=['Reciprocal mixed scalar-metric constraints and incident metric acting on generated background scalar.',
                    'Higher/Born incident field interaction and return of this new wave into photons/free matter.',
                    'Complete physical background, EOS/derivative/spatial/interpolation/nonlinear errors, static and observational comparison.'],
        settings=SETTINGS,gates=dict(quadrature=.002,cadence=.02,port=1e-12,boundary=1e-12,independent=1e-10),
        budgets=CAPS,total_seconds=sum(CAPS.values()),CPU_threads=1,virtual_GiB=3,new_EOS_calls=0,new_fluid_steps=0,
        forecast='Measure and save the first fine background field. Twice its cost for each remaining field plus10s must fit180s. Reuse that field. Measure the first fine charge integral before remaining controls within120s.',
        stop='Preserve raw failures; do not expand clocks, spatial degree, physical interval, repeated fluid trajectories or acceptance thresholds.',
        bindings={str(p):sha(p) for p in files}))
    write(OUT/'symbolic.json',symbolic());(OUT/'registered-producer.py').write_bytes(Path(__file__).read_bytes())


def initialize():
    infinity.prior.initialize();wave.OUT=FIELDS;wave.base.OUT=FIELDS
    retained.constraints.OUT=FIELDS


def background():
    initialize();d=current.source(128)
    np.savez_compressed(FIELDS/'source-128.npz',**d)
    coarse={k:(v[::2] if v.ndim and len(v)==len(d['t']) else v) for k,v in d.items()}
    np.savez_compressed(FIELDS/'source-128-17.npz',**coarse)
    model=wave.Response();fn=FunctionType(wave.Response.run.__code__,dict(wave.Response.run.__globals__,OUT=FIELDS))
    start=time.monotonic();rows=[fn(model,128,8)];first=time.monotonic()-start
    forecast=4*first+10;write(FIELDS/'pilot.json',dict(first_seconds=first,remaining_upper_seconds=forecast,eligible=forecast<CAPS['background']-first))
    assert forecast<CAPS['background']-first
    rows += [fn(model,n,q) for n,q in SETTINGS[1:]]
    write(FIELDS/'result.json',dict(classification='Counterexample candidate',fields_constructed=True,rows=rows,
        scope='The newly combined saved current source only; no physical background re-evolution or continuum certificate.'))


class Lapse(retained.metric.Lapse):
    def boundary(self,d,steps,order):
        t=d['t'];end=t[-1];mu,mw,ids,solution,invariant=self.rays(order,end)
        bg=np.load(infinity.prior.EV/'accepted-ports-128.npz');lum=bg['angular_luminosity'];h=float(bg['h'])
        gx,gw=leg.leggauss(order);gamma=1-1/np.sqrt(2);histories=[]
        for path in PACKETS:
            p=np.load(path);hp=end/int(path.stem.split('-')[1]);stage=p['accepted_angular_times']
            histories.append((stage,LD(hp)*np.tile([1-gamma,gamma],len(stage)//2)[:,None]*p['accepted_angular_luminosity']))
        photon=[];energy=[]
        for now in t:
            count=int(round(now/h))
            if not count:photon.append(0.);energy.append(0.);continue
            tt=((np.arange(count)[:,None]+(gx+1)/2)*h).ravel()
            weights=((h*gw/2)[None,:,None]*lum[:count,None,ids]*mu*mw).reshape(len(tt),len(mu)).astype(LD)
            for stage,packets in histories:
                take=stage<=now+16*np.finfo(float).eps*end
                tt=np.r_[tt,stage[take]];weights=np.concatenate([weights,packets[take][:,ids]*mu*mw])
            state=solution.sol(np.maximum((now-tt)/end,0.)).reshape(3,len(mu),len(tt)).transpose(0,2,1)
            rp=self.r0*state[0];mp=state[1];a=self.bg.metric(rp.ravel()/self.model.m.R)[3].reshape(rp.shape)
            photon.append(float(LD(G)/LD(C)**4*np.sum(weights*(state[2]/self.r0-mp*mp/(rp*a)),dtype=LD)))
            energy.append(float(np.sum(weights,dtype=LD)))
        photon=np.array(photon);energy=np.array(energy)
        error=float(np.max(abs(energy-d['outer_cumulative_energy_erg']))/max(np.max(abs(energy)),1.));assert error<1e-12,error
        z=(gx+1)/2;_,N,b,_,_=self.bg.metric(self.r0/(z*self.model.m.R))
        return photon,energy,dict(ray_invariant=invariant,emitted_energy_relative=error,
            unit_ADM_lapse_kernel_per_cm=float(np.sum(gw/(2*self.r0*N*b**1.5))),all_current_packet_histories_applied=True)


def coefficients():
    initialize();model=Lapse();fn=FunctionType(retained.metric.Lapse.run.__code__,dict(retained.metric.Lapse.run.__globals__,OUT=METRIC))
    rows=[]
    for n,qorder in SETTINGS:
        row=fn(model,n,qorder);d=dict(np.load(FIELDS/f'source-{n}.npz'));field=np.load(FIELDS/f'fields-{n}-g{qorder}.npz')
        metric=np.load(METRIC/f'metric-{n}-g{qorder}.npz');m=model.response;m.setup(d,qorder)
        q=wave.base.flow.initial.Quadrature(d['edges'],qorder);r=m.r;z=m.z;a=z['lapse']*np.sqrt(z['b'])
        f=np.array([np.interp(r,field['radius_E'],v) for v in field['delta_phi']])
        lam=r*z['Phi']*f+m.J/(r*z['b'])
        stress=np.asarray(d['metric_stress_erg']/d['volume']*LD(G)/LD(C)**4,float)[:,m.ids]
        zr=(2/(r*z['b'])+4*np.pi*r*z['A4']/z['b']*(3*z['Pg']-z['Eg']-z['Kg']-z['Er']+z['R4']))*lam
        zr+=4*np.pi*r*z['A4']/z['b']*(z['alpha']*(7*z['Pg']-z['Eg']-3*z['Kg'])*f-stress)
        centers=retained.constraints.centers(m,d,field,qorder)
        outer_lam=centers['delta_lambda'][:,-1];outer=metric['delta_nu_faces'][:,-1]-outer_lam
        z0=m.coeff(np.array([model.r0]));a0=float(z0['lapse'][0]*np.sqrt(z0['b'][0]))
        radiation=metric['outer_photon_lapse']+G/C**4*metric['emitted_energy_erg']/(model.r0*a0)
        scalar=metric['outer_scalar_lapse']-z0['Phi'][0]*field['U'][:,-1]
        residual=metric['outer_ADM_residual_lapse']-metric['asymptotic_mass_residual_cm']/(model.r0*a0)
        external=radiation+scalar+residual
        interface_error=float(np.max(abs(outer-external))/max(np.max(abs(outer)),1e-290))
        assert interface_error<1e-12,interface_error
        jac=1/(np.exp(-2*z['phi']**2)*(1+z['alpha']*r*z['Phi']))
        zet=[]
        for boundary,gradient in zip(outer,zr):
            integ,faces=q.integrate(gradient*jac)
            zet.append((boundary-(faces[-1]-integ)).ravel())
        zet=np.array(zet);vg=z['lapse']**2*(2*z['mass']/r**3+4*np.pi*z['A4']*(z['Pg']+z['Pr']-z['Eg']-z['Er']))
        density=-2*zet*vg-a*a*zr/r
        # ponytail: retain the existing node/center field interpolation. Its
        # continuum error stays open; no new reconstruction is invented here.
        np.savez_compressed(OUT/f'operator-{n}-g{qorder}.npz',t=d['t'],x=m.x-m.xfaces[-1],
            xfaces=m.xfaces-m.xfaces[-1],radius=r,zeta=zet,zeta_prime=zr,volume_density=density,
            outer_zeta=outer,external_zeta=external,outer_radiation_zeta=radiation,
            outer_scalar_zeta=scalar,outer_ADM_residual_zeta=residual,interface_relative=interface_error,
            M=float(d['M_cm']),r0=model.r0,inverse=m.inverse)
        rows.append(dict(row,interface_relative=interface_error))
    write(OUT/'coefficients.json',dict(classification='Counterexample candidate',constructed=True,rows=rows,
        constraint_identity=read(OUT/'symbolic.json')['passed'],actual_current_background_used=True,full_goal_complete=False))


def integrate(z,clock,spatial_extra=0):
    """Exact time/pulse cuts for the declared linear-time nodal polynomial."""
    t=z['t'];faces=z['xfaces'];order=z['inverse'].shape[-1];mid=(faces[:-1]+faces[1:])/2;half=np.diff(faces)/2
    co=np.einsum('cij,tcj->tci',z['inverse'],z['volume_density'].reshape(len(t),-1,order))
    gx,gw=leg.leggauss((order+10)//2+spatial_extra);tx,tw=leg.leggauss(5);D=t[-1]/2;eta=infinity.incident.ETA
    boundary=[];exterior_boundary=[];bulk=[]
    for u in clock:
        if u==0:boundary.append(0.);exterior_boundary.append(0.);bulk.append(0.);continue
        cuts=np.unique(np.r_[faces,-C*t,C*(t-u),C*(D-t),C*(D-u)/2,-C*u/2])
        cuts=cuts[(cuts>=max(faces[0],-C*u/2))&(cuts<=0)]
        lo=cuts[:-1];hi=cuts[1:];ids=np.clip(np.searchsorted(faces,(lo+hi)/2,side='right')-1,0,len(mid)-1)
        xx=((lo+hi)[:,None]/2+(hi-lo)[:,None]*gx/2).ravel();ww=((hi-lo)[:,None]*gw/2).ravel();ids=np.repeat(ids,len(gx))
        vv=leg.legvander((xx-mid[ids])/half[ids],order-1)
        values=np.einsum('tij,ij->ti',co[:,ids],vv)
        lower=np.maximum(t[:-1,None],-xx/C);upper=np.minimum(t[1:,None],np.minimum(D-xx/C,u+xx/C))
        width=np.maximum(upper-lower,0);qt=lower[...,None]+width[...,None]*(tx+1)/2
        v=values[:-1,:,None]+(qt-t[:-1,None,None])/np.diff(t)[:,None,None]*(values[1:]-values[:-1])[:,:,None]
        pulse=eta*float(z['r0'])*infinity.incident.pulse((qt+xx[None,:,None]/C)/D)
        bulk.append(C/2*np.sum(ww[None,:]*width/2*np.sum(tw*v*pulse,axis=-1),dtype=LD))
        lo=t[:-1];hi=np.minimum(t[1:],min(u,D));width=np.maximum(hi-lo,0);qt=lo[:,None]+width[:,None]*(tx+1)/2
        value=z['outer_zeta'][:-1,None]+(qt-lo[:,None])/np.diff(t)[:,None]*np.diff(z['outer_zeta'])[:,None]
        ux=eta*float(z['r0'])/(C*D)*infinity.incident.pulse(qt/D,1)
        boundary.append(C/2*np.sum(width[:,None]*tw/2*value*ux,dtype=LD))
        value=z['external_zeta'][:-1,None]+(qt-lo[:,None])/np.diff(t)[:,None]*np.diff(z['external_zeta'])[:,None]
        exterior_boundary.append(-C/2*np.sum(width[:,None]*tw/2*value*ux,dtype=LD))
    return np.asarray(bulk,LD),np.asarray(boundary,LD),np.asarray(exterior_boundary,LD)


def readout():
    clock=np.load(infinity.OUT/'charge-parts.npz')['t'];paths=[];rows=[];start=time.monotonic()
    for i,(n,q) in enumerate(SETTINGS):
        z=np.load(OUT/f'operator-{n}-g{q}.npz');mark=time.monotonic();bulk,boundary,external=integrate(z,clock)
        charge=-bulk/LD(z['M']);edge=-boundary/LD(z['M']);ext=-external/LD(z['M'])
        error=float(np.max(abs(edge+ext))/max(np.max(abs(edge)),LD('1e-290')));assert error<1e-12,error
        np.savez_compressed(OUT/f'charge-{n}-g{q}.npz',t=clock,interior_volume=charge,
            interior_domain_boundary=edge,matched_exterior_domain_boundary=ext,paired_total=charge+edge+ext)
        paths.append(charge+edge+ext);rows.append(dict(source=n,order=q,seconds=time.monotonic()-mark,boundary_relative=error))
        if not i:
            forecast=4*rows[-1]['seconds']+10
            write(OUT/'readout-pilot.json',dict(first_seconds=rows[-1]['seconds'],remaining_upper_seconds=forecast,eligible=forecast<CAPS['readout']-(time.monotonic()-start)))
            assert forecast<CAPS['readout']-(time.monotonic()-start)
    norm=max(float(np.max(abs(paths[0]))),1e-290)
    controls={name:float(np.max(abs(v-paths[0]))/norm) for name,v in zip(['quadrature','cadence'],paths[1:])}
    old=np.load(infinity.OUT/'charge-parts.npz');e0=np.load(infinity.BACKGROUND)['epsilon'].astype(LD)
    correction=paths[0]/(1-e0);combined=old['body']+correction
    np.savez_compressed(OUT/'applied-charge.npz',t=clock,old_body=old['body'],interior_operator_return=correction,
        body_plus_selected_return=combined,direct=old['direct'])
    result=dict(classification='Counterexample candidate',passed=controls['quadrature']<.002 and controls['cadence']<.02,
        controls=controls,endpoint_interior_operator_return=float(correction[-1]),endpoint_previous_body=float(old['body'][-1]),
        endpoint_combined_selected=float(combined[-1]),return_over_previous_body=float(correction[-1]/old['body'][-1]),
        current_interior_metric_applied_to_incident_primary=True,domain_boundary_paired=True,
        full_reciprocal_mixed_constraints=False,new_wave_applied_to_material_photons=False,full_goal_complete=False,
        rows=rows,seconds=time.monotonic()-start)
    write(OUT/'result.json',result);print(json.dumps(result),flush=True);assert result['passed'],result


if __name__=='__main__':
    action=sys.argv[1];assert action in ['prepare','background','coefficients','readout']
    receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));infinity.incident.native.deadline(CAPS[action])
    started=time.monotonic();cpu=time.process_time();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        globals()[action]()
    except Exception as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(action=action,seconds=time.monotonic()-started,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
