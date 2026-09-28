"""Counterexample candidate: reconcile saved mass and reference-energy ports.

Reuse the actual finite histories. No new fluid steps or fitted mass offset.
"""
from pathlib import Path
import json, resource, sys, time
import numpy as np
import sympy as sp
import bound_native_generated_scalar as previous
import solve_native_incident_reciprocal as coupled

OUT=Path('native-mass-energy166-work')
background=previous.background; inf=previous.inf
read,write,sha=previous.read,previous.write,previous.sha
C=inf.C; G=inf.G; LD=np.longdouble
CAPS=dict(prepare=30,ledger=90,work=120,audit=30)
TOTAL=sum(CAPS.values())
CAPS['cell']=60  # Reassigned from the unused90-second ledger allocation.
CAPS['commutator']=60


def symbolic():
    e,n,v,u,nu,zeta=sp.symbols('e n v u nu zeta')
    h=sp.symbols('h')
    # Occupation in this code represents count per REFERENCE phase measure.
    number_flux=(n+h*v)*(1+h*zeta)
    physical=number_flux*e*(1+h*(u+nu))
    ref=number_flux*e
    assert sp.expand(sp.diff(physical-ref,h).subs(h,0)-n*e*(u+nu))==0
    return dict(classification='Proven',passed=True,
        identity='For reference-cell photon counts and Eref=a0*epsilon, the instantaneous coordinate-energy port is Lphys=Lref+(u+nu)*L0 at first order. The speed term is already in Lref; no second area/volume factor is added.',
        boundary='This is instantaneous -p_t at the SAME radius. In a time-dependent metric it is not a conserved Killing energy along the whole ray; propagation work and scalar/matter flux remain necessary.',
        scope='Algebra for the declared coordinates, not a numerical mass-conservation verdict.')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    files=[Path(__file__),Path(coupled.__file__),Path(inf.incident.__file__),
        Path(background.__file__),Path(background.retained.constraints.__file__),
        Path(coupled.prior.base.old.matter.branch.base.__file__),
        previous.OUT/'result.json',previous.OUT/'audit.json',
        previous.previous.matched.OUT/'mass-fine.npz',
        inf.prior.EV/'coupled-128.npz',inf.prior.EV/'accepted-ports-128.npz']
    for n in [64,128]:
        files += [inf.BEFORE/f'gr/source-{n}-reference-128.npz',
                  inf.BEFORE/f'photons-precise/steps-{n}-reference-128.npz',
                  inf.BEFORE/f'material-analytic/steps-{n}-reference-128.npz',
                  inf.SELF/f'metric/metric-{n}-g8.npz']
    files += [inf.incident.FIELDS/f'born-g{q}.npz' for q in [4,8]]
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='b9ae9ab14',
        claim='Reconcile the same-interface mass constraint with actual saved photon/material reference-energy histories, canonical geometry subtraction, cell lapse weighting, photon frequency work and the physical boundary energy conversion.',
        decision='Identify whether a source/port conversion or a stored work/discretization mismatch accounts for the paired mass residual. Apply only a derived correction, never a fitted residual subtraction. If additional evolution is needed, first specify its missing equation.',
        controls=dict(algebra=1e-12,time=.02,quadrature=.002,mass_closure_fraction=.002),
        budget=dict(actions=CAPS,total_seconds=TOTAL,CPU_threads=1,virtual_GiB=3,new_fluid_steps=0,new_EOS_roots=0,new_rays=0),
        measured_basis='Previous generated-scalar input/setup and all actions6.68s; saved archive inspection1.9s. Two saved paths,17 background states and explicit quadratures capped at240s plus prepare/audit60s. No long integration forecast.',
        stop='No automatic clock/mesh/horizon enlargement or weakened gate. Preserve all historical verdicts and every inspected component. A bookkeeping identity is not physical closure.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)}))
    write(OUT/'symbolic.json',symbolic())


def weights(model,d,order):
    q=background.wave.base.flow.initial.Quadrature(d['edges'],order)
    _,_,a,B,_=model.geo(q.r.ravel()-model.model.m.RJ)
    w=q.r*q.r*B.reshape(q.r.shape)
    return q.h*((w*a.reshape(q.r.shape))@q.w)/(q.h*(w@q.w))


def ledger():
    background.initialize();m=background.wave.Response();rows=[]
    for n in [64,128]:
        d=np.load(inf.BEFORE/f'gr/source-{n}-reference-128.npz')
        z=np.load(inf.BEFORE/f'material-analytic/steps-{n}-reference-128.npz')
        p=np.load(inf.BEFORE/f'photons-precise/steps-{n}-reference-128.npz')
        met=np.load(inf.SELF/f'metric/metric-{n}-g8.npz')
        energy=d['gas_nonrest_energy_erg']+d['photon_energy_erg']+d['baryon_g'].astype(LD)*LD(d['cx'])*LD(C)**2
        measured=met['asymptotic_mass_residual_cm'].astype(LD)*LD(C)**4/LD(G)
        a=d['a'].astype(LD);mean=weights(m,d,8).astype(LD)
        ports=-d['inner_cumulative_energy_erg']+d['outer_cumulative_energy_erg']
        center=(energy*a).sum(1,dtype=LD)+ports
        mean_delta=(energy*(mean-a)).sum(1,dtype=LD)
        mean_mass=center+mean_delta
        ids=[np.argmin(abs(z['t']-t)) for t in d['t']]
        rest=LD(m.model.m.a0)*LD(d['cx'])*LD(C)**2
        state=z['history_scaled'][ids]
        material=(state[:,2]+rest*state[:,0])*LD(z['amplitude'])
        photon=p['moments'][:,0].astype(LD)
        stored=(material+photon).sum(1,dtype=LD)+ports
        canonical=(energy*a-material-photon).sum(1,dtype=LD)
        identity=float(np.max(abs(mean_mass-measured))/max(np.max(abs(measured)),LD('1e-290')))
        expanded=float(np.max(abs(stored+canonical+mean_delta-measured))/max(np.max(abs(measured)),LD('1e-290')))
        q4=weights(m,d,4).astype(LD)
        spatial=float(np.max(abs((energy*(mean-q4)).sum(1,dtype=LD)))/max(np.max(abs(measured)),LD('1e-290')))
        np.savez_compressed(OUT/f'ledger-{n}.npz',t=d['t'],reference_state_and_ports_erg=stored,
            canonical_counterterm_erg=canonical,cell_lapse_weighting_erg=mean_delta,
            reconstructed_port_erg=mean_mass,stored_GR_port_erg=measured,
            center_lapse=a,mean_lapse=mean,cell_energy=energy,
            material_reference_erg=material,photon_reference_erg=photon,
            radial_ports_erg=ports,collision_energy_erg=p['collision_transfer'][:,:,0])
        rows.append(dict(steps=n,identity_relative=identity,expanded_relative=expanded,quadrature_relative=spatial,
            reference_state_and_ports_erg=float(stored[-1]),canonical_counterterm_erg=float(canonical[-1]),
            cell_lapse_weighting_erg=float(mean_delta[-1]),constraint_port_erg=float(measured[-1]),
            photon_reference_frequency_work_erg=float(p['escape'][2]),spectral_ghost_energy_erg=float(p['escape'][1]),
            material_ledger_energy_erg=float((z['ledger_scaled'][2]+rest*z['ledger_scaled'][0])*LD(z['amplitude'])),
            material_discard_energy_erg=float((z['discard_scaled'][2]+rest*z['discard_scaled'][0])*LD(z['amplitude']))))
    result=dict(classification='Counterexample candidate',rows=rows,
        algebraic_reconstruction_passed=max(max(v['identity_relative'],v['expanded_relative']) for v in rows)<1e-12,
        physical_mass_closure_verified=False,full_goal_complete=False)
    write(OUT/'ledger.json',result);print(json.dumps(result),flush=True)
    assert result['algebraic_reconstruction_passed'],result


def work():
    inf.prior.initialize();driver=inf.incident.Driver(8);model=driver.model
    data=np.load(inf.prior.EV/'coupled-128.npz')
    t=np.load(inf.BEFORE/'gr/source-128-reference-128.npz')['t']
    ids=[np.argmin(abs(data['snapshot_t']-v)) for v in t]
    assert np.max(abs(data['snapshot_t'][ids]-t))<1e-18
    I=np.concatenate([data['snapshot_bulk_I'][ids],data['snapshot_I'][ids].sum(1)],axis=1)
    b=model.bulk;W=np.r_[b.W,model.W];area=np.r_[b.area[:-1],model.area]
    Eangle=np.sum(I*(b.d['num']*b.d['Einf']),axis=3,dtype=LD)*LD(4*np.pi)*W[None,:,None]*b.w
    outer=LD(4*np.pi*C)*area[-1]*np.sum(I[:,-1]*(b.w*b.mu)[None,:,None]*(b.d['num']*b.d['Einf']),axis=(1,2),dtype=LD)
    variance=b.mu2-b.mu*b.mu
    assert np.max(abs(variance-np.diff(b.edges_mu)**2/12))<1e-15
    rows=[];gamma=1-1/np.sqrt(2)
    for n in [64,128]:
        p=np.load(inf.BEFORE/f'photons-precise/steps-{n}-reference-128.npz');stage=p['accepted_angular_times']
        missing=[];original=[];correct=[];converted=[]
        for now in stage:
            j=np.clip(np.searchsorted(t,now,side='left')-1,0,len(t)-2);f=(now-t[j])/(t[j+1]-t[j])
            Er=(1-f)*Eangle[j]+f*Eangle[j+1]
            U,Ut,Ux=driver.wave(now,driver.xc);lamt=driver.zc['Phi']*Ut
            missing.append(-np.sum(Er*variance*lamt[:,None],dtype=LD))
            original.append(-np.sum(Er*b.mu**2*lamt[:,None],dtype=LD))
            correct.append(-np.sum(Er*b.mu2*lamt[:,None],dtype=LD))
            u0=driver.wave(now,np.array([0.]))[0][0];uo=driver.wave(now,driver.xout)[0]
            ext=np.sum(driver.eq.h*((2*driver.zout['Phi']*uo/(driver.rout*driver.zout['b'])).reshape(driver.eq.r.shape)@driver.eq.w),dtype=LD)
            ell=(driver.z0['alpha'][0]/driver.r0+driver.z0['Phi'][0])*u0-ext
            converted.append(((1-f)*outer[j]+f*outer[j+1])*ell)
        weights=LD(t[-1]/n)*np.tile([1-gamma,gamma],n)
        cols={k:np.cumsum(np.asarray(v,LD)*weights,dtype=LD) for k,v in dict(missing_angular_work=missing,old_radial_work=original,correct_radial_work=correct,outer_coordinate_energy_conversion=converted).items()}
        np.savez_compressed(OUT/f'work-{n}.npz',stage_t=stage,mu=b.mu,mu2=b.mu2,variance=variance,**cols)
        rows.append(dict(steps=n,**{k+'_erg':float(v[-1]) for k,v in cols.items()}))
    mass=np.load(previous.previous.matched.OUT/'mass-fine.npz')
    old=mass['matched_homogeneous_port_cm'][-1]*LD(C)**4/LD(G)
    result=dict(classification='Counterexample candidate',rows=rows,paired_mass_residual_erg=float(old),
        angular_work_corrected_residual_erg=float(old+cols['missing_angular_work'][-1]),
        no_fluid_or_transport_replay=True,physical_mass_closure_verified=False,full_goal_complete=False)
    write(OUT/'work.json',result);print(json.dumps(result),flush=True)


def cell():
    write(OUT/'cell-plan.json',dict(classification='Counterexample candidate',
        reason='Measured cell lapse weighting and angular-moment work are too small to explain the0.307erg paired residual. Test the spatial geometry counterterm at the same finite-volume interface before any replay.',
        prediction='The initial canonical energy removed by the GR operator is evaluated at quadrature points, whereas the exported add-back uses cell centers. Compare them without fitting any correction.',
        budget_seconds=60,original_total_seconds=TOTAL,new_fluid_steps=0,new_roots=0,
        bindings={str(p):sha(p) for p in [Path(__file__),OUT/'registered-producer.py',OUT/'work-producer.py',OUT/'ledger.json',OUT/'work.json']}))
    coupled.initialize(2);mat=coupled.Material(128,128);driver=mat.driver
    d=np.load(inf.BEFORE/'gr/source-128-reference-128.npz');old=np.load(OUT/'ledger-128.npz')
    model=background.wave.Response();rows=[]
    states=[mat.point(k) for k in range(17)]
    data=np.load(inf.prior.EV/'coupled-128.npz');ids=[np.argmin(abs(data['snapshot_t']-v)) for v in d['t']]
    I=np.concatenate([data['snapshot_bulk_I'][ids],data['snapshot_I'][ids].sum(1)],axis=1)
    b=mat.model.bulk;W=np.r_[b.W,mat.model.W]
    Er=np.sum(I*(b.d['num']*b.d['Einf'])*b.w[None,None,:,None],axis=(2,3),dtype=LD)*LD(4*np.pi)*W
    bg=np.array([v['Q'][2]+mat.rest*v['Q'][0] for v in states],LD)+Er
    for qorder in [4,8]:
        model.setup(d,qorder);z=model.z;r=model.r
        q=background.wave.base.flow.initial.Quadrature(d['edges'],qorder)
        _,_,a,B,_=model.geo(q.r.ravel()-model.model.m.RJ)
        vol=LD(4*np.pi)*(q.h[:,None]*q.w*q.r*q.r*B.reshape(q.r.shape)).astype(LD)
        v=vol.sum(1,dtype=LD);measure=(vol*a.reshape(q.r.shape)).ravel()
        out=[];init=[];bulk=[]
        for k,t in enumerate(d['t']):
            U=driver.wave(t,driver.optical(r))[0];u=z['alpha']*U/r;lam=z['Phi']*U;s=3*u+lam
            canonical=((z['Eg']+z['Pg'])*s+4*z['Er']*u+(z['Er']+z['Pr'])*lam)*LD(C)**4/LD(G)
            state=-np.repeat(bg[k]/d['a']/v,qorder)*s
            ca=np.sum(measure*canonical,dtype=LD);st=np.sum(measure*state,dtype=LD)
            init.append(ca);bulk.append(st);out.append(ca+st)
        center_counter=(old['cell_energy']*old['mean_lapse']-(old['material_reference_erg']+old['photon_reference_erg'])*old['mean_lapse']/old['center_lapse']).sum(1,dtype=LD)
        delta=np.asarray(out)-center_counter
        np.savez_compressed(OUT/f'cell-g{qorder}.npz',t=d['t'],initial_canonical_erg=init,
            physical_volume_erg=bulk,counterterm_erg=out,center_counterterm_erg=center_counter,
            spatial_counterterm_shift_erg=delta,proper_volume=v,stored_volume=d['volume'])
        rows.append(dict(order=qorder,spatial_counterterm_shift_erg=float(delta[-1]),
            continuous_counterterm_erg=float(out[-1]),center_counterterm_erg=float(center_counter[-1]),
            volume_relative=float(np.max(abs(v/d['volume']-1)))))
    result=dict(classification='Counterexample candidate',rows=rows,
        physical_mass_closure_verified=False,full_goal_complete=False)
    write(OUT/'cell.json',result);print(json.dumps(result),flush=True)


def commutator():
    write(OUT/'commutator-plan.json',dict(classification='Counterexample candidate',
        claim='Measure the actual mismatch between cell-centered continuum redshift work and the discrete shared-face energy flux product rule.',
        prediction='For reference photon counts the spatial frequency work must cancel sum ell_i*dE_i/dt plus the physical boundary-energy conversion. The continuum point derivative need not be the adjoint of the existing upwind finite-volume flux.',
        stop='No replay or alteration of saved trajectories. A measured source correction is not its coupled-response charge.',
        budget_seconds=60,total_previous_budget_seconds=TOTAL,
        bindings={str(p):sha(p) for p in [Path(__file__),OUT/'work.json',OUT/'cell.json']}))
    inf.prior.initialize();driver=inf.incident.Driver(8);model=driver.model
    d=np.load(inf.prior.EV/'coupled-128.npz');t=np.load(inf.BEFORE/'gr/source-128-reference-128.npz')['t']
    ids=[np.argmin(abs(d['snapshot_t']-v)) for v in t]
    I=np.concatenate([d['snapshot_bulk_I'][ids],d['snapshot_I'][ids].sum(1)],axis=1)
    b=model.bulk;W=np.r_[b.W,model.W];area=np.r_[b.area[:-1],model.area]
    density=np.sum(I*(b.d['num']*b.d['Einf']),axis=3,dtype=LD)*LD(4*np.pi)*b.w
    cc=driver.zc['lapse']*np.sqrt(driver.zc['b']);gamma=1-1/np.sqrt(2);rows=[]
    for n in [64,128]:
        p=np.load(inf.BEFORE/f'photons-precise/steps-{n}-reference-128.npz');stage=p['accepted_angular_times']
        point=[];face=[];checks=[]
        for now in stage:
            j=np.clip(np.searchsorted(t,now,side='left')-1,0,len(t)-2);f=(now-t[j])/(t[j+1]-t[j])
            q=(1-f)*density[j]+f*density[j+1];field=driver.at(now)
            ell=field['delta_log_lapse'];prime=field['delta_nu_prime']+field['delta_u_prime']
            flux=LD(C)*area[1:-1]*np.sum(np.where((b.mu>0)[None,:],q[:-1],q[1:])*b.mu,axis=1,dtype=LD)
            direct=-np.sum(q*W[:,None]*C*cc[:,None]*b.mu*prime[:,None],dtype=LD)
            discrete=-np.sum(np.diff(ell)*flux,dtype=LD)
            # Independent summation-by-parts evaluation, with zero end flux;
            # endpoint conversion is kept separately in work.json.
            rate=-np.diff(np.r_[LD(0),flux,LD(0)])
            checks.append(float(abs(discrete+np.sum(ell*rate,dtype=LD))/max(abs(discrete),LD('1e-290'))))
            point.append(direct);face.append(discrete)
        weights=LD(t[-1]/n)*np.tile([1-gamma,gamma],n)
        point=np.cumsum(weights*np.asarray(point,LD),dtype=LD)
        face=np.cumsum(weights*np.asarray(face,LD),dtype=LD)
        np.savez_compressed(OUT/f'commutator-{n}.npz',stage_t=stage,continuum_point_work_erg=point,
            discrete_face_work_erg=face,missing_product_rule_work_erg=face-point)
        rows.append(dict(steps=n,continuum_point_work_erg=float(point[-1]),discrete_face_work_erg=float(face[-1]),
            missing_product_rule_work_erg=float(face[-1]-point[-1]),summation_by_parts_relative=max(checks)))
    old=read(OUT/'work.json')['paired_mass_residual_erg']
    result=dict(classification='Counterexample candidate',rows=rows,
        paired_residual_before_erg=old,paired_residual_after_source_work_only_erg=old+rows[-1]['missing_product_rule_work_erg'],
        source_product_rule_passed=max(v['summation_by_parts_relative'] for v in rows)<1e-12,
        correction_applied_to_coupled_transport=False,physical_mass_closure_verified=False,full_goal_complete=False)
    write(OUT/'commutator.json',result);print(json.dumps(result),flush=True)
    assert result['source_product_rule_passed']


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS
    receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():
                target=OUT/'registered-producer.py' if p==str(Path(__file__)) and (OUT/'registered-producer.py').exists() else p
                assert sha(target)==h,p
        globals()[action]()
    except Exception as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(action=action,seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
