"""Counterexample candidate: reciprocal matter/photon response to incident GR.

Same retained finite background, external packet and clocks as Phase155.
Solve the existing 17-knot waveform equations; no nonlinear/continuum claim.
"""
from pathlib import Path
from types import FunctionType
import json,resource,sys,time
import numpy as np
import sympy as sp
import solve_native_incident_material as prior
import def_native_matter_photon_feedback as feedback

base=prior.base;AMP=base.AMP;LD=np.longdouble
OUT=Path('native-incident-reciprocal156-work');read,write,sha=prior.read,prior.write,prior.sha
CAPS=dict(prepare=30,check=90,photon_pilot=120,photon_production=1100,
          material_pilot=60,material_production=400,residual=60,compact=180)
MAX_SWEEPS=3
CAPS['photon_repair']=79
CAPS['material_repair']=50
CAPS['material_deep_pilot']=45


def paths(sweep):
    if sweep==0:return prior.lift.PHOTON,prior.MATERIAL
    folder=OUT/f'sweep-{sweep}'
    return folder/'photons-precise',folder/('material-probed' if sweep==1 else 'material-analytic')


def transfers(m,path,n):
    with np.load(path/f'steps-{n}-reference-128.npz') as p:
        ids=[int(np.argmin(abs(p['t']-t))) for t in m.t]
        assert np.max(abs(p['t'][ids]-m.t))<1e-18
        c=p['collision_transfer'][ids]
        m.transfer=np.stack([np.zeros_like(c[:,:,0]),p['moments'][ids,3]/m.a,c[:,:,0],c[:,:,1]],axis=1)/AMP


def initialize(sweep):
    global Response,Material,PHOTON,MATERIAL
    assert 1<=sweep<=MAX_SWEEPS
    PHOTON,MATERIAL=paths(sweep);previous_photon,previous_material=paths(sweep-1)
    prior.initialize();Photon=base.Response;Matter=base.Material
    class Reciprocal(Photon):
        def __init__(self,n):
            super().__init__(n);self.material=Matter(128,n)
            transfers(self.material,previous_photon,n)
            with np.load(previous_material/f'steps-{n}-reference-128.npz') as d:
                ids=[int(np.argmin(abs(d['t']-t))) for t in self.t]
                assert np.max(abs(d['t'][ids]-self.t))<1e-18
                motion=d['history_scaled'][ids].copy()
            offset=(self.a.astype(LD)-self.model.m.a0)*self.model.cx*LD(base.C)**2*motion[:,0]
            self.energy_offset=np.asarray(offset,float);motion[:,2]-=offset
            self.motion=np.asarray(motion,float)
            self.mechanical=np.asarray(motion[:,[2,3]]-self.material.transfer[:,[2,3]],float)
            self.xi=np.array([self.material.model.mech.xi@np.r_[0.,-np.cumsum(z[0,:self.nb])] for z in self.motion])
            self.velocity_jet_error=0.;self.mapping_error=0.;self.map_checks=[]
        point=feedback.Response.point
        def collision(self,c,x,g,source=False):
            p,q,e,b=feedback.mono.Response.collision(self,c,x,g,source)
            return p,q+c['mechanical'] if source else q,e,b
        def local(self,t):
            self.set_stage(t);c=feedback.Response.local(self,t)
            p,_,e,b=self.collision(c,self.lift(t)[0],np.zeros((self.n,2)))
            c['q']+=p;c['qb']+=b;c['qe']+=e
            return c
    # Reuse the accepted lifted, actual-stage SDIRK owner. The gas unknown is
    # now total Etilde/H, while only collision transfers go to the free fluid.
    runner=(prior.OUT/'expanded-lifted-run.py').read_text()
    runner=runner.replace('rtol=1e-12','rtol=1e-14')
    def change(a,b):
        nonlocal runner
        assert runner.count(a)==1,(a,runner.count(a));runner=runner.replace(a,b)
    change('impulse=np.zeros(self.n);','impulse=np.zeros(self.n);transfer=np.zeros((self.n,2));')
    change("pressure=np.einsum('nj,nj->n',pmap,g)*self.volume",
           "pressure=(np.einsum('nj,nj->n',pmap,g)+self.local(t)['pressure_source'])*self.volume")
    change('transfer_history.append(g*np.stack([self.eu,self.nu],axis=-1)*AMPLITUDE)',
           'transfer_history.append(transfer.copy()*AMPLITUDE)')
    change('g[:,0]*self.eu,g[:,1]*self.nu,impulse',
           'g[:,0]*self.eu+np.array([np.interp(t,self.t,v) for v in self.energy_offset.T]),g[:,1]*self.nu,impulse')
    change("transfer_history=list(z['collision_transfer']);", "transfer_history=list(z['collision_transfer']);transfer=z['collision_transfer'][-1]/AMPLITUDE;")
    change('t=k*h;c=self.local(t+gamma*h);inverse=', 't=k*h;c=self.local(t+gamma*h);mechanical=c["mechanical"].copy();inverse=')
    change("qg=self.gas(c['q'],c['qb'],c['qe'])", "qg=self.gas(c['q'],c['qb'],c['qe'])+mechanical")
    change('c=self.local(t+h);inverse=', 'c=self.local(t+h);c["mechanical"]=mechanical;inverse=')
    change('qg=self.gas(c["q"],c["qb"],c["qe"])','qg=self.gas(c["q"],c["qb"],c["qe"])+mechanical')
    anchor='p2,g2,es2,_=self.collision(c,z,gz,True)'
    change(anchor,anchor+'\n        transfer+=h*((1-gamma)*g1+gamma*g2-mechanical)*np.stack([self.eu,self.nu],axis=-1)\n        ledger+=h*np.array([-np.sum(mechanical[:,1]*self.nu),np.sum(mechanical[:,0]*self.eu)])')
    change('port_history=port_history,transfer_history=transfer_history,','port_history=port_history,transfer_history=transfer_history,transfer=transfer,')
    change('collision_transfer=transfer_history,accepted_angular_times=',
           'collision_transfer=transfer_history,energy_offset_reference=AMPLITUDE*self.energy_offset,energy_offset_t=self.t,accepted_angular_times=')
    change("    row['passed']=", "    actual=g[:,0]*self.eu+np.array([np.interp(count*h,self.t,v) for v in self.energy_offset.T])\n    row.update(endpoint_material_nonrest_energy_erg=row['endpoint_material_reference_energy_erg'],endpoint_material_reference_energy_erg=float(np.sum(actual)*AMPLITUDE),endpoint_material_energy_L1_erg=float(np.sum(abs(actual))*AMPLITUDE),velocity_jet_relative=self.velocity_jet_error,primitive_mapping_relative=self.mapping_error,conserved_material_input_applied=True)\n    row['passed']=self.velocity_jet_error<1e-4 and self.mapping_error<1e-10 and ")
    namespace=dict(Photon.run.__globals__,OUT=PHOTON,write=write)
    exec(compile(runner,__file__,'exec'),namespace);Reciprocal.run=namespace['run'];Response=Reciprocal
    class Returned(Matter):
        def __init__(self,reference=128,steps=128):
            super().__init__(reference,steps);transfers(self,PHOTON,steps)
        # Arithmetic probes only: same physical eta. The initial common probe
        # under-resolved deep changes when the new atmospheric response grew.
        def rhs(self,t,z,probe=1.):return super().rhs(t,z,probe*64)
        run_owner=FunctionType(Matter.run_owner.__code__,dict(Matter.run_owner.__globals__,OUT=MATERIAL),argdefs=Matter.run_owner.__defaults__)
        def run(self,steps,label,limit=None,restart=None):
            self.run_owner(steps,label,limit,restart);return read(MATERIAL/f'{label}.json')
    if sweep>=2:
        import def_native_incident_deep_tangent as deep
        Returned.deep_tangent=deep.deep_tangent
        Returned.directional_rhs,source=deep.rhs(base.old.aligned)
        (OUT/'expanded-deep-tangent-rhs.py').write_text(source)
    Material=Returned
    (OUT/f'sweep-{sweep}/expanded-run.py').write_text(runner)


def prepare():
    assert not OUT.exists();OUT.mkdir()
    assert read(prior.OUT/'audit.json')['passed']
    for i in range(1,MAX_SWEEPS+1):
        for p in paths(i):p.mkdir(parents=True)
    B,E,a,a0,c,T,M=sp.symbols('B E a a0 c T M')
    Et=E-(a-a0)*c*B
    assert sp.expand(Et+(a-a0)*c*B-E)==0
    assert sp.expand((T+M)-M-T)==0
    write(OUT/'plan.json',dict(classification='Counterexample candidate',before_checkpoint='22f56148dc20d0e82f3c41ae8d28e22bf13ae959',
        claim='Return actual driven baryon/momentum/inventory and noncollisional Etilde/H to simultaneous photons and total gas Etilde/H; return only the new collision transfers to free matter. Test reciprocal finite waveform residual, then read its compact charge.',
        decision='Whether the original one-way incident response survives reciprocal material-photon coupling. Stop once both paths satisfy declared residual; otherwise at three sweeps or nondecreasing residual.',
        equations='Etilde=Eref-(a_ref-a_surface)*cx*c^2*B. Mechanical cumulative M=Etilde/H from prior free matter minus its collision-only transfer. New gas evolves C(new photons,new Etilde/H;prior B,S,xi)+dM/dt. Free matter receives only C. Exact driver/lift remains unchanged.',
        interpolation='Same17 canonical waveform knots; mechanical derivative uses the left interval at closing stages. A finite waveform residual is not a continuum contraction or error certificate.',
        budget=dict(per_action_seconds=CAPS,max_sweeps=MAX_SWEEPS,total_action_seconds=5600,CPU_threads=1,virtual_GiB=3),
        forecast='Previous full lifted photons233.78s and material166.75s. New moving coefficient maps add17 points per clock. Equal-horizon4/8 photon and two-step material prefixes measure admission;2x forecast with previous late-step floor must fit1100s/400s. Later contention and Krylov cost remain extrapolated.',
        gates=dict(time=.02,reciprocal=.002,energy_H=.002,conservation=1e-8,linear=1e-12,mapping=1e-10,velocity=1e-4,material_probe=.002,branch=.01),
        stop='Any original gate, nondecreasing full material residual after the second sweep, three sweeps or5600s. No automatic refinement, new waveform, amplitude, horizon or background replay.',
        scope='First-order retained constitutive maps on the corrected evolving background, prescribed external GR only. No native uniform derivatives, full nonlinear star, self-GR fixed point, complete infinity, orbital or static-nonabsorption claim.',
        symbolic=dict(classification='Proven',passed=True,scope='Invertible conserved energy coordinate and exact collision/mechanical partition only.'),
        bindings={str(p):sha(p) for p in [Path(__file__),Path(prior.__file__),Path(prior.lift.__file__),Path(base.__file__),Path(feedback.__file__),prior.OUT/'audit.json',prior.OUT/'expanded-lifted-run.py']}))


def check(sweep):
    initialize(sweep);m=Response(128);rows=[]
    for k in [0,8,16]:
        c=m.local(m.t[k]);z=np.zeros_like(m.I[0]);g=np.zeros((m.n,2))
        p,q,e,b=m.collision(c,z,g,True);paired=m.gas(p,b,e)
        err=float(np.max(abs(q-paired-c['mechanical']))/max(np.max(abs(q)),1.))
        energy=float(abs(np.sum(p*m.Eweight)+e[1].sum()+np.sum(paired[:,0]*m.eu))/max(np.sum(abs(p)*m.Eweight),np.sum(abs(paired[:,0])*m.eu),1.))
        species=float(abs(np.sum(b*m.Nweight)-np.sum(paired[:,1]*m.nu))/max(np.sum(abs(b)*m.Nweight),1.))
        partition=float(np.max(abs(m.motion[k,[2,3]]-m.mechanical[k]-m.material.transfer[k,[2,3]]))/max(np.max(abs(m.motion[k,[2,3]])),1.))
        rows.append(dict(k=k,mechanical_partition=partition,paired_gas=err,energy=energy,species=species))
        assert max(err,energy,species,partition)<1e-10,rows[-1]
    assert m.velocity_jet_error<1e-4 and m.mapping_error<1e-10
    write(OUT/f'sweep-{sweep}/check.json',dict(classification='Counterexample candidate',passed=True,rows=rows,velocity_jet=m.velocity_jet_error,mapping=m.mapping_error))


def photon(sweep,pilot):
    assert read(OUT/f'sweep-{sweep}/check.json')['passed'];initialize(sweep);start=time.monotonic();rows=[]
    if not pilot:assert read(PHOTON/'pilot.json')['eligible']
    for n in [64,128]:
        label=f'pilot-{n}' if pilot else f'steps-{n}-reference-128'
        if pilot and (PHOTON/f'{label}.json').exists():
            row=read(PHOTON/f'{label}.json');assert row['passed'];rows.append(row);continue
        m=Response(n)
        rows.append(m.run(n,label,n//16 if pilot else None,None if pilot else f'pilot-{n}'))
    def history(n):return np.load(PHOTON/(f'pilot-{n}.npz' if pilot else f'steps-{n}-reference-128.npz'))['moments'][:,[0,1,2,3,5,6]]
    a,b=history(64),history(128)
    errors=(np.max(np.sum(abs(a-b),axis=2),axis=0)/np.maximum(np.max(np.sum(abs(b),axis=2),axis=0),1e-290)).tolist()
    result=dict(classification='Counterexample candidate',passed=all(r['passed'] for r in rows) and max(errors)<.02,rows=rows,time_comparison=errors,seconds=time.monotonic()-start)
    if pilot:
        forecasts=[]
        for row in rows:
            old=read(prior.lift.PHOTON/f"steps-{row['steps']}-reference-128.json")
            cost=max(row['stepping_seconds']/row['new_steps'],old['stepping_seconds']/old['new_steps'])
            forecasts.append(17*row['operator_point_seconds']/row['operator_points']+(row['steps']-row['completed_steps'])*cost+20)
        result.update(upper_remaining_seconds=2*sum(forecasts),eligible=result['passed'] and 2*sum(forecasts)<CAPS['photon_production'])
    write(PHOTON/('pilot.json' if pilot else 'result.json'),result);print(json.dumps(result),flush=True)
    assert result['eligible' if pilot else 'passed'],result


def material(sweep,pilot):
    initialize(sweep);assert read(PHOTON/'result.json')['passed'];start=time.monotonic();rows=[]
    if not pilot:assert read(MATERIAL/'pilot.json')['eligible']
    for n in [64,128]:
        label=f'pilot-{n}' if pilot else f'steps-{n}-reference-128'
        if pilot and (MATERIAL/f'{label}.json').exists():
            row=read(MATERIAL/f'{label}.json');assert row['passed'];rows.append(row);continue
        start_one=time.monotonic();m=Material(128,n)
        row=m.run(n,label,2 if pilot else None,None if pilot else f'pilot-{n}')
        row.update(worker_wall_seconds=time.monotonic()-start_one,physical_branch_ratio=m.physical_branch_ratio,maximum_owner_error=max(v['owner_error'] for v in m.cache.values()))
        if not pilot:
            z=np.load(MATERIAL/f'{label}.npz')['delta_scaled'];rr=[m.rhs(m.t[-1],z,p)[0] for p in [.5,1.,2.]]
            row['endpoint_probe_half_nominal_double']=[(np.sum(abs(r-rr[1]),axis=1)/np.maximum(np.sum(abs(rr[1]),axis=1),1.)).astype(float).tolist() for r in [rr[0],rr[2]]]
            row['passed']=row['passed'] and np.max(row['endpoint_probe_half_nominal_double'])<.002
        row['nominal_arithmetic_probe']=512
        row['passed']=row['passed'] and row['physical_branch_ratio']<.01 and row['maximum_owner_error']<1e-8
        write(MATERIAL/f'{label}.json',row);rows.append(row);assert row['passed'],row
    result=dict(classification='Counterexample candidate',passed=True,rows=rows,seconds=time.monotonic()-start)
    if pilot:
        forecasts=[]
        for r in rows:
            old=read(prior.MATERIAL/f"steps-{r['steps']}-reference-128.json")
            forecasts.append(old['raw_owner_calls']*r['seconds']/r['raw_owner_calls']+r['worker_wall_seconds']-r['seconds']+5)
        result.update(upper_remaining_seconds=2*sum(forecasts),eligible=2*sum(forecasts)<CAPS['material_production'])
    write(MATERIAL/('pilot.json' if pilot else 'production.json'),result);print(json.dumps(result),flush=True)
    assert result['eligible' if pilot else 'passed'],result


def relative(a,b):
    return (np.max(np.sum(abs(a-b),axis=2),axis=0)/np.maximum(np.max(np.sum(abs(b),axis=2),axis=0),1e-290)).astype(float).tolist()


def residual(sweep):
    initialize(sweep);rows=[];histories=[]
    for n in [64,128]:
        d=np.load(MATERIAL/f'steps-{n}-reference-128.npz');p=np.load(PHOTON/f'steps-{n}-reference-128.npz')
        old=np.load(paths(sweep-1)[1]/f'steps-{n}-reference-128.npz')
        assert np.array_equal(d['t'],old['t'])
        ids=[int(np.argmin(abs(d['t']-t))) for t in p['t']]
        assert np.max(abs(d['t'][ids]-p['t']))<1e-18
        actual=d['history_scaled'];state=relative(old['history_scaled'],actual)
        ph=p['moments'][:,[1,2]]/AMP;mat=actual[ids][:,[2,3]]
        mismatch=relative(ph,mat);histories.append(actual[ids])
        balance=float(np.max(abs(np.sum(actual,axis=2,dtype=LD)+d['discards_scaled']-d['ledgers_scaled'])/np.maximum(d['norms_scaled'],1.)))
        gamma=1-1/np.sqrt(2);h=d['t'][-1]/n
        flux=(p['accepted_angular_luminosity']@(np.arange(1,8,2)/32)).reshape(n,2)
        angular=float(abs(h*np.sum(flux*[1-gamma,gamma])-p['radial_ports'][-1,1,1])/max(h*np.sum(abs(flux)),1e-290))
        row=dict(steps=n,material_waveform_residual=state,photon_material_energy_H_residual=mismatch,material_balance=balance,angular_port=angular)
        rows.append(row);assert balance<1e-8 and angular<1e-12,row
    timing=relative(*histories);maximum=max(max(r['material_waveform_residual']+r['photon_material_energy_H_residual']) for r in rows)
    before=read(OUT/f'sweep-{sweep-1}/residual.json')['maximum_residual'] if sweep>1 else None
    result=dict(classification='Counterexample candidate',rows=rows,time_comparison=timing,maximum_residual=maximum,previous_residual=before,
        finite_reciprocal_residual_passed=maximum<.002,passed=max(timing)<.02,
        next_sweep_allowed=maximum>=.002 and sweep<MAX_SWEEPS and (before is None or maximum<before),
        full_nonlinear_fixed_point=False,self_GR_returned=False,full_goal_complete=False)
    write(OUT/f'sweep-{sweep}/residual.json',result);print(json.dumps(result),flush=True);assert result['passed'],result


if __name__=='__main__':
    action=sys.argv[1];sweep=int(sys.argv[2]) if len(sys.argv)>2 else 0
    assert action in CAPS
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));base.drive.native.deadline(CAPS[action])
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        if action=='prepare':prepare()
        else:
            plan=read(OUT/'plan.json')
            for p,h in plan['bindings'].items():
                target=OUT/'registered-producer.py' if Path(p).name==Path(__file__).name else p
                assert sha(target)==h,p
            repair=read(OUT/'deep-tangent-plan.json');assert sha(__file__)==repair['source_sha256']
            assert sha('verification/def_native_incident_deep_tangent.py')==repair['deep_source_sha256']
            assert sha(OUT/'probed-producer.py')==read(OUT/'material-repair-plan.json')['source_sha256']
            assert sha(OUT/'linear-producer.py')==read(OUT/'linear-repair-plan.json')['source_sha256']
            if sweep>1:assert read(OUT/f'sweep-{sweep-1}/residual.json')['next_sweep_allowed']
            receipts=[read(p) for p in OUT.rglob('*-receipt.json')]
            assert sum(p['seconds'] for p in receipts)+CAPS[action]<=plan['budget']['total_action_seconds']
            if sweep>=2 and action.startswith('material_'):assert read(OUT/'deep-tangent-check.json')['passed']
            if action.startswith(('photon_','material_')):globals()[action.split('_')[0]](sweep,action.endswith(('pilot','repair')))
            else:globals()[action](sweep)
    except Exception as exc:error=repr(exc);raise
    finally:
        receipt=OUT/f'{action}-{sweep}-receipt.json';assert not receipt.exists()
        write(receipt,dict(seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,error=error,source_sha256=sha(__file__)))
