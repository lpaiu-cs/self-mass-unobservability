"""Evolve both sides of the old atmosphere interface with the same reactions.

Counterexample candidate: a finite reactive interior patch on the saved metric;
not a replacement for the entire stellar interior or evolved photon transport.
"""
from pathlib import Path
from types import FunctionType
import json
import signal
import sys
import time
import numpy as np
import def_native_reactive_flow as old

ex=old.exchange
C=old.C
OUT=ex.OUT.parent/'def-native-reactive-interface'
write=old.write
thermo=old.old.previous.old.bank.prior

# Reuse the exact background construction with only the inner domain edge
# moved. The measured atmosphere is still[-200,1200]m on the same cell sizes.
base_source=old.replace(thermo.init_source,'np.linspace(-20000,120000,n+1)','np.linspace(-40000,120000,n+1)')
base_ns=dict(thermo.init_namespace);exec(compile(base_source,__file__,'exec'),base_ns)


class InteriorBase(thermo.Flow):
    __init__=base_ns['__init__']


rhs_source=old.replace(old.rhs_source,'np.r_[self.eos.y0,y],np.r_[y,self.eos.y0]',
    'np.r_[self.incoming_y(y),y],np.r_[y,self.eos.y0]')
rhs_source=old.replace(rhs_source,'return rate,ledger,dt', '''j=self.join
    local_photon=float(np.sum(energy[j:]*m.vol[j:])-np.sum(work[j:]))
    local=np.array([C*(flux[0,j]-flux[0,-1]),C*(flux[2,j]-flux[2,-1])+local_photon,
        local_photon,C*(flux[3,j]-flux[3,-1]),np.sum(species[j:]*m.vol[j:])])
    ledger=np.r_[ledger,local]
    self.maximum_inner_characteristic=max(self.maximum_inner_characteristic,float(max(C*m.a[:j]/m.B[:j]*(abs(v[:j])+cs[:j]))))
    return rate,ledger,dt''')
ns=dict(old.ns,OUT=OUT);exec(compile(rhs_source,__file__,'exec'),ns)


class Flow(old.Flow):
    rhs=ns['rhs']
    primitive=FunctionType(old.Flow.primitive.__code__,dict(old.Flow.primitive.__globals__,OUT=OUT))

    def __init__(self,atmosphere_cells,boundary='copy',temperature='inherited'):
        assert atmosphere_cells%7==0 and boundary in ['copy','saved']
        self.atmosphere_cells=atmosphere_cells;self.join=atmosphere_cells//7;self.n=n=atmosphere_cells+self.join;self.boundary=boundary
        self.base=m=InteriorBase(n);self.eos=ex.EOS()
        rho,v,sigma=m.primitive(m.initial);theta=np.log(m.eos(rho,sigma)[3])
        inherited_temperature=theta.copy()
        if temperature=='physical':
            r=m.R+m.x/m.As;actual=np.interp(np.minimum(r,m.R),m.env['r'],np.log(m.env['T']));theta=np.where(rho>0,actual,theta)
        self.initial_temperature=theta.copy();self.temperature_model=temperature
        self.initial_temperature_change=float(max(abs(np.expm1(theta[rho>0]-inherited_temperature[rho>0]))))
        leftT=float(m.eos(np.array([m.left[0]]),np.array([m.left[2]]))[3][0]);m.left=(m.left[0],0.,np.log(leftT))
        if temperature=='physical':m.left=(m.left[0],0.,float(theta[0]))
        self.initial=self.conserved(rho,v,theta,np.full(n,self.eos.y0),m.a)[0]
        self.seed=theta.copy();self.max_recovery=0.;self.max_optical=0.;self.minimum_sigma=float(min(theta));self.maximum_sigma=float(max(theta));self.scalar_roots=0
        self.max_absorption=0.;self.max_speed=0.;self.min_y=self.eos.y0;self.maximum_inner_characteristic=0.
        self.background_cell=np.array([rho,np.zeros(n),theta])
        env=m.env;d=m.eos.d;r=m.R+m.xf/m.As
        density=np.exp(np.interp(np.minimum(r,m.R),env['r'],np.log(env['rho'])))/m.eos.rho0
        entropy=np.interp(np.minimum(r,m.R),d['initial_r'],d['initial_sigma'])
        theta=np.log(m.eos(density,entropy)[3]);vacuum_theta=float(np.log(m.eos(np.array([0.]),np.array([0.]))[3][0]))
        if temperature=='physical':theta=np.interp(np.minimum(r,m.R),env['r'],np.log(env['T']))
        self.background_left=np.array([density,np.zeros(n+1),theta]);self.background_right=self.background_left.copy()
        self.background_left[:,m.xf>0]=np.array([0.,0.,vacuum_theta])[:,None]
        self.background_right[:,m.xf>=0]=np.array([0.,0.,vacuum_theta])[:,None]
        assert m.xf[self.join]==-20000 and abs(m.dx-140000/atmosphere_cells)<1e-10

    def incoming_y(self,y):return y[0] if self.boundary=='copy' else self.eos.y0

    def reconstruct(self,V,t):
        m=self.base;delta=V[:3]-self.background_cell
        if self.boundary=='copy':ghost_delta=delta[:,0].copy()
        else:ghost_delta=np.array([0.,np.interp(t,m.hist_t,m.hist_v),0.])
        ext=np.column_stack([ghost_delta,delta,np.zeros(3)]);slope=np.zeros_like(ext)
        slope[:,1:-1]=thermo.task.minmod(ext[:,1:-1]-ext[:,:-2],ext[:,2:]-ext[:,1:-1])
        L=self.background_left+ext[:,:-1]+slope[:,:-1]/2;R=self.background_right+ext[:,1:]-slope[:,1:]/2
        if self.boundary=='copy':ghost=self.background_left[:,0]+ghost_delta
        else:ghost=np.array(m.left);ghost[1]=np.interp(t,m.hist_t,m.hist_v)
        ext=np.column_stack([ghost,V[:3],np.zeros(3)]);slope=np.zeros_like(ext)
        slope[:,1:-1]=thermo.task.minmod(ext[:,1:-1]-ext[:,:-2],ext[:,2:]-ext[:,1:-1])
        left=ext[:,:-1]+slope[:,:-1]/2;right=ext[:,1:]-slope[:,1:]/2;left[:,:2]=L[:,:2];right[:,:2]=R[:,:2]
        assert min(left[0])>=-1e-13 and min(right[0])>=-1e-13
        left[0]=np.maximum(left[0],0);right[0]=np.maximum(right[0],0)
        y=V[3];ext=np.r_[self.incoming_y(y),y,self.eos.y0];slope=np.zeros_like(ext)
        slope[1:-1]=thermo.task.minmod(ext[1:-1]-ext[:-2],ext[2:]-ext[1:-1])
        return np.vstack([left,ext[:-1]+slope[:-1]/2]),np.vstack([right,ext[1:]-slope[1:]/2])

    def run(self,label):
        assert not (OUT/f'{label}.json').exists();start=time.monotonic();m=self.base;j=self.join;U=self.initial.copy();t=0.;end=float(m.hist_t[-1]);ledger=np.zeros(10);discard=np.zeros((2,4));steps=0;next_dump=0.;history=[];local_history=[];snapshots=[]
        while t<end:
            k,l,dt=self.rhs(U,t);dt=min(dt,end-t);trial=U+dt*k
            for attempt in range(12):
                k2,l2,second_dt=self.rhs(trial,t+dt)
                if dt<=second_dt*(1+1e-12):break
                dt=min(dt/2,second_dt);trial=U+dt*k
            else:raise AssertionError('SSP stage time step')
            nxt=(U+trial+dt*k2)/2;tiny=nxt[0]<self.eos.floor
            discard[0]+=np.sum(nxt[:,tiny]*m.vol[tiny],axis=1)
            local_tiny=tiny.copy();local_tiny[:j]=False;discard[1]+=np.sum(nxt[:,local_tiny]*m.vol[local_tiny],axis=1)
            nxt[:,tiny]=0;U=nxt;ledger+=dt*(l+l2)/2;t+=dt;steps+=1
            assert steps<20000
            if t>=next_dump or t>=end:
                outside=np.sum(U[0,m.x>=0]*m.vol[m.x>=0]);energy=np.sum((U[2]-self.initial[2])*m.vol);local_energy=np.sum((U[2,j:]-self.initial[2,j:])*m.vol[j:])
                history.append([t,outside,energy,*ledger[:5],*discard[0]]);local_history.append([t,outside,local_energy,*ledger[5:],*discard[1]]);snapshots.append(U.copy());next_dump+=end/32
                np.savez_compressed(OUT/f'{label}-progress.npz',U=U,initial=self.initial,time_seconds=t,ledger=ledger,conserved_discard=discard,history=history,local_history=local_history,snapshots=snapshots,steps=steps,seed=self.seed)
        rho,v,lt,y=self.primitive(U);p,u,gamma,T,kap=self.eos(rho,lt);self.eos.y=self.eos.y0;ip,iu,*_=self.eos(self.initial[0],self.initial_temperature);self.eos.y=y
        trace=-self.eos.cx*U[0]*v*v/(1+np.sqrt(1-v*v))+rho*u-3*p-(self.initial[0]*iu-3*ip)
        entropy=(self.eos.evaluate(rho,lt)[-1]-float(self.eos.d['s0']))/self.eos.sunit;scale=4*np.pi*m.RJ**2*self.eos.rho0
        rows=[]
        for side,lo in enumerate([0,j]):
            h=np.array([history,local_history][side]);port=ledger[5*side:5*side+5];lost=discard[side];volume=m.vol[lo:]
            baryon=abs(np.sum((U[0,lo:]-self.initial[0,lo:])*volume)+lost[0]-port[0])/np.sum(self.initial[0,lo:]*volume)
            energy_residual=np.sum((U[2,lo:]-self.initial[2,lo:])*volume)+lost[2]-port[1]
            response=max(np.sum(volume*(rho[lo:]*v[lo:]**2+p[lo:])),abs(h[-1,2]),1e-100)
            energy=abs(energy_residual)/response
            species=abs(np.sum((U[3,lo:]-self.initial[3,lo:])*volume)+lost[3]-port[3]-port[4])/np.sum(self.initial[3,lo:]*volume)
            row=dict(classification='Counterexample candidate',passed=bool(baryon<1e-10 and energy<1e-8 and species<1e-9 and self.max_recovery<1e-8),
                cells=self.n-lo,actual_evolved_cells=self.n,inner_edge_cm=float(m.xf[lo]),steps=steps,seconds=time.monotonic()-start,
                baryon_ledger_relative=float(baryon),energy_ledger_over_response=float(energy),species_ledger_relative=float(species),maximum_primitive_energy_relative=self.max_recovery,
                gas_outside_original_radius_g=float(h[-1,1]*scale),integrated_trace_energy_erg=float(np.sum(trace[lo:]*volume)*scale*C*C),neutral_inventory=float(np.sum(U[3,lo:]*volume)),
                inner_baryon_into_layer_g=float(port[0]*scale),inner_Killing_nonrest_energy_into_layer_erg=float((port[1]-port[2])*scale*C*C),photon_energy_into_gas_erg=float(port[2]*scale*C*C),
                dilute_baryon_g=float(lost[0]*scale),dilute_Killing_nonrest_energy_erg=float(lost[2]*scale*C*C),
                maximum_absorption_energy_optical_depth=self.max_absorption,maximum_scattering_optical_depth=self.max_optical,
                maximum_speed_over_c=self.max_speed,maximum_inner_characteristic_cm_s=self.maximum_inner_characteristic,
                inner_characteristic_distance_cm=self.maximum_inner_characteristic*end,buffer_width_cm=20000.,minimum_logT=self.minimum_sigma,maximum_logT=self.maximum_sigma,
                finite_H_reactions=True,actual_reactive_interface=bool(side==1),far_boundary=self.boundary,temperature_model=self.temperature_model,initial_temperature_change=self.initial_temperature_change,full_GR_scalar_feedback=False,physical_chemistry_closed=False,final_charge_solved=False,full_goal_complete=False)
            tag=('full-' if side==0 else '')+label
            np.savez_compressed(OUT/f'{tag}.npz',U=U[:,lo:],initial=self.initial[:,lo:],x_cm=m.x[lo:],volume=volume,rho=rho[lo:]*self.eos.rho0,velocity_cm_s=v[lo:]*C,logT=lt[lo:],sigma=entropy[lo:],T=T[lo:],pressure=p[lo:]*self.eos.rho0*C*C,
                history=h,snapshots=np.array(snapshots)[:,:,lo:],conserved_discard=lost,initial_logT=self.initial_temperature[lo:])
            write(OUT/f'{tag}.json',row);rows.append(row);assert row['passed'],row
        print(label,json.dumps(rows[-1]),flush=True);return rows[-1]


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='b46d23f',
        claim='Replace the frozen thermal/chemical boundary at-200m by an actual interface between reactive interior and atmosphere cells using one conservative face flux.',
        decision='Test whether the original7.7739percent trace failure is removed by matching the physical equations across that interface. Retain the failed original model and gates.',
        domain_change='Explicit physical-model reassessment: add only the200m interior buffer needed to place the old interface between evolved cells. Measured atmosphere remains[-200,1200]m with exactly the original224/448/896 cell widths. Total grids256/512/1024; no automatic further extension.',
        model='Saved native rho,T and metric in the buffer, same finite H reactions and prescribed thin photon bath. Other species fixed. No artificial reaction heating correction, incoming temperature history or fitted boundary flux. Copy the far ghost perturbation; test the old frozen ghost once as a far-boundary control.',
        reuse='The entire Phase107 constitutive/rate bank and original failed flow; no new EOS table or long whole-star evolution.',
        budget=dict(native_calls=100,native_seconds=30,flow_seconds=200,CPU_threads=1,memory_GB=2,paths=['pilot224','448','896','448 saved-ghost control']),
        forecast='Phase107448 took12.42s and896 took38.86s. Added14.3percent cells at unchanged dx predicts about14-20s and45-65s; new interior recovery and boundary behavior remain unmeasured. Measure448 before fine dispatch.',
        gates=dict(native=.002,baryon=1e-10,energy=1e-8,species=1e-9,primitive=1e-8,mass_grid=.02,trace_grid=.02,neutral_grid=.02,far_boundary_trace=.001,far_boundary_interface=.001),
        stop='Do not enlarge density/temperature/composition support, grids, buffer or horizon on failure. Stop on cost/gate/support failure; a local buffer does not complete the whole stellar interior, photon transport or final charge.',
        inputs={str(p):ex.old.cold.sha(p) for p in [Path(__file__),Path(old.__file__),Path(ex.__file__),ex.OUT/'repaired-bank.npz',old.OUT/'result.json',ex.OUT/'trace-budget-decomposition.json']}))
    (OUT/'extended-background-source.py').write_text(base_source);(OUT/'extended-rhs-source.py').write_text(rhs_source)


def preflight():
    assert not (OUT/'preflight.json').exists();start=time.monotonic();signal.alarm(30);f=Flow(448);m=f.base;native=ex.Native(cap=100);checks=[]
    rho=f.initial[0];t=f.initial_temperature;p,u,g,T,kap=f.eos(rho,t)
    for i in [0,16,32,63,64,65]:
        a=native.state(float(np.log(rho[i])),float(t[i]),f.eos.y0);errors=[abs(p[i]*f.eos.rho0*C*C/a['raw'][1]-1),abs(u[i]*C*C/a['raw'][2]-1)]
        checks.append(dict(cell=i,x_cm=float(m.x[i]),rho=float(rho[i]*f.eos.rho0),T=float(T[i]),relative=list(map(float,errors))))
    assert max(max(q['relative']) for q in checks)<.002
    prior=old.Flow(448);same=bool(np.array_equal(f.initial[:,f.join:],prior.initial));assert same,'Original atmosphere initialization changed'
    assert np.array_equal(m.x[f.join:],prior.base.x) and np.array_equal(m.vol[f.join:],prior.base.vol)
    result=dict(classification='Counterexample candidate',passed=True,original_atmosphere_bitwise=True,native_controls=checks,EOS_calls=native.ion.calls,seconds=time.monotonic()-start,maximum_log_density=float(np.log(max(rho))),saved_EOS_log_density_ceiling=float(f.eos.x[-1]))
    np.savez_compressed(OUT/'native-controls.npz',**{k:np.array([z[k] for z in native.ion.states]) for k in native.ion.states[0]})
    write(OUT/'preflight.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


def run():
    assert json.loads((OUT/'preflight.json').read_text())['passed'];assert not (OUT/'result.json').exists();start=time.monotonic();signal.alarm(40)
    pilot=Flow(224).run('pilot-224');signal.alarm(max(1,int(60-(time.monotonic()-start))))
    coarse=Flow(448).run('cells-448');elapsed=time.monotonic()-start;forecast=7*coarse['seconds']+8
    write(OUT/'measured-budget.json',dict(pilot_seconds=pilot['seconds'],coarse_seconds=coarse['seconds'],elapsed_seconds=elapsed,forecast_fine_and_control_seconds=forecast,remaining_seconds=200-elapsed,assumption='Fine5.6x measured coarse plus1.4x for the same-grid far-boundary control and8s setup.'))
    assert forecast<200-elapsed,'Remaining physical-interface budget';signal.alarm(max(1,int(200-elapsed)))
    fine=Flow(896).run('cells-896');control=Flow(448,'saved').run('boundary-448')
    errors={k:abs(coarse[k]/fine[k]-1) for k in ['gas_outside_original_radius_g','integrated_trace_energy_erg','neutral_inventory']}
    trace_boundary=abs(control['integrated_trace_energy_erg']-coarse['integrated_trace_energy_erg'])/abs(coarse['integrated_trace_energy_erg'])
    # Scale the small boundary mass/energy terms by their trace consequence,
    # not by a nearly zero individual port.
    e=ex.EOS();u0=(e.u0[0]+(e.u0[1]-e.u0[0])*(e.y0-e.ys[0])/np.diff(e.ys)[0])
    boundary_interface=max(abs(2*u0*(control['inner_baryon_into_layer_g']-coarse['inner_baryon_into_layer_g'])),abs(control['inner_Killing_nonrest_energy_into_layer_erg']-coarse['inner_Killing_nonrest_energy_into_layer_erg']))/abs(coarse['integrated_trace_energy_erg'])
    passed=bool(max(errors.values())<.02 and trace_boundary<.001 and boundary_interface<.001)
    result=dict(classification='Counterexample candidate',passed=passed,refinement=errors,far_boundary_trace_relative=float(trace_boundary),far_boundary_interface_trace_scale=float(boundary_interface),seconds=time.monotonic()-start,
        actual_reactive_interface_evolved=True,original_failed_verdict_preserved=True,full_stellar_interior=False,full_photon_transport=False,full_GR_scalar_feedback=False,physical_chemistry_closed=False,final_charge_solved=False,full_goal_complete=False)
    write(OUT/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True);assert passed


PHYSICAL=OUT/'physical'


class PhysicalFlow(Flow):
    def __init__(self,n,boundary='copy'):super().__init__(n,boundary,'physical')
    run=FunctionType(Flow.run.__code__,dict(Flow.run.__globals__,OUT=PHYSICAL))
    primitive=FunctionType(Flow.primitive.__code__,dict(Flow.primitive.__globals__,OUT=PHYSICAL))


def physical_prepare():
    assert not PHYSICAL.exists();PHYSICAL.mkdir();start=time.monotonic();signal.alarm(20)
    f=PhysicalFlow(448);m=f.base;active=f.initial[0]>0;actual=np.interp(np.minimum(m.R+m.x/m.As,m.R),m.env['r'],np.log(m.env['T']))
    mismatch=float(max(abs(f.initial_temperature[active]-actual[active])));assert mismatch==0
    inherited=Flow(448);old_error=float(max(abs(np.expm1(inherited.initial_temperature[active]-actual[active]))))
    assert old_error>.002
    write(PHYSICAL/'plan.json',dict(classification='Counterexample candidate',checkpoint='b46d23f',
        prior_numeric_verdict=True,prior_physical_temperature_profile_passed=False,maximum_inherited_temperature_relative=old_error,
        root_cause='The inherited initial-entropy samples end at the old-200m boundary. Interpolating them farther inward clamps entropy, then the isentropic EOS raises T with rho. That is not the saved almost-isothermal envelope temperature.',
        correction='Initialize gas logT directly from the saved physical envelope at every occupied cell and every gas face, including the deeper ghost. Retain actual rho, geometry, EOS, rates and conservative evolution. Preserve the first numerical interface pass as a different initial-value problem.',
        scope='Actual reactive interior patch, not whole-star physical closure. All remaining radiation and chemistry premises remain explicit.',
        paths=['448 physical','896 physical','448 physical saved-ghost control'],remaining_flow_seconds=110,prior_flow_seconds=json.loads((OUT/'result.json').read_text())['seconds'],original_flow_budget=200,
        remaining_native_calls=70,remaining_native_seconds=20,unchanged_gates=json.loads((OUT/'plan.json').read_text())['gates'],
        stop='No more background, grid, duration, support or gate change on failure; assess the actual saved profile.',
        bindings={str(p):ex.old.cold.sha(p) for p in [Path(__file__),OUT/'inherited-temperature-source.py',OUT/'result.json',Path(thermo.task.old.prior.OUT)/'final-envelope.npz',ex.OUT/'repaired-bank.npz']}))
    native=ex.Native(cap=70);p,u,g,T,kap=f.eos(f.initial[0],f.initial_temperature);checks=[]
    for i in [0,16,32,63,64,65]:
        a=native.state(float(np.log(f.initial[0,i])),float(f.initial_temperature[i]),f.eos.y0)
        errors=[abs(p[i]*f.eos.rho0*C*C/a['raw'][1]-1),abs(u[i]*C*C/a['raw'][2]-1)]
        checks.append(dict(cell=i,x_cm=float(m.x[i]),T=float(T[i]),relative=list(map(float,errors))))
    assert max(max(z['relative']) for z in checks)<.002
    np.savez_compressed(PHYSICAL/'native-controls.npz',**{k:np.array([z[k] for z in native.ion.states]) for k in native.ion.states[0]})
    write(PHYSICAL/'preflight.json',dict(classification='Counterexample candidate',passed=True,saved_temperature_log_error=mismatch,maximum_inherited_temperature_relative=old_error,native_controls=checks,EOS_calls=native.ion.calls,seconds=time.monotonic()-start))
    signal.alarm(0);print((PHYSICAL/'preflight.json').read_text(),flush=True)


def physical_run():
    assert json.loads((PHYSICAL/'preflight.json').read_text())['passed'];assert not (PHYSICAL/'result.json').exists();start=time.monotonic();signal.alarm(30)
    coarse=PhysicalFlow(448).run('cells-448');elapsed=time.monotonic()-start;forecast=6.4*coarse['seconds']+5
    write(PHYSICAL/'measured-budget.json',dict(coarse_seconds=coarse['seconds'],forecast_fine_and_control_seconds=forecast,remaining_seconds=110-elapsed,
        assumption='Existing same-dx producer measured fine/coarse3.13 and control/coarse0.85; allow5x plus1.4x coarse and5s setup for corrected initial temperatures.'))
    assert forecast<110-elapsed,'Corrected temperature remaining budget';signal.alarm(max(1,int(110-elapsed)))
    fine=PhysicalFlow(896).run('cells-896');control=PhysicalFlow(448,'saved').run('boundary-448')
    errors={k:abs(coarse[k]/fine[k]-1) for k in ['gas_outside_original_radius_g','integrated_trace_energy_erg','neutral_inventory']}
    trace_boundary=abs(control['integrated_trace_energy_erg']-coarse['integrated_trace_energy_erg'])/abs(coarse['integrated_trace_energy_erg'])
    e=ex.EOS();u0=(e.u0[0]+(e.u0[1]-e.u0[0])*(e.y0-e.ys[0])/np.diff(e.ys)[0])
    interface_boundary=max(abs(2*u0*(control['inner_baryon_into_layer_g']-coarse['inner_baryon_into_layer_g'])),abs(control['inner_Killing_nonrest_energy_into_layer_erg']-coarse['inner_Killing_nonrest_energy_into_layer_erg']))/abs(coarse['integrated_trace_energy_erg'])
    passed=bool(max(errors.values())<.02 and trace_boundary<.001 and interface_boundary<.001)
    result=dict(classification='Counterexample candidate',passed=passed,refinement=errors,far_boundary_trace_relative=float(trace_boundary),far_boundary_interface_trace_scale=float(interface_boundary),seconds=time.monotonic()-start,
        saved_physical_temperature_applied=True,actual_reactive_interface_evolved=True,old_trace_failure_preserved=True,inherited_temperature_numeric_pass_preserved=True,
        full_stellar_interior=False,full_photon_transport=False,full_GR_scalar_feedback=False,physical_chemistry_closed=False,final_charge_solved=False,full_goal_complete=False)
    write(PHYSICAL/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True);assert passed


if __name__=='__main__':globals()[sys.argv[1]]()
