"""Read the updated actual coupled paths into the same retarded GR charge."""
from pathlib import Path
from types import FunctionType
import json,signal,time,shutil
import numpy as np
import def_native_refined_thermochemistry as run
import verify_native_stage_energy_charge as old

OUT=run.GR;flow=run.flow;write=flow.write;sha=flow.sha;C=flow.C


def endpoint():
    """Independent same-native endpoint after the new actual evolution."""
    model=run.Coupled();b=model.bulk;e=b.eos;z=np.load(run.EV/'checkpoint-128.npz');assert int(z['completed'])==128
    model.Pi=z['Pi'];model.h=z['h'];model.j=z['j'];model.mass=model.mass0-model.h[1:]+model.h[:-1];model.set_material(flow.old.END)
    theta=z['theta'];eta=z['eta'];p,u,*_=e.gas(theta,eta);ab,em,*_=e.radiation(theta,eta)
    native=run.chem.old.Native(cap=200);rows=[];new=[];base=e.base.radiation(theta+e.theta0,eta)
    df=e.spectral['frequency'][0]*(1-e.f)+e.spectral['frequency'][1]*e.f
    for j in range(b.n):
        run.chem.setup(native,e.base.d,j)
        state=native.state(float(np.log1p(e.density_shift[j])+np.log1p(e.x[j])),
            float(np.log(e.base.d['T'][j])+e.theta0[j]+theta[j]),float(e.base.d['y0'][j]*(1+eta[j])))
        raw=state['raw'];pp=raw[1]-e.inventory[0,j]*e.xi[j];uu=raw[2]-e.inventory[1,j]*e.xi[j]
        chi,ee=run.chem.prior.coefficients(native,state,b.d['Einf']/b.d['a'][j]);coeff=np.array([chi+ee,ee])
        for k in [0,1]:coeff[k]-=base[k][j]*e.rinventory[k,j]*e.xi[j]*(1+df[j,k]*e.frequency_shift[j])
        new.append(coeff)
        rows.append(dict(cell=j,pressure_relative=float(abs(p[j]/pp-1)),energy_relative=float(abs(z['u'][j]/uu-1)),
            population_error=float(state['population_error']),raw=raw.tolist(),fraction=state['fraction'].tolist(),affinity=float(state['affinity'])))
    nr=np.array(new);I=z['xb']*b.scale;fac=b.factor[:,None,None];number=b.volume[:,None,None]*b.w[None,:,None]*b.d['num']/b.d['a'][:,None,None]**3
    # Opposite H transfer is paired with the exact same bound-free source.
    beta=model.velocity();boost=-beta[:,None,None]*b.mu[None,:,None]
    spectral=e.spectral['frequency'][0]*(1-e.f)+e.spectral['frequency'][1]*e.f
    def bound(a,ee):
        return fac*(ee[:,None,:]-(a-ee)[:,None,:]*I)+fac*boost*((ee*(1+spectral[:,1]))[:,None,:]*(1+I)-(a*(1+spectral[:,0]))[:,None,:]*I)
    current=bound(ab,em);delta=bound(nr[:,0]-ab,nr[:,1]-em)
    Hdef=float(np.sum(abs(np.sum(delta*number,axis=(1,2))))/max(np.sum(abs(np.sum(current*number,axis=(1,2)))),1.))
    errors=[]
    for k,field in [(0,I),(1,np.ones_like(I)),(1,I)]:
        weight=number*b.d['Einf']*fac*field
        errors.append(float(np.sum(abs(nr[:,k]-[ab,em][k])[:,None,:]*weight)/max(np.sum(nr[:,k,None,:]*weight),1.)))
    thermo=max(max(r['pressure_relative'],r['energy_relative']) for r in rows)
    return dict(classification='Counterexample candidate',passed=thermo<.002 and max(errors)<.002 and Hdef<.02,
        maximum_constitutive_relative=thermo,positive_rate_relative=errors,net_neutral_exchange_relative=Hdef,
        native_calls=native.ion.calls,rows=rows,uniform_EOS_error_bound=False,scope='Actual new fine endpoint only; earlier323-state controls are on the previous saved path.')


def conserved_inventory():
    """Check the inherited baryon gate directly from both terminal states."""
    model=run.Coupled();f=model.flow;m=model.m;rows=[]
    for steps in [64,128]:
        z=np.load(run.EV/f'checkpoint-{steps}.npz');assert int(z['completed'])==steps
        baryon=abs(np.sum(-np.diff(z['h']),dtype=np.longdouble)+
            (np.sum((z['U'][0]-f.initial[0])*m.vol,dtype=np.longdouble)+z['discard'][0])*model.gas_scale/C**2)
        initial=np.sum(f.initial[0]*m.vol)*model.gas_scale/C**2
        rows.append(dict(steps=steps,baryon_relative_to_initial_atmosphere=float(baryon/initial),
            baryon_relative_to_actual_port=float(baryon/max(abs(z['ledger'][0]*model.gas_scale/C**2),1.))))
    return dict(classification='Counterexample candidate',passed=max(r['baryon_relative_to_initial_atmosphere'] for r in rows)<1e-10,
        rows=rows,nonlinear_stage_residual_upper_from_unchanged_owner=2e-11,
        scope='Each successful inherited implicit stage exits only below2e-11 or raises; its observed maximum was not separately stored. Whole photon/neutral port history was not captured and is not certified here.')


def main():
    assert not (OUT/'result.json').exists();assert json.loads((run.EV/'result.json').read_text())['passed']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        claim='Apply actual new native-thermal coupled matter/photon sources and accepted stage ports to retarded charge; recompute its source-dependent conditional bounds and independently check the new native endpoint.',
        reuse='Same corrected initial GR operator, same64/128 clocks and531 cells. No new evolution. Reuse coefficient geometry only; no old source norm or interval is inherited.',
        gates=dict(time=.02,quadrature=.002,independent=1e-9,energy=1e-8,baryon=1e-10,native=.002,net_H=.02),
        budget=dict(seconds=60,native_calls=200,CPU_threads=1,memory_GB=3),
        limits='Current initial first-variation GR operator, restricted finite H/Thomson model and declared floor condition. Initial thermodynamic coefficients in the GR operator retain their prior representation; not nonlinear GR or a certified full EOS/transport fixed point.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(run.__file__),Path(old.__file__),run.EV/'source-64.npz',run.EV/'source-128.npz',run.EV/'result.json',run.OUT/'bank.npz',run.OUT/'controls/result.json']}))
    start=time.monotonic();signal.signal(signal.SIGALRM,flow.old.optical.timeout);signal.alarm(60)
    native=endpoint();write(OUT/'native-endpoint.json',native);inventory=conserved_inventory();write(OUT/'conserved-inventory.json',inventory)
    # Preserve successful terminal state as an artifact, not just a rolling
    # untracked checkpoint. This avoids replay for the next native audit.
    for steps in [64,128]:shutil.copyfile(run.EV/f'checkpoint-{steps}.npz',run.EV/f'final-{steps}.npz')
    m=old.gr.Response();paths={}
    for steps,order in [(128,8),(64,8),(128,4)]:
        d=dict(np.load(run.EV/f'source-{steps}.npz'));wave=old.read(m,d,order);paths[steps,order]=wave
        np.savez_compressed(OUT/f'wave-{steps}-g{order}.npz',**wave)
    d=dict(np.load(run.EV/'source-128.npz'));fine=paths[128,8];norm=max(np.max(abs(fine['free_scalar'])),1e-300)
    time_error=float(max(abs(fine['free_scalar'][::2]-paths[64,8]['free_scalar']))/norm)
    quadrature=float(max(abs(fine['free_scalar']-paths[128,4]['free_scalar']))/norm)
    direct,residual=old.independent.direct(m,d,8);agreement=abs(direct/fine['direct_scalar'][-1]-1)
    enclose=FunctionType(old.enclosure.__code__,dict(vars(old),OUT=OUT))
    bound=enclose(m,d,np.load(run.EV/'history-128.npz'));end=float(fine['free_scalar'][-1]);error=bound['total_conditional_normalized_bound']
    interval=[float(np.nextafter(end-error,-np.inf)),float(np.nextafter(end+error,np.inf))]
    previous=np.load(old.OUT/'wave-128-g8.npz');delta=fine['free_scalar']-previous['free_scalar']
    old_end=float(previous['free_scalar'][-1]);deeper=json.loads((run.failed.prior.OUT/'result.json').read_text())
    result=dict(classification='Counterexample candidate',passed=bool(native['passed'] and inventory['passed'] and time_error<.02 and quadrature<.002 and agreement<1e-9 and bound['remaining_source_port_relative']<1e-8 and interval[0]>0),
        actual_refined_thermal_EOS_evolution_completed=True,actual_updated_sources_applied_to_GR=True,
        endpoint_free_scalar=end,previous_endpoint_free_scalar=old_end,endpoint_relative_change=end/old_end-1,
        waveform_change_relative=float(max(abs(delta))/norm),pressure_only_prediction=deeper['native_pressure_endpoint'],
        time_relative=time_error,quadrature_relative=quadrature,independent_direct_relative=agreement,
        retarded_inverse_residual_cm=residual,conditional_scalar_interval=interval,bound=bound,
        native_endpoint_passed=bool(native['passed']),conserved_inventory=inventory,seconds=time.monotonic()-start,
        full_source_error_enclosed=False,atmosphere_native_inverse_audited=False,full_floor_feedback_enclosed=False,
        uniform_EOS_derivative_bound=False,coupled_fixed_point_verified=False,nonlinear_GR=False,final_charge_solved=False,full_goal_complete=False)
    write(OUT/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
