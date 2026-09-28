"""Native endpoint and moving-scattering checks of the actual coupled run."""
from types import SimpleNamespace
import json
import signal
import time
import numpy as np
import def_native_two_way_atmosphere as task


def main():
    out=task.OUT;assert not (out/'audit.json').exists();start=time.monotonic();signal.signal(signal.SIGALRM,task.optical.timeout);signal.alarm(20)
    verdict=json.loads((out/'result.json').read_text());rows=[v for v in verdict['paths'] if v['completed_steps']==v['steps']]
    assert rows,'No completed coupled path to audit'
    row=rows[-1];cells=row['atmosphere_cells'];steps=row['steps'];v=np.load(out/f'cells-{cells}-steps-{steps}.npz')
    flow=task.Flow(cells);m=flow.base;U=v['U'];rho,velocity,lt,y=flow.primitive(U);p,u,g,T,kap=flow.eos(rho,lt)
    tab=task.optical.Spectrum();native=tab.native;photons=v['I'].sum(0);q=photons.shape[1];mu=(np.arange(q)+.5)*2/q-1;w=np.ones(q)/q
    d=np.load(task.deep.prior.prior.OUT/'bank-16-8.npz');E=d['Einf'];number=d['num'];active=np.flatnonzero(rho>=flow.eos.floor)
    ids=sorted(set([int(active[np.argmin(lt[active])]),int(active[np.argmax(abs(velocity[active]))])]))
    checks=[]
    for j in ids:
        state=native.state(float(np.log(rho[j])),float(lt[j]),float(y[j]));D=(1-velocity[j]*mu)/np.sqrt(1-velocity[j]**2);energy=D[:,None]*E[None,:]/m.a[j]
        chi,em=task.optical.photons.coefficients(native,state,energy.ravel());ab=(chi+em).reshape(q,-1);em=em.reshape(q,-1)
        aa,ee=tab.coefficients(np.array([rho[j]*flow.eos.rho0]),np.array([lt[j]]),np.array([y[j]]),energy[None]);aa,ee=aa[0],ee[0]
        errors=[]
        for exact,estimate,field in [(ab,aa,photons[j]),(em,ee,1+photons[j])]:
            for weight in [number,number*E]:
                weights=w[:,None]*D[:,None]*weight[None,:]*field
                errors.append(float(np.sum(abs(exact-estimate)*weights)/max(float(np.sum(exact*weights)),1e-250)))
        constitutive=float(max(abs(p[j]*flow.eos.rho0*task.C**2/state['raw'][1]-1),abs(u[j]*task.C**2/state['raw'][2]-1)))
        checks.append(dict(cell=j,rho=float(rho[j]*flow.eos.rho0),T=float(T[j]),y=float(y[j]),velocity_over_c=float(velocity[j]),constitutive=constitutive,
            actual_spectrum_rate_relative=max(errors),native_population_error=state['population_error']))
    fields=v['snapshot_I'].sum(1);J=np.einsum('tiqf,q->tif',fields,w);H=np.einsum('tiqf,q->tif',fields,w*mu);K=np.einsum('tiqf,q->tif',fields,w*(mu*mu+(2/q)**2/12))
    variance=np.divide(J*K-H*H,J*J,out=np.zeros_like(J),where=J>1e-140)
    angular=bool(fields.min()>=0 and variance.min()>=-1e-12 and np.all(K<=J*(1+1e-14)))
    # Use the actual outgoing field with a controlled nonzero electron speed.
    # An elastic moving kernel must satisfy Q=v*F including frequency exits.
    probe=object.__new__(task.Coupled);probe.q=q;probe.freq=len(E);probe.mu=mu;probe.w=w;probe.number=number;probe.E=E
    probe.bulk=SimpleNamespace(edges_mu=np.linspace(-1,1,q+1));probe.scatter_number_error=0.
    j=ids[-1];vv=np.array([.001]);field=photons[j:j+1]
    collision,escape=probe.scattering(field,np.array([rho[j]*flow.eos.rho0]),vv,np.array([kap[j]]),m.a[j:j+1])
    energy=float(np.einsum('iqf,q,f->i',collision,w,number*E)[0]/m.a[j]**3+escape[0,1])
    momentum=float(np.einsum('iqf,q,f->i',collision,w*mu,number*E)[0]/m.a[j]**3+escape[0,2])
    elastic=abs(energy-vv[0]*momentum)/max(abs(energy),abs(vv[0]*momentum),1e-250)
    result=dict(classification='Counterexample candidate',passed=bool(angular and elastic<1e-8 and probe.scatter_number_error<1e-12 and max(max(z['constitutive'],z['actual_spectrum_rate_relative']) for z in checks)<.002),
        audited_path=f'cells-{cells}-steps-{steps}',native_checks=checks,minimum_photon_occupation=float(fields.min()),minimum_variance=float(variance.min()),
        elastic_moving_energy_momentum_relative=float(elastic),scattering_number_relative=probe.scatter_number_error,native_calls=native.ion.calls,
        seconds=time.monotonic()-start,source_sha256=task.sha(__file__),full_coupled_comparison_passed=verdict['passed'],final_charge_solved=False)
    task.write(out/'audit.json',result);signal.alarm(0);print(json.dumps(result),flush=True);assert result['passed']


def cold_root():
    """Same-equation density-start repair on the exact failed cold cell."""
    import sys
    from scipy.optimize import brentq
    out=task.OUT;assert not (out/'cold-density-repair.json').exists()
    task.write(out/'cold-density-repair-plan.json',dict(classification='Counterexample candidate',
        failure='Native status125 comes from eos_calc trial mass-density underflow. It occurs inside the density-variable initializer for a finite requested density, before the constraint solver can return the desired cold state.',
        repair='Route status125 through the existing electron-variable solve used for status4. Solve exactly the same density equation with unchanged monotonicity and residual gates; do not substitute an ideal EOS or change any table support.',
        decision='A native root below240K confirms the missing constitutive support. Only after that repair may a separate bounded support-and-resume plan be considered. No fluid replay in this control.',
        additional_native_call_cap=40,wall_seconds=20,previous_probe_call_upper_bound=16,
        reason_for_reassessment='The original cold probe terminated at160K before recording its counter; retain its16call upper bound. Register this bounded initializer repair separately from the original1100call phase budget.',
        source_sha256=task.sha(__file__),ion_source_sha256=task.sha(task.interface.ex.old.cold.old.__file__)))
    signal.signal(signal.SIGALRM,task.optical.timeout);signal.alarm(20);start=time.monotonic()
    f=task.Flow(896);m=f.base;d=np.load(out/'root-failure-state.npz');i=416;U=d['U'][:,i];D=U[0]
    tau=(U[2]-(m.a[i]-m.a0)*f.eos.cx*D)/m.a[i];y=U[3]/D;native=task.optical.ex.Native(cap=40);rows=[];failure=None;rootT=None
    try:
        def residual(T):
            # Four pressure/momentum iterations match the actual scalar inverse.
            pp=0.
            for k in range(4):
                v=U[1]/(f.eos.cx*D+tau+pp);r=np.sqrt(1-v*v);rho=D*r;W=1/r;wm=v*v/(r*(1+r))
                state=native.state(float(np.log(rho)),float(np.log(T)),float(y));nextp=state['raw'][1]/(f.eos.rho0*task.C**2)
                if k and abs(nextp-pp)<max(abs(nextp),1e-100)*1e-12:pp=nextp;break
                pp=nextp
            u=state['raw'][2]/task.C**2;res=(f.eos.cx*D*wm+(rho*u+pp)*W*W-pp-tau)/tau
            rows.append(dict(T=float(T),relative_energy_residual=float(res),population_error=state['population_error']))
            return float(res)
        a=residual(240.);b=residual(160.);assert a*b<0,('Native positive cold bracket',a,b)
        rootT=float(brentq(residual,160.,240.,xtol=1e-5));final=residual(rootT)
        assert abs(final)<2e-11
    except Exception as exc:failure=repr(exc)
    result=dict(classification='Counterexample candidate',passed=failure is None,failure=failure,cold_cell=i,
        native_root_temperature_K=rootT,checks=rows,native_calls=native.ion.calls,
        density_fallbacks=native.ion.density_fallbacks,seconds=time.monotonic()-start,
        atmosphere_table_support_unchanged=True,full_fine_path_completed=False)
    task.write(out/'cold-density-repair.json',result);signal.alarm(0);print(json.dumps(result),flush=True);assert result['passed']


if __name__=='__main__':
    import sys
    globals()[sys.argv[1] if len(sys.argv)>1 else 'main']()
