"""Apply the existing heat input to the complete canonical weak GR system."""
from pathlib import Path
import argparse
import json
import resource
import signal
import time
import numpy as np
import def_gr_energy_fem as fem
import def_conduction_radau as radau

OUT=fem.OUT/'evolution';write=fem.write
digest=fem.base.prior.digest
FIELDS=fem.base.FIELDS
BANK=fem.base.prior.OUT


def readout(model,u,v,t):
    heat=model.heat;p=model.original;r=model.grid
    flux,energy=heat.faces(t);q=u-model.H@energy;qt=v-model.H@(flux*p.radiation.geometry.tc)
    nf=model.surface_index+1;z=q[model.indices[:nf,0]];vz=qt[model.indices[:nf,0]]
    psi=np.r_[q[model.indices[:-1,1]],0.]
    N,a=p.radiation.geometry.metric(r[:nf])
    cv=np.interp(p.native,r[:nf],a/N*r[:nf]*vz*fem.base.task.h.gr.C)
    cf=np.interp(p.native,r,psi)
    row=dict(tau=float(t),velocity_mass_RMS_m_s=float(np.sqrt(p.weights@(cv*cv))),
        scalar_mass_RMS=float(np.sqrt(p.weights@(cf*cf))),
        **{name:float(np.sqrt(p.weights[m]@cv[m]**2/p.weights[m].sum())) for name,m in zip(FIELDS[2:],p.masks)})
    return row,q,qt,cv,cf


def solve(model,steps,label):
    start=time.monotonic();dt=1/steps;step=radau.Step(model.K,-model.M,dt)
    u=np.zeros(model.size);v=u.copy();history=[];native_v=[];native_psi=[]
    forcing=lambda t:model.load@model.heat.faces(t)[1]
    step_start=time.monotonic()
    for j in range(steps+1):
        row,q,qt,cv,cf=readout(model,u,v,j*dt);history.append(row);native_v.append(cv);native_psi.append(cf)
        if j<steps:u,v=step.advance(u,v,j*dt,forcing)
    flux,energy=model.heat.faces(1.)
    balance=float(abs(np.sum(-np.diff(energy),dtype=np.longdouble))/max(abs(energy).max(),1e-100))
    result=dict(classification='Counterexample candidate',steps=steps,outer=model.outer,coarse_grid=model.coarse,
        linear_residual=step.error,heat_telescoping=balance,seconds=time.monotonic()-start,
        step_seconds=time.monotonic()-step_start,history=history)
    np.savez_compressed(OUT/(label+'.npz'),grid=model.grid,indices=model.indices,q=q,qt=qt,u=u,ut=v,
        native_radius=model.original.native,native_velocity=np.array(native_v),native_scalar=np.array(native_psi),
        weights=model.original.weights,masks=np.array(model.original.masks),heat_energy=energy,heat_flux=flux)
    write(OUT/(label+'.json'),result);print('EVOLUTION',label,result['seconds'],row,flush=True)
    return result


def prepare():
    assert not OUT.exists();OUT.mkdir();signal.alarm(60)
    paths=[Path(__file__),Path(fem.__file__),Path(fem.canonical.__file__),Path(radau.__file__),
        BANK/'fine-bank.npz',BANK/'coarse-bank.npz',fem.OUT/'symbolic.json',fem.OUT/'assembly.json']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='29fe9acd',
        claim='Apply the same 4012-face heat input and exact microscopic startup to actual fluid/scalar evolution from the complete continuous canonical GR weak form; test the previously failed component histories.',
        changed='Full consistent weak spatial stiffness, inertia, physical pressure and source. Same original nodes, horizon, frozen background and heat input. No partial mass correction and no modal clipping.',
        decision='Run the original16/32/64 time counts. Stop if any of four components fails. Only if all pass, run fixed coefficient, outer3R and every-other-node spatial contrasts. Never refine on failure.',
        gates=dict(time_relative=.02,time_order=1.5,coefficient=.02,outer=.002,spatial_contrast=.02,linear_residual=1e-9,heat_balance=2e-13),
        budget=dict(pilot_steps=8,production_steps=304,maximum_paths=6,hard_seconds=240,CPU_threads=1,memory_GB=3,new_EOS_calls=0,automatic_expansion=False),
        scope='A candidate discretization of the original continuous linear GR equations. Changing the spatial scheme requires its own validation. No claim of full nonlinear, photon, outer-current, charge or observational closure.',
        bindings={str(p):digest(p) for p in paths}))
    write(OUT/'control.json',radau.control())
    start=time.monotonic();model=fem.Model(BANK/'fine-bank.npz');setup=time.monotonic()-start
    pilot=solve(model,8,'pilot');forecast=1.5*(4*setup+304/8*pilot['seconds'])
    write(OUT/'pilot-budget.json',dict(classification='Counterexample candidate',setup_seconds=setup,
        solve_seconds=pilot['seconds'],forecast_seconds=forecast,assumption='Four setups and304 steps scaled from8 measured steps with50 percent margin. Coefficient, outer and coarser-grid cost unmeasured.'))
    print('FORECAST',forecast,flush=True);signal.alarm(0)


def run():
    assert not (OUT/'result.json').exists();plan=json.loads((OUT/'plan.json').read_text())
    for p,s in plan['bindings'].items():assert digest(Path(p))==s,p
    assert json.loads((OUT/'pilot-budget.json').read_text())['forecast_seconds']<240
    signal.alarm(240);resource.setrlimit(resource.RLIMIT_AS,(int(3e9),int(3e9)));start=time.monotonic()
    model=fem.Model(BANK/'fine-bank.npz');cases={str(n):solve(model,n,'fine-'+str(n)) for n in [16,32,64]}
    comparisons={}
    for field in FIELDS:
        a,b,c=[np.array([row[field] for row in cases[str(n)]['history']]) for n in [16,32,64]]
        norm=max(abs(c).max(),1e-100);d1=np.max(abs(a-b[::2]))/norm;d2=np.max(abs(b-c[::2]))/norm
        comparisons[field]=dict(time_previous=float(d1),time_last=float(d2),order=float(np.log2(d1/d2)))
    time_passed=all(row['time_last']<.02 and row['order']>1.5 for row in comparisons.values())
    if time_passed:
        for label,bank,outer,coarse in [('coefficient','coarse-bank.npz',2,False),('outer','fine-bank.npz',3,False),('spatial','fine-bank.npz',2,True)]:
            other=fem.Model(BANK/bank,outer,coarse);cases[label]=solve(other,64,label+'-64')
            for field in FIELDS:
                c=np.array([row[field] for row in cases['64']['history']]);d=np.array([row[field] for row in cases[label]['history']])
                comparisons[field][label]=float(np.max(abs(c-d))/max(abs(c).max(),1e-100))
    passed=time_passed and all(row.get('coefficient',1)<.02 and row.get('outer',1)<.002 and row.get('spatial',1)<.02 for row in comparisons.values())
    passed=passed and all(row['linear_residual']<1e-9 and row['heat_telescoping']<2e-13 for row in cases.values())
    result=dict(classification='Counterexample candidate',actual_GR_fluid_scalar_evolved=True,
        time_passed=time_passed,passed=passed,comparisons=comparisons,endpoint=cases['64']['history'][-1],
        paths=list(cases),seconds=time.monotonic()-start,memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
        original_failure_resolved=False,physical_pressure_metric_reconstruction_checked=False,
        full_dynamic_charge_solved=False,photon_transport=False,nonlinear_feedback=False,
        limitation='Even a numerical pass still requires full pressure/metric reconstruction and spatial consistency audit before declaring the earlier failure repaired. The coarse coefficient bank previously failed its own local tau gate.')
    write(OUT/'result.json',result);signal.alarm(0);print('RESULT',json.dumps(result),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','run']);globals()[parser.parse_args().action]()
