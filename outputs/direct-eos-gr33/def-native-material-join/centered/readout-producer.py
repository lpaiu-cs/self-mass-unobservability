"""Counterexample candidate: remove the diagnosed coarse acoustic mass diffusion.

Same physical shared HLL face, native EOS, photon equations and finite volumes.
The interior central acoustic flux is restricted to this smooth short horizon;
it is not a shock method or an unconditional RK2 stability theorem.
"""
from pathlib import Path
from types import FunctionType,MethodType
import inspect
import json
import signal
import sys
import textwrap
import time
import numpy as np
import sympy as sp
import def_native_material_join as base
import verify_native_material_join as verify

old=base.old;previous=base.previous;C=base.C;write=base.write;sha=base.sha
OUT=base.OUT/'centered';replace=base.replace
method=textwrap.dedent(inspect.getsource(base.Coupled.material_rhs))
method=replace(method,'+(pl-pr)/(2*Z)','')
method=replace(method,'+Z*(vl-vr)/2','')
ns=dict(vars(base),OUT=OUT);exec(compile(method,__file__,'exec'),ns)


class Coupled(base.Coupled):
    material_rhs=ns['material_rhs']
    run=FunctionType(base.Coupled.run.__code__,dict(base.Coupled.run.__globals__,OUT=OUT),argdefs=base.Coupled.run.__defaults__)


read_ns=dict(verify.namespace,OUT=OUT,run=sys.modules[__name__])
for name in ['source','readout','audit']:
    fn=verify.namespace[name];read_ns[name]=FunctionType(fn.__code__,read_ns,argdefs=fn.__defaults__)


def prepare():
    assert not OUT.exists();OUT.mkdir()
    diagnostic=json.loads((base.OUT/'diagnostic.json').read_text());assert diagnostic['diffusive_rest_over_direct']>.99
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='b024c099c',
        failure='Shared face evolution completed and passed time/energy/baryon checks, but the new coarse upwind acoustic pressure-jump mass flux alone reconstructs100.147percent of the large negative direct charge. The original upwind trajectory, strict neutral-sum failure and producer remain intact.',
        repair='Remove only the artificial interior pressure-jump mass flux and velocity-jump momentum dissipation. Use central acoustic face velocity/pressure on the same native16 volumes. Retain the actual nonlinear shared-face HLL flux, simultaneous SSP material coupling, native thermochemistry and photons. No higher grid or longer history.',
        decision='Does the same free material interface yield a finite direct residual after the diagnosed artificial baryon transport is removed? Keep physical sign unresolved if the original2percent time/readout gates fail.',
        stability='Central acoustic semidiscretization has a skew-adjoint constant-coefficient flat limit. SSP RK2 is not unconditionally energy stable; measure its finite-horizon flat modal amplification on the actual deep acoustic scale, alongside nonlinear support/conservation pilots. No shock or long-time claim.',
        resource_reassessment='The original450s dispatch completed in about325s. A new numerical equation is required; do not silently reuse the remaining budget or rerun the failed upwind equation. Add at most400s for exactly the two corrected64/128 paths, after their two-step pilot forecasts pass1.6x remainder plus10s. Reuse all native banks. Combined production caps850s; actual earlier cost retained.',
        budget=dict(pilot_seconds=30,production_seconds=400,saved_readout_seconds=45,CPU_threads=1,memory_GB=3,new_native_bank_calls=0,paths=[64,128]),
        gates=dict(energy=1e-8,baryon=1e-10,source=1e-9,time_direct=.02,time_total=.02,time_trace=.02,native=.002,flat_modal_growth=.000001),
        stop='Stop on original positivity/support/conservation/time gates, measured forecast or hard alarm. No automatic third path, refinement or gate relaxation.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(base.__file__),Path(verify.__file__),base.OUT/'diagnostic.json',base.OUT/'production.json',base.OUT/'result.json']}))
    (OUT/'expanded-material-rhs.py').write_text(method)
    plan=json.loads((base.OUT/'readout-plan.json').read_text());plan.update(
        claim='Actual centered-interior/shared-HLL material source and retarded charge; preserve the rejected upwind result separately.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(base.__file__),Path(verify.__file__)]})
    write(OUT/'readout-plan.json',plan);(OUT/'readout-producer.py').write_bytes(Path(__file__).read_bytes())


def pilot():
    assert not (OUT/'pilot.json').exists();start=time.monotonic();signal.signal(signal.SIGALRM,old.optical.timeout);signal.alarm(30)
    m=Coupled();sound=np.sqrt(m.face_K/(m.cx*m.face_rho));rate=float(max(sound)/min(np.diff(m.edge)))
    # Independent periodic flat acoustic operator and its energy weight.
    D=np.array([[0,1,-1],[-1,0,1],[1,-1,0]],float)/2;L=np.block([[np.zeros((3,3)),-D],[-D,np.zeros((3,3))]])
    assert np.array_equal(L+L.T,np.zeros((6,6)))
    h=old.END/64;growth=float(np.expm1(4*64*np.log1p((2*rate*h)**4/4)/2));assert growth<1e-6
    x=sp.symbols('x',real=True);assert sp.expand((1-x*x/2)**2+x*x-1-x**4/4)==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='Periodic flat constant-coefficient central acoustic operator is skew-adjoint. RK2 squared modal amplification is1+x^4/4, not exactly1. No nonlinear or stratified unconditional stability theorem.'))
    rows=[]
    for steps in [64,128]:
        row=(m if steps==64 else Coupled()).run(steps,f'pilot-{steps}',2);rows.append(row)
        if not row['passed']:break
    forecast=None if len(rows)!=2 or not all(r['passed'] for r in rows) else sum(r['seconds']/2*(r['steps']-2) for r in rows)
    upper=None if forecast is None else 1.6*forecast+10
    write(OUT/'pilot.json',dict(classification='Counterexample candidate',paths=rows,forecast_seconds=forecast,upper_seconds=upper,
        eligible=bool(upper is not None and upper<400),flat_finite_horizon_modal_growth= growth,seconds=time.monotonic()-start));signal.alarm(0)


def production():
    assert json.loads((OUT/'pilot.json').read_text())['eligible'];assert not (OUT/'production.json').exists()
    for path,value in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert sha(path)==value,path
    start=time.monotonic();signal.signal(signal.SIGALRM,old.optical.timeout);signal.alarm(400);rows=[]
    for steps in [64,128]:
        row=Coupled().run(steps,f'coupled-{steps}',restart=f'pilot-{steps}');rows.append(row)
        if not row['passed']:break
    write(OUT/'production.json',dict(classification='Counterexample candidate',passed=bool(len(rows)==2 and all(r['passed'] for r in rows)),paths=rows,seconds=time.monotonic()-start,
        shared_material_face=True,coarse_interior_acoustic_dissipation_removed=True,full_GR_feedback=False,final_charge_solved=False));signal.alarm(0)


def readout():read_ns['readout']()


def audit():
    read_ns['audit']()
    start=time.monotonic();m=Coupled();z=np.load(OUT/'coupled-128.npz');m.h=z['h'];m.Pi=z['Pi'];m.mass=m.mass0-m.h[1:]+m.h[:-1];m.set_material(old.END)
    m.recover_material(m.material_state(z['u'],z['eta']),old.END,z['theta']);f=m.flow;V=f.primitive(z['U']);_,R=f.reconstruct(V,old.END);flux=m.join_flux(R[:,0],V[3,0]);port=m.mflux.copy()
    # Independently retain all endpoint fluxes and difference in extended
    # precision. Do not subtract huge rounded divergences to infer a tiny port.
    src=replace(method,'return [mass,force,-np.diff(energy),-np.diff(neutral)]','return mass,force,energy,neutral')
    scope=dict(vars(base));exec(compile(src,__file__,'exec'),scope)
    mass,force,energy,neutral=MethodType(scope['material_rhs'],m)(z['u'],z['theta'],z['eta']);rows=[]
    for key,F in [('mass',mass),('energy',energy),('neutral',neutral)]:
        LD=np.longdouble;d=-np.diff(F.astype(LD));relative=float(abs(d.sum(dtype=LD)+LD(F[-1])-LD(F[0]))/max(abs(F[-1]),1.))
        rows.append(dict(quantity=key,extended_precision_port_relative=relative))
    converted=flux*np.array([C,C*C,C**3,C*f.eos.nH])*4*np.pi*m.m.RJ**2*f.eos.rho0
    conversion=float(max(abs(converted-port)/np.maximum(abs(port),1.)))
    passed=max(conversion,max(r['extended_precision_port_relative'] for r in rows))<1e-8
    write(OUT/'join-audit.json',dict(classification='Counterexample candidate',passed=bool(passed),endpoint_flux_checks=rows,
        conversion_relative=conversion,integrated_shared_mass_g=float(z['scalar_join_mass']),integrated_shared_momentum_g_cm_s=float(z['scalar_join_momentum']),
        integrated_shared_Killing_energy_erg=float(z['scalar_join_energy']),integrated_shared_neutral=float(z['scalar_join_neutral']),
        whole_neutral_reaction_trajectory_audited=False,physical_spatial_reconstruction_certified=False,seconds=time.monotonic()-start))
    assert passed


if __name__=='__main__':globals()[sys.argv[1]]()
