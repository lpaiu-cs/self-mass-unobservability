"""Counterexample candidate: local radial comparison of actual coupled paths."""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import inspect
import json
import signal
import sys
import time
import numpy as np
import sympy as sp
import def_native_boundary_layer as run
import verify_native_material_join as inherited

OUT=run.OUT;write=run.write;sha=run.sha
runtime=SimpleNamespace(**dict(vars(run),previous=run.feedback))
ns=dict(run.prior.read_ns,OUT=OUT,run=runtime)
for name in ['source','readout','audit']:
    f=run.prior.read_ns[name];ns[name]=FunctionType(f.__code__,ns,argdefs=f.__defaults__)


def prepare():
    assert not (OUT/'readout-plan.json').exists()
    plan=json.loads((run.prior.OUT/'readout-plan.json').read_text())
    plan.update(claim='Export actual19-cell material/photon sources; use exactly the same observer clock and retarded operator, then compare with the completed16-cell reference.',
        limits='The last-cell local radial comparison is empirical. The other fifteen cells, spectral/angular discretization, first-order inventory/density and initial/dynamic metric closure remain unverified globally.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(run.__file__),Path(inherited.__file__),run.prior.OUT/'source-128.npz',run.prior.OUT/'wave-128.npz']})
    write(OUT/'readout-plan.json',plan);(OUT/'readout-producer.py').write_bytes(Path(__file__).read_bytes())
    # Reused analytical identities apply at every internal radial partition.
    F,h=sp.symbols('F h');assert sp.expand(-h*F+h*F)==0
    v=sp.symbols('v0:5');assert sp.expand(sum(v[i]-v[i+1] for i in range(4))-(v[0]-v[-1]))==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        scope='A common face flux cancels when adjacent native cells are summed. Refining a cell adds internal cancelling faces; the two external ports remain. This identity is not a spatial-error theorem.'))


def readout():ns['readout']()


def radial():
    assert not (OUT/'radial.json').exists();start=time.monotonic();signal.signal(signal.SIGALRM,run.old.optical.timeout);signal.alarm(15)
    fine=np.load(OUT/'wave-128.npz');coarse=np.load(run.prior.OUT/'wave-128.npz');assert np.array_equal(fine['t'],coarse['t'])
    sf=np.load(OUT/'source-128.npz');sc=np.load(run.prior.OUT/'source-128.npz');assert np.array_equal(sf['t'],sc['t'])
    def relative(a,b):return float(np.max(abs(a-b))/max(float(np.max(abs(a))),1e-300))
    rows=dict(space_direct=relative(fine['direct_relative'],coarse['direct_relative']),
        space_total=relative(fine['direct_plus_photon_mass'],coarse['direct_plus_photon_mass']),
        space_surface_spectrum=relative(sf['luminosity_per_mu'],sc['luminosity_per_mu']))
    result=json.loads((OUT/'result.json').read_text());passed=result['passed'] and all(v<.02 for v in rows.values())
    z=np.load(OUT/'coupled-128.npz');zc=np.load(run.prior.OUT/'coupled-128.npz')
    ports={key:dict(refined=float(z['scalar_'+key]),original=float(zc['scalar_'+key]),relative_to_refined=float(abs(z['scalar_'+key]-zc['scalar_'+key])/max(abs(z['scalar_'+key]),1.))) for key in ['join_mass','join_momentum','join_energy','join_neutral']}
    row=dict(classification='Counterexample candidate',passed=bool(passed),controls=rows,shared_material_ports=ports,
        original_cells=16,refined_cells=19,nearest_center_face_gap_m=300.,observer_clock_bitwise_identical=True,
        endpoint_direct_original=float(coarse['direct_relative'][-1]),endpoint_direct_refined=float(fine['direct_relative'][-1]),
        endpoint_total_original=float(coarse['direct_plus_photon_mass'][-1]),endpoint_total_refined=float(fine['direct_plus_photon_mass'][-1]),
        local_last_cell_radial_comparison_only=True,global_radial_error_certified=False,space_time_error_interaction_certified=False,
        seconds=time.monotonic()-start,new_fluid_steps=0,final_charge_solved=False)
    write(OUT/'radial.json',row);signal.alarm(0);print(json.dumps(row),flush=True)


def audit():
    ns['audit']()
    old=run.prior.previous;original=np.load(old.chem.prior.OUT/'bank-16-8.npz');new=np.load(OUT/'geometry.npz')
    for key in ['r','a','B','rho','T','phi','raw','y0','thermo','target']:assert np.array_equal(original[key][:15],new[key][:15]),key
    for oldpath,newpath,keys,axis in [
        (run.angular.OUT/'thermal-support/bank.npz',OUT/'thermal-support/bank.npz',['raw','rates'],0),
        (old.motion.OUT/'native.npz',OUT/'native.npz',['raw','pr','ur','pt','ut','K'],1),
        (old.OUT/'bank.npz',OUT/'bank.npz',['density','frequency'],1)]:
        a=np.load(oldpath);b=np.load(newpath)
        for key in keys:
            aa,bb=(a[key][:15],b[key][:15]) if axis==0 else (a[key][:,:15],b[key][:,:15])
            assert np.array_equal(aa,bb),(str(oldpath),key)
    m=run.Coupled();d=m.bulk.d;z=np.load(OUT/'coupled-128.npz');m.h=z['h'];m.Pi=z['Pi'];m.mass=m.mass0-m.h[1:]+m.h[:-1];m.set_material(run.old.END)
    m.recover_material(m.material_state(z['u'],z['eta']),run.old.END,z['theta']);f=m.flow;V=f.primitive(z['U']);_,R=f.reconstruct(V,run.old.END);flux=m.join_flux(R[:,0],V[3,0])
    scale=4*np.pi*m.m.RJ**2*f.eos.rho0
    check=flux*np.array([run.C,run.C**2,run.C**3,run.C*f.eos.nH])*scale
    error=float(np.max(abs(check-m.mflux)/np.maximum(abs(m.mflux),1.)));assert error<1e-12
    assert m.bulk.area[-1]==m.area[0] and np.max(abs(d['edges'][:16]-original['edges'][:16]))==0
    row=dict(classification='Counterexample candidate',passed=True,original15_native_banks_bitwise=True,shared_photon_face_area_identical=True,
        shared_material_flux_units_relative=error,new_cells_filled_by_actual_native_EOS=True,
        nonlinear_temperature_endpoint_K=(d['T'][-4:]*np.exp(z['theta'][-4:])).tolist(),
        relative_density_endpoint=m.bulk.eos.x[-4:].tolist(),nearest_center_to_shared_face_m=float((m.m.rf[0]-d['r'][-1])/100),
        complete_neutral_trajectory_audit=False,global_radial_error_certified=False,full_GR_feedback=False,final_charge_solved=False)
    write(OUT/'local-audit.json',row);print(json.dumps(row),flush=True)


def gr():
    assert json.loads((OUT/'radial.json').read_text())['passed']
    FunctionType(inherited.gr.__code__,dict(vars(inherited),run=runtime,OUT=OUT))()


if __name__=='__main__':globals()[sys.argv[1]]()
