"""Independent radiation on the current GR background without LTE double count."""
from pathlib import Path
import argparse
import json
import time
import numpy as np
import sympy as sp
import gr_radiation_eos_split as old
import def_ionic_structure_transport as model

ex=model.ex;h=model.h;OUT=model.OUT.parent/'def-photon-matter-split'
CELLS=[0,500,1722,2600,4122,5734]


def symbolic():
    e,P,E,q,pr=sp.symbols('e P E q pr');eq=sp.symbols('E_LTE')
    g=sp.diag(-1,1,1,1);u=sp.Matrix([1,0,0,0])
    reference=(e+P)*(u*u.T)+P*g
    equilibrium=4*eq/3*(u*u.T)+eq*g/3
    photons=sp.Matrix([[E,q,0,0],[q,pr,0,0],[0,0,(E-pr)/2,0],[0,0,0,(E-pr)/2]])
    matter=(e-eq+P-eq/3)*(u*u.T)+(P-eq/3)*g
    assert sp.simplify(reference-equilibrium-matter)==sp.zeros(4)
    assert sp.trace(g*photons)==0 and sp.trace(g*equilibrium)==0
    assert sp.simplify(reference+(photons-equilibrium)-matter-photons)==sp.zeros(4)
    return dict(classification='Proven',passed=True,
        split='T_total=T_native_LTE+[R-R_LTE]=T_matter+R. Subtract the complete LTE radiation tensor, not energy alone.',
        evolution='For interaction four-force G, div R=-G and div T_matter=G. Equivalently div T_native_LTE=-div(R-R_LTE), in addition to the declared scalar/matter forces. Initial R=R_LTE keeps the accepted mechanical background unchanged.',
        scalar_trace='For massless radiation, trace R=trace R_LTE=0, so trace T_matter=trace T_native_LTE at the same rho,T,X. Radiation still changes matter temperature, metric and motion; this is not a no-signal statement.',
        scope='Vacuum-dispersion photons with the native LTE radiation free energy. In-medium polarization/dispersion and complete physical EOS errors are not supplied by this algebra.')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    paths=[Path(__file__),Path(old.__file__),old.BRIDGE,old.model.LIB,model.base.thermal.OUT/'coefficients.npz']
    ex.write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='0f6b09ed',
        claim='Provide matter-only pressure, energy and thermal derivatives for independent photon stress on the current saved GR states; retain the same native EOS constants and initial total stress.',
        bindings={str(p.relative_to(h.ROOT)) if p.is_relative_to(h.ROOT) else str(p):h.digest(p) for p in paths},
        cells=CELLS,gates=dict(pressure_derivative_relative=1e-7,energy_erg_g_or_ULP='max(2,32 ulp)',capacity_positive=True),
        budget=dict(native_EOS_calls=6,hard_seconds=30,CPU_workers=1,new_stellar_steps=0,automatic_expansion=False),
        input_convention='Convert stored pressure-mode rho_P,rho_T,u_P,u_T to density-mode derivatives before subtracting the same compiled LTE photon free energy. Do not interpret pressure-mode columns as density-mode.',
        limits='Exact split within the native vacuum-photon EOS candidate. No plasma-dispersion correction, opacity model certificate or completed independent photon evolution.'))
    ex.write(OUT/'symbolic.json',symbolic())


def run():
    assert not (OUT/'result.json').exists();plan=json.loads((OUT/'plan.json').read_text())
    for p,sha in plan['bindings'].items():assert h.digest(h.ROOT/p)==sha,p
    start=time.monotonic();data,physical=model.base.inputs();gas=old.GasEOS()
    density_mode=np.copy(data['raw']);density_mode[:,[5,6,9,10]]=data['thermo'][:,:4]
    split=old.split(density_mode,physical['T'],gas.a_rad);checks=[]
    assert np.min(split['P_gas'])>0 and np.min(split['cvT_gas'])>0
    for i in CELLS:
        native=gas(2,float(data['lnd'][i]),float(data['lnT'][i]),data['X'][i])
        mapping=dict(P_gas=1,chi_rho_gas=5,chi_T_gas=6,du_dlnrho_gas=9,cvT_gas=10)
        errors={k:float(abs(native[j]-split[k][i])/max(abs(split[k][i]),1)) for k,j in mapping.items()}
        budget=max(2.,32*np.spacing(abs(native[2]+native[1]/native[0])))
        energy_score=float(abs(native[2]-split['u_gas'][i])/budget)
        entropy_score=float(physical['T'][i]*abs(native[3]-split['s_gas'][i])/budget)
        checks.append(dict(index=i,passed=max(errors.values())<1e-7 and max(energy_score,entropy_score)<1,
            relative_errors=errors,energy_score=energy_score,entropy_score=entropy_score))
    np.savez_compressed(OUT/'matter.npz',**split,rho_B=np.exp(data['lnd']),T=physical['T'],compiled_a_rad_cgs=gas.a_rad,
        compiled_c_cm_s=gas.c_light,initial_flux=np.zeros(len(physical['T'])),initial_radial_pressure=split['pressure_radiation'])
    result=dict(classification='Counterexample candidate',passed=all(r['passed'] for r in checks),checks=checks,
        full_grid_algebraic_cells=len(physical['T']),native_control_cells=len(CELLS),seconds=time.monotonic()-start,
        matter_pressure_capacity_positive=True,independent_photon_evolution_completed=False,whole_star_heat_closed=False,full_dynamic_charge_solved=False)
    ex.write(OUT/'result.json',result);print('RESULT',result,flush=True);assert result['seconds']<30


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','run']);globals()[parser.parse_args().action]()
