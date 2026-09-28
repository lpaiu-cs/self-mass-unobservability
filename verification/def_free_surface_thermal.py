"""Same-background thermal coefficients and actual native reaction sources.

Pressure-mode EOS derivatives are transformed, never read as density-mode
columns. Source values use the actual saved electron abundance/degeneracy;
the native derivative columns are deliberately not consumed.
"""
from pathlib import Path
import argparse
import json
import time
import numpy as np
import sympy as sp
import def_free_surface_response_normalized as surface
import gr_two_carrier_evolution as two
import source_retry as retry

h=surface.h
g=retry.g
OUT=g.OUT/'def-free-surface-thermal'
CACHE=Path('/home/lpaiu/work/def-free-surface-thermal50')


def inputs():
    bg=dict(np.load(h.OUT/'absolute-shoot/background-0.001.npz'))
    raw=np.load(h.OUT/'absolute-shoot/native-audit.npz')['raw']
    data,_=h.inputs();body=h.Structure(.001)
    state=dict(lnd=np.log(raw[:,0]),lnT=bg['thermo'][:,1],X=bg['X'],dm=bg['dm'],
        lnR=np.log(bg['states'][:,0]*body.R*100),L=np.zeros(len(raw)))
    return bg,raw,state


def transformed(raw):
    rp,rt=raw[:,7],raw[:,8]
    cr=1/rp;ct=-rt/rp
    ur=raw[:,9]/rp;ut=raw[:,10]-raw[:,9]*rt/rp
    ad=(raw[:,1]/raw[:,0]-ur)/ut
    return np.c_[cr,ct,ur,ut,ad,raw[:,10]-raw[:,1]/raw[:,0]*rt]


def symbolic():
    rp,rt,up,ut=sp.symbols('rp rt up ut',nonzero=True)
    dr,dt=sp.symbols('dr dt');dp=(dr-rt*dt)/rp
    assert sp.expand(up*dp+ut*dt-up/rp*dr-(ut-up*rt/rp)*dt)==0
    A,N,K,a,T,Tr,ar,nr=sp.symbols('A N K a T Tr ar nr',positive=True)
    flux=-N**2*A**2*K/a*(A*Tr+A*T*ar+A*T*nr)
    assert sp.simplify(flux+N*A**2*K/a*(A*N*(Tr+T*(ar+nr))))==0
    return dict(classification='Proven',passed=True,
        EOS='chi_rho=1/rho_P; chi_T=-rho_T/rho_P; u_lnrho=u_lnP/rho_P; cv*T=u_lnT_at_P-u_lnP*rho_T/rho_P.',
        transport='L_infinity=-4*pi*r_E^2*N*A^2*K_J/a_E * d(A*N*T_J)/dr_E. This is the declared LTE diffusion law, not a thin-atmosphere closure.',
        source='Native eta derivative arguments are not used as physical response derivatives; value independence is checked by replaying all source states with changed derivative arguments.')


def prepare():
    assert not OUT.exists() and not CACHE.exists();OUT.mkdir();CACHE.mkdir()
    files=[Path(__file__),h.OUT/'absolute-shoot/background-0.001.npz',
        h.OUT/'absolute-shoot/native-audit.npz',h.OUT/'absolute-shoot/lapse.npz',
        surface.OUT/'background.npz',Path(two.__file__),Path(two.radiative.tables.__file__),
        retry.OUT/'plan.json',Path(g.s.v.__file__),g.c.fresh.BINARY]
    h.write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='b3552307',
        bindings={str(p.relative_to(h.ROOT)) if p.is_relative_to(h.ROOT) else str(p):h.digest(p) for p in files},
        symbolic=symbolic(),
        claim='Connect finite-temperature thermal coefficients and all 26-species native rates on the same nonzero mechanical background; use them for conservative thermal/reactive increments and free-surface mass/charge response.',
        actual_sources='Current saved native rho,T,X,free-electron abundance,eta. Value-only native evaluation; returned derivative matrices are discarded. Replay with eta derivative arguments changed from zero to (+1,-1) must leave rate/heat/loss outputs bitwise unchanged.',
        thermal_boundary='Photospheric outgoing luminosity remains an explicit parameter; diffusion is used only as an interior constitutive candidate. No Rosseland-to-Planck identification and no cold opacity extrapolation.',
        gates=dict(EOS_chain_relative=1e-7,Gamma_identity_relative=1e-7,source_derivative_argument_values_bitwise=True,baryon_relative=1e-10),
        budget=dict(native_EOS_controls=6,source_pilot_cells=32,full_source_batches=2,source_timeout_seconds=120,coefficient_timeout_seconds=60,CPU_threads=1,GPU=False,automatic_expansion=False),
        cost_basis='Earlier identical 5735-cell native source trace took 8.50 s; changed-state cost remains to be measured with the 32-cell pilot. Reuse all saved native thermodynamics.',
        full_dynamic_charge_solved=False))


def coefficients():
    assert not (OUT/'coefficients.npz').exists();began=time.monotonic()
    bg,raw,state=inputs();thermo=transformed(raw);plan=json.loads((OUT/'plan.json').read_text())
    gamma=thermo[:,0]+thermo[:,1]*thermo[:,4]
    identity=float(np.max(abs(gamma/raw[:,4]-1)));assert identity<plan['gates']['Gamma_identity_relative']
    eos=h.molecular.model.EOS();control=[]
    for i in [0,600,1175,2300,4000,5734]:
        b=eos(2,state['lnd'][i],state['lnT'][i],state['X'][i])
        calculated=thermo[i,:4];actual=b[[5,6,9,10]]
        error=float(np.max(abs(actual-calculated)/np.maximum(abs(actual),1)))
        control.append(dict(cell=i,relative_chain_error=error));assert error<plan['gates']['EOS_chain_relative'],control[-1]
    op=two.radiative.tables.Opacity()
    opacity=np.array([two.opacity_parts(op,row) for row in zip(state['lnd'],state['lnT'],state['X'])])
    lap=np.load(h.OUT/'absolute-shoot/lapse.npz');body=h.Structure(.001)
    phi=.001*(1+body.mu*bg['states'][:,3]);A=np.exp(-2*phi*phi);N=np.exp(lap['nu_mid'])
    r=bg['states'][:,0]*body.R*100;rf=bg['faces'][:,0]*body.R*100
    mass=bg['states'][:,1]*body.B*100;metric=1/np.sqrt(1-2*mass/r)
    T=np.exp(state['lnT']);rho=raw[:,0];K=16*5.670400e-5*T[:,None]**3/(3*rho[:,None]*opacity[:,[0,3]])
    theta=A*N*T;gradient=np.diff(theta)/np.diff(r)
    # Outward positive luminosity, shell order outer to inner.
    fac=4*np.pi*rf[1:-1]**2*np.sqrt(N[:-1]*N[1:])*np.exp(np.log(A[:-1]*A[1:]))/np.sqrt(metric[:-1]*metric[1:])
    L=-fac[:,None]*np.sqrt(K[:-1]*K[1:])*gradient[:,None]
    heat=(np.r_[L.sum(1),0]-np.r_[0,L.sum(1)])/(state['dm']*A*N)
    cap=thermo[:,3];cp=thermo[:,5]
    causal=4*np.pi*rf[0]**2*N[0]**2*A[0]**4*4*5.670400e-5*T[0]**4
    np.savez_compressed(OUT/'coefficients.npz',thermo=thermo,opacity=opacity,K=K,A=A,N=N,metric=metric,
        radius_cm=r,faces_cm=rf,luminosity=L,zero_outer_heat_erg_g_s=heat,
        outer_LTE_c_energy_luminosity=causal,raw=raw,**state)
    defect=float(abs(np.sum(state['dm']*A*N*heat,dtype=np.longdouble))/max(np.sum(abs(L)),1))
    assert defect<1e-13
    h.write(OUT/'coefficients.json',dict(classification='Counterexample candidate',cells=len(raw),
        EOS_controls=control,Gamma_identity_relative=identity,closed_face_energy_telescoping=defect,
        capacity_positive=bool(np.min(cap)>0 and np.min(cp)>0),
        maximum_frozen_isochoric_logT_rate=float(np.max(abs(heat/cap))),
        outer_LTE_c_energy_luminosity_erg_s=float(causal),
        seconds=time.monotonic()-began,native_EOS_calls=6,
        surface_luminosity_prescribed=False,physical_radiation_closure=False,full_dynamic_charge_solved=False))
    print('COEFFICIENTS',identity,defect,'seconds',time.monotonic()-began,flush=True)


def configure():
    retry.configure();g.s.OUT=OUT;g.s.CACHE=CACHE
    g.c.native.OUT=OUT;g.c.native.CACHE=CACHE;g.c.native.context()
    shim=g.s.v.shim;shim.OUT=OUT;shim.CACHE=CACHE
    shim.save=lambda name,value:h.write(OUT/name,value)


def source(label,indices,derivative_arguments):
    configure();bg,raw,state=inputs();data={k:v[indices] for k,v in state.items()}
    aux=np.c_[raw[indices,13]/raw[indices,0],raw[indices,12],
        np.full(len(indices),derivative_arguments[0]),np.full(len(indices),derivative_arguments[1])]
    g.s.auxiliary=lambda _:aux
    began=time.monotonic();native,corrected=g.s.evaluate(label,data)
    assert np.array_equal(native['aux_used'],aux)
    assert np.max(abs(np.log(native['rho'])-data['lnd']))<1e-12
    assert np.max(abs(np.log(native['T'])-data['lnT']))<1e-12
    assert np.max(abs(native['X']-data['X']))<1e-12
    return native,corrected,time.monotonic()-began


def pilot():
    assert not (OUT/'source-pilot.json').exists()
    selected=np.unique(np.linspace(0,5734,32).astype(int))
    _,_,elapsed=source('pilot',selected,[0,0])
    forecast=4*8.50+2*elapsed
    h.write(OUT/'source-pilot.json',dict(classification='Counterexample candidate',cells=selected.tolist(),seconds=elapsed,
        full_two_batch_forecast_seconds=forecast,basis='Four times the prior full-grid trace plus two current pilot overheads; changed full-grid rates remain unmeasured.',within_budget=forecast<110))
    print('SOURCE PILOT',elapsed,'forecast',forecast,flush=True)


def sources():
    assert not (OUT/'sources.json').exists()
    assert json.loads((OUT/'source-pilot.json').read_text())['within_budget']
    rows=np.arange(5735);a,x,first=source('actual',rows,[0,0]);b,y,second=source('argument-control',rows,[1,-1])
    equal={key:np.array_equal(x[key],y[key]) for key in ['dxdt','heat','neutrino']}
    assert all(equal.values()),equal
    _,profile=g.c.mesa(OUT/'actual-profile.data.gz');loss=profile['non_nuc_neu']
    coeff=np.load(OUT/'coefficients.npz');dm=coeff['dm'];A=coeff['A'];N=coeff['N'];R=x['dxdt']
    rest=(g.c.W/g.c.A-1)*(h.gr.C*100)**2
    # Total atomic rest energy plus native internal energy is conserved; do not
    # add nuclear Q again. A later finite EOS update accounts for u(X) as well.
    heating=-R@rest-x['neutrino']-loss
    baryon=float(np.max(abs(R.sum(1))/np.maximum(abs(R).sum(1),1e-100)))
    assert baryon<1e-10,baryon
    np.savez_compressed(OUT/'sources.npz',dxdt=R,neutrino=x['neutrino'],thermal_neutrino=loss,
        rest_to_internal_heating=heating,heat_Q_reference=x['heat'])
    h.write(OUT/'sources.json',dict(classification='Counterexample candidate',cells=len(R),
        rate_values_bitwise_independent_of_derivative_arguments=equal,baryon_source_relative=baryon,
        proper_heating_range_erg_g_s=[float(heating.min()),float(heating.max())],
        redshifted_net_power_erg_s=float(dm@(A*A*N*N*heating)),
        maximum_species_rate=float(abs(R).max()),seconds=first+second,
        internal_energy_composition_derivative_included=False,returned_Jacobian_used=False,
        finite_thermal_reactive_step=False,full_dynamic_charge_solved=False))
    print('SOURCES',len(R),'seconds',first+second,'max heating',heating.max(),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','coefficients','pilot','sources'])
    globals()[parser.parse_args().action]()
