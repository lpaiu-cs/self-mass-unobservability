"""Unclamped stellar response versus its own static surface impedance.

Reuse Phase83 full Dirichlet-to-Neumann data. Transport through vacuum to the
physical stellar surface before freezing the static comparator: freezing at
an arbitrary computational outer radius would include that vacuum as matter.
"""
from pathlib import Path
import argparse
import json
import signal
import time
import numpy as np
import mpmath as mp
import sympy as sp
from scipy.integrate import solve_ivp
import def_orbital_charge_fem as prior

ROOT=prior.exterior.s.ROOT
OUT=prior.OUT.parent/'def-full-orbital-exterior'
write=prior.write;ext=prior.exterior


def symbolic():
    Z0,dZ,Zo,drive,h=sp.symbols('Z0 dZ Zo drive h',nonzero=True)
    D=Z0-Zo
    assert sp.simplify((drive/(D+dZ)-drive/D)/h+drive*dZ/(h*D*(D+dZ)))==0
    q,M,dq,dM=sp.symbols('q M dq dM',nonzero=True);eps=sp.symbols('eps')
    assert sp.diff((q+eps*dq)/(M+eps*dM),eps).subs(eps,0)==dq/M-q*dM/M**2
    # The standing scattering solution has unit modulus for real interior Z.
    x,y,z=sp.symbols('x y z',real=True)
    S=-(z-x+sp.I*y)/(z-x-sp.I*y)
    assert sp.simplify(S*sp.conjugate(S)-1)==0
    return dict(classification='Proven',passed=True,
        identities='Same incoming: outgoing contrast=-drive*dZ/((Z0-Zout)*(Z0+dZ-Zout)*h). d(Q/M)=dQ/M-Q*dM/M^2. Real lossless interior impedance gives unitary one-channel scattering.',
        scope='Linear algebra and normalization identities, not thermal equilibrium or full binary matching.')


def background():
    saved=np.load(prior.patch.base.surface.OUT/'background.npz')
    rs=float(saved['r'][-1]);mu=float(saved['m'][-1]/rs);flux=float(saved['v'][-1]*rs)
    M,K,bg=ext.background(mu,flux,2e-13)
    return mu,flux,M,K,bg


def transport(z,omega,outer,bg,K,tol):
    # Carry (psi,psi') at unit outer Dirichlet amplitude back to the surface.
    def rhs(r,y):
        m,F=bg.sol(1/r);b=1-2*m/r;Phi=K/(F*r*r)
        return [y[1],-(2/r+2*m/(r*r*b))*y[1]-(omega*omega/(F*F)+2*Phi*Phi/b)*y[0]]
    sol=solve_ivp(rhs,[outer,1],[1.,z/outer],method='DOP853',rtol=tol,atol=tol*1e-5)
    assert sol.success
    return sol.y[1,-1]/sol.y[0,-1],sol.y[0,-1],sol.nfev


def static_match(mu,flux,Z,amplitude):
    # Independent exact Just differential, including the scalar contribution
    # to ADM mass and fixed-areal-radius metric constraint.
    mp.mp.dps=70;mu,q,Z,a=[mp.mpf(str(v)) for v in [mu,flux,Z,amplitude]]
    E=ext.q.exterior.exact
    mass=lambda u,v:u+v*v*E(u,v)[0]
    charge=lambda u,v:v*E(u,v)[1]
    far=lambda u,v:v*E(u,v)[2]
    du=(1-2*mu)*q*a;dv=Z*a
    dM=mp.diff(lambda u:mass(u,q),mu)*du+mp.diff(lambda v:mass(mu,v),q)*dv
    dQ=mp.diff(lambda u:charge(u,q),mu)*du+mp.diff(lambda v:charge(mu,v),q)*dv
    dphi=a+mp.diff(lambda u:far(u,q),mu)*du+mp.diff(lambda v:far(mu,v),q)*dv
    alpha=charge(mu,q)/mass(mu,q)
    bare=dQ/(mass(mu,q)*dphi);norm=-alpha*dM/(mass(mu,q)*dphi)
    identity=dM/(mass(mu,q)*dphi)-alpha
    assert abs(identity)<mp.mpf('1e-50')
    return dict(classification='Counterexample candidate',alpha_background=float(alpha),
        incident_static_calibration=float(dphi),mass_log_derivative=float(dM/(mass(mu,q)*dphi)),
        unnormalized_charge_gain=float(bare),mass_normalization_gain=float(norm),normalized_charge_gain=float(bare+norm),
        exact_mass_work_identity_error=mp.nstr(abs(identity),10),
        scope='Adiabatic zero-frequency family selected by the saved GR constraint; no complete thermal-equilibrium family or finite-frequency ADM-charge identity.')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        claim='Compute the full unclamped outgoing response against the same stars static surface impedance, and exact static ADM charge/mass differential.',
        comparator='Freeze full interior Z(0) at physical surface r/R=1; retain exact frequency-dependent vacuum propagation and identical incoming amplitudes. This is the same adiabatic inventory/entropy/composition mechanical reference, not a thermal equilibrium star or a whole static EFT.',
        method='Reuse all16 Phase83 impedances. Back-propagate exact vacuum scalar+metric ODE to surface, compare tolerances, then use direct rational impedance difference. No new EOS, FEM solve or time evolution.',
        gates=dict(surface_transport_relative=1e-6,static_incident_calibration=1e-8,current=2e-8,spatial=.02,outer=.002,quadrature=.002,unitarity=1e-12),
        budget=dict(total_seconds=120,new_EOS_calls=0,new_FEM_solves=0,new_time_steps=0,CPU_threads=1),
        decision='Stop if radius/tolerance controls fail; do not extend the frequency set or orbital horizon. Report full radiative response and finite-radius mass separately from an EOS-matched binary sensitivity.',
        bindings={str(p):prior.go.task.digest(p) for p in [Path(__file__),Path(ext.__file__),Path(ext.q.exterior.__file__),prior.OUT/'result.json',prior.OUT/'audit.json',prior.OUT/'pilot.json',prior.BENCH,
            prior.patch.base.surface.OUT/'background.npz']}))
    write(OUT/'symbolic.json',symbolic())


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for p,h in plan['bindings'].items():assert prior.go.task.digest(Path(p))==h,p
    assert not (OUT/'result.json').exists();signal.alarm(120);start=time.monotonic()
    old=json.loads((prior.OUT/'result.json').read_text());pilot=json.loads((prior.OUT/'pilot.json').read_text())
    mu,flux,M,K,bg=background();F=float(bg.y[1,-1]);ratio=pilot['R_m']/pilot['ADM_geom_m'];rows={};maximum=0.
    waves={n:ext.outgoing(mu,flux,old['cases']['p4'][n]['omega_R_over_c'],tol=2e-13) for n in range(4)}
    for label,prior_rows in old['cases'].items():
        outer=3 if label=='outer' else 2;surface=[]
        for row in prior_rows:
            z=row['Z_clamped'][0]+row['Z_difference'][0];w=row['omega_R_over_c']
            a,amp,_=transport(z,w,outer,bg,K,2e-12)
            b,amp2,_=transport(z,w,outer,bg,K,2e-13)
            err=abs(a-b)/max(abs(b),1e-100);maximum=max(maximum,err)
            surface.append((b,amp2))
        z0=surface[0][0];case=[]
        for n,(z,amp) in enumerate(surface):
            w=prior_rows[n]['omega_R_over_c'];wave=waves[n];h=wave['h'];Zo=wave['impedance']
            drive=np.exp(-1j*w)/(F*h);D=z0-Zo;dz=z-z0
            af=drive/(z-Zo);astatic=drive/D
            delta=-drive*dz/(D*(D+dz));tail=delta/h;charge=-ratio*tail
            # Quasilocal Misner-Sharp mass from the exact vacuum constraint.
            # The radiative ADM mass contrast vanishes at infinity; this finite
            # surface mass includes the exterior scalar field's energy budget.
            mass_ratio=(1-2*mu)*flux*delta/mu
            scattering=-(z-np.conj(Zo))/(z-Zo)*np.conj(h)/h
            pair=lambda v:[float(v.real),float(v.imag)]
            case.append(dict(harmonic=n,surface_impedance=float(z),surface_impedance_difference=float(dz),
                full_surface_scalar=pair(af),static_surface_scalar=pair(astatic),
                outgoing_full_minus_static=pair(tail),radiative_charge_gain=pair(charge),
                surface_quasilocal_mass_fraction_gain=pair(mass_ratio),
                S_outgoing_over_incoming=pair(scattering),unitarity_error=float(abs(abs(scattering)-1)),
                original_outer_solution_surface_relative=float(abs(af-complex(*prior_rows[n]['boundary_amplitude'])*amp)/abs(af)),
                incident_amplitude=old['rows'][n]['drive_amplitude'],
                delta_alpha_radiative_over_phi0=float(abs(charge)*old['rows'][n]['drive_amplitude']/.001),
                surface_mass_fraction_amplitude=float(abs(mass_ratio)*old['rows'][n]['drive_amplitude'])))
        rows[label]=case
    matching=static_match(mu,flux,rows['p4'][0]['surface_impedance'],rows['p4'][0]['full_surface_scalar'][0])
    comparisons=[]
    for n in [1,2,3]:
        ref=complex(*rows['p4'][n]['radiative_charge_gain']);item=dict(harmonic=n)
        for label in ['p2','outer','quadrature']:item[label]=abs(complex(*rows[label][n]['radiative_charge_gain'])-ref)/abs(ref)
        comparisons.append(item)
    g=plan['gates'];passed=maximum<g['surface_transport_relative'] and abs(matching['incident_static_calibration']-1)<g['static_incident_calibration']
    passed=passed and all(r['p2']<g['spatial'] and r['outer']<g['outer'] and r['quadrature']<g['quadrature'] for r in comparisons)
    passed=passed and max(r['unitarity_error'] for case in rows.values() for r in case)<g['unitarity']
    passed=passed and max(w['current_relative_error'] for w in waves.values())<g['current']
    result=dict(classification='Counterexample candidate',passed=passed,rows=rows,comparisons=comparisons,static_matching=matching,
        maximum_transport_relative=maximum,seconds=time.monotonic()-start,
        whole_adiabatic_interior_response=True,clamped_subtraction=False,static_comparator_at_physical_surface=True,
        full_thermal_background=False,full_binary_sensitivity_matched=False,observational_nuisance_applied=False,full_goal_complete=False)
    write(OUT/'result.json',result);print('RESULT',json.dumps({k:v for k,v in result.items() if k!='rows'}),flush=True)
    print('READOUTS',json.dumps(rows['p4']),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
