"""EOS-matched Maxwell/Vlasov density spectrum for collective photon scattering.

Positive velocity quadrature turns the collisionless response into a symmetric
diagonal-plus-rank-one oscillator matrix. Integrate its positive spectral atoms
over photon cells instead of sampling potentially unresolved plasma resonances.
"""
from pathlib import Path
import argparse
import ctypes
import json
import signal
import time
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import atomic_mass
from scipy.linalg import eigh
from scipy.special import roots_hermitenorm, wofz
import def_photon_finite_jump as previous
import eos_species_inventory as inventory_reader

old=previous.old;base=old.base;OUT=previous.OUT.parent/'def-photon-collective'


def prepare():
    assert not (OUT/'plan.json').exists();OUT.mkdir(exist_ok=True)
    (OUT/'inventory-source.py').write_bytes(Path(__file__).read_bytes())
    paths=[OUT/'inventory-source.py',Path(inventory_reader.__file__),old.previous.matter.old.BRIDGE,
        old.OUT/'bank.npz',previous.OUT/'manifest.json']
    old.ex.write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='9b1fce346',
        claim='Construct collective electron/ion photon scattering from the SAME gas EOS populations and connect it to the previously accepted material/photon response. Preserve unresolved bound-electron/opacity/atmosphere/GR requirements.',
        first_decision='Determine free/charged/neutral inventories and plasma scales; compare the independently derived RPA static scattering with the saved scalar opacity without fitting it.',
        method='Read existing native ionization module immediately after the current molecular gas-only call. Reuse the prior inventory checker, not its older EOS library. Velocity quadrature defines a positive Hamiltonian Vlasov spectral measure; no arbitrary Lorentz broadening of unresolved collective poles.',
        sources=['https://docs.plasmapy.org/en/stable/_modules/plasmapy/diagnostics/thomson.html',
            'https://arxiv.org/html/1601.01005'],
        budget=dict(native_EOS_calls=4,new_opacity_queries=0,inventory_seconds=60,
            spectral_pilot_seconds=120,CPU_workers=1,GPU=False,new_stellar_steps=0,
            production='Assess measured cost and error of the spectral measure BEFORE evolving. No automatic grid/duration expansion.'),
        gates=dict(inventory=1e-10,charge=1e-10,reported_ions=1e-12,saved_ne=1e-12,
            algebra_static=1e-9,algebra_second_moment=1e-9),
        limits='Collisionless weak-coupling Maxwell/RPA reference. Atomic form factors, polarizability, nonideal corrections and material/GR recoil must not be silently declared complete.',
        bindings={p.relative_to(old.h.ROOT).as_posix() if p.is_relative_to(old.h.ROOT) else str(p):old.h.digest(p) for p in paths}))
    print('PREPARED',OUT,flush=True)


def inventory():
    assert not (OUT/'inventory.npz').exists();signal.alarm(60);start=time.monotonic()
    d,p=base.inputs();gas=old.previous.matter.old.GasEOS();lib=gas.gas_lib
    gas.inventory_lib=lib;native=lib.ionization_inventory
    def capture(mode,value,t,eps,out,info):
        raw=np.full(24,np.nan);native(mode,value,t,eps,raw,info)
        out[:]=raw[:22];out[20]=0.;gas.molecules=raw[22:].copy()
    gas.call=capture
    snapshots=[];rows=[]
    for T in [float(p['T'][0]),.0015*1.602176634e-16/base.k,.002*1.602176634e-16/base.k,float(p['T'][0])]:
        snap=inventory_reader.InventoryEOS.snapshot(gas,float(d['lnd'][0]),float(np.log(T)),d['X'][0])
        row=inventory_reader.check(snap,d['X'][0],float(d['lnd'][0]),gas)
        assert not row['missing_nonzero_elements'] and row['inventory_error']<1e-10 and row['charge_error']<1e-10 and row['reported_fraction_error']<1e-12
        rows.append(dict(T=T,**row));snapshots.append(snap)
    assert all(np.array_equal(snapshots[0][key],snapshots[-1][key]) for key in snapshots[0])
    s=snapshots[0];ym=(d['X'][0]/inventory_reader.g.c.A)@gas.mapping;cx=float(ym@gas.weights);eps=ym/cx
    factor=np.exp(d['lnd'][0])*cx*base.Avogadro*1e6
    ni=s['number_fractions']*factor;charges=np.arange(29)
    ne=float(s['eos'][13]*base.Avogadro*1e6);bank=np.load(old.OUT/'bank.npz')
    assert abs(ne/float(bank['ne'])-1)<1e-12
    # Species with equal element mass share the Maxwell susceptibility; sum Z^2
    # weights exactly instead of replacing the mixture by an effective ion.
    strength=(ni*charges**2).sum(axis=1)/ne
    mass=gas.weights*atomic_mass;active=strength>0
    strength=strength[active];mass=mass[active]
    molecular_ion=factor*eps[0]*s['molecular_H_fractions'][1]/2
    if molecular_ion>0:strength=np.r_[strength,molecular_ion/ne];mass=np.r_[mass,2*gas.weights[0]*atomic_mass]
    T=float(bank['T']);ld=np.sqrt(base.k*T/(4*np.pi*base.E2*ne));sigma=np.sqrt(base.k*T/base.m_e)
    wp=sigma/ld;up=base.hbar*wp/(base.k*T);wave=base.k*T/(base.hbar*base.c)
    radius=(3/(4*np.pi*ne))**(1/3);thermal=np.sqrt(2*np.pi*base.hbar**2/(base.m_e*base.k*T))
    np.savez_compressed(OUT/'inventory.npz',ni=ni,charges=charges,neutral_masses=gas.weights*atomic_mass,
        ion_strength=strength,ion_mass=mass,ne=ne,T=T,lambda_D=ld,sigma_e=sigma,plasma_frequency=wp,
        plasma_u=up,wave_lambda=wave*ld,full_EOS=s['eos'],molecular_H_fractions=s['molecular_H_fractions'],
        native_temperature_number_fractions=np.array([v['number_fractions'] for v in snapshots[:3]]))
    info=dict(classification='Counterexample candidate',passed=True,rows=rows,history_bitwise=True,
        free_electrons_m3=ne,neutral_m3=ni[:,0].tolist(),ion_strength=strength.tolist(),
        mass_groups=len(mass),effective_Z_squared_weight=float(strength.sum()),
        lambda_D_m=ld,plasma_u=up,wave_lambda=wave*ld,coupling_parameter=float(base.E2/(radius*base.k*T)),
        degeneracy_parameter=float(ne*thermal**3/2),seconds=time.monotonic()-start,
        actual_EOS_library=True,native_opacity_populations_identified=False)
    old.ex.write(OUT/'inventory.json',info);print('INVENTORY',info,flush=True);signal.alarm(0)


def static():
    s=np.load(OUT/'inventory.npz');bank=np.load(old.OUT/'bank.npz');mu,w=leggauss(128)
    strength=float(s['ion_strength'].sum());u=bank['u'];x2=(u[:,None]*float(s['wave_lambda']))**2*2*(1-mu)
    S=(x2+strength)/(x2+1+strength);pred=(S*(3/8*(1+mu**2)*w)).sum(1)
    actual=bank['rate_s']/float(bank['rate_e']);rows=[]
    for target in [.001,.01,.1,.3,1,3,5,10,20,50]:
        j=int(np.argmin(abs(u-target)));rows.append(dict(u=float(u[j]),RPA=float(pred[j]),native=float(actual[j]),relative=float(pred[j]/actual[j]-1)))
    np.savez_compressed(OUT/'static.npz',u=u,RPA_scattering_ratio=pred,native_scattering_ratio=actual)
    old.ex.write(OUT/'static.json',dict(classification='Counterexample candidate',rows=rows,
        low_q_limit=strength/(1+strength),fit_performed=False,
        interpretation='Independent EOS-population RPA sum rule, integrated over dipole angle. A comparison with a total scalar opacity is not identification of its atomic components.'))
    print('STATIC',rows,flush=True)


def modes(x,order,s):
    """Positive RPA measure at x=q*lambda_D, phase velocity in sigma_e units."""
    velocity,weights=roots_hermitenorm(order)
    keep=velocity>0;velocity=velocity[keep];weights=weights[keep]*2/np.sqrt(2*np.pi)
    speeds=np.r_[1.,np.sqrt(base.m_e/s['ion_mass'])]
    strength=np.r_[1.,s['ion_strength']]
    v=(speeds[:,None]*velocity).ravel();w=np.tile(weights,len(speeds))
    b=np.repeat(np.sqrt(strength),len(velocity))*np.sqrt(w)/x
    b[:len(velocity)]*=-1;g=v*b
    # Rotate the charge direction before adding the large Coulomb eigenvalue.
    # This avoids subtracting large entries to recover weak ion modes at small q.
    h=g/np.linalg.norm(g);h[0]+=np.copysign(1.,h[0]);h/=np.linalg.norm(h)
    hv=h*v*v;A=np.diag(v*v)-2*np.outer(h,hv)-2*np.outer(hv,h)+4*(h@hv)*np.outer(h,h)
    A[0,0]+=g@g
    lam,U=eigh(A,driver='evr');assert np.all(lam>0)
    electron=np.zeros(len(v));electron[:len(velocity)]=velocity*np.sqrt(weights)
    electron-=2*h*(h@electron)
    W=(electron@U)**2/lam
    return np.sqrt(lam),W


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','inventory','static']);globals()[p.parse_args().action]()
