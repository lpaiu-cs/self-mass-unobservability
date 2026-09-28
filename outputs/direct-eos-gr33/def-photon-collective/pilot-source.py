"""EOS-matched Maxwell/Vlasov density spectrum for collective photon scattering.

Positive velocity quadrature turns the collisionless response into a symmetric
diagonal-plus-rank-one oscillator matrix. Integrate its positive spectral atoms
over photon cells instead of sampling potentially unresolved plasma resonances.
"""
from pathlib import Path
from types import FunctionType
import argparse
import ctypes
import json
import signal
import time
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import atomic_mass
from scipy.linalg import svd
from scipy.special import roots_hermitenorm, ndtr
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
    b[:len(velocity)]*=-1
    # Whiten the electrostatic Gibbs covariance BEFORE spectral decomposition.
    # Singular values of H^(1/2) D avoid squaring the condition number; direct
    # observation weights avoid division by inaccurate tiny eigenvalues.
    h=b/np.linalg.norm(b);h[0]+=np.copysign(1.,h[0]);h/=np.linalg.norm(h)
    Q=np.eye(len(v))-2*np.outer(h,h);scale=np.sqrt(1+b@b)
    F=Q*v;F[0]*=scale
    U,singular,_=svd(F,full_matrices=False,lapack_driver='gesvd')
    electron=np.zeros(len(v));electron[:len(velocity)]=np.sqrt(weights)
    electron=Q@electron;electron[0]/=scale
    W=(electron@U)**2
    return singular,W


def spectrum(order,points):
    target=OUT/f'spectrum-{order}-{points}.npz'
    if target.exists():return dict(np.load(target))
    start=time.monotonic();s=np.load(OUT/'inventory.npz');grid=np.geomspace(1e-7,1e3,points)
    velocities=[];weights=[];errors=[]
    for x in grid:
        v,w=modes(x,order,s);idx=np.argsort(v);v=v[idx];w=w[idx]/2
        exact=(x*x+s['ion_strength'].sum())/(x*x+1+s['ion_strength'].sum())
        errors.append([float(abs(2*w.sum()/exact-1)),float(abs(2*w@(v*v)-1))])
        velocities.append(v);weights.append(w)
    velocities=np.array(velocities);weights=np.array(weights)
    assert np.max(errors)<1e-9,('spectral identities',np.max(errors))
    table=dict(grid=grid,velocities=velocities,weights=weights,
        cumulative_mass=np.c_[np.zeros(points),np.cumsum(weights,axis=1)],
        cumulative_first=np.c_[np.zeros(points),np.cumsum(weights*velocities,axis=1)])
    np.savez_compressed(target,**table)
    old.ex.write(target.with_suffix('.json'),dict(classification='Counterexample candidate',passed=True,
        velocity_order=order,points=points,seconds=time.monotonic()-start,
        static_and_second_moment_max=np.max(errors,axis=0).tolist(),
        formula='Gibbs-whitened Vlasov oscillators: B=H^(1/2) D^2 H^(1/2), H=I+b b^T. SVD avoids forming the squared spectrum; weights are squared electron-density projections. No fitted normalization.'))
    return table


def overlap(table,node,a,b,c,d):
    """Integral of a trapezoid against the positive half of a spectral measure."""
    if table is None:
        cdf=lambda x:ndtr(x)-.5
        mom=lambda x:(1-np.exp(-x*x/2))/np.sqrt(2*np.pi)
        A,B,C,D=[cdf(x) for x in [a,b,c,d]];MA,MB,MC,MD=[mom(x) for x in [a,b,c,d]]
    else:
        v=table['velocities'][node];mass=table['cumulative_mass'][node];first=table['cumulative_first'][node]
        ia,ib,ic,id=[np.searchsorted(v,x,side='right') for x in [a,b,c,d]]
        A,B,C,D=[mass[x] for x in [ia,ib,ic,id]];MA,MB,MC,MD=[first[x] for x in [ia,ib,ic,id]]
    return (MB-MA)-a*(B-A)+(b-a)*(C-B)+d*(D-C)-(MD-MC)


def cell_density(table,x,scale,difference,hi,hj):
    """Exact rectangular cell convolution of each discrete spectral atom.

    q is fixed at the two cell representatives. Only the slow geometry and
    equilibrium prefactors use those representatives; no point sampling of poles.
    """
    total=hi+hj;flat=abs(hi-hj)
    a=(difference-total)/scale;b=(difference-flat)/scale
    c=(difference+flat)/scale;d=(difference+total)/scale
    same=difference==0
    # For the diagonal use twice the positive-half triangle. Off-diagonal
    # disjoint source cells have no negative-frequency overlap.
    a=np.maximum(a,0);b=np.maximum(b,0)
    if table is None:values=overlap(None,0,a,b,c,d)
    else:
        lg=np.log(table['grid']);q=np.log(x);node=np.searchsorted(lg,q)-1
        assert np.all((node>=0)&(node<len(lg)-1)),('q table range',float(x.min()),float(x.max()))
        fraction=(q-lg[node])/(lg[node+1]-lg[node]);values=np.zeros_like(x)
        for n in np.unique(node):
            ids=node==n
            for shift,weight in [(0,1-fraction[ids]),(1,fraction[ids])]:
                values[ids]+=weight*overlap(table,int(n+shift),a[ids],b[ids],c[ids],d[ids])
    # Roundoff in cumulative first moments can leave negative ulps in an empty
    # interval. Reject anything beyond a floating-point cancellation allowance.
    allowance=64*np.finfo(float).eps*(1+d)
    assert np.all(values>=-allowance),float(np.min(values/allowance))
    values=np.maximum(values,0);values[same]*=2
    return values*scale/(4*hi*hj)


def coefficients(bank,order,nangle,nfreq=1):
    assert nfreq==1,'Cell overlap is integrated explicitly; use spectrum/grid comparisons'
    start=time.monotonic();u=bank['u'][bank['u']<=60];C=bank['Ci'][:len(u)];N=len(u)
    s=np.load(OUT/'inventory.npz');theta=float(bank['theta']);width=12*np.sqrt(theta)
    # Include the collective plasma-frequency displacement in the pair window.
    stop=np.searchsorted(u,u*np.exp(width)+4*float(s['plasma_u']),side='right')
    counts=stop-np.arange(N);i=np.repeat(np.arange(N,dtype=np.int32),counts)
    j=np.concatenate([np.arange(a,b,dtype=np.int32) for a,b in enumerate(stop)])
    table=None if bank.get('collective_free',False) else spectrum(int(bank['velocity_order']),int(bank['q_points']))
    B=15*float(bank['arad'])*float(bank['T'])**3/np.pi**4
    occupation=1/np.expm1(u);measure=C/(B*u**4*occupation*(1+occupation))
    half=np.diff(bank['edges_u'][:N+1])/2
    nodes,weights=leggauss(nangle);t=(nodes+1)/2;mu=1-2*t*t;aw=2*weights*t
    P=np.polynomial.legendre.legvander(mu,order-1);alpha=np.zeros((len(i),order))
    for begin in range(0,len(i),40000):
        end=min(begin+40000,len(i));ii=i[begin:end];jj=j[begin:end]
        v=u[ii];z=u[jj];delta=z-v
        quantum=np.ones_like(delta);nonzero=delta>0;quantum[nonzero]=delta[nonzero]/(2*np.sinh(delta[nonzero]/2))
        common=B*float(bank['rate_e'])*measure[ii]*measure[jj]*v*z*quantum*np.sqrt(
            occupation[ii]*(1+occupation[ii])*occupation[jj]*(1+occupation[jj]))
        for cosine,weight,pl in zip(mu,aw,P):
            g=np.sqrt(delta*delta+2*v*z*(1-cosine));scale=g*np.sqrt(theta)
            density=cell_density(table,g*float(s['wave_lambda']),scale,delta,half[ii],half[jj])
            val=common*(3/8*(1+cosine*cosine))*weight*density
            alpha[begin:end]+=val[:,None]*pl
    assert np.all(alpha[:,0]>=0) and np.max(abs(alpha)-alpha[:,0,None])<1e-12
    return i,j,alpha,dict(seconds=time.monotonic()-start,pairs=len(i),frequency_cells=N,nangle=nangle,
        collective=table is not None,velocity_order=bank.get('velocity_order'),q_points=bank.get('q_points'),
        model='Maxwell collisionless multi-ion RPA plus detailed-balance FDT factor; exact rectangular-cell convolution of positive velocity spectral atoms. Vacuum photon kinematics; no native scalar-rate normalization.')


class Operator(previous.Operator):
    # Reuse the frozen sparse coupled solver without modifying Phase62 source.
    __init__=FunctionType(previous.Operator.__init__.__code__,dict(vars(previous),coefficients=coefficients),
        argdefs=previous.Operator.__init__.__defaults__,closure=previous.Operator.__init__.__closure__)


def pilot():
    assert not (OUT/'pilot.json').exists();signal.alarm(120);start=time.monotonic()
    b=dict(np.load(old.OUT/'bank.npz'));b.update(velocity_order=64,q_points=65)
    op=Operator(b,4,1,16);build=time.monotonic()-start;r=op.evolve(float(b['H']/b['c']),16)
    result=dict(classification='Counterexample candidate',seconds=time.monotonic()-start,
        build_seconds=build,coefficients=op.info,balance=r['balance'],energy_residual=r['energy_equation_residual'],
        entropy_growth=r['entropy_growth'],solver_residual=op.max_solver_error,iterations=op.max_iterations)
    old.ex.write(OUT/'pilot.json',result);print('PILOT',result,flush=True);signal.alarm(0)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','inventory','static','pilot']);globals()[p.parse_args().action]()
