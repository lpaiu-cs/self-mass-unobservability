"""Conservative null-ray transport of the saved native neutrino emissivity.

The emitter is the same free-surface background. No moment closure, diffusion
speed or immediate escape is substituted for null propagation. The background
metric is frozen; this is a radiation/source operator for subsequent coupling,
not a self-gravitating stellar evolution.
"""
from pathlib import Path
import argparse
import json
import time
import numpy as np
import sympy as sp
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss
import def_free_surface_response_normalized as surface
import def_free_surface_thermal as thermal

h=surface.h
OUT=thermal.OUT.parent/'def-causal-neutrinos'


def symbolic():
    r,mu,c=sp.symbols('r mu c',positive=True)
    N=sp.Function('N')(r);a=sp.Function('a')(r)
    nr=sp.diff(N,r)/N
    weight=a*r*r/N**3
    radial=c*r*r/N**2*mu
    angular=c*r*r/N**2*(1-mu*mu)*(1/r-nr)
    assert sp.simplify(sp.diff(radial,r)+sp.diff(angular,mu))==0
    assert sp.simplify(sp.diff(r*r/N**2,r)/2-r*r/N**2*(1/r-nr))==0
    qJ,A,rho,dm,eps=sp.symbols('qJ A rho dm eps',positive=True)
    assert sp.simplify((4*sp.pi*a*r*r*A**5*N**2*rho*eps)/(4*sp.pi*a*r*r*A**3*rho)-A*A*N*N*eps)==0
    return dict(classification='Proven',passed=True,
        invariant='W=N^4 I_E; conserved redshifted energy measure proportional to a*r^2/N^3 dr dmu dW.',
        characteristics='dr/dt=c*N*mu/a; dmu/dt=c*N/a*(1-mu^2)*(1/r-dlnN/dr); impact b=r*sqrt(1-mu^2)/N is constant.',
        emission='Einstein proper-volume emissivity q_E=A^5*rho_J*epsilon_nu. Redshifted energy injection is dm*A^2*N^2*epsilon_nu per coordinate second.',
        source_cancellation='Use the identical emitted-energy increment as minus matter and plus radiation; the sum vanishes cell by cell. Photons already in the LTE EOS are not added again.',
        ray_inventory='For cumulative emission H, segment energy is H(t-tau_enter)-H(t-tau_exit); escaped energy is H(t-tau_escape). Adjacent segments telescope for arbitrary causal source histories.',
        limitation='Collisionless massless neutrinos in a prescribed static metric. Finite neutrino mass, collisions, photon transfer, fluid response and time-dependent metric remain separate.')


class Geometry:
    def __init__(self,flat=False):
        d=np.load(surface.OUT/'background.npz');self.R=float(d['R']*d['r'][-1]);self.tc=self.R/h.gr.C
        self.r=d['r']/d['r'][-1];self.N=d['N'];self.m=d['m']/d['r'][-1]
        self.flat=flat
        # Exact exterior Just background follows the existing vacuum solver.
        if not flat:
            mu=self.m[-1];flux=float(d['r'][-1]*d['v'][-1])
            _,_,self.exterior=surface.old.exterior.background(mu,flux,2e-12)

    def metric(self,r):
        r=np.asarray(r,float)
        if self.flat:return np.ones_like(r),np.ones_like(r)
        N=np.interp(r,self.r,self.N);m=np.interp(r,self.r,self.m)
        central=r<self.r[1];m=np.where(central,self.m[1]*(r/self.r[1])**3,m)
        outside=r>1
        if np.any(outside):
            values=self.exterior.sol(1/r[outside]);m=np.asarray(m).copy();N=np.asarray(N).copy()
            m[outside]=values[0];N[outside]=values[1]/np.sqrt(1-2*m[outside]/r[outside])
        b=1-2*np.divide(m,r,out=np.zeros_like(r),where=r>0)
        return N,1/np.sqrt(b)


def emitters(n,flat=False):
    edges=np.linspace(0,1,n+1)
    if flat:
        # Exact uniform-volume shell weights and the volume centroid.
        power=np.diff(edges**3)
        radius=.75*np.diff(edges**4)/power
        return edges,radius,power
    d=np.load(thermal.OUT/'coefficients.npz');s=np.load(thermal.OUT/'sources.npz')
    geom=Geometry();r=d['radius_cm']/(100*geom.R)
    power=d['dm']*d['A']**2*d['N']**2*(s['neutrino']+s['thermal_neutrino'])
    assert np.min(power)>=0 and np.max(r)<1
    bin_id=np.minimum(np.searchsorted(edges,r,side='right')-1,n-1)
    sums=np.bincount(bin_id,weights=power,minlength=n)
    centres=np.bincount(bin_id,weights=power*r,minlength=n)
    centres=np.divide(centres,sums,out=(edges[:-1]+edges[1:])/2,where=sums>0)
    assert abs(sums.sum()/power.sum()-1)<1e-14
    return edges,centres,sums


def build(n,angles,order=16,flat=False):
    start=time.monotonic();geometry=Geometry(flat)
    _,radius,power=emitters(n,flat);edges=np.linspace(0,2,2*n+1)
    nodes,weights=leggauss(angles);mu=nodes[nodes>0];weights=weights[nodes>0]/2
    gx,gw=leggauss(order)
    segments=[];surface_crossings=[];escape=[];minimum_margin=1.
    for emitter,(re,P) in enumerate(zip(radius,power)):
        if P==0:continue
        Ne=float(geometry.metric(np.array([re]))[0][0])
        for cosine,angular_weight in zip(mu,weights):
            impact=re*np.sqrt(1-cosine*cosine)/Ne
            turning=brentq(lambda r:r/float(geometry.metric(np.array([r]))[0][0])-impact,0,re,xtol=1e-15)
            points=np.unique(np.r_[turning,re,edges[edges>turning]])
            z=np.sqrt(np.maximum(points-turning,0));lo,hi=z[:-1],z[1:]
            zz=(lo[:,None]+hi[:,None])/2+(hi-lo)[:,None]*gx/2
            rr=turning+zz*zz;N,a=geometry.metric(rr)
            D=-np.expm1(2*(np.log(impact)+np.log(N)-np.log(rr)))
            assert D.min()>0,(emitter,cosine,D.min())
            duration=(hi-lo)/2*np.sum(gw*2*zz*a/(N*np.sqrt(D)),axis=1)
            travel=np.r_[0,np.cumsum(duration)]
            ie=int(np.flatnonzero(points==re)[0]);Ts=float(np.interp(1.,points,travel));Te=travel[ie]
            rmid=(points[:-1]+points[1:])/2;Nm,_=geometry.metric(rmid)
            abs_mu=np.sqrt(np.maximum(0,1-(impact*Nm/rmid)**2))
            bins=np.minimum(np.searchsorted(edges,rmid,side='right')-1,len(edges)-2)
            # Both signs of isotropic emission; inward rays turn before leaving.
            for sign in [1,-1]:
                if sign==1:
                    selected=np.arange(ie,len(duration));tin=travel[selected]-Te;tout=travel[selected+1]-Te
                    for j,t0,t1 in zip(selected,tin,tout):segments.append((bins[j],emitter,angular_weight,t0,t1,abs_mu[j]))
                    tau_surface=Ts-Te;tau_outer=travel[-1]-Te
                else:
                    for j in range(ie):segments.append((bins[j],emitter,angular_weight,Te-travel[j+1],Te-travel[j],-abs_mu[j]))
                    for j in range(len(duration)):segments.append((bins[j],emitter,angular_weight,Te+travel[j],Te+travel[j+1],abs_mu[j]))
                    tau_surface=Ts+Te;tau_outer=travel[-1]+Te
                minimum_margin=min(minimum_margin,tau_surface-(1-re))
                surface_crossings.append((emitter,angular_weight,tau_surface))
                escape.append((emitter,angular_weight,tau_outer))
    assert minimum_margin>-1e-10
    result=dict(edges=edges,emitter_radius=radius,emitter_power=power,segments=np.array(segments),
        surface_crossings=np.array(surface_crossings),outer_crossings=np.array(escape),R=geometry.R,tc=geometry.tc)
    return result,dict(seconds=time.monotonic()-start,segments=len(segments),rays=len(escape),minimum_radial_light_cone_margin=minimum_margin)


def cumulative(t,history):
    t=np.maximum(t,0)
    if history=='constant':return t
    # Compact triangular source pulse; tests history convolution and turnarounds.
    return np.where(t<1,t*t/2,np.where(t<2,2*t-t*t/2-1,1.))


def rate(t,history):
    if history=='constant':return (t>=0).astype(float)
    return np.where((t>=0)&(t<1),t,np.where((t>=1)&(t<2),2-t,0.))


def moments(ray,times,history='constant'):
    seg=ray['segments'];bins=seg[:,0].astype(int);ids=seg[:,1].astype(int)
    power=ray['emitter_power'];P=power.sum();w=power[ids]*seg[:,2]
    results=[];energies=[];stresses=[]
    for t in times:
        energy=w*(cumulative(t-seg[:,3],history)-cumulative(t-seg[:,4],history))
        E=np.bincount(bins,weights=energy,minlength=len(ray['edges'])-1)
        F=np.bincount(bins,weights=energy*seg[:,5],minlength=len(E))
        Pr=np.bincount(bins,weights=energy*seg[:,5]**2,minlength=len(E))
        boundary=[]
        for name in ['surface_crossings','outer_crossings']:
            q=ray[name];pp=power[q[:,0].astype(int)]*q[:,1]
            boundary.append((float(pp@cumulative(t-q[:,2],history)),float(pp@rate(t-q[:,2],history))))
        emitted=P*float(cumulative(t,history));escaped,L=boundary[1]
        defect=abs(E.sum()+escaped-emitted)/max(P,emitted,1e-100)
        cone=max(float(np.max(abs(F)-E)),float(np.max(Pr-E)),float(-np.min(Pr)))/max(P,1e-100)
        trace=float(np.max(abs(-E+Pr+2*((E-Pr)/2))))/max(P,1e-100)
        assert defect<2e-13 and cone<2e-14 and trace<2e-14,(defect,cone,trace)
        # This debit is precisely the energy credited to radiation at emission,
        # not the delayed surface/outer flux; material EOS evolution is separate.
        matter_debit=-power*float(cumulative(t,history))
        results.append(dict(time_crossings=float(t),time_seconds=float(t*ray['tc']),
            emitted_redshifted_erg=float(emitted*ray['tc']),matter_source_debit_erg=float(matter_debit.sum()*ray['tc']),
            radiation_inventory_redshifted_erg=float(E.sum()*ray['tc']),outer_escaped_erg=float(escaped*ray['tc']),
            surface_escaped_erg=float(boundary[0][0]*ray['tc']),surface_luminosity_erg_s=boundary[0][1],
            outer_luminosity_erg_s=L,energy_balance_relative=defect,trace_relative=trace,cone_excess=cone))
        energies.append(E*ray['tc']);stresses.append(np.array([E,F,Pr])*ray['tc'])
    return results,np.array(energies),np.array(stresses)


def prepare():
    assert not OUT.exists();OUT.mkdir()
    paths=[Path(__file__),surface.OUT/'background.npz',thermal.OUT/'coefficients.npz',thermal.OUT/'sources.npz',
        h.ROOT/'verification/gr_radiation_metric_equations.py',h.ROOT/'verification/def_radiative_exterior.py']
    h.write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='08c631fe',
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths},symbolic=symbolic(),
        claim='Complete a finite-time causal, angle-resolved neutrino energy/stress operator on the same free-surface background, with shared matter debit and radiation credit. Remove instantaneous escape from the dynamic mass budget.',
        source='Actual Phase50 native nuclear plus thermal neutrino emissivity, isotropic in the initially resting fluid. Constant initial emissivity is a frozen-source response, not a nonlinear reaction evolution.',
        approximation='Collisionless massless neutrinos, prescribed initial GR/DEF metric, no incoming radiation and initially empty neutrino field. Retain the original LTE photons in the EOS; this operator adds neutrinos only.',
        controls='Exact flat uniform-sphere flight-distance distribution; compact triangular source history; positive energy, massless trace, moment cone, radial light-cone bound and total emitted=in-flight+escaped energy.',
        grid_pairs=[[48,24],[96,48]],quadrature_order=16,quadrature_control_order=24,
        times_crossings=np.linspace(0,3.5,71).tolist(),
        gates=dict(balance_relative=2e-13,flat_escaped_energy_relative=.01,flat_flux_time_L1_relative=.02,
                   spatial_escaped_energy_relative=.02,angular_radial_moment_L1_relative=.03,quadrature_escaped_relative=1e-7),
        budget=dict(pilot_cells=12,pilot_angles=12,hard_pilot_seconds=30,production_hard_seconds=180,
                    workers=1,BLAS_threads=1,native_calls=0,GPU=False,estimated_peak_memory_MB=256,maximum_transport_builds=5,automatic_expansion=False),
        scope='Finite-time radiation/source subproblem. No claim of evolved fluid, updated Einstein geometry, physical photon atmosphere, neutrino-opacity certification or full dynamic-charge completion.'))
    ray,measurement=build(12,12)
    # Pessimistic work scaling includes source, angle and radial segment counts.
    measurement['production_forecast_seconds']=measurement['seconds']*((96/12)**2*(48/12)+(48/12)**2*(24/12)+2*(48/12)**2*(24/12)*24/16)*1.5
    measurement['array_bytes']=int(sum(v.nbytes for v in ray.values() if isinstance(v,np.ndarray)))
    h.write(OUT/'pilot.json',measurement);np.savez_compressed(OUT/'pilot.npz',**ray)
    print(json.dumps(measurement),flush=True)


def run():
    assert not (OUT/'result.json').exists();start=time.monotonic()
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,sha in plan['bindings'].items():assert h.digest(h.ROOT/rel)==sha,rel
    pilot=json.loads((OUT/'pilot.json').read_text());assert pilot['production_forecast_seconds']<165,pilot
    times=np.array(plan['times_crossings']);cases={};measurements={}
    for label,n,angles,order,flat in [('coarse',48,24,16,False),('fine',96,48,16,False),
                                     ('quadrature',48,24,24,False),('flat',48,24,16,True)]:
        ray,measurement=build(n,angles,order,flat);rows,E,stress=moments(ray,times)
        pulse,pE,pStress=moments(ray,times,'pulse')
        np.savez_compressed(OUT/(label+'.npz'),**ray,times=times,energy=E,stress_inventory=stress,pulse_energy=pE,pulse_stress_inventory=pStress)
        h.write(OUT/(label+'.json'),dict(classification='Counterexample candidate',measurement=measurement,constant=rows,pulse=pulse))
        cases[label]=(rows,E,stress);measurements[label]=measurement
        print(label,measurement,flush=True);assert time.monotonic()-start<170
    flat=cases['flat'][0];x=times;analytic=np.where(x<2,3*x*x/8-x**4/64,x-.75)
    flat_ray=np.load(OUT/'flat.npz');got=np.array([r['surface_escaped_erg'] for r in flat])/float(flat_ray['tc'])
    flat_energy=float(np.max(abs(got-analytic))/np.max(analytic))
    expected_flux=np.where(x<2,3*x/4-x**3/16,1.)
    flat_flux=float(np.trapezoid(abs(np.array([r['surface_luminosity_erg_s'] for r in flat])-expected_flux),x)/np.trapezoid(expected_flux,x))
    def boundary_difference(a,b):
        aa=np.array([r['surface_escaped_erg'] for r in cases[a][0]]);bb=np.array([r['surface_escaped_erg'] for r in cases[b][0]])
        return float(np.max(abs(aa-bb))/max(abs(bb).max(),1e-100))
    space=boundary_difference('coarse','fine');quad=boundary_difference('quadrature','coarse')
    fine_stress=cases['fine'][2].reshape(len(times),3,96,2).sum(axis=-1)
    coarse_stress=cases['coarse'][2]
    moment=float(np.max(np.sum(abs(fine_stress-coarse_stress),axis=-1))/np.max(np.sum(abs(fine_stress),axis=-1)))
    state=dict(classification='Counterexample candidate',finite_time_neutrino_transport=True,
        flat_energy_relative=flat_energy,flat_flux_time_L1_relative=flat_flux,
        spatial_surface_escaped_relative=space,quadrature_surface_escaped_relative=quad,moment_inventory_relative=moment,
        passed=flat_energy<.01 and flat_flux<.02 and space<.02 and quad<1e-7 and moment<.03,
        seconds=time.monotonic()-start,measurements=measurements,
        radiation_energy_momentum_available=True,matter_EOS_time_evolved=False,Einstein_metric_time_evolved=False,
        physical_neutrino_opacity_certified=False,physical_photon_atmosphere=False,full_dynamic_charge_solved=False)
    h.write(OUT/'result.json',state);print(json.dumps(state),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
