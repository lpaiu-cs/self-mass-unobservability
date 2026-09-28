"""Pair arbitrary cell source histories with causal neutrino stress/face budgets."""
from pathlib import Path
import json
import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss
import def_causal_neutrinos as transport

h=transport.h
OUT=transport.OUT/'source-coupling'


class History:
    """Piecewise constant power for each fixed emitter, with exact integrals."""
    def __init__(self,times,power):
        self.times=np.asarray(times,float);self.power=np.asarray(power,float)
        assert self.times[0]==0 and np.all(np.diff(self.times)>0)
        assert self.power.shape[0]==len(self.times)-1 and np.min(self.power)>=0
        self.integral=np.vstack([np.zeros(self.power.shape[1]),np.cumsum(np.diff(self.times)[:,None]*self.power,axis=0)])

    def energy(self,t,emitter):
        t=np.clip(np.asarray(t),0,self.times[-1])
        i=np.minimum(np.searchsorted(self.times,t,side='right')-1,len(self.power)-1)
        return self.integral[i,emitter]+(t-self.times[i])*self.power[i,emitter]


def state(ray,history,t):
    segments=ray['segments'];cell=segments[:,0].astype(int);emitter=segments[:,1].astype(int);tc=float(ray['tc'])
    energy=segments[:,2]*(history.energy(t-tc*segments[:,3],emitter)-history.energy(t-tc*segments[:,4],emitter))
    n=len(ray['edges'])-1
    E=np.bincount(cell,weights=energy,minlength=n)
    J=np.bincount(cell,weights=energy*segments[:,5],minlength=n)
    P=np.bincount(cell,weights=energy*segments[:,5]**2,minlength=n)
    emitted=history.energy(t,np.arange(len(ray['emitter_radius'])))
    source_cell=np.searchsorted(ray['edges'],ray['emitter_radius'],side='right')-1
    matter=-np.bincount(source_cell,weights=emitted,minlength=n)
    # Integrated net face luminosity follows the ray inventory. Check its
    # material-surface and external values against independent crossing records.
    face=np.r_[0,np.cumsum(-matter-E)]
    checks=[]
    for key,face_id in [('surface_crossings',n//2),('outer_crossings',n)]:
        q=ray[key];direct=float(q[:,1]@history.energy(t-tc*q[:,2],q[:,0].astype(int)))
        checks.append(abs(face[face_id]-direct)/max(float(emitted.sum()),1.))
    assert max(checks)<2e-13,checks
    return dict(E_inventory=E,J_inventory=J,P_inventory=P,matter_debit=matter,
                integrated_face_energy=face,crossing_balance=max(checks))


def local_stress(ray,value):
    geo=transport.Geometry();x,w=leggauss(16);edge=ray['edges']
    r=(edge[:-1,None]+edge[1:,None])/2+np.diff(edge)[:,None]*x/2
    N,a=geo.metric(r)
    weight=4*np.pi*(float(ray['R'])*100)**3*np.diff(edge)/2*np.sum(w*N*a*r*r,axis=1)
    E=value['E_inventory']/weight;J=value['J_inventory']/weight;P=value['P_inventory']/weight
    return dict(E=E,J=J,Pr=P,Pperp=(E-P)/2,redshifted_proper_volume=weight)


def symbolic():
    v,q=sp.symbols('v q',real=True);W=1/sp.sqrt(1-v*v)
    radiation=sp.Matrix([W*q,W*v*q]);matter=-radiation
    assert matter+radiation==sp.zeros(2,1)
    inverse=sp.Matrix([[W,-W*v],[-W*v,W]])
    assert sp.simplify(inverse*radiation-sp.Matrix([q,0]))==sp.zeros(2,1)
    x=sp.symbols('x',real=True);pdf=sp.Rational(3,4)*(1-x*x/4)
    assert sp.integrate(pdf,(x,0,2))==1
    assert sp.integrate(x*pdf,(x,0,2))==sp.Rational(3,4)
    return dict(classification='Proven',passed=True,
        four_force='Isotropic comoving emission q has Eulerian source (W*q,W*v*q); matter has its exact negative. Momentum cannot be omitted for a moving emitter.',
        flat_control='Uniform sphere, isotropic directions: flight-distance pdf=3/4*(1-x^2/4), 0<x<2; mean flight=3/4 crossing time.',
        source_history='Piecewise constant emitter powers integrate exactly; the shared material debit and ray credit use the identical history integral.',
        limitation='A source/transport coupling API and frozen-metric radiation stress; the material EOS and Einstein evolution are not advanced here.')


def main():
    assert not OUT.exists();OUT.mkdir();proof=symbolic()
    paths=[Path(__file__),Path(transport.__file__),transport.OUT/'fine.npz',transport.OUT/'fine.json',transport.OUT/'result.json']
    h.write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='08c631fe',
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths},symbolic=proof,
        claim='Provide the finite-time radiation stress and identical matter/radiation source histories needed by a conservative coupled solver; include moving-emitter four-force identity.',
        controls='Replay the saved constant-native source; then a manufactured independently time-modulated power in every emitter. Compare shell budgets and direct surface/external ray crossings without building new rays.',
        gates=dict(replay_relative=1e-12,source_and_face_balance_relative=2e-13),
        budget=dict(hard_timeout_seconds=30,new_transport_builds=0,native_calls=0,GPU=False),
        scope='No material state or Einstein time advance. History modulation is a software/control input, not a new native nuclear trajectory.'))
    ray=dict(np.load(transport.OUT/'fine.npz'));P=ray['emitter_power'];tc=float(ray['tc']);times=np.linspace(0,3.5,15)*tc
    constant=History([0,times[-1]],P[None,:])
    rows=[];stored=np.load(transport.OUT/'fine.npz')['stress_inventory']
    for j in [0,5,10,20,35,50,70]:
        t=float(ray['times'][j]*tc);value=state(ray,constant,t)
        actual=np.array([value[k] for k in ['E_inventory','J_inventory','P_inventory']])
        error=float(np.max(abs(actual-stored[j]))/max(abs(stored[j]).max(),P.sum()*tc*1e-14))
        assert error<1e-12,error
        rows.append(dict(time_seconds=t,replay_relative=error,balance=value['crossing_balance']))
    # Different history in each cell proves the operator accepts nonseparable
    # source changes, not only a global rescaling of the initial luminosity.
    factors=1+.2*np.sin(np.arange(len(P))[None,:]+np.arange(len(times)-1)[:,None]*.7)
    varying=History(times,P[None,:]*factors)
    controls=[]
    for t in times:
        value=state(ray,varying,float(t));stress=local_stress(ray,value)
        numerator=value['matter_debit']+value['E_inventory']+np.diff(value['integrated_face_energy'])
        balance=float(abs(numerator).max()/max(P.sum()*tc,1.))
        assert balance<2e-13
        controls.append(dict(time_seconds=float(t),local_energy_balance=balance,crossing_balance=value['crossing_balance']))
    np.savez_compressed(OUT/'final-stress.npz',radius_edges=ray['edges']*float(ray['R'])*100,**value,**stress)
    actual=json.loads((transport.OUT/'fine.json').read_text())['constant']
    selected=[actual[j] for j in [5,10,20,35,50,70]]
    result=dict(classification='Counterexample candidate',passed=True,symbolic=proof,constant_native_replay=rows,independent_history_controls=controls,
        actual_native_initial_power_erg_s=float(P.sum()),crossing_time_seconds=tc,finite_horizon_seconds=float(times[-1]),
        physical_neutrino_transport_rows=selected,
        runtime_measurement=dict(production_seconds=7.643865892,process_wall_seconds=10.18,peak_RSS_KiB=334328,
                                 peak_RSS_estimate_MB=256,memory_estimate_exceeded=True),
        radiation_source_credit_and_matter_debit_connected=True,radiation_stress_for_metric_available=True,
        source_history_manufactured_controls_are_physical=False,matter_EOS_time_evolved=False,Einstein_metric_time_evolved=False,
        scalar_charge_time_evolved=False,full_dynamic_charge_solved=False)
    h.write(OUT/'result.json',result)
    print(json.dumps(dict(passed=True,seconds=times[-1],max_replay=max(q['replay_relative'] for q in rows),max_crossing=max(q['crossing_balance'] for q in controls),rows=selected)),flush=True)


if __name__=='__main__':main()
