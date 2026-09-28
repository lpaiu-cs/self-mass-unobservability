"""Conservative matter/radius derivatives and a saved-trajectory charge readout.

Counterexample candidate: the readout is a quasistatic scalar projection of a
moving material snapshot, not the radiative charge of the fast pulse itself.
The exterior wall and the original momentum reference subtraction are retained.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from types import SimpleNamespace
import argparse
import json
import time
import numpy as np
import sympy as sp
import def_normalized_charge as q
import def_resolved_scalar_pulse_run as production

s,e,ld=q.s,q.e,q.ld
OUT=e.g.OUT/'def-fluid-charge'


def geometry(star,xi):
    rf=star.rf*(1+np.r_[0,xi]);f=star.fraction
    r=((1-f)*rf[:-1]**3+f*rf[1:]**3)**(ld(1)/3)
    volume=4*np.pi/3*np.diff(rf**3)
    return SimpleNamespace(n=star.n,r=r,rf=rf,volume=volume,fraction=f,beta=star.beta,
        distance=np.r_[ld(1),np.diff(r),rf[-1]-r[-1]],
        scalar_area=4*np.pi*np.r_[ld(0),r[:-1]*r[1:],r[-1]*rf[-1]],
        face_weight=np.zeros(star.n+1))


def readout(g,E,R,trace):
    """Leading scalar charge on instantaneous constrained matter stress/geometry."""
    mf=np.r_[E[0]*0,np.cumsum(e.GRAV*g.volume*E)]
    mass=mf[:-1]+g.fraction*np.diff(mf)
    b=1-2*mass/g.r;bf=np.r_[E[0]*0+1,1-2*mf[1:]/g.rf[1:]]
    H=q.lapse(g,E,R,b,bf,np.zeros(g.n+1));mu=mf[-1]/g.rf[-1]
    G=-bf[-1]*np.log1p(-2*mu)/(2*mu)
    # q.scalar's pressure argument encodes trace only; lapse uses radial stress.
    w,J,d,c,res=q.scalar(g,E,(trace+E)/3,bf,H,mu,G)
    return J/mu,dict(a=1/np.sqrt(b),b=b,bf=bf,m=mass,mf=mf,H=H,w=w,
        scalar_residual=float(np.max(abs(res))/max(abs(J),ld('1e-30'))))


def differential(star,base,eta,theta,xi,velocity):
    """Fixed-composition tangent: B changes, moving faces, T and radial velocity."""
    f=star.fraction;rf=star.rf;r=star.r;V=star.volume
    xf=np.r_[xi[0]*0,xi]
    dlogV=3*np.diff(rf**3*xf)/np.diff(rf**3)
    dr=((1-f)*rf[:-1]**3*xf[:-1]+f*rf[1:]**3*xf[1:])/r**2
    er=base['E']+base['rho']*base['raw'][:,9]
    et=base['rho']*base['raw'][:,10]
    k=base['a']**2/r
    F=eta-dlogV+base['a']**2*base['m']*dr/r**2
    dm=np.zeros(star.n+1,dtype=np.result_type(eta,theta,xi,velocity,base['E']))
    for i in range(star.n):
        gv=e.GRAV*V[i]
        forcing=er[i]*F[i]+et[i]*theta[i]+2*base['Q'][i]*velocity[i]+base['E'][i]*dlogV[i]
        dm[i+1]=((1-gv*er[i]*k[i]*(1-f[i]))*dm[i]+gv*forcing)/(1+gv*er[i]*k[i]*f[i])
    dlnrho=F-k*((1-f)*dm[:-1]+f*dm[1:])
    deps=er*dlnrho+et*theta
    dp=base['P']*(base['raw'][:,5]*dlnrho+base['raw'][:,6]*theta)
    dE=deps+2*base['Q']*velocity;dR=dp+2*base['Q']*velocity
    dtrace=-deps+3*dp
    value,_=readout(geometry(star,xi),base['E']+dE,base['R']+dR,base['trace']+dtrace)
    return value


def native_snapshot(star,template,eta=None,theta=None,xi=None,velocity=None):
    """Preserve declared Eulerian shell baryons, T_J, X, v and Q while re-matching.

    Leading phi_infinity -> 0 projection: finite-field effects are not certified.
    For a stored dynamical state, remove its frame factors before evaluating EOS.
    """
    zero=np.zeros(star.n,dtype=ld)
    eta=zero if eta is None else eta;theta=zero if theta is None else theta
    xi=zero if xi is None else xi;velocity=zero if velocity is None else velocity
    g=geometry(star,xi);delta=template['delta'].copy()
    la=template.get('logA',zero)
    delta[:,0]-=3*la;delta[:,1]+=theta-la;delta[:,2]+=velocity
    # The conserved B has no conformal rest-mass factor; Q transforms as A^4.
    for j in [3,4]:delta[:,j]=(star.base[:,j]+delta[:,j])*np.exp(-4*la)-star.base[:,j]
    target=template['a']*template['rho']*template['W']*star.volume*(1+eta)
    targetT=star.reference['T']*np.exp(delta[:,1])
    rows=[]
    for iteration in range(24):
        raw=star.native_aux(delta,zero);z=star.fluid(delta,zero,raw)
        value,met=readout(g,z['E'],z['R'],z['trace'])
        rho=target/(met['a']*z['W']*g.volume)
        nxt=np.log(rho/star.reference['rho'])
        err=float(abs(nxt-delta[:,0]).max());rows.append(err)
        if err<2e-17:break
        delta[:,0]=nxt
    else:raise RuntimeError(('Native projection cap',rows))
    baryon=float(abs(met['a']*z['rho']*z['W']*g.volume/target-1).max())
    assert baryon<2e-13 and met['scalar_residual']<2e-13
    assert float(abs(z['T']/targetT-1).max())<2e-17
    z.update(met,delta=delta)
    return value,z,dict(baryon_relative_residual=baryon,scalar_residual=met['scalar_residual'],iterations=len(rows))


def symbolic():
    a,rho,W,V,B=sp.symbols('a rho W V B',positive=True)
    assert sp.diff(sp.log(B/(a*W*V)),V)==-1/V
    eps,P,Q,v=sp.symbols('eps P Q v',real=True)
    E=(eps+P*v*v+2*Q*v)/(1-v*v)
    R=(eps*v*v+P+2*Q*v)/(1-v*v)
    assert sp.simplify(-E+R+2*P+eps-3*P)==0
    assert sp.diff(E,v).subs(v,0)==2*Q
    m,r,dm,dr=sp.symbols('m r dm dr')
    h=sp.Symbol('h');lna=-sp.log(1-2*(m+h*dm)/(r+h*dr))/2
    assert sp.simplify(sp.diff(lna,h).subs(h,0)-(dm/r-m*dr/r**2)/(1-2*m/r))==0
    return dict(classification='Proven',passed=True,scope='B=a rho W V; Lorentz/heat trace identity; first velocity stress and moving-radius metric differential. No stationary background or orbital signal proof.')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    files=[Path(__file__),Path(q.__file__),Path(s.__file__),Path(s.old.__file__),Path(s.m.__file__),
        s.OUT/'initial.npz',q.OUT/'zero.npz',q.OUT/'gradient.npz',s.OUT/'production/manifest.json',
        e.g.OUT/'gr-normalized-charge-milestone-manifest.json']
    e.write(OUT/'plan.json',dict(classification='Counterexample candidate',symbolic=symbolic(),
        bindings={p.relative_to(s.ROOT).as_posix():e.digest(p) for p in files},
        claim='Connect shell baryon redistribution, moving radii, radial velocity/heat stress and temperature to a scalar exterior readout. Apply a same-inventory quasistatic projection to all nine saved pulse endpoints, without rerunning their evolution.',
        model='Leading phi_infinity->0 scalar readout; instantaneous radial constraints, fixed stored composition in differential controls; full evolving composition in native endpoint projections. No radiative finite-frequency exterior or stationary orbit.',
        control_steps=[.001,.0003],control_directions=['interior_radius','homology_radius','baryon_redistribution','velocity'],
        gates=dict(native=2e-13,derivative_relative=.005,derivative_absolute=2e-12,minimum_time_order=.8,
            maximum_relative_time_difference=.1,minimum_signal_over_twice_time_difference=10,minimum_charge_response=1e-12),
        decision='Report charge response convergence or its failure without enlarging trajectories. Keep necessary amplification separate from a dynamical gain bound. A frozen snapshot charge is not the radiative charge of the fast pulse.',
        budget=dict(hard_timeout_seconds=180,workers=4,blas_threads=1,maximum_native_snapshots=27,
            maximum_snapshot_iterations=24,new_evolution_steps=0,automatic_expansion=False,
            estimate_seconds=[20,100],basis='Phase41 15 native snapshots took 19.08s; new geometry and moving-matter iteration rate is unmeasured.')))


def run():
    start=time.monotonic();plan=json.loads((OUT/'plan.json').read_text())
    for rel,sha in plan['bindings'].items():assert e.digest(s.ROOT/rel)==sha,rel
    assert not (OUT/'result.json').exists();rows=[]
    with ProcessPoolExecutor(max_workers=4,initializer=s.old.imported.initial.original.worker_init) as pool:
        star,initial=s.initialize(pool,-4,ld('.001'),ld(1));initial['delta']=np.zeros_like(star.base)
        baseline,base,checks=native_snapshot(star,initial)
        np.savez_compressed(OUT/'baseline.npz',**{k:v for k,v in base.items() if isinstance(v,np.ndarray)})
        gradients={};z=np.zeros(star.n,dtype=ld)
        for slot,name in enumerate(['baryon','temperature','radius','velocity']):
            values=[]
            for i in range(star.n):
                args=[np.zeros(star.n,dtype=np.clongdouble) for _ in range(4)];args[slot][i]=1e-24j
                values.append(float(differential(star,base,*args).imag/1e-24))
            gradients[name]=np.asarray(values)
        previous=np.load(q.OUT/'gradient.npz')['gradient']
        assert np.max(abs(previous-gradients['temperature']))<1e-10*abs(previous).sum()
        # Independent radial shapes; no duplicate sign/uniform controls.
        shape=np.sin(np.pi*star.rf[1:]/star.rf[-1]);shape[-1]=0
        eta=np.cos(np.pi*star.r/star.rf[-1]);B=initial['a']*initial['rho']*star.volume
        eta-=np.sum(B*eta)/np.sum(B)
        vel=np.sin(np.pi*star.r/star.rf[-1])
        directions=[('interior_radius',[z,z,shape,z]),('homology_radius',[z,z,np.ones(star.n),z]),
            ('baryon_redistribution',[eta,z,z,z]),('velocity',[z,z,z,vel])]
        for name,args in directions:
            expected=sum(gradients[key]@arg for key,arg in zip(gradients,args));differences=[]
            for h in plan['control_steps']:
                pair=[]
                for sign in [-1,1]:
                    value,state,check=native_snapshot(star,initial,*[sign*ld(h)*arg for arg in args]);pair.append(value)
                    np.savez_compressed(OUT/f'{name}-{h}-{sign}.npz',**{k:v for k,v in state.items() if isinstance(v,np.ndarray)})
                    rows.append(dict(name=f'{name}-{h}-{sign}',charge=float(value),**check))
                differences.append(float((pair[1]-pair[0])/(2*ld(h))))
            error=max(abs(v-expected) for v in differences)
            passed=error<plan['gates']['derivative_absolute']+plan['gates']['derivative_relative']*abs(expected)
            rows.append(dict(name=name,expected=float(expected),differences=differences,error=float(error),passed=bool(passed)))
        endpoints=[]
        for n in [48,96,192]:
            for label in ['driven','undriven','decoupled']:
                path=production.OUT/f'{label}-{n}'/f'step-{n:03d}.npz'
                template=dict(np.load(path))
                value,state,check=native_snapshot(star,template)
                # Conserved cumulative baryons identify material-surface motion
                # on the fixed Eulerian grid. No outer-wall expansion is inferred.
                masses=template['B']*star.volume
                before=np.r_[ld(0),np.cumsum(B)];after=np.r_[ld(0),np.cumsum(masses)]
                radii=np.interp(np.asarray(before,float),np.asarray(after,float),np.asarray(star.rf**3,float))**(1/3)
                displacement=radii-np.asarray(star.rf,float);displacement[[0,-1]]=0
                theta=template['delta'][:,1]-template['logA']
                thermal=float(gradients['temperature']@theta)
                entry=dict(label=label,steps=n,charge=float(value),thermal_only_linear_increment=thermal,
                    maximum_material_face_displacement_over_radius=float(abs(displacement).max()/star.rf[-1]),
                    source_sha256=e.digest(path),**check)
                endpoints.append(entry)
                np.savez_compressed(OUT/f'{label}-{n}.npz',**{k:v for k,v in state.items() if isinstance(v,np.ndarray)},material_face_displacement=displacement)
                print('ENDPOINT',json.dumps(entry),flush=True)
        calls=star.pool.evaluations
    readouts=[]
    for name,other in [('driven_minus_undriven','undriven'),('material_coupling_contrast','decoupled')]:
        values=[];thermal=[]
        for n in [48,96,192]:
            a=next(v for v in endpoints if v['steps']==n and v['label']=='driven')
            b=next(v for v in endpoints if v['steps']==n and v['label']==other)
            values.append(a['charge']-b['charge']);thermal.append(a['thermal_only_linear_increment']-b['thermal_only_linear_increment'])
        diff=np.abs(np.diff(values));amp=abs(values[-1]);order=float(np.log2(diff[0]/diff[1]))
        ratio=float(amp/(2*diff[1]));relative=float(diff[1]/amp)
        passed=amp>=1e-12 and order>=.8 and relative<=.1 and ratio>=10
        readouts.append(dict(name=name,charge_differences=values,thermal_only=thermal,time_differences=diff.tolist(),
            empirical_time_order=order,signal_over_twice_time_difference=ratio,relative_time_difference=relative,passed=bool(passed)))
    np.savez_compressed(OUT/'gradients.npz',**gradients)
    result=dict(classification='Counterexample candidate',controls_passed=all(r.get('passed',True) for r in rows),
        charge_time_gates_passed=all(r['passed'] for r in readouts),baseline_charge=float(baseline),controls=rows,
        gradient_l1={k:float(abs(v).sum()) for k,v in gradients.items()},endpoints=endpoints,readouts=readouts,
        native_calls=calls,seconds=time.monotonic()-start,new_evolution_steps=0,physical_orbit_solved=False,
        free_surface_evolved=False,fluid_gain_bounded=False,radiative_exterior_solved=False,symbolic=symbolic())
    e.write(OUT/'result.json',result)
    e.write(OUT/'manifest.json',dict(sha256={p.relative_to(s.ROOT).as_posix():e.digest(p) for p in OUT.iterdir() if p.is_file() and p.name not in ['manifest.json','run.log']}))
    print('FINAL',json.dumps({k:v for k,v in result.items() if k not in ['controls','endpoints']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
