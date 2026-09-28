"""Replay saved RK mass additions and retain small shell increments separately."""
import json,sys
from fractions import Fraction as F
from types import FunctionType,MethodType,SimpleNamespace
import numpy as np
import gr_direct_cell_integrals as cells

g=cells.g;OUT=g.OUT/'gr-mass-increment-replay'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def rational(x):
    return dict(numerator=str(x.numerator),denominator=str(x.denominator),display=float(x))


def prepare():
    assert not OUT.exists();OUT.mkdir();cells.verify()
    paths=[g.ROOT/'verification/gr_mass_increment_replay.py',g.ROOT/'verification/common_eos.py',
        g.ROOT/'verification/baryon_entropy.py',g.ROOT/'verification/audit_structured_enthalpy.py',
        g.OUT/'initial-state-17-4.npz',g.OUT/'initial-structure-17-4.json',g.OUT/'reference-state.npz',
        g.OUT/'initial-adiabats-17.npz',cells.OUT/'manifest.json']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='f161e59',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},cells=[0,1,2],
        baseline='Original RK4 state/RHS and subdivision 4, with each proposed and stored mass addition recorded. Compare against original step and saved face products bitwise.',
        repair='Same RK4 formula and binary64 RHS evaluation, accumulate small three-coordinate shifts in longdouble from a fixed cell starting point. Retain shell mass increments without subtracting two total masses. No native EOS precision increase is claimed.',
        repair_subdivisions=[4,8],direct_mass_relative_tolerance=1e-5,
        limits='Post-diagnostic arithmetic repair on three outer cells. Original source, saved full-star state and time driver remain unchanged. Finite refinement and local mass comparisons do not certify full-star geometry, EOS physics, continuum errors or time evolution.'))


def inputs():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    return plan,dict(np.load(g.OUT/'initial-state-17-4.npz'))


def step(solver,a,b,start,i,shifted):
    n=max(solver.sub,int(np.ceil(abs(b-a)/(.12/solver.sub))));h=(b-a)/n
    base=start.copy();delta=np.zeros(3,dtype=np.longdouble);y=start.copy();records=[]
    for k in range(n):
        x=a+k*h
        if shifted:
            # ponytail: native RHS remains binary64; small increments are the
            # retained state. Certified RHS/ODE errors are a separate task.
            y=base.astype(np.longdouble)+delta
            rhs=lambda xx,yy:solver.rhs(xx,np.asarray(yy,dtype=float),i,solver.mat.B,True)
        else:rhs=lambda xx,yy:solver.rhs(xx,yy,i,solver.mat.B,True)
        k1=rhs(x,y);k2=rhs(x+h/2,y+h*k1/2)
        k3=rhs(x+h/2,y+h*k2/2);k4=rhs(x+h,y+h*k3)
        dy=h*(k1+2*k2+2*k3+k4)/6
        if shifted:delta+=dy
        else:
            before=float(y[1]);y=y+dy;records.append([before,float(dy[1]),float(y[1])])
    return (np.asarray(base.astype(np.longdouble)+delta,dtype=float) if shifted else y),delta,np.array(records)


def arithmetic(trace):
    increments=sum((F(float(r[1])) for r in trace),F(0))
    actual=F(float(trace[-1,2]))-F(float(trace[0,0]))
    errors=sum((F(float(after))-F(float(before))-F(float(dy)) for before,dy,after in trace),F(0))
    assert actual-increments==errors
    assert np.array_equal(trace[1:,0],trace[:-1,2])
    return dict(classification='Proven',scope='Exact rational identity of the saved binary64 mass additions only.',
        steps=len(trace),unchanged_mass_steps=int(np.sum(trace[:,0]==trace[:,2])),
        proposed_below_half_ulp=int(np.sum(abs(trace[:,1])<np.spacing(trace[:,0])/2)),
        proposed_sum=rational(increments),stored_sum=rational(actual),rounding_sum=rational(errors),
        missing_mass_fraction_of_proposed=float(errors/(-increments)))


def run():
    plan,state=inputs();reference=dict(np.load(g.OUT/'reference-state.npz'))
    stats=dict(calls=0,evaluations=0,maximum_score=0.,label='mass-increment-replay')
    inverse=FunctionType(g.audit.strict_invert.__code__,dict(g.audit.strict_invert.__globals__,
        ROOT_STATS=stats,e=SimpleNamespace(OUT=OUT,save=save)))
    solver=object.__new__(g.c.Structure)
    initialize=FunctionType(g.c.Structure.__init__.__code__,dict(g.c.Structure.__init__.__globals__,OUT=g.OUT,EOS=g.EOS))
    initialize(solver,'initial',reference,17,4)
    local_state=FunctionType(g.c.Structure.state.__code__,dict(g.c.Structure.state.__globals__,be=SimpleNamespace(invert=inverse)))
    solver.state=MethodType(local_state,solver);m=solver.mat;B=m.B
    parameters=json.loads((g.OUT/'initial-structure-17-4.json').read_text())['parameters']
    rs=m.R*np.exp(parameters[1]);ms=g.c.gr.TARGET*np.exp(parameters[2])
    p,e,b=solver.state(m.lp[0],0);f=1-2*ms/rs
    rB=np.sqrt(f)/(4*np.pi*rs*rs*b);massB=e/b*np.sqrt(f)
    lpB=-(e+p)*(ms+4*np.pi*rs**3*p)/(4*np.pi*rs**4*b*np.sqrt(f)*p)
    w0=min(m.dm[0]*1e-6,1e-8/abs(lpB*B));seed=massB*B*w0
    start=np.array([(rs-rB*B*w0)/m.R,(ms-seed)/B,m.lp[0]-lpB*B*w0])
    assert ms==state['mass_faces_geom'][0] and rs==state['radius_faces_m'][0]
    records=[];y=start.copy()
    for i in plan['cells']:
        a=np.log(w0 if i==0 else m.outer[i]);b=np.log(m.outer[i+1])
        original=solver.step(a,b,y.copy(),i,B,True)
        replay,_,trace=step(solver,a,b,y,i,False)
        assert np.array_equal(original,replay),(i,'original-step replay')
        assert replay[0]*m.R==state['radius_faces_m'][i+1],(i,'saved radius')
        assert replay[1]*B==state['mass_faces_geom'][i+1],(i,'saved mass')
        np.savez_compressed(OUT/f'cell-{i}-baseline.npz',mass_additions=trace,start=y,end=replay)
        row=arithmetic(trace);row['cell']=i
        integrated=-sum((F(float(dy)) for dy in trace[:,1]),F(0))*F(B)+(F(seed) if i==0 else 0)
        saved=F(float(state['mass_faces_geom'][i]))-F(float(state['mass_faces_geom'][i+1]))
        # Include seed, normalized-state conversion and final product rounding
        # separately; the complete equality is exact for the saved numbers.
        seed_or_start=F(float(state['mass_faces_geom'][i]))-F(float(y[1]))*F(B)-(F(seed) if i==0 else 0)
        step_error=(F(float(y[1]))-F(float(replay[1]))+sum((F(float(dy)) for dy in trace[:,1]),F(0)))*F(B)
        product_error=F(float(replay[1]))*F(B)-F(float(state['mass_faces_geom'][i+1]))
        assert saved-integrated==seed_or_start+step_error+product_error
        row.update(proposed_plus_seed_mass_cm=rational(integrated*100),saved_difference_cm=rational(saved*100),
            seed_or_start_conversion_cm=rational(seed_or_start*100),step_addition_error_cm=rational(step_error*100),
            endpoint_product_error_cm=rational(product_error*100),full_difference_identity_passed=True)
        records.append(row);y=replay;print('MASS ADDITION REPLAY',i,row['missing_mass_fraction_of_proposed'],flush=True)
    direct={row['cell']:row['paths'][-1]['integrated_geometric_mass_cm'] for row in json.loads((cells.OUT/'result.json').read_text())['records']}
    paths=[]
    for sub in plan['repair_subdivisions']:
        solver.sub=sub;y=start.copy();rows=[]
        for i in plan['cells']:
            a=np.log(w0 if i==0 else m.outer[i]);b=np.log(m.outer[i+1])
            y,delta,_=step(solver,a,b,y,i,True)
            shell=float(-delta[1]*np.longdouble(B)*100+(np.longdouble(seed)*100 if i==0 else 0))
            relative=shell/direct[i]-1
            row=dict(cell=i,shell_geometric_mass_cm=shell,direct_EOS_mass_relative_difference=relative,
                direct_mass_passed=bool(abs(relative)<plan['direct_mass_relative_tolerance']),
                radius_endpoint_difference_m=float(y[0]*m.R-state['radius_faces_m'][i+1]))
            rows.append(row);np.savez_compressed(OUT/f'cell-{i}-shifted-{sub}.npz',end=y,shift=delta)
            print('SHIFTED MASS',sub,row,flush=True)
        paths.append(dict(subdivision=sub,rows=rows))
    # Positive control: repeatedly adding a sub-ULP mass change to unity loses
    # every update; storing the shift resolves the declared cumulative change.
    plain=1.;shift=np.longdouble(0);tiny=1e-18
    for _ in range(1000):plain-=tiny;shift-=tiny
    assert plain==1. and abs(float(shift)/(1000*tiny)+1)<1e-14
    save('result.json',dict(classification='Counterexample candidate',completed=True,baseline_records=records,
        shifted_paths=paths,entropy_root_statistics=stats,positive_control_passed=True,
        longdouble_mantissa_bits=np.finfo(np.longdouble).nmant,
        repair_direct_mass_passed=all(r['direct_mass_passed'] for path in paths for r in path['rows']),
        full_star_recomputed=False,continuous_error_certified=False,full_GR_evolution=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify()


def verify():
    plan,_=inputs()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    result=json.loads((OUT/'result.json').read_text());assert result['completed'] and result['positive_control_passed']
    for i,row in zip(plan['cells'],result['baseline_records']):
        recomputed=arithmetic(np.load(OUT/f'cell-{i}-baseline.npz')['mass_additions'])
        assert all(row[k]==v for k,v in recomputed.items())
    print('PASS exact saved RK mass-addition replay; finite shifted-cell comparison only',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
