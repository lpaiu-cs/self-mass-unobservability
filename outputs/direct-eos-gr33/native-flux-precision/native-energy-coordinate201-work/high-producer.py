"""Read the saved rejected pair; test the native energy-coordinate round trip."""
from pathlib import Path
from types import FunctionType
import json,os,sys,time,numpy as np
from decimal import localcontext
import continue_precise_native as prior
OUT=Path('native-energy-coordinate201-work');OLD=prior.OUT
roundoff=len(sys.argv)>1 and sys.argv[1]=='roundoff'
high=len(sys.argv)>1 and sys.argv[1]=='high'
assert not (OUT/('high-result.json' if high else 'roundoff-result.json' if roundoff else 'result.json')).exists();OUT.mkdir(exist_ok=True)
for s in [0,1]:
    for part in ['photons','material']:(OUT/f'sweep-{s}/{part}').mkdir(parents=True,exist_ok=True)
for src in list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
    dst=OUT/src.relative_to(OLD);dst.parent.mkdir(parents=True,exist_ok=True)
    if not dst.exists():os.link(src,dst)
    assert prior.sha(dst)==prior.sha(src)
prior.joint.previous.original.inf.incident.native.deadline(120)
start=time.monotonic();init=FunctionType(prior.initialize.__code__,dict(prior.initialize.__globals__,OUT=OUT));init(False)
m=prior.owner.Model(64);engine=prior.joint.previous.engine;LD=prior.LD
z=dict(np.load(OLD/'last-accepted-64.npz'));f=dict(np.load(OLD/'failed-linear-64.npz'))
t,h=z['next_time'][()],z['next_step'][()];pairs=[m.unpack(row)[1] for row in f['solution'].reshape(2,-1)]
old=np.array([m.native(t+c*h,g) for c,g in zip(prior.joint.C,pairs)])
if high:
    import precise_native_baryon as native
    states=[];rates=[];timing=[];scale=np.linalg.norm(f['rhs'])
    for digits in [40,70]:
        mark=time.monotonic()
        with native.mp.workdps(digits):
            tangent=native.build(prior.joint,prior.owner)
            values=np.array([native.native_B(m,t+c*h,g,tangent) for c,g in zip(prior.joint.C,pairs)])
            defect=native.hp(np.array(pairs)[:,:,2])-native.hp(z['g'][:,2])-native.hp(h)*(native.hp(prior.joint.A)@values)
            cast=lambda a:np.array([[LD(str(v)) for v in row] for row in a])
            rates.append(cast(values));states.append(cast(defect));timing.append(time.monotonic()-mark)
    result=dict(classification='Counterexample candidate',seconds=time.monotonic()-start,
        digits=[40,70],evaluation_seconds=timing,actual_B=[float(np.linalg.norm(v)/scale) for v in states],
        precision_change=float(np.linalg.norm(states[0]-states[1])/scale),
        original_gate=1e-12,new_physical_steps=0,actual_stage_accepted=False,final_charge_conclusion='unadjudicated')
    np.savez_compressed(OUT/'high-comparison.npz',old_native=old,high_B=rates,high_defect=states)
    prior.write(OUT/'high-result.json',result);print(json.dumps(result,indent=2));sys.exit()
if roundoff:
    guides=[m.unpack(row)[1] for row in f['guess'].reshape(2,-1)];selected=[];maps=[]
    for c,guide,g in zip(prior.joint.C,guides,pairs):
        tangent,reset=prior.joint.selected_tangent();m.native(t+c*h,guide,tangent=tangent)
        reset(False);selected.append(m.native(t+c*h,g,tangent=tangent))
        maps.append(m.jacobian(t+c*h,guide))
    exact=prior.precision.exact;rows=prior.precision.baryon_rows;A=prior.joint.A;scale=np.linalg.norm(f['rhs'])
    measure=lambda a:float(np.linalg.norm(a)/scale)
    with localcontext() as ctx:
        ctx.prec=80
        products=[rows(J,g) for (J,_),g in zip(maps,pairs)]
        anchors=[rows(J,g) for (J,_),g in zip(maps,guides)]
        predicted=np.array([[LD(str(exact(base[k,2])+products[j][k]-anchors[j][k])) for k in range(m.n)] for j,(_,base) in enumerate(maps)])
        precise=np.array([[LD(str(exact(pairs[i][k,2])-exact(z['g'][k,2])-exact(h)*sum(exact(A[i,j])*exact(old[j,k,2]) for j in range(2)))) for k in range(m.n)] for i in range(2)])
    raw=np.array(pairs)[:,:,2]-z['g'][:,2]-h*(A@old[:,:,2])
    selected=np.array(selected)
    result=dict(classification='Counterexample candidate',seconds=time.monotonic()-start,
        actual_B=measure(raw),exact_outer_sum_B=measure(precise),outer_arithmetic_error=measure(raw-precise),
        branch_selected_native_exactly_matches_actual=bool(np.array_equal(selected,old)),
        actual_minus_selected_native_B=measure(h*A@(old[:,:,2]-selected[:,:,2])),
        actual_minus_affine_prediction_B=measure(h*A@(old[:,:,2]-predicted)),
        original_gate=1e-12,new_physical_steps=0,actual_stage_accepted=False,final_charge_conclusion='unadjudicated')
    np.savez_compressed(OUT/'roundoff-comparison.npz',actual_native=old,selected_native=selected,predicted_B=predicted,raw_B=raw,exact_B=precise)
    prior.write(OUT/'roundoff-result.json',result);print(json.dumps(result,indent=2));sys.exit()
mark='et=z.astype(LD).copy();et[2]-=(m.a.astype(LD)-m.model.m.a0)*m.model.cx*C*C*et[0]'
assert engine.source.count(mark)==1
source=engine.source.replace(mark,'et=z.astype(LD).copy();et[2]=m.native_etilde')
ns=dict(engine.tangent.__globals__);exec(compile(source,__file__,'exec'),ns)
new=[];energy=[]
for c,g in zip(prior.joint.C,pairs):
    exact=g[:,0]*m.eu;q=m.conserved(g)
    reconstructed=q[2]-(m.material.a.astype(LD)-m.material.model.m.a0)*m.material.model.cx*engine.C*engine.C*q[0]
    m.material.native_etilde=exact
    new.append(m.native(t+c*h,g,tangent=ns['tangent']))
    error=reconstructed-exact
    energy.append(dict(cell261_direct=str(exact[261]),cell261_roundtrip_error=str(error[261]),
        cell261_relative=float(abs(error[261])/max(abs(exact[261]),LD('1e-290')))))
new=np.array(new);scale=np.linalg.norm(f['rhs']);A=prior.joint.A
gas=np.array(pairs);initial=z['g']
measure=lambda x:float(np.linalg.norm(x)/scale)
oldB=gas[:,:,2]-initial[:,2]-h*(A@old[:,:,2]);newB=gas[:,:,2]-initial[:,2]-h*(A@new[:,:,2])
result=dict(classification='Counterexample candidate',seconds=time.monotonic()-start,
    energy=energy,old_actual_B=measure(oldB),direct_energy_actual_B=measure(newB),
    B_stage_change=measure(h*A@(new[:,:,2]-old[:,:,2])),original_gate=1e-12,
    new_physical_steps=0,actual_stage_accepted=False,final_charge_conclusion='unadjudicated')
np.savez_compressed(OUT/'coordinate-comparison.npz',old_native=old,direct_energy_native=new,old_B=oldB,new_B=newB)
prior.write(OUT/'result.json',result);print(json.dumps(result,indent=2))
