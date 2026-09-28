"""Bounded, reproducible closure of Phase130 artifact and local pairing checks."""
from pathlib import Path
import inspect, json, signal, time
import numpy as np
import sympy as s
import def_native_refined_thermochemistry as run

start=time.monotonic(); signal.signal(signal.SIGALRM,run.flow.old.optical.timeout); signal.alarm(40)
sha=run.sha; read=lambda p:json.loads(p.read_text())
plans=[run.failed.OUT/'plan.json',run.OUT/'plan.json',run.OUT/'execution-plan.json',run.GR/'plan.json']
bindings=[]
for plan in plans:
    for name,digest in read(plan)['bindings'].items():
        p=Path(name)
        if plan==run.OUT/'plan.json' and p.name=='def_native_refined_thermochemistry.py':
            p=run.OUT/'bank-producer.py'
        assert sha(p)==digest,(str(plan),str(p))
        bindings.append(dict(plan=str(plan),file=str(p),sha256=digest))
assert sha(run.failed.OUT/'failed-audit-producer.py')==sha(run.failed.__file__)
assert not read(run.failed.OUT/'result.json')['passed']
for p in [run.OUT/'controls/result.json',run.EV/'result.json',run.GR/'result.json',run.GR/'native-endpoint.json']:
    assert read(p)['passed'],str(p)
old=np.load(run.flow.previous.OUT/'thermal-support/bank.npz');new=np.load(run.OUT/'bank.npz')
for i,dt in enumerate(old['offsets']):
    j=np.flatnonzero(new['offsets']==dt); assert len(j)==1
    for key in ['raw','rates']: assert np.array_equal(old[key][:,:,i],new[key][:,:,j[0]])
assert np.array_equal(old['ratios'],new['ratios'])

a,e,n,da,de,dn,I,K,B,A,E,S=s.symbols('a e n da de dn I K B A E S')
def F(a,e,n):return K*(e-(a-e)*I)+B*(e*E*(1+I)-a*A*I)+n*S
identity=s.expand(F(a+da,e+de,n+dn)-F(a,e,n)-F(da,de,dn))
assert identity==0
run.write(run.OUT/'symbolic-affine-collision.json',dict(classification='Proven',passed=True,
    identity=str(identity),scope='At fixed state, spectral slopes and angular interpolation, collision is linear in absorption, emission and electron density. Signed differences are algebraic inputs, not physical rates.'))

model=run.Coupled();b=model.bulk;d=b.d;rows=[]
for steps in [64,128]:
    p=run.EV/f'final-{steps}.npz';assert sha(p)==sha(run.EV/f'checkpoint-{steps}.npz')
    z=np.load(p);assert int(z['completed'])==steps
    model.Pi=z['Pi'];model.h=z['h'];model.j=z['j']
    model.mass=model.mass0-model.h[1:]+model.h[:-1];model.set_material(run.flow.old.END)
    en=b.volume[:,None]*d['num']*d['Einf']*b.scale/d['a'][:,None]**3
    num=b.volume[:,None]*d['num']*b.scale/d['a'][:,None]**3
    paired_E=b.gas_weight[:,None]*b.en_scaled
    paired_N=(model.mass*d['thermo'][:,4]*d['y0'])[:,None]*b.num_scaled
    er=float(np.max(abs(paired_E/en-1)));nr=float(np.max(abs(paired_N/num-1)))
    assert max(er,nr)<1e-10,(steps,er,nr)
    rows.append(dict(steps=steps,paired_energy_relative=er,paired_H_relative=nr,final_sha256=sha(p)))
owner=inspect.getsource(b.implicit)
assert 'if err<2e-11:break' in owner and "else:raise AssertionError(('Unsplit Newton cap',err))" in owner
modules={}
import sys
for module in list(sys.modules.values()):
    name=getattr(module,'__file__',None)
    if name and Path(name).suffix=='.py' and Path(name).resolve().is_relative_to(Path.cwd()/'verification'):
        modules[str(Path(name).resolve().relative_to(Path.cwd()))]=sha(name)
result=dict(classification='Counterexample candidate',passed=True,bindings_checked=len(bindings),bindings=bindings,
    preserved_old_native_states=190,old_failed_verdict_preserved=True,paired_source_rows=rows,
    symbolic_linear_collision_passed=True,implicit_residual_gate_code_enforced=2e-11,
    measured_implicit_maximum_residual=None,global_H_port_ledger_certified=False,
    verification_module_sha256=modules,seconds=time.monotonic()-start)
run.write(run.OUT/'final-audit.json',result);signal.alarm(0)
print(json.dumps({k:v for k,v in result.items() if k not in ['bindings','verification_module_sha256']}))
