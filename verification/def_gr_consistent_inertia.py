"""Counterexample candidate: repair underintegrated principal wave inertia.

Static K, source, physical grid and time counts stay fixed. The uniform-medium
limit is the positive consistent linear-element mass, not an arbitrary filter.
Full variable-coefficient GR energy and pressure reconstruction remain checks.
"""
from pathlib import Path
from types import SimpleNamespace
import argparse
import json
import signal
import time
import resource
import numpy as np
import sympy as sp
from scipy.sparse import coo_matrix
import def_gr_transient_radau as method

task=method.prior.old.task
OUT=method.OUT/'consistent-inertia';write=method.write;digest=method.digest


def assemble(bg,fn):
    K,D,forcing=task.coupled.assemble(bg,fn)
    ops,B,_,_=task.coupled.operators(bg.mid,fn)
    dx=np.diff(bg.grid);correction=dx[:,None,None]**2/12*np.einsum('nij,njk->nik',ops[:,:,:4],B)
    correction[:,[1,3],:]=0
    rows=[];cols=[];values=[];n=len(dx)
    for shift,sign in [(0,1),(4,-1)]:
        rows.extend((np.arange(4*n).reshape(n,4,1)+np.zeros((1,1,4),int)).ravel())
        cols.extend((4*np.arange(n)[:,None,None]+np.arange(4)[None,None,:]+shift+np.zeros((1,4,1),int)).ravel())
        values.extend((sign*correction).ravel())
    delta=coo_matrix((values,(rows,cols)),shape=D.shape).tocsc();delta.eliminate_zeros()
    return K,D+delta,forcing


def solve(n,bank,outer,label):
    coupled=SimpleNamespace(**dict(vars(task.coupled),assemble=assemble))
    runtime=SimpleNamespace(**dict(vars(task),coupled=coupled))
    ns=dict(vars(method.prior.old),OUT=OUT,task=runtime,Step=method.Step)
    exec(compile((method.prior.OUT/'solver-source.py').read_text(),str(method.prior.OUT/'solver-source.py'),'exec'),ns)
    return ns['solve'](n,bank,False,outer,label)


def prepare():
    assert not OUT.exists();OUT.mkdir()
    paths=[Path(__file__),Path(method.__file__),method.OUT/'result.json',method.OUT/'midpoint-mode-diagnostic.json',
        method.prior.OUT/'fine-bank.npz',method.prior.OUT/'coarse-bank.npz',method.prior.OUT/'solver-source.py']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        reassessment='Extended stage arithmetic leaves the failed response unchanged; fifth-order stepping still fails. The actual sampled modes show midpoint kinetic cancellation down to0.0312. Address the spatial principal inertia before another time refinement.',
        claim='Add the O(dx^2) constitutive acceleration correction whose constant-coefficient elimination yields the exact positive linear-element mass. Reapply this candidate to the actual GR response with unchanged input, static K,source,time counts and accuracy gates.',
        decision='A passing time response is only candidate discretization progress. Before physical adoption, require variable-coefficient energy/constraint and corrected material-pressure reconstruction. Do not call full GR/physical closure or silently use the modified auxiliary pressure as the old thermodynamic pressure.',
        correction='D_q,left += dx^2*(A*B)_q/12; D_q,right -= same, q rows0,2 only. Constant medium gives lambda=-12*sin(theta/2)^2/(dx^2*(2+cos(theta))) instead of -4*tan(theta/2)^2/dx^2.',
        budget=dict(paths=5,steps=240,hard_seconds=120,cpu_threads=1,memory_GB=3,new_EOS_calls=0,new_collision_states=0),
        forecast='Measured five three-stage paths50.7s. Same-size sparse correction forecast45-90s, hard120s. No additional mesh/time refinement or another candidate after failure.',
        gates=dict(time_relative=.02,time_order=1.5,coefficient_relative=.02,outer_relative=.002,linear_residual=1e-9,heat_balance=2e-13),
        bindings={str(p):digest(p) for p in paths}))
    x,y,h=sp.symbols('x y h',real=True);s=sp.symbols('s',positive=True)
    exact=sp.integrate(((1-s)*x+s*y)**2,(s,0,1));mid=((x+y)/2)**2
    assert sp.simplify(exact-mid-(x-y)**2/12)==0
    tan2=sp.symbols('t',nonnegative=True);lam=-4*tan2/(h*h*(1+tan2/3))
    assert sp.simplify(sp.limit(lam,tan2,sp.oo)+12/h**2)==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        identity='Integral of a linear field squared minus midpoint square equals (right-left)^2/12. Positive consistent mass eigenvalues are h/2,h/6; midpoint mass has a checkerboard null. Constant-wave frequency is bounded by sqrt(12)*c/h after the correction.',
        boundary='Uniform principal wave subsystem only; full variable GR operator energy, boundaries, auxiliary pressure and source-consistent reconstruction are not certified.'))


def run():
    assert not (OUT/'result.json').exists();signal.alarm(120);begin=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(int(3e9),int(3e9)))
    for p,h in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert digest(Path(p))==h,p
    cases={}
    for n in [16,32,64]:cases[str(n)]=solve(n,method.prior.OUT/'fine-bank.npz',2,'heat-'+str(n))
    cases['coarse']=solve(64,method.prior.OUT/'coarse-bank.npz',2,'coefficient-64')
    cases['outer']=solve(64,method.prior.OUT/'fine-bank.npz',3,'outer-64')
    comparisons={}
    for field in ['velocity_mass_RMS_m_s','scalar_mass_RMS','old_interface_velocity_RMS_m_s','new_interface_velocity_RMS_m_s']:
        series=lambda name:np.array([x[field] for x in cases[name]['history']])
        a,b,c=series('16'),series('32'),series('64');norm=max(abs(c).max(),1e-100)
        e1=float(max(abs(a-b[::2]))/norm);e2=float(max(abs(b-c[::2]))/norm)
        comparisons[field]=dict(time_previous=e1,time_last=e2,order=float(np.log2(e1/e2)),
            coefficients=float(max(abs(c-series('coarse')))/norm),outer=float(max(abs(c-series('outer')))/norm))
    balance=max(r['heat_telescoping'] for r in cases.values());residual=max(r['linear_residual'] for r in cases.values())
    passed=all(x['time_last']<.02 and x['order']>1.5 and x['coefficients']<.02 and x['outer']<.002 for x in comparisons.values()) and balance<2e-13 and residual<1e-9
    result=dict(classification='Counterexample candidate',actual_GR_candidate_evolved=True,numerical_gates_passed=passed,
        comparisons=comparisons,max_heat_balance=balance,max_linear_residual=residual,
        seconds=time.monotonic()-begin,memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
        full_variable_GR_energy_certified=False,modified_auxiliary_pressure_reconstructed=False,
        original_physical_failure_resolved=False,full_dynamic_charge_solved=False)
    write(OUT/'result.json',result);signal.alarm(0);print('INERTIA',json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
