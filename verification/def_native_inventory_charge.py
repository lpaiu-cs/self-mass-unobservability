"""Same direct scalar readout, now with the actual fixed-inventory outflow."""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import json
import signal
import time
import numpy as np
import def_native_inventory_flow as flow
import def_native_release_charge as read

OUT=flow.OUT/'charge'


def adapter(n):
    f=flow.DiluteFlow(n);m=f.base;m.eos=f.eos;m.primitive=f.primitive
    return m


def main():
    assert not OUT.exists();OUT.mkdir();begin=time.monotonic();signal.alarm(45)
    assert json.loads((flow.OUT/'dilute/result.json').read_text())['passed']
    assert json.loads((flow.OUT/'dilute/audit.json').read_text())['passed']
    flow.write(OUT/'plan.json',dict(classification='Counterexample candidate',seconds=45,EOS_calls=30,new_fluid_steps=0,
        claim='Apply the completed frozen-inventory conservative histories to the exact same direct retarded Green and local acoustic readout used for LTE. Compare on448/896 grids and against the saved same-grid896 LTE wave.',
        gates=dict(wave_grid=.02,native_bulk_gamma=.002),
        boundary='Uniform initial chemical inventories and no finite reactions; prescribed nonlinear gas with local acoustic bulk. Full scalar/metric feedback, absorptive photons and a physical tail are not certified. A surviving direct component is not the final charge.',
        inputs={str(p):flow.cold.sha(p) for p in [Path(__file__),Path(flow.__file__),Path(read.__file__),flow.OUT/'dilute-eos.npz',flow.OUT/'dilute/cells-448.npz',flow.OUT/'dilute/cells-896.npz']}))
    m=adapter(896);x=-20000.;r=m.R+x/m.As
    rho=float(np.exp(np.interp(r,m.env['r'],np.log(m.env['rho']))));T=float(np.exp(np.interp(r,m.env['r'],np.log(m.env['T']))))
    native=flow.Native(cap=30);xx=np.log(rho/m.eos.rho0);lt=np.log(T);raw=native(xx,lt);h=1e-4
    ar=(native(xx+h,lt)-native(xx-h,lt))/(2*h);at=(native(xx,lt+h)-native(xx,lt-h))/(2*h)
    gamma=ar[1]/raw[1]+at[1]/raw[1]*(raw[1]/raw[0]-ar[2])/at[2]
    native_biased_gamma=float(raw[4]);raw[4]=gamma
    p,u,g,_,_=m.eos(np.array([rho/m.eos.rho0]),np.array([lt]));error=abs(float(g[0])/gamma-1);assert error<.002
    np.savez_compressed(OUT/'native-bulk.npz',rho=rho,T=T,raw=raw[None],native_biased_gamma=native_biased_gamma,fixed_gamma=gamma)
    flow.write(OUT/'bulk-thermo.json',dict(classification='Counterexample candidate',passed=True,EOS_calls=native.ion.calls,
        fixed_gamma=float(gamma),interpolated_gamma_relative=float(error),equilibrium_derivative_reused=False))
    ns=dict(vars(read),OUT=OUT,prior=SimpleNamespace(Flow=adapter,OUT=flow.OUT/'dilute'))
    reader=FunctionType(read.readout.__code__,ns,argdefs=read.readout.__defaults__)
    old=np.load(flow.previous.OUT/'charge/cells-896-linear-g12.npz');times=old['u_seconds']
    fine,fr=reader(896,'linear',times,raw[None]);coarse,cr=reader(448,'linear',times,raw[None])
    peak=max(abs(fine));error=float(max(abs(fine-coarse))/peak)
    change=float(max(abs(fine-old['normalized_charge']))/max(abs(old['normalized_charge'])))
    result=dict(classification='Counterexample candidate',passed=error<.02,wave_grid_relative=error,
        fixed_endpoint=float(fine[-1]),LTE_same_grid_endpoint=float(old['normalized_charge'][-1]),
        relative_change_from_LTE=change,endpoint_ratio=float(fine[-1]/old['normalized_charge'][-1]),
        same_sign=bool(fine[-1]*old['normalized_charge'][-1]>0),endpoint_components_cm=fr['endpoint_components_cm'],
        seconds=time.monotonic()-begin,actual_inventory_flow_applied=True,retarded_direct_component_computed=True,
        finite_chemistry=False,full_GR_scalar_feedback=False,final_charge_solved=False,full_goal_complete=False)
    flow.write(OUT/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True);assert result['passed']


if __name__=='__main__':main()
