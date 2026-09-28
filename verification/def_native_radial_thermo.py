"""Use the saved native isentropic derivatives in the radial gas closure."""
from pathlib import Path
from types import FunctionType
import argparse
import inspect
import json
import signal
import time
import textwrap
import numpy as np
from scipy.interpolate import CubicHermiteSpline
import def_native_radial_release as task

OUT=task.OUT


class EOS(task.EOS):
    def __init__(self):
        super().__init__();raw=self.d['raw'];x=self.x
        self.logp=CubicHermiteSpline(x,np.log(raw[:,:,1]/(self.rho0*task.C**2)),raw[:,:,4],axis=1)
        self.gamma=self.logp.derivative()
        self.u=CubicHermiteSpline(x,raw[:,:,2]/task.C**2,raw[:,:,1]/raw[:,:,0]/task.C**2,axis=1)
        ad=(raw[:,:,1]/raw[:,:,0]-raw[:,:,9])/raw[:,:,10]
        self.logT=CubicHermiteSpline(x,np.log(self.d['T']),ad,axis=1)

    def __call__(self,rho,sigma):
        # The inherited call enforces the unchanged density/entropy domain
        # and supplies the unchanged native electron-scattering inventory.
        _,_,_,_,kap=super().__call__(rho,sigma)
        xx=np.log(np.maximum(rho,self.floor));z=np.clip(sigma,self.sigma[0],0)
        fraction=(z-self.sigma[0])/(self.sigma[1]-self.sigma[0]);ids=np.minimum(fraction.astype(int),1);f=fraction-ids;cols=np.arange(len(rho))
        def blend(fn):
            a=fn(xx);return (1-f)*a[ids,cols]+f*a[ids+1,cols]
        active=rho>=self.floor
        p=np.exp(blend(self.logp))*active;u=blend(self.u)*active;gamma=blend(self.gamma);T=np.exp(blend(self.logT))
        assert np.all(gamma>1) and np.all(p>=0),'Native Hermite thermodynamic stability'
        return p,u,gamma,T,kap


init_source=textwrap.dedent(inspect.getsource(task.Flow.__init__))
anchor="np.interp(ri,previous['radius'],v)"
assert init_source.count(anchor)==1
init_source=init_source.replace(anchor,"np.interp(ri,previous['radius'],np.asarray(v,float))")
init_namespace=dict(vars(task),EOS=EOS)
exec(compile(init_source,__file__,'exec'),init_namespace)


class Flow(task.Flow):
    __init__=init_namespace['__init__']


def audit():
    assert not (OUT/'thermo-control.json').exists()
    failure=json.loads((OUT/'eos.json').read_text());assert not failure['passed']
    task.write(OUT/'thermo-reassessment.json',dict(classification='Counterexample candidate',
        previous_maximum_relative=max(max(r['relative']) for r in failure['controls']),
        decision='Preserve the failed independent PCHIP gamma interpolation. Use dlnP/dlnrho=Gamma1, du/dlnrho=P/rho and dlnT/dlnrho=(P/rho-u_lnrho)/cvT from the SAME native isentropes in cubic Hermite interpolation. Gamma is now the derivative of the same pressure representation. No new grid, density/entropy domain or relaxed gate.',
        additional_native_call_cap=50,additional_seconds=20,total_call_cap_including_prior_failures=1270,
        controls='Repeat the original eight checks and add four unseen density/entropy controls near the chemical transition.',
        bindings={str(p.relative_to(task.old.ROOT)):task.old.photons.digest(p) for p in [Path(__file__),Path(task.__file__),OUT/'eos.npz',OUT/'eos.json']}))
    start=time.monotonic();signal.alarm(20);eos=EOS();d=eos.d;fan=task.prior.Fan(call_cap=50,reuse=True);rows=[]
    tests=[(x,f) for x in [-.037,-.7,-3.3,-10.3] for f in [.25,.75]]+[(x,f) for x,f in zip([-.14,-.52,-.86,-1.07],[.4,.6,.4,.6])]
    for xx,frac in tests:
        sig=float(d['sigma'][0])*frac;target=float(d['s0'])+float(d['sunit'])*sig
        p,u,gamma,T,kap=eos(np.array([np.exp(xx)]),np.array([sig]));guess=np.log(T[0])
        for _ in range(8):
            raw=fan.call(np.log(eos.rho0)+xx,guess);error=(raw[3]-target)*np.exp(guess)/raw[10]
            if abs(error)<2e-12:break
            guess-=error
        else:raise AssertionError('Native control root')
        errors=[float(abs(p[0]*eos.rho0*task.C**2/raw[1]-1)),float(abs(u[0]*task.C**2/raw[2]-1)),float(abs(gamma[0]/raw[4]-1)),float(abs(T[0]/np.exp(guess)-1))]
        rows.append(dict(log_density=xx,entropy_fraction=frac,relative=errors,new_control=len(rows)>=8))
        task.write(OUT/'thermo-control-progress.json',dict(rows=rows,native_calls=fan.calls))
    passed=max(max(r['relative']) for r in rows)<.002
    data=dict(classification='Counterexample candidate',passed=passed,controls=rows,native_calls=fan.calls,
        total_native_calls=failure['total_native_calls']+fan.calls,seconds=time.monotonic()-start,
        maximum_relative=max(max(r['relative']) for r in rows),same_native_states=True,whole_physical_EOS_certified=False)
    task.write(OUT/'thermo-control.json',data);signal.alarm(0);print(json.dumps(data),flush=True);assert passed


# Reuse the frozen flow driver and its unchanged gates. The failed EOS audit
# remains failed; only the new derivative-consistent producer can use its own
# new audit. The replacement is checked and its source is saved for review.
source=inspect.getsource(task.run)
anchor="(OUT/'eos.json').read_text()"
assert source.count(anchor)==1
source=source.replace(anchor,"(OUT/'thermo-control.json').read_text()")
source=source.replace('signal.alarm(240)','signal.alarm(235)')
namespace=dict(vars(task),Flow=Flow)
exec(compile(source,__file__,'exec'),namespace)


def run():
    spec=json.loads((OUT/'thermo-reassessment.json').read_text())
    for p,h in spec['bindings'].items():
        actual=OUT/'registered-thermo.py' if p=='verification/def_native_radial_thermo.py' else task.old.ROOT/p
        assert task.old.photons.digest(actual)==h,p
    task.write(OUT/'initialization-repair.json',dict(classification='Counterexample candidate',
        failed_before_first_fluid_step=True,error='numpy.interp rejects the saved float128 boundary-velocity array',
        correction='Convert only the imported boundary-velocity interpolation ordinate to binary64, the declared flow arithmetic; preserve all stored bulk fields and the unchanged EOS, fluid equations and gates.',
        failed_process_wall_seconds=2.4481423,reserved_repair_seconds=5,remaining_production_cap_seconds=235,
        source_sha256=task.old.photons.digest(Path(__file__)),prior_source_sha256=task.old.photons.digest(OUT/'registered-thermo.py')))
    (OUT/'reused-flow-init.py').write_text(init_source,encoding='utf-8')
    (OUT/'reused-thermo-driver.py').write_text(source,encoding='utf-8')
    namespace['run']()


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['audit','run']);globals()[parser.parse_args().action]()
