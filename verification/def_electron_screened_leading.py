"""Leading screened direct+exchange probabilities, omitting their interference.

This is the small-transfer approximation specified after Eq. 12 of
Shternin-Yakovlev 2006. Do not combine differently oriented retarded
propagators into an unvalidated coherent exchange amplitude.
"""
from pathlib import Path
import json
import time
import numpy as np
from scipy.interpolate import CubicSpline
import def_electron_dynamic_screening as screen
import def_electron_screening_resume as recovery

exchange=screen.exchange
OUT=exchange.OUT/'leading-screening'


def amplitude(state,grid,values):
    spline=CubicSpline(grid,values,extrapolate=False)
    def probability(out1,in1,out2,in2):
        q=out1-in1;q2=np.sum(q*q,axis=1)
        Eo=np.sqrt(1+np.sum(out1*out1,axis=1));Ei=np.sqrt(1+np.sum(in1*in1,axis=1))
        omega=np.sum(q*(out1+in1),axis=1)/(Eo+Ei)
        phase=omega/np.sqrt(q2);assert np.max(abs(phase))<=grid[-1]
        Pi=spline(abs(phase));Pi=Pi.real+1j*np.sign(phase)[:,None]*Pi.imag
        J1=screen.currents(out1,in1);J2=screen.currents(out2,in2)
        charge=np.einsum('nai,nbj->nabij',J1[:,0],J2[:,0])
        vector=np.einsum('nkai,nkbj->nabij',J1[:,1:],J2[:,1:])
        term=charge/(q2+Pi[:,0])[:,None,None,None,None]-(vector-(omega*omega/q2)[:,None,None,None,None]*charge)/(q2-omega*omega+Pi[:,1])[:,None,None,None,None]
        return np.sum(abs(term)**2,axis=(1,2,3,4))
    def value(p1,p2,p3,p4,unused_qs2):
        return (4*np.pi*exchange.model.alpha)**2/4*(probability(p3,p1,p4,p2)+probability(p4,p1,p3,p2))
    return value


def main():
    assert not OUT.exists();OUT.mkdir();old=screen.OUT;h=screen.h
    parent=json.loads((old/'plan.json').read_text())
    for p,sha in parent['bindings'].items():assert h.digest(h.ROOT/p)==sha,p
    paths=[Path(__file__),Path(recovery.__file__),old/'plan.json',old/'resume-plan.json']+list(old.glob('polarization-*.npz'))
    parent['bindings'].update({p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths})
    parent['preserved_failure']='The coherent dynamically screened extension failed inverse-event probability equality: max relative 0.846846290632167 on the fixed control. Ward and static spin-trace controls passed. No production ran. Its retarded-channel exchange interference prescription was not validated and is rejected.'
    parent['declared_change']='Use the source small-transfer approximation: sum direct and exchange probabilities and omit their mutual interference, retaining longitudinal/transverse interference within each channel. This changes the approximation explicitly; it is not a correction proving the rejected full coherent model. Missing finite-q recoil and exchange-interference errors remain.'
    exchange.write(OUT/'plan.json',parent)
    for path in old.glob('polarization-*.npz'):(OUT/path.name).write_bytes(path.read_bytes())
    screen.OUT=OUT;screen.make_amplitude=amplitude
    records=[]
    for index in parent['cells']:
        state=exchange.equilibrium(index);d=np.load(OUT/f'polarization-{index}.npz');grid,values=d['phase'],d['polarization']
        sample=(np.linspace(0,len(grid)-2,16).astype(int)+.5)*grid[-1]/(len(grid)-1)
        exact=np.array([screen.polarization(a,state) for a in sample]);estimate=CubicSpline(grid,values)(sample)
        direct=exchange.model.kinetic(np.array([index]),256)['conductivity_SI'][0]
        projected=exchange.transfer(state,np.zeros_like(state['G']),7)['K']
        records.append(dict(cell=index,relative_static_compressibility=float(abs(values[0,0]/state['qs2']-1)),withheld_interpolation_relative=float(np.max(abs(exact-estimate))/state['qs2']),EI_direct_SI_comparison=float(abs(projected/direct-1))))
    state=exchange.equilibrium(parent['cells'][1]);d=np.load(OUT/f"polarization-{state['index']}.npz")
    amp=amplitude(state,d['phase'],d['polarization']);vertex=recovery.controls(state,amp)
    exchange.amplitude=amp;began=time.monotonic();exchange.bracket(state,12,exchange.SEEDS[0]);elapsed=time.monotonic()-began
    exchange.write(OUT/'preparation.json',dict(records=records,vertex_checks=vertex,pilot_seconds=elapsed,
        production_linear_forecast_seconds=elapsed*384,table_seconds=None,table_cost='Reused prior completed tables.',estimate='Same 4096-event pilot and 120s process cap; no event or grid expansion.'))
    print('PREPARATION',json.dumps(records),'VERTEX',vertex,'FORECAST',elapsed*384,flush=True)
    screen.run()


if __name__=='__main__':main()
