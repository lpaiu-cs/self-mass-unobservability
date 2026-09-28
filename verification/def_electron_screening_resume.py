"""Reuse completed polarization tables after the control-only norm failure."""
from pathlib import Path
import json
import time
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.stats import qmc
import def_electron_dynamic_screening as screen


def controls(state,amp):
    exchange=screen.exchange
    checks=dict(Ward_relative=0.,reverse_relative=0.,spin_trace_relative=0.)
    def audited(p1,p2,p3,p4,qs2):
        for out,inp in [(p3,p1),(p4,p2),(p4,p1),(p3,p2)]:
            J=screen.currents(out,inp);q=out-inp
            eout=np.sqrt(1+np.sum(out*out,axis=1));ein=np.sqrt(1+np.sum(inp*inp,axis=1))
            omega=np.sum(q*(out+inp),axis=1)/(eout+ein)
            residue=np.einsum('nk,nkab->nab',q,J[:,1:])-omega[:,None,None]*J[:,0]
            scale=np.linalg.norm(q,axis=1)*np.sqrt(np.sum(abs(J[:,1:])**2,axis=(1,2,3)))
            checks['Ward_relative']=max(checks['Ward_relative'],float(np.max(np.linalg.norm(residue,axis=(1,2))/np.maximum(scale,1e-30))))
        value=amp(p1,p2,p3,p4,qs2);reverse=amp(p3,p4,p1,p2,qs2)
        checks['reverse_relative']=float(np.max(abs(value/reverse-1)))
        a=screen.currents(p3,p1)[:,0];b=screen.currents(p4,p2)[:,0]
        d=screen.currents(p4,p1)[:,0];e=screen.currents(p3,p2)[:,0]
        qd=np.sum((p3-p1)**2,axis=1)+qs2;qe=np.sum((p4-p1)**2,axis=1)+qs2
        M=np.einsum('nai,nbj->nabij',a,b)/qd[:,None,None,None,None]-np.einsum('nbi,naj->nabij',d,e)/qe[:,None,None,None,None]
        trace=(4*np.pi*exchange.model.alpha)**2/4*np.sum(abs(M)**2,axis=(1,2,3,4))
        checks['spin_trace_relative']=float(np.max(abs(trace/screen.STATIC_AMPLITUDE(p1,p2,p3,p4,qs2)-1)))
        return value
    exchange.amplitude=audited
    exchange.events(qmc.Sobol(5,scramble=True,seed=5598).random_base2(10),state)
    exchange.amplitude=amp
    assert max(checks.values())<2e-12,checks
    return checks


def main():
    out=screen.OUT;exchange=screen.exchange;h=screen.h
    assert not (out/'preparation.json').exists()
    plan=json.loads((out/'plan.json').read_text())
    for p,sha in plan['bindings'].items():assert h.digest(h.ROOT/p)==sha,p
    paths=[Path(__file__),Path(screen.__file__),out/'plan.json']+list(out.glob('polarization-*.npz'))
    exchange.write(out/'resume-plan.json',dict(classification='Counterexample candidate',
        failure='Three polarization tables completed; the Ward-control scale called numpy.linalg.norm with three axes and raised ValueError: Improper number of dimensions to norm. No dynamic-screening production ran.',
        change='Use sqrt(sum(abs(J)^2)) for that diagnostic Frobenius norm. Preserve the original source, polarization arrays, physics, grids, event seeds and gates. No table rerun.',
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths}))
    records=[]
    for index in plan['cells']:
        state=exchange.equilibrium(index);d=np.load(out/f'polarization-{index}.npz');grid,values=d['phase'],d['polarization']
        sample=(np.linspace(0,len(grid)-2,16).astype(int)+.5)*grid[-1]/(len(grid)-1)
        exact=np.array([screen.polarization(a,state) for a in sample]);estimate=CubicSpline(grid,values)(sample)
        direct=exchange.model.kinetic(np.array([index]),256)['conductivity_SI'][0]
        projected=exchange.transfer(state,np.zeros_like(state['G']),7)['K']
        records.append(dict(cell=index,relative_static_compressibility=float(abs(values[0,0]/state['qs2']-1)),
            withheld_interpolation_relative=float(np.max(abs(exact-estimate))/state['qs2']),EI_direct_SI_comparison=float(abs(projected/direct-1))))
    state=exchange.equilibrium(plan['cells'][1]);d=np.load(out/f"polarization-{state['index']}.npz")
    amp=screen.make_amplitude(state,d['phase'],d['polarization']);vertex_checks=controls(state,amp)
    exchange.amplitude=amp;began=time.monotonic();exchange.bracket(state,12,exchange.SEEDS[0]);elapsed=time.monotonic()-began
    exchange.write(out/'preparation.json',dict(records=records,vertex_checks=vertex_checks,pilot_seconds=elapsed,
        production_linear_forecast_seconds=elapsed*384,table_seconds=None,table_cost='Reused all three completed tables; failed preparation did not persist its original walltime.',
        estimate='Same 4096-event pilot scaling and 120s production cap.'))
    print('CONTROL',json.dumps(records),'VERTEX',vertex_checks,'FORECAST',elapsed*384,flush=True)
    screen.run()


if __name__=='__main__':main()
