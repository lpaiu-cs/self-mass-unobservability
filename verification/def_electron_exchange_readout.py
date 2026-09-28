"""Low-momentum scattering and finite-basis heat-memory readout."""
from pathlib import Path
import json
import time
import numpy as np
import sympy as sp
from scipy.stats import qmc
import def_electron_screened_fermi as resolved

exchange=resolved.exchange
OUT=exchange.OUT/'collision-readout'


def main():
    h=exchange.h;source=resolved.OUT
    assert not OUT.exists();assert json.loads((source/'result.json').read_text())['numerical_gates_passed'];OUT.mkdir()
    plan=json.loads((source/'plan.json').read_text());paths=[Path(__file__),Path(resolved.__file__),source/'plan.json',source/'result.json',exchange.model.base.BENCHMARK]
    paths+=list(source.glob('cell-*.npz'))+list(source.glob('polarization-*.npz'))
    exchange.write(OUT/'plan.json',dict(classification='Counterexample candidate',
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths},cells=plan['cells'],
        purpose='Check the actual linearized ee out-scattering rate at zero and small momentum, then read the constrained finite-basis transfer. Do not infer a continuum inverse-moment bound from a few positive loss rates.',
        momenta_over_pF=[0.,1e-4,.01],seeds=[5511,5512],powers=[13,15],
        gates=dict(loss_quadrature_relative=.03,zero_momentum_limit_relative=.001),
        budget=dict(hard_seconds=60,CPU_workers=1,events=3*3*2*2**15,native_calls=0,stellar_steps=0,automatic_expansion=False)))
    s,l=sp.symbols('s lambda',positive=True)
    assert sp.integrate(sp.exp(-l*s),(s,0,sp.oo))==1/l
    exchange.write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        finite_basis='For K(s)=sum_j w_j/(lambda_j+s), w_j>=0 and lambda_j>0, tau_eff=sum(w_j/lambda_j^2)/sum(w_j/lambda_j). The zero-current unit step rises as 1-sum[(w_j/lambda_j)exp(-lambda_j*t)]/sum(w_j/lambda_j).',
        harmonic_bound='abs(K(i*omega)/K(0)-1) <= abs(omega)*tau_eff; follows termwise from abs[1/(lambda+i*omega)-1/lambda] <= abs(omega)/lambda^2.',
        limitation='A statement about the fixed Galerkin operator only. It does not bound spatial heat diffusion, global stellar response, omitted collision physics, or the continuum inverse moment.'))
    coefficient=np.load(exchange.model.base.thermal.OUT/'coefficients.npz')
    benchmark=json.loads(exchange.model.base.BENCHMARK.read_text());coordinate=3*2*np.pi/benchmark['leading_drive']['period_seconds']
    rows=[];start=time.monotonic()
    for index in plan['cells']:
        state=exchange.equilibrium(index);data=np.load(source/f'polarization-{index}.npz')
        exchange.amplitude=resolved.leading.amplitude(state,data['phase'],data['polarization'])
        rates=[];coarse=[]
        for seed in [5511,5512]:
            samples=qmc.Sobol(5,scramble=True,seed=seed).random_base2(15)
            each=[];smaller=[]
            for fraction in [0.,1e-4,.01]:
                total=0.;initial=0.
                for first in range(0,len(samples),4096):
                    w,_,diag=exchange.events(samples[first:first+4096],state,fixed_p=fraction*state['xF'])
                    assert diag['energy']<2e-13 and diag['momentum']<2e-13 and diag['detailed_balance_log']<2e-10
                    total+=w.sum()
                    if first<2**13:initial+=w.sum()
                each.append(total/2**15*exchange.RATE);smaller.append(initial/2**13*exchange.RATE)
            rates.append(each);coarse.append(smaller)
        mean=np.mean(rates,axis=0);qc=np.mean(coarse,axis=0);quadrature=float(np.max(abs(qc/mean-1)))
        limit=float(abs(mean[1]/mean[0]-1));assert np.min(mean)>0
        result=np.load(source/f'cell-{index}.npz');poles=result['fine_poles'];weights=result['fine_weights']
        normalized=weights/poles;normalized/=normalized.sum();tau=float(result['fine_tau'])
        factors=np.array([0.,1.,5.,20.]);step=1-np.sum(normalized[:,None]*np.exp(-poles[:,None]*tau*factors),axis=0)
        proper=coordinate/(coefficient['A'][index]*coefficient['N'][index])
        row=dict(cell=index,out_scattering_rates=mean.tolist(),loss_quadrature_relative=quadrature,zero_momentum_limit_relative=limit,
            finite_basis_tau_seconds=tau,finite_basis_third_harmonic_memory_upper=proper*tau,
            unit_step_times_over_tau=factors.tolist(),unit_step_normalized_flux=step.tolist())
        rows.append(row);np.savez_compressed(OUT/f'cell-{index}.npz',loss_rates=np.array(rates),coarse_loss_rates=np.array(coarse),step_times_seconds=tau*factors,step_normalized_flux=step)
        print('CELL',index,'LOSS0',mean[0],'LIMIT',limit,'QUADRATURE',quadrature,flush=True)
    result=dict(classification='Counterexample candidate',rows=rows,
        numerical_gates_passed=all(r['loss_quadrature_relative']<.03 and r['zero_momentum_limit_relative']<.001 for r in rows),
        low_momentum_out_scattering_no_longer_vanishes=True,continuum_first_memory_certified=False,
        elapsed_seconds=time.monotonic()-start,new_native_calls=0,new_stellar_steps=0,full_dynamic_charge_solved=False)
    exchange.write(OUT/'result.json',result)


if __name__=='__main__':main()
