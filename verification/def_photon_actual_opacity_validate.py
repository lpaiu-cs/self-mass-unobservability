"""Audit direct atomic spectra and couple their declared LTE absorption candidate.

Native population/emission failures are independent of the numerical response
verdict. A positive LTE opacity operator cannot certify those physical inputs.
"""
from pathlib import Path
import argparse
import json
import resource
import signal
import time
import numpy as np
import sympy as sp
from scipy.sparse import diags
import def_photon_actual_opacity as a
import def_photon_collective_continuum as c
from def_photon_finite_jump_validate import algebra,diagnostics

OUT=a.OUT;old=a.previous.old


def spectrum(name):
    path=OUT/name
    assert json.loads((path/'run.json').read_text())['returncode']==0
    raw=np.loadtxt(path/'fort.88')
    raw=raw[np.argsort(raw[:,0])]
    frequency,idx=np.unique(raw[:,0],return_index=True)
    # Repeated batch endpoints can differ in approximate line-wing selection.
    # Keep the first, report the duplicate spread rather than hiding it.
    count=np.diff(np.r_[idx,len(raw)])
    lo=np.minimum.reduceat(raw[:,1],idx);hi=np.maximum.reduceat(raw[:,1],idx)
    assert np.all(np.isfinite(raw)) and np.all(raw[:,1]>0)
    return raw[idx],dict(samples=len(raw),unique_samples=len(idx),duplicate_count=int((count-1).sum()),
        duplicate_opacity_relative=float(np.max(hi/lo-1)),
        native_Kirchhoff_max=float(np.max(abs(raw[:,2]/raw[:,1]-1))),
        frequency_Hz=[float(frequency[0]),float(frequency[-1])])


def populations():
    s=np.load(a.previous.OUT/'inventory.npz');b=np.load(old.OUT/'bank.npz')
    Z=list(a.previous.inventory_reader.g.d.CHARGES);native={};implicit={}
    text=(OUT/'balanced-fine/fort.28').read_text().splitlines()
    for line in text:
        p=line.split()
        if p[0]=='STAGE':native[(int(p[1]),int(p[2]))]=float(p[3])
        if p[0]=='ION':implicit[(int(p[1]),int(p[2]))]=float(p[5])
        if p[0]=='STATE':state=[float(v) for v in p[1:]]
    assert np.allclose(state,[float(b['T']),float(b['rho']),float(b['ne'])/1e6],rtol=1e-13,atol=0)
    rows=[]
    for (z,q),value in native.items():
        expected=float(s['ni'][Z.index(z),q]/1e6)
        rows.append(dict(Z=z,charge=q,native_cm3=value,EOS_cm3=expected,relative=value/expected-1))
    charge=sum(q*v for (z,q),v in native.items())
    # F has no explicit atom. Its implicit ion inventory still carries charge.
    charge+=sum(q*v for (z,q),v in implicit.items() if z==9)
    return dict(classification='Counterexample candidate',rows=rows,T_rho_ne_matched=True,
        explicit_charge_relative=charge/(float(b['ne'])/1e6)-1,
        maximum_stage_relative=max(abs(r['relative']) for r in rows),
        passed=bool(max(abs(r['relative']) for r in rows)<.001),
        limitation='Same total elemental ratios and T/rho/ne do not imply the same ion/level populations. Native atomic mass and partition/occupation conventions differ; no population rescaling was performed.')


def bank():
    assert not (OUT/'bank.npz').exists()
    b=dict(np.load(old.OUT/'bank.npz'));T=float(b['T']);h=2*np.pi*old.base.hbar
    infrared,irinfo=spectrum('infrared');reports={'infrared':irinfo};rates={}
    keep=b['u']<=60;N=int(keep.sum());edges=b['edges_u'][:N+1]
    x,w=np.polynomial.legendre.leggauss(8)
    capacity0=15*float(b['arad'])*T**3/np.pi**4
    for label,name in [('coarse','balanced-coarse'),('fine','balanced-fine')]:
        raw,info=spectrum(name);reports[label]=info
        # Each half-open band owns its samples. No blending of different T.
        raw=np.vstack([infrared[infrared[:,0]<raw[0,0]],raw])
        u=h*raw[:,0]/(old.base.k*T);chi=raw[:,1]
        # The piecewise-linear spectral interpolant is a declared numerical
        # representation. Split integrals at EVERY source knot before mapping
        # onto the existing photon cells, so narrow lines are not point sampled.
        mesh=np.unique(np.r_[edges,u[(u>edges[0])&(u<edges[-1])]])
        midpoint=(mesh[:-1]+mesh[1:])/2;half=np.diff(mesh)/2
        q=midpoint[:,None]+half[:,None]*x
        cq=np.interp(q,u,chi)
        density=capacity0*q**4*np.exp(-q)/(-np.expm1(-q))**2
        mass=(density*w*half[:,None]).sum(axis=1)
        amount=(density*cq*w*half[:,None]).sum(axis=1)
        ids=np.searchsorted(edges,midpoint,side='right')-1
        C=np.bincount(ids,mass,minlength=N);Q=np.bincount(ids,amount,minlength=N)
        cap_error=float(np.max(abs(C/b['Ci'][:N]-1)))
        assert cap_error<1e-8,cap_error
        rate=np.array(b['rate_a']);rate[:N]=float(b['c'])*Q/C
        assert np.all(np.isfinite(rate)) and np.all(rate>0)
        rates[label]=rate
        reports[label].update(capacity_relative=cap_error,spectral_segments=len(midpoint),
            u_range=[float(u[0]),float(u[-1])],below_source_capacity_fraction=float(np.sum(mass[midpoint<u[0]])/C.sum()),
            mapping='Exact source-knot splits; positive 8-point integration of piecewise-linear volume opacity with Planck derivative. Constant lowest opacity below first source frequency. No mean or relaxation-rate fitting.')
    np.savez_compressed(OUT/'bank.npz',**b,rate_a_fine=rates['fine'],rate_a_coarse=rates['coarse'])
    p=populations();a.write(OUT/'populations.json',p)
    a.write(OUT/'spectral-audit.json',dict(classification='Counterexample candidate',spectra=reports,
        direct_temperature=True,populations_passed=p['passed'],
        native_emission_passed=max(r['native_Kirchhoff_max'] for r in reports.values())<1e-6,
        physical_opacity_certified=False,
        scope='Conditional LTE absorption candidate on the inherited finite photon bins. Emission is defined by Kirchhoff in that candidate; this does not reclassify the native emission residual as passing. Subcell transport/line redistribution, missing atomic continua and mass/partition consistency are not bounded by the source-grid comparison.'))
    # A reproducible analytic check for the actual source correction.
    U,x,B=sp.symbols('U x B',positive=True)
    chi=U*(1-sp.exp(-x));j=U*B*sp.exp(-x)
    assert sp.simplify(j/chi-B/(sp.exp(x)-1))==0
    a.write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        statement='For an LTE transition with common positive unstimulated opacity U_nu, chi_nu=U_nu(1-exp(-h nu/kT)) and j_nu=U_nu(2h nu^3/c^2)exp(-h nu/kT) obey j_nu=chi_nu B_nu. Replacing both factors by batch-center values does not preserve this identity at other frequencies.',
        limit='The identity assumes LTE populations of the SAME microscopic model. It supplies neither cross-section accuracy nor an EOS-population reconciliation.'))
    print('BANK',reports,'EOS charge',p['explicit_charge_relative'],flush=True)


def set_absorption(op,rate):
    rate=rate[:len(op.u)];delta=rate-op.rate_a
    op.L=op.L+diags(np.repeat(delta,op.order),format='csc')
    op.P=op.L+1j*op.kc*op.V
    op.q[::op.order]=op.jump_q-rate*np.sqrt(op.Ci/op.Cm)
    op.aa=op.jump_aa+float(rate@op.Ci/op.Cm)
    op.rate_a=rate.copy();op.max_iterations=0;op.max_solver_error=0.


def run():
    assert not (OUT/'response.json').exists();signal.alarm(480)
    plan=json.loads((OUT/'response-plan.json').read_text())
    for name,value in plan['bindings'].items():assert a.digest(old.h.ROOT/name)==value,name
    resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)))
    start=time.monotonic();b=dict(np.load(OUT/'bank.npz'));b['rate_a']=b['rate_a_fine']
    b.update(velocity_order=256,q_points=257)
    op=c.Operator(b,8,1,24);duration=float(b['H']/b['c']);Cm=float(b['Cm'])
    rows=[];runs=[]
    for n in [16,32,64]:
        r=op.evolve(duration,n);runs.append(r)
        row=dict(name=f'fine-{n}',algebra=algebra(op,b),**diagnostics(op,r));rows.append(row)
        a.write(OUT/f'fine-{n}.json',row)
        print('RESPONSE',n,row['balance'],row['solver_residual'],flush=True)
    errors=[old.compare(r,s,op) for r,s in zip(runs[:-1],runs[1:])]
    order=float(np.log2(errors[0]/errors[1]));fine=runs[-1]
    np.savez_compressed(OUT/'fine-response.npz',T=fine['T'],E=fine['E'],u=op.u,Ci=op.Ci,Cm=Cm,duration=duration)
    set_absorption(op,b['rate_a_coarse']);coarse=op.evolve(duration,64)
    b['rate_a']=b['rate_a_coarse']
    row=dict(name='coarse-64',algebra=algebra(op,b),**diagnostics(op,coarse));rows.append(row)
    a.write(OUT/'coarse-64.json',row)
    np.savez_compressed(OUT/'coarse-response.npz',T=coarse['T'],E=coarse['E'],u=op.u,Ci=op.Ci,Cm=Cm,duration=duration)
    grid=old.compare(fine,coarse,op)
    prior=np.load(c.OUT/'q-grid.npz')
    assert np.array_equal(op.u,prior['u']) and np.array_equal(op.Ci,prior['Ci'])
    difference=float(np.sqrt(abs(fine['T']-prior['T'])**2+np.linalg.norm(fine['E']-prior['E'].ravel())**2)/np.sqrt(Cm))
    passed=bool(order>=1.8 and errors[-1]<.001 and grid<.001 and all(
        r['balance']<1e-9 and r['energy_residual']<1e-9 and r['solver_residual']<1e-11
        and r['entropy_growth']<1e-10 and max(r['algebra']['energy_number_null_relative'])<1e-10 for r in rows))
    result=dict(classification='Counterexample candidate',numerical_response_passed=passed,
        physical_opacity_certified=False,populations_passed=False,native_emission_passed=False,
        source_grid_difference_initial=grid,time_order=order,time_differences_initial=errors,
        new_atomic_model_difference_initial=difference,
        temperature_perturbation_K=float((fine['T']/np.sqrt(Cm)).real),
        temperature_difference_K=float(((fine['T']-prior['T'])/np.sqrt(Cm)).real),
        rows=rows,coefficients=op.info,seconds=time.monotonic()-start,
        memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
        decision='The direct-T atomic/LTE candidate is coupled, but cannot replace the prior opacity as a certified same-EOS input. Its difference includes the atomic model, populations, lines and spectral representation; it is NOT a measured TOPS temperature-interpolation error.',
        atmosphere_closed=False,GR_feedback_complete=False,full_dynamic_charge_solved=False)
    a.write(OUT/'response.json',result);print('RESULT',result,flush=True);signal.alarm(0)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['bank','run']);globals()[parser.parse_args().action]()
