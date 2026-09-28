"""Finite-T long-wavelength longitudinal/transverse collision kernel.

Reuses the relativistic conserving events and collision brackets. The
semiclassical polarization averages the angular response over the actual
Fermi distribution. Quantum recoil / full finite-q RPA remains excluded.
"""
from pathlib import Path
import argparse
import json
import time
import resource
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline
from scipy.special import expit,spence
from scipy.stats import qmc
import def_electron_energy_exchange as exchange

h=exchange.h
OUT=exchange.OUT/'dynamic-screening'
PAULI=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)
IDENTITY=np.eye(2)
STATIC_AMPLITUDE=exchange.amplitude


def polarization(a,state):
    """Retarded-sign convention matches chi_t=+i*pi*a/(4*v) at small a."""
    theta=state['theta'];eta=state['eta'];top=eta+60
    def density(x):
        E=1+theta*x;p=np.sqrt(theta*x*(2+theta*x));f=expit(eta-x)
        return f*(1-f),E,p
    i0=quad(lambda x:density(x)[0]*density(x)[1]*density(x)[2],0,top,epsabs=1e-12,epsrel=2e-10,limit=150)[0]
    factor=4*exchange.model.alpha/np.pi
    if a==0:return factor*i0+0j,0j
    Ea=1/np.sqrt(1-a*a);xa=(Ea-1)/theta
    def integrand(x,component):
        weight,E,p=density(x);v=p/E
        logarithm=np.log((a+v)/abs(a-v)) if v!=a else 0.
        polynomial=E*E if component==0 else (1-a*a)*E*E-1
        return weight*polynomial*logarithm
    points=[xa] if 0<xa<top else None
    j1=quad(integrand,0,top,args=(0,),points=points,epsabs=1e-11,epsrel=2e-9,limit=200)[0]
    j2=quad(integrand,0,top,args=(1,),points=points,epsabs=1e-11,epsrel=2e-9,limit=200)[0]
    z=eta-xa;F=expit(z);logarithm=np.logaddexp(0,z)
    fermione=-spence(1+np.exp(z))
    E2tail=Ea*Ea*F+2*theta*Ea*logarithm+2*theta*theta*fermione
    longitudinal=factor*(i0-a*j1/2-1j*np.pi*a*E2tail/2)
    transverse=factor*(a*a*i0/2+a*j2/4+1j*np.pi*a*((1-a*a)*E2tail-F)/4)
    return longitudinal,transverse


def table(state,count=257):
    max_energy=1+2*state['theta']*(state['eta']+40)
    maximum=np.sqrt(1-1/max_energy**2)*(1+1e-12)
    grid=np.linspace(0,maximum,count)
    values=np.array([polarization(a,state) for a in grid])
    return grid,values


def currents(out,inp):
    eo=np.sqrt(1+np.sum(out*out,axis=1));ei=np.sqrt(1+np.sum(inp*inp,axis=1))
    scale=np.sqrt((eo+1)*(ei+1))
    a,b=exchange.density_vertex(out,inp)
    J=np.empty((len(inp),4,2,2),complex)
    J[:,0]=a[:,None,None]*IDENTITY+1j*np.einsum('nk,kab->nab',b,PAULI)
    vp=inp/(ei+1)[:,None];vo=out/(eo+1)[:,None]
    for j in range(3):
        unit=np.zeros(3);unit[j]=1
        b=np.cross(unit,vp-vo)*scale[:,None]
        a=(vp[:,j]+vo[:,j])*scale
        J[:,j+1]=a[:,None,None]*IDENTITY+1j*np.einsum('nk,kab->nab',b,PAULI)
    return J


def make_amplitude(state,grid,values):
    spline=CubicSpline(grid,values,extrapolate=False)
    def channel(out1,in1,out2,in2):
        q=out1-in1;q2=np.sum(q*q,axis=1)
        Eo=np.sqrt(1+np.sum(out1*out1,axis=1));Ei=np.sqrt(1+np.sum(in1*in1,axis=1))
        # Stable energy difference, especially for near-forward scattering.
        omega=np.sum((out1-in1)*(out1+in1),axis=1)/(Eo+Ei)
        phase=omega/np.sqrt(q2);assert np.max(abs(phase))<=grid[-1]
        Pi=spline(abs(phase));Pi=Pi.real+1j*np.sign(phase)[:,None]*Pi.imag
        J1=currents(out1,in1);J2=currents(out2,in2)
        charge=np.einsum('nai,nbj->nabij',J1[:,0],J2[:,0])
        vector=np.einsum('nkai,nkbj->nabij',J1[:,1:],J2[:,1:])
        longitudinal=charge/(q2+Pi[:,0])[:,None,None,None,None]
        transverse=(vector-(omega*omega/q2)[:,None,None,None,None]*charge)/(q2-omega*omega+Pi[:,1])[:,None,None,None,None]
        return longitudinal-transverse
    def result(p1,p2,p3,p4,unused_qs2):
        direct=channel(p3,p1,p4,p2)
        swap=channel(p4,p1,p3,p2).swapaxes(1,2)
        M=direct-swap
        value=(4*np.pi*exchange.model.alpha)**2/4*np.sum(abs(M)**2,axis=(1,2,3,4))
        assert np.all(np.isfinite(value)) and np.min(value)>0
        return value
    return result


def controls(state,amp):
    checks=dict(Ward_relative=0.,reverse_relative=0.,spin_trace_relative=0.)
    def audited(p1,p2,p3,p4,qs2):
        for out,inp in [(p3,p1),(p4,p2),(p4,p1),(p3,p2)]:
            J=currents(out,inp);q=out-inp
            eout=np.sqrt(1+np.sum(out*out,axis=1));ein=np.sqrt(1+np.sum(inp*inp,axis=1))
            omega=np.sum(q*(out+inp),axis=1)/(eout+ein)
            residue=np.einsum('nk,nkab->nab',q,J[:,1:])-omega[:,None,None]*J[:,0]
            scale=np.linalg.norm(q,axis=1)*np.linalg.norm(J[:,1:],axis=(1,2,3))
            checks['Ward_relative']=max(checks['Ward_relative'],float(np.max(np.linalg.norm(residue,axis=(1,2))/np.maximum(scale,1e-30))))
        value=amp(p1,p2,p3,p4,qs2);reverse=amp(p3,p4,p1,p2,qs2)
        checks['reverse_relative']=float(np.max(abs(value/reverse-1)))
        a=currents(p3,p1)[:,0];b=currents(p4,p2)[:,0]
        d=currents(p4,p1)[:,0];e=currents(p3,p2)[:,0]
        qd=np.sum((p3-p1)**2,axis=1)+qs2;qe=np.sum((p4-p1)**2,axis=1)+qs2
        M=np.einsum('nai,nbj->nabij',a,b)/qd[:,None,None,None,None]-np.einsum('nbi,naj->nabij',d,e)/qe[:,None,None,None,None]
        trace=(4*np.pi*exchange.model.alpha)**2/4*np.sum(abs(M)**2,axis=(1,2,3,4))
        checks['spin_trace_relative']=float(np.max(abs(trace/STATIC_AMPLITUDE(p1,p2,p3,p4,qs2)-1)))
        return value
    exchange.amplitude=audited
    exchange.events(qmc.Sobol(5,scramble=True,seed=5598).random_base2(10),state)
    exchange.amplitude=amp
    assert max(checks.values())<2e-12,checks
    return checks


def prepare():
    assert not OUT.exists();OUT.mkdir()
    parent=json.loads((exchange.OUT/'plan.json').read_text())
    paths=[Path(__file__),Path(exchange.__file__),exchange.OUT/'plan.json',exchange.OUT/'result.json']
    plan=dict(classification='Counterexample candidate',
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths},cells=parent['cells'],seeds=parent['seeds'],powers=parent['powers'],basis_sizes=parent['basis_sizes'],
        source='Shternin-Yakovlev 2006 equations 5-9. Finite-temperature semiclassical extension explicitly averages the same angular functions over the ideal Fermi energy derivative; does not supply finite-q quantum recoil.',
        claim='Replace the static density-only ee interaction by longitudinal and transverse dynamically screened Dirac currents in the same conserving event integral, with paired samples and unchanged physical states.',
        gates=parent['gates']|dict(polarization_interpolation_relative=1e-5,polarization_table_bracket_relative=.002,Ward_identity_relative=2e-12,EI_direct_SI_comparison=.02),
        budget=dict(table_nodes=257,withheld_points=16,pilot_power=12,production_hard_seconds=120,CPU_workers=1,native_calls=0,stellar_steps=0,automatic_expansion=False))
    exchange.write(OUT/'plan.json',plan)
    records=[];start=time.monotonic()
    for index in plan['cells']:
        state=exchange.equilibrium(index);grid,values=table(state)
        np.savez_compressed(OUT/f'polarization-{index}.npz',phase=grid,polarization=values)
        sample=(np.linspace(0,len(grid)-2,16).astype(int)+.5)*grid[-1]/(len(grid)-1)
        exact=np.array([polarization(a,state) for a in sample]);estimate=CubicSpline(grid,values)(sample)
        error=float(np.max(abs(exact-estimate))/state['qs2'])
        direct=exchange.model.kinetic(np.array([index]),256)['conductivity_SI'][0]
        projected=exchange.transfer(state,np.zeros_like(state['G']),7)['K']
        SI_error=float(abs(projected/direct-1))
        records.append(dict(cell=index,relative_static_compressibility=abs(values[0,0]/state['qs2']-1),withheld_interpolation_relative=error,EI_direct_SI_comparison=SI_error))
    state=exchange.equilibrium(plan['cells'][1]);data=np.load(OUT/f"polarization-{state['index']}.npz")
    exchange.amplitude=make_amplitude(state,data['phase'],data['polarization'])
    vertex_checks=controls(state,exchange.amplitude)
    began=time.monotonic();exchange.bracket(state,12,exchange.SEEDS[0]);elapsed=time.monotonic()-began
    exchange.write(OUT/'preparation.json',dict(records=records,vertex_checks=vertex_checks,table_seconds=began-start,pilot_seconds=elapsed,
        production_linear_forecast_seconds=elapsed*3*4*32,estimate='Paired same-event count, linear pilot scaling; 120s process cap unchanged.'))
    print('PREPARE',json.dumps(records,default=lambda v:float(v.real)),'FORECAST',elapsed*384,flush=True)


def run():
    assert not (OUT/'result.json').exists();start=time.monotonic();plan=json.loads((OUT/'plan.json').read_text());prep=json.loads((OUT/'preparation.json').read_text())
    for p,sha in plan['bindings'].items():assert h.digest(h.ROOT/p)==sha,p
    assert prep['production_linear_forecast_seconds']<100
    assert max(row['withheld_interpolation_relative'] for row in prep['records'])<1e-5
    assert max(row['EI_direct_SI_comparison'] for row in prep['records'])<.02
    rows=[]
    for index in plan['cells']:
        state=exchange.equilibrium(index);data=np.load(OUT/f'polarization-{index}.npz')
        grid,values=data['phase'],data['polarization']
        exchange.amplitude=make_amplitude(state,grid,values)
        fine=[];coarse=[];checks=[]
        for seed in exchange.SEEDS:
            a,b,d=exchange.bracket(state,17,seed,15);fine.append(a);coarse.append(b);checks.append(d)
        EE=np.mean(fine,axis=0);small=np.mean(coarse,axis=0)
        results=[exchange.transfer(state,EE,n) for n in [3,5,7]]
        qc=exchange.transfer(state,small,7);individual=[exchange.transfer(state,e,7) for e in fine]
        errors={key:dict(quadrature=abs(qc[key]/results[-1][key]-1),basis=abs(results[-2][key]/results[-1][key]-1),scramble_relative_std=float(np.std([v[key] for v in individual],ddof=1)/results[-1][key])) for key in ['K','tau']}
        # Independent table-size control uses the same already registered seed.
        exchange.amplitude=make_amplitude(state,grid[::2],values[::2])
        low_table,_,_=exchange.bracket(state,15,exchange.SEEDS[0])
        table_error=float(np.linalg.norm(low_table-coarse[0])/np.linalg.norm(coarse[0]))
        previous=np.load(exchange.OUT/f'cell-{index}.npz');static_K=float(previous['fine_K'])
        np.savez_compressed(OUT/f'cell-{index}.npz',**state,EE=EE,coarse_EE=small,replicate_EE=np.array(fine),**{'fine_'+key:v for key,v in results[-1].items()})
        row=dict(cell=index,K=results[-1]['K'],tau=results[-1]['tau'],min_projected_rate=results[-1]['min_rate'],errors=errors,table_bracket_relative=table_error,native_ratio=results[-1]['K']/state['native_K'],static_screening_K_ratio=results[-1]['K']/static_K,diagnostics=checks)
        rows.append(row);print('CELL',index,'NATIVE',row['native_ratio'],'STATIC',row['static_screening_K_ratio'],'ERRORS',json.dumps(errors),flush=True)
    passed=all(row['table_bracket_relative']<.002 and all(e['quadrature']<.03 and e['basis']<.02 for e in row['errors'].values()) for row in rows)
    result=dict(classification='Counterexample candidate',rows=rows,numerical_gates_passed=passed,
        native_compatibility_passed=all(abs(r['native_ratio']-1)<.2 for r in rows),
        elapsed_seconds=time.monotonic()-start,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        finite_q_quantum_polarization_included=False,continuum_inverse_moment_certified=False,whole_star_closure_replaced=False,full_dynamic_charge_solved=False)
    exchange.write(OUT/'result.json',result);print('RESULT',passed,'SECONDS',result['elapsed_seconds'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
