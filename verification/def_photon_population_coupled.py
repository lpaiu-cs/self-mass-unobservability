"""Counterexample candidate: finite ion-stage kinetics in the saved photon solve.

Internal levels and unretained chemical coordinates remain equilibrated. The
native constrained susceptibility, not an ideal-gas guess, supplies latent heat.
"""
from pathlib import Path
from types import FunctionType, SimpleNamespace
import argparse
import ctypes
import difflib
import json
import re
import shutil
import signal
import subprocess
import time
import numpy as np
from scipy import sparse
from scipy.linalg import lu_factor, lu_solve
from scipy.sparse.linalg import splu
import def_photon_hhe_coupled as old

OUT=old.OUT.parent/'def-photon-population-coupled'
CACHE=old.CACHE.parent/'photon-population-coupled'
LIB=CACHE/'libfree_eos_population.so'
write=old.write;digest=old.digest


def prepare():
    assert not (OUT/'build.json').exists();OUT.mkdir(exist_ok=True);CACHE.mkdir(exist_ok=True)
    if (OUT/'plan.json').exists():
        write(OUT/'setup-first-failure.json',dict(error='eos_calc.f90 is an included original source, not a copied cache file',native_builds=0,EOS_calls=0,coupled_paths=0))
    paths=[Path(__file__),old.LIB,old.CACHE/'gas.so',old.OUT/'bank.npz',
        old.OUT/'channels.npz',old.old.a.OUT/'eos-state.npz',old.old.a.OUT/'catalog.json',
        old.old.RATES/'cross-sections.npz',old.old.RATES/'rates.npz',
        old.OUT/'fine-64.npz',old.OUT/'response.json']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='aab855f9',
        claim='Evolve retained ion-stage populations, photons and material heat together on the original interval. Obtain the constrained native chemical susceptibility by conjugate external fields; remove its relaxational heat capacity from the thermal reservoir.',
        decision='Compare the actual finite-population path against the saved instantaneous-LTE path. No claim of full non-LTE levels, collision completeness, original complete-opacity repair or stellar GR completion.',
        model='18 retained ionization edges; higher-charge cumulative populations are conjugate to dimensionless species chemical fields. Other chemical coordinates and internal levels equilibrate conditionally. Native susceptibility and equilibrium temperature tangent define the local quadratic free energy. Bound-free events carry one ionization, one photon and reciprocal thermal/chemical energy. No invented relaxation time or rate multiplier.',
        budget=dict(build_seconds=60,EOS_calls=75,input_seconds=90,coupled_paths=5,coupled_seconds=480,cpu_threads=1,memory_GB=4,new_stellar_steps=0),
        forecast='Phase69 six paths218s,1.59GB; EOS call0.106s and build3.95s. New small reservoir Schur block is unmeasured; forecast180-420s, cap480s. Reuse saved LTE and zero-opacity controls; fixed24134 cells,8 angles,H/c duration,16/32/64 steps plus two input contrasts. No expansion after failure.',
        gates=dict(native_unbiased_bitwise=True,susceptibility_symmetry=1e-5,susceptibility_step=1e-5,positive_frozen_capacity=True,time_order=1.8,time_error=.001,input_error=.001,energy_balance=1e-9,solver=1e-11,entropy_growth=1e-10),
        bindings={str(p):digest(p) for p in paths}))
    source=old.old.a.old.p.SOURCE
    for p in old.CACHE.glob('*.mod'):shutil.copyfile(p,CACHE/p.name)
    owner=(old.CACHE/'mod_eos_calc.f90').read_text()
    owner=owner.replace('  public eos_calc','  public eos_calc\n  real(fp_kind), save, public :: population_field(318)=0._fp_kind')
    (CACHE/'mod_eos_calc.f90').write_text(owner)
    original=(source/'eos_calc.f90').read_text();s=original
    anchor='  ! Local variables';assert anchor in s
    s=s.replace(anchor,'  real(fp_kind) :: population_dv(size(dv))\n\n'+anchor,1)
    # Constant external fields change values, not derivatives with respect to
    # the EOS Newton variables. All callers share this lowest-level hook.
    start=s.index('  nionsp2 = ')
    s=s[:start]+re.sub(r'\bdv\b','population_dv',s[start:])
    start=s.index('  nionsp2 = ')
    s=s[:start]+'  population_dv=dv+population_field(:size(dv))\n'+s[start:]
    (CACHE/'eos_calc.f90').write_text(s)
    shutil.copyfile(source/'ionize.f90',CACHE/'ionize.f90')
    (OUT/'source.patch').write_text(''.join(difflib.unified_diff(original.splitlines(True),s.splitlines(True),fromfile='phase69/eos_calc.f90',tofile='phase70/eos_calc.f90')))
    begin=time.monotonic();commands=[
        ['gfortran','-cpp','-O2','-fPIC','-fcheck=all','-ffree-line-length-none','-I'+str(CACHE),'-I'+str(old.old.a.old.PARENT),'-I'+str(source),'-c','mod_eos_calc.f90','-o','mod_eos_calc.o']]
    objects=json.loads((old.old.a.OUT/'build.json').read_text())['reused_objects']
    objects=[p for p in objects if not p.endswith('/mod_eos_calc.o')]
    commands.append(['gfortran','-shared','-Wl,-soname,'+LIB.name,'mod_eos_calc.o',str(old.CACHE/'mod_excitation.o'),str(old.old.a.CACHE/'mod_free_eos_detailed.o'),str(old.old.a.CACHE/'mod_free_eos.o'),*objects,'-llapack','-lblas','-o',str(LIB)])
    bridge=old.old.a.old.p.a.OUT.parent/'gr-radiation-eos-split/gas-bridge.f90'
    commands.append(['gfortran','-O2','-fPIC','-shared','-I'+str(CACHE),'-I'+str(old.old.a.old.PARENT),str(bridge),'-L'+str(CACHE),'-Wl,-rpath,'+str(CACHE),'-lfree_eos_population','-o','gas.so'])
    receipts=[]
    for i,cmd in enumerate(commands):
        r=subprocess.run(cmd,cwd=CACHE,capture_output=True,text=True,timeout=60)
        (OUT/f'build-{i}.log').write_text(r.stdout+r.stderr)
        receipts.append(dict(command=cmd,returncode=r.returncode));assert r.returncode==0,r.stderr
    assert time.monotonic()-begin<60
    write(OUT/'build.json',dict(seconds=time.monotonic()-begin,receipts=receipts,
        bindings={str(p):digest(p) for p in [LIB,CACHE/'gas.so',CACHE/'mod_eos_calc.f90',CACHE/'eos_calc.f90']}))
    print('BUILT',time.monotonic()-begin,flush=True)


def thermo():
    assert not (OUT/'thermo.json').exists();signal.alarm(90);begin=time.monotonic()
    prior=old.old.a.old.p.a.previous;module=prior.old.previous.matter.old
    gas=object.__new__(module.GasEOS);init=module.GasEOS.__init__
    FunctionType(init.__code__,dict(init.__globals__,BRIDGE=CACHE/'gas.so'),closure=init.__closure__)(gas)
    gas.inventory_lib=gas.gas_lib;native=gas.gas_lib.ionization_inventory
    def call(mode,value,t,eps,out,info):
        raw=np.full(24,np.nan);native(mode,value,t,eps,raw,info)
        out[:]=raw[:22];out[20]=0.;gas.molecules=raw[22:].copy()
    gas.call=call
    fields=np.ctypeslib.as_array((ctypes.c_double*318).in_dll(gas.gas_lib,'__mod_eos_calc_MOD_population_field'))
    saved=dict(np.load(old.old.a.OUT/'eos-state.npz'));catalog=json.loads((old.old.a.OUT/'catalog.json').read_text())
    edges=[(r['element_index']-1,r['charge'],r['ion_index']-1) for r in catalog['rows']]+[(0,0,0),(1,1,2)]
    st=dict(np.load(old.old.RATES/'station.npz'));T=float(st['T']);rho=float(st['rho']);conversion=rho*float(st['mass_scale'])*6.02214076e23
    Z=[1,2,6,7,8,10,11,12,13,14,15,16,17,18,20,22,24,25,26,28,3,4,5,9]
    assert sum(Z)==316 and all(Z[r['element_index']-1]==r['Z'] for r in catalog['rows'])
    bias=np.zeros((318,len(edges)))
    for j,(e,q,i) in enumerate(edges):bias[i:i+Z[e]-q,j]=1.
    d,_=prior.base.inputs();r=float(d['lnd'][0]);X=d['X'][0];snapshots=[]
    def snapshot(source):
        fields[:]=bias@source
        s=prior.inventory_reader.InventoryEOS.snapshot(gas,r,float(np.log(T)),X)
        snapshots.append(s);return s
    def coords(s):return np.array([s[e,q+1:].sum() for e,q,i in edges])*conversion
    base=snapshot(np.zeros(len(edges)))
    same=all(np.array_equal(base[k],saved[k][0]) for k in base);assert same
    responses=[]
    for step in [1e-3,5e-4]:
        columns=[]
        for j in range(len(edges)):
            source=np.zeros(len(edges));source[j]=step
            plus=snapshot(source);minus=snapshot(-source)
            columns.append((coords(plus['number_fractions'])-coords(minus['number_fractions']))/(2*step))
        responses.append(np.array(columns).T)
    fields[:]=0
    replay=snapshot(np.zeros(len(edges)));assert all(np.array_equal(base[k],replay[k]) for k in base)
    a=[]
    for minus,plus,step in [(3,4,2e-4),(7,8,1e-4)]:
        a.append((coords(saved['number_fractions'][plus])-coords(saved['number_fractions'][minus]))/(2*step*T))
    S=responses[-1];scale=np.sqrt(np.diag(S));normalized=S/scale[:,None]/scale[None,:]
    symmetry=float(np.max(abs(normalized-normalized.T)))
    step_error=float(np.max(abs((responses[0]-S)/scale[:,None]/scale[None,:])))
    # Symmetrize only after a quantified native Maxwell-symmetry test.
    S=(S+S.T)/2;M=S/scale[:,None]/scale[None,:];ev=np.linalg.eigvalsh(M)
    k=float(catalog['constants'][-1]);a=np.array(a);H=k*T*np.linalg.inv(M)/scale[:,None]/scale[None,:]
    U=np.linalg.cholesky(T*H).T;g=U@a[-1]
    Ceq=float(np.load(old.OUT/'bank.npz')['Cm']);Cf=Ceq-g@g
    result=dict(classification='Counterexample candidate',same_EOS_bitwise=same,
        symmetry=symmetry,step_error=step_error,min_normalized_eigenvalue=float(ev[0]),
        equilibrium_heat_capacity=Ceq,frozen_heat_capacity=float(Cf),chemical_heat_capacity=float(g@g),
        temperature_tangent_difference=float(np.linalg.norm(U@(a[0]-a[1]))/np.sqrt(Ceq)),
        native_EOS_calls=len(snapshots),seconds=time.monotonic()-begin,
        limitation='Constrained ion-stage Hessian; intra-stage levels and all unretained coordinates still equilibrate. Biased EOS energy/entropy are not used as physical values; only its native abundance response is used.')
    write(OUT/'thermo.json',result)
    np.savez_compressed(OUT/'thermo.npz',S=S,raw_response=responses,a=a,U=U,g=g,Cf=Cf,Ceq=Ceq,edges=edges,bias=bias,coordinates=coords(base['number_fractions']),temperature=T)
    np.savez_compressed(OUT/'field-states.npz',**{key:np.array([s[key] for s in snapshots]) for key in base})
    assert symmetry<1e-5 and step_error<1e-5 and ev[0]>0 and Cf>0
    signal.alarm(0);print('THERMO',result,flush=True)


def inputs():
    assert not (OUT/'moments.npz').exists();signal.alarm(90);begin=time.monotonic()
    b=dict(np.load(old.OUT/'bank.npz'));data=dict(np.load(old.old.RATES/'rates.npz'))
    cross=dict(np.load(old.old.RATES/'cross-sections.npz'));hh=dict(np.load(old.OUT/'channels.npz'))
    catalog=json.loads((old.old.a.OUT/'catalog.json').read_text());erg,ryd,c2,c,k=catalog['constants'];h=erg/c
    T=float(b['T']);n=int((b['u']<=60).sum());edges=b['edges_u'][:n+1]
    capacity0=15*float(b['arad'])*T**3/np.pi**4
    atom=ctypes.CDLL(str(old.CACHE/'atomic.so'));fn=atom.hydrogen_cross
    array=np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS')
    fn.argtypes=[ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_double,array,array,array];fn.restype=None
    native=np.loadtxt(old.old.rates.OUT/'fort.93')
    def integrate(stride,order):
        moments=np.zeros((3,n,18));gx,gw=np.polynomial.legendre.leggauss(order)
        def add(stage,mesh,amplitude):
            if len(mesh)<2:return
            mid=(mesh[:-1]+mesh[1:])/2;half=np.diff(mesh)/2;u=mid[:,None]+half[:,None]*gx
            density=capacity0*u**4*np.exp(-u)/(-np.expm1(-u))**2
            weight=c*density*amplitude(u)*(-np.expm1(-u))*half[:,None]*gw
            ids=np.searchsorted(edges,mid,side='right')-1
            for power in range(3):np.add.at(moments[power,:,stage],ids,(weight/(k*T*u)**power).sum(1))
        query=cross['query'];meta=cross['metadata'];forward=data['bf_forward']
        for level in np.unique(query[:,0]):
            ids=np.flatnonzero(query[:,0]==level);ids=ids[meta[ids,3]>=0]
            ids=ids[np.unique(np.r_[np.arange(0,len(ids),stride),len(ids)-1])]
            u=meta[ids,2];amp=forward[ids];stage=int(meta[ids[0],0])
            if u[0]>=edges[-1]:continue
            upper=min(edges[-1],u[-1]);mesh=np.unique(np.r_[u[(u>=edges[0])&(u<=upper)],edges[(edges>=u[0])&(edges<=upper)],upper])
            add(stage,mesh,lambda uq:np.interp(uq,u,amp))
        for j,z,ground in [(0,1,1),(1,2,35)]:
            for level in range(10):
                threshold=hh['binding_cm_inverse'][j,level]*erg/(k*T)
                mesh=np.unique(np.r_[edges[edges>threshold],max(edges[0],threshold)])
                def amplitude(u):
                    freq=np.ascontiguousarray(((native[ground-1,3]/(level+1)**2+(u-threshold)*k*T)/old.old.rates.NATIVE_H).ravel())
                    bf=np.zeros_like(freq);ff=np.zeros_like(freq);fn(len(freq),z,level+1,T,freq,bf,ff)
                    return hh['populations'][j,level]*bf.reshape(u.shape)
                add(16+j,mesh,amplitude)
        return moments
    fine=integrate(1,8);coarse=integrate(2,4)
    # Metal source thinning formerly kept Gauss8. Its smooth subinterval
    # quadrature difference is now included in this stricter combined contrast.
    assert np.all(fine>=0) and np.all(coarse>=0)
    bf_rate=fine[0].sum(1)/b['Ci'][:n]
    remaining=b['rate_a_fine'][:n]-bf_rate
    assert remaining.min()>=-1e-10*max(b['rate_a_fine'])
    hhe_error=float(np.linalg.norm(fine[0,:,16:].sum(1)-c*hh['bf8'])/np.linalg.norm(c*hh['bf8']))
    assert hhe_error<1e-12
    np.savez_compressed(OUT/'moments.npz',fine=fine,coarse=coarse)
    write(OUT/'input.json',dict(classification='Counterexample candidate',seconds=time.monotonic()-begin,
        HHe_previous_input_relative=hhe_error,retained_ionization_edges=18,
        bound_free_photon_energy_moments=3,missing_collisional_rates=True,internal_levels_LTE=True,
        bindings={str(p):digest(p) for p in [OUT/'moments.npz',OUT/'thermo.npz',OUT/'thermo.json',old.OUT/'bank.npz']}))
    signal.alarm(0);print('INPUT',time.monotonic()-begin,flush=True)


def photon_operator(b):
    v=old.old.v;directory=old.OUT
    m=SimpleNamespace(**dict(vars(v.c.m),OUT=directory))
    spec=FunctionType(v.c.spectrum.__code__,dict(vars(v.c),OUT=directory,m=m))
    coef=FunctionType(v.c.raw_coefficients.__code__,dict(v.c.raw_coefficients.__globals__,OUT=directory,spectrum=spec),argdefs=v.c.raw_coefficients.__defaults__)
    class Operator(v.c.Operator):
        __init__=FunctionType(v.c.Operator.__init__.__code__,dict(v.c.Operator.__init__.__globals__,coefficients=coef),argdefs=v.c.Operator.__init__.__defaults__,closure=v.c.Operator.__init__.__closure__)
    return Operator(b,8,1,24)


def coupled_blocks(op,thermo,moments):
    G0,G1,G2=moments;Cf=float(thermo['Cf']);U=thermo['U'];a=thermo['a'][-1]
    latent=U.T@(U@a);s1=G1.sum(0);s2=G2.sum(0)
    A=np.zeros((19,19));A[0,0]=op.aa+((-2*latent*s1+latent**2*s2).sum())/Cf
    A[1:,1:]=(U*s2)@U.T;A[0,1:]=U@(s1-latent*s2)/np.sqrt(Cf);A[1:,0]=A[0,1:]
    # Only the isotropic photon moment exchanges population with the bath.
    B=np.empty((len(op.u),19));B[:,0]=op.q[::op.order]+G1@latent/np.sqrt(Cf*op.Ci)
    B[:,1:]=-(G1@U.T)/np.sqrt(op.Ci)[:,None]
    weights=np.r_[np.sqrt(Cf),U@a]
    # Integrated Cauchy-Schwarz checks nonnegative event covariance.
    gram=G0*G2-G1**2
    assert gram.min()>=-1e-12*max(float((G0*G2).max()),1.)
    null_m=A@weights+B.T@np.sqrt(op.Ci)
    null_p=op.L@op.energy+op.off(op.energy)
    null_p[::op.order]+=B@weights
    scale=max(np.linalg.norm(A@weights),np.linalg.norm(B.T@np.sqrt(op.Ci)),1.)
    assert np.linalg.norm(null_m)/scale<1e-10
    assert np.linalg.norm(null_p)/max(np.linalg.norm(op.L@op.energy),1.)<1e-10
    return A,B,weights


def evolve(op,A,B,weights,duration,steps):
    dt=duration/steps;gamma=1-1/np.sqrt(2);h=gamma*dt
    factor=splu(sparse.eye(op.size,format='csc')+h*op.P)
    fullB=np.zeros((op.size,19),complex);fullB[::op.order]=B
    R=factor.solve(h*fullB);schur=lu_factor(np.eye(19)+h*A-h*B.T@R[::op.order])
    del fullB
    contraction=h*op.off_bound;assert contraction<.5
    maxres=0.;maxiter=0
    # Keep the sparse frequency/angular operation and the small reservoir
    # matrix separate; never assemble a dense full collision matrix.
    def derivative(x,E):
        dx=-A@x-B.T@E[::op.order];de=-op.P@E-op.off(E)
        de[::op.order]-=B@x;return dx,de
    def solve(x,E):
        nonlocal maxres,maxiter
        def diagonal(erhs):
            y=factor.solve(erhs);z=lu_solve(schur,x-h*B.T@y[::op.order]);return z,y-R@z
        z,y=diagonal(E)
        for iteration in range(16):
            zz,yy=diagonal(E-h*op.off(y))
            delta=np.linalg.norm(yy-y);den=max(np.sqrt(np.vdot(yy,yy).real+np.vdot(zz,zz).real),1e-30)
            z,y=zz,yy
            if contraction/(1-contraction)*delta/den<1e-12:break
        else:raise AssertionError('finite-population fixed point failed')
        dx,de=derivative(z,y)
        residual=np.sqrt(np.linalg.norm(z-h*dx-x)**2+np.linalg.norm(y-h*de-E)**2)/den
        maxres=max(maxres,float(residual));maxiter=max(maxiter,iteration+1);assert residual<1e-11
        return z,y
    x=weights.astype(complex);E=np.zeros(op.size,complex)
    initial=float(weights@weights);integrated_flux=0j;balance=0.;entropy=0.;norm=initial
    def flux(e):return -1j*op.kc*(np.sqrt(op.Ci)@e.reshape(-1,op.order)[:,1])/np.sqrt(3)
    for step in range(steps):
        x1,e1=solve(x,E);f1,g1=derivative(x1,e1)
        x2,e2=solve(x+(1-gamma)*dt*f1,E+(1-gamma)*dt*g1)
        integrated_flux+=dt*((1-gamma)*flux(e1)+gamma*flux(e2))
        x,E=x2,e2
        energy=weights@x+op.energy@E
        balance=max(balance,float(abs(energy-initial-integrated_flux)/initial))
        nextnorm=float(np.vdot(x,x).real+np.vdot(E,E).real)
        entropy=max(entropy,(nextnorm-norm)/initial);norm=nextnorm
    return dict(x=x,E=E,balance=balance,entropy_growth=entropy,solver_residual=maxres,
        max_solver_iterations=maxiter,norm=norm,steps=steps)


def run():
    import resource
    assert not (OUT/'response.json').exists();signal.alarm(480);begin=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)))
    thermo=dict(np.load(OUT/'thermo.npz'));moments=dict(np.load(OUT/'moments.npz'))
    b=dict(np.load(old.OUT/'bank.npz'));b['Cm']=thermo['Cf'];b.update(velocity_order=256,q_points=257)
    op=photon_operator(b);duration=float(b['H']/b['c']);Ceq=float(thermo['Ceq']);rows=[]
    def path(name,steps,which='fine'):
        old.old.v.set_absorption(op,b['rate_a_'+which])
        A,B,weights=coupled_blocks(op,thermo,moments['coarse' if which=='coarse' else 'fine'])
        result=evolve(op,A,B,weights,duration,steps)
        record={k:v for k,v in result.items() if k not in ['x','E']};record['name']=name;rows.append(record)
        write(OUT/(name+'.json'),record)
        np.savez_compressed(OUT/(name+'.npz'),x=result['x'],E=result['E'],weights=weights,Cf=thermo['Cf'],Ceq=Ceq,duration=duration)
        print('POPULATION',name,'seconds',time.monotonic()-begin,'balance',result['balance'],flush=True)
        return result
    def difference(x,y):return float(np.sqrt(np.linalg.norm(x['x']-y['x'])**2+np.linalg.norm(x['E']-y['E'])**2)/np.sqrt(Ceq))
    solutions=[path('fine-'+str(n),n) for n in [16,32,64]]
    errors=[difference(x,y) for x,y in zip(solutions[:-1],solutions[1:])];order=float(np.log2(errors[0]/errors[1]))
    fine=solutions[-1];contrasts={}
    for name in ['coarse','quadrature4']:contrasts[name]=difference(fine,path(name+'-64',64,name))
    previous=dict(np.load(old.OUT/'fine-64.npz'));deltaT=previous['T']/np.sqrt(Ceq)
    oldstate=dict(x=np.r_[np.sqrt(thermo['Cf']),thermo['g']]*deltaT,E=previous['E'])
    effect=difference(fine,oldstate);temperature=fine['x'][0]/np.sqrt(thermo['Cf'])
    lag=fine['x'][1:]-thermo['g']*temperature
    tail=json.loads((old.OUT/'input.json').read_text())['line_integration'][-1]['omitted_wing_response_bound']
    tail_rescale=(1+op.Ci.sum()/float(thermo['Cf']))/(1+op.Ci.sum()/Ceq)
    tail*=tail_rescale
    passed=bool(order>=1.8 and errors[-1]<.001 and contrasts['coarse']+2*tail<.001 and contrasts['quadrature4']<.0001
        and all(r['balance']<1e-9 and r['solver_residual']<1e-11 and r['entropy_growth']<1e-10 for r in rows))
    result=dict(classification='Counterexample candidate',actual_coupled_evolution=True,numerical_response_passed=passed,
        time_order=order,time_differences_initial=errors,input_differences_initial=contrasts,
        source_plus_two_tails=contrasts['coarse']+2*tail,temperature_perturbation_K=float(temperature.real),
        line_tail_capacity_rescale=float(tail_rescale),
        LTE_temperature_perturbation_K=float(deltaT.real),finite_vs_LTE_difference=effect,
        population_disequilibrium_norm=float(np.linalg.norm(lag)/np.sqrt(Ceq)),
        duration_seconds=duration,rows=rows,seconds=time.monotonic()-begin,
        memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
        physical_opacity_certified=False,original_complete_input_failure_resolved=False,
        full_GR_photon_feedback_evolved=False,full_dynamic_charge_solved=False,
        boundary='Finite retained ion-stage radiative kinetics only. Internal levels, unretained chemistry and thermal velocities remain conditionally equilibrated. No electron-impact rates, complete redistribution, missing opacity, atmosphere or moving stellar GR.')
    write(OUT/'response.json',result);signal.alarm(0);print('FINISHED',result,flush=True)


def symbolic():
    import sympy as sp
    heat,latent,cf,rootchi,photon=sp.symbols('e l C U P',positive=True)
    v=sp.Matrix([(heat-latent)/sp.sqrt(cf),rootchi,-heat/photon])
    energy=sp.Matrix([sp.sqrt(cf),latent/rootchi,photon])
    assert sp.simplify(v.dot(energy))==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        identity='Each photoionization transfers e-l to heat, l to stored chemical energy, and removes e from photons. Its whitened rank-one positive collision bracket conserves their total. Cf=Ceq-T a^T H a avoids double counting; equilibrium lifting has energy/norm Ceq.',
        limitation='Conditional on a positive chemical Hessian and positive Cf; no proof of microscopic or full nonlinear physical completeness.'))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','thermo','inputs','run','symbolic'])
    globals()[parser.parse_args().action]()
