"""Conservative finite-T electron-electron collision brackets, not a BGK rate.

Natural units hbar=m_e=c=1 inside the collision integral. Exact relativistic
two-body kinematics and spin-averaged direct/exchange density vertices;
static longitudinal screening only. Transverse/dynamic screening is absent.
The finite Galerkin calculation is not a continuum spectral-gap certificate.
"""
from pathlib import Path
import argparse
import json
import time
import resource
import urllib.request
import numpy as np
import sympy as sp
from scipy.constants import hbar,m_e,c,k,elementary_charge,epsilon_0
from scipy.special import expit
from scipy.stats import qmc
from scipy.linalg import cholesky,solve_triangular,null_space,eigh
from numpy.polynomial.legendre import legvander
import def_electron_response_regular as regular

model=regular.model
model.angular=regular.angular
h=model.h
OUT=model.base.OUT.parent/'def-electron-energy-exchange'
RATE=m_e*c*c/hbar
KUNIT=(m_e*c/hbar)**3*k*c*c/RATE
DEGREE=6
SEEDS=[5501,5502,5503,5504]


def write(path,value):
    # numpy scalar diagnostics are converted only at the JSON boundary.
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False,
        default=lambda v:v.item() if isinstance(v,np.generic) else v.tolist())+'\n')


def equilibrium(index,order=256):
    row=model.kinetic(np.array([index]),order)
    d,data=model.base.inputs();theta=k*data['T'][index]/(m_e*c*c)
    eta=float(row['eta'][0]);x=row['x'][0];p=np.sqrt(theta*x*(2+theta*x));E=np.sqrt(1+p*p)
    # W dx = Nunit/m_e * p^3/(3*pi^2*E) f(1-f) dx.
    # G has p^4 dp/(3*pi^2), with dp/dx=theta*E/p.
    measure=row['weights'][0]/row['tau'][0]/(m_e*c/hbar)**3*m_e*theta*E*E
    basis=legvander((x-eta)/4,DEGREE)
    G=np.einsum('ni,nj,n->ij',basis,basis,measure)
    EI=np.einsum('ni,nj,n->ij',basis,basis,measure/(row['tau'][0]*RATE))
    B=np.c_[basis.T@(measure/E),basis.T@(measure*(x-eta)/E)]
    # Compressibility of the same ideal finite-T electron distribution.
    raw=model.distribution(np.array([data['T'][index]]),np.array([data['ne'][index]]),order)
    qs2=hbar*hbar*elementary_charge**2/epsilon_0*raw[5][0]/(m_e*c)**2
    return dict(index=int(index),theta=theta,eta=eta,qs2=qs2,G=G,EI=EI,B=B,
        native_K=float(row['native_conductivity_SI'][0]),published_nuee=float(row['ee_rate'][0]),
        ne=float(data['ne'][index]),T=float(data['T'][index]),xF=float(data['xF'][index]))


def energy(u,state):
    """Exact truncated logistic proposal on kinetic energy x in [0,eta+40]."""
    f0=expit(state['eta']);fmax=expit(-40.);normal=f0-fmax
    occupancy=f0-u*normal
    x=state['eta']+np.log1p(-occupancy)-np.log(occupancy)
    return x,occupancy,normal


def density_vertex(out,inp):
    eo=np.sqrt(1+np.sum(out*out,axis=1));ei=np.sqrt(1+np.sum(inp*inp,axis=1))
    scale=np.sqrt((eo+1)*(ei+1))
    scalar=scale+np.sum(out*inp,axis=1)/scale
    vector=np.cross(out,inp)/scale[:,None]
    return scalar,vector


def multiply(a,b):
    # (a0 I+i a.sigma)(b0 I+i b.sigma).
    return a[0]*b[0]-np.sum(a[1]*b[1],axis=1),a[0][:,None]*b[1]+b[0][:,None]*a[1]-np.cross(a[1],b[1])


def amplitude(p1,p2,p3,p4,qs2):
    a=density_vertex(p3,p1);b=density_vertex(p4,p2)
    d=density_vertex(p4,p1);e=density_vertex(p3,p2)
    norm=lambda j:2*(j[0]*j[0]+np.sum(j[1]*j[1],axis=1))
    interference=2*multiply(multiply(multiply(a,(d[0],-d[1])),b),(e[0],-e[1]))[0]
    direct=np.sum((p3-p1)**2,axis=1)+qs2;exchange=np.sum((p4-p1)**2,axis=1)+qs2
    total=norm(a)*norm(b)/direct**2+norm(d)*norm(e)/exchange**2-2*interference/(direct*exchange)
    assert np.min(total)>0
    return (4*np.pi*model.alpha)**2/4*total


def events(sample,state,fixed_p=None):
    theta=state['theta'];eta=state['eta'];n=len(sample)
    x1,f1,norm=energy(sample[:,0],state);x2,f2,_=energy(sample[:,1],state)
    p1mag=np.sqrt(theta*x1*(2+theta*x1));p2mag=np.sqrt(theta*x2*(2+theta*x2))
    if fixed_p is not None:
        p1mag=np.full(n,fixed_p);x1=p1mag*p1mag/(np.sqrt(1+p1mag*p1mag)+1)/theta;f1=expit(eta-x1)
    mu=2*sample[:,2]-1;p1=np.c_[np.zeros(n),np.zeros(n),p1mag]
    p2=np.c_[p2mag*np.sqrt(1-mu*mu),np.zeros(n),p2mag*mu]
    E1=1+theta*x1;E2=1+theta*x2
    dot=E1*E2-p1mag*p2mag*mu;s=2+2*dot;root=np.sqrt(s);gamma=(E1+E2)/root
    beta=(p1+p2)/(E1+E2)[:,None]
    pstar=p1+((gamma*gamma/(gamma+1))*np.sum(beta*p1,axis=1)-gamma*E1)[:,None]*beta
    kmag=np.sqrt(np.sum(pstar*pstar,axis=1));direction=pstar/kmag[:,None]
    delta=state['qs2']/(2*kmag*kmag)
    branch=sample[:,3]<.5;u=np.where(branch,2*sample[:,3],2*sample[:,3]-1)
    y=2*u*delta/(delta+2*(1-u));cosine=np.where(branch,1-y,y-1)
    pdf=.25*delta*(delta+2)*(1/(1-cosine+delta)**2+1/(1+cosine+delta)**2)
    angle=2*np.pi*sample[:,4];other=np.c_[-direction[:,2],np.zeros(n),direction[:,0]]
    transverse=np.sin(angle)[:,None]*other
    transverse[:,1]+=np.cos(angle)
    p3star=kmag[:,None]*(cosine[:,None]*direction+np.sqrt(np.maximum(0,1-cosine*cosine))[:,None]*transverse)
    bp=np.sum(beta*p3star,axis=1);Estar=root/2
    p3=p3star+((gamma*gamma/(gamma+1))*bp+gamma*Estar)[:,None]*beta
    p4=-p3star+(-(gamma*gamma/(gamma+1))*bp+gamma*Estar)[:,None]*beta
    E3=np.sqrt(1+np.sum(p3*p3,axis=1));E4=np.sqrt(1+np.sum(p4*p4,axis=1))
    z=np.array([x1-eta,x2-eta,(E3-1)/theta-eta,(E4-1)/theta-eta]).T
    loghole=-np.logaddexp(0,-z);logf=-np.logaddexp(0,z)
    detailed=logf[:,0]+logf[:,1]+loghole[:,2]+loghole[:,3]-logf[:,2]-logf[:,3]-loghole[:,0]-loghole[:,1]
    # Full final solid angle is used: 1/2 prevents identical-final double counting.
    cross=amplitude(p1,p2,p3,p4,state['qs2'])/(128*np.pi**2*s)
    speed=np.sqrt((dot-1)*(dot+1))/(E1*E2)
    diagnostics=dict(energy=float(np.max(abs(E1+E2-E3-E4))),
        momentum=float(np.max(np.linalg.norm(p1+p2-p3-p4,axis=1))),
        detailed_balance_log=float(np.max(abs(detailed))))
    if fixed_p is not None:
        weight=2/np.pi*theta*p2mag*E2*norm*np.exp(loghole[:,2]+loghole[:,3]-loghole[:,0]-loghole[:,1])*speed*cross/pdf
        return weight,None,diagnostics
    # g^2/4=1 for spin multiplicity g=2. Global orientation averages the
    # vector collision bracket to its dot product divided by three.
    weight=theta**2*norm**2*p1mag*E1*p2mag*E2*np.exp(loghole[:,2]+loghole[:,3]-loghole[:,0]-loghole[:,1])*speed*cross/(6*np.pi**3*pdf)
    polynomials=legvander(z/4,DEGREE)
    change=(p1[:,None,:]*polynomials[:,0,:,None]+p2[:,None,:]*polynomials[:,1,:,None]
            -p3[:,None,:]*polynomials[:,2,:,None]-p4[:,None,:]*polynomials[:,3,:,None])
    assert diagnostics['energy']<2e-13 and diagnostics['momentum']<2e-13 and diagnostics['detailed_balance_log']<2e-10
    # This column is an exact collision invariant; remove its roundoff only
    # after checking the event-level conserved momentum above.
    change[:,0,:]=0
    return weight,change,diagnostics


def bracket(state,power,seed,coarse_power=None):
    samples=qmc.Sobol(5,scramble=True,seed=seed).random_base2(power)
    matrix=np.zeros((DEGREE+1,DEGREE+1));coarse=np.zeros_like(matrix);checks=dict(energy=0.,momentum=0.,detailed_balance_log=0.)
    for start in range(0,len(samples),4096):
        w,delta,diag=events(samples[start:start+4096],state)
        value=np.einsum('n,nia,nja->ij',w,delta,delta,optimize=True)
        matrix+=value
        if coarse_power is not None and start<2**coarse_power:coarse+=value
        for key in checks:checks[key]=max(checks[key],diag[key])
    return matrix/len(samples),coarse/2**coarse_power if coarse_power else None,checks


def transfer(state,EE,size):
    G=state['G'][:size,:size];C=state['EI'][:size,:size]+EE[:size,:size];B=state['B'][:size]
    L=cholesky(G,lower=True)
    whiten=lambda a:solve_triangular(L,solve_triangular(L,a,lower=True).T,lower=True).T
    D=whiten(C);b=solve_triangular(L,B,lower=True)
    Q=null_space(b[:,0][None,:]);rates,vectors=eigh(Q.T@D@Q)
    assert np.min(rates)>0
    weights=(vectors.T@Q.T@b[:,1])**2
    static=np.sum(weights/rates);memory=np.sum(weights/rates**2)/static/RATE
    frequencies=np.array([.01,1.,100.])/memory
    response=np.array([np.sum(weights/(rates+1j*w/RATE))/static for w in frequencies])
    return dict(K=KUNIT*static,tau=memory,min_rate=float(rates.min()*RATE),
        poles=rates*RATE,weights=weights,frequency_response=response)


def symbolic():
    a,b,c0=sp.symbols('a b c',positive=True);d=a*b/c0
    f=lambda v:1/(1+v)
    assert sp.factor(f(a)*f(b)*(1-f(c0))*(1-f(d))-f(c0)*f(d)*(1-f(a))*(1-f(b)))==0
    return dict(classification='Proven',passed=True,
        detailed_balance='For a_i=exp((E_i-mu)/T), E1+E2=E3+E4 implies a1*a2=a3*a4 and f1*f2*(1-f3)*(1-f4)=f3*f4*(1-f1)*(1-f2).',
        bracket='C_ee[i,j]=(g^2/4) integral d3p1 d3p2/(2pi)^6 f1 f2 (1-f3)(1-f4) v_M d_sigma Delta_psi_i dot Delta_psi_j /3, g=2; d_sigma has identical-final factor 1/2 and spin average 1/4.',
        invariants='Each event conserves 1,E,p exactly. Positive event weights make the linearized collision form symmetric positive semidefinite. Momentum column is an exact null vector, not a relaxation term.',
        constrained_response='After G-whitening and zero-current projection, K(s)=KUNIT sum w_j/(lambda_j+s/RATE), w_j>=0, lambda_j>0 for the computed finite basis. This establishes finite-system passivity and a decaying homogeneous response.',
        limitation='No continuum gap or physical error bound follows from a finite Galerkin gap.')


def controls(state):
    momentum=1e-4;angle=.7;qs2=.3*momentum**2
    p1=np.array([[0.,0.,momentum]]);p2=-p1
    p3=np.array([[momentum*np.sin(angle),0.,momentum*np.cos(angle)]]);p4=-p3
    direct=np.sum((p3-p1)**2)+qs2;exchange=np.sum((p4-p1)**2)+qs2
    expected=16*(4*np.pi*model.alpha)**2*(direct**-2+exchange**-2-1/(direct*exchange))
    actual=amplitude(p1,p2,p3,p4,qs2)[0];nr_error=abs(actual/expected-1)
    assert nr_error<1e-7
    sample=qmc.Sobol(5,scramble=True,seed=5599).random_base2(10)
    w,delta,checks=events(sample,state)
    matrix=np.einsum('n,nia,nja->ij',w,delta,delta)/len(sample)
    assert np.array_equal(matrix[0],np.zeros(DEGREE+1))
    assert np.linalg.eigvalsh(matrix)[0]>-1e-12*np.trace(matrix)
    # Recover a constant-nu conductivity independently from its direct
    # thermoelectric moments. This catches natural-unit/measure mistakes.
    C=np.copy(state['G']);G=state['G'];B=state['B']
    R=np.linalg.solve(G,B);direct_K=KUNIT*(B[:,1]@R[:,1]-(B[:,0]@R[:,1])**2/(B[:,0]@R[:,0]))
    modified=dict(state,EI=C)
    reconstructed=transfer(modified,np.zeros_like(C),7)['K']
    assert abs(reconstructed/direct_K-1)<1e-12
    return dict(classification='Counterexample candidate',NR_spin_exchange_relative=nr_error,
        constant_rate_projection_relative=float(abs(reconstructed/direct_K-1)),
        exact_momentum_null=True,event_checks=checks,passed=True)


def prepare():
    assert not OUT.exists();OUT.mkdir()
    url='https://arxiv.org/html/astro-ph/0608371'
    with urllib.request.urlopen(url,timeout=20) as reply:(OUT/'shternin-yakovlev2006.html').write_bytes(reply.read())
    previous=regular.OUT/'response.npz';indices=np.load(previous)['indices']
    selected=indices[np.array([0,len(indices)//2,-1])]
    paths=[Path(__file__),Path(model.__file__),Path(regular.__file__),Path(model.base.__file__),model.base.thermal.OUT/'coefficients.npz',previous,OUT/'shternin-yakovlev2006.html']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='72662fc5',
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths},source=url,
        claim='Construct and solve an energy-exchanging electron collision operator on actual core states; test whether the previous low-momentum elastic obstruction survives explicit Pauli-blocked ee scattering.',
        model='Exact relativistic dispersion and two-body kinematics; spin-averaged direct/exchange static longitudinal screened Born density interaction in the plasma rest frame; ideal Fermi electrons; previous correlated-ion elastic kernel retained. No transverse/dynamic screening, fitted relaxation floor or fitted conductivity prefactor.',
        comparison='Published 1999 DC ee rate and native conductivity are independent comparisons; neither calibrates this collision operator. Finite basis passivity and sampled scattering do not certify the full plasma or full stellar response.',
        cells=selected.tolist(),seeds=SEEDS,powers=[15,17],basis_sizes=[3,5,7],energy_domain='0 <= kinetic energy/(kT) <= eta+40; FD thermal proposal with exact density/Jacobian.',
        gates=dict(event_energy=2e-13,event_momentum=2e-13,detailed_balance_log=2e-10,quadrature_relative=.03,basis_relative=.02,native_K_compatibility=.2,equilibrium_quadrature=1e-6),
        budget=dict(pilot_cells=1,pilot_power=12,production_hard_seconds=120,CPU_workers=1,native_calls=0,stellar_steps=0,automatic_expansion=False)))
    state=equilibrium(selected[1]);start=time.monotonic();bracket(state,12,SEEDS[0]);elapsed=time.monotonic()-start
    write(OUT/'pilot.json',dict(seconds=elapsed,production_linear_forecast_seconds=elapsed*3*4*2**(17-12),
        basis='Measured 4096-event collision bracket; linear event-count assumption, input and output overhead excluded. Stop at 120 seconds.'))
    write(OUT/'symbolic.json',symbolic())
    write(OUT/'controls.json',controls(state))
    print('PILOT',elapsed,'FORECAST',elapsed*3*4*32,flush=True)


def run():
    assert not (OUT/'result.json').exists();start=time.monotonic();plan=json.loads((OUT/'plan.json').read_text())
    for p,sha in plan['bindings'].items():assert h.digest(h.ROOT/p)==sha,p
    assert json.loads((OUT/'pilot.json').read_text())['production_linear_forecast_seconds']<100
    rows=[]
    for index in plan['cells']:
        state=equilibrium(index);check=equilibrium(index,512)
        eqerr=max(float(np.linalg.norm(state[key]-check[key])/np.linalg.norm(check[key])) for key in ['G','EI','B'])
        fine=[];coarse=[];diagnostics=[]
        for seed in SEEDS:
            a,b,diag=bracket(state,17,seed,15);fine.append(a);coarse.append(b);diagnostics.append(diag)
        EE=np.mean(fine,axis=0);small=np.mean(coarse,axis=0)
        values=[transfer(state,EE,size) for size in [3,5,7]]
        qc=transfer(state,small,7);individual=[transfer(state,a,7) for a in fine]
        errors={key:dict(quadrature=abs(qc[key]/values[-1][key]-1),basis=abs(values[-2][key]/values[-1][key]-1),scramble_relative_std=float(np.std([v[key] for v in individual],ddof=1)/values[-1][key])) for key in ['K','tau']}
        np.savez_compressed(OUT/f'cell-{index}.npz',**state,EE=EE,coarse_EE=small,replicate_EE=np.array(fine),replicate_coarse_EE=np.array(coarse),**{'fine_'+key:v for key,v in values[-1].items()})
        rows.append(dict(cell=index,equilibrium_quadrature=eqerr,errors=errors,K=values[-1]['K'],tau=values[-1]['tau'],native_ratio=values[-1]['K']/state['native_K'],min_projected_collision_rate=values[-1]['min_rate'],events=4*2**17,diagnostics=diagnostics))
        print('CELL',index,'K_RATIO',rows[-1]['native_ratio'],'ERRORS',json.dumps(errors),flush=True)
    numerical=all(row['equilibrium_quadrature']<1e-6 and all(e['quadrature']<.03 and e['basis']<.02 for e in row['errors'].values()) for row in rows)
    result=dict(classification='Counterexample candidate',rows=rows,numerical_gates_passed=numerical,
        native_compatibility_passed=all(abs(r['native_ratio']-1)<.2 for r in rows),
        elapsed_seconds=time.monotonic()-start,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        full_physical_collision_error_certified=False,whole_star_closure_replaced=False,full_dynamic_charge_solved=False)
    write(OUT/'result.json',result);print('RESULT',numerical,'SECONDS',result['elapsed_seconds'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
