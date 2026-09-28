"""Multi-charge Yukawa HNC structure for finite-temperature electron transport.

No degenerate-fit cutoff, ion-charge averaging, or conductivity calibration.
Classical ions, static compressibility screening, HNC without a bridge term,
and elastic Born scattering remain declared physical approximations.
"""
from pathlib import Path
import argparse
import json
import time
import resource
import urllib.request
import numpy as np
from scipy.fft import dst
from scipy.special import erfc
from scipy.interpolate import PchipInterpolator,CubicSpline
from numpy.polynomial.legendre import leggauss,legvander
import def_conduction_core as previous

ex=previous.ex
base=ex.model.base
h=ex.h
OUT=previous.OUT.parent/'def-ionic-structure-transport'
CELLS=[1722,2000,2500,2800,3000,3500,4000,4122,4928,5734]


class Transform:
    def __init__(self,n,length):
        self.n=n;self.length=length;self.dr=length/(n+1);self.dk=np.pi/length
        self.r=np.arange(1,n+1)*self.dr;self.q=np.arange(1,n+1)*self.dk

    def forward(self,f):
        shape=(-1,)+(1,)*(f.ndim-1)
        return 2*np.pi*self.dr/self.q.reshape(shape)*dst(self.r.reshape(shape)*f,type=1,axis=0)

    def inverse(self,F):
        shape=(-1,)+(1,)*(F.ndim-1)
        return self.dk/(4*np.pi**2*self.r.reshape(shape))*dst(self.q.reshape(shape)*F,type=1,axis=0)


def parameters(index):
    d,x=base.inputs();ne=x['ne'][index];T=x['T'][index];ae=(3/(4*np.pi*ne))**(1/3)
    charges=np.unique(base.thermal.g.c.Z)
    densities=np.array([x['ion'][index,base.thermal.g.c.Z==z].sum() for z in charges])*ae**3
    keep=densities>0;charges=charges[keep];densities=densities[keep]
    distribution=ex.model.distribution(np.array([T]),np.array([ne]),256)
    kappa2=base.elementary_charge**2/base.epsilon_0*distribution[5][0]*ae**2
    gamma=base.E2/(base.k*T*ae)
    return dict(index=index,ae=ae,T=T,ne=ne,charges=charges,densities=densities,gamma=gamma,kappa=np.sqrt(kappa2),
                T_over_TF=float(x['theta'][index]),charge_deficit=float(x['ionization_deficit'][index]))


def hnc(state,n=1024,length=32.,strength=1.):
    start=time.monotonic();tr=Transform(n,length);r=tr.r;q=tr.q;kap=state['kappa'];Z=state['charges'];rho=state['densities']
    charge=Z[:,None]*Z[None,:]*state['gamma']*strength;root=np.sqrt(rho);m=len(Z)
    potential=np.exp(-kap*r)/r
    # Gaussian-smoothed Yukawa tail has an exact Fourier transform. The
    # remaining direct correlation is regular at the Coulomb singularity.
    long=np.exp(kap*kap/4)/(2*r)*(np.exp(-kap*r)*erfc(kap/2-r)-np.exp(kap*r)*erfc(kap/2+r))
    long_q=4*np.pi*np.exp(-q*q/4)/(q*q+kap*kap)
    V=potential[:,None,None]*charge;Vl=long[:,None,None]*charge
    gamma=np.zeros((n,m,m));trace=[];eye=np.eye(m)
    for iteration in range(600):
        exponent=-V+gamma
        assert np.max(exponent)<50,('closure overflow',iteration,np.max(exponent))
        short=np.expm1(exponent)-gamma+Vl
        C=tr.forward(short)-long_q[:,None,None]*charge
        C=(C+C.swapaxes(1,2))/2
        matrix=eye-C*root[None,:,None]*root[None,None,:]
        Cr=C*root[None,None,:]
        indirect_q=Cr@np.linalg.solve(matrix,Cr.swapaxes(1,2))
        new=tr.inverse(indirect_q);error=float(np.max(abs(new-gamma)))
        if iteration%20==0:trace.append([iteration,error])
        if error<1e-8:break
        gamma+=.2*(new-gamma)
    else:raise RuntimeError(('HNC iteration limit',state['index'],n,error))
    # q=0 evaluated from the same regular integral, not an extrapolated fit.
    C0=4*np.pi*tr.dr*np.einsum('n,nij->ij',r*r,short)-4*np.pi/kap**2*charge
    allC=np.concatenate([C0[None],C]);M=eye-allC*root[None,:,None]*root[None,None,:]
    minimum=float(np.linalg.eigvalsh(M).min());assert minimum>0
    z=root*Z;structure=np.einsum('i,ni->n',z,np.linalg.solve(M,np.broadcast_to(z,(len(M),m))[...,None])[...,0])/(z@z)
    assert np.min(structure)>0
    return dict(q=np.r_[0,q],Scharge=structure,r=r,log_g=exponent,
        iterations=iteration+1,residual=error,minimum_inverse_structure_eigenvalue=minimum,trace=np.array(trace),seconds=time.monotonic()-start,
        gamma=state['gamma'],kappa=kap,charges=Z,densities=rho)


def kinetic(state,structure,order=256,angular_order=64):
    T=state['T'];ne=state['ne'];eta,x,p,fp,dx,_,_=ex.model.distribution(np.array([T]),np.array([ne]),order)
    eta=float(eta[0]);x=x[0];p=p[0];fp=fp[0];dx=dx[0];theta=base.k*T/(base.m_e*base.c**2);E=np.sqrt(1+p*p)
    physical=base.m_e*base.c*p;velocity=base.c*p/E
    qe2=state['kappa']**2/state['ae']**2
    u=base.hbar**2*qe2/(4*physical**2);top=np.log1p(1/u);gx,gw=leggauss(angular_order)
    y=top[:,None]*(gx+1)/2;transfer=u[:,None]*np.expm1(y)
    q=2*physical[:,None]/base.hbar*state['ae']*np.sqrt(transfer)
    interp=PchipInterpolator(structure['q'],structure['Scharge'],extrapolate=False)
    S=interp(np.minimum(q,structure['q'][-1]));S=np.where(q>structure['q'][-1],1.,S)
    coulomb=top/4*np.sum(gw*(-np.expm1(-y))*(1-(velocity/base.c)[:,None]**2*transfer)*S,axis=1)
    charge2=np.sum(state['densities']*state['charges']**2)/state['ae']**3
    nu=4*np.pi*base.E2**2*charge2/(physical**2*velocity)*coulomb
    measure=theta*E*p**3*fp*dx/(3*np.pi**2);basis=legvander((x-eta)/4,ex.DEGREE)
    G=np.einsum('ni,nj,n->ij',basis,basis,measure)
    EI=np.einsum('ni,nj,n->ij',basis,basis,measure*nu/ex.RATE)
    B=np.c_[basis.T@(measure/E),basis.T@(measure*(x-eta)/E)]
    original=ex.equilibrium(state['index'],order)
    gerror=float(np.linalg.norm(G-original['G'])/np.linalg.norm(G));berror=float(np.linalg.norm(B-original['B'])/np.linalg.norm(B))
    assert max(gerror,berror)<1e-12
    original.update(G=G,EI=EI,B=B)
    weight=physical**3/(3*np.pi**2*base.hbar**3*base.m_e*E)*fp*dx/nu
    mean=np.sum(weight*(x-eta))/np.sum(weight)
    diagonal=base.k*base.k*T*np.sum(weight*(x-eta-mean)**2)
    return original,dict(coulomb_min=float(coulomb.min()),direct_EI_K=float(diagonal),
        tail_structure_last=float(structure['Scharge'][-1]),G_relative=gerror,B_relative=berror)


def controls():
    tr=Transform(512,32);f=np.exp(-tr.r**2);exact=np.pi**1.5*np.exp(-tr.q**2/4)
    transformed=tr.forward(f);error=float(np.max(abs(transformed-exact))/exact.max())
    roundtrip=float(np.max(abs(tr.inverse(transformed)-f)));assert max(error,roundtrip)<1e-12
    state=parameters(3000);zero=hnc(state,128,32,strength=0.);assert np.array_equal(zero['Scharge'],np.ones(129))
    weak=hnc(state,512,32,strength=1e-4);q=weak['q'];qion=4*np.pi*state['gamma']*1e-4*np.sum(state['densities']*state['charges']**2)
    reference=(q*q+state['kappa']**2)/(q*q+state['kappa']**2+qion)
    weak_error=float(np.max(abs(weak['Scharge']-reference)));assert weak_error<1e-5
    return dict(classification='Counterexample candidate',passed=True,Gaussian_transform_error=error,roundtrip_error=roundtrip,zero_coupling_exact=True,weak_RPA_absolute_error=weak_error)


def prepare():
    assert not OUT.exists();OUT.mkdir()
    urls={'shaffer-starrett2020.html':'https://arxiv.org/html/2002.00928','blouin2020.html':'https://arxiv.org/html/2006.16390','ionic-correlations2017.html':'https://arxiv.org/html/1707.01509'}
    for name,url in urls.items():
        with urllib.request.urlopen(url,timeout=20) as response:(OUT/name).write_bytes(response.read())
    paths=[Path(__file__),Path(ex.__file__),Path(base.__file__),base.thermal.OUT/'coefficients.npz',previous.OUT/'manifest.json']+[OUT/name for name in urls]
    ex.write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='2ad1201f',bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths},
        sources=urls,cells=CELLS,
        claim='Remove the strong-degeneracy cutoff from the ion structure used by the microscopic conduction operator, using the actual multi-charge mixture and finite-T electron compressibility. Construct the paired ee+ei response in the partially degenerate and classical fully ionized layers.',
        model='Classical multi-species Yukawa HNC with no bridge function, ideal-electron compressibility screening, static ion structure and elastic Born/Mott electron-ion scattering. Every present nuclear charge group is retained. Same Pauli-blocked dynamically screened leading ee collision bracket as Phase55. This is a declared microscopic candidate, not the full average-atom mean-force model of the references.',
        limits='No neutral/partially ionized opacity, ion dynamics, full finite-q electron screening, non-Born phase shifts, coherent ee-channel interference or continuum physical error certificate. Do not graft a new hard zero-flux edge onto a stellar path.',
        grids=[[512,32],[1024,32]],box_controls=dict(cells=[1722,3000,5734],grid=[1536,48]),
        gates=dict(HNC_residual=1e-8,positive_structure=True,EI_grid_relative=.005,EI_box_relative=.005,energy_angular_relative=.0001,ee_quadrature=.03,basis=.02),
        budget=dict(pilot_states=1,HNC_iteration_limit=600,production_hard_seconds=180,workers=1,BLAS_threads=1,native_calls=0,new_stellar_steps=0,automatic_expansion=False)))
    pilot()


def pilot():
    assert not (OUT/'pilot.json').exists()
    ex.write(OUT/'controls.json',controls());s=parameters(4122);a=hnc(s,512,32);b=hnc(s,1024,32)
    for name,values in [('coarse',a),('fine',b)]:np.savez_compressed(OUT/f'pilot-{name}.npz',**values)
    ex.write(OUT/'pilot.json',dict(HNC_coarse_seconds=a['seconds'],HNC_fine_seconds=b['seconds'],iterations=[a['iterations'],b['iterations']],
        forecast_seconds=(a['seconds']+b['seconds'])*10+b['seconds']*1.5*3+7*6,
        assumption='Measured two HNC grids at the old cutoff; ten states plus three longer boxes. Seven new ee states at 6 seconds each based on Phase56. Import, source download and output overhead excluded; stop at 180 seconds.'))
    print('PILOT',a['seconds'],b['seconds'],'ITERATIONS',a['iterations'],b['iterations'],flush=True)


def symbolic():
    import sympy as sp
    q,e,i=sp.symbols('q e i',positive=True)
    C=-4*sp.pi/(q*q+e);rho=i/(4*sp.pi)
    assert sp.factor(1/(1-rho*C)-(q*q+e)/(q*q+e+i))==0
    return dict(classification='Proven',passed=True,
        OZ='For D=diag(sqrt(n_s)), M=I-D C D and H=C+C D M^-1 D C. If M is positive definite, S_Z=(DZ)^T M^-1(DZ)/sum(n_s Z_s^2)>0.',
        weak_limit='At weak ionic coupling the rank-one Yukawa direct correlation yields S_Z(q)=(q^2+kappa_e^2)/(q^2+kappa_e^2+kappa_i^2).',
        collision='The positive charge structure times the positive screened elastic angular kernel gives a positive EI Gram matrix. Add the conserving positive EE bracket before zero-current projection, not independent fitted thermal resistivities.',
        limit='These algebraic identities neither prove HNC accuracy nor certify a continuum collision gap.')


def run():
    assert not (OUT/'result.json').exists()
    plan=json.loads((OUT/'plan.json').read_text())
    for p,sha in plan['bindings'].items():assert h.digest(h.ROOT/p)==sha,p
    assert json.loads((OUT/'pilot.json').read_text())['forecast_seconds']<150
    ex.write(OUT/'symbolic.json',symbolic());start=time.monotonic();rows=[]
    for index in CELLS:
        state=parameters(index);assert abs(state['charge_deficit'])<1e-4
        structures=[]
        for n,L in plan['grids']:
            path=OUT/f'hnc-{index}-{n}-{L}.npz'
            if index==4122:
                structure=dict(np.load(OUT/f"pilot-{'coarse' if n==512 else 'fine'}.npz"))
            else:structure=hnc(state,n,L)
            np.savez_compressed(path,**structure);structures.append(structure)
        fine,diag=kinetic(state,structures[1]);zero=np.zeros_like(fine['EI']);ei=ex.transfer(fine,zero,7)
        comparisons={}
        def compare(name,other):
            answer=ex.transfer(other,zero,7)
            comparisons[name]={k:float(abs(answer[k]/ei[k]-1)) for k in ['K','tau']}
        compare('structure_grid',kinetic(state,structures[0])[0])
        compare('energy',kinetic(state,structures[1],128)[0])
        compare('angle',kinetic(state,structures[1],256,128)[0])
        if index in plan['box_controls']['cells']:
            n,L=plan['box_controls']['grid'];larger=hnc(state,n,L)
            np.savez_compressed(OUT/f'hnc-{index}-{n}-{L}.npz',**larger)
            compare('box',kinetic(state,larger)[0])
        reused=index in [4122,4928,5734]
        if reused:
            old=np.load(previous.OUT/f'cell-{index}.npz')
            assert fine['T']==old['T'] and fine['ne']==old['ne']
            EE=old['EE'];coarse=old['coarse_EE'];replicates=old['replicate_EE'];table_error=0.
            event_checks={'reuse':True};extra={}
        else:
            m=previous.microscopic
            maximum=np.sqrt(1-1/(1+2*fine['theta']*(fine['eta']+40))**2)*(1+1e-12)
            grid=m.nodes(fine,maximum);values=np.array([m.screen.polarization(a,fine) for a in grid])
            mids=(grid[1:]+grid[:-1])/2;vf=fine['xF']/np.sqrt(1+fine['xF']**2)
            selected=np.unique(np.r_[np.linspace(0,255,16).astype(int),np.argsort(abs(mids-vf))[:16]])
            reference=np.array([m.screen.polarization(a,fine) for a in mids[selected]])
            table_error=float(np.max(abs(reference-CubicSpline(grid,values)(mids[selected])))/fine['qs2'])
            ex.amplitude=m.leading.amplitude(fine,grid,values);replicates=[];coarse_rows=[];event_checks=[]
            for seed in ex.SEEDS:
                a,b,checks=ex.bracket(fine,17,seed,15);replicates.append(a);coarse_rows.append(b);event_checks.append(checks)
            EE=np.mean(replicates,axis=0);coarse=np.mean(coarse_rows,axis=0);extra=dict(phase=grid,polarization=values)
        answer=ex.transfer(fine,EE,7);small=ex.transfer(fine,coarse,7);low=ex.transfer(fine,EE,5)
        errors={k:dict(quadrature=float(abs(small[k]/answer[k]-1)),basis=float(abs(low[k]/answer[k]-1))) for k in ['K','tau']}
        np.savez_compressed(OUT/f'cell-{index}.npz',**fine,EE=EE,coarse_EE=coarse,replicate_EE=replicates,
            **extra,**{'fine_'+k:v for k,v in answer.items()})
        numerical=all(v<.005 for key in ['structure_grid','box'] if key in comparisons for v in comparisons[key].values())
        numerical=numerical and all(v<1e-4 for key in ['energy','angle'] for v in comparisons[key].values())
        numerical=numerical and all(v['quadrature']<.03 and v['basis']<.02 for v in errors.values()) and table_error<1e-5
        row=dict(index=index,passed=numerical,T_over_TF=state['T_over_TF'],charge_deficit=state['charge_deficit'],
            HNC_iterations=[int(s['iterations']) for s in structures],comparisons=comparisons,errors=errors,
            K=answer['K'],tau=answer['tau'],native_ratio=answer['K']/fine['native_K'],EI_only_K=ei['K'],
            direct_EI_relative=float(abs(ei['K']/diag['direct_EI_K']-1)),table_error=table_error,
            diagnostic=diag,reused_EE=reused,event_checks=event_checks)
        rows.append(row);ex.write(OUT/'progress.json',rows)
        print('CELL',index,'PASS',numerical,'K',row['K'],'RATIO',row['native_ratio'],'ERRORS',errors,flush=True)
        assert time.monotonic()-start<175,'production time budget'
    result=dict(classification='Counterexample candidate',passed=all(r['passed'] for r in rows),checks=rows,
        seconds=time.monotonic()-start,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        fully_ionized_sample_operator_constructed=True,physical_error_certified=False,
        whole_star_heat_closed=False,full_dynamic_charge_solved=False,new_stellar_steps=0)
    ex.write(OUT/'result.json',result);print('RESULT',result['passed'],'SECONDS',result['seconds'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','pilot','run']);globals()[p.parse_args().action]()
