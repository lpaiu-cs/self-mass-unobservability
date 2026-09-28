"""Finite-jump, angle-dependent free-electron reference for the frozen LTE tangent.

Sazonov & Sunyaev (2000), Eq. 7, with measured detailed-balance symmetrization.
This replaces BOTH native elastic scattering and Kompaneets; it is not a
reconstruction of ATOMIC's collective/bound-electron scattering mixture.
"""
from pathlib import Path
import argparse
import gc
import json
import resource
import signal
import time
import numpy as np
from numpy.polynomial.legendre import leggauss, legvander
from scipy.integrate import quad_vec
from scipy.sparse import csr_matrix, diags, eye
from scipy.sparse.linalg import splu
import def_photon_spatial_coupling as old

OUT=old.OUT.parent/'def-photon-finite-jump'
SOURCE='https://arxiv.org/html/astro-ph/9910280v2'


def kernel(u,v,mu,theta):
    """2*pi*K per outgoing u per scattering cosine, source Eq. 7a-d."""
    d=1-mu;g=np.sqrt((v-u)**2+2*u*v*d)
    e=np.sqrt(2*d)/g*(v-u+theta*u*v*d)
    exponent=e*e/(4*d*theta)
    bracket=(1+mu*mu+(1/8-mu-63*mu*mu/8+5*mu**3)*theta
        -mu*(1+mu)/2*e*e-3*(1+mu*mu)/(32*d*d)*e**4/theta
        +mu*(1-mu*mu)*e*theta*u+(1+mu*mu)/(8*d)*e**3*u
        +d*d*theta*theta*u*v)
    # The asymptotic polynomial is not used in its exponentially remote tail.
    assert np.all(bracket[exponent<=64]>=0) if np.ndim(bracket) else bracket>=0 or exponent>64
    return 3/16*np.sqrt(2/np.pi/theta)*v/(u*g)*np.maximum(bracket,0)*np.exp(-exponent)*(exponent<=64)


def moments(theta):
    """Independent adaptive outgoing-frequency integrals of the source kernel."""
    x,w=leggauss(48);t=(x+1)/2;mu=1-2*t*t;angular_w=w*2*t
    rows=[]
    for u in [.1,1.,5.,15.,50.]:
        def f(z):
            v=u*np.exp(2*np.sqrt(theta)*z);dv=2*np.sqrt(theta)*v
            k=kernel(u,v,mu,theta)*angular_w*dv
            return np.array([k.sum(),k.sum()*(v-u)/u,k.sum()*((v-u)/u)**2])
        val,err=quad_vec(f,-9,9,epsabs=1e-12,epsrel=1e-10,points=[0.])
        rows.append(dict(u=u,moments=val.tolist(),quadrature_error=float(err),
            normalization_error=float(abs(val[0]-(1-2*theta*u))),
            drift_relative=float(abs(val[1]/(theta*(4-u))-1)),
            diffusion_relative=float(abs(val[2]/(2*theta)-1))))
    return rows


def prepare():
    assert not (OUT/'plan.json').exists();OUT.mkdir(exist_ok=True)
    paths=[Path(__file__),Path(old.__file__),old.OUT/'bank.npz',old.OUT/'result.json',old.OUT/'mode-1.npz']
    old.ex.write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='892ea0265',
        claim='Remove the small-jump expansion from the actual frozen surface photon/material response, retaining joint scattering-angle/frequency redistribution, induced scattering, recoil and the same absorptive EOS tangent.',
        decision='Compare the finite-jump response against diffusion and elastic controls at the same free-electron density, before committing to atmosphere/GR evolution with an unvalidated diffusion approximation.',
        source=SOURCE,model='Eq.7 free Maxwellian electrons; unpolarized intensity. Pairwise geometric symmetrization of the equilibrium transition flux. Both native scalar elastic scattering and Kompaneets are replaced, not added. Native absorption, matter heat capacity, temperature, density and spatial scale unchanged.',
        cells='Same 24134 fixed-opacity cells through u=60. Exact LTE bin capacities define positive quadrature weights; bin representative energies define finite photon number. Compare one/two quadrature nodes inside the SAME source cells, not a new opacity fit.',
        angular='Normalized Legendre modes. K_l=integral K P_l. Positive angular quadrature gives |K_l|<=K_0. Streaming is the existing symmetric Galerkin operator.',
        range='Pairs within |log(v/u)|<=12 sqrt(theta), about six widest Doppler sigmas; exponent>64 omitted. This cutoff is tested against 16 sqrt(theta) on continuum moments, not claimed as a rigorous physical tail envelope.',
        pilot=dict(seconds=60,angular_modes=4,scattering_quadrature=16,steps=16),
        production=dict(seconds=240,angular_modes=[8,12],scattering_quadrature=[24,48],steps=[16,32,64],frequency_quadrature=[1,2],kH=1),
        gates=dict(time_order=1.8,time_difference_initial=.001,angular_difference_initial=.001,
            frequency_quadrature_difference_initial=.001,scattering_quadrature_difference_initial=.001,
            energy_balance=1e-9,energy_equation_residual=1e-9,entropy_growth=1e-10,
            moment_normalization=1e-7,moment_relative=.005,heating_relative=.005,solver_relative=1e-11),
        budget=dict(CPU_workers=1,GPU=False,new_EOS_calls=0,new_physical_queries=0,new_stellar_steps=0,
            memory_cap_GB=4,automatic_expansion=False,forecast='Production assessed after the four-mode pilot; stop on cap or failed algebra/solver checks.'),
        limits='Free-electron reference, not native collective/Rayleigh/bound-electron certification. Actual-temperature opacity, atmosphere, moving matter/GR feedback and full charge remain incomplete.',
        bindings={p.relative_to(old.h.ROOT).as_posix():old.h.digest(p) for p in paths}))
    print('PREPARED',OUT,flush=True)


def coefficients(bank,order,nangle,nfreq=1,width=12):
    start=time.monotonic();u=bank['u'][bank['u']<=60];C=bank['Ci'][:len(u)];N=len(u)
    theta=float(bank['theta']);B=15*float(bank['arad'])*float(bank['T'])**3/np.pi**4
    stop=np.searchsorted(u,u*np.exp(width*np.sqrt(theta)),side='right')
    counts=stop-np.arange(N);i=np.repeat(np.arange(N,dtype=np.int32),counts)
    j=np.concatenate([np.arange(a,b,dtype=np.int32) for a,b in enumerate(stop)])
    x,w=leggauss(nangle);t=(x+1)/2;mu=1-2*t*t;aw=2*w*t;P=legvander(mu,order-1)
    # Exact bin capacity fixes the positive integration measure, not a mean fit.
    gx,gw=leggauss(nfreq);gw=gw/2;edges=bank['edges_u'][:N+1]
    nodes=u[:,None]+np.diff(edges)[:,None]*gx/2
    occ=1/np.expm1(nodes);f=nodes**4*occ*(1+occ)
    weights=C[:,None]/(B*(f@gw))[:,None]*gw
    a=np.zeros((len(i),order));max_balance=0.;weighted_balance=0.;total=0.
    for begin in range(0,len(i),40000):
        end=min(begin+40000,len(i));ii=i[begin:end];jj=j[begin:end];block=a[begin:end]
        for q in range(nfreq):
            for r in range(nfreq):
                v=nodes[ii,q];z=nodes[jj,r];nv=occ[ii,q];nz=occ[jj,r]
                weight=B*float(bank['rate_e'])*weights[ii,q]*weights[jj,r]
                for cosine,ww,pl in zip(mu,aw,P):
                    forward=v*v*nv*(1+nz)*kernel(v,z,cosine,theta)
                    reverse=z*z*nz*(1+nv)*kernel(z,v,cosine,theta)
                    geom=np.sqrt(forward*reverse);val=weight*ww*geom
                    block+=val[:,None]*pl
                    diff=abs(forward-reverse);scale=forward+reverse
                    valid=scale>1e-200
                    max_balance=max(max_balance,float(np.max(diff[valid]/scale[valid],initial=0)))
                    weighted_balance+=float(np.sum(weight*ww*diff));total+=float(np.sum(weight*ww*scale))
    assert np.all(a[:,0]>=0) and np.max(abs(a)-a[:,0,None])<1e-12
    info=dict(seconds=time.monotonic()-start,pairs=len(i),frequency_cells=N,nangle=nangle,nfreq=nfreq,
        balance_relative_max=max_balance,balance_relative_weighted=weighted_balance/total,
        interpretation='Geometric equilibrium-flux symmetrization; measured change, not exact relativistic certification.')
    return i,j,a,info


class Operator(old.Operator):
    def __init__(self,bank,order,kH,nangle=24,nfreq=1):
        super().__init__(dict(bank,rate_s=np.zeros_like(bank['rate_s'])),order,kH,compton=False)
        i,j,a,self.info=coefficients(bank,order,nangle,nfreq);same=i==j;other=~same
        u=self.u;root=np.sqrt(self.Ci);Cm=self.Cm;N=len(u);cols=np.arange(N)*order
        ip=i[other];jp=j[other];a0=a[other,0];du=u[jp]-u[ip]
        d=np.bincount(ip,a0*u[ip]**2/self.Ci[ip],minlength=N)+np.bincount(jp,a0*u[jp]**2/self.Ci[jp],minlength=N)
        diagonal=np.repeat(d[:,None],order,axis=1)
        diagonal+=(a[same,0,None]-a[same])*u[:,None]**2/self.Ci[:,None]
        q0=(np.bincount(ip,a0*du*u[ip]/root[ip],minlength=N)
            -np.bincount(jp,a0*du*u[jp]/root[jp],minlength=N))/np.sqrt(Cm)
        self.q[cols]+=q0;self.aa+=float(a0@(du*du)/Cm)
        self.jump_q=q0;self.jump_aa=float(a0@(du*du)/Cm)
        self.diagonal=diagonal;self.gain=[]
        rows=np.r_[ip,jp];columns=np.r_[jp,ip]
        for ell in range(order):
            v=-a[other,ell]*u[ip]*u[jp]/(root[ip]*root[jp])
            self.gain.append(csr_matrix((np.r_[v,v],(rows,columns)),shape=(N,N)))
        self.L=self.L+diags(diagonal.ravel(),format='csc');self.P=self.L+1j*self.kc*self.V
        self.off_bound=max(float(np.max(np.asarray(abs(g).sum(axis=1)))) for g in self.gain)
        self.schur=float(self.aa-np.sum(self.q[cols]**2/self.L.diagonal()[cols]))
        assert self.schur>=-1e-10
        self.max_iterations=0;self.max_solver_error=0.
        # Store physical photon collision force through the actual operator.
        self.rate_s=np.full(N,float(bank['rate_e']))

    def off(self,E):
        z=E.reshape(-1,self.order);return np.column_stack([g@z[:,ell] for ell,g in enumerate(self.gain)]).ravel()

    def rhs(self,T,E):
        a,z=super().rhs(T,E);return a,z-self.off(E)

    def solver(self,dt):
        factor=splu(eye(self.size,format='csc')+dt*self.P);q=dt*self.q
        response=factor.solve(q.astype(complex));denom=1+dt*self.aa-q@response
        contraction=dt*self.off_bound
        assert contraction<.5,('fixed-point budget',contraction)
        def base(T,E):
            y=factor.solve(E);x=(T-q@y)/denom;return x,y-response*x
        def solve(T,E):
            U,W=base(T,E);scale=max(np.sqrt(abs(T)**2+np.vdot(E,E).real),1e-100)
            for iteration in range(1,17):
                V,Z=base(T,E-dt*self.off(W))
                delta=np.sqrt(abs(V-U)**2+np.linalg.norm(Z-W)**2)
                error=contraction/(1-contraction)*delta/scale
                U,W=V,Z
                if error<1e-12:break
            else:raise RuntimeError('finite-jump solve did not converge')
            f,g=self.rhs(U,W);res=np.sqrt(abs(U-dt*f-T)**2+np.linalg.norm(W-dt*g-E)**2)/scale
            assert res<1e-11,res
            self.max_iterations=max(self.max_iterations,iteration);self.max_solver_error=max(self.max_solver_error,float(res))
            return U,W
        return solve


def pilot():
    assert not (OUT/'pilot.json').exists();signal.alarm(60);start=time.monotonic()
    b=dict(np.load(old.OUT/'bank.npz'));op=Operator(b,4,1,16);build=time.monotonic()-start
    r=op.evolve(float(b['H']/b['c']),16)
    result=dict(classification='Counterexample candidate',build_seconds=build,seconds=time.monotonic()-start,
        peak_memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
        coefficients=op.info,off_norm_bound=op.off_bound,preconditioner_schur=op.schur,
        iterations=op.max_iterations,solver_residual=op.max_solver_error,
        balance=r['balance'],energy_residual=r['energy_equation_residual'],entropy_growth=r['entropy_growth'])
    old.ex.write(OUT/'pilot.json',result);print('PILOT',result,flush=True);signal.alarm(0)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','pilot']);globals()[p.parse_args().action]()
