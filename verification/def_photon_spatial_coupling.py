"""Coupled spectral/angular LTE tangent with material recoil and spatial transport.

The declared candidate uses dipole elastic scattering and the angle-averaged,
nonrelativistic free-electron Kompaneets operator. The native scalar scattering
column does NOT identify this gain kernel or certify its physical error.
"""
from pathlib import Path
from types import FunctionType
import argparse
import json
import signal
import time
import urllib.request
import numpy as np
import sympy as sp
from scipy.sparse import coo_matrix, diags, eye
from scipy.sparse.linalg import splu
import def_photon_native_measure as previous

ex=previous.ex;h=previous.h;base=previous.matter.model.base
OUT=previous.OUT.parent/'def-photon-spatial-coupling'
SOURCE='https://arxiv.org/html/2209.06240'


def prepare():
    assert not (OUT/'plan.json').exists();OUT.mkdir(exist_ok=True)
    d,p=base.inputs();old=json.loads((previous.native.OUT/'requests.json').read_text())[2]
    T=float(base.k*p['T'][0]/1.602176634e-16)
    row=dict(old,name='phActual61',fields=dict(old['fields'],mixname='phActual61',temps=f'{T:.15g}',datype='cont'))
    ex.write(OUT/'requests.json',[row]);(OUT/'lanl-tops-form.html').write_bytes((previous.native.OUT/'lanl-tops-form.html').read_bytes())
    paths=[Path(__file__),Path(previous.__file__),previous.OUT/'manifest.json',previous.matter.old.BRIDGE,
           base.thermal.OUT/'coefficients.npz',OUT/'requests.json',OUT/'lanl-tops-form.html']
    ex.write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='250385ee1',
        claim='Connect material temperature, resolved absorption, angular scattering, photon frequency redistribution and spatial transport in one unsplit linear operator at the actual saved surface state.',
        previous_turn='Progress: native photon grid, conservative local absorptive tangent and failed-window replay accepted in Phase60.',
        sources=[SOURCE,'https://arxiv.org/html/1601.01005'],
        actual_temperature_keV=T,actual_temperature_K=float(p['T'][0]),
        model='Vacuum photons; fixed local state and Fourier spatial modes. Native absorption and scalar scattering rates. Elastic unpolarized dipole phase eigenvalues 1,0,1/10,0,...; this is an explicit leading-order candidate, not a recovered ATOMIC differential kernel. Angle-averaged Kompaneets energy redistribution uses the same material EOS free-electron density and exchanges energy with that material. Stimulated scattering and photon-number conservation retained. No fitted timescale or opacity normalization.',
        frequency='Retain native source cells through u=60 and their exact Planck derivative capacities; report excluded capacity. Adjacent frequency exchange uses a positive detailed-balance bracket. Tail and frequency-discretization checks remain distinct from physical cross-section uncertainty.',
        spatial='Legendre angular projection of exact local Fourier streaming. Wavenumbers 0,1/H_T,10/H_T with H_T measured between the two actual outer cells in proper distance. These are local component modes, not the inhomogeneous atmosphere or a derived orbital drive.',
        gates=dict(native_means_relative=.001,energy_number_algebra=1e-10,linear_residual=1e-9,
                   time_order=1.8,time_difference_initial=.001,angular_difference_initial=.001,tail_capacity_fraction=1e-17),
        budget=dict(new_physical_queries=1,network_seconds=45,native_EOS_calls=1,
                    pilot_seconds=30,production_seconds=180,CPU_workers=1,GPU=False,new_stellar_steps=0,automatic_expansion=False),
        limits='This candidate does not close native angular/energy redistribution, polarization, bound-electron recoil, collective corrections to Compton redistribution, nonlinear opacity updates, the atmosphere, GR feedback or full dynamic charge.',
        bindings={p.relative_to(h.ROOT).as_posix() if p.is_relative_to(h.ROOT) else str(p):h.digest(p) for p in paths}))


def fetch():
    # TOPS does not interpolate temperature. Enforce its published input domain.
    from bs4 import BeautifulSoup
    soup=BeautifulSoup((OUT/'lanl-tops-form.html').read_bytes(),'html.parser')
    available=[float(o.get('value',o.get_text(strip=True))) for o in soup.find('select',{'name':'tlow'}).find_all('option')]
    requested=[float(r['fields']['temps']) for r in json.loads((OUT/'requests.json').read_text())]
    assert all(T in available for T in requested),'TOPS requires database temperatures; use a declared interpolation model'
    assert not (OUT/'retrieval.json').exists();start=time.monotonic();signal.alarm(45)
    try:
        with urllib.request.urlopen(SOURCE,timeout=15) as r:(OUT/'jiang2022.html').write_bytes(r.read())
        env=dict(vars(previous.native.retrieval),OUT=OUT)
        env['save']=lambda name,value:ex.write(OUT/name,value)
        fn=FunctionType(previous.native.retrieval.fetch.__code__,env,argdefs=previous.native.retrieval.fetch.__defaults__)
        result=fn(json.loads((OUT/'requests.json').read_text())[0])
        ex.write(OUT/'retrieval.json',dict(result=result,seconds=time.monotonic()-start))
        print('RETRIEVAL',result,flush=True)
    finally:signal.alarm(0)


def bank():
    assert not (OUT/'bank.npz').exists();row=json.loads((OUT/'requests.json').read_text())[0]
    T=float(row['fields']['temps']);sources=[np.load(previous.OUT/(name+'-bank.npz')) for name in ['surf1000T0015','surf1000T002']]
    # Exact common rational lattice for these two source edge grids. Merging
    # floating almost-duplicates would otherwise introduce artificial tiny bins.
    integer=[np.rint(s['edges_keV']*3200000).astype(np.int64) for s in sources]
    assert all(np.max(abs(i/3200000-s['edges_keV']))<1e-12 for i,s in zip(integer,sources))
    energy_edges=np.unique(np.concatenate(integer))/3200000
    centers=(energy_edges[:-1]+energy_edges[1:])/2;edges=energy_edges/T;u=centers/T
    weights=previous.loss.weights(edges,1);fraction=np.log(T/.0015)/np.log(.002/.0015)
    assert 0<fraction<1
    fields={}
    for key in ['absorption','scattering']:
        samples=[s[key][np.clip(np.searchsorted(s['edges_keV'],centers,side='right')-1,0,len(s[key])-1)] for s in sources]
        fields[key]=np.exp((1-fraction)*np.log(samples[0])+fraction*np.log(samples[1]))
        fields[key+'_source_states']=np.array(samples)
    ka=fields['absorption'];ks=fields['scattering']
    means=np.array([weights[:,0]@ka,1/(weights[:,1]@(1/(ka+ks)))])
    d,p=base.inputs();gas=previous.matter.old.GasEOS();T_K=float(p['T'][0])
    eos=gas(2,float(d['lnd'][0]),float(d['lnT'][0]),d['X'][0]);Cm=np.exp(d['lnd'][0])*eos[10]/T_K
    Ci=4*gas.a_rad*T_K**3*weights[:,1];rho=float(row['fields']['dens'])
    ne=float(base.Avogadro*eos[13]*1e6);sigma=8*np.pi/3*(base.E2/(base.m_e*base.c**2))**2
    rate_e=ne*sigma*base.c;theta=base.k*T_K/(base.m_e*base.c**2)
    distance=abs(d['radius_cm'][1]-d['radius_cm'][0])*float(d['A'][0]*d['metric'][0])
    scale=float(distance/abs(d['lnT'][1]-d['lnT'][0]))
    np.savez_compressed(OUT/'bank.npz',u=u,edges_u=edges,weights=weights,Ci=Ci,Cm=Cm,
        **fields,rate_a=gas.c_light*rho*ka,rate_s=gas.c_light*rho*ks,
        rate_e=rate_e,theta=theta,rate_C=rate_e*theta,arad=gas.a_rad,T=T_K,ne=ne,
        c=gas.c_light,H=scale,local_A=float(d['A'][0]),local_N=float(d['N'][0]),rho=rho)
    ex.write(OUT/'bank.json',dict(classification='Counterexample candidate',passed=True,
        actual_material_state=True,actual_opacity_directly_queried=False,T_K=T_K,T_keV=T,
        opacity_interpolation='Positive log T/log opacity interpolation at fixed physical photon energy on the exact union of both frozen source-cell partitions. Endpoint values extended as in Phase60; no mean fitting.',
        interpolation_fraction=float(fraction),frequency_cells=len(u),computed_means=means.tolist(),free_electrons_m3=ne,
        Thomson_rate_s=rate_e,Compton_rate_s=rate_e*theta,theta=theta,H_T_cm=scale,
        native_differential_scattering_kernel_identified=False,actual_temperature_interpolation_error_certified=False))
    print('BANK',json.loads((OUT/'bank.json').read_text()),flush=True)


def symbolic():
    a,b,C,c,d=sp.symbols('a b C c d',positive=True)
    column=sp.Matrix([-(b-a)/sp.sqrt(C),-a/sp.sqrt(c),b/sp.sqrt(d)])
    energy=sp.Matrix([sp.sqrt(C),sp.sqrt(c),sp.sqrt(d)])
    number=sp.Matrix([0,sp.sqrt(c)/a,sp.sqrt(d)/b])
    assert sp.simplify(column.dot(energy))==0 and sp.simplify(column.dot(number))==0
    mu,v=sp.symbols('mu v',real=True);P2=lambda x:(3*x*x-1)/2
    phase=1+P2(mu)*P2(v)/2
    assert sp.integrate(phase,(v,-1,1))/2==1
    assert sp.integrate(phase*v,(v,-1,1))==0
    assert sp.simplify(sp.integrate(phase*P2(v),(v,-1,1))/2-P2(mu)/10)==0
    return dict(classification='Proven',passed=True,
        structure='Whitened variables z=(sqrt(Cm)*dT, E_l/sqrt(Ci)); z_dot=-(L+i*k*c*V)z. Absorption and Kompaneets are sums of positive outer products; elastic dipole losses are nonnegative; V is real symmetric. Thus the energy norm cannot grow.',
        Kompaneets='For adjacent source cells i,j use b=(-(u_j-u_i)/sqrt(Cm),-u_i/sqrt(Ci),u_j/sqrt(Cj)) and alpha=15*a*T^3/pi^4*rate_C*u_face^4*n0*(1+n0)/(u_j-u_i). The edge generator -alpha*b*b^T preserves material-plus-photon energy and photon number, has the LTE temperature and Bose chemical-potential null modes, and includes induced scattering.',
        origin='Linearizing n_x+n+n^2 about n0=(exp(x)-1)^-1 gives n0(1+n0)*d[delta_n/(n0(1+n0))-x*delta_T/T]/dx. The positive edge bracket discretizes this expression. Exact bin capacities retain the prior radiation heat capacity; the point energy assigned to a cell defines the finite photon-number invariant.',
        spatial='Angular multiplication by mu is symmetric in normalized Legendre modes. Its conservative energy moment is d(Cm*T+sum E0)/dt=-i*k*c*sum E1/sqrt(3). Photon momentum includes a collision force that must be deposited in a moving material solver; this component calculation records the force and does not evolve matter velocity.',
        limitations='Conditional properties of the declared discrete candidate, not physical certification of an ATOMIC gain kernel or of Kompaneets truncation.')


class Operator:
    def __init__(self,bank,order,kH=0.,compton=True,dipole=True):
        keep=bank['u']<=60.;u=bank['u'][keep];Ci=bank['Ci'][keep];Cm=float(bank['Cm'])
        self.u=u;self.Ci=Ci;self.Cm=Cm;self.order=order;self.size=len(u)*order
        self.rate_a=bank['rate_a'][keep];self.rate_s=bank['rate_s'][keep]
        self.kc=kH*float(bank['c']/bank['H']);self.tail=float(bank['Ci'][~keep].sum()/bank['Ci'].sum())
        root=np.sqrt(Ci);cols=np.arange(len(u))*order;self.energy=np.zeros(self.size);self.energy[cols]=root
        self.number=np.zeros(self.size);self.number[cols]=root/u
        loss=np.repeat(self.rate_a,order).reshape(-1,order)
        angular=np.ones(order);angular[0]=0
        if dipole and order>2:angular[2]=.9
        loss+=self.rate_s[:,None]*angular
        L=diags(loss.ravel(),format='csc');q=np.zeros(self.size);q[cols]=-self.rate_a*root/np.sqrt(Cm)
        aa=float(self.rate_a@Ci/Cm)
        if compton:
            face=bank['edges_u'][1:len(u)];n=np.exp(-face)/(-np.expm1(-face))
            alpha=15*float(bank['arad'])*float(bank['T'])**3/np.pi**4*float(bank['rate_C'])*face**4*n*(1+n)/np.diff(u)
            b0=-np.diff(u)/np.sqrt(Cm);left=-u[:-1]/root[:-1];right=u[1:]/root[1:]
            idx=np.arange(len(u)-1);B=coo_matrix((np.r_[left,right],(np.r_[cols[:-1],cols[1:]],np.r_[idx,idx])),shape=(self.size,len(idx))).tocsc()
            L=L+(B@diags(alpha)@B.T).tocsc();q+=np.asarray(B@(alpha*b0)).ravel();aa+=float(alpha@(b0*b0))
            self.K=B;self.alpha=alpha;self.b0=b0
        else:self.K=None
        ids=np.arange(self.size).reshape(-1,order);ell=np.arange(order-1)
        off=(ell+1)/np.sqrt((2*ell+1)*(2*ell+3));values=np.tile(off,len(u))
        V=coo_matrix((np.r_[values,values],(np.r_[ids[:,:-1].ravel(),ids[:,1:].ravel()],np.r_[ids[:,1:].ravel(),ids[:,:-1].ravel()])),shape=(self.size,self.size)).tocsc()
        self.L=L;self.q=q;self.aa=aa;self.V=V;self.P=L+1j*self.kc*V

    def solver(self,dt):
        factor=splu(eye(self.size,format='csc')+dt*self.P);q=dt*self.q
        response=factor.solve(q.astype(complex));denom=1+dt*self.aa-q@response
        def solve(T,E):
            y=factor.solve(E);x=(T-q@y)/denom
            return x,y-response*x
        return solve

    def rhs(self,T,E):return -self.aa*T-self.q@E,-self.q*T-self.P@E

    def evolve(self,duration,steps):
        gamma=1-1/np.sqrt(2);dt=duration/steps;solve=self.solver(gamma*dt)
        T=complex(np.sqrt(self.Cm));E=np.zeros(self.size,complex);initial=self.Cm
        flux=0j;balance=0.;growth=0.;last=initial;residual=0.
        for _ in range(steps):
            U,W=solve(T,E);V,Z=solve(T+(1-gamma)/gamma*(U-T),E+(1-gamma)/gamma*(W-E))
            f1=self.rhs(U,W);f2=self.rhs(V,Z)
            flux+=dt*((1-gamma)*(f1[0]*np.sqrt(self.Cm)+self.energy@f1[1])+gamma*(f2[0]*np.sqrt(self.Cm)+self.energy@f2[1]))
            energy=np.sqrt(self.Cm)*V+self.energy@Z
            balance=max(balance,float(abs(energy-self.Cm-flux)/self.Cm))
            for temp,rad,(dtemp,drad) in [(U,W,f1),(V,Z,f2)]:
                streaming=-1j*self.kc*np.sum(self.Ci**.5*rad.reshape(-1,self.order)[:,1])/np.sqrt(3) if self.order>1 else 0j
                residual=max(residual,float(abs(np.sqrt(self.Cm)*dtemp+self.energy@drad-streaming)/(1+np.sqrt(self.Cm)*abs(dtemp)+np.linalg.norm(drad)*np.sqrt(self.Ci.sum()))))
            norm=float(abs(V)**2+np.vdot(Z,Z).real);growth=max(growth,(norm-last)/initial);last=norm
            T,E=V,Z
        return dict(T=T,E=E,balance=balance,entropy_growth=growth,energy_equation_residual=residual)


def compare(a,b,op):
    return float(np.sqrt(abs(a['T']-b['T'])**2+np.linalg.norm(a['E']-b['E'])**2)/np.sqrt(op.Cm))


def pilot():
    assert not (OUT/'pilot.json').exists();start=time.monotonic();b=np.load(OUT/'bank.npz')
    op=Operator(b,6,1);duration=float(b['H']/b['c']);r=op.evolve(duration,16)
    result=dict(seconds=time.monotonic()-start,size=op.size,frequency_cells=len(op.u),duration_seconds=duration,
        energy_balance=r['balance'],entropy_growth=r['entropy_growth'],energy_equation_residual=r['energy_equation_residual'],
        production_forecast='Measured 6 angular modes at 16 steps. Budget includes 6/12 mode comparisons at 16/32/64 steps, 3 spatial modes and isolated redistribution controls. 180 seconds hard cap; no automatic enlargement.')
    ex.write(OUT/'pilot.json',result);print('PILOT',result,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','fetch','bank','pilot'])
    globals()[p.parse_args().action]()
