"""First-order core heat-current contribution to the existing GR time path.

The expansion parameter multiplies the transport channel, as it multiplies
the frozen reaction source in Phase52. Background gradients drive the poles;
temperature/metric feedback into that channel is second order in this
expansion. This is not a completed nonlinear thermal star or an envelope.
"""
from pathlib import Path
import argparse
import json
import time
import resource
import numpy as np
import sympy as sp
from scipy.sparse.linalg import splu
import def_reactive_cauchy as coupled
import def_conduction_core as core

h=coupled.h
OUT=core.OUT/'coupled'


def symbolic():
    t=sp.symbols('t');y=sp.Function('y')(t);H=sp.Function('H')(t);K,D,F=sp.symbols('K D F')
    assert sp.expand(K*(y-H)-D*sp.diff(y-H,t,2)-D*sp.diff(H,t,2)-F-(K*y-D*sp.diff(y,t,2)-F-K*H))==0
    N,a,w,A,Q=sp.symbols('N a w A Q',positive=True)
    xdot=-N*Q/(a*A**4*w)
    assert sp.simplify(A**4*w*a/N*xdot+Q)==0
    rate,tau,L=sp.symbols('rate tau L',positive=True)
    integral=L*(tau+(sp.exp(-rate*tau)-1)/rate)
    assert sp.simplify(sp.diff(integral,tau)-L*(1-sp.exp(-rate*tau)))==0
    return dict(classification='Proven',passed=True,
        heat_momentum='At first order around zero flux, radial momentum includes partial_tau Q_E plus A^4*w*a/N*xi_tautau. Q_E=G*Rstar^2*q_E/c^5. The pressure equation gains -a*partial_tau Q_E/(N*A^4*p).',
        lift='H_zeta=N/(a*A^4*w*r)*integral Q_E d_tau; H_f=r*Phi*H_zeta; other components zero. y=u-H yields K*u-D*u_tt=F+K*H for the discrete source D*H_tt. Initial y,y_t are zero because all pole fluxes start at zero.',
        pole_integral='Each frozen-gradient pole has L_j(t)=L_j0*(1-exp(-lambda_j*A*N*t)). Its exact time integral is used for material energy and J; there is no microsecond-scale numerical stepping or arbitrary relaxation floor.',
        mass='The same internal face energy is debited from one cell and credited to its neighbor; N*a*J=-G*1e-7*E_face/(c^4*Rstar). Total enclosed core heat energy is zero at the exterior.',
        scope='First order in the strength of a frozen-background core transport channel. The discrete momentum source uses the existing midpoint inertia matrix acting on the kinematic lift; its continuum limit is the local heat momentum term. Nonlinear feedback, envelope/photon currents and spatial convergence are not proved.')


class Heat:
    def __init__(self,geometry,bank):
        self.geometry=geometry;self.d=np.load(coupled.thermal.OUT/'coefficients.npz');d=self.d
        self.edges=np.r_[d['faces_cm'][::-1]/(100*geometry.R),3.]
        self.rnative=d['radius_cm'][::-1]/(100*geometry.R)
        # The same redshift-weighted proper volumes used by the radiation debit.
        self.integral=coupled.Radiation.integral.__get__(self)
        self.volumes=4*np.pi*(100*geometry.R)**3*self.integral(self.edges[:-1],self.edges[1:])
        b=np.load(bank);face=b['faces'];theta=d['A']*d['N']*np.exp(d['lnT'])
        gradient=np.diff(theta)/np.diff(d['radius_cm'])
        N=np.sqrt(d['N'][face-1]*d['N'][face]);A=np.sqrt(d['A'][face-1]*d['A'][face]);a=np.sqrt(d['metric'][face-1]*d['metric'][face])
        fac=-4*np.pi*d['faces_cm'][face]**2*N*A*A/a*gradient[face-1]*1e5
        self.amplitude=fac[:,None]*b['mode_K_SI']
        self.rates=b['poles_proper_s']*(A*N)[:,None]
        self.face_ids=len(d['faces_cm'])-1-face
        self.cache={}

    def faces(self,tau):
        seconds=tau*self.geometry.tc;x=self.rates*seconds
        flux=np.zeros(len(self.edges));energy=flux.copy()
        flux[self.face_ids]=np.sum(self.amplitude*(-np.expm1(-x)),axis=1)
        energy[self.face_ids]=np.sum(self.amplitude*(x+np.expm1(-x))/self.rates,axis=1)
        return flux,energy

    def projection(self,bg):
        key=id(bg['r'])
        if key not in self.cache:
            r=bg['r'];ids=np.clip(np.searchsorted(self.edges,r,side='right')-1,0,len(self.volumes)-1)
            fraction=4*np.pi*(100*self.geometry.R)**3*self.integral(self.edges[ids],np.minimum(r,self.edges[ids+1]))/self.volumes[ids]
            self.cache[key]=(ids,np.clip(fraction,0,1))
        return self.cache[key]

    def source(self,tau,bg):
        _,energy=self.faces(tau);increment=-np.diff(energy);ids,fraction=self.projection(bg)
        prefix=np.r_[0,np.cumsum(increment)];enclosed=prefix[ids]+fraction*increment[ids]
        r=bg['r'];N,a=self.geometry.metric(r);A4=np.exp(-8*bg['phi']**2)
        geo=h.gr.G*.1*self.geometry.R**2/h.gr.C**4
        loss=-increment[ids]/self.volumes[ids]*geo/A4
        interp=lambda v:np.interp(r,self.rnative,v[::-1])
        rho=interp(self.d['raw'][:,0]);U=-loss/(rho*geo)
        rr=interp(self.d['raw'][:,8]/self.d['thermo'][:,5])*U
        J=enclosed*h.gr.G*1e-7/(h.gr.C**4*self.geometry.R*N*a)
        return np.array([rr,loss,np.zeros(len(r)),np.zeros(len(r)),J])

    def lift(self,tau,bg,derivative=False):
        flux,energy=self.faces(tau);quantity=flux*self.geometry.tc if derivative else energy
        r=bg['r'];N,a=self.geometry.metric(r);inside=(r<1)&(bg['p']>0)
        integral=np.interp(r,self.edges,quantity)*h.gr.G*1e-7/(4*np.pi*h.gr.C**4*self.geometry.R*N*N*r*r)
        A4=np.exp(-8*bg['phi']**2);z=np.zeros(len(r))
        z[inside]=N[inside]*integral[inside]/(a[inside]*A4[inside]*(bg['e'][inside]+bg['p'][inside])*r[inside])
        return np.c_[z,np.zeros(len(r)),r*bg['v']*z,np.zeros(len(r))].ravel()


class Radiation(coupled.Radiation):
    def __init__(self,ray,bank,neutrino):
        super().__init__(ray);self.heat=Heat(self.geometry,bank);self.neutrino=neutrino

    def source(self,tau,bg,projection):
        result=self.heat.source(tau,bg)
        if self.neutrino:result+=super().source(tau,bg,projection)
        return result


def solve(steps,bank,neutrino=True,outer=2,label=None):
    start=time.monotonic();ray=dict(np.load(coupled.OUT/'fine-rays.npz'));radiation=Radiation(ray,bank,neutrino)
    bg=coupled.Background(radiation,outer);fn,_=coupled.reactive.symbolic();K,D,forcing=coupled.assemble(bg,fn)
    dt=1/steps;matrix=K-4*D/dt**2;scale=np.asarray(abs(matrix).sum(1)).ravel()
    lu=splu(matrix.multiply((1/scale)[:,None]).tocsc());extended=matrix.astype(np.longdouble)
    u=np.zeros(matrix.shape[0]);velocity=u.copy();acceleration=u.copy();history=[];max_residual=0.
    native_r=bg.native['radius_cm'][::-1]/(100*radiation.geometry.R);weights=bg.native['dm'][::-1];weights/=weights.sum()
    def output(tau):
        y=(u-radiation.heat.lift(tau,bg.nodes)).reshape(-1,4)
        v=(velocity-radiation.heat.lift(tau,bg.nodes,True)).reshape(-1,4)
        r=bg.grid;N,a=radiation.geometry.metric(r);physical_v=a/N*r*v[:,0]*h.gr.C;scalar=y[:,2]-r*y[:,0]*bg.nodes['v']
        cv=np.interp(native_r,r,physical_v);cf=np.interp(native_r,r,scalar)
        return dict(tau=tau,velocity_mass_RMS_m_s=float(np.sqrt(weights@(cv*cv))),scalar_mass_RMS=float(np.sqrt(weights@(cf*cf)))),y,v
    history.append(output(0.)[0])
    for j in range(1,steps+1):
        tau=j*dt;lift=radiation.heat.lift(tau,bg.nodes)
        rhs=forcing(tau)+K@lift-D@(4*u/dt**2+4*velocity/dt+acceleration)
        answer=lu.solve(rhs/scale)
        for _ in range(2):
            defect=rhs.astype(np.longdouble)-extended@answer.astype(np.longdouble)
            answer+=lu.solve(np.asarray(defect/scale,float))
        residual=float(np.max(abs(matrix@answer-rhs)/(abs(matrix)@abs(answer)+abs(rhs)+1e-100)));max_residual=max(max_residual,residual)
        new_acceleration=4*(answer-u-dt*velocity)/dt**2-acceleration
        velocity+=dt*(acceleration+new_acceleration)/2;acceleration=new_acceleration;u=answer
        assert np.all(np.isfinite(u)) and max_residual<1e-9
        history.append(output(tau)[0])
    _,y,v=output(1.);flux,energy=radiation.heat.faces(1.);increment=-np.diff(energy)
    balance=float(abs(np.sum(increment,dtype=np.longdouble))/max(abs(energy).max(),1e-100))
    heat_source=radiation.heat.source(1.,bg.nodes);rr,loss,E,P,J=radiation.source(1.,bg.nodes,bg.node_projection)
    rho=np.interp(bg.grid,radiation.heat.rnative,bg.native['raw'][::-1,0]);geo=h.gr.G*.1*radiation.geometry.R**2/h.gr.C**4
    cvT=np.interp(bg.grid,radiation.heat.rnative,bg.native['thermo'][::-1,3])
    logT=bg.nodes['adiabatic_T_rho']*y[:,1]/bg.nodes['gamma']-heat_source[1]/(rho*geo*cvT)
    if neutrino:logT+=bg.nodes['theta_rate']*radiation.geometry.tc
    # Derive the mass and radial metric from the SAME heat/neutrino debit.
    r=bg.grid;xi=r*y[:,0];b=1-2*bg.nodes['m']/r;A4=np.exp(-8*bg.nodes['phi']**2)
    Dm=r*r*b*bg.nodes['v']*y[:,2]-(4*np.pi*r*r*A4*bg.nodes['p']+r*r*b*bg.nodes['v']**2/2)*xi+J
    dl=(Dm/r-bg.nodes['m']*xi/r**2)/b
    row=dict(classification='Counterexample candidate',steps=steps,outer=outer,neutrino=neutrino,seconds=time.monotonic()-start,
        linear_residual=max_residual,heat_telescoping=balance,history=history,max_heat_luminosity_erg_s=float(abs(flux).max()),
        maximum_heat_redistributed_erg=float(abs(energy).max()),max_abs_Delta_log_T=float(abs(logT[r<1]).max()),
        exterior_heat_mass_source=float(abs(heat_source[4,r>=1]).max()))
    if label:
        np.savez_compressed(OUT/(label+'.npz'),grid=bg.grid,response=y,velocity=v,Delta_m_over_R=Dm,Delta_lambda=dl,Delta_log_T=logT,
            heat_face_r=radiation.heat.edges,heat_luminosity_erg_s=flux,heat_cumulative_energy_erg=energy)
        core.ex.write(OUT/(label+'.json'),row)
    print('PATH',label,steps,'SECONDS',row['seconds'],'END',history[-1],flush=True)
    return row


def prepare():
    assert not OUT.exists();OUT.mkdir();assert json.loads((core.OUT/'result.json').read_text())['passed']
    paths=[Path(__file__),Path(coupled.__file__),Path(core.__file__),core.OUT/'result.json',core.OUT/'bank.npz',core.OUT/'coarse-bank.npz',
           coupled.OUT/'fine-rays.npz',coupled.OUT/'fine-64.npz',coupled.OUT/'fine-64.json',coupled.surface.OUT/'background.npz',coupled.thermal.OUT/'coefficients.npz',coupled.thermal.OUT/'sources.npz']
    core.ex.write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='cca3717e',bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths},
        claim='Add the microscopic core heat-current contribution to the actual free-surface fluid/scalar/metric time path, including exact pole startup, conservative material energy, GR mass and heat momentum.',
        expansion='First order in the strength of the initially zero heat channel, with its physical microscopic rate bank and initial background gradient. Perturbed temperature/metric/acceleration in the heat drive multiply the channel parameter and are next order; the finite-amplitude feedback error is not bounded.',
        domain='Internal faces of the original degenerate core only; omitted outer currents, photons, physical surface and nonlinear feedback remain unclosed. This is a component of the full linear source response, not an assertion of a physically insulated core.',
        horizon_seconds=.23080495568542375,steps=[16,32,64],additional_paths=['heat-only-64','coarse-bank-64','outer-64'],
        gates=dict(linear_residual=1e-9,heat_balance=2e-13,time_relative=.02,time_order=1.5,coefficient_relative=.02,outer_relative=.002,superposition=2e-8),
        budget=dict(pilot_steps=8,production_steps=304,maximum_paths=6,hard_seconds=180,CPU_workers=1,BLAS_threads=1,native_calls=0,automatic_expansion=False)))
    core.ex.write(OUT/'symbolic.json',symbolic());row=solve(8,core.OUT/'bank.npz',label='pilot')
    forecast=row['seconds']*304/8*1.2
    core.ex.write(OUT/'pilot-budget.json',dict(seconds=row['seconds'],forecast_seconds=forecast,assumption='Conservative proportional scaling of measured eight-step construction+solve to 304 steps, despite factorization reuse per path. Stop at 180 seconds.'))
    print('FORECAST',forecast,flush=True)


def run():
    assert not (OUT/'result.json').exists();plan=json.loads((OUT/'plan.json').read_text())
    for p,s in plan['bindings'].items():assert h.digest(h.ROOT/p)==s,p
    assert json.loads((OUT/'pilot-budget.json').read_text())['forecast_seconds']<170
    start=time.monotonic();cases={}
    for steps in [16,32,64]:cases[str(steps)]=solve(steps,core.OUT/'bank.npz',label=f'combined-{steps}')
    cases['heat']=solve(64,core.OUT/'bank.npz',False,label='heat-only-64')
    cases['coarse']=solve(64,core.OUT/'coarse-bank.npz',label='coarse-bank-64')
    cases['outer']=solve(64,core.OUT/'bank.npz',outer=3,label='outer-64')
    comparisons={}
    for field in ['velocity_mass_RMS_m_s','scalar_mass_RMS']:
        series=lambda name:np.array([x[field] for x in cases[name]['history']])
        a,b,c=series('16'),series('32'),series('64');norm=max(abs(c).max(),1e-100)
        d1=np.max(abs(a-b[::2]))/norm;d2=np.max(abs(b-c[::2]))/norm
        comparisons[field]=dict(time_previous=float(d1),time_last=float(d2),order=float(np.log2(d1/d2)),
            coefficients=float(np.max(abs(c-series('coarse')))/norm),outer=float(np.max(abs(c-series('outer')))/norm))
    combined=np.load(OUT/'combined-64.npz');only=np.load(OUT/'heat-only-64.npz');old=np.load(coupled.OUT/'fine-64.npz')
    superposition={k:float(np.max(abs(combined[k]-only[k]-old[k]))/max(abs(combined[k]).max(),1e-100)) for k in ['response','velocity']}
    passed=all(v['time_last']<.02 and v['order']>1.5 and v['coefficients']<.02 and v['outer']<.002 for v in comparisons.values())
    passed=passed and max(superposition.values())<2e-8 and max(x['heat_telescoping'] for x in cases.values())<2e-13
    result=dict(classification='Counterexample candidate',coupled_core_heat_time_path_completed=True,passed=passed,comparisons=comparisons,superposition=superposition,
        seconds=time.monotonic()-start,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        heat_only_endpoint=cases['heat']['history'][-1],combined_endpoint=cases['64']['history'][-1],
        full_temperature_feedback=False,physical_core_interface_closed=False,spatial_error_certified=False,photon_transport=False,full_dynamic_charge_solved=False)
    core.ex.write(OUT/'result.json',result);print('RESULT',json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
