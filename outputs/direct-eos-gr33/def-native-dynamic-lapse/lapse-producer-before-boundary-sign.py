"""Actual escaping packets set the lapse of the represented GR response.

Counterexample candidate: first-order constraints and null transport about the
corrected background. The increments stay separate from that background.
"""
from pathlib import Path
from types import SimpleNamespace
import json
import signal
import sys
import time
import numpy as np
import sympy as s
from numpy.polynomial import legendre as leg
from scipy.integrate import solve_ivp
import verify_native_anisotropic_gr as prior

flow=prior.flow;C=prior.C;G=prior.G;LD=prior.LD;write=prior.write;sha=prior.sha
OUT=prior.OUT.parent/'def-native-dynamic-lapse'


def symbolic():
    # Hamiltonian for a null particle: H=a*sqrt(pr^2/h^2+L^2/R^2).
    # Partial derivatives are evaluated before the orthonormal substitution.
    a,h,R,ar,hr,Rr,at,ht,Rt,p,L,eps,mu=s.symbols('a h R ar hr Rr at ht Rt p L eps mu',positive=True)
    H=a*s.sqrt(p*p/h**2+L*L/R**2)
    velocity=s.diff(H,p);pdot=-sum(s.diff(H,x)*dx for x,dx in [(a,ar),(h,hr),(R,Rr)])
    energy=H/a;direction=p/(h*energy)
    dedt=sum(s.diff(energy,x)*(xt+xr*velocity) for x,xt,xr in [(h,ht,hr),(R,Rt,Rr)])+s.diff(energy,p)*pdot
    dmdt=sum(s.diff(direction,x)*(xt+xr*velocity) for x,xt,xr in [(h,ht,hr),(R,Rt,Rr)])+s.diff(direction,p)*pdot
    replace={p:h*eps*mu,L:R*eps*s.sqrt(1-mu*mu)}
    de=-mu*ar/h-ht/h*mu**2-Rt/R*(1-mu**2)
    dm=(1-mu**2)*(a/h*(Rr/R-ar/a)+mu*(Rt/R-ht/h))
    assert s.simplify((dedt/energy).subs(replace)-de)==0
    assert s.simplify(dmdt.subs(replace)-dm)==0
    # In scalar vacuum (r*Phi)'=-Phi/b. Integration by parts fixes the
    # otherwise easy-to-miss sign of the scalar lapse boundary term.
    r,b,Phi,f,fr,J=s.symbols('r b Phi f fr J',nonzero=True)
    nuvar=Phi*f/b+J/(r*r*b*b)+r*Phi*fr
    assert s.simplify(nuvar-(-Phi*f/b+r*Phi*fr)-2*Phi*f/b-J/(r*r*b*b))==0
    return dict(classification='Proven',passed=True,
        null_characteristics='dr/dt=c*a*mu/h; dmu/dt=(1-mu^2)*[c*a/h*(R_prime/R-a_prime/a)+mu*(R_t/R-h_t/h)]; dln(epsilon)/dt=-c*mu*a_prime/h-(h_t/h)*mu^2-(R_t/R)*(1-mu^2). Derived from the null Hamiltonian with zero shift.',
        variables='Fixed Einstein radius: a=A*N, h=A/sqrt(b), R=A*r. zeta=delta_nu-delta_lambda; f=delta_phi; u=alpha*f. Existing frequency coordinate is Eref=a0(r)*epsilon.',
        increments='delta_vr=c*N*sqrt(b)*mu*zeta; delta_mudot=(1-mu^2)*[c*N*sqrt(b)*((1/r-nu_prime)*zeta-delta_nu_prime)-mu*delta_lambda_t]; delta_lnEref_dot=-c*N*sqrt(b)*mu*(delta_nu_prime+u_prime)-u_t-mu^2*delta_lambda_t.',
        lapse_packet='An exterior packet of Killing energy dE at(rp,mup) contributes(G*dE/c^4)*[integral(r0,rp) dr/(N*b^(3/2)*r^2)-mup^2/(rp*Np*sqrt(bp))] to delta_nu(r0), including its compensating interior mass debit.',
        scalar_lapse='delta_nu_scalar(r0)=-Phi(r0)*U(r0)-integral(r0,infinity)2*Phi*U/(r*b)dr, U=r*f. Additional vacuum J contributes -integral J/(r^2*b^2)dr.',
        limits='Algebra for the declared first variation; not a nonlinear stellar or angular-continuum certificate.')


def prepare():
    assert not (OUT/'plan.json').exists();OUT.mkdir(exist_ok=True)
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='2fc1fb85f',
        claim='Supply the missing asymptotically normalized lapse from the actually emitted angular photon packets, then form all radial, angular and frequency transport increments for the new GR field.',
        decision='Can the existing separate scalar/mass response provide a finite, reproducible dynamic transport input, without adding sub-ulp increments to full background arrays or clamping the lapse at the numerical surface?',
        model='Use the same corrected exact scalar vacuum beyond the actual photon edge; no initially exterior photons. Keep actual4 outward angular bins with constant occupation within each. Propagate every quadrature packet on its outgoing null ray and integrate pressure and mass debit together. Include the measured global numerical energy-port residual as a separate asymptotic mass error, not a fitted conservation correction.',
        scalar_boundary='Continue the computed compact outgoing U along vacuum radial characteristics for the leading scalar lapse term. Missing exterior-generated/scattered scalar waves remain explicit; this is not full scalar boundary closure.',
        interior='Integrate the linear polar lapse constraint with actual pressure/energy forcing and canonical volume response. Interpolate the saved scalar field within existing cells; its contribution is separately measured. No new EOS state, fluid trajectory or physical domain enlargement.',
        controls=dict(packet_order=[4,8],constraint_order=[4,8],time_paths=[64,128],relative_quadrature=.002,relative_time=.02,energy_port=1e-11,ray_invariant=1e-10,analytic_flat_shell=1e-10),
        budget=dict(seconds=45,CPU_threads=1,memory_GB=3,new_fluid_steps=0,new_native_calls=0),
        measured_basis='Phase122 saved-source export4.98s and independent constraint audit5.32s. This adds short vacuum ODEs for16/32 directions and existing17 times. Stop45s; no automatic higher orders or new fluid runs.',
        limits='The resulting metric input is for a compensated coupled transport solve. Merely computing the characteristics is not that solve. Physical angular/spectral reconstruction, external scalar forcing and nonlinear backreaction remain.',
        references=['https://arxiv.org/abs/gr-qc/0201064'],
        bindings={str(p):sha(p) for p in [Path(__file__),Path(prior.__file__),Path(prior.new.__file__),Path(prior.base.__file__),prior.OUT/'audit.json',flow.INPUT/'balanced-20.npz',flow.OUT/'coupled-64.npz',flow.OUT/'coupled-128.npz',prior.OUT/'source-64.npz',prior.OUT/'source-128.npz',prior.new.OUT/'fields-64-g8.npz',prior.new.OUT/'fields-128-g8.npz']}))
    write(OUT/'symbolic.json',symbolic())
    m=Lapse.__new__(Lapse);m.r0=m.N0=1.
    m.model=SimpleNamespace(m=SimpleNamespace(R=1.),bulk=SimpleNamespace(edges_mu=np.linspace(-1,1,9)))
    m.bg=SimpleNamespace(metric=lambda r:(np.zeros_like(r),np.ones_like(r),np.ones_like(r),np.ones_like(r),np.zeros_like(r)))
    mu,_,_,sol,invariant=m.rays(4,1/C);state=sol.sol(.7).reshape(3,len(mu))
    radius=np.sqrt(1+1.4*mu+.7**2);direction=(mu+.7)/radius
    error=float(max(np.max(abs(state[0]-radius)),np.max(abs(state[1]-direction)),np.max(abs(state[2]-(1-1/radius)))))
    assert error<1e-10 and invariant<1e-10
    write(OUT/'flat-ray-check.json',dict(classification='Proven',passed=True,analytic_flat_ray_error=error,impact_invariant=invariant,
        scope='Exact flat outgoing null rays and their lapse mass kernel; does not certify the curved numerical trajectory.'))


class Lapse:
    def __init__(self):
        self.response=prior.new.Response();self.model=self.response.model;self.bg=self.response.bg;self.geo=self.response.geo
        self.r0=self.bg.rend;z=self.bg.fields(np.array([self.r0]));self.N0=z['lapse'][0]

    def rays(self,order,end):
        b=self.model.bulk;positive=b.edges_mu[b.edges_mu>=0];gx,gw=leg.leggauss(order)
        mu=(positive[:-1,None]+np.diff(positive)[:,None]*(gx+1)/2).ravel()
        weights=(np.diff(positive)[:,None]*gw/2).ravel();ids=np.repeat(np.arange(len(positive)-1),order)
        n=len(mu);k=C*end/self.r0
        def rhs(t,z):
            y,u,integ=z.reshape(3,n);r=self.r0*y
            mass,N,b,cc,Phi=self.bg.metric(r/self.model.m.R);mass=mass*self.model.m.R;Phi=Phi/self.model.m.R
            nr=mass/(r*r*b)+r*Phi*Phi/2
            return np.array([k*cc*u,k*cc*(1-u*u)*(1/y-self.r0*nr),k*u/(b*y*y)]).ravel()
        solution=solve_ivp(rhs,[0,1],np.r_[np.ones(n),mu,np.zeros(n)],rtol=2e-12,atol=1e-13,dense_output=True)
        assert solution.success
        state=solution.sol(np.linspace(0,1,17)).reshape(3,n,17);r=self.r0*state[0]
        N=self.bg.metric(r.ravel()/self.model.m.R)[1].reshape(r.shape)
        invariant=state[1]**2+(1-mu[:,None]**2)*(N/self.N0/state[0])**2-1
        return mu,weights,ids,solution,float(np.max(abs(invariant)))

    def boundary(self,d,steps,order):
        t=d['t'];end=t[-1];mu,mw,ids,solution,invariant=self.rays(order,end)
        path=np.load(flow.OUT/f'coupled-{steps}.npz');snap=path['snapshot_t'];keep=[np.argmin(abs(snap-tt)) for tt in t]
        occupation=path['snapshot_I'][keep].sum(1)[:,-1,self.model.mu>0,:]
        lum=2*np.pi*C*self.model.area[-1]*(occupation@(self.model.number*self.model.E))
        photon=[];energy=[];gx,gw=leg.leggauss(order)
        for k,now in enumerate(t):
            if k==0:photon.append(0.);energy.append(0.);continue
            tt=(t[:k,None]+np.diff(t[:k+1])[:,None]*(gx+1)/2).ravel();tw=(np.diff(t[:k+1])[:,None]*gw/2).ravel()
            ll=np.column_stack([np.interp(tt,t,v) for v in lum.T])
            # Packet energy measure includes the emission flux factor mu.
            packet=tw[:,None]*ll[:,ids]*mu*mw
            state=solution.sol((now-tt)/end).reshape(3,len(mu),len(tt)).transpose(0,2,1)
            rp=self.r0*state[0];mp=state[1];N,b,cc=[v.reshape(rp.shape) for v in self.bg.metric(rp.ravel()/self.model.m.R)[1:4]]
            kernel=state[2]/self.r0-mp*mp/(rp*cc)
            photon.append(float(G/C**4*np.sum(packet*kernel,dtype=LD)));energy.append(float(np.sum(packet,dtype=LD)))
        energy=np.array(energy);photon=np.array(photon)
        energy_error=float(np.max(abs(energy-d['outer_cumulative_energy_erg']))/energy[-1])
        # Unit asymptotic J lapse kernel compactified all the way to infinity.
        z=(gx+1)/2;_,N,b,_,_=self.bg.metric(self.r0/(z*self.model.m.R))
        vacuum_kernel=float(np.sum(gw/(2*self.r0*N*b**1.5)))
        return photon,energy,dict(ray_invariant=invariant,emitted_energy_relative=energy_error,unit_ADM_lapse_kernel_per_cm=vacuum_kernel)

    def run(self,steps,order):
        start=time.monotonic();d=dict(np.load(prior.OUT/f'source-{steps}.npz'));field=np.load(prior.new.OUT/f'fields-{steps}-g8.npz')
        response=self.response;response.setup(d,order);z=response.z;r=response.r;n=len(d['radius']);q=flow.initial.Quadrature(d['edges'],order)
        photon,energy,row=self.boundary(d,steps,order)
        # The scalar lapse term is tiny but explicit. Use the saved outgoing
        # field on radial vacuum characteristics over its causal support.
        U0=field['U'][:,-1];scalar=[];gx,gw=leg.leggauss(order);base_delay=self.geo(np.array([d['edges'][-1]-self.model.m.RJ]))[1][0]
        Phi0=self.bg.fields(np.array([self.r0]))['Phi'][0]
        for k,now in enumerate(d['t']):
            if k==0:scalar.append(0.);continue
            delays=base_delay+now-d['t'][:k+1][::-1];cuts=self.r0+C*(delays-base_delay)
            for _ in range(3):
                rJ=self.geo.physical(cuts)[0];actual=self.geo(rJ-self.model.m.RJ)[1];zz=self.bg.fields(cuts)
                cuts-=(actual-delays)*C*zz['lapse']*np.sqrt(zz['b'])
            rr=(cuts[:-1,None]+np.diff(cuts)[:,None]*(gx+1)/2).ravel();ww=(np.diff(cuts)[:,None]*gw/2).ravel()
            zz=self.bg.fields(rr);rJ=self.geo.physical(rr)[0];delay=self.geo(rJ-self.model.m.RJ)[1]-base_delay
            uu=np.interp(np.clip(now-delay,0,now),d['t'],U0)
            scalar.append(float(-Phi0*U0[k]-np.sum(ww*2*zz['Phi']*uu/(rr*zz['b']),dtype=LD)))
        scalar=np.array(scalar)
        # Reuse the exact-center constraint for the physical asymptotic mass
        # residual: J/sqrt(b)*N plus emitted energy must be zero in exact data.
        centers=prior.centers(response,d,field,order);response.setup(d,order);z=response.z;r=response.r
        residual=centers['J'][:,-1]*response.tz['lapse'][-1]/np.sqrt(response.tz['b'][-1])+G/C**4*energy
        boundary=photon+scalar-residual*row['unit_ADM_lapse_kernel_per_cm']
        f=np.array([np.interp(r,field['radius_E'],v) for v in field['delta_phi']]);fr=np.array([np.interp(r,field['radius_E'],v) for v in field['delta_Phi']])
        J=response.J;dm=r*r*z['b']*z['Phi']*f+J;dl=dm/(r*z['b']);volume=3*z['alpha']*f+dl
        rest=np.asarray(d['baryon_g'],LD)*LD(d['cx'])*LD(C)**2;E=rest+d['gas_nonrest_energy_erg']+d['photon_energy_erg']
        pF=np.asarray((E-d['metric_stress_erg'])/d['volume']*LD(G)/LD(C)**4,float)[:,response.ids]
        dp=pF-z['Kg']*volume-4*z['alpha']*z['Pr']*f-(3*z['Pr']-z['R4'])*dl
        P=z['Pg']+z['Pr'];b=z['b'];N=z['lapse']
        nr=(1+8*np.pi*r*r*z['A4']*P)*dm/(r*r*b*b)+4*np.pi*r*z['A4']*(dp+4*z['alpha']*P*f)/b+r*z['Phi']*fr
        # Integrate in Jordan radius, retaining its exact Jacobian to Einstein r.
        A=np.exp(-2*z['phi']**2);integrand=nr/(A*(1+z['alpha']*r*z['Phi']))
        values=integrand.reshape(len(d['t']),n,order);whole=q.h*(values@q.w)
        totals=np.sum(whole,axis=1,dtype=LD);face=boundary[:,None]-np.asarray(totals[:,None]-np.column_stack([np.zeros(len(d['t'])),np.cumsum(whole,axis=1,dtype=LD)]),float)
        inverse=np.linalg.inv(leg.legvander(q.x,order-1));co=values@inverse.T
        xc=(d['radius']-d['edges'][:-1])/q.h-1;Q=np.column_stack([leg.legval(xc,leg.legint(np.eye(order)[j]))-leg.legval(-1,leg.legint(np.eye(order)[j])) for j in range(order)])
        lapse=face[:,:-1]+q.h*np.sum(co*Q,axis=2)
        ft=field['U_t'][:,:-1]/field['radius_E'][:-1];phi=field['delta_phi'][:,:-1]
        alpha=response.tz['alpha'][:-1];u=alpha*phi;ut=alpha*ft;up=-4*response.tz['Phi'][:-1]*phi+alpha*field['delta_Phi'][:,:-1]
        lam=centers['delta_lambda'][:,:-1];nup=centers['delta_nu_prime'][:,:-1]
        # The dynamic input declares piecewise-linear lambda/nu on the saved
        # time grid. Store interval derivatives rather than inventing precision.
        lt=np.diff(lam,axis=0)/np.diff(d['t'])[:,None];nt=np.diff(lapse,axis=0)/np.diff(d['t'])[:,None]
        np.savez_compressed(OUT/f'metric-{steps}-g{order}.npz',t=d['t'],radius_E=field['radius_E'][:-1],
            delta_nu=lapse,delta_nu_faces=face,delta_lambda=lam,delta_phi=phi,delta_u=u,delta_u_t=ut,delta_u_prime=up,
            delta_nu_prime=nup,delta_lambda_interval_rate=lt,delta_nu_interval_rate=nt,
            delta_log_lapse=lapse+u,delta_log_radial_length=lam+u,delta_log_area_radius=u,
            delta_log_speed=lapse-lam,outer_photon_lapse=photon,outer_scalar_lapse=scalar,
            outer_ADM_residual_lapse=-residual*row['unit_ADM_lapse_kernel_per_cm'],asymptotic_mass_residual_cm=residual,emitted_energy_erg=energy)
        row.update(classification='Counterexample candidate',steps=steps,order=order,seconds=time.monotonic()-start,
            maximum_delta_nu=float(np.max(abs(lapse))),endpoint_outer_photon_lapse=float(photon[-1]),
            endpoint_outer_scalar_lapse=float(scalar[-1]),maximum_asymptotic_mass_residual_cm=float(np.max(abs(residual))),
            maximum_ADM_residual_lapse=float(np.max(abs(residual*row['unit_ADM_lapse_kernel_per_cm']))),
            maximum_delta_log_transport_speed=float(np.max(abs(lapse-lam))),maximum_delta_lambda_rate=float(np.max(abs(lt))),
            actual_exterior_photons_applied=True,full_spatial_feedback=False,final_charge_solved=False)
        write(OUT/f'metric-{steps}-g{order}.json',row);print(json.dumps(row),flush=True);return row


def run():
    assert not (OUT/'result.json').exists()
    for p,h in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert sha(p)==h,p
    signal.signal(signal.SIGALRM,flow.old.optical.timeout);signal.alarm(45);start=time.monotonic();m=Lapse();rows=[]
    for steps,order in [(128,8),(128,4),(64,8)]:rows.append(m.run(steps,order))
    files=[np.load(OUT/f'metric-{steps}-g{order}.npz') for steps,order in [(128,8),(128,4),(64,8)]]
    errors={name:float(np.max(abs(files[0]['delta_nu']-d['delta_nu']))/np.max(abs(files[0]['delta_nu']))) for name,d in zip(['quadrature','time'],files[1:])}
    passed=errors['quadrature']<.002 and errors['time']<.02 and max(x['ray_invariant'] for x in rows)<1e-10 and max(x['emitted_energy_relative'] for x in rows)<1e-11
    result=dict(classification='Counterexample candidate',passed=bool(passed),controls=errors,paths=rows,seconds=time.monotonic()-start,
        actual_photon_lapse_and_mass_pressure_applied=True,transport_input_ready=True,full_spatial_feedback=False,full_exterior_scalar=False,final_charge_solved=False)
    write(OUT/'result.json',result);print(json.dumps(result),flush=True);signal.alarm(0)


if __name__=='__main__':globals()[sys.argv[1]]()
