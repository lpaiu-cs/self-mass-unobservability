"""Simultaneous lapse, grey heat-flux geometry and moving-photon feedback.

Counterexample candidate: native frozen-coefficient tangent model with the
prescribed finite-pressure support. The support's gravitational stress tensor
and a nonlinear atmosphere are not supplied by this transport correction.
"""
from pathlib import Path
import argparse
import inspect
import json
import resource
import signal
import textwrap
import time
import numpy as np
from scipy.sparse import diags
from scipy.interpolate import CubicHermiteSpline
import def_native_balanced_flux as prior
import def_native_moving_rays as moving

old=prior.old
OUT=prior.OUT.parent/'def-native-metric-transport'
write=old.write


stage_source=textwrap.dedent(inspect.getsource(prior.Model.stage))
for a,b in [
    ('def step(state,t):','def solve_once(state,t,correction,photons):'),
    ('self.photon_force(float(t))','photons'),
    ('t*self.LE0+LE@pred-dr','t*self.LE0+LE@pred-dr+correction'),
    ('targetflux=Lq@q+t*self.LE0+LE@e','targetflux=Lq@q+t*self.LE0+LE@e+correction'),
    ('return step','return self.close_stage(solve_once)')]:
    assert stage_source.count(a)==1,a
    stage_source=stage_source.replace(a,b)
namespace=dict(prior.Model.stage.__globals__)
exec(compile(stage_source,__file__,'exec'),namespace)

evolve_source=textwrap.dedent(inspect.getsource(prior.Model.evolve))
anchor='self.history_t=[0.];self.history_e=[0.];self.history_d=[0.]'
assert evolve_source.count(anchor)==1
evolve_source=evolve_source.replace(anchor,anchor+'\n    self.history_z=[0.];self.history_zd=[0.];self.max_closure_error=0.;self.max_closure_iterations=0')
anchor='self.history_d.append(float(state[3][-1]))'
assert evolve_source.count(anchor)==1
evolve_source=evolve_source.replace(anchor,anchor+'\n            self.history_z.append(float((self.surfaceV[0]@state[0])[0]));self.history_zd.append(float((self.surfaceV[0]@state[1])[0]))')
anchor='rows.append(row);temp.append(T);vel.append(velocity);scalar.append(s);surface.append(surf)'
assert evolve_source.count(anchor)==1
evolve_source=evolve_source.replace(anchor,'self.closure(state,t)\n        row.update(surface_lapse=float(self.last_lapse_surface),maximum_lapse=float(max(abs(self.last_lapse_native))))\n        '+anchor)
anchor='moving_surface_solved=False,final_charge_solved=False,full_goal_complete=False)'
assert evolve_source.count(anchor)==1
evolve_source=evolve_source.replace(anchor,'max_closure_relative=self.max_closure_error,max_closure_iterations=self.max_closure_iterations,\n        '+anchor)
evolve_namespace=dict(prior.Model.evolve.__globals__);evolve_namespace['OUT']=OUT
exec(compile(evolve_source,__file__,'exec'),evolve_namespace)


class Model(prior.Model):
    stage=namespace['stage']
    evolve=evolve_namespace['evolve']

    def __init__(self):
        start=time.monotonic();super().__init__()
        mask=self.points['r']<1;p={k:v[mask] for k,v in self.points.items()}
        r=p['r'];b=1-2*p['m']/r;A=np.exp(-8*p['phi']**2);alpha=-4*p['phi'];w=p['e']+p['p']
        gamma=p['gamma'];v=p['v'];ad=np.interp(r,self.native,self.thermo[:,4]);V,D=self.evaluation(r)
        nr=p['m']/(r*r*b)+4*np.pi*r*A*p['p']/b+r*v*v/2
        g=nr+alpha*v
        dlz=r*r*v*v/2-4*np.pi*r*r*A*p['p']/b-p['m']/(r*b)
        mass=-diags(4*np.pi*r**3*A*w)@V[0]+diags(r*r*b*v)@V[1]
        pressure=-diags(gamma*p['p']*r)@D[0]+diags(-gamma*p['p']*(3+dlz+3*alpha*r*v)+r*w*g)@V[0]
        pressure-=diags(gamma*p['p']*(r*v+3*alpha))@V[1]
        cm=(1+8*np.pi*r*r*A*p['p'])/(r*r*b*b);cp=4*np.pi*r*A/b
        self.nuprime_q=(diags(cm)@mass+diags(r*v)@D[1]+diags(cp)@(pressure+diags(4*alpha*p['p'])@V[1])).tocsc()
        _,loss,J=self.sources(p);cj=cm-cp*gamma*p['p']/(r*b)
        self.nuprime_e=(diags(cj)@J-diags(cp*ad)@loss).tocsc()
        _,loss0,J0=self.source_values(p,self.f0);self.nuprime0=cj*J0-cp*ad*loss0
        self.lapse_weights=self.weights[mask].astype(np.longdouble)
        self.lapse_native_cut=np.searchsorted(r,self.native)
        self.lapse_face_cut=np.searchsorted(r,self.edges[1:-1])
        self.rays_motion=moving.Rays(self.bg,self.points['r'][self.outside])
        ext={k:v[self.outside] for k,v in self.points.items()};re=ext['r'];be=1-2*ext['m']/re
        Ve,_=self.evaluation(re);ps=self.bg.sample(np.array([1.]));bs=1-2*ps['m'][0]
        self.lapse_q_surface=ps['v'][0]*self.surfaceV[1].toarray().ravel()-2*(Ve[1].T@(self.weights[self.outside]*ext['v']/be))
        self.lapse_q_surface=self.lapse_q_surface.astype(np.longdouble)
        pp=self.bg.sample(self.native);self.alpha_native=-4*pp['phi']
        face=self.bg.sample(self.edges[1:-1]);rf=face['r'];bf=1-2*face['m']/rf;Af=np.exp(-8*face['phi']**2)
        self.faceV,_=self.evaluation(rf);self.face_alpha=-4*face['phi'];self.face_rb=rf*bf
        self.face_mass_q=-diags(4*np.pi*rf**3*Af*(face['e']+face['p']))@self.faceV[0]+diags(rf*rf*bf*face['v'])@self.faceV[1]
        self.face_J=-old.G/(old.C**4*self.bg.R)*np.sqrt(bf)/face['N']
        nsprime=ps['m'][0]/bs+4*np.pi*np.exp(-8*ps['phi'][0]**2)*ps['p'][0]/bs+ps['v'][0]**2/2
        self.alpha_surface=-4*ps['phi'][0]
        self.surface_geometry=2+4*self.alpha_surface*ps['v'][0]+2*nsprime
        self.setup_seconds=time.monotonic()-start

    def radiation(self,state,t):
        q,v,e,d=state;times=self.history_t
        energy=self.history_e;flux=self.history_d;z=self.history_z;zd=self.history_zd
        if t>times[-1]:
            times=times+[float(t)];energy=energy+[float(e[-1])];flux=flux+[float(d[-1])]
            z=z+[float((self.surfaceV[0]@q)[0])];zd=zd+[float((self.surfaceV[0]@v)[0])]
        rays=self.rays_motion;u=t-rays.delay;active=u>0;clipped=np.maximum(u,0)
        if len(times)>1:
            heat=CubicHermiteSpline(times,energy,flux);move=CubicHermiteSpline(times,z,zd)
            H=(self.f0[-1]*clipped+heat(clipped))*active
            L=(self.f0[-1]+heat(clipped,1))*active
            hm,epm,pm=rays.moments(float(t),move)
        else:
            H=np.zeros_like(u);L=np.zeros_like(u);hm=epm=pm=np.zeros(len(rays.r))
        factor=old.G/(old.C**4*self.bg.R);r=rays.r;mu=rays.mu
        J=-factor*(H@rays.weights+self.f0[-1]*hm)*np.sqrt(rays.b)/rays.n
        P=factor*((L*mu)@rays.weights+self.f0[-1]*pm)/(4*np.pi*r*r*rays.n**2)
        EP=factor*((L*(1/mu-mu))@rays.weights+self.f0[-1]*epm)/(4*np.pi*r*r*rays.n**2)
        force=self.photon_test@(2*rays.c*rays.v*J/rays.b**2-4*np.pi*r**3*rays.c*rays.v*EP/rays.b)
        lapse=-np.sum(self.weights[self.outside]*(J/(r*r*rays.b**2)+4*np.pi*r*P/rays.b),dtype=np.longdouble)
        return force,lapse

    def closure(self,state,t):
        q,v,e,d=state;photons,lapse_photons=self.radiation(state,float(t))
        ns=self.lapse_q_surface@q+lapse_photons
        derivative=self.nuprime_q@q+t*self.nuprime0+self.nuprime_e@e
        integral=np.r_[np.longdouble(0),np.cumsum(self.lapse_weights*derivative,dtype=np.longdouble)]
        nu=ns-(integral[-1]-integral[self.lapse_native_cut])
        nf=ns-(integral[-1]-integral[self.lapse_face_cut])
        psi=self.nativeV[1]@q;pf=self.faceV[1]@q
        dm=self.face_mass_q@q+self.face_J*(t*self.f0[:-1]+e[:-1])
        internal=self.thermal_operator@(nu+self.alpha_native*psi)
        internal+=self.f0[:-1]*(nf+2*self.face_alpha*pf-dm/self.face_rb)
        surface=self.f0[-1]*(self.surface_geometry*(self.surfaceV[0]@q)[0]+4*self.alpha_surface*(self.surfaceV[1]@q)[0]+2*ns)
        self.last_lapse_surface=ns;self.last_lapse_native=nu
        return np.r_[internal,surface],photons

    def close_stage(self,solve_once):
        def step(state,t):
            guess=state;correction,photons=self.closure(guess,t)
            for iteration in range(1,7):
                answer=solve_once(state,t,correction,photons)
                new_c,new_p=self.closure(answer,t)
                # Assess the added closure itself, not its tiny ratio to the
                # old baseline flux. No residual-floor success is substituted.
                err=max(float(max(abs(new_c-correction))/max(max(abs(new_c)),1e-100)),
                        float(max(abs(new_p-photons))/max(max(abs(new_p)),1e-100)))
                if err<1e-8:
                    self.max_closure_error=max(self.max_closure_error,err)
                    self.max_closure_iterations=max(self.max_closure_iterations,iteration)
                    return answer
                correction,photons=new_c,new_p
            raise AssertionError(('Coupled metric/transport closure did not converge',float(t),err))
        return step


def symbolic():
    import sympy as s
    r,m,p,A,phi,psi,psip,J,dp,alpha=s.symbols('r m p A phi psi psip J dp alpha',nonzero=True)
    b=1-2*m/r;dm=r*r*b*phi*psi+J
    nu=m/(r*r*b)+4*s.pi*r*A*p/b+r*phi**2/2
    variation=s.diff(nu,m)*dm+s.diff(nu,p)*dp+s.diff(nu,A)*4*alpha*A*psi+s.diff(nu,phi)*psip
    target=(1+8*s.pi*r*r*A*p)*dm/(r*r*b*b)+4*s.pi*r*A*(dp+4*alpha*p*psi)/b+r*phi*psip
    assert s.simplify(variation-target)==0
    phip=-(2/r+2*m/(r*r*b))*phi
    assert s.simplify(phi+r*phip-phi/b+2*phi/b)==0
    u,Pg,Pr,E,F=s.symbols('u Pg Pr E F')
    # Energy and normal momentum supplied by the explicit pressure support.
    assert s.expand((F+u*(Pg+Pr))-(F+u*Pr)-u*Pg)==0
    assert s.expand((Pg+Pr+u*F)-(Pr+u*F)-Pg)==0
    return dict(classification='Proven',passed=True,
        lapse='delta_nu_prime=(1+8*pi*r^2*A4*p)*delta_m/(r^2*b^2)+4*pi*r*A4*(delta_p+4*alpha*p*psi)/b+r*Phi*psi_prime. Exterior integration gives nu_s=Phi_s*psi_s-2*integral(Phi*psi/b)-integral[J/(r^2*b^2)+4*pi*r*Pr/b].',
        transport='delta L contains op*(delta lnT+alpha*psi+delta nu) plus L0*(delta nu+2*alpha*psi-delta m/(r*b)). Surface delta L/L0=4*Delta lnT+(2+4*alpha*Phi_s+2*nu_s_prime)*zeta+4*alpha*psi_s+2*delta nu_s.',
        support='Prescribed gas pressure supplies normal momentum Pg and moving-surface energy flux v*Pg. These local matching identities do not supply its gravitational stress tensor.',
        scope='Algebra of the declared tangent equations and fixed transport coefficients. Not an isolated nonlinear atmosphere or a global final-charge theorem.')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='075c83dcb',
        claim='Connect the same material/scalar/metric solution to bulk redshifted heat transport, surface area/lapse luminosity and moving-photon feedback inside every coupled time stage.',
        changed='Reconstruct radial lapse from the same pressure/mass/scalar constraints; include its volume and surface heat-flux terms. Use the actual current stage emission/motion in causal photons and iterate all added terms to1e-8.',
        preserved='Native EOS/baryons, composition, C1 source, p4 consistent mass, coefficients, horizon, finite-pressure mechanical boundary, original numerical gates. No new background roots or native EOS calls.',
        boundary='The externally maintained Pg is prescribed; its required work and momentum are stated. Its gravitational stress tensor and the full radiation-corrected scalar junction remain unclosed. No claim of an isolated stellar surface.',
        paths=['p4-8 pilot','p4-32','p4-64'],gates=dict(temperature=.02,velocity=.03,scalar=.03,closure=1e-8),
        budget=dict(total_seconds=240,CPU_threads=1,memory_GB=5,maximum_closure_iterations=6),
        decision='Stop on failed closure/tangent/residual or measured budget. No degree, mesh, period, profile or threshold expansion. Compare the saved same-order Phase96 path without relabeling its unresolved spatial error.',
        symbolic=symbolic(),bindings={str(p.relative_to(old.ROOT)):old.photons.digest(p) for p in [Path(__file__),Path(prior.__file__),Path(moving.__file__),prior.OUT/'p4-64.npz',old.OUT/'inputs.npz',old.OUT/'coefficients.npz']}))
    (OUT/'reused-stage.py').write_text(stage_source,encoding='utf-8');(OUT/'reused-evolution.py').write_text(evolve_source,encoding='utf-8')


def run():
    assert not (OUT/'result.json').exists();spec=json.loads((OUT/'plan.json').read_text())
    for p,h in spec['bindings'].items():assert old.photons.digest(old.ROOT/p)==h,p
    signal.alarm(230);resource.setrlimit(resource.RLIMIT_AS,(int(5e9),int(5e9)));start=time.monotonic()
    m=Model();horizon=json.loads((old.OUT/'plan.json').read_text())['horizon_seconds']
    pilot=m.evolve(horizon,8,'p4-8');forecast=time.monotonic()-start+1.5*(pilot['seconds']*96/8+10)
    write(OUT/'pilot-budget.json',dict(classification='Counterexample candidate',forecast_seconds=forecast,setup_seconds=m.setup_seconds,evolution8_seconds=pilot['seconds'],
        assumption='Actual simultaneous-closure p4 pilot scaled to96 further steps with50percent margin; startup/import is outside internal timing.'))
    assert forecast<230,'Measured closure path exceeds budget; no automatic expansion'
    for n in [32,64]:m.evolve(horizon,n,f'p4-{n}')
    a=np.load(OUT/'p4-32.npz');b=np.load(OUT/'p4-64.npz');c=np.load(prior.OUT/'p4-64.npz');errors={};change={}
    for f in ['temperature','velocity','scalar']:
        scale=max(np.max(abs(b[f])),1e-100)
        errors[f]=float(np.max(abs(a[f]-b[f][::2]))/scale);change[f]=float(np.max(abs(c[f]-b[f]))/scale)
    write(OUT/'result.json',dict(classification='Counterexample candidate',time_passed=all(errors[f]<spec['gates'][f] for f in errors),
        time_relative=errors,old_new_same_p4_relative=change,seconds=time.monotonic()-start,
        actual_bulk_metric_heat_feedback=True,actual_surface_metric_luminosity_feedback=True,actual_moving_ray_feedback=True,
        external_support_gravity_closed=False,radiation_corrected_scalar_junction_closed=False,spatial_failure_resolved=False,final_charge_solved=False,full_goal_complete=False))
    signal.alarm(0)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','run']);globals()[parser.parse_args().action]()
