"""Conservative stress packet and charge-source ledger for the native fan.

Counterexample candidate: the nonlinear local gas fan is mapped to the leading
thin-layer GR constraints. This is not a global dynamical charge solution.
"""
from pathlib import Path
import json
import signal
import time
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.integrate import cumulative_simpson
import def_native_vacuum_release as task

OUT=task.OUT


def identities():
    import sympy as s
    rho,h,p,v,c,x,e0=s.symbols('rho h p v c x e0',real=True)
    W2=1/(1-v*v);E=rho*h*W2-p;S=rho*h*W2*v
    radial=S*v+p;trace=rho*h-4*p
    assert s.simplify(E-radial-2*p-trace)==0
    assert s.simplify(trace-E+S*v+3*p)==0
    xi=(v-c)/(1-v*c)
    assert s.simplify((xi*E-S)-(-rho*h*c/(1-v*c)-xi*p))==0
    # Self-similar conservation F'=xi*U' gives (xi*U-F)'=U.
    U,F=s.Function('U')(x),s.Function('F')(x)
    assert s.simplify(s.diff(x*U-F,x).subs(s.diff(F,x),x*s.diff(U,x))-U)==0
    return dict(classification='Proven',passed=True,
        energy_primitive='K_E(xi)=xi*(E-E0*H(-xi))-S. dK_E/dxi=E-E0*H(-xi). It vanishes at both ends of a complete conservative vacuum fan.',
        trace_integral='Integral delta(e-3p) dxi = -Integral[S*v+3*(p-p0*H(-xi))] dxi, if total lab energy is conserved. Never multiply a raw baryon quadrature residual by c^2 and call it charge.',
        general_weight='Integral w*delta(trace) = [w*K_E] - Integral w_prime*K_E - Integral w*(S*v+3*delta p).',
        mass_constraint='In the local thin-layer limit, delta m_geom=(G/c^4)*A_s*sqrt(b_s)*area_s*c*tau*K_E; its time derivative obeys the radial energy-flux constraint.',
        scope='Conservation identities and leading local constraint map, not an error-certified spherical solution or Bondi scalar charge.')


def synthetic():
    # An independent exact gamma-law SR invariant, with appreciable v/c,
    # checks the same stress quadrature used for the native, much slower gas.
    x=np.linspace(0,-18,145);rho=np.exp(x);g=5/3;p=.01*rho**g;u=p/(rho*(g-1));h=1+u+p/rho
    cs=np.sqrt(g*p/(rho*h));Y=2/np.sqrt(g-1)*(np.arctanh(cs[0]/np.sqrt(g-1))-np.arctanh(cs/np.sqrt(g-1)))
    v=np.tanh(Y);xi=(v-cs)/(1-v*cs);raw=np.zeros((len(x),21));raw[:,0]=rho;raw[:,1]=p;raw[:,2]=u
    d=dict(xi=xi,log_density_ratio=x,raw=raw,velocity_over_c=v,enthalpy=h)
    moments=task.integrate_fan(d)
    errors=dict(baryon=abs(moments['baryon'])/(xi[-1]-xi[0]),
        energy=abs(moments['nonrest_energy'])/(.01*(xi[-1]-xi[0])),momentum=abs(moments['momentum']/.01-1))
    assert max(errors.values())<.002,errors
    return errors


def packet(d,env,R,c,t):
    xi=d['xi'];rho=d['raw'][:,0];p=d['raw'][:,1];u=d['raw'][:,2];v=d['velocity_over_c'];cs=d['cs']
    h=d['enthalpy'];E=rho*h/(1-v*v)-p;S=rho*h*v/(1-v*v);left=xi<0
    KE=-rho*h*cs/(1-v*cs)-xi*p-xi*E[0]*left
    A=float(env['A'][-1]);N=float(env['N'][-1]);b=float(env['b'][-1]);tau=A*N*t
    area=4*np.pi*(A*R)**2;scale=c*tau
    radius=R+np.sqrt(b)/A*scale*xi
    deltaPr=S*v+p-p[0]*left
    dm=task.old.G/c**4*A*np.sqrt(b)*area*scale*KE
    # Fixed surface coefficients, not a fictitious exact spherical EOS solve.
    lapse_prime=dm/(R*R*b*b)+4*np.pi*task.old.G/c**4*R*A**4/b*deltaPr
    integral=cumulative_simpson(lapse_prime,x=radius,initial=0)
    lapse=integral-integral[-1]
    np.savez_compressed(OUT/'gr-source-packet.npz',xi=xi,radius_cm=radius,proper_time_s=tau,
        local_density=rho,local_temperature=d['T'],local_velocity_cm_s=v*c,
        lab_energy=E,lab_momentum=S,radial_pressure=S*v+p,tangential_pressure=p,
        delta_trace=rho*(h-p/rho)-3*p-(rho[0]*(h[0]-p[0]/rho[0])-3*p[0])*left,
        conservative_energy_primitive=KE,delta_m_geom_cm=dm,delta_lapse=lapse)
    return dict(maximum_abs_delta_m_geom_cm=float(max(abs(dm))),
        maximum_abs_delta_lapse=float(max(abs(lapse))),
        energy_primitive_at_resolved_tail=float(KE[-1]),
        full_tail_mass_constraint_closed=False,
        scope='The tail endpoint lapse is set to zero only for this resolved source packet. The unqueried tail, variable coefficients, spherical fluid backreaction and photon exchange remain unclosed.')


def main():
    assert not (OUT/'audit.json').exists();assert not (OUT/'audit-plan.json').exists()
    task.write(OUT/'audit-plan.json',dict(classification='Counterexample candidate',
        claim='Construct the conservative nonlinear stress input for the actual stellar radius and its leading GR mass/lapse response; adjudicate the resolved integrated trace without rest-mass quadrature contamination.',
        gates=dict(exact_gamma_law_stress_moments=.002,native_momentum_relative=.002,conservative_trace_refinement=.01),
        budget_seconds=30,native_calls=0,new_global_evolutions=0,
        tail_premises='In addition to the registered sound-speed envelope, assume0<=u<=u_cut and0<=p/rho<=p_cut/rho_cut below the last native state. This explicit premise supplies a stress bound, not a native cold-EOS certification.',
        bindings={str(p.relative_to(task.old.ROOT)):task.old.photons.digest(p) for p in [Path(__file__),Path(task.__file__),OUT/'coarse.npz',OUT/'fine.npz',OUT/'result.json']}))
    start=time.monotonic();signal.alarm(30);symbolic=identities();control=synthetic()
    d=np.load(OUT/'fine.npz');a=np.load(OUT/'coarse.npz');result=json.loads((OUT/'result.json').read_text())
    mf=task.integrate_fan(d);mc=task.integrate_fan(a)
    # This subtraction removes only the energy integral by the exact
    # conservation identity; it is not a post-hoc fit of a spurious mass mode.
    tf=mf['trace_minus_baryon_rest']-mf['nonrest_energy'];tc=mc['trace_minus_baryon_rest']-mc['nonrest_energy']
    refinement=abs(tf-tc)/abs(tf);momentum_error=abs(mf['momentum']/d['raw'][0,1]-1)
    env=np.load(task.old.prior.OUT/'final-envelope.npz');R=float(env['r'][-1]);A=float(env['A'][-1]);N=float(env['N'][-1]);alpha=-4*float(env['phi'][-1])
    c=task.old.C;area=4*np.pi*(A*R)**2;length=c*result['local_proper_time_seconds']
    model=task.old.Background();adm=model.M*R
    rho,p,u=d['raw'][-1,:3];assert u>=0
    vmax=np.tanh(d['rapidity'][-1]+result['conditional_tail_rapidity_upper']);span=vmax-d['xi'][-1]
    hrest=d['enthalpy'][0]-d['raw'][0,2]-d['raw'][0,1]/d['raw'][0,0]
    tail=(rho*(hrest+u+p/rho)*vmax*vmax/(1-vmax*vmax)+3*p)*span
    conversion=alpha*A*N*task.old.G/c**4*area*length/adm
    source=conversion*tf;tail_bound=abs(conversion)*tail
    naive=conversion*(mf['trace_minus_baryon_rest']+hrest*mf['baryon'])
    radial=packet(d,env,R,c,result['coordinate_horizon_seconds'])
    passed=refinement<.01 and momentum_error<.002
    task.write(OUT/'audit.json',dict(classification='Counterexample candidate',passed=bool(passed),
        symbolic=symbolic,exact_gamma_law_relative=control,native_momentum_relative=momentum_error,
        conservative_trace_refinement_relative=refinement,
        conservative_resolved_trace_integral=tf,conditional_tail_trace_integral_upper=float(tail),
        resolved_trace_energy_erg=area*length*tf,
        instantaneous_normalized_scalar_source_proxy=source,
        conditional_tail_proxy_absolute_upper=tail_bound,
        naive_uncorrected_source_proxy=naive,
        numerical_scope='The source proxy freezes alpha*A*N and area at the actual surface. It is a conservative local source moment, NOT an evolved asymptotic charge; coarse/fine agreement is numerical evidence, not a rigorous discretization enclosure.',
        tail_scope='Tail interval is conditional on the explicit sound-speed and thermodynamic premises; no native EOS calls below rho/rho0=exp(-18).',
        gr_source_packet=radial,seconds=time.monotonic()-start,
        full_spherical_GR_evolution=False,photon_matter_exchange_solved=False,
        final_charge_solved=False,full_goal_complete=False))
    signal.alarm(0);print((OUT/'audit.json').read_text(),flush=True)


if __name__=='__main__':main()
