"""Conditional gas free-surface enclosure without extrapolating a failed EOS.

The theorem assumes the declared molecular gas branch persists to vacuum;
it does not certify real-atmosphere radiation, condensation, or plasma error.
"""
from pathlib import Path
import json
import time
import numpy as np
import mpmath as mp
import sympy as sp
from scipy.integrate import solve_ivp
import def_native_atmosphere_density as atmosphere

h=atmosphere.h
OUT=atmosphere.parent.OUT/'surface-enclosure'


def symbolic():
    # Fixed entropy/composition: dh=dP/rho, hydro: P'=-(e+P)(nu'+beta*phi*phi').
    H,rho,nu_r,beta,phi,v=sp.symbols('H rho nu_r beta phi v',nonzero=True)
    assert sp.simplify(-(rho*H)*(nu_r+beta*phi*v)/(rho*H)+nu_r+beta*phi*v)==0
    z,alpha,phib=sp.symbols('z alpha phib')
    assert sp.expand(z+beta*((phib+alpha*z)**2-phib**2)/2-((1+beta*alpha*phib)*z+beta*alpha**2*z**2/2))==0
    x,mu,q0,N0,ph0=sp.symbols('x mu q0 N0 ph0',positive=True)
    a,q,f,n=sp.symbols('a q f n',real=True);N=N0*sp.exp(mu*n);b=1-2*mu*a/x
    F=sp.Matrix([q0*q0*q*q/(2*mu*N*N*x*x),0,q0*q/(ph0*N*sp.sqrt(b)*x*x),
        a/(x*x*b)+q0*q0*q*q/(2*mu*N*N*b*x**3)])
    J=F.jacobian([a,q,f,n])
    assert all(J[i,2]==0 for i in range(4))
    # These derivative expressions are bounded explicitly in main's state box.
    assert sp.simplify(J[3,0]-(1/(x*x*b*b)+q0*q0*q*q/(N*N*b*b*x**4)))==0
    return dict(classification='Proven',passed=True,first_integral=True,Just_quadratic=True,
        vacuum_jacobian=sp.sstr(J),scope='Algebraic identities; numerical input and global gas-branch premises remain separate.')


def just_surface(ctx,mu,qr,phib,H,R):
    # Source Just map, expressed through the lapse increment z. Avoid subtracting
    # nearly equal radius/enthalpy values and the small quadratic root.
    beta=ctx.mpf(-4);b=1-2*mu;d=qr*qr+2*mu/b;alpha=2*qr/d
    invlam=ctx.sqrt(1+alpha*alpha);lam=1/invlam;t=invlam*d/(d+2)
    taub=-ctx.log((1+t)/(1-t));k=1+beta*alpha*phib
    z=2*H/(k+ctx.sqrt(k*k+2*beta*alpha*alpha*H))
    tauf=taub+2*z/lam
    ratio=ctx.exp((1-lam)*(tauf-taub)/2)*(-ctx.expm1(taub))/(-ctx.expm1(tauf))
    ff=ctx.exp(tauf);bf=ff*(1+(1-lam)*ctx.expm1(-tauf)/2)**2
    return R*ratio,z,phib+alpha*z,R*ratio*(1-bf)/2


def main():
    assert not OUT.exists();OUT.mkdir();start=time.monotonic();mp.mp.dps=75;mp.iv.dps=65
    paths=[Path(__file__),atmosphere.OUT/'table.json',atmosphere.OUT/'result.json',
        atmosphere.parent.BACKGROUND/'background-0.001.npz',atmosphere.parent.BACKGROUND/'lapse.npz',
        h.OLD/'molecular-state-17-8.npz',h.molecular.model.OUT/'direct_ion_bridge.f90',
        h.ROOT/'verification/gr_scalar_nonlinear_exterior.py']
    h.write(OUT/'plan.json',dict(classification='Counterexample candidate',
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths},
        claim='Close the mechanical material-vacuum radius/junction for an explicitly assumed dilute gas completion using h*A*N and a uniform tail self-gravity bound. Do not extend native NaN states or fit a polytrope.',
        premises=['Same fixed entropy and composition extend smoothly to P=rho=T=0; no phase transition or external radiation bath.',
            'e+P=rho*h, h>=h0>0, 0<=3P<=e, and P decreases outward.',
            'Dilute ground state is H2 plus neutral atoms of the declared FreeEOS species; nonideal terms vanish. The saved bridge subtracts the H2 ground energy per neutral gram, which must be multiplied by CX for baryon-gram h0.',
            'Static DEF beta=-4 equations, no heat flux, no surface layer. This is mechanical Cauchy data, not a stationary thermal star.'],
        gates=dict(conditional_surface_relative_halfwidth=1e-8,added_baryon_fraction_upper=1e-12,vacuum_ODE_radius_difference_m=.001),
        budget=dict(hard_timeout_seconds=30,native_calls=0,evolution_steps=0,automatic_expansion=False),
        method='Uniform weighted-state Gronwall bound around an exact Just vacuum test-atmosphere; bound all omitted material sources by the base pressure. No native curve interpolation enters this surface location.'))
    data,_=h.inputs();body=h.Structure(.001);bg=np.load(atmosphere.parent.BACKGROUND/'background-0.001.npz');lap=np.load(atmosphere.parent.BACKGROUND/'lapse.npz')
    y=bg['faces'][0];R=mp.mpf(float(y[0]*body.R));M=mp.mpf(float(y[1]*body.B));mu=M/R
    phib=mp.mpf(float(.001*(1+body.mu*y[3])));qr=mp.mpf(float(.001*body.mu*y[4]*float(R)/body.R))
    nub=mp.mpf(float(lap['nu_faces'][0]));Nb=mp.exp(nub);q0=Nb*mp.sqrt(1-2*mu)*qr
    a=np.array(json.loads((atmosphere.OUT/'table.json').read_text())['rows'][0]['raw']);cx=mp.mpf(float(data['CX'][0]));c=mp.mpf(float(h.gr.C*100))
    h0=cx-cx*mp.mpf(float(a[11]))/(c*c)
    hb=cx+(mp.mpf(float(a[2]))+mp.mpf(float(a[1]))/mp.mpf(float(a[0])))/(c*c)
    H=mp.log1p((hb-h0)/h0);P=mp.mpf(float(h.gr.G))*mp.mpf(float(a[1]))*mp.mpf('.1')/mp.mpf(float(h.gr.C))**4
    # Global state box for both full and vacuum solutions, until either surface.
    xmax=1/(1-H/mu);dx=xmax-1;Amax=mp.mpf(1);bmin=1-2*mu*mp.mpf('1.01')
    Nmin=Nb*mp.exp(-mp.mpf('.1')*mu);Nmax=Nb*mp.exp(mp.mpf('.1')*mu);phmax=phib*mp.mpf('1.01')
    Sa=4*mp.pi*R*R*xmax**4*P/(mu*mu)
    Sq=4*mp.pi*4*phmax*Nmax*R*R*xmax**4*P/(mp.sqrt(bmin)*mu*abs(q0))
    Sn=4*mp.pi*R*R*xmax*P*dx/(mu*bmin)
    Bupper=4*mp.pi*R**3*xmax**4*P/(mp.sqrt(bmin)*mp.mpf(float(body.B))*mu*h0)
    # Absolute row sums of the symbolic vacuum Jacobian throughout the box.
    qmax=mp.mpf('1.01');Lrows=[q0*q0/(mu*Nmin*Nmin)*(qmax+mu*qmax*qmax),0,
        abs(q0)/phib/Nmin/mp.sqrt(bmin)*(mu*qmax/bmin+1+mu*qmax),
        1/bmin**2+q0*q0*qmax*qmax/(Nmin*Nmin*bmin*bmin)+q0*q0*qmax/(mu*Nmin*Nmin*bmin)+q0*q0*qmax*qmax/(Nmin*Nmin*bmin)]
    L=max(Lrows);delta=mp.exp(L*dx)*(Sa+Sq+Sn)
    error_psi=(mu+4*phmax*phib)*delta;radius_bound=R*error_psi*xmax*xmax/mu
    # Closing bootstrap: neither full nor vacuum trajectory can leave the box.
    excursions=[q0*q0*qmax*qmax*dx/(2*mu*Nmin*Nmin)+Sa,Sq,
        abs(q0)*qmax*dx/(phib*Nmin*mp.sqrt(bmin)),
        (mp.mpf('1.01')/bmin+q0*q0*qmax*qmax/(2*mu*Nmin*Nmin*bmin))*dx+Sn]
    assert all(v<limit for v,limit in zip(excursions,[.01,.01,.01,.1])),excursions
    rv,z,pf,mf=just_surface(mp,mu,qr,phib,H,R)
    # Interval evaluation of the closed expression, treating saved inputs as
    # exact binary64 data; not an enclosure of uncertain physical inputs.
    iv=mp.iv
    mui=iv.mpf([mu,mu]);qri=iv.mpf([qr,qr]);pi=iv.mpf([phib,phib]);Hi=iv.mpf([H,H]);Ri=iv.mpf([R,R])
    rvi=just_surface(iv,mui,qri,pi,Hi,Ri)[0]
    lo=np.nextafter(float(rvi.a)-float(radius_bound),-np.inf);hi=np.nextafter(float(rvi.b)+float(radius_bound),np.inf)
    # Independent vacuum IVP, increments avoid background cancellation.
    def rhs(x,v):
        ma=float(mu)+v[0];ph=float(phib)+v[1];nu=float(nub)+v[2];N=np.exp(nu);b=1-2*ma/x
        return [float(q0*q0)/(2*N*N*x*x),float(q0)/(N*np.sqrt(b)*x*x),ma/(x*x*b)+float(q0*q0)/(2*N*N*b*x**3)]
    def event(x,v):return v[2]-2*((float(phib)+v[1])**2-float(phib)**2)-float(H)
    event.terminal=True;event.direction=1
    sol=solve_ivp(rhs,(1,float(xmax)*1.00001),[0,0,0],method='DOP853',rtol=2e-12,atol=[1e-28,1e-23,1e-23],events=event,max_step=.001)
    assert sol.success and len(sol.t_events[0])==1
    independent=float(R)*sol.t_events[0][0];independent_error=abs(independent-float(rv))
    warm=json.loads((atmosphere.OUT/'result.json').read_text())['rows'][-1]
    # Evaluate exact first integral at the saved warm endpoint to expose curve
    # interpolation drift, rather than pretending ODE tolerance is EOS error.
    endpoint=np.load(atmosphere.OUT/'atmosphere-1e-12.npz')['state'][:,-1]
    coldrow=json.loads((atmosphere.OUT/'table.json').read_text())['rows'][-1];aa=np.array(coldrow['raw'])
    hw=cx+(mp.mpf(float(aa[2]))+mp.mpf(float(aa[1]))/mp.mpf(float(aa[0])))/(c*c)
    Hwarm=mp.log(hb/hw);expected_warm=just_surface(mp,mu,qr,phib,Hwarm,R)[0]
    record=dict(classification='Counterexample candidate',conditional_mechanical_surface_passed=bool(float(radius_bound/R)<1e-8 and float(Bupper)<1e-12 and independent_error<.001),
        symbolic=symbolic(),native_calls=0,seconds=time.monotonic()-start,
        photosphere_radius_m=float(R),conditional_surface_radius_m=float(rv),thickness_m=float(rv-R),
        conditional_radius_interval_m=[lo,hi],material_source_radius_halfwidth_m=float(radius_bound),
        total_added_baryon_fraction_upper=float(Bupper),material_mass_fraction_upper=float(Sa),
        scalar_flux_fraction_change_upper=float(Sq),weighted_Gronwall_error=float(delta),vacuum_Lipschitz_bound=float(L),bootstrap_excursions=list(map(float,excursions)),
        molecular_ground_enthalpy=float(h0),base_enthalpy=float(hb),log_enthalpy_drop=float(H),
        vacuum_reference_surface=dict(phi=float(pf),mass_geom_m=float(mf),log_lapse=float(nub+z),scalar_flux_Q_geom_m=float(q0*R)),
        independent_vacuum_ODE_radius_difference_m=independent_error,
        warm_native_segment_first_integral_radius_discrepancy_m=float(warm['radius_m']-expected_warm),
        warm_segment_endpoint_first_integral_defect=float(endpoint[4]-2*((float(phib)+.001*endpoint[2])**2-float(phib)**2)-Hwarm),
        no_surface_layer_junction='P=rho=T=0 and continuous m,phi,Q,N in the assumed gas limit. Lagrangian material condition Delta P=0; scalar matching uses the saved moving-surface adapter.',
        thermal_Cauchy_condition='Zero initial material heat flux at the vacuum endpoint. No radiative equilibrium or orbit-long stationarity inferred.',
        certified_scope='Analytic conditional material-source bound plus interval evaluation at frozen input data. Not an EOS existence/domain proof, physical EOS certification, or dynamic fluid resolvent bound.',
        native_cold_branch_globally_certified=False,physical_radiative_atmosphere=False,thermal_stationarity=False,dynamic_exterior_coupled=False,full_dynamic_charge_solved=False)
    h.write(OUT/'result.json',record);assert record['conditional_mechanical_surface_passed'];print(json.dumps(record),flush=True)


if __name__=='__main__':main()
