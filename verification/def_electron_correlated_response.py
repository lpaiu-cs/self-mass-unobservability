"""Published correlated-ion potential on the frozen electron-response cohort.

The potential is fitted in the degenerate limit. Holding its background
parameters fixed off the Fermi surface is an explicit finite-T extension.
The published ee rate is used only for a DC comparison, never as a full
frequency-dependent electron-electron collision operator.
"""
from pathlib import Path
import argparse
import json
import time
import urllib.request
import resource
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import hbar,m_e,c,k,epsilon_0,elementary_charge,atomic_mass
from scipy.special import expit,exp1
from scipy.integrate import quad
import def_electron_collision_response as base

h=base.h
OUT=base.OUT/'correlated-ions'
alpha=base.E2/(hbar*c)


def distribution(T,ne,order):
    """Four panels at the same total node count as the preserved first trial."""
    t=k*T/(m_e*c*c);xf=hbar*(3*np.pi**2*ne)**(1/3)/(m_e*c)
    eta=xf*xf/(np.sqrt(1+xf*xf)+1)/t
    gx,gw=leggauss(order//4)
    gx=np.concatenate([(gx+1)/8+j/4 for j in range(4)])
    gw=np.tile(gw/8,4)
    top=np.sqrt(eta+60);u=top[:,None]*gx;x=u*u
    dx=2*top[:,None]*gw*u;pp=np.sqrt(t[:,None]*x*(2+t[:,None]*x))
    density=(m_e*c/hbar)**3/np.pi**2*t[:,None]*(1+t[:,None]*x)*pp
    for _ in range(12):
        f=expit(eta[:,None]-x);slope=f*(1-f)
        number=np.sum(density*f*dx,axis=1);derivative=np.sum(density*slope*dx,axis=1)
        error=number/ne-1
        if np.max(abs(error))<2e-13:break
        eta-=(number-ne)/derivative
    assert np.max(abs(error))<2e-11
    return eta,x,pp,slope,dx,derivative/(k*T),float(np.max(abs(error)))


def angular(u,w,v2):
    """Positive integral; analytic form except its cancellation-prone domain."""
    u,w,v2=np.broadcast_arrays(u,w,v2);out=np.empty(u.shape)
    regular=(u<1)&(w>.1)&(u*w<100)
    a,b,v=u[regular],w[regular],v2[regular]
    ew=-np.expm1(-b);ln=np.log1p(1/a)
    e=np.exp(a*b)*(exp1(a*b)-exp1((a+1)*b))
    l1=.5*(ln+a/(a+1)*ew-(1+a*b)*e)
    l2=.5*(1-ew/b-a*a/(a+1)*ew-2*a*ln+a*(2+a*b)*e)
    out[regular]=l1-v*l2
    if np.any(~regular):
        a,b,v=u[~regular],w[~regular],v2[~regular]
        nodes,weights=leggauss(32)
        top=np.log1p(1/a);z=top[:,None]*(nodes+1)/2
        momentum=a[:,None]*np.expm1(z)
        out[~regular]=top/4*np.sum(weights*(-np.expm1(-z))*(1-v[:,None]*momentum)*(-np.expm1(-b[:,None]*momentum)),axis=1)
    assert np.min(out)>0
    return out


def kinetic(indices,order):
    d,data=base.inputs();T=data['T'][indices];ne=data['ne'][indices]
    ions=data['ion'][indices];Z=base.thermal.g.c.Z;A=base.thermal.g.c.A
    eta,x,pdim,fp,dx,dndmu,error=distribution(T,ne,order)
    xf=data['xF'][indices];pf=m_e*c*xf;vf=c*xf/np.sqrt(1+xf*xf)
    tf=hbar*hbar*(elementary_charge**2/epsilon_0*dndmu)
    momentum=m_e*c*pdim;energy=m_e*c*c*np.sqrt(1+pdim*pdim);velocity=momentum*c*c/energy
    total=np.zeros_like(x);nuF=np.zeros(len(T));gamma_all=[];eta_all=[];inelastic=[]
    for j,zj in enumerate(Z):
        if not np.any(ions[:,j]>0):continue
        ai=(3*zj/(4*np.pi*ne))**(1/3)
        gamma=zj*zj*base.E2/(k*T*ai);qd2=3*gamma/(ai*ai)
        beta=np.pi*alpha*zj*vf/c
        qs2=(qd2*(1+.06*gamma)*np.exp(-np.sqrt(gamma))+tf/hbar**2)*np.exp(-beta)
        wc=13*(1+beta/3)/qd2
        tp=hbar/k*np.sqrt(elementary_charge**2/epsilon_0*zj*ne/(A[j]*atomic_mass))
        ion_eta=T/tp;eta0=.19/zj**(1/6)
        G=ion_eta/np.sqrt(ion_eta*ion_eta+eta0*eta0)*(1+.122*beta*beta)
        Gk=G+.0105*(1-1/zj)*(1+(vf/c)**3*beta)*ion_eta/(ion_eta*ion_eta+.0081)**1.5
        alpha0=4*(pf/hbar)**2*ai*ai/(3*gamma*ion_eta)
        D=np.exp(-alpha0*2.8*np.exp(-9.1*ion_eta)/4)
        u=hbar*hbar*qs2[:,None]/(4*momentum*momentum)
        w=wc[:,None]*4*(momentum/hbar)**2
        term=angular(u,w,(velocity/c)**2)*(G*D)[:,None]
        total+=ions[:,j,None]*zj*zj*term
        LF=angular(hbar*hbar*qs2/(4*pf*pf),wc*4*(pf/hbar)**2,(vf/c)**2)
        nuF+=ions[:,j]*zj*zj*LF*Gk*D
        gamma_all.append(gamma);eta_all.append(ion_eta);inelastic.append(Gk/G-1)
    nu=4*np.pi*base.E2**2*total/(momentum*momentum*velocity)
    nuF*=4*np.pi*base.E2**2/(pf*pf*vf)
    tau=1/nu;z=x-eta[:,None]
    W=momentum**3/(3*np.pi**2*hbar**3)*c*c/energy*fp*dx
    weights=W*tau;moments=np.array([np.sum(weights*z**j,axis=1) for j in range(3)])
    mean=moments[1]/moments[0];variance=np.sum(weights*(z-mean[:,None])**2,axis=1)
    conductivity=k*k*T*variance
    memory=np.sum(weights*tau*(z-mean[:,None])**2,axis=1)/variance
    # Published ee thermal rate at the Fermi surface: DC comparison only.
    mass=m_e*np.sqrt(1+xf*xf);ktf2=4*alpha/np.pi*(pf/hbar)**2*c/vf
    tpe=hbar/k*np.sqrt(elementary_charge**2/epsilon_0*ne/mass);y=np.sqrt(3)*tpe/T
    J=(1+6/(5*xf*xf)+2/(5*xf**4))*(y**3/(3*(1+.07414*y)**3)*np.log1p((2.81-.81*(vf/c)**2)/y)+np.pi**5/6*y**4/(13.91+y)**4)
    nuee=3*alpha**2*(k*T)**2/(2*np.pi**3*hbar*mass*c*c)*(2*pf/(hbar*np.sqrt(ktf2)))**3*J
    degenerate_prefactor=np.pi**2*k*k*T*ne/(3*mass)
    return dict(indices=indices,order=order,eta=eta,x=x,z=z,weights=weights,tau=tau,
        moments=moments,variance=variance,conductivity_SI=conductivity,tau_eff=memory,
        density_relative_error=error,native_conductivity_SI=d['K'][indices,1]*1e-5,
        DC_degenerate_ei_conductivity=degenerate_prefactor/nuF,
        DC_degenerate_ei_ee_conductivity=degenerate_prefactor/(nuF+nuee),
        ee_rate=nuee,ei_rate_Fermi=nuF,ion_gamma=np.array(gamma_all),T_over_Tp=np.array(eta_all),
        inelastic_relative_rate_correction=np.array(inelastic))


def prepare():
    assert not OUT.exists();OUT.mkdir()
    url='https://arxiv.org/html/astro-ph/9903127'
    with urllib.request.urlopen(url,timeout=20) as reply:(OUT/'potekhin-et-al1999.html').write_bytes(reply.read())
    paths=[Path(__file__),Path(base.__file__),base.OUT/'plan.json',base.OUT/'result.json',base.thermal.OUT/'coefficients.npz',OUT/'potekhin-et-al1999.html']
    h.write(OUT/'plan.json',dict(classification='Counterexample candidate',
        sources={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths},
        source_url=url,equations='25-31 and classical mixture sum for the correlated-ion potential; 32-33 for the DC-only ee rate.',
        decision='Resolve the factor-2.6 static mismatch by the published physical terms, without fitting the saved conductivity; do not promote a DC ee rate into an energy-resolved dynamic operator.',
        inherited_failures='Base numerical 1e-6 criterion failed by 1.57e-6; static compatibility failed. Both remain immutable.',
        quadrature_change='Same 128/256 total nodes, four panels in sqrt(kinetic energy) to resolve the Fermi window. No automatic resolution increase.',
        model='Published correlated screened effective potential with mixture sum; background potential parameters fixed off the Fermi surface and elastic G_sigma. G_kappa and ee rates enter only the separate degenerate DC comparison. No fitted prefactors.',
        gates=dict(angular_relative=1e-6,energy_relative=1e-6,static_compatibility_relative=.2),
        budget=dict(pilot_cells=16,orders=[128,256],hard_seconds=60,CPU_workers=1,native_calls=0,stellar_time_steps=0,automatic_expansion=False),
        limits='Finite-T extension, mixture prescription and source fits are approximations, not a certified physical error bound; no ee dynamic operator, outer envelope or photon solution.'))
    _,data=base.inputs();indices=np.flatnonzero(data['selected']);sample=indices[np.linspace(0,len(indices)-1,16).astype(int)]
    began=time.monotonic();pilot=kinetic(sample,128);elapsed=time.monotonic()-began
    h.write(OUT/'pilot.json',dict(seconds=elapsed,forecast_seconds=elapsed*len(indices)/len(sample)*3,
        estimate='Linear cell/node scaling for the two orders; excludes import, output and source-retrieval overhead. Hard limit remains 60s.'))
    errors=[]
    for u in [.001,.1,.999,1.,10.,1e6]:
        for w in [1e-8,.01,.1,1.,20.,1e5]:
            for v in [0.,.2,.9]:
                top=np.log1p(1/u)
                f=lambda t:.5*(-np.expm1(-t))*(1-v*u*np.expm1(t))*(-np.expm1(-w*u*np.expm1(t)))
                exact=quad(f,0,top,epsabs=1e-28,epsrel=2e-11,limit=200)[0]
                actual=angular(np.array([u]),np.array([w]),np.array([v]))[0]
                errors.append(abs(actual/exact-1))
    h.write(OUT/'angular-controls.json',dict(classification='Counterexample candidate',max_relative_error=max(errors),passed=max(errors)<1e-6))
    print('PILOT',elapsed,'FORECAST',elapsed*len(indices)/len(sample)*3,'ANGULAR',max(errors),flush=True)


def run():
    assert not (OUT/'result.json').exists();start=time.monotonic()
    plan=json.loads((OUT/'plan.json').read_text())
    for p,sha in plan['sources'].items():assert h.digest(h.ROOT/p)==sha,p
    assert json.loads((OUT/'pilot.json').read_text())['forecast_seconds']<50
    assert json.loads((OUT/'angular-controls.json').read_text())['passed']
    _,data=base.inputs();indices=np.flatnonzero(data['selected'])
    coarse=kinetic(indices,128);fine=kinetic(indices,256)
    differences={q:float(np.max(abs(coarse[q]/fine[q]-1))) for q in ['conductivity_SI','tau_eff']}
    controls=[]
    for factor in [.01,1.,100.]:
        a=base.response(coarse,factor/coarse['tau_eff']);b=base.response(fine,factor/fine['tau_eff'])
        controls.append(float(np.max(abs(a-b))));assert np.min((1+b).real)>0
    ratios={q:[float(np.min(fine[q]/fine['native_conductivity_SI'])),float(np.max(fine[q]/fine['native_conductivity_SI']))] for q in ['conductivity_SI','DC_degenerate_ei_conductivity','DC_degenerate_ei_ee_conductivity']}
    passed=max(differences.values())<1e-6 and max(controls)<1e-6
    np.savez_compressed(OUT/'response.npz',**fine)
    result=dict(classification='Counterexample candidate',cohort_cells=len(indices),
        numerical_gates_passed=passed,quadrature_relative=differences,frequency_quadrature_absolute=controls,
        static_conductivity_ratios=ratios,DC_ee_to_ei_rate_range=[float(np.min(fine['ee_rate']/fine['ei_rate_Fermi'])),float(np.max(fine['ee_rate']/fine['ei_rate_Fermi']))],
        DC_combined_compatibility_passed=max(abs(v-1) for v in ratios['DC_degenerate_ei_ee_conductivity'])<.2,
        T_over_Tp_range=[float(fine['T_over_Tp'].min()),float(fine['T_over_Tp'].max())],
        max_inelastic_relative_rate=float(fine['inelastic_relative_rate_correction'].max()),
        tau_eff_ei_seconds_range=[float(fine['tau_eff'].min()),float(fine['tau_eff'].max())],
        elapsed_seconds=time.monotonic()-start,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        full_ee_collision_operator_identified=False,physical_error_certified=False,whole_star_closure_replaced=False,full_dynamic_charge_solved=False)
    h.write(OUT/'result.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
