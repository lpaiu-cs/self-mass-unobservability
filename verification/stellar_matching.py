"""Specified SLy/DEF stellar equilibrium and scalar response; no inferred EOS.

G=c=1, lengths in metres. The scalar background is zero. Scalar perturbations
therefore decouple from fluid/metric perturbations at first order. Frequencies
use exp(-i omega t). See notes/REQUEST13_STELLAR_DERIVATION.md.
"""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq, root

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'outputs/research-remediation'
GM_SUN = 1.3271244e20
C = 299792458.
MSUN = GM_SUN/C**2
TARGET_MASS = 1.4378144085


class Star:
    def __init__(self, logpc, rtol=2e-9, r0=.5, interpolation='pchip'):
        path = OUT/'sources/LALSimNeutronStarEOS_SLY.dat'
        table = np.loadtxt(path)
        assert table.shape[1] == 2 and np.all(np.diff(table,axis=0)>0)
        self.eos = PchipInterpolator(np.log(table[:,0]),np.log(table[:,1]),extrapolate=False)
        self.interpolation=interpolation
        if interpolation=='lal':
            import lalsimulation as ls
            self.native_eos=ls.SimNeutronStarEOSFromFile(str(path))
            self.native_energy=ls.SimNeutronStarEOSEnergyDensityOfPressureGeometerized
        self.pmin,self.pmax = table[0,0],table[-1,0]
        self.rtol = rtol
        self.r0 = r0
        self.pc = np.exp(logpc)
        assert self.pmin < self.pc < self.pmax
        ec = self.energy(self.pc)

        def rhs(r,y):
            m,lp,nu = y
            p = np.exp(max(lp,np.log(self.pmin)))
            e = self.energy(p)
            f = 1-2*m/r
            assert f > 0
            nup = (m+4*np.pi*r**3*p)/(r*r*f)
            return [4*np.pi*r*r*e,-(e+p)/p*nup,nup]

        def surface(r,y):
            return y[1]-np.log(self.pmin)
        surface.terminal=True
        surface.direction=-1
        p0=self.pc-2*np.pi/3*(ec+self.pc)*(ec+3*self.pc)*r0**2
        self.sol=solve_ivp(rhs,(r0,40000.),[4*np.pi/3*ec*r0**3,np.log(p0),0.],
            method='DOP853',rtol=rtol,atol=[1e-8,1e-11,1e-12],events=surface,dense_output=True,max_step=200.)
        assert self.sol.success and len(self.sol.t_events[0]) == 1
        self.radius=float(self.sol.t[-1]); self.mass=float(self.sol.y[0,-1])
        self.nu_shift=.5*np.log(1-2*self.mass/self.radius)-self.sol.y[2,-1]

    def energy(self,p):
        if self.interpolation=='lal':
            assert self.pmin <= p <= self.pmax
            return self.native_energy(p,self.native_eos)
        value=self.eos(np.log(p))
        assert np.isfinite(value), 'EOS extrapolation forbidden'
        return float(np.exp(value))

    def interior(self,r,beta,w):
        m,lp,nu=self.sol.sol(r)
        p=np.exp(max(lp,np.log(self.pmin))); e=self.energy(p)
        f=1-2*m/r
        nup=(m+4*np.pi*r**3*p)/(r*r*f)
        lamp=(4*np.pi*r*e-m/r**2)/f
        return 2/r+nup-lamp,(w*w*np.exp(-2*(nu+self.nu_shift))+4*np.pi*beta*(-e+3*p))/f

    def scalar(self,beta,w=0.):
        def rhs(r,y):
            a,b=self.interior(r,beta,w)
            return [y[1],-a*y[1]-b*y[0]]
        ec=self.energy(self.pc)
        b=w*w*np.exp(-2*self.nu_shift)+4*np.pi*beta*(-ec+3*self.pc)
        r0=self.r0
        sol=solve_ivp(rhs,(r0,self.radius),[1-b*r0*r0/6,-b*r0/3],method='DOP853',
            rtol=self.rtol,atol=1e-12,max_step=100.)
        assert sol.success
        return sol.y[:,-1]

    def static(self,beta):
        psi,dpsi=self.scalar(beta)
        r=self.radius; m=self.mass; f=1-2*m/r
        q=-r*r*f*dpsi
        phi=psi+q/(2*m)*np.log(f)
        return dict(beta=beta,central_normalization_phi_infinity=float(phi),
                    charge_per_central_scalar_m=float(q),susceptibility_m=float(q/phi))

    def scattering(self,beta,omega_s,far_cycles=30.,relative_response=False):
        # Exterior equation for u=r*phi in tortoise radius x. Integrate in r:
        # f^2 u''+f f' u'+(w^2-f*2M/r^3)u=0.
        w=omega_s/C; m=self.mass; r=self.radius
        far=far_cycles/w
        def rhs(t,y):
            # log radius avoids resolving the entire exterior on a stellar mesh.
            rr=np.exp(t); f=1-2*m/rr
            u,v=y  # v=du/d(log r)
            return [v,(1-2*m/(rr*f))*v-((w*rr/f)**2-2*m/(rr*f))*u]
        ffar=1-2*m/far
        # exp(i*w*x) outgoing, with an irrelevant overall phase removed.
        sol=solve_ivp(rhs,(np.log(far),np.log(r)),[1+0j,1j*w*far/ffar],method='DOP853',
            rtol=min(self.rtol,2e-11),atol=2e-13)
        assert sol.success
        out,logder=sol.y[:,-1]
        dout=logder/r
        phi,dphi=self.scalar(beta,w)
        u=r*phi; du=phi+r*dphi
        if relative_response:
            phi0,dphi0=self.scalar(0.,w)
            y=du/u; y0=(phi0+r*dphi0)/(r*phi0)
            # Exact exterior Wronskian is 2 i w/f(R). Cancel it analytically
            # before dividing S/S0-1 by 2 i w, avoiding small-frequency loss.
            response=(y-y0)/((1-2*m/r)*(dout-out*y)*(np.conj(out)*y0-np.conj(dout)))
            return response
        # Aout/Ain from a 2x2 Wronskian; incoming is conjugate outgoing.
        ratio=(np.conj(out)*du-np.conj(dout)*u)/(dout*u-out*du)
        # Restore the removed far phase to use exp(i*w*r*) at infinity.
        xfar=far+2*m*np.log(far/(2*m)-1)
        return -ratio*np.exp(-2j*w*xfar)

    def outgoing_pole(self,beta,guess,angle=1.8,length=2e6):
        # Exterior complex radial contour. The outgoing solution decays at the
        # remote end; Riccati integration avoids exponential amplitude overflow.
        m=self.mass; rr=self.radius; direction=np.exp(1j*angle)
        def mismatch(x):
            omega=complex(*x); w=omega/C
            far=rr+direction*length; ffar=1-2*m/far
            def rhs(s,y):
                r=rr+direction*s; f=1-2*m/r
                v=f*2*m/r**3
                return direction*(-y*y-2*m/(r*r*f)*y-(w*w-v)/(f*f))
            ext=solve_ivp(rhs,(length,0.),[1j*w/ffar],method='DOP853',rtol=2e-10,atol=1e-13)
            assert ext.success
            phi,dphi=self.scalar(beta,w)
            z=rr*((phi+rr*dphi)/(rr*phi)-ext.y[0,-1])
            return [z.real,z.imag]
        fit=root(mismatch,guess,tol=1e-8)
        error=float(np.linalg.norm(mismatch(fit.x)))
        return dict(omega_real_per_s=float(fit.x[0]),omega_imag_per_s=float(fit.x[1]),
                    angle=angle,outer_length_m=length,success=bool(fit.success),dimensionless_mismatch=error)


def prototype():
    rows=[]
    for x in [.9,1.18,1.2,1.3,1.34,1.35,1.4,1.5]:
        dx=x-1.3474
        c2=dx*(-.3796+.3294*x-.09174*x*x)
        c4=4.270-7.804*x+5.545*x*x-1.337*x**3
        if dx < 0:
            wr=dx*dx*(.4465-.9698*x+.6670*x*x)
            wi=dx*(-.20769+.21629*x-.11647*x*x)
            k=c2; q=0.
        else:
            wr=dx*dx*(18.753-20.205*x+5.374*x*x)
            wi=dx*(1.5439-1.2407*x+.20706*x*x)
            q=np.sqrt(-6*c2/c4); k=-2*c2
        inertia=k/(wr*wr+wi*wi)
        # This damping reproduces the complex pole, rather than assuming Gamma=1.
        gamma=2*inertia*wi
        assert k>0 and inertia>0 and wi>0
        rows.append(dict(baryon_mass_solar=x,c2=c2,c4=c4,charge_solar=q,
            omega_real_per_solar_mass=wr,omega_damping_per_solar_mass=wi,
            inertia_solar_mass=inertia,pole_matched_damping=gamma,
            decay_seconds=MSUN/C/wi,low_frequency_tau_seconds=gamma/k*MSUN/C))
    control=next(r for r in rows if r['baryon_mass_solar']==1.2)
    assert np.isclose(control['c2'],.01716,rtol=.002)
    assert np.isclose(control['c4'],.5791,rtol=.002)
    assert np.isclose(control['inertia_solar_mass'],53.7,rtol=.03)
    return rows


def main():
    rows=prototype()
    results=[]
    for interpolation,rtol,r0 in [('pchip',2e-8,1.),('pchip',2e-10,.25),('lal',2e-8,1.),('lal',2e-10,.25)]:
        pc=brentq(lambda lp:Star(lp,rtol,r0,interpolation).mass/MSUN-TARGET_MASS,
                  np.log(1e-11),np.log(2e-10),xtol=2e-10)
        star=Star(pc,rtol,r0,interpolation)
        assert abs(star.mass/MSUN-TARGET_MASS)<1e-6
        control=star.static(0.)
        assert abs(control['susceptibility_m'])<1e-9
        stat=star.static(-4.)
        result=dict(interpolation=interpolation,rtol=rtol,r0_m=r0,central_pressure_geom=star.pc,
            gravitational_mass_solar=star.mass/MSUN,radius_m=star.radius,
            compactness=star.mass/star.radius,static=stat)
        results.append(result)
        print('STAR',result,flush=True)
    # Finite-frequency comparison, independently subtracting the beta=0 geometry.
    freq=[]
    for far in [30.,60.]:
        for w in [1.,3.,10.,30.,100.,300.]:
            s=star.scattering(-4.,w,far); s0=star.scattering(0.,w,far)
            raw_response=(s-s0)/(2j*w/C)
            # Remove the beta=0 Schwarzschild scattering phase. Subtracting S0
            # alone retains a frequency-dependent propagation phase in H.
            response=(s/s0-1)/(2j*w/C)
            stable=star.scattering(-4.,w,far,True)
            assert abs(stable-response)/abs(stable)<2e-8
            assert abs(abs(s/s0)-1)<1e-9, 'Elastic flux conservation'
            freq.append(dict(omega_per_s=w,far_cycles=far,
                susceptibility_real_m=response.real,susceptibility_imag_m=response.imag,
                raw_response_real_m=raw_response.real,raw_response_imag_m=raw_response.imag,
                inferred_inverse_real_per_m=(1/response).real,
                inferred_damping=-(1/response).imag/(w/C)))
            print('RESPONSE',freq[-1],flush=True)
    # Independent public implementation uses a different table interpolant.
    import lal
    import lalsimulation as lalsim
    eos=lalsim.SimNeutronStarEOSFromFile(str(OUT/'sources/LALSimNeutronStarEOS_SLY.dat'))
    radius,mass,love=lalsim.SimNeutronStarTOVODEIntegrate(star.pc*lal.C_SI**4/lal.G_SI,eos)
    independent=dict(radius_m=radius,gravitational_mass_solar=mass/lal.MSUN_SI,love_k2=love,
        radius_relative_difference=radius/star.radius-1,mass_relative_difference=mass/lal.MSUN_SI/(star.mass/MSUN)-1)
    print('LAL_CROSSCHECK',independent,flush=True)
    dense=[]
    for size in [4096,16384]:
        # Sample the SAME native p->e continuum. This diagnoses the enthalpy
        # construction, rather than changing the EOS or selecting another mass.
        ps=np.geomspace(star.pmin,star.pmax,size)
        es=np.array([star.energy(p) for p in ps])
        path=OUT/f'SLy-native-dense-{size}.dat'
        np.savetxt(path,np.column_stack([ps,es]))
        refined=lalsim.SimNeutronStarEOSFromFile(str(path))
        rr,mm,kk=lalsim.SimNeutronStarTOVODEIntegrateWithTolerance(star.pc*lal.C_SI**4/lal.G_SI,refined,1e-9)
        row=dict(points=size,radius_m=rr,gravitational_mass_solar=mm/lal.MSUN_SI,
                 radius_relative_difference=rr/star.radius-1,mass_relative_difference=mm/lal.MSUN_SI/(star.mass/MSUN)-1)
        dense.append(row)
        print('LAL_DENSE',row,flush=True)
    poles=[]
    for angle,length in [(1.8,2e6),(2.,2e6),(1.8,4e6)]:
        pole=star.outgoing_pole(-4.,[1500.,-6500.],angle,length)
        poles.append(pole)
        print('POLE',pole,flush=True)
    orbital=[]
    for period in [1.6293990080893948,327.25512703377643]:
        w=2*np.pi/(period*86400)
        response=star.scattering(-4.,w,60.,True)
        # For this conservative one-channel scattering problem flux conservation
        # fixes the imaginary inverse exactly; use it to remove cancellation noise.
        k=(1/response).real
        flux_response=1/(k-1j*w/C)
        orbital.append(dict(period_days=period,real_m=flux_response.real,imag_m=flux_response.imag,
            raw_complex_imag_m=response.imag,phase_radians=float(np.angle(flux_response))))
    data=dict(status='Counterexample candidate',theory='massless DEF, beta=-4, phi_infinity=0, unscalarized branch',
        scope='Specified-EOS structure and linear scalar scattering; not companion matching or a J0337 gravity limit',
        eos_sha256=hashlib.sha256((OUT/'sources/LALSimNeutronStarEOS_SLY.dat').read_bytes()).hexdigest(),
        prototype=rows,stars=results,frequency_response=freq,independent_LAL=independent,independent_LAL_refinement=dense,poles=poles,orbital_response=orbital)
    (OUT/'stellar-matching.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')


if __name__=='__main__': main()
