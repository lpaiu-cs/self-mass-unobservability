"""Specified finite-background DEF stars and coupled companion response (Request 17)."""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy import constants
from scipy.integrate import cumulative_simpson, solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq, least_squares

from stellar_matching import Star, C, MSUN, GM_SUN

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/nonzero-drive17'
BETA=-4.

def save(name,data):
    OUT.mkdir(exist_ok=True)
    (OUT/name).write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+'\n')

def electron_gas(x):
    """P, total energy, ion rest-energy in SI; stable nonrelativistic series."""
    x=np.asarray(x);me=constants.m_e;hbar=constants.hbar
    pref=me**4*C**5/(24*np.pi**2*hbar**3)
    n=(me*C/hbar)**3*x**3/(3*np.pi**2)
    rest=2*constants.atomic_mass*C**2*n
    pressure=pref*(x*(2*x*x-3)*np.sqrt(1+x*x)+3*np.arcsinh(x))
    kinetic=3*pref*(x*(1+2*x*x)*np.sqrt(1+x*x)-np.arcsinh(x))-me*C*C*n
    from scipy.special import binom
    ps=8*pref*sum(binom(-.5,k)*x**(2*k+5)/(2*k+5) for k in range(8))
    es=24*pref*sum(binom(.5,k)*x**(2*k+3)/(2*k+3) for k in range(1,9))
    return np.where(x<.05,ps,pressure),rest+np.where(x<.05,es,kinetic),rest

class EOS:
    def __init__(self,kind,size):
        self.kind=kind;self.size=size
        if kind=='SLy':
            import lalsimulation as ls
            path=ROOT/'outputs/research-remediation/sources/LALSimNeutronStarEOS_SLY.dat'
            table=np.loadtxt(path);self.native=ls.SimNeutronStarEOSFromFile(str(path))
            self.energy=lambda p:ls.SimNeutronStarEOSEnergyDensityOfPressureGeometerized(p,self.native)
            p=np.geomspace(table[0,0],table[-1,0],size);e=np.array([self.energy(x) for x in p])
            # d log(n) = d epsilon/(epsilon+p); n=(epsilon+p) exp(-h),
            # h=integral dp/(epsilon+p). Its overall normalization cancels at fixed N.
            h=cumulative_simpson(p/(e+p),x=np.log(p),initial=0.)
            rest=(e+p)*np.exp(-h)
            self.rmax=40000.;self.step=150.;self.r0=.25
        else:
            x=np.geomspace(1e-4,30.,size);p,e,rest=electron_gas(x)
            conversion=constants.G/C**4;p*=conversion;e*=conversion;rest*=conversion
            en=PchipInterpolator(np.log(p),np.log(e),extrapolate=False)
            self.energy=lambda pressure:float(np.exp(en(np.log(pressure))))
            self.rmax=2e8;self.step=1e5;self.r0=25.
            h=cumulative_simpson(p/(e+p),x=np.log(p),initial=0.)
        self.pmin=float(p[0]);self.pmax=float(p[-1])
        self.rest=PchipInterpolator(np.log(p),np.log(rest),extrapolate=False)
        self.h_from_lp=PchipInterpolator(np.log(p),h,extrapolate=False)
        self.lp_from_h=PchipInterpolator(h,np.log(p),extrapolate=False)

class ScalarStar(Star):
    """Reuse the frozen linear-response methods only on this class's phi=0 branch."""
    def __init__(self,eos,logpc,phic=0.,rtol=2e-9,r0_factor=1.):
        self.eos17=eos;self.pc=np.exp(logpc);self.pmin=eos.pmin;self.pmax=eos.pmax
        self.rtol=rtol;self.r0=eos.r0*r0_factor;self.energy=eos.energy
        assert self.pmin<self.pc<self.pmax
        r0=self.r0;ec=self.energy(self.pc);ac=BETA*phic;a4=np.exp(2*BETA*phic**2)
        scalar2=4*np.pi/3*a4*ac*(ec-3*self.pc)
        pressure2=-(ec+self.pc)*(4*np.pi*a4*(ec/3+self.pc)+ac*scalar2)
        p0=self.pc+pressure2*r0*r0/2
        nb0=4*np.pi/3*r0**3*np.exp(1.5*BETA*phic**2+eos.rest(logpc))
        def rhs(r,y):
            m,h,nu,phi,psi,nb=y
            lp=float(eos.lp_from_h(max(h,0.)))
            assert np.isfinite(lp) and lp<np.log(eos.pmax), 'EOS domain exceeded'
            p=max(np.exp(lp),eos.pmin);e=self.energy(p)
            f=1-2*m/r;assert f>0
            a4=np.exp(2*BETA*phi**2);alpha=BETA*phi
            dm=4*np.pi*r*r*a4*e+.5*r*r*f*psi*psi
            dn=(m+4*np.pi*r**3*a4*p)/(r*r*f)+.5*r*psi*psi
            dl=(dm/r-m/r**2)/f
            return [dm,-(dn+alpha*psi),dn,psi,
                    -(2/r+dn-dl)*psi+4*np.pi*a4*alpha*(e-3*p)/f,
                    4*np.pi*r*r*np.exp(1.5*BETA*phi**2+eos.rest(np.log(p)))/np.sqrt(f)]
        def surface(r,y):return y[1]
        surface.terminal=True;surface.direction=-1
        self.full=solve_ivp(rhs,(r0,eos.rmax),[4*np.pi/3*r0**3*a4*ec,eos.h_from_lp(np.log(p0)),0.,phic+scalar2*r0*r0/2,scalar2*r0,nb0],
          rtol=rtol,atol=[1e-10,1e-12,1e-13,1e-15,1e-20,1e-10],method='DOP853',events=surface,max_step=eos.step,dense_output=True)
        assert self.full.success and len(self.full.t_events[0])==1
        self.radius=float(self.full.t[-1]);m,lp,nu,phi,psi,nb=self.full.y[:,-1];R=self.radius
        def exterior(u,y):
            m,phi,z,nu=y;f=1-2*m*u/R
            return [-.5*f*z*z/R,-z/R,2*m*z/(R*f),-m/(R*f)-.5*z*z*u/R**2]
        ext=solve_ivp(exterior,(1.,0.),[m,phi,R*R*psi,nu],method='DOP853',rtol=min(rtol,1e-11),atol=[1e-11,1e-16,1e-12,1e-13])
        assert ext.success
        self.mass,self.phi_infinity,z,self.nu_infinity=map(float,ext.y[:,-1]);self.charge=-z;self.baryon=float(nb)
        self.nu_shift=-self.nu_infinity
        # Existing linear methods expect exactly the GR m,logp,nu profile.
        from types import SimpleNamespace
        self.sol=SimpleNamespace(sol=lambda r:np.array([self.full.sol(r)[0],eos.lp_from_h(np.maximum(self.full.sol(r)[1],0.)),self.full.sol(r)[2]]))
    def record(self):
        return dict(gravitational_mass_solar=self.mass/MSUN,radius_m=self.radius,central_pressure_geom=self.pc,
          phi_infinity=self.phi_infinity,charge_m=self.charge,baryon_normalization_m=self.baryon)
    def scalar(self,beta,w=0.):
        assert abs(self.phi_infinity)<1e-20, 'Linear scalar solver requires the zero branch'
        def rhs(r,y):
            a,b=self.interior(r,beta,w)
            return [y[1],-a*y[1]-b*y[0]]
        ec=self.energy(self.pc);b=w*w*np.exp(-2*self.nu_shift)+4*np.pi*beta*(-ec+3*self.pc);r0=self.r0
        sol=solve_ivp(rhs,(r0,self.radius),[1-b*r0*r0/6,-b*r0/3],method='DOP853',rtol=self.rtol,
          atol=1e-12,max_step=self.eos17.step)
        assert sol.success
        return sol.y[:,-1]

def masses():
    # Frozen 4-body IVP masses in geometric metres; first three form the triple.
    from remaining_audit import hex_fraction
    lines=(ROOT/'outputs/validated-variational/ivp.hex').read_text().splitlines()
    vals=[float(hex_fraction(v)) for line in lines[2:] for v in line.split()]
    # The export is 24 phase variables then four dimensionless masses.
    head=[float(hex_fraction(v)) for v in lines[1].split()]
    return np.array(vals[24:28])*head[3]/MSUN

def kepler_geometry():
    from remaining_audit import hex_fraction
    lines=(ROOT/'outputs/validated-variational/ivp.hex').read_text().splitlines()
    h=np.array([float(hex_fraction(v)) for v in lines[1].split()])
    z=np.array([float(hex_fraction(v)) for line in lines[2:] for v in line.split()])
    pos=z[:12].reshape(4,3)*h[3];vel=z[12:24].reshape(4,3)*h[3]/h[4];m=masses()[:3]
    inner=sum(m[:2]);pcenter=(pos[0]*m[0]+pos[1]*m[1])/inner;vcenter=(vel[0]*m[0]+vel[1]*m[1])/inner
    result=[]
    for r,v,mt in [(pos[0]-pos[1],vel[0]-vel[1],inner),(pcenter-pos[2],vcenter-vel[2],sum(m))]:
        mu=GM_SUN*mt;rr=np.linalg.norm(r);hh=np.cross(r,v)
        aa=1/(2/rr-v.dot(v)/mu);ev=np.cross(v,hh)/mu-r/rr;ecc=np.linalg.norm(ev)
        ex=ev/ecc;ey=np.cross(hh/np.linalg.norm(hh),ex)
        assert aa>0 and 0<ecc<1
        result.append(dict(a_m=aa,ecc=ecc,ex=ex.tolist(),ey=ey.tolist(),period_days=2*np.pi*np.sqrt(aa**3/mu)/86400))
    return result

def coupled():
    data=json.loads((OUT/'stars.json').read_text());fine=[r for r in data['records'] if r['precision']=='fine']
    chi=np.array([r['linear_static']['susceptibility_m'] for r in fine]);m=masses()[:3]*MSUN
    geometry=kepler_geometry();pairs=[(0,1),(0,2),(1,2)];summaries=[]
    for size in [64,128]:
        a,b=np.meshgrid(2*np.pi*np.arange(size)/size,2*np.pi*np.arange(size)/size,indexing='ij')
        vectors=[]
        for mean,g in zip([a,b],geometry):
            E=mean.copy()
            for _ in range(12):E-=(E-g['ecc']*np.sin(E)-mean)/(1-g['ecc']*np.cos(E))
            vectors.append(g['a_m']*((np.cos(E)-g['ecc'])[...,None]*g['ex']+np.sqrt(1-g['ecc']**2)*np.sin(E)[...,None]*g['ey']))
        r,R=vectors;f=m[1]/sum(m[:2]);outer_fraction=m[2]/sum(m)
        positions=np.stack([f*r+outer_fraction*R,-(1-f)*r+outer_fraction*R,-sum(m[:2])/sum(m)*R],axis=-2)
        distances=np.stack([np.linalg.norm(positions[...,i,:]-positions[...,j,:],axis=-1) for i,j in pairs],axis=-1)
        K=np.zeros((size,size,3,3))
        for k,(i,j) in enumerate(pairs):K[...,i,j]=K[...,j,i]=1/distances[...,k]
        L=np.diag(1/chi)-K;q=np.linalg.solve(L,np.ones((size,size,3,1)))[...,0] # q/phi
        assert np.max(abs(np.einsum('...ij,...j->...i',L,q)-1))<1e-12
        Delta=np.stack([q[...,i]*q[...,j]/(m[i]*m[j]) for i,j in pairs],axis=-1) # Delta/phi^2
        frozen=np.stack([q[...,0]*chi[j]/(m[0]*m[j]) for j in [1,2]],axis=-1)
        fft=np.fft.fft2(Delta,axes=(0,1))/size**2
        oldfft=np.fft.fft2(frozen,axes=(0,1))/size**2
        rows=[]
        for i,j in [(1,0),(0,1),(1,-1)]:
            amp=2*abs(fft[i,j])*1e-10;oldamp=2*abs(oldfft[i,j])*1e-10
            rows.append(dict(carrier=[i,j],all_pairs_amplitude=amp.tolist(),frozen_companion_amplitude=oldamp.tolist(),ratio_to_frozen_companions=(amp[:2]/oldamp).tolist()))
        summaries.append(dict(grid=size,min_separations_m=distances.min(axis=(0,1)).tolist(),max_separations_m=distances.max(axis=(0,1)).tolist(),
          mean_charge_at_phi_1e5_m=(q.mean(axis=(0,1))*1e-5).tolist(),mean_Delta_at_phi_1e5=(Delta.mean(axis=(0,1))*1e-10).tolist(),
          max_charge_feedback_fraction=np.max(abs(q/chi-1),axis=(0,1)).tolist(),carriers=rows))
        if size==128:np.savez(OUT/'coupled-grid.npz',q_over_phi=q,Delta_over_phi2=Delta,distances_m=distances)
    save('coupled.json',dict(classification='Counterexample candidate',scope='Leading small-background scalar monopoles on specified 3-body Kepler geometry; no live timing or dynamical damping fit',
      chi_m=chi.tolist(),masses_geom_m=m.tolist(),geometry=geometry,pairs=pairs,grids=summaries,finite_background_nonlinear_terms_in_grid=False))

def stars():
    m=masses();print('MASSES',m,flush=True);records=[]
    for precision,size,rtol,r0 in [('base',4096,2e-9,1.),('fine',16384,2e-11,.5)]:
        for k,kind in enumerate(['SLy','WD','WD']):
            eos=EOS(kind,size)
            bracket=(np.log(1e-11),np.log(2e-10)) if kind=='SLy' else tuple(np.log(electron_gas(np.array([.03,3.]))[0]*constants.G/C**4))
            pc=brentq(lambda lp:ScalarStar(eos,lp,rtol=rtol,r0_factor=r0).mass/MSUN-m[k],*bracket,xtol=2e-11)
            zero=ScalarStar(eos,pc,rtol=rtol,r0_factor=r0);stat=zero.static(BETA)
            assert abs(zero.static(0.)['susceptibility_m'])<1e-7
            responses=[];normalization=stat['central_normalization_phi_infinity']
            for field in [1e-5,-1e-5,3e-5,1e-4]:
                def residual(z):
                    obj=ScalarStar(eos,pc+z[0],field*z[1],rtol=rtol,r0_factor=r0)
                    return [np.log(obj.baryon/zero.baryon),obj.phi_infinity/field-1]
                def jac(z):
                    h=1e-4
                    return np.column_stack([(np.array(residual(z+np.eye(2)[i]*h))-residual(z-np.eye(2)[i]*h))/(2*h) for i in range(2)])
                fitted=least_squares(residual,[0.,1/normalization],jac=jac,bounds=([-.01,.5/normalization],[.01,2/normalization]),
                  xtol=1e-10,ftol=1e-10,gtol=1e-12,max_nfev=20)
                err=float(np.max(abs(np.array(residual(fitted.x)))))
                assert err<1e-7,(kind,precision,field,fitted.message,err)
                obj=ScalarStar(eos,pc+fitted.x[0],field*fitted.x[1],rtol=rtol,r0_factor=r0)
                row=dict(target_phi=field,**obj.record(),baryon_relative_error=obj.baryon/zero.baryon-1,
                  q_over_phi_m=obj.charge/obj.phi_infinity,shooting_error=err,solver_success=bool(fitted.success))
                responses.append(row);print('FIELD',precision,k,row,flush=True)
            records.append(dict(precision=precision,index=k,eos=kind,eos_points=size,rtol=rtol,r0_factor=r0,zero=zero.record(),linear_static=stat,finite_background=responses))
            save('stars-partial.json',dict(classification='Counterexample candidate',records=records))
    save('stars.json',dict(classification='Counterexample candidate',theory='DEF beta=-4; specified SLy and cold ideal mu_e=2 electron EOS',masses_solar=m.tolist(),records=records))

if __name__=='__main__':{'stars':stars,'coupled':coupled}[sys.argv[1]]()
