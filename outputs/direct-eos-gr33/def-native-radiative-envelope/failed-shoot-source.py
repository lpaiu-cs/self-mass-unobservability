"""Native EOS envelope matched to an explicit Eddington outgoing boundary.

Counterexample candidate: quasistatic grey transport, not a static spacetime
with nonzero heat flux, non-LTE atmosphere, or fixed-inventory replacement star.
The original inner P,T,r,m,phi,phi',N and composition are kept at the junction.
"""
from pathlib import Path
import argparse
import hashlib
import json
import resource
import signal
import time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import def_free_surface_thermal as prior

ROOT = Path(__file__).resolve().parents[1]
OUT = prior.OUT.parent/'def-native-radiative-envelope'
write = prior.h.write


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def symbolic():
    import sympy as s
    eg, pg, pr, gravity, force = s.symbols('eg pg pr g force')
    gas = -(eg+pg)*gravity+force
    temperature = -gravity-force/(4*pr)
    assert s.expand(gas+4*pr*temperature+(eg+pg+4*pr)*gravity) == 0
    return dict(classification='Proven', passed=True,
        identity='Pgas prime=-(e_gas+Pgas)*g+rho*kappa*A*a*F/c; lnT prime=-g-rho*kappa*A*a*F/(4*Prad*c). Their sum gives Ptotal prime=-(etotal+Ptotal)*g.',
        assumptions='Isotropic E=3Prad=a_rad*T^4, grey flux transport, g=(ln(A*N)) prime; primes use Einstein radius. Hydrostatic/thermal constraints only, not the time-radial Einstein equation.')


class Envelope:
    def __init__(self):
        self.d = dict(np.load(prior.OUT/'coefficients.npz'))
        self.bg = dict(np.load(prior.surface.OUT/'background.npz'))
        self.eos = prior.h.molecular.model.EOS()
        self.opacity = prior.two.radiative.tables.Opacity()
        inputs, _ = prior.h.inputs()
        self.i = 132
        d, i = self.d, self.i
        self.X = d['X'][i]
        assert np.max(abs(d['X'][:i+1]-self.X)) == 0
        self.cx = float(inputs['CX'][i])
        self.c = prior.h.gr.C*100
        self.G = prior.h.gr.G*1000
        self.arad = float(np.load(prior.OUT.parent/'def-photon-matter-split/matter.npz')['compiled_a_rad_cgs'])
        self.sigma = self.arad*self.c/4
        self.R = float(d['radius_cm'][i])
        rb = self.R/(float(self.bg['R'])*100)
        self.m = float(np.interp(rb,self.bg['r'],self.bg['m'])*float(self.bg['R'])*100)
        self.phi = float(np.interp(rb,self.bg['r'],self.bg['phi']))
        self.v = float(np.interp(rb,self.bg['r'],self.bg['v'])/(float(self.bg['R'])*100))
        self.nu = float(np.log(d['N'][i]))
        self.T = float(np.exp(d['lnT'][i]))
        self.Pg = float(d['raw'][i,1]-self.arad*self.T**4/3)
        self.Lref = float(d['luminosity'][i,0])
        self.total_baryon = float(d['dm'].sum())
        self.calls = 0
        self.max_calls = 60000
        self.rows = []
        self.maximum_integrations = 30

    def state(self,x,y,L):
        self.calls += 1
        assert self.calls < self.max_calls, 'Registered native call budget'
        r = self.R*(1+y[0]); m = self.m*(1+y[1])
        T = np.exp(y[2]); phi = self.phi+.001*y[3]
        v = self.v+.001*y[4]/self.R; N = np.exp(self.nu+y[5])
        A = np.exp(-2*phi*phi); b = 1-2*m/r
        pg = np.exp(x); pr = self.arad*T**4/3; pt = pg+pr
        assert T > 4000 and b > 0 and r > 0, ('EOS/domain',T,r,b)
        a = self.eos(1,float(np.log(pt)),float(y[2]),self.X)
        rho = a[0]; en = rho*(self.cx*self.c**2+a[2]); eg = en-3*pr
        assert rho > 0 and eg > 0 and abs(a[1]/pt-1)<1e-7
        kap = prior.two.opacity_parts(self.opacity,(np.log(rho),y[2],self.X))[0]
        pgeo = self.G*pt/self.c**4; egeo = self.G*en/self.c**4
        nr = m/(r*r*b)+4*np.pi*r*A**4*pgeo/b+r*v*v/2
        grav = nr-4*phi*v
        F = L/(4*np.pi*r*r*N*N*A**4)
        optical = rho*kap*A/np.sqrt(b)
        force = optical*F/self.c
        pprime = -(eg+pg)*grav+force
        assert pprime < 0, ('Gas pressure inversion/Eddington failure',force/((eg+pg)*grav))
        tprime = -grav-force/(4*pr)
        rp,rt=a[7:9]; up,ut=a[9:11]
        cpT=ut-pt/rho*rt
        nabla_ad=-pt/rho*rt/cpT
        nabla=tprime/(-(en+pt)*grav/pt)
        return dict(r=r,m=m,T=T,phi=phi,v=v,N=N,A=A,b=b,Pgas=pg,Prad=pr,
            Ptotal=pt,rho=rho,energy=en,opacity=kap,F=F,g=grav,nr=nr,
            gas_pressure_prime=pprime,logT_prime=tprime,optical=optical,
            pgeo=pgeo,egeo=egeo,radiative_acceleration_fraction=force/((eg+pg)*grav),
            nabla=nabla,nabla_ad=nabla_ad,entropy=float(a[3]),
            opacity_logR=float(np.log10(rho)-3*np.log10(T)+18))

    def integrate(self,L,cut=1.,tol=2e-5,save=None):
        assert len(self.rows)<self.maximum_integrations, 'Registered integration cap'
        before=self.calls; start=time.monotonic()
        y0=np.zeros(8); y0[2]=np.log(self.T)
        def rhs(x,y):
            z=self.state(x,y,L)
            r,m,phi,v,A,b=[z[k] for k in ['r','m','phi','v','A','b']]
            drdx=z['Pgas']/z['gas_pressure_prime']
            mr=4*np.pi*r*r*A**4*z['egeo']+r*r*b*v*v/2
            vr=4*np.pi/b*(-4*phi*A**4*(z['egeo']-3*z['pgeo'])+r*v*A**4*(z['egeo']-z['pgeo']))-2*(r-m)/(r*r*b)*v
            return drdx*np.array([1/self.R,mr/self.m,z['logT_prime'],v/.001,
                vr*self.R/.001,z['nr'],4*np.pi*r*r*A**3*z['rho']/np.sqrt(b)/self.total_baryon,-z['optical']])
        def floor(x,y):
            r=self.R*(1+y[0]); phi=self.phi+.001*y[3]; N=np.exp(self.nu+y[5])
            F=L/(4*np.pi*r*r*N*N*np.exp(-8*phi*phi))
            return np.log(2*self.sigma*np.exp(4*y[2])/F)-np.log(.75)
        floor.terminal=True; floor.direction=-1
        sol=solve_ivp(rhs,(np.log(self.Pg),np.log(cut)),y0,rtol=tol,
            atol=np.array([1e-10,1e-19,1e-8,1e-15,1e-15,1e-15,1e-20,1e-7]),
            max_step=.25,events=floor,dense_output=True)
        assert sol.success,sol.message
        z=self.state(sol.t[-1],sol.y[:,-1],L)
        mismatch=float(np.log(2*self.sigma*z['T']**4/z['F']))
        if sol.status==1:
            mismatch-=float(sol.t[-1]-np.log(cut))
        row=dict(Linfinity=L,cut_pressure=cut,tolerance=tol,boundary_log_mismatch=mismatch,
            reached_cut=sol.status==0,native_calls=self.calls-before,seconds=time.monotonic()-start,
            endpoint=z,added_baryon_fraction=float(sol.y[6,-1]))
        self.rows.append(row)
        write(OUT/'progress.json',dict(classification='Counterexample candidate',integrations=self.rows,native_calls=self.calls))
        if save:
            assert sol.status==0 and abs(mismatch)<1e-5
            grid=np.linspace(np.log(self.Pg),np.log(cut),401); states=sol.sol(grid)
            values=[self.state(x,y,L) for x,y in zip(grid,states.T)]
            np.savez_compressed(OUT/(save+'.npz'),logPgas=grid,states=states,
                **{key:np.array([a[key] for a in values]) for key in values[0]},Linfinity=L)
            row.update(maximum_radiative_acceleration_fraction=max(a['radiative_acceleration_fraction'] for a in values),
                maximum_nabla_minus_adiabatic=max(a['nabla']-a['nabla_ad'] for a in values),
                opacity_logR_range=[min(a['opacity_logR'] for a in values),max(a['opacity_logR'] for a in values)],
                tau_base=float(-states[7,-1]),
                base_flux_mismatch_from_original=L/self.Lref-1)
            write(OUT/(save+'.json'),dict(classification='Counterexample candidate',**row))
        return row

    def shoot(self,cut,tol,label,previous=None):
        history=[]
        cached={}
        if previous is None:
            pilot=json.loads((OUT/'pilot.json').read_text())['row']
            cached[0.]=pilot['boundary_log_mismatch']
        def objective(logratio):
            if logratio in cached:return cached[logratio]
            assert len(history)<24, 'Registered shooting evaluation cap'
            row=self.integrate(self.Lref*np.exp(logratio),cut,tol)
            history.append([logratio,row['boundary_log_mismatch']])
            return row['boundary_log_mismatch']
        if previous is None:
            lo,hi=(0.,np.log(2)) if cached[0.]>0 else (np.log(.25),0.)
        else:
            center=np.log(previous['Linfinity']/self.Lref)
            value=objective(center)
            lo,hi=(center,center+.0021) if value>0 else (center-.0021,center)
        fit=brentq(objective,lo,hi,xtol=2e-10,rtol=2e-10,maxiter=22)
        row=self.integrate(self.Lref*np.exp(fit),cut,tol,save=label)
        row['shooting_evaluations']=len(history)
        return row


def prepare():
    assert not OUT.exists();OUT.mkdir()
    paths=[Path(__file__),Path(prior.__file__),prior.OUT/'coefficients.npz',
        prior.surface.OUT/'background.npz',prior.OUT.parent/'def-photon-matter-split/matter.npz',
        ROOT/'verification/opacity_tables.py',ROOT/'verification/direct_ion_eos.py']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='5e753c069',
        claim='Construct a native EOS, state-dependent opacity envelope with outward luminosity and simultaneous hydrostatic/radiative constraints, joined at native cell132. Test the old nonstationary outer temperature/flux against the new matched boundary.',
        model='Quasistatic DEF hydrostatic constraints plus grey Eddington E=3P=aT^4 and constant redshifted luminosity. Separate gas radiation force and Tolman-corrected thermal gradient. At low gas pressure impose E=2F/c (no incoming half-isotropic Eddington boundary). No Rosseland-to-Planck substitution.',
        physical_limits='Grey angular closure, LTE and constant-luminosity envelope are explicit assumptions. Nonzero flux forbids an exactly static metric: record required mass loss and heat-stress scale. No stationary spacetime, exact spectral atmosphere, fixed-inventory whole-star replacement, convection closure or new observable claim.',
        junction='Fix original r,m,phi,phi prime,N,total P,T and uniform outer composition at cell132. Solve only L to meet the exterior temperature boundary. Compare resulting L with the old diffusive inner L; do not silently declare those matched. Added/replaced exterior baryons are recorded, not forced equal.',
        paths=['Pgas cutoff1dyn/cm2 tolerance2e-5','same cutoff tolerance2e-6','cutoff0.25dyn/cm2 tolerance2e-6'],
        gates=dict(boundary_log=1e-5,luminosity_tolerance_relative=.002,temperature_tolerance_relative=.001,
            luminosity_cutoff_relative=.002,positive_outward_flux=True),
        budget=dict(pilot_seconds=120,total_seconds=600,CPU_threads=1,memory_GB=4,
            native_calls=60000,maximum_shoot_evaluations=24,automatic_expansion=False),
        forecast='EOS finite controls in Phase89 about1.21s for35 calls; new envelope step rate unmeasured. Measure one full native integration before shooting.',
        references=['https://docs.mesastar.org/en/24.03.1/atm/t-tau.html'],
        bindings={str(p):digest(p) for p in paths}))
    write(OUT/'symbolic.json',symbolic())
    resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)));signal.alarm(120)
    start=time.monotonic(); p=Envelope(); row=p.integrate(p.Lref)
    elapsed=time.monotonic()-start
    forecast=elapsed+1.5*row['seconds']*60+30
    write(OUT/'pilot.json',dict(classification='Counterexample candidate',row=row,seconds=elapsed,
        forecast_seconds=forecast,forecast_assumption='60 full integrations plus50percent margin and30s; actual shoot counts and refined-cost increase unmeasured.',
        base_cell=p.i,base_r_cm=p.R,base_T=p.T,base_gas_pressure=p.Pg,
        original_Linfinity=p.Lref,original_outer_baryon_fraction=float(p.d['dm'][:p.i+1].sum()/p.total_baryon)))
    print('ENVELOPE PILOT',elapsed,forecast,row['boundary_log_mismatch'],flush=True)


def run():
    plan=json.loads((OUT/'plan.json').read_text());pilot=json.loads((OUT/'pilot.json').read_text())
    budget=json.loads((OUT/'reuse-budget.json').read_text())
    assert not (OUT/'result.json').exists() and budget['forecast_seconds']<600
    for name,value in plan['bindings'].items(): assert digest(name)==value,name
    signal.alarm(int(600-pilot['seconds']));resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)))
    start=time.monotonic();p=Envelope()
    rows=[]
    for cut,tol,label in [(1.,2e-5,'coarse'),(1.,2e-6,'fine'),(.25,2e-6,'thin')]:
        rows.append(p.shoot(cut,tol,label,rows[-1] if rows else None))
    a,b,c=rows; tolerance=abs(a['Linfinity']/b['Linfinity']-1);cutoff=abs(c['Linfinity']/b['Linfinity']-1)
    coarse=np.load(OUT/'coarse.npz');fine=np.load(OUT/'fine.npz')
    temperature=float(max(abs(coarse['T']/fine['T']-1)))
    period=json.loads((OUT.parent/'def-orbital-charge-fem/pilot.json').read_text())['row']['omega_R_over_c']
    # The original orbital frequency uses its own free-surface radius.
    orbit=json.loads((OUT.parent/'def-orbital-charge-fem/pilot.json').read_text())
    R=float(np.load(prior.surface.OUT/'background.npz')['r'][-1]*np.load(prior.surface.OUT/'background.npz')['R'])*100
    orbit_seconds=2*np.pi*R/(p.c*period)
    mass_geom=float(orbit['ADM_geom_m'])*100
    mass_gram=mass_geom*p.c*p.c/p.G
    write(OUT/'result.json',dict(classification='Counterexample candidate',
        passed=tolerance<.002 and cutoff<.002 and temperature<.001,
        rows=rows,luminosity_tolerance_relative=tolerance,temperature_tolerance_relative=temperature,
        luminosity_cutoff_relative=cutoff,required_Linfinity=b['Linfinity'],
        required_minus_original_L_relative=b['Linfinity']/p.Lref-1,
        necessary_mass_loss_gram_s=b['Linfinity']/p.c**2,
        original_orbit_seconds=orbit_seconds,energy_over_original_ADM_per_orbit=b['Linfinity']*orbit_seconds/(mass_gram*p.c**2),
        native_calls=p.calls,seconds=time.monotonic()-start,total_compute_seconds=pilot['seconds']+time.monotonic()-start,
        quasi_static_grey_envelope_solved=True,full_time_radial_Einstein_equation_solved=False,
        unchanged_interior_luminosity_matched=False,same_inventory_whole_star=False,
        spectral_atmosphere=False,full_goal_complete=False))
    print('ENVELOPE RESULT',(OUT/'result.json').read_text(),flush=True)


def reuse_budget():
    assert not (OUT/'result.json').exists() and not (OUT/'reuse-budget.json').exists()
    pilot=json.loads((OUT/'pilot.json').read_text());plan=json.loads((OUT/'plan.json').read_text())
    estimate=pilot['seconds']+1.5*30*pilot['row']['seconds']+60
    assert estimate<600
    plan['bindings'][str(Path(__file__))]=digest(Path(__file__))
    plan['resource_revision']='Preserve initial940s forecast and pilot. Reuse pilot objective and previous shooting roots; cap production at30 integrations including final profiles. Original cases/gates and600s wall cap unchanged.'
    plan['budget']['maximum_production_integrations']=30
    write(OUT/'plan.json',plan)
    write(OUT/'reuse-budget.json',dict(classification='Counterexample candidate',forecast_seconds=estimate,
        initial_forecast_rejected=True,reused_pilot=True,reused_shooting_roots=True,
        assumption='30 integrations at measured pilot speed plus50percent and60s for final401-point native readouts; refined/cutoff speeds unmeasured. No automatic expansion.',
        unchanged_physical_model_and_gates=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','reuse_budget','run'])
    globals()[parser.parse_args().action]()
