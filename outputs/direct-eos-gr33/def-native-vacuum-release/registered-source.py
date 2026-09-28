"""Native gas-EOS rarefaction at the finite-pressure stellar cut.

Counterexample candidate: local planar, adiabatic chemical-equilibrium gas
release with independent photons. The unqueried dilute tail is enclosed under
an explicit sound-speed envelope; no cold EOS continuation is fabricated.
"""
from pathlib import Path
import argparse
import json
import signal
import time
import numpy as np
from scipy.optimize import brentq
from scipy.interpolate import PchipInterpolator
from scipy.integrate import cumulative_simpson, simpson
import def_native_metric_transport as prior
import gr_radiation_eos_split as split

old=prior.old
OUT=prior.OUT.parent/'def-native-vacuum-release'
write=old.write


def symbolic():
    import sympy as s
    D,e,p,v=s.symbols('D e p v',real=True);W=1/s.sqrt(1-v*v)
    energy=(e+p)*W*W-p;radial=(e+p)*W*W*v*v+p
    assert s.simplify(energy-radial-2*p-(e-3*p))==0
    a0,a,g=s.symbols('a0 a g',positive=True)
    invariant=2/s.sqrt(g-1)*s.atanh(a/s.sqrt(g-1))
    # Along a gamma-law isentrope da/dlnrho=a*(g-1-a^2)/2.
    assert s.simplify(s.diff(invariant,a)*a*(g-1-a*a)/2-a)==0
    return dict(classification='Proven',passed=True,
        native='At fixed gas entropy, dlnT/dlnrho=(P/rho-u_lnrho)/cvT, cs^2=Gamma1*P/(e+P). The left rarefaction has rapidity Y=integral(cs*dlnrho) from the current density to the initial density and xi=(tanh(Y)-cs)/(1-tanh(Y)*cs).',
        trace='E_lab-P_radial-2*P_tangential=e_gas-3*P_gas. LTE photons must be subtracted before using the gas entropy and sound speed.',
        boundary='A stress-free gas-vacuum material interface requires Pgas=0. Releasing Pgas>0 produces a nonlinear fan, not a zero-amplitude perturbation of a supported surface.',
        tail='If cs(rho)<=cs_cut*(rho/rho_cut)^(delta/2) for0<rho<=rho_cut, the remaining rapidity is<=2*cs_cut/delta. This is a declared tail premise, not a native EOS proof.',
        scope='Relativistic isentropic simple-wave identities. Planar geometry, equilibrium gas chemistry and negligible external gas heating are separate applicability conditions.')


class Fan:
    def __init__(self):
        self.env=dict(np.load(old.prior.OUT/'final-envelope.npz'))
        self.core=dict(np.load(old.prior.OUT/'final-core.npz'))
        self.eos=split.GasEOS();self.total=old.prior.envelope.Envelope()
        self.rho=float(self.env['rho'][-1]);self.T=float(self.env['T'][-1]);self.X=self.core['X'][0]
        self.cx=self.total.cx;self.c=old.C;self.calls=0
        self.base=self.call(np.log(self.rho),np.log(self.T));self.entropy=float(self.base[3])
        assert abs(self.base[1]/self.env['Pgas'][-1]-1)<1e-7

    def call(self,lr,lt):
        self.calls+=1;assert self.calls<=1600,'Native call budget'
        row=self.eos(2,float(lr),float(lt),self.X)
        assert np.all(np.isfinite(row)) and row[1]>0 and row[10]>0 and row[4]>1,('EOS domain/stability',lr,lt,row[:11])
        return row

    def table(self,step,label):
        start=time.monotonic();density=np.linspace(0,-18,int(round(18/step))+1)
        rows=[];temperatures=[];last=np.log(self.T);guess=last
        for index,x in enumerate(density):
            if index==0:a=self.base;lt=last
            else:
                lastrow=rows[-1];ad=(lastrow[1]/lastrow[0]-lastrow[9])/lastrow[10]
                guess=last+ad*(x-density[index-1]);cache={}
                def error(lt):
                    if lt not in cache:cache[lt]=self.call(np.log(self.rho)+x,lt)
                    return (cache[lt][3]-self.entropy)*np.exp(lt)/cache[lt][10]
                lo,hi=guess-.15,last+.01
                elo,ehi=error(lo),error(hi)
                assert elo*ehi<0,('Entropy bracket failed',x,elo,ehi)
                lt=brentq(error,lo,hi,xtol=2e-12,rtol=2e-13);a=cache.get(lt)
                if a is None:a=self.call(np.log(self.rho)+x,lt)
            rows.append(a);temperatures.append(np.exp(lt));last=lt
            # Stop at the declared domain floor instead of extending a failed
            # cold EOS. The unused density interval is never filled in.
            if temperatures[-1]<250:break
        rows=np.array(rows);T=np.array(temperatures);density=density[:len(rows)]
        enthalpy=self.cx*self.c**2+rows[:,2]+rows[:,1]/rows[:,0]
        energy=rows[:,0]*(self.cx*self.c**2+rows[:,2])
        ad=(rows[:,1]/rows[:,0]-rows[:,9])/rows[:,10]
        gamma=rows[:,5]+rows[:,6]*ad
        chain=np.max(abs(gamma/rows[:,4]-1));assert chain<1e-4,chain
        cs=np.sqrt(gamma*rows[:,1]/(energy+rows[:,1]));assert max(cs)<1
        rapidity=cumulative_simpson(cs,x=-density,initial=0)
        velocity=np.tanh(rapidity);xi=(velocity-cs)/(1-velocity*cs)
        assert np.all(np.diff(xi)>0),('Non-convex EOS/simple-wave failure',np.min(np.diff(xi)))
        entropy_error=max(abs((rows[:,3]-self.entropy)*T/rows[:,10]));assert entropy_error<1e-8,entropy_error
        np.savez_compressed(OUT/f'{label}.npz',log_density_ratio=density,T=T,raw=rows,energy=energy,enthalpy=enthalpy,cs=cs,rapidity=rapidity,velocity_over_c=velocity,xi=xi)
        data=dict(classification='Counterexample candidate',native_states=len(rows),native_calls=self.calls,seconds=time.monotonic()-start,
            native_gamma_chain_relative=float(chain),native_entropy_scaled_error=float(entropy_error),minimum_T=float(T[-1]),minimum_density_ratio=float(np.exp(density[-1])),
            minimum_gamma1=float(min(gamma)),maximum_gamma1=float(max(gamma)),head_sound_speed_m_s=float(cs[0]*self.c/100),resolved_maximum_speed_m_s=float(velocity[-1]*self.c/100))
        write(OUT/f'{label}.json',data);return data


def integrate_fan(d,order=64):
    # Integrate all resolved stresses in characteristic speed, not a uniform
    # radial grid that could miss a rarefaction narrower than one native cell.
    x=d['xi'];z,w=np.polynomial.legendre.leggauss(order)
    boundaries=np.sort(np.unique(np.r_[x,0.]))
    boundaries=boundaries[(boundaries>=x[0])&(boundaries<=x[-1])]
    at=(boundaries[:-1,None]+np.diff(boundaries)[:,None]*(z+1)/2).ravel();weights=(np.diff(boundaries)[:,None]*w/2).ravel()
    lr=PchipInterpolator(x,d['log_density_ratio'])(at)
    pressure=np.exp(PchipInterpolator(x,np.log(d['raw'][:,1]))(at))
    specific=PchipInterpolator(x,d['raw'][:,2])(at)
    speed=PchipInterpolator(x,d['velocity_over_c'])(at)
    rho=d['raw'][0,0]*np.exp(lr);W=1/np.sqrt(1-speed**2);Wminus=speed**2/(np.sqrt(1-speed**2)*(1+np.sqrt(1-speed**2)))
    hrest=float(d['enthalpy'][0]-d['raw'][0,2]-d['raw'][0,1]/d['raw'][0,0])
    energy=rho*(hrest+specific);left=at<0
    # Subtract the conserved baryon rest term before measuring the small
    # non-rest energy and trace ledgers.
    baryon=rho*W-rho[0]*0-d['raw'][0,0]*left
    nonrest=rho*hrest*W*Wminus+(rho*specific+pressure)*W*W-pressure-d['raw'][0,0]*d['raw'][0,2]*left
    trace_minus_rest=-rho*hrest*Wminus+rho*specific-3*pressure-(d['raw'][0,0]*d['raw'][0,2]-3*d['raw'][0,1])*left
    mom=(energy+pressure)*W*W*speed
    return dict(baryon=float(weights@baryon),nonrest_energy=float(weights@nonrest),trace_minus_baryon_rest=float(weights@trace_minus_rest),
        momentum=float(weights@mom),pressure_change=float(weights@(pressure-d['raw'][0,1]*left)),
        matter_outside=float(weights@(rho*W*(at>=0))),positive_trace_envelope=float(weights@(abs(baryon)*hrest+abs(trace_minus_rest))))


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='75ee1384d',
        claim='Replace the assumption of a supported almost-stationary gas cut by its actual native gas-EOS nonlinear vacuum rarefaction, and quantify the resulting matter/stress input for the GR connection.',
        preserved='Current native surface density, temperature, composition, chemical energy convention and same EOS with its compiled LTE radiation removed. No refit polytrope or frozen-ion substitute.',
        physical_conditions='Local planar Riemann release into gas vacuum, gas entropy conserved, equilibrium native chemistry, independent photons. Curvature, gravity, heat/ionization kinetics and full stellar feedback are not asserted negligible without separate assessment.',
        table=dict(log_density_interval=[0,-18],coarse_step=.25,fine_step=.125,minimum_temperature_K=250,native_call_cap=1600),
        tail='No native extrapolation: assume cs(rho)<=cs_cut*(rho/rho_cut)^0.05 below the last native state. This yields an explicit residual rapidity/mass envelope; the premise remains conditional.',
        gates=dict(native_entropy=1e-8,native_gamma_chain=1e-4,characteristic_speed_refinement=.005,conserved_baryon_relative=.002,conserved_nonrest_energy_relative=.005),
        budget=dict(total_seconds=180,CPU_threads=1,new_whole_star_integrations=0),
        decision='One coarse and one fine native isentrope only, with measured pilot budget. Stop on EOS failure/nonconvex fan or budget; no density/temperature extension, extra grids or relaxed gates.',
        reference='https://doi.org/10.12942/lrr-1999-3',reference_scope='Native isentropic Riemann-invariant construction follows the relativistic simple-wave method; this source does not validate our EOS or stellar assumptions.',
        symbolic=symbolic(),bindings={str(p.relative_to(old.ROOT)):old.photons.digest(p) for p in [Path(__file__),Path(split.__file__),old.prior.OUT/'final-envelope.npz',old.prior.OUT/'final-core.npz',prior.OUT/'p4-64.npz']}))


def run():
    assert not (OUT/'result.json').exists();plan=json.loads((OUT/'plan.json').read_text())
    for p,h in plan['bindings'].items():assert old.photons.digest(old.ROOT/p)==h,p
    signal.alarm(170);start=time.monotonic();fan=Fan()
    first=fan.table(.25,'coarse');forecast=time.monotonic()-start+2.5*(first['seconds']+3)
    write(OUT/'pilot-budget.json',dict(classification='Counterexample candidate',forecast_seconds=forecast,native_calls=fan.calls,
        assumption='Measured coarse native table, twice as many density intervals and25percent margin for the fine table, plus bookkeeping.'))
    assert forecast<170,'Measured native table budget exceeded'
    second=fan.table(.125,'fine');a=np.load(OUT/'coarse.npz');b=np.load(OUT/'fine.npz')
    speed_error=float(max(abs(np.interp(a['log_density_ratio'][::-1],b['log_density_ratio'][::-1],b['velocity_over_c'][::-1])-a['velocity_over_c'][::-1]))/max(b['velocity_over_c']))
    moments=integrate_fan(b);check=integrate_fan(b,96)
    for k in moments:assert abs(moments[k]-check[k])<1e-8*max(abs(moments[k]),1e-100),(k,moments[k],check[k])
    rho0=fan.rho;energy_scale=rho0*max(abs(b['raw'][:,2]-b['raw'][0,2]))+fan.base[1]
    span=b['xi'][-1]-b['xi'][0];baryon_error=abs(moments['baryon'])/(rho0*span)
    nonrest_error=abs(moments['nonrest_energy'])/(energy_scale*span)
    tail_rapidity=20*float(b['cs'][-1]);vacuum_upper=np.tanh(float(b['rapidity'][-1])+tail_rapidity)
    tail_span=vacuum_upper-float(b['xi'][-1]);tail_mass=float(b['raw'][-1,0])*tail_span/np.sqrt(1-vacuum_upper**2)
    env=fan.env;R=float(env['r'][-1]);A=float(env['A'][-1]);N=float(env['N'][-1]);bs=float(env['b'][-1])
    horizon=json.loads((old.OUT/'plan.json').read_text())['horizon_seconds'];proper_time=A*N*horizon
    length_scale=fan.c*proper_time
    head=length_scale*float(b['xi'][0]);tail_interval=length_scale*np.array([float(b['velocity_over_c'][-1]),vacuum_upper])
    area=4*np.pi*(A*R)**2;outside=area*length_scale*moments['matter_outside'];tail_mass_bound=area*length_scale*tail_mass
    gravity=fan.c**2*float(env['g'][-1])*np.sqrt(bs)/A
    params=dict(gravity_impulse_over_head_sound_speed=gravity*proper_time/(b['cs'][0]*fan.c),
        maximum_fan_span_over_radius=float((tail_interval[1]-head)/(A*R)),
        initial_density_scale_height_fraction=float(abs(head)*abs(env['gas_pressure_prime'][-1])/env['Pgas'][-1]),
        initial_Rosseland_optical_depth_scale=float(env['rho'][-1]*env['opacity'][-1]*(tail_interval[1]-head)),
        optical_scope='A Rosseland value is not an absorption/heating or kinetic-rate bound.')
    passed=speed_error<.005 and baryon_error<.002 and nonrest_error<.005
    write(OUT/'result.json',dict(classification='Counterexample candidate',native_local_fan_numerics_passed=bool(passed),
        seconds=time.monotonic()-start,native_calls=fan.calls,speed_refinement_relative=speed_error,
        resolved_baryon_conservation_relative=baryon_error,resolved_nonrest_energy_conservation_relative=nonrest_error,
        moments=moments,conditional_tail_rapidity_upper=tail_rapidity,conditional_tail_baryon_g_upper=tail_mass_bound,
        coordinate_horizon_seconds=horizon,local_proper_time_seconds=proper_time,
        inward_head_distance_m=float(head/100),vacuum_front_distance_m_interval=(tail_interval/100).tolist(),
        resolved_maximum_velocity_m_s=float(b['velocity_over_c'][-1]*fan.c/100),gas_outside_original_cut_g=outside,
        local_approximation_parameters=params,
        equilibrium_chemistry_during_release_certified=False,unqueried_native_cold_tail_certified=False,
        full_spherical_GR_fan_coupled=False,final_charge_solved=False,full_goal_complete=False))
    signal.alarm(0)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','run']);globals()[parser.parse_args().action]()
