"""Finite native H reactions in the existing conservative spherical flow."""
from pathlib import Path
from types import SimpleNamespace
import inspect
import json
import signal
import sys
import textwrap
import time
import numpy as np
import sympy as s
import def_native_hydrogen_exchange as exchange

old=exchange.old
OUT=exchange.OUT/'flow'
C=old.C
write=old.write


def replace(source,before,after):
    assert source.count(before)==1,(before,source.count(before))
    return source.replace(before,after)


primitive_source=textwrap.dedent(inspect.getsource(old.previous.old.Flow.primitive))
primitive_source=replace(primitive_source,'m=self.base;D=', 'm=self.base;D=') # Bind the inherited implementation explicitly.
primitive_source=replace(primitive_source,'sigma=self.seed.copy();p=',
    'all_y=np.divide(U[3],D,out=np.full_like(D,self.eos.y0),where=active);self.eos.y=all_y\n    sigma=self.seed.copy();p=')
primitive_source=replace(primitive_source,'def residual(theta):','def residual(theta):\n                self.eos.y=np.array([all_y[i]])')
primitive_source=replace(primitive_source,'for _ in range(4):\n            v=', 'self.eos.y=all_y\n        for _ in range(4):\n            v=')
primitive_source=replace(primitive_source,'return np.array([rho,v,sigma])','self.eos.y=all_y\n    return np.array([rho,v,sigma,all_y])')
ns=dict(vars(old.previous.old),OUT=OUT,temperature=SimpleNamespace(parent=SimpleNamespace(EOS=exchange.EOS)))
exec(compile(primitive_source,__file__,'exec'),ns)

rhs_source=textwrap.dedent(old.previous.rhs_source)
rhs_source=replace(rhs_source,'rho,v,sigma=V','rho,v,sigma,y=V')
rhs_source=replace(rhs_source,'self.conserved(rho,v,sigma,m.a)','self.conserved(rho,v,sigma,y,m.a)')
rhs_source=replace(rhs_source,'rate=-C*np.diff(flux)/m.vol', '''# A passive species must follow the same mass flux. Reconstructing rho
    # and y independently does not preserve positivity of their product.
    # ponytail: first-order donor y; retain only if the registered grid comparison passes.
    donor=np.where(flux[0]>=0,np.r_[self.eos.y0,y],np.r_[y,self.eos.y0])
    flux[3]=flux[0]*donor
    rate=-C*np.diff(flux)/m.vol''')
rhs_source=replace(rhs_source,'return rate,ledger,dt', '''self.eos.y=y
    reaction=self.eos.reactions(rho,sigma)
    dilution=(1-mu)/2
    absorbed=reaction[:,0]*dilution[:,None]
    emitted=reaction[:,1]+reaction[:,2]*dilution[:,None]
    net=absorbed-emitted
    species=-m.a*rho*net[:,0]
    Q=rho*self.eos.nH*net[:,1]/C**2
    # Declared bath approximation: evaluate photons in the local static
    # frame, use its hemisphere mean for absorption; isotropic spontaneous
    # emission. A velocity-frame error estimate is recorded, not hidden.
    force_energy=rho*self.eos.nH*(absorbed[:,1]-dilution*reaction[:,2,1])*(1+mu)/2/C**2
    momentum=m.a*W*(force_energy+v*Q)
    energy=m.a*m.a*W*(Q+v*force_energy)
    rate[1]+=momentum;rate[2]+=energy;rate[3]+=species
    transfer=float(np.sum(energy*m.vol));ledger[1]+=transfer;ledger[2]+=transfer
    ledger=np.r_[ledger,C*(flux[3,0]-flux[3,-1]),np.sum(species*m.vol)]
    stiffness=dilution*reaction[:,0,0]/y+(reaction[:,1,0]+dilution*reaction[:,2,0])/(1-y)
    dt=min(dt,.2/max(float(max(stiffness)),1e-100))
    for amount,change in [(U[3],rate[3]),(U[0]-U[3],rate[0]-rate[3])]:
        losses=(change<0)&(U[0]>=self.eos.floor)
        if np.any(losses):dt=min(dt,.8*float(min(amount[losses]/(-change[losses]))))
    optical=float(np.sum(m.B*m.dx*rho*self.eos.nH*absorbed[:,1]/(C**3*Fbase)))
    self.max_absorption=max(self.max_absorption,optical)
    self.max_speed=max(self.max_speed,float(max(abs(v))))
    self.min_y=min(self.min_y,float(min(y[rho>=self.eos.floor])))
    return rate,ledger,dt''')
exec(compile(rhs_source,__file__,'exec'),ns)

run_source=textwrap.dedent(inspect.getsource(old.previous.old.Flow.run))
run_source=replace(run_source,'ledger=np.zeros(3);discard=np.zeros(3)','ledger=np.zeros(5);discard=np.zeros(4)')
run_source=replace(run_source,'k2,l2,_=self.rhs(trial,t+dt);nxt=(U+trial+dt*k2)/2', '''for attempt in range(12):
            k2,l2,second_dt=self.rhs(trial,t+dt)
            if dt<=second_dt*(1+1e-12):break
            dt=min(dt/2,second_dt);trial=U+dt*k
        else:raise AssertionError('SSP stage time step')
        nxt=(U+trial+dt*k2)/2''')
run_source=replace(run_source,'rho,v,sigma=self.primitive(U);','rho,v,sigma,y=self.primitive(U);')
run_source=replace(run_source,'ip,iu,*_=self.eos(ir,np.log(iT))','self.eos.y=self.eos.y0;ip,iu,*_=self.eos(ir,np.log(iT));self.eos.y=y')
run_source=replace(run_source,'passed=baryon<1e-10', '''species_residual=np.sum((U[3]-self.initial[3])*m.vol)+discard[3]-ledger[3]-ledger[4]
    species_error=abs(species_residual)/max(np.sum(self.initial[3]*m.vol),1e-100)
    passed=species_error<1e-9 and self.max_absorption<.001 and baryon<1e-10''')
run_source=replace(run_source,'scattering_work_into_gas_erg=', 'photon_energy_into_gas_erg=')
run_source=replace(run_source,'full_GR_scalar_feedback=False,', '''neutral_fraction_minimum=self.min_y,maximum_speed_over_c=self.max_speed,
        species_ledger_relative=float(species_error),maximum_absorption_energy_optical_depth=self.max_absorption,
        finite_H_reactions=True,photon_Killing_energy_paired=True,
        static_frame_bath_approximation=True,full_GR_scalar_feedback=False,''')
exec(compile(run_source,__file__,'exec'),ns)


class Flow(old.DiluteFlow):
    primitive=ns['primitive'];rhs=ns['rhs'];run=ns['run']

    def __init__(self,n):
        baseline=old.DiluteFlow(n);self.__dict__.update(baseline.__dict__);m=self.base;rho=self.initial[0].copy();self.eos=exchange.EOS()
        self.initial=self.conserved(rho,np.zeros(n),self.seed,np.full(n,self.eos.y0),m.a)[0]
        self.max_absorption=0.;self.max_speed=0.;self.min_y=self.eos.y0

    def conserved(self,rho,v,lt,y,a):
        self.eos.y=y
        U,F,thermo=old.previous.old.Flow.conserved(self,rho,v,lt,a)
        return np.vstack([U,U[0]*y]),np.vstack([F,F[0]*y]),thermo

    def reconstruct(self,V,t):
        left,right=old.previous.Flow.reconstruct(self,V[:3],t)
        y=V[3];ext=np.r_[self.eos.y0,y,self.eos.y0];slope=np.zeros_like(ext)
        slope[1:-1]=old.previous.old.bank.task.minmod(ext[1:-1]-ext[:-2],ext[2:]-ext[1:-1])
        return np.vstack([left,ext[:-1]+slope[:-1]/2]),np.vstack([right,ext[1:]-slope[1:]/2])


def symbolic():
    A,B,lam,u,n=s.symbols('A B lambda u n',positive=True)
    balance=A*(n-s.exp(lam-u)*(1+n))
    assert s.simplify(balance.subs({lam:0,n:1/(s.exp(u)-1)}))==0
    a,W,v,Q,F=s.symbols('a W v Q F',real=True)
    gas=a*a*W*(Q+v*F);photon=-gas
    assert s.expand(gas+photon)==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        scope='Pointwise reciprocal thermodynamic closure balances at zero H affinity in a same-temperature Planck bath. The declared Lorentz-source Killing-energy transfer cancels its paired photon ledger exactly. No assertion of a microscopic inverse or evolved photon transport.'))


def run():
    assert not (OUT/'result.json').exists();assert json.loads((exchange.OUT/'repaired-bank.json').read_text())['passed']
    assert (OUT/'first-flow-source.py').exists()
    (OUT/'pilot-224-progress.npz').rename(OUT/'first-pilot-progress.npz')
    for name in ['primitive','rhs','run']:(OUT/f'expanded-{name}.py').rename(OUT/f'first-expanded-{name}.py')
    write(OUT/'species-transport-repair.json',dict(classification='Counterexample candidate',
        failure='Independent primitive rho and neutral-fraction reconstruction produced negative neutral H in the first pilot, after the saved0.548163ms state. Preserve the original source and checkpoint.',
        repair='Transport D*y with the exact same mass flux and an upwind cell neutral fraction. Limit the combined transport/reaction Euler decrement of each species to80percent and enforce the admissible second-stage step in SSP RK2. No abundance clipping or change to reaction rates, EOS, physical domain, horizon or gates.',
        approximation='The added passive-species advection is first order. Judge the reactive output and abundance through the same448/896 comparison.',
        reserved_failed_seconds=10,remaining_flow_seconds=170,source_sha256=old.cold.sha(__file__)))
    start=time.monotonic();signal.alarm(170);symbolic()
    for name,src in [('primitive',primitive_source),('rhs',rhs_source),('run',run_source)]: (OUT/f'expanded-{name}.py').write_text(src)
    write(OUT/'binding.json',dict(classification='Counterexample candidate',source_sha256=old.cold.sha(__file__),exchange_sha256=old.cold.sha(exchange.__file__),bank_sha256=old.cold.sha(exchange.OUT/'repaired-bank.npz'),seconds=180,
        source_boundary='The reaction changes D*y, fluid momentum and total gas energy. The opposite photon Killing-energy transfer is accumulated. The optically thin Planck bath is prescribed; no frequency-dependent photon transport, omitted chemical channels or complete physical charge is claimed.'))
    pilot=Flow(224).run('pilot-224');forecast=20*pilot['seconds']*1.4
    write(OUT/'measured-budget.json',dict(pilot_seconds=pilot['seconds'],forecast_remaining_seconds=forecast,remaining_seconds=180-(time.monotonic()-start),assumption='CFL cell-squared scaling plus40percent, unresolved fine-grid inversion cost.'))
    assert forecast<180-(time.monotonic()-start),'Measured evolution budget'
    rows=[Flow(n).run(f'cells-{n}') for n in [448,896]]
    errors={k:abs(rows[0][k]/rows[1][k]-1) for k in ['gas_outside_original_radius_g','integrated_trace_energy_erg']}
    frozen=json.loads((old.OUT/'dilute/cells-896.json').read_text());change={k:rows[-1][k]/frozen[k]-1 for k in errors}
    result=dict(classification='Counterexample candidate',passed=bool(max(errors.values())<.02),refinement=errors,relative_change_from_frozen=change,seconds=time.monotonic()-start,finite_H_reactions=True,physical_chemistry_closed=False,full_GR_scalar_feedback=False,final_charge_solved=False,full_goal_complete=False)
    write(OUT/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True);assert result['passed']


def production():
    assert not (OUT/'result.json').exists();assert json.loads((OUT/'pilot-224.json').read_text())['passed']
    start=time.monotonic();signal.alarm(40)
    write(OUT/'budget-reassessment.json',dict(classification='Counterexample candidate',
        stopped='The completed224-cell reactive pilot took8.09s and562 steps. Naive cell-squared scaling forecasts226s and correctly stopped the production dispatch.',
        decision='The coarse pilot was limited by the measured reaction time step, not only the hydrodynamic CFL. Reuse it without rerunning. Measure only the already registered448 path under40s, then permit896 only if its measured cell-squared forecast fits the original remaining155s. No additional path or increased aggregate budget.',
        reserved_previous_flow_seconds=25,remaining_seconds=155,coarse_seconds=40,fine_forecast='5.6 times measured448 time plus3s setup',source_sha256=old.cold.sha(__file__)))
    coarse=Flow(448).run('cells-448');elapsed=time.monotonic()-start;forecast=5.6*coarse['seconds']+3
    write(OUT/'production-budget.json',dict(coarse_seconds=coarse['seconds'],forecast_fine_seconds=forecast,remaining_seconds=155-elapsed))
    assert forecast<155-elapsed,'Measured fine-grid remaining budget'
    signal.alarm(max(1,int(155-elapsed)));fine=Flow(896).run('cells-896')
    errors={k:abs(coarse[k]/fine[k]-1) for k in ['gas_outside_original_radius_g','integrated_trace_energy_erg']}
    frozen=json.loads((old.OUT/'dilute/cells-896.json').read_text());change={k:fine[k]/frozen[k]-1 for k in errors}
    result=dict(classification='Counterexample candidate',passed=bool(max(errors.values())<.02),refinement=errors,relative_change_from_frozen=change,seconds=time.monotonic()-start,finite_H_reactions=True,physical_chemistry_closed=False,full_GR_scalar_feedback=False,final_charge_solved=False,full_goal_complete=False)
    write(OUT/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True);assert result['passed']


if __name__=='__main__':globals()[sys.argv[1]]()
