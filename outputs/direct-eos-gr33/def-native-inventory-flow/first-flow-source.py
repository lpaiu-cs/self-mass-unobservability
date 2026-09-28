"""Actual conservative release with native fixed chemical inventories."""
from pathlib import Path
from types import FunctionType, SimpleNamespace
import json
import inspect
import re
import signal
import sys
import time
import textwrap
import numpy as np
from scipy.interpolate import RectBivariateSpline
import def_native_cold_population as cold
import def_native_metric_release as previous

OUT=cold.OUT.parent/'def-native-inventory-flow'
C=previous.C
write=cold.write

# Freeze populated stages as before, but bound initially negligible stages
# when they become populated away from the original adiabat. Their initial
# underflow zeros are never inverted into a finite target population.
constrain_source=textwrap.dedent(inspect.getsource(cold.old.Ions.constrain))
constrain_source=constrain_source.replace(
    "ratios=np.log(np.maximum(target[e,:z+1],1e-290)/np.maximum(current[e,:z+1],1e-290))",
    "desired=np.where(active[e,:z+1],target[e,:z+1],target[e].sum()*1e-18)\n            ratios=np.log(np.maximum(desired,1e-290)/np.maximum(current[e,:z+1],1e-290))")
constrain_source=constrain_source.replace('take=active[e,1:z+1]',
    'take=active[e,1:z+1]|(current[e,1:z+1]>target[e].sum()*1e-18)')
constrain_namespace=dict(vars(cold.old));exec(compile(constrain_source,__file__,'exec'),constrain_namespace)


class InventoryIons(cold.StableIons):
    constrain=constrain_namespace['constrain']


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='d4c274345',
        claim='Apply the repaired same-native fixed-inventory constitutive law to the actual conservative spherical gas release with paired photon scattering work.',
        decision='Compare actual mass/trace histories against the retained LTE producer, before assigning the old LTE wave a physical interpretation. This frozen-reaction control does not complete finite chemistry or absorptive photon closure.',
        model='Uniform initial surface ionic inventories, all populated24-element stages frozen; internal levels and trace stages remain conditionally equilibrated. Use the saved physical rho,T profile and metric. Quantify the initialization change caused by imposing the surface inventory across the local layer.',
        native_grid=dict(log_density=[-6,.02],log_density_nodes=9,temperature_K=[240,20000],temperature_nodes=9),
        evolution=dict(paths=[224,448,896],horizon_seconds=.0034344311179287023,domain_m=[-200,1200],CFL=.35),
        budget=dict(native_calls=2000,native_seconds=60,flow_seconds=120,CPU_threads=1,memory_GB=2),
        forecast='Phase105149 constrained evaluations2.40s.81 grid nodes plus12 controls expected10-35s; new off-adiabat states and conservative inversions are unmeasured. Measure224-cell pilot before either production path.',
        gates=dict(native_constitutive=.002,first_law=1e-4,initial_pressure=.002,primitive_energy=1e-8,baryon=1e-10,energy=1e-8,mass_refinement=.02,trace_refinement=.02),
        stop='No domain/grid/horizon enlargement, no rate borrowing from the different Phase66 EOS, no equilibrium thermal derivative reused at fixed chemical inventory, no final physical charge claim.',
        bindings={str(p):cold.sha(p) for p in [Path(__file__),cold.OUT/'adiabat.npz',cold.OUT/'audit.json',cold.CACHE/'stable/gas.so',cold.CACHE/'libfree_eos_native_cold_stable.so',previous.OUT/'cells-896.npz',previous.OUT/'cells-1792.npz']}))


class Native:
    def __init__(self,cap=2000):
        self.ion=InventoryIons(cap);self.fan=self.ion.fan;self.prefix=np.load(cold.OUT/'adiabat.npz')
        self.lr=np.log(self.fan.rho);self.base=self.ion.snapshot(self.lr,np.log(self.fan.T),np.zeros(318))
        self.target=self.base['number_fractions'];raw=(cold.old.SOURCE/'mod_ionization_data.f90').read_text()
        body=raw[raw.index('monatomic_ip(nions) ='):];body=body[:body.index(']')]
        body='\n'.join(line.split('!')[0] for line in body.splitlines())
        ip=np.array([float(v.replace('d','e')) for v in re.findall(r'([-+]?[0-9]+\.[0-9]*(?:[edED][-+]?[0-9]+)?)_fp_kind',body)])
        assert len(ip)==316
        c2=float(json.loads((cold.old.native.OUT.parent/'def-photon-shared-atomic/catalog.json').read_text())['constants'][2])
        self.energy=np.concatenate([np.cumsum(ip[self.ion.starts[e]:self.ion.starts[e+1]]) for e in range(24)])*c2
        self.charges=np.concatenate([np.arange(1,z+1) for z in self.ion.Z])
        self.active=np.concatenate([self.target[e,1:z+1]>self.target[e].sum()*1e-18 for e,z in enumerate(self.ion.Z)])
        self.diss=float(re.search(r'h2diss = ([0-9.]+)_fp_kind',raw).group(1))*c2
        self.molion=self.diss+ip[0]*c2-21375.95*c2

    def __call__(self,x,lt):
        d=self.prefix;j=int(np.argmin(abs(np.log(d['T'])-lt)))
        last=np.log(d['T'][j]);dr=x-d['log_density_ratio'][j];dt=lt-last;inv=np.exp(-lt)-np.exp(-last)
        fields=d['fields'][j].copy()
        fields[:316]+=np.where(self.active,self.energy*inv+self.charges*(dr-1.5*dt),0)
        fields[316]+=-self.diss*inv-dr+1.5*dt;fields[317]+=self.molion*inv+dr-1.5*dt
        a,_,_=self.ion.constrain(self.lr+x,lt,self.target,fields,target_molecules=self.base['molecular_H_fractions'],tolerance=1e-12)
        return a['eos']


def bank():
    assert not (OUT/'eos.npz').exists();start=time.monotonic();signal.alarm(58);native=Native(cap=1959)
    (OUT/'constrain-source.py').write_text(constrain_source)
    write(OUT/'bank-repair-plan.json',dict(classification='Counterexample candidate',
        failure='At low density and20000K, initially negligible high charge stages become populated. The earlier active-stage-only constraint then cannot preserve total target inventories.',
        intervention='Freeze initially populated stages; suppress newly populated initially negligible stages to their declared1e-18 elemental fraction ceiling. Keep the original full-vector1e-12 acceptance check and16iteration cap. This explicitly extends the trace-stage closure off the original adiabat; it is not a claim that the old unconstrained trace model worked there.',
        reserved_prior_calls=41,remaining_calls=1959,reserved_prior_seconds=2,remaining_seconds=58,
        reuse='Recover the eight accepted first density-column states from their saved native rows.',source_sha256=cold.sha(__file__)))
    x=np.linspace(-6,.02,9);lt=np.linspace(np.log(240.),np.log(20000.),9);raw=np.zeros((9,9,21));done=np.zeros((9,9),bool)
    saved=np.load(OUT/'first-bank-native-states.npz')
    for j,t in enumerate(lt[:8]):
        ids=np.flatnonzero((abs(saved['lrho']-(native.lr+x[0]))<1e-12)&(abs(saved['logT']-t)<1e-12))
        assert len(ids)>0;raw[0,j]=saved['eos'][ids[-1]];done[0,j]=True
    try:
        for i,xx in enumerate(x):
            for j,t in enumerate(lt):
                if not done[i,j]:raw[i,j]=native(float(xx),float(t));done[i,j]=True
            np.savez_compressed(OUT/'bank-progress.npz',x=x,lt=lt,raw=raw,done=done,EOS_calls=native.ion.calls)
        old=previous.Flow(224).eos
        np.savez_compressed(OUT/'eos.npz',x=x,lt=lt,raw=raw,rho0=native.fan.rho,cx=native.fan.cx,
            s0=native.base['eos'][3],sunit=old.sunit,target=native.target,molecules=native.base['molecular_H_fractions'])
        eos=EOS();checks=[]
        for xx,tt in [(-.01,13000),(-.11,14000),(-.37,9000),(-1.1,6900),(-2.3,3500),(-3.1,1900),(-4.2,850),(-5.3,390),(-5.9,260),(-3.7,12000),(-5.5,7000),(-.8,18000)]:
            t=np.log(tt);a=native(xx,t);p,u,gamma,T,kap,cv,s=eos.evaluate(np.array([np.exp(xx)]),np.array([t]))
            h=1e-4;rp=native(xx+h,t);rm=native(xx-h,t);tp=native(xx,t+h);tm=native(xx,t-h)
            ar=(rp-rm)/(2*h);at=(tp-tm)/(2*h);g=ar[1]/a[1]+at[1]/a[1]*(a[1]/a[0]-ar[2])/at[2]
            errors=[float(abs(p[0]*eos.rho0*C*C/a[1]-1)),float(abs(u[0]*C*C/a[2]-1)),float(abs(gamma[0]/g-1)),float(abs(cv[0]*C*C/at[2]-1))]
            laws=[float((tt*at[3]-at[2])/at[2]),float((tt*ar[3]-ar[2]+a[1]/a[0])/(a[1]/a[0]))]
            checks.append(dict(x=xx,T=tt,relative=errors,first_law=laws))
        result=dict(classification='Counterexample candidate',passed=bool(max(max(c['relative']) for c in checks)<.002 and max(max(abs(np.array(c['first_law']))) for c in checks)<1e-4),
            checks=checks,EOS_calls=native.ion.calls,seconds=time.monotonic()-start,source_sha256=cold.sha(__file__))
        write(OUT/'eos.json',result);print(json.dumps(result),flush=True);assert result['passed']
    except Exception as exc:
        write(OUT/'bank-failure.json',dict(error=repr(exc),EOS_calls=native.ion.calls,seconds=time.monotonic()-start,states=int(done.sum())))
        raise
    finally:
        np.savez_compressed(OUT/'bank-native-states.npz',**{k:np.array([s[k] for s in native.ion.states]) for k in native.ion.states[0]})
    signal.alarm(0)


def recover_bank():
    assert not (OUT/'eos.json').exists();saved=np.load(OUT/'bank-native-states.npz');lr=np.log(float(np.load(OUT/'eos.npz')['rho0']))
    failure=json.loads((OUT/'second-bank-failure.json').read_text());assert failure['states']==81
    def native(x,t):
        ids=np.flatnonzero((abs(saved['lrho']-lr-x)<1e-12)&(abs(saved['logT']-t)<1e-12))
        assert len(ids)>0,('Saved constitutive audit state absent',x,t)
        return saved['eos'][ids[-1]]
    native.ion=SimpleNamespace(calls=failure['EOS_calls'])
    # Recover only the completed result serialization from saved exact inputs.
    source=inspect.getsource(bank);section=source[source.index('        eos=EOS();checks=[]'):source.index('    except Exception as exc:')]
    section=textwrap.dedent(section)
    (OUT/'recovered-bank-report.py').write_text(section)
    ns=dict(globals(),native=native,start=time.monotonic()-failure['seconds'])
    exec(compile(section,str(OUT/'recovered-bank-report.py'),'exec'),ns)
    write(OUT/'bank-report-recovery.json',dict(prior_error=failure['error'],new_native_calls=0,source_sha256=cold.sha(__file__),passed=True))


class EOS:
    def __init__(self):
        self.d=np.load(OUT/'eos.npz');d=self.d;self.x=d['x'];self.lt=d['lt'];self.rho0=float(d['rho0']);self.cx=float(d['cx']);self.sunit=float(d['sunit']);self.floor=np.exp(self.x[0]);self.top=np.exp(self.x[-1]);self.fan=SimpleNamespace(calls=0)
        a=d['raw'];rho=self.rho0*np.exp(self.x[:,None]);T=np.exp(self.lt[None,:])
        self.R=float(a[-1,-1,1]/a[-1,-1,0]/T[0,-1]);self.u0=float(a[-1,-1,2]-1.5*self.R*T[0,-1])
        # Interpolate the small nonideal remainder, retaining exact ideal T
        # dependence. No equilibrium EOS derivative enters fixed chemistry.
        self.f=[RectBivariateSpline(self.x,self.lt,v) for v in [np.log(a[:,:,1]/rho/T),a[:,:,2]-self.u0-1.5*self.R*T,a[:,:,3],a[:,:,13]/a[:,:,0]]]

    def limits(self,rho):return np.full_like(rho,self.lt[0]),np.full_like(rho,self.lt[-1])

    def evaluate(self,rho,lt):
        active=rho>=self.floor;assert max(rho)<=self.top*(1+1e-9),'Fixed-inventory density support'
        assert np.all((lt[active]>=self.lt[0]-1e-12)&(lt[active]<=self.lt[-1]+1e-12)),'Fixed-inventory temperature support'
        x=np.log(np.maximum(rho,self.floor));t=np.clip(lt,*self.lt[[0,-1]]);T=np.exp(t);lp,uf,sf,kf=self.f
        p=np.exp(lp.ev(x,t))*self.rho0*rho*T
        u=self.u0+1.5*self.R*T+uf.ev(x,t);cv=1.5*self.R*T+uf.ev(x,t,dy=1);ur=uf.ev(x,t,dx=1)
        chit=1+lp.ev(x,t,dy=1);chir=1+lp.ev(x,t,dx=1)
        gamma=chir+chit*(p/(self.rho0*np.maximum(rho,self.floor))-ur)/cv
        assert np.all(cv[active]>0) and np.all(gamma[active]>1),'Fixed-inventory stability'
        return p/(self.rho0*C*C)*active,u/C**2*active,gamma,T,kf.ev(x,t)*6.6524587321e-25/1.66053906660e-24*active,cv/C**2,sf.ev(x,t)

    def __call__(self,rho,lt):return self.evaluate(rho,lt)[:5]


namespace=dict(previous.old.Flow.primitive.__globals__,OUT=OUT,
    temperature=SimpleNamespace(parent=SimpleNamespace(EOS=EOS)))


class Flow(previous.Flow):
    primitive=FunctionType(previous.old.Flow.primitive.__code__,namespace)
    run=FunctionType(previous.old.Flow.run.__code__,dict(previous.old.Flow.run.__globals__,OUT=OUT))

    def __init__(self,n):
        super().__init__(n);m=self.base;rho=self.initial[0].copy();theta=self.seed.copy();old=self.eos
        p0,u0,*_=old(rho,theta);self.eos=EOS();p,u,*_=self.eos(rho,theta);active=rho>0
        self.initialization=dict(pressure_relative=float(max(abs(p[active]/p0[active]-1))),internal_energy_relative=float(max(abs(u[active]/u0[active]-1))))
        assert self.initialization['pressure_relative']<.002,'Physical initial profile changed'
        self.initial=self.conserved(rho,np.zeros(n),theta,m.a)[0]
        self.background_cell=np.array([rho,np.zeros(n),theta])
        write(OUT/f'initial-{n}.json',dict(classification='Counterexample candidate',**self.initialization,
            boundary='The physical rho,T profile is retained; imposing one surface ion inventory changes its internal energy. This is a frozen-chemistry control, not exact original nonuniform species advection.'))


def run():
    assert json.loads((OUT/'eos.json').read_text())['passed'];assert not (OUT/'result.json').exists();start=time.monotonic();signal.alarm(120)
    pilot=Flow(224).run('pilot-224');forecast=pilot['seconds']*(4+16)*1.4
    write(OUT/'measured-budget.json',dict(pilot_seconds=pilot['seconds'],forecast_remaining_seconds=forecast,remaining_seconds=120-(time.monotonic()-start),assumption='Explicit CFL cell-squared scaling plus40percent; fine inversions are unmeasured.'))
    assert forecast<120-(time.monotonic()-start),'Evolution budget'
    rows=[Flow(n).run(f'cells-{n}') for n in [448,896]]
    errors={k:abs(rows[0][k]/rows[1][k]-1) for k in ['gas_outside_original_radius_g','integrated_trace_energy_erg']}
    result=dict(classification='Counterexample candidate',passed=max(errors.values())<.02,refinement=errors,seconds=time.monotonic()-start,
        actual_frozen_inventory_flow=True,finite_reactions=False,physical_chemistry_closed=False,final_charge_solved=False,full_goal_complete=False)
    write(OUT/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True);assert result['passed']


if __name__=='__main__':globals()[sys.argv[1]]()
