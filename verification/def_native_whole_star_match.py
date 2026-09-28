"""Fixed material core plus native radiative envelope: six-variable match.

Counterexample candidate: quasistatic grey Cauchy background, not thermal
equilibrium or a zero-pressure material/vacuum junction. The core entropy is
retained; the reconstructed envelope entropy is allowed to change.
"""
from pathlib import Path
import argparse
import json
import signal
import time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
import def_native_radiative_envelope as envelope

h=envelope.prior.h
ROOT=envelope.ROOT
OUT=envelope.OUT.parent/'def-native-whole-star-match'
write=envelope.write


class Match:
    def __init__(self,sub=4,tol=2e-6,cut=1.):
        self.core=h.Structure(.001,sub)
        self.env=envelope.Envelope()
        self.tol,self.cut=tol,cut
        self.i=132
        self.native_fraction=float(self.core.outer[self.i]+self.core.f[self.i]/2)
        old=json.loads((OUT.parent/'def-native-atmosphere/density-coordinate/result.json').read_text())
        self.old_atmosphere_fraction=old['rows'][-1]['added_baryon_fraction']
        self.target=self.native_fraction+self.old_atmosphere_fraction
        self.history=[]
        self.last=None

    def initial(self):
        s=self.core
        old=np.load(h.OUT/'absolute-shoot/background-0.001.npz')['parameters']
        d=np.load(envelope.OUT/'fine.npz')
        return np.array([old[0],np.log(d['r'][-1]/(s.R*100)),
            np.log(d['m'][-1]/(s.M*100)),old[3],d['r'][-1]*d['v'][-1]/(.001*s.mu),np.log(float(d['Linfinity']))])

    def outer(self,x,record=False):
        s,p=self.core,self.env
        _,lr,lm,_,j,ll=x
        R=s.R*100*np.exp(lr);M=s.M*100*np.exp(lm);q=.001*s.mu*j
        h.mp.mp.dps=50
        exact=h.exterior.exact(h.mp.mpf(M/R),h.mp.mpf(q))
        phi=.001-q*float(exact[2])
        nu=float(h.mp.log1p(-2*h.mp.mpf(M/R))/2-h.mp.mpf(q)**2*exact[3])
        p.R,p.m,p.phi,p.v,p.nu=R,M,phi,q/R,nu
        L=np.exp(ll);A=np.exp(-2*phi*phi)
        p.T=(L/(8*np.pi*R*R*np.exp(2*nu)*A**4*p.sigma))**.25
        y0=np.zeros(8);y0[2]=np.log(p.T)
        def rhs(logpg,y):
            z=p.state(logpg,y,L)
            r,m,phi,v,A,b=[z[k] for k in ['r','m','phi','v','A','b']]
            dr=z['Pgas']/z['gas_pressure_prime']
            mr=4*np.pi*r*r*A**4*z['egeo']+r*r*b*v*v/2
            vr=4*np.pi*A**4/b*(-4*phi*(z['egeo']-3*z['pgeo'])+r*v*(z['egeo']-z['pgeo']))-2*(r-m)/(r*r*b)*v
            return dr*np.array([1/p.R,mr/p.m,z['logT_prime'],v/.001,
                vr*p.R/.001,z['nr'],4*np.pi*r*r*A**3*z['rho']/np.sqrt(b)/p.total_baryon,-z['optical']])
        def baryons(logpg,y):return y[6]+self.target
        baryons.terminal=True;baryons.direction=-1
        sol=solve_ivp(rhs,(np.log(self.cut),np.log(2e6)),y0,method=envelope.DomainRK45,
            rtol=self.tol,atol=[1e-13,1e-23,1e-12,1e-19,1e-19,1e-19,1e-26,1e-10],
            max_step=.3,events=baryons,dense_output=True)
        assert sol.success and sol.status==1,(sol.message,sol.y[6,-1],self.target)
        z=p.state(sol.t[-1],sol.y[:,-1],L)
        state=np.array([z['r']/(s.R*100),z['m']/(s.B*100),np.log(z['Ptotal']),
            -j*float(exact[2])+sol.y[3,-1]/s.mu,z['v']*s.R*100/(.001*s.mu)])
        target_temperature=s.state(state[2],self.i)[3][1]
        thermal=np.log(z['T'])-target_temperature
        saved=None
        if record:
            grid=np.linspace(sol.t[-1],np.log(self.cut),401);states=sol.sol(grid)
            values=[p.state(a,b,L) for a,b in zip(grid,states.T)]
            saved=dict(logPgas=grid,states=states,**{k:np.array([v[k] for v in values]) for k in values[0]},
                Linfinity=L,baryon_outside_base=self.target+states[6])
        return state,float(thermal),z,saved

    def branches(self,x,record=False):
        s=self.core
        pc,_,_,uc,_,_=x
        p,en,rho,_=s.state(pc,len(s.f)-1)
        phic=.001*(1+s.mu*uc);A=np.exp(-2*phic**2)
        pE,enE=A**4*p,A**4*en;r0=1.
        Phi=4*np.pi/3*s.beta*(1+s.mu*uc)*(enE-3*pE)*r0
        q0=4*np.pi*A**3*rho*r0**3/(3*s.B)
        drop=2*np.pi/3*(enE+3*pE)*r0*r0+s.beta*phic*.001*Phi*r0/2
        yc=np.array([r0/s.R,4*np.pi/3*enE*r0**3/s.B,
            pc-(en+p)/p*drop,uc+Phi*r0/(2*s.mu),Phi*s.R/s.mu])
        inner=[(q0,*yc)]
        for i in range(len(s.f)-1,s.split-1,-1):
            low=q0 if i==len(s.f)-1 else s.inner[i+1]
            yc=s.step(np.log(low),np.log(s.inner[i]),yc,i,s.B,False)
            if record:inner.append((s.inner[i],*yc))
        yo,thermal,z,saved=self.outer(x,record)
        outer=[(self.native_fraction,*yo)]
        for i in range(self.i,s.split):
            low=self.native_fraction if i==self.i else s.outer[i]
            yo=s.step(np.log(low),np.log(s.outer[i+1]),yo,i,s.B,True)
            if record:outer.append((s.outer[i+1],*yo))
        return np.r_[yc-yo,thermal],np.array(inner),np.array(outer),z,saved

    def objective(self,x):
        assert len(self.history)<45,'Registered objective call cap'
        start=time.monotonic();error,*_=self.branches(x)
        self.history.append(dict(parameters=x.tolist(),residual=error.tolist(),seconds=time.monotonic()-start))
        write(OUT/'progress.json',dict(classification='Counterexample candidate',history=self.history,native_calls=self.env.calls))
        print('MATCH',len(self.history),max(abs(error)),error.tolist(),flush=True)
        return error

    def save(self,x,label):
        err,inner,outer,z,saved=self.branches(x,True)
        s=self.core;n=len(s.f)
        # Keep the original baryon cell labels. The first retained midpoint is
        # exactly the envelope junction; native cells0..131 are replaced.
        values=[];thermo=[]
        for i in range(self.i,n):
            if i==self.i:
                value=outer[0,1:]
            else:
                outside=i<s.split
                if outside:
                    item=outer[i-self.i];begin=item[0];value=item[1:]
                    mid=(s.outer[i]+s.outer[i+1])/2
                else:
                    item=inner[n-1-i];begin=item[0];value=item[1:]
                    mid=(s.inner[i]+s.inner[i+1])/2
                value=s.step(np.log(begin),np.log(mid),value,i,s.B,outside)
            values.append(value);thermo.append(s.state(value[2],i)[3])
        np.savez_compressed(OUT/f'{label}-core.npz',states=values,thermo=thermo,
            inner=inner,outer=outer,parameters=x,original_indices=np.arange(self.i,n),
            dm=np.r_[s.data['dm'][self.i]/2,s.data['dm'][self.i+1:]],X=s.data['X'][self.i:],
            entropy=s.ref[self.i:,3],R_scale_m=s.R,B_scale_m=s.B,mu_scale=s.mu)
        np.savez_compressed(OUT/f'{label}-envelope.npz',**saved)
        q=.001*s.mu*x[4];R=s.R*np.exp(x[1]);M=s.M*np.exp(x[2])
        ext=h.exterior.exact(h.mp.mpf(M/R),h.mp.mpf(q))
        adm=M+R*q*q*float(ext[0]);charge=-R*q*float(ext[1])
        row=dict(classification='Counterexample candidate',parameters=x.tolist(),residual=err.tolist(),
            radius_m=R,ADM_geom_m=adm,scalar_charge_geom_m=charge,alpha=-charge/adm,
            L_infinity=float(np.exp(x[5])),base=z,envelope_baryon_fraction=self.target,
            native_outer_replaced_fraction=self.native_fraction,old_native_atmosphere_fraction=self.old_atmosphere_fraction,
            retained_core_baryon_g=float(s.data['dm'][self.i]/2+s.data['dm'][self.i+1:].sum()),
            reconstructed_envelope_baryon_g=self.target*self.env.total_baryon,
            whole_material_g=(1+self.old_atmosphere_fraction)*self.env.total_baryon,
            core_entropy_and_isotope_template_retained=True,thermal_stationarity=False,
            full_dynamic_charge_solved=False,full_goal_complete=False)
        write(OUT/f'{label}.json',row)
        return row


def prepare():
    assert not OUT.exists();OUT.mkdir()
    files=[Path(__file__),Path(envelope.__file__),Path(h.__file__),
        h.OUT/'absolute-shoot/background-0.001.npz',h.OUT/'extended-table.npz',
        h.OLD/'molecular-state-17-8.npz',h.OLD/'molecular-adiabats-17.npz',
        envelope.OUT/'fine.npz',OUT.parent/'def-native-atmosphere/density-coordinate/result.json']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='1a36074b8',
        claim='Solve a whole fixed-inventory mechanical background with an actual native radiative envelope and regular scalar centre/asymptotic field. Eliminate artificial envelope addition or scaling of all original baryons.',
        material='Retain native cell132 inner half plus cells133..5734, their isotope fractions and entropy. Replace the remaining native cells and the previously integrated native atmosphere by a grey envelope containing exactly the same total baryons and uniform outer composition. The older unresolved cold tail remains separately bounded.',
        unknowns='central logP, surface logR/logM, central scaled scalar, surface scalar flux, logL; five mechanical junction conditions plus native entropy/temperature continuity at the fixed material junction.',
        thermal='Solve Eddington radiation transport in the envelope. Core entropy remains its initial profile, so the core need not be in thermal equilibrium. Determine a single common base face flux and explicitly compute residual core heating rather than deleting a flux mismatch.',
        gates=dict(junction_maximum=1e-8,base_logT=1e-7,inventory_relative_to_envelope=1e-8,
            refinement_logL=1e-4,refinement_logR=1e-6,refinement_alpha_relative=1e-7,
            native_core_logrho=1e-7,native_core_entropy_over_cv=1e-7),
        budget=dict(total_seconds=600,pilot_seconds=60,maximum_objective_calls=45,native_calls=60000,
            BLAS_threads=1,automatic_expansion=False,old_long_evolutions=0),
        production='One six-variable root using existing core RK4/sub4 and native envelope tolerance2e-6; one fixed-parameter core sub8/envelope2e-7 refinement, no automatic new fit.',
        decision='If material and thermal junctions pass, use the saved whole-star background for revised coupled response; do not count a finite-pressure grey surface as a complete physical atmosphere.',
        bindings={str(p.relative_to(ROOT)):envelope.digest(p) for p in files},
        symbolic=envelope.symbolic()))


def pilot():
    assert not (OUT/'pilot.json').exists();signal.alarm(60);start=time.monotonic()
    model=Match();x=model.initial();before=time.monotonic();error,*_=model.branches(x)
    elapsed=time.monotonic()-before
    write(OUT/'pilot.json',dict(classification='Counterexample candidate',parameters=x.tolist(),
        residual=error.tolist(),seconds=time.monotonic()-start,objective_seconds=elapsed,
        native_calls=model.env.calls,old_atmosphere_fraction=model.old_atmosphere_fraction,
        target_envelope_fraction=model.target,
        forecast_45_objectives_seconds=45*elapsed*1.25+60))
    signal.alarm(0);print((OUT/'pilot.json').read_text(),flush=True)


def run():
    assert not (OUT/'result.json').exists()
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,sha in plan['bindings'].items():assert envelope.digest(ROOT/rel)==sha,rel
    pilot=json.loads((OUT/'pilot.json').read_text())
    assert pilot['forecast_45_objectives_seconds']+pilot['seconds']<600,'Reassess budget, no automatic expansion'
    write(OUT/'budget.json',dict(classification='Counterexample candidate',pilot=pilot['seconds'],
        forecast_seconds=pilot['forecast_45_objectives_seconds'],upper_bound_seconds=600,
        assumption='Same cached core and native atmosphere speed plus25percent, extra60s for saved profiles and a fixed-parameter refinement; unmeasured paths may differ.',proceed=True))
    signal.alarm(int(600-pilot['seconds']));start=time.monotonic();model=Match()
    steps=np.array([2e-6,2e-6,2e-6,.001,.0001,2e-5])
    def jac(x):
        columns=[]
        for i,step in enumerate(steps):
            v=np.zeros(6);v[i]=step
            columns.append((model.objective(x+v)-model.objective(x-v))/(2*step))
        return np.array(columns).T
    fit=root(model.objective,np.array(pilot['parameters']),jac=jac,method='hybr',
        options=dict(xtol=1e-10,maxfev=30))
    value=model.save(fit.x,'matched')
    passed=max(abs(np.array(value['residual'])))<1e-8
    assert passed,(fit.message,value['residual'])
    fine=Match(sub=8,tol=2e-7);refined=fine.save(fit.x,'refined')
    errors=dict(junction=float(max(abs(np.array(refined['residual'])))),
        base_logT=abs(refined['residual'][-1]),
        radius_log=abs(np.log(refined['base']['r']/value['base']['r'])),
        material=float(abs(refined['envelope_baryon_fraction']/value['envelope_baryon_fraction']-1)))
    # Fixed-parameter residual checks independently test the matched root.
    accepted=errors['junction']<1e-8 and errors['base_logT']<1e-7 and errors['material']<1e-8
    result=dict(classification='Counterexample candidate',passed=bool(accepted),
        parameters=fit.x.tolist(),root_message=str(fit.message),objective_calls=len(model.history),
        native_calls=model.env.calls+fine.env.calls,matched=value,refined=refined,
        refinement=errors,seconds=time.monotonic()-start,total_compute_seconds=pilot['seconds']+time.monotonic()-start,
        whole_fixed_material_background_matched=bool(accepted),thermal_stationarity=False,
        finite_pressure_surface=True,full_goal_complete=False)
    write(OUT/'result.json',result);signal.alarm(0)
    print('FINAL',json.dumps(result),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','pilot','run'])
    globals()[parser.parse_args().action]()
