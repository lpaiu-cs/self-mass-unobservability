"""Conservative gas evolution using the registered native energy EOS."""
from pathlib import Path
import argparse
import json
import signal
import time
import numpy as np
import def_native_energy_release as bank
import def_native_energy_temperature as temperature

OUT=bank.OUT
C=bank.C
write=bank.write


class Flow:
    def __init__(self,n):
        self.base=m=bank.prior.Flow(n);self.eos=temperature.EOS();self.n=n
        rho,v,entropy=m.primitive(m.initial);sigma=np.log(m.eos(rho,entropy)[3])
        leftT=float(m.eos(np.array([m.left[0]]),np.array([m.left[2]]))[3][0]);m.left=(m.left[0],0.,np.log(leftT))
        self.initial=self.conserved(rho,v,sigma,m.a)[0]
        self.seed=sigma.copy();self.max_recovery=0.;self.max_optical=0.;self.minimum_sigma=float(min(sigma));self.maximum_sigma=float(max(sigma))

    def conserved(self,rho,v,sigma,a):
        p,u,gamma,T,kap=self.eos(rho,sigma);root=np.sqrt(1-v*v);W=1/root;wm=v*v/(root*(1+root));D=rho*W
        h=self.eos.cx+u+p/np.maximum(rho,self.eos.floor);S=rho*h*W*W*v
        K=a*(self.eos.cx*D*wm+(rho*u+p)*W*W-p)+(a-self.base.a0)*self.eos.cx*D
        U=np.array([D,S,K]);F=np.array([D*v,S*v+p,(K+a*p)*v])
        cs=np.sqrt(gamma*p/np.maximum(rho*h,1e-100))
        assert max(cs)<1,'Causal sound speed'
        return U,F,(p,u,gamma,T,kap,cs)

    def primitive(self,U):
        m=self.base;D=np.maximum(U[0],0);active=D>=self.eos.floor
        tau=(U[2]-(m.a-m.a0)*self.eos.cx*D)/m.a
        sigma=self.seed.copy();p=np.zeros_like(D)
        for iteration in range(16):
            v=np.divide(U[1],self.eos.cx*D+tau+p,out=np.zeros_like(D),where=active)
            assert max(abs(v))<1,'Subluminal conservative inverse'
            root=np.sqrt(1-v*v);rho=D*root;W=1/root;wm=v*v/(root*(1+root))
            lo,hi=self.eos.limits(rho);sigma=np.clip(sigma,lo,hi)
            p,u,gamma,T,kap,cvT,entropy=self.eos.evaluate(rho,sigma)
            recovered=self.eos.cx*D*wm+(rho*u+p)*W*W-p
            error=recovered-tau
            relative=np.max(np.where(active,abs(error)/np.maximum(abs(tau),1e-100),0))
            if relative<2e-11:break
            derivative=np.maximum(rho*W*W*cvT,1e-100)
            correction=np.where(active,error/derivative,0)
            # sigma is the reconstructed log temperature in this producer.
            # Domain-bracketed iterates do not clamp the final physical state:
            # a root that cannot satisfy energy within the native bank fails.
            sigma=np.clip(sigma-np.clip(correction,-.5,.5),lo,hi)
        else:
            np.savez_compressed(OUT/'primitive-failure.npz',U=U,rho=rho,sigma=sigma,error=error,tau=tau,relative=relative)
            raise AssertionError(('Conservative primitive root',relative,float(min(sigma[active])),float(max(sigma[active]))))
        self.max_recovery=max(self.max_recovery,float(relative));self.seed=sigma.copy()
        self.minimum_sigma=min(self.minimum_sigma,float(min(sigma[active])));self.maximum_sigma=max(self.maximum_sigma,float(max(sigma[active])))
        return np.array([rho,v,sigma])

    def rhs(self,U,t):
        m=self.base;V=self.primitive(U);rho,v,sigma=V
        _,_,thermo=self.conserved(rho,v,sigma,m.a);p,u,gamma,T,kap,cs=thermo
        ghost=np.array(m.left);ghost[1]=np.interp(t,m.hist_t,m.hist_v)
        ext=np.column_stack([ghost,V,np.zeros(3)]);slope=np.zeros_like(ext)
        slope[:,1:-1]=bank.task.minmod(ext[:,1:-1]-ext[:,:-2],ext[:,2:]-ext[:,1:-1])
        L=ext[:,:-1]+slope[:,:-1]/2;R=ext[:,1:]-slope[:,1:]/2
        UL,FL,tl=self.conserved(*L,m.af);UR,FR,tr=self.conserved(*R,m.af)
        sl=np.minimum(0,np.minimum((L[1]-tl[-1])/(1-L[1]*tl[-1]),(R[1]-tr[-1])/(1-R[1]*tr[-1])))
        sr=np.maximum(0,np.maximum((L[1]+tl[-1])/(1+L[1]*tl[-1]),(R[1]+tr[-1])/(1+R[1]*tr[-1])))
        flux=np.divide(sr*FL-sl*FR+sl*sr*(UR-UL),sr-sl,out=np.zeros_like(FL),where=sr>sl)*m.af*m.area
        rate=-C*np.diff(flux)/m.vol
        W=1/np.sqrt(1-v*v);h=self.eos.cx+u+p/np.maximum(rho,self.eos.floor);E=rho*h*W*W-p
        rate[1]+=C*m.dx/m.vol*(-(m.r/m.RJ)**2*E*m.ap+2*m.a*m.r/m.RJ**2*p)
        area=(m.r/m.RJ)**2;Fbase=m.F0/(area*(m.a/m.a0)**2)
        mu=np.sqrt(np.maximum(0,1-(m.RJ/m.r)**2*(m.a/m.a0)**2));Eg=2*Fbase/(1+mu);Pg=2*Fbase*(1+mu+mu*mu)/(3*(1+mu))
        work=np.zeros(self.n)
        for _ in range(2):
            F=Fbase+np.r_[0,np.cumsum(work[:-1])]/(m.a*m.a*area*C)
            fcom=self.eos.rho0*rho*kap*W*W*((1+v*v)*F-v*(Eg+Pg))
            force=m.a*C*W*fcom;work=-m.vol*m.a*m.a*C*W*v*fcom
        rate[1]+=force;rate[2]-=work/m.vol
        self.max_optical=max(self.max_optical,float(np.sum(m.vol/area*self.eos.rho0*rho*kap)))
        ledger=np.array([C*(flux[0,0]-flux[0,-1]),C*(flux[2,0]-flux[2,-1])-sum(work),-sum(work)])
        dt=.35*m.dx/np.max(C*m.a/m.B*(abs(v)+cs+1e-100))
        return rate,ledger,dt

    def run(self,label):
        assert not (OUT/f'{label}.json').exists();m=self.base;start=time.monotonic();U=self.initial.copy();t=0.;end=float(m.hist_t[-1]);ledger=np.zeros(3);discard=np.zeros(3);steps=0;next_dump=0.;history=[];snapshots=[]
        while t<end:
            k,l,dt=self.rhs(U,t);dt=min(dt,end-t);trial=U+dt*k
            assert min(trial[0])>-1e-13,'Density positivity'
            k2,l2,_=self.rhs(trial,t+dt);nxt=(U+trial+dt*k2)/2
            tiny=nxt[0]<self.eos.floor;discard+=np.sum(nxt[:,tiny]*m.vol[tiny],axis=1);nxt[:,tiny]=0
            U=nxt;ledger+=dt*(l+l2)/2;t+=dt;steps+=1
            assert steps<20000,'Evolution step cap'
            if t>=next_dump or t>=end:
                outside=np.sum(U[0,m.x>=0]*m.vol[m.x>=0]);energy=np.sum((U[2]-self.initial[2])*m.vol)
                history.append([t,outside,energy,*ledger,*discard]);snapshots.append(U.copy());next_dump+=end/32
        rho,v,sigma=self.primitive(U);p,u,gamma,T,kap=self.eos(rho,sigma)
        ir,iv,isent=self.base.primitive(self.base.initial);iT=self.base.eos(ir,isent)[3];ip,iu,*_=self.eos(ir,np.log(iT))
        trace=-self.eos.cx*U[0]*v*v/(1+np.sqrt(1-v*v))+rho*u-3*p-(self.base.initial[0]*iu-3*ip)
        scale=4*np.pi*m.RJ**2*self.eos.rho0
        energy_residual=np.sum((U[2]-self.initial[2])*m.vol)+discard[2]-ledger[1]
        response=max(np.sum(m.vol*(rho*v*v+p)),abs(history[-1][2]),1e-100)
        baryon=abs(np.sum((U[0]-self.initial[0])*m.vol)+discard[0]-ledger[0])/np.sum(self.initial[0]*m.vol)
        energy_error=abs(energy_residual)/response
        entropy=(self.eos.evaluate(rho,sigma)[-1]-float(self.eos.d['s0']))/self.eos.sunit
        np.savez_compressed(OUT/f'{label}.npz',U=U,initial=self.initial,x_cm=m.x,volume=m.vol,rho=rho*self.eos.rho0,velocity_cm_s=v*C,logT=sigma,sigma=entropy,T=T,pressure=p*self.eos.rho0*C*C,
            history=history,snapshots=snapshots,conserved_discard=discard,initial_entropy=self.base.primitive(self.base.initial)[2])
        passed=baryon<1e-10 and energy_error<1e-8 and self.max_recovery<1e-8
        row=dict(classification='Counterexample candidate',passed=bool(passed),cells=self.n,steps=steps,seconds=time.monotonic()-start,
            baryon_ledger_relative=float(baryon),energy_ledger_over_response=float(energy_error),energy_ledger_residual_erg=float(energy_residual*scale*C*C),
            maximum_primitive_energy_relative=self.max_recovery,minimum_logT=self.minimum_sigma,maximum_logT=self.maximum_sigma,
            final_minimum_sigma=float(min(entropy[rho>=self.eos.floor])),final_maximum_sigma=float(max(entropy[rho>=self.eos.floor])),
            gas_outside_original_radius_g=float(history[-1][1]*scale),integrated_trace_energy_erg=float(np.sum(trace*m.vol)*scale*C*C),
            dilute_baryon_g=float(discard[0]*scale),dilute_Killing_nonrest_energy_erg=float(discard[2]*scale*C*C),
            inner_baryon_into_layer_g=float(ledger[0]*scale),inner_Killing_nonrest_energy_into_layer_erg=float((ledger[1]-ledger[2])*scale*C*C),
            scattering_work_into_gas_erg=float(ledger[2]*scale*C*C),maximum_scattering_optical_depth=self.max_optical,
            full_GR_scalar_feedback=False,physical_chemistry_closed=False,final_charge_solved=False,full_goal_complete=False)
        write(OUT/f'{label}.json',row);print(label,json.dumps(row),flush=True);assert passed
        return row


def run():
    assert not (OUT/'result.json').exists();assert json.loads((OUT/'temperature-audit.json').read_text())['passed'];start=time.monotonic();signal.alarm(240)
    write(OUT/'flow-source-binding.json',dict(source_sha256=bank.old.photons.digest(Path(__file__)),EOS_source_sha256=bank.old.photons.digest(Path(temperature.__file__)),eos_sha256=bank.old.photons.digest(OUT/'temperature-columns.npz'),primitive_variable='log temperature; sigma variable name in the producer refers to logT, entropy is reconstructed separately from the native table')))
    pilot=Flow(224).run('pilot-224');forecast=pilot['seconds']*(1+16+64)*1.3
    write(OUT/'flow-budget.json',dict(pilot_seconds=pilot['seconds'],forecast_seconds=forecast,assumption='Cell-squared explicit-CFL cost at the same horizon,30percent margin. Finer-grid recovery cost unmeasured.'))
    assert forecast<235,'Conservative flow forecast exceeds budget'
    rows=[Flow(n).run(f'cells-{n}') for n in [896,1792]]
    errors={k:abs(rows[0][k]/rows[1][k]-1) for k in ['gas_outside_original_radius_g','integrated_trace_energy_erg']}
    result=dict(classification='Counterexample candidate',passed=all(r['passed'] for r in rows) and max(errors.values())<.02,refinement=errors,
        seconds=time.monotonic()-start,conservative_energy_evolved=True,energy_defect_added_as_heat=False,full_GR_scalar_feedback=False,final_charge_solved=False,full_goal_complete=False)
    write(OUT/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['run']);globals()[p.parse_args().action]()
