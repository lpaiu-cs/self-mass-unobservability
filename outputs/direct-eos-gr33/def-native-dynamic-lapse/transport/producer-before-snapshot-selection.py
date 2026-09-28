"""Compensated geodesic response on the actual shared photon mesh.

This evolves the geometric transport subsystem. Coupling its increments into
the collision/thermochemical and fluid subsystems remains a separate task.
"""
from pathlib import Path
import json
import signal
import sys
import time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu
import def_native_dynamic_lapse as metric

flow=metric.flow;C=metric.C;LD=metric.LD;write=metric.write;sha=metric.sha
OUT=metric.OUT/'transport';GAMMA=1-1/np.sqrt(2)


def prepare():
    assert not OUT.exists();OUT.mkdir()
    assert json.loads((metric.OUT/'result.json').read_text())['passed']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        claim='Apply the completed lapse/mass/scalar input to the actual shared531-cell,8-angle,152-frequency photon transport operator and integrate its separate distribution increment over the full original horizon.',
        equations='For conserved phase-space packet count N, d(delta_N)/dt=L0*delta_N-div(delta_characteristic_velocity*N0). Use the original common radial/angular operator for L0. The signed frequency generator transfers neighboring energy-node packets and retains the exact number/energy of two ghost exits.',
        source='Actual saved coupled angular spectra drive the perturbation, with the derived radial-speed, angular-drift and frequency-work coefficients. Interpolate stored fields and spectra linearly in time; use interval derivatives of the same lambda and conformal factor, not an inconsistent background-plus-increment sum.',
        boundary='Hold the prescribed incoming packet flux at the deep edge fixed; no incoming external photons. The perturbed outgoing ports are evolved. This is an explicit drive boundary, not a closure of the deeper star.',
        subsystem='The streaming subsystem is solved here. The stored background already includes actual absorption/emission/scattering/heat/H and hydrodynamics, but their additional response to delta_N is not yet in L0. Do not call this full coupled feedback or use it as the final charge.',
        conservation='Count and reference-frequency energy changes equal the same SDIRK-weighted source, radial ports and spectral ghosts. Frequency drift work is a separate metric-work term, never artificial gas heating.',
        controls=dict(paths=[[64,128],[128,128],[128,64]],time=.02,background_time=.02,conservation=1e-10,frequency_moments=1e-12),
        budget=dict(seconds=60,CPU_threads=1,memory_GB=3,new_native_calls=0,new_nonlinear_fluid_steps=0),
        forecast='Original linear band/sparse photon stages are already used by the full solver. These paths omit native roots and nonlinear Newton loops; plan only the three named linear paths, halt60s and preserve any prefix on failure. No automatic extra mesh or paths.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(metric.__file__),metric.OUT/'result.json',metric.OUT/'metric-128-g8.npz',metric.OUT/'metric-64-g8.npz',flow.OUT/'coupled-128.npz',flow.OUT/'coupled-64.npz']}))


class Transport:
    def __init__(self,reference):
        self.model=m=flow.Coupled();b=m.bulk;self.reference=reference;self.n=b.n+m.n;self.q=m.q;self.nf=m.freq
        self.W=np.r_[b.W,m.W];self.area=np.r_[b.area[:-1],m.area];self.mu=b.mu;self.edges_mu=b.edges_mu;self.w=b.w;self.E=b.d['Einf'];self.num=b.d['num']
        self.A=flow.old.transport(self.W,self.area,self.mu,self.edges_mu)
        self.weights=4*np.pi*self.W[:,None,None]*self.w[None,:,None]*self.num[None,None,:]
        data=np.load(flow.OUT/f'coupled-{reference}.npz');self.t=data['snapshot_t']
        self.I=np.concatenate([data['snapshot_bulk_I'],data['snapshot_I'].sum(1)],axis=1)
        self.g=dict(np.load(metric.OUT/'corrected'/f'metric-{reference}-g8.npz'))
        assert len(self.t)==len(self.g['t']) and np.max(abs(self.t-self.g['t']))<1e-18
        self.t=self.g['t'].copy()
        self.r=self.g['radius_E'];bg=m.m.bg.fields(self.r);self.cc=bg['lapse']*np.sqrt(bg['b']);self.nuprime=bg['nu_prime']
        self.extended=np.r_[self.E[0]**2/self.E[1],self.E,self.E[-1]**2/self.E[-2]]

    def frequency(self,I,omega):
        count=I*self.num;pos=omega>=0
        gap=np.where(pos[:,:,None],self.extended[2:]-self.E,self.E-self.extended[:-2])
        moved=count*abs(omega[:,:,None])*self.E/gap
        up=np.where(pos[:,:,None],moved,0.);down=np.where(pos[:,:,None],0.,moved)
        change=-moved.copy();change[:,:,1:]+=up[:,:,:-1];change[:,:,:-1]+=down[:,:,1:]
        escapeN=down[:,:,0]+up[:,:,-1];escapeE=down[:,:,0]*self.extended[0]+up[:,:,-1]*self.extended[-1]
        work=(count*self.E*omega[:,:,None]).sum(2)
        expectedN=change.sum(2)+escapeN;expectedE=(change*self.E).sum(2)+escapeE-work
        number_error=float(np.max(abs(expectedN))/max(np.max(moved.sum(2)),1e-300))
        energy_error=float(np.max(abs(expectedE))/max(np.max(abs(work)),np.max(moved@self.E),1e-300))
        return change/self.num,escapeN,escapeE,work,max(number_error,energy_error)

    def source(self,t):
        j=min(np.searchsorted(self.t,t,side='left')-1,len(self.t)-2);j=max(j,0);h=self.t[j+1]-self.t[j];a=(t-self.t[j])/h
        blend=lambda name:(1-a)*self.g[name][j]+a*self.g[name][j+1]
        I=(1-a)*self.I[j]+a*self.I[j+1];zeta=blend('delta_log_speed');nr=blend('delta_nu_prime');up=blend('delta_u_prime')
        lt=self.g['delta_lambda_interval_rate'][j];ut=(self.g['delta_u'][j+1]-self.g['delta_u'][j])/h
        # The shared face has exactly one perturbation flux. At the two open
        # boundaries the imposed incoming packet perturbation is zero.
        zface=np.interp(np.r_[self.model.bulk.d['edges'][:-1],self.model.m.rf],np.r_[self.model.bulk.d['r'],self.model.m.r],zeta)
        rad=np.zeros((self.n+1,self.q,self.nf));positive=self.mu>0
        rad[1:-1,positive]=I[:-1,positive];rad[1:-1,~positive]=I[1:,~positive]
        rad[-1,positive]=I[-1,positive];rad[0,~positive]=I[0,~positive]
        rad*=C*self.area[:,None,None]*zface[:,None,None]*self.mu[None,:,None]
        source=-np.diff(rad,axis=0)/self.W[:,None,None]
        mu=self.edges_mu[1:-1]
        drift=(1-mu*mu)[None,:]*(C*self.cc[:,None]*((1/self.r-self.nuprime)*zeta-nr)[:,None]-mu[None,:]*lt[:,None])
        angular=np.zeros((self.n,self.q+1,self.nf));angular[:,1:-1]=drift[:,:,None]*I[:,:-1]
        source-=np.diff(angular,axis=1)/np.diff(self.edges_mu)[None,:,None]
        omega=-C*self.cc[:,None]*self.mu[None,:]*(nr+up)[:,None]-ut[:,None]-self.mu[None,:]**2*lt[:,None]
        frequency,en,ee,work,error=self.frequency(I,omega);source+=frequency
        spectralN=float(np.sum(4*np.pi*self.W[:,None]*self.w*en,dtype=LD))
        spectralE=float(np.sum(4*np.pi*self.W[:,None]*self.w*ee,dtype=LD))
        metricwork=float(np.sum(4*np.pi*self.W[:,None]*self.w*work,dtype=LD))
        port=4*np.pi*np.sum((rad[0]-rad[-1])*self.w[:,None]*self.num,axis=0)
        return source,np.array([spectralN,spectralE,metricwork,float(port.sum()-spectralN),float(port@self.E+metricwork-spectralE)]),error

    def moments(self,x):return np.array([np.sum(x*self.weights,dtype=LD),np.sum(x*self.weights*self.E,dtype=LD)],float)

    def port(self,x):
        inner=np.where((self.mu<0)[:,None],x[0],0.);outer=np.where((self.mu>0)[:,None],x[-1],0.)
        packet=4*np.pi*C*np.sum((self.area[0]*inner-self.area[-1]*outer)*(self.w*self.mu)[:,None]*self.num,axis=0)
        return np.array([packet.sum(),packet@self.E])

    def run(self,steps):
        start=time.monotonic();h=self.t[-1]/steps;lu=splu(sparse.eye(self.n*self.q,format='csc')-GAMMA*h*self.A)
        def solve(x):return lu.solve(x.reshape(self.n*self.q,self.nf)).reshape(x.shape)
        x=np.zeros_like(self.I[0]);ledger=np.zeros(2);ghost=np.zeros(3);balance=0.;moment_error=0.;times=[0.];hist=[np.zeros((4,self.n))];source_norm=0.
        for k in range(steps):
            s1,l1,e1=self.source((k+GAMMA)*h);y=solve(x+GAMMA*h*s1)
            first=(self.A@y.reshape(self.n*self.q,self.nf)).reshape(x.shape)+s1
            s2,l2,e2=self.source((k+1)*h);nxt=solve(x+(1-GAMMA)*h*first+GAMMA*h*s2)
            second=(self.A@nxt.reshape(self.n*self.q,self.nf)).reshape(x.shape)+s2
            impulse=h*((1-GAMMA)*(self.port(y)+l1[3:])+GAMMA*(self.port(nxt)+l2[3:]));ledger+=impulse
            ghost+=h*((1-GAMMA)*l1[:3]+GAMMA*l2[:3]);moment_error=max(moment_error,e1,e2)
            source_norm+=h*((1-GAMMA)*np.sum(abs(s1)*self.weights*self.E)+GAMMA*np.sum(abs(s2)*self.weights*self.E))
            x=nxt;balance=max(balance,float(np.max(abs(self.moments(x)-ledger)/np.maximum(abs(ledger),1.))))
            if (k+1)%(steps//16)==0:
                times.append((k+1)*h);hist.append(np.array([np.sum(x*self.weights,axis=(1,2)),np.sum(x*self.weights*self.E,axis=(1,2)),np.sum(x*self.weights*self.E*self.model.bulk.mu2[None,:,None],axis=(1,2)),np.sum(abs(x)*self.weights*self.E,axis=(1,2))]))
        label=f'steps-{steps}-reference-{self.reference}'
        np.savez_compressed(OUT/f'{label}.npz',t=times,moments=hist,delta_count_scaled_occupation=x,number_energy_ledger=ledger,
            spectral_number_energy_and_metric_work=ghost,positive_source_energy_norm=source_norm,radius_E=self.r)
        row=dict(classification='Counterexample candidate',steps=steps,reference=self.reference,seconds=time.monotonic()-start,
            conservation_relative=balance,frequency_moment_relative=moment_error,endpoint_packet_number=float(ledger[0]),endpoint_reference_energy_erg=float(ledger[1]),
            endpoint_spectral_number=float(ghost[0]),endpoint_spectral_energy_erg=float(ghost[1]),integrated_metric_frequency_work_erg=float(ghost[2]),
            endpoint_absolute_energy_erg=float(hist[-1][3].sum()),positive_source_energy_norm_erg=float(source_norm),
            actual_shared_geodesic_transport_evolved=True,collisional_thermochemical_response=False,hydrodynamic_response=False,full_GR_feedback=False,final_charge_solved=False)
        write(OUT/f'{label}.json',row);print(json.dumps(row),flush=True);return row


def run():
    assert not (OUT/'result.json').exists()
    for p,h in json.loads((OUT/'execution-plan.json').read_text())['bindings'].items():assert sha(p)==h,p
    signal.signal(signal.SIGALRM,flow.old.optical.timeout);signal.alarm(52);start=time.monotonic();m=Transport(128)
    rows=[m.run(64),m.run(128),Transport(64).run(128)]
    a=np.load(OUT/'steps-128-reference-128.npz')['moments'];comparisons={}
    for key,name in [('time','steps-64-reference-128'),('background_time','steps-128-reference-64')]:
        other=np.load(OUT/f'{name}.npz')['moments'];comparisons[key]=float(np.max(np.sum(abs(a[:,:3]-other[:,:3]),axis=2)/np.maximum(np.max(np.sum(abs(a[:,:3]),axis=2),axis=0),1.)))
    passed=max(comparisons.values())<.02 and max(r['conservation_relative'] for r in rows)<1e-10 and max(r['frequency_moment_relative'] for r in rows)<1e-12
    result=dict(classification='Counterexample candidate',passed=bool(passed),comparisons=comparisons,paths=rows,seconds=time.monotonic()-start,
        actual_metric_transport_subsystem_completed=True,full_coupled_metric_feedback=False,final_charge_solved=False)
    write(OUT/'result.json',result);print(json.dumps(result),flush=True);signal.alarm(0)


if __name__=='__main__':globals()[sys.argv[1]]()
