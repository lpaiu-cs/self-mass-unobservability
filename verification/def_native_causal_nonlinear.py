"""Counterexample candidate: nonlinear heat/H/photon exchange on saved radii."""
from pathlib import Path
import json
import signal
import sys
import time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu
import def_native_causal_photons as prior

old=prior.old
OUT=prior.OUT/'nonlinear'
C=prior.C
write=prior.write


def prepare():
    assert not OUT.exists();OUT.mkdir()
    failed=json.loads((prior.OUT/'pilot.json').read_text());assert not failed['passed']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        failure='The actual16cell spectral/heat/H pilot conserves energy but gives2.55percent logT change and61.7percent neutral change, exceeding the original tangent gates. Preserve it; no spatial refinement can repair that physical approximation.',
        repair='Replace tangent thermochemistry by same-native nonlinear H inventories and gas energy, with nonlinear stimulated/emitted photon rates. Reuse the identical16cells,152frequencies,metric,Cauchy field and time interval. The temperature and chemical response are solved implicitly together with photons.',
        support=dict(logT_offset=[-.06,.06],neutral_ratio=[.25,2.5],temperature_nodes=3,chemical_planes=2),
        gates=dict(native_constitutive=.002,native_rate=.002,newton_scaled=1e-9,energy=1e-8,time_refinement=.02,positive_populations=True),
        budget=dict(new_native_calls=1200,new_native_seconds=45,nonlinear_seconds=120,paths_time_steps=[64,128],cells=16,CPU_threads=1,memory_GB=3),
        resource_reassessment='Native full16cell bank and process startup took7.10s before a NumPy-bool report serialization failure. The saved coefficients are reused. Its exact public-call counter was not serialized; its3500call cap is an upper bound, not an observed count. This new nonlinear constitutive task has its own explicit1200call cap; combined worst-case4700, actual first count unknown. No completed native bank or fluid history is recomputed.',
        prior_operational_failures=['NumPy bool report serialization; bank NPZ survived. Future bank NPZ now saves call count and audit errors before report serialization.','The first angular precondition used the half-isotropic surface ratio0.5 as an interior maximum. Correct P1 positive-moment condition is(H/J)^2<=1/3; measured maximum0.502277 is allowed. No numerical trajectory was run under the mistaken check.'],
        forecast='Pilot factor0.051s; three Newton factors per step suggest29.4s for64+128steps; state-dependent factor cost and iteration count unmeasured. Hard120s alarm, maximum10 Newton iterations and no automatic extra path.',
        limits='No new spatial convergence claim from16cells; P1 angular closure, prescribed initial spectrum, frozen mechanical density/metric, other ionic inventories frozen and reciprocal rather than microscopic reverse rates remain. Final charge is not solved.',
        bindings={str(p):prior.sha(p) for p in [Path(__file__),Path(prior.__file__),prior.OUT/'initial-producer.py',prior.OUT/'bank-16-8.npz',prior.OUT/'pilot.json']}))


def setup(native,d,j):
    rho,T=d['rho'][j],d['T'][j];eq=native.ion.snapshot(np.log(rho),np.log(T),np.zeros(318));target=eq['number_fractions']
    native.lr=np.log(rho);native.base=eq;native.target=target;native.epsH=float(target[0,:2].sum());native.y0=float(target[0,0]/native.epsH)
    native.nH=native.epsH*native.fan.cx*old.NA
    native.active=np.concatenate([target[e,1:z+1]>target[e].sum()*1e-18 for e,z in enumerate(native.ion.Z)])
    native.prefix=dict(T=np.array([T]),fields=np.zeros((1,318)),log_density_ratio=np.zeros(1))
    assert abs(native.y0/d['y0'][j]-1)<1e-10


def bank():
    start=time.monotonic();signal.alarm(45);assert not (OUT/'bank.npz').exists()
    d=np.load(prior.OUT/'bank-16-8.npz');native=old.Native(cap=1200);n=len(d['r']);m=len(d['Einf'])
    offsets=np.array([-.06,0.,.06]);ratios=np.array([.25,2.5]);raw=np.zeros((n,2,3,21));rates=np.zeros((n,2,3,2,m))
    for j in range(n):
        setup(native,d,j)
        for iy,factor in enumerate(ratios):
            for it,theta in enumerate(offsets):
                y=native.y0*factor;s=native.state(0.,np.log(d['T'][j])+theta,y);k,e=prior.coefficients(native,s,d['Einf']/d['a'][j])
                raw[j,iy,it]=s['raw'];rates[j,iy,it]=[(k+e)/y,e/(1-y)]
    np.savez_compressed(OUT/'bank.npz',raw=raw,rates=rates,offsets=offsets,ratios=ratios,native_calls=native.ion.calls)
    eos=Table();checks=[]
    for j in range(n):
        setup(native,d,j);theta=.017 if j%2 else -.023;eta=.4 if j%2 else -.15
        s=native.state(0.,np.log(d['T'][j])+theta,native.y0*(1+eta));k,e=prior.coefficients(native,s,d['Einf']/d['a'][j]);a=k+e
        t=np.full(n,theta);y=np.full(n,eta);p,u,*_=eos.gas(t,y);aa,ee,*_=eos.radiation(t,y)
        mask=(a>0)&(e>0);rate=float(max(np.max(abs(aa[j,mask]/a[mask]-1)),np.max(abs(ee[j,mask]/e[mask]-1))))
        checks.append(dict(cell=j,constitutive=float(max(abs(p[j]/s['raw'][1]-1),abs(u[j]/s['raw'][2]-1))),rate=rate))
    # All saved equilibrium states independently check the two chemical planes.
    p,u,*_=eos.gas(np.zeros(n),np.zeros(n));aa,ee,*_=eos.radiation(np.zeros(n),np.zeros(n))
    bb=1/np.expm1(d['Einf'][None,:]/(d['a']*prior.K*d['T'])[:,None])
    balance=float(np.max(abs(ee-(aa-ee)*bb)/(ee+abs((aa-ee)*bb)+1e-100)))
    baseline=float(max(np.max(abs(p/d['raw'][:,1]-1)),np.max(abs(u/d['raw'][:,2]-1))))
    result=dict(classification='Counterexample candidate',passed=bool(max(a['rate'] for a in checks)<.002 and max(a['constitutive'] for a in checks)<.002 and baseline<.002),
        checks=checks,baseline_constitutive_relative=baseline,unadjusted_LTE_balance_relative=balance,
        native_calls=native.ion.calls,seconds=time.monotonic()-start,
        note='The raw interpolation imbalance is recorded. Evolution evolves this declared interpolant directly; it is not removed by subtracting a source residual.',source_sha256=prior.sha(__file__))
    write(OUT/'bank.json',result);signal.alarm(0);print(json.dumps(result),flush=True);assert result['passed']


class Table:
    def __init__(self):
        self.d=d=np.load(prior.OUT/'bank-16-8.npz');v=np.load(OUT/'bank.npz');self.raw=v['raw'];self.rates=v['rates'];self.ratios=v['ratios'];self.n=len(d['r'])
        h=.06
        def poly(v):return np.stack([v[:,:,1],(v[:,:,2]-v[:,:,0])/(2*h),(v[:,:,2]+v[:,:,0]-2*v[:,:,1])/(2*h*h)],axis=2)
        rho=d['rho'][:,None];T=d['T'][:,None]
        self.R=self.raw[:,:,1,1]/(rho*T);self.u0=self.raw[:,:,1,2]-1.5*self.R*T
        temps=T[:,:,None]*np.exp(v['offsets'])[None,None,:]
        self.uc=poly(self.raw[:,:,:,2]-self.u0[:,:,None]-1.5*self.R[:,:,None]*temps)
        self.pc=poly(np.log(self.raw[:,:,:,1]/(rho[:,:,None]*temps)))
        self.nec=poly(self.raw[:,:,:,13]/1.66053906660e-24)
        self.mask=self.rates[:,:,1]>0
        logs=np.log(np.maximum(self.rates,1e-300))
        # Keep the exact photon Boltzmann exponential out of interpolation.
        # Otherwise a tiny relative temperature interpolation error is
        # multiplied by h*nu/kT in the Wien tail, even with perfect EOS data.
        logs[:,:,:,1]+=d['Einf'][None,None,None,:]/(d['a'][:,None,None,None]*prior.K*temps[:,:,:,None])
        self.rc=poly(logs)

    def polynomial(self,co,t):
        tt=t.reshape((self.n,1)+(1,)*(co.ndim-3));value=co[:,:,0]+tt*(co[:,:,1]+tt*co[:,:,2]);der=co[:,:,1]+2*tt*co[:,:,2]
        return value,der

    def mix(self,values,eta):
        w=((1+eta-self.ratios[0])/np.diff(self.ratios)[0]).reshape((self.n,)+(1,)*(values.ndim-2))
        return values[:,0]*(1-w)+values[:,1]*w,(values[:,1]-values[:,0])/np.diff(self.ratios)[0]

    def gas(self,theta,eta):
        T=self.d['T'][:,None]*np.exp(theta[:,None]);q,qt=self.polynomial(self.uc,theta)
        u=self.u0+1.5*self.R*T+q;ut=1.5*self.R*T+qt
        lp,lpt=self.polynomial(self.pc,theta);p=np.exp(lp)*self.d['rho'][:,None]*T
        en,ent=self.polynomial(self.nec,theta)
        uu,uy=self.mix(u,eta);pu,py=self.mix(p,eta);uuT,_=self.mix(ut,eta);pt,_=self.mix(p*(1+lpt),eta)
        ne,ney=self.mix(en,eta);net,_=self.mix(ent,eta)
        return pu,uu,uuT,uy,pt,py,ne,net,ney

    def radiation(self,theta,eta):
        lp,lpt=self.polynomial(self.rc,theta)
        boltz=self.d['Einf'][None,:]/((self.d['a']*prior.K*self.d['T']*np.exp(theta))[:,None])
        lp[:,:,1]-=boltz[:,None,:];lpt[:,:,1]+=boltz[:,None,:]
        rates=np.exp(lp)*self.mask;ratesT=rates*lpt
        v,vy=self.mix(rates,eta);vt,_=self.mix(ratesT,eta);y=self.d['y0']*(1+eta);dy=self.d['y0']
        a=v[:,0]*y[:,None];e=v[:,1]*(1-y[:,None]);at=vt[:,0]*y[:,None];et=vt[:,1]*(1-y[:,None])
        ay=vy[:,0]*y[:,None]+v[:,0]*dy[:,None];ey=vy[:,1]*(1-y[:,None])-v[:,1]*dy[:,None]
        return a,e,at,et,ay,ey

    def temperature(self,energy,eta):
        theta=np.zeros(self.n)
        for _ in range(8):
            _,u,ut,*_=self.gas(theta,eta);step=(u-energy)/ut;theta-=step
            if max(abs(step))<1e-12:break
        else:raise AssertionError('Conservative gas temperature root')
        assert max(abs(theta))<=.06 and np.all((eta>=-.75)&(eta<=1.5)),('Nonlinear support',min(theta),max(theta),min(eta),max(eta))
        return theta


class Model(prior.Model):
    def __init__(self):
        super().__init__('bank-16-8');self.eos=Table();d=self.d;n=self.n;m=self.m
        self.u0=self.eos.gas(np.zeros(n),np.zeros(n))[1];self.uscale=d['thermo'][:,0]
        self.num=d['num'][None,:]*self.scale/(d['a']**3*d['rho']*d['thermo'][:,4]*d['y0'])[:,None]
        self.energy=-d['num'][None,:]*d['Einf']*self.scale/(d['a']**4*d['rho']*self.uscale)[:,None]
        # Reuse the exact transport entries; remove all material/collision
        # entries from the tangent matrix rather than duplicating streaming.
        A=self.A.tolil();A[:,self.it]=0.;A[:,self.iy]=0.;A[self.it,:]=0.;A[self.iy,:]=0.
        I=self.I;F=self.F
        for i in range(n):A[I[i],I[i]]=0.
        for i in range(n+1):A[F[i],F[i]]=0.
        self.transport=A.tocsc();initial=np.zeros(self.size);initial[:self.nj]=(self.J/self.scale).ravel();initial[self.nj:self.nj+self.nh]=(self.H/self.scale).ravel()
        self.forcing=self.transport@initial

    def rhs(self,state,jacobian=True):
        d=self.d;n=self.n;I=self.I;F=self.F;eta=state[self.iy];theta=self.eos.temperature(self.u0+self.uscale*state[self.it],eta)
        p,u,ut,uy,pt,py,ne,net,ney=self.eos.gas(theta,eta)
        ab,em,at,et,ay,ey=self.eos.radiation(theta,eta);J=state[:self.nj].reshape(n,self.m)*self.scale+self.J
        HH=state[self.nj:self.nj+self.nh].reshape(n+1,self.m)+self.H/self.scale
        fac=d['a']*C*d['rho']*d['thermo'][:,4]
        collision=fac[:,None]*(em*(1+J)-ab*J)/self.scale
        base=self.transport@state+self.forcing;base[I]+=collision;base[self.it]=(self.energy*collision).sum(1);base[self.iy]=(self.num*collision).sum(1)
        opacity=d['rho'][:,None]*d['thermo'][:,4,None]*(ab-em)+ne[:,None]*6.6524587321e-25
        rate=d['face_a'][1:-1,None]*C*(opacity[:-1]+opacity[1:])/2;base[F[1:-1]]-=rate*HH[1:-1]
        if not jacobian:return base,theta,p,u
        rows=[];cols=[];vals=[]
        def add(row,col,value):
            rr,cc,vv=np.broadcast_arrays(row,col,value);rows.extend(rr.ravel());cols.extend(cc.ravel());vals.extend(vv.ravel())
        cj=fac[:,None]*(em-ab);ct=fac[:,None]*(et*(1+J)-at*J)/self.scale;cy=fac[:,None]*(ey*(1+J)-ay*J)/self.scale
        cw=ct*(self.uscale/ut)[:,None];cy-=ct*(uy/ut)[:,None]
        add(I,I,cj);add(I,self.it[:,None],cw);add(I,self.iy[:,None],cy)
        for index,weight in [(self.it,self.energy),(self.iy,self.num)]:
            add(index[:,None],I,weight*cj);add(index,self.it,(weight*cw).sum(1));add(index,self.iy,(weight*cy).sum(1))
        add(F[1:-1],F[1:-1],-rate)
        ot=d['rho'][:,None]*d['thermo'][:,4,None]*(at-et)+net[:,None]*6.6524587321e-25
        oy=d['rho'][:,None]*d['thermo'][:,4,None]*(ay-ey)+ney[:,None]*6.6524587321e-25
        ow=ot*(self.uscale/ut)[:,None];oy-=ot*(uy/ut)[:,None]
        for index,der in [(self.it,ow),(self.iy,oy)]:
            for shift in [0,1]:
                sel=slice(shift,n-1+shift);v=-d['face_a'][1:-1,None]*C*der[sel]/2*HH[1:-1]
                add(F[1:-1],index[sel,None],v)
        jac=self.transport+sparse.coo_matrix((vals,(rows,cols)),shape=(self.size,self.size)).tocsc()
        return base,theta,p,u,jac

    def run(self,steps,label):
        assert not (OUT/(label+'.npz')).exists();start=time.monotonic();h=prior.END/steps;state=np.zeros(self.size);history=[];flux_energy=0.;balance=0.;max_newton=0
        d=self.d;port=4*np.pi*C*d['edges']**2/d['face_a']**2;q=d['num']*d['Einf']*self.scale
        scale_state=np.r_[np.maximum((self.J/self.scale).ravel(),1e-8),np.maximum((self.H/self.scale).ravel(),1e-8),np.ones(2*self.n)]
        def power(s):return float(port[0]*(q@(s[self.inner]+self.H[0]/self.scale))-port[-1]*(q@(s[self.outer]+self.H[-1]/self.scale)))
        def energy(s):return float(self.energyJ@s[:self.nj]+(self.volume*d['rho']*d['a']*self.uscale)@s[self.it])
        def record(s):
            _,t,p,u=self.rhs(s,False);tr=d['rho']*(u-self.u0)-3*(p-self.eos.gas(np.zeros(self.n),np.zeros(self.n))[0])
            history.append((t.copy(),s[self.iy].copy(),tr))
        record(state);snapshots=[state.copy()];failed=None
        try:
            for step in range(steps):
                f0=self.rhs(state,False)[0];new=state.copy()
                for it in range(10):
                    rhs,theta,p,u,jac=self.rhs(new);res=new-state-h*(f0+rhs)/2;res[self.outer]=new[self.outer]-.5*new[self.I[-1]]
                    err=float(max(abs(res)/scale_state))
                    if err<1e-9:break
                    left=(sparse.eye(self.size,format='csc')-h*jac/2).tolil()
                    for f,j in zip(self.outer,self.I[-1]):left.rows[f]=[int(j),int(f)];left.data[f]=[-.5,1.]
                    delta=splu(left.tocsc()).solve(-res)
                    # Native table bounds govern Newton trials as well as
                    # accepted states. Backtracking never clips a population.
                    for cut in range(12):
                        trial=new+delta*2.**(-cut)
                        try:
                            rr=self.rhs(trial,False)[0];check=trial-state-h*(f0+rr)/2;check[self.outer]=trial[self.outer]-.5*trial[self.I[-1]]
                        except AssertionError:continue
                        if max(abs(check)/scale_state)<err:break
                    else:raise AssertionError('Nonlinear line search')
                    new=trial
                else:raise AssertionError('Newton iteration cap')
                max_newton=max(max_newton,it);flux_energy+=h*(power(state)+power(new))/2
                balance=max(balance,abs(energy(new)-flux_energy));state=new;record(state);snapshots.append(state.copy())
        except Exception as exc:failed=repr(exc)
        hist=np.array(history);trace=hist[:,2]@(self.volume*d['a']);ss=np.array(snapshots);photon=ss[:,:self.nj].reshape(-1,self.n,self.m)*self.scale+self.J
        scale=max(abs(flux_energy),max(abs(trace)),1.)
        np.savez_compressed(OUT/(label+'.npz'),times=np.arange(len(history))*h,theta=hist[:,0],eta=hist[:,1],trace_density=hist[:,2],trace_energy=trace,state=ss,J=photon)
        result=dict(classification='Counterexample candidate',passed=bool(failed is None and balance/scale<1e-8 and photon.min()>=-1e-10*self.J.max()),failure=failed,
            steps=steps,completed_steps=len(history)-1,seconds=time.monotonic()-start,energy_balance_relative=float(balance/scale),maximum_Newton_iterations=max_newton,
            maximum_logT_change=float(max(abs(hist[:,0]).ravel())),maximum_relative_neutral_change=float(max(abs(hist[:,1]).ravel())),minimum_photon_occupation=float(photon.min()),trace_endpoint_erg=float(trace[-1]),
            nonlinear_native_chemistry_and_heat_evolved=failed is None,spatial_convergence=False,full_GR_feedback=False,final_charge_solved=False,source_sha256=prior.sha(__file__))
        write(OUT/(label+'.json'),result);print(json.dumps(result),flush=True);return result


def run():
    assert json.loads((OUT/'repaired-bank.json').read_text())['passed'];signal.alarm(120);start=time.monotonic();m=Model();rows=[]
    for n in [64,128]:
        row=m.run(n,f'steps-{n}');rows.append(row)
        if not row['passed']:break
    result=dict(classification='Counterexample candidate',passed=False,paths=rows,seconds=time.monotonic()-start)
    if len(rows)==2 and all(r['passed'] for r in rows):
        a=np.load(OUT/'steps-64.npz');b=np.load(OUT/'steps-128.npz');err=float(np.max(abs(a['trace_energy']-b['trace_energy'][::2]))/np.max(abs(b['trace_energy'])))
        result.update(passed=err<.02,time_trace_relative=err)
    write(OUT/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


def repair():
    assert not (OUT/'repaired-bank.json').exists();start=time.monotonic();signal.alarm(20)
    failed=json.loads((OUT/'bank.json').read_text());assert not failed['passed']
    write(OUT/'interpolation-reassessment.json',dict(classification='Counterexample candidate',
        failure='Maximum frequencywise reverse-rate interpolation error0.004064 exceeds the original0.002 gate. The exact exp(-photon_energy/kT) was interpolated in log-temperature.',
        repair='Factor the known Boltzmann exponential analytically and interpolate only its prefactor. Reuse all96 native states and keep all physical supports and gates. No extra table node, frequency or fluid cell.',
        remaining_native_calls=1200-failed['native_calls'],remaining_native_seconds=45-failed['seconds'],source_sha256=prior.sha(__file__)))
    d=np.load(prior.OUT/'bank-16-8.npz');native=old.Native(cap=1200-failed['native_calls']);eos=Table();checks=[];saved=[]
    for j in range(len(d['r'])):
        setup(native,d,j);theta=.017 if j%2 else -.023;eta=.4 if j%2 else -.15
        s=native.state(0.,np.log(d['T'][j])+theta,native.y0*(1+eta));k,e=prior.coefficients(native,s,d['Einf']/d['a'][j]);a=k+e
        t=np.full(len(d['r']),theta);y=np.full(len(d['r']),eta);p,u,*_=eos.gas(t,y);aa,ee,*_=eos.radiation(t,y);mask=(a>0)&(e>0)
        checks.append(dict(cell=j,constitutive=float(max(abs(p[j]/s['raw'][1]-1),abs(u[j]/s['raw'][2]-1))),rate=float(max(np.max(abs(aa[j,mask]/a[mask]-1)),np.max(abs(ee[j,mask]/e[mask]-1))))))
        saved.append(dict(theta=theta,eta=eta,raw=s['raw'],absorption=a,emission=e))
    np.savez_compressed(OUT/'native-controls.npz',**{k:np.array([a[k] for a in saved]) for k in saved[0]},native_calls=native.ion.calls)
    result=dict(classification='Counterexample candidate',passed=bool(max(a['rate'] for a in checks)<.002 and max(a['constitutive'] for a in checks)<.002),checks=checks,
        new_native_calls=native.ion.calls,total_native_calls=failed['native_calls']+native.ion.calls,seconds=time.monotonic()-start,source_sha256=prior.sha(__file__))
    write(OUT/'repaired-bank.json',result);signal.alarm(0);print(json.dumps(result),flush=True);assert result['passed']


if __name__=='__main__':globals()[sys.argv[1]]()
