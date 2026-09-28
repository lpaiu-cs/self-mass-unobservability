"""Put saved native atmospheric states through the actual collision owner."""
from pathlib import Path
import json,signal,time
import numpy as np
import def_native_atmosphere_inverse as audit

run=audit.run;OUT=audit.OUT;C=run.C;write=run.write;sha=run.sha


def main():
    assert not (OUT/'collision-result.json').exists()
    write(OUT/'collision-plan.json',dict(classification='Counterexample candidate',
        claim='Use all223 conserved native fine endpoint states with the actual angular photons and atmospheric moving collision owner. Compare positive channels and net H,energy,momentum transfers, including the changed primitive and electrons.',
        reuse='Saved native roots, fractions and affinities; no new native state evaluations or time steps. Hydro is zeroed for both evaluations to isolate the implemented local collision exchange.',
        gates=dict(positive_channels=.002,net_exchange=.02,paired_energy=1e-10),
        budget=dict(seconds=40,CPU_threads=1,memory_GB=3,new_native_states=0),
        limits='Actual fine endpoint only. No complete atmospheric history, spectral continuum, derivative or feedback certificate.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(audit.__file__),OUT/'production-samples.json',OUT/'production.json',run.EV/'final-128.npz']}))
    start=time.monotonic();signal.signal(signal.SIGALRM,run.flow.old.optical.timeout);signal.alarm(40)
    model=run.Coupled();f=model.flow;m=model.m;z=np.load(run.EV/'final-128.npz');U=z['U'];I=z['I'].sum(0)
    V=f.primitive(U);kap=f.eos(V[0],V[2])[4];ids=np.flatnonzero(U[0]>=f.eos.floor)
    rows={r['cell']:r for r in json.loads((OUT/'production-samples.json').read_text()) if r['steps']==128}
    assert set(rows)==set(map(int,ids));newV=V.copy();newkap=kap.copy()
    for j in ids:
        r=rows[int(j)];newV[:,j]=[r['rho'],r['v'],r['lt'],r['y']]
        newkap[j]=r['raw'][13]/r['raw'][0]*6.6524587321e-25/1.66053906660e-24
    hydro=f.hydro;coeff=model.spectrum.coefficients;cross=model.spectrum.cross
    fraction=np.array([rows[int(j)]['fraction'] for j in ids]);affinity=np.array([rows[int(j)]['affinity'] for j in ids])
    def native_coeff(rho,lt,y,energy):
        assert len(rho)==len(ids)
        ab=np.zeros_like(energy);em=np.zeros_like(energy);T=np.exp(lt)[:,None,None]
        for k in range(10):
            sigma=cross(energy,k+1);pop=y*fraction[:,k]
            ab+=pop[:,None,None]*sigma
            good=pop>0;logpop=np.full_like(pop,-np.inf);logpop[good]=np.log(pop[good])
            em+=np.exp(logpop[:,None,None]+affinity[:,None,None]-energy/(run.chem.old.K*T))*sigma
        return ab,em
    def evaluate(primitive,opacity,coefficient):
        f.hydro=lambda state,t:(np.zeros_like(state),np.zeros(3),1e100,primitive,opacity)
        model.spectrum.coefficients=coefficient
        return model.local_rhs(U,I,run.flow.old.END)
    try:before=evaluate(V,kap,coeff);after=evaluate(newV,newkap,native_coeff)
    finally:f.hydro=hydro;model.spectrum.coefficients=coeff
    volume=4*np.pi*m.RJ*m.RJ*m.vol;unit=f.eos.rho0*C*C*volume
    net=[]
    for k,label in [(3,'neutral_H'),(2,'Killing_energy'),(1,'radial_momentum')]:
        weight=volume*f.eos.rho0*(f.eos.nH if k==3 else C*C)
        old=before[0][k]*weight;new=after[0][k]*weight
        net.append(dict(quantity=label,relative=float(np.sum(abs(new-old))/max(np.sum(abs(old)),1.))))
    energy=[]
    for label,rhs in [('table',before),('native',after)]:
        photon=np.sum(rhs[1]*model.energy_weight,dtype=np.longdouble);gas=np.sum(rhs[0][2]*unit,dtype=np.longdouble);escape=rhs[2][-1]
        energy.append(dict(owner=label,relative=float(abs(photon+gas+escape)/max(abs(photon),abs(gas),1.))))
    channels=[];a=m.a[ids];field=I[ids];number=volume[ids,None,None]*model.w[None,:,None]*model.number/a[:,None,None]**3
    # Both complete physical primitives are retained in the rate comparison.
    def positive(primitive,fn):
        rho,v,lt,y=primitive[:,ids];D=(1-v[:,None]*model.mu[None,:])/np.sqrt(1-v*v)[:,None]
        E=model.E[None,None,:]*D[:,:,None]/a[:,None,None];ab,em=fn(rho*f.eos.rho0,lt,y,E)
        factor=a[:,None,None]*C*(rho*f.eos.rho0*f.eos.nH)[:,None,None]*D[:,:,None]
        return [factor*ab*field,factor*em,factor*em*field]
    old=positive(V,coeff);new=positive(newV,native_coeff)
    for label,o,n in zip(['absorption','spontaneous','stimulated'],old,new):
        channels.append(dict(channel=label,number_relative=float(np.sum(abs(n-o)*number)/max(np.sum(n*number),1.)),
            energy_relative=float(np.sum(abs(n-o)*number*model.E)/max(np.sum(n*number*model.E),1.))))
    passed=max(r['relative'] for r in net)<.02 and max(r['relative'] for r in energy)<1e-10 and max(max(r['number_relative'],r['energy_relative']) for r in channels)<.002
    result=dict(classification='Counterexample candidate',passed=bool(passed),native_cells=len(ids),net_exchange=net,
        positive_channels=channels,paired_energy=energy,seconds=time.monotonic()-start,new_native_state_evaluations=0,
        actual_collision_owner_used=True,full_atmosphere_history_audited=False,actual_atmosphere_EOS_re_evolution=False,final_charge_solved=False)
    write(OUT/'collision-result.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
