"""Native EOS audit and a usable core/envelope background with one base flux."""
from pathlib import Path
import json
import signal
import time
import numpy as np
import sympy as sp
from scipy.integrate import simpson
import def_native_whole_star_match as task


def main():
    out=task.OUT
    assert not (out/'audit.json').exists()
    result=json.loads((out/'final-result.json').read_text());assert result['passed']
    remaining=int(600-result['total_compute_seconds'])
    assert remaining>=30,'Audit forecast no longer fits remaining original budget'
    task.write(out/'audit-plan.json',dict(classification='Counterexample candidate',
        claim='Check native core thermodynamics, independently integrate envelope baryons and luminosity, retain isotope inventory and construct the full finite-pressure background/lapse. Put one luminosity on the shared material face and retain the resulting core heat imbalance.',
        gates=dict(native_logrho=1e-7,native_entropy_over_cv=1e-7,envelope_mass_relative=1e-5,
            luminosity_relative=1e-5,isotope_inventory_relative=1e-13),
        budget=dict(hard_seconds=remaining,native_calls=430,new_roots=0,new_evolution_paths=0,
            forecast_seconds_range=[15,30],total_cap_seconds=600),
        bindings={str(p.relative_to(task.ROOT)):task.envelope.digest(p) for p in
            [Path(__file__),Path(task.__file__),out/'final-result.json',out/'final-core.npz',out/'final-envelope.npz']}))
    start=time.monotonic();signal.alarm(remaining)
    Li,Lb,Lo,Qc,Qe=sp.symbols('Li Lb Lo Qc Qe')
    assert sp.expand((Li-Lb+Qc)+(Lb-Lo+Qe)-(Li-Lo+Qc+Qe))==0
    core=np.load(out/'final-core.npz');env=np.load(out/'final-envelope.npz')
    model=task.Match();s=model.core;p=model.env;eos=p.eos
    original=np.load(task.envelope.prior.OUT/'coefficients.npz')
    states=core['states'].copy();thermo=core['thermo'].copy()
    # The junction is a cell boundary; give the retained half of cell132 its
    # actual midpoint instead of treating the saved boundary as a cell centre.
    first=model.native_fraction;end=s.outer[model.i+1]
    states[0]=s.step(np.log(first),np.log((first+end)/2),core['outer'][0,1:],model.i,s.B,True)
    thermo[0]=s.state(states[0,2],model.i)[3]
    ids=np.unique(np.r_[0,1,2,3,10,30,100,np.linspace(200,len(states)-1,12).astype(int)])
    audits=[];native={}
    for j in ids:
        i=j+model.i;raw=eos(1,float(states[j,2]),float(thermo[j,1]),s.data['X'][i]);native[int(j)]=raw
        cv=(raw[10]-raw[9]*raw[8]/raw[7])/np.exp(thermo[j,1])
        err=[abs(np.log(raw[0])-thermo[j,0]),abs(raw[3]-s.ref[i,3])/max(abs(cv),1)]
        assert max(err)<1e-7,(i,err)
        audits.append(dict(cell=int(i),logrho_error=float(err[0]),entropy_over_cv=float(err[1])))
    # Integrate the exact cellwise isentropic first integral from the newly
    # solved envelope lapse, keeping lapse continuous at composition jumps.
    faces=np.r_[core['outer'][:,1:],core['inner'][-2:0:-1,1:],
        [[0.,0.,core['parameters'][0],core['parameters'][3],0.]]]
    assert len(faces)==len(states)+1,(len(faces),len(states))
    def phi(y):return .001*(1+s.mu*y[3])
    def enthalpy(y,i):
        v=s.state(y[2],i)[3]
        return np.longdouble(s.data['CX'][i])+np.longdouble(v[2])+np.exp(np.longdouble(y[2])-np.longdouble(v[0]))/np.longdouble(p.c)**2
    nu=float(np.log(env['N'][0]));nu_faces=[nu];nu_mid=[];gammas=[]
    for j,i in enumerate(range(model.i,len(s.f))):
        a,b=faces[j],faces[j+1];mid=states[j]
        ha,hb,hm=enthalpy(a,i),enthalpy(b,i),enthalpy(mid,i)
        nu_mid.append(nu+2*(phi(mid)**2-phi(a)**2)-float(np.log1p((hm-ha)/ha)))
        nu+=2*(phi(b)**2-phi(a)**2)-float(np.log1p((hb-ha)/ha))
        nu_faces.append(nu)
        # Derivative of the same bounded EOS interpolant used in the solve.
        dx=mid[2]-s.lp[i]
        xx,co=s.extra[i] if i in s.extra else (s.offset,s.coef[:,:,i])
        k=min(len(xx)-2,max(0,int(np.searchsorted(xx,dx)-1)));t=dx-xx[k]
        derivative=(3*co[0,k]*t+2*co[1,k])*t+co[2,k]
        gammas.append(1/derivative[0])
    gamma_error=max(abs(gammas[j]/raw[4]-1) for j,raw in native.items())
    rho=np.exp(thermo[:,0]);pressure=np.exp(states[:,2]);T=np.exp(thermo[:,1]);cx=s.data['CX'][model.i:]
    energy=rho*p.c**2*(cx+thermo[:,2]);r=states[:,0]*s.R*100;m=states[:,1]*s.B*100
    ph=.001*(1+s.mu*states[:,3]);v=.001*s.mu*states[:,4]/(s.R*100)
    A=np.exp(-2*ph**2);N=np.exp(nu_mid)
    raws=[]
    for P,temperature in zip(env['Ptotal'],env['T']):raws.append(eos(1,float(np.log(P)),float(np.log(temperature)),p.X))
    raws=np.array(raws)
    mass=simpson(4*np.pi*env['r']**2*env['A']**3*env['rho']/np.sqrt(env['b']),x=env['r'])/p.total_baryon
    mass_error=abs(mass/model.target-1);assert mass_error<1e-5,mass_error
    theta=env['A']*env['N']*env['T'];integral=simpson(env['optical']*env['N']**2/env['r']**2,x=env['r'])
    L=4*np.pi*p.arad*p.c/3*(theta[0]**4-theta[-1]**4)/integral
    luminosity_error=abs(L/float(env['Linfinity'])-1);assert luminosity_error<1e-5
    # Baryon and isotope equality is assessed on the tiny replaced component,
    # as well as the whole star, so total-mass rounding cannot hide a mismatch.
    old_isotopes=np.sum(s.data['dm'][:,None]*s.data['X'],axis=0)+model.old_atmosphere_fraction*p.total_baryon*p.X
    new_isotopes=np.sum(core['dm'][:,None]*core['X'],axis=0)+model.target*p.total_baryon*p.X
    isotope_error=float(np.max(abs(new_isotopes-old_isotopes)/np.maximum(abs(old_isotopes),1)))
    assert isotope_error<1e-13,isotope_error
    oldatm=json.loads((out.parent/'def-native-atmosphere/density-coordinate/result.json').read_text())['rows'][-1]
    oldbound=json.loads((task.envelope.prior.surface.OUT/'surface-enclosure-corrected/result.json').read_text())
    # Same hydrostatic pressure-column bound applied above the last native
    # atmosphere point; the unresolved cold branch premise remains explicit.
    baseP=json.loads((out.parent/'def-native-atmosphere/density-coordinate/table.json').read_text())['rows'][0]['raw'][1]
    tail_bound=oldbound['total_added_baryon_fraction_upper']*oldatm['endpoint_pressure_dyn_cm2']/baseP
    # Shared base flux is the envelope's L. Adjacent core radiative face flux
    # is computed at its actual native temperatures; their difference is heat.
    kap=[task.envelope.prior.two.opacity_parts(p.opacity,(float(np.log(rho[j])),float(np.log(T[j])),s.data['X'][model.i+j]))[0] for j in [0,1]]
    K=4*p.arad*p.c*T[:2]**3/(3*rho[:2]*kap)
    face=faces[1];rf=face[0]*s.R*100;mf=face[1]*s.B*100;af=np.exp(-2*phi(face)**2);nf=np.exp(nu_faces[1])
    flux=-4*np.pi*rf*rf*nf*af*af*np.sqrt(K[0]*K[1])*np.sqrt(1-2*mf/rf)*(A[1]*N[1]*T[1]-A[0]*N[0]*T[0])/(r[1]-r[0])
    baseL=float(env['Linfinity']);heat=float((flux-baseL)/(core['dm'][0]*A[0]*N[0]))
    cv=native[0][10]-native[0][9]*native[0][8]/native[0][7]
    temperature_rate=heat/cv
    oldcenter=.001*(1+s.mu*core['parameters'][3])
    center=s.state(core['parameters'][0],len(s.f)-1)
    def merged(a,b,central):return np.r_[central,a[::-1],b]
    np.savez_compressed(out/'background.npz',
        radius_cm=merged(r,env['r'],0),mass_geom_cm=merged(m,env['m'],0),
        phi=merged(ph,env['phi'],oldcenter),phi_prime_cm=merged(v,env['v'],0),
        lapse=merged(N,env['N'],np.exp(nu_faces[-1])),
        pressure_cgs=merged(pressure,env['Ptotal'],np.exp(core['parameters'][0])),
        energy_cgs=merged(energy,env['energy'],center[1]*p.c**4/(p.G*1e4)),
        density_cgs=merged(rho,env['rho'],np.exp(center[3][0])),
        temperature_K=merged(T,env['T'],np.exp(center[3][1])),
        gamma1=merged(np.array(gammas),raws[:,4],native[len(states)-1][4]),
        core_faces=faces,core_logN_faces=nu_faces,core_logN_mid=nu_mid,
        core_states=states,core_thermo=thermo,envelope_native_raw=raws,
        junction_index=len(states)+1,base_luminosity=baseL,core_neighbour_face_luminosity=flux)
    seconds=time.monotonic()-start
    audit=dict(classification='Counterexample candidate',passed=True,native_core_controls=audits,
        symbolic=dict(classification='Proven',passed=True,identity='Shared junction luminosity cancels in the sum of adjacent material energy equations; no zero-volume interface heat reservoir is introduced.'),
        envelope_baryon_integral_relative=mass_error,luminosity_integral_relative=luminosity_error,
        whole_isotope_relative=isotope_error,core_gamma_interpolant_native_relative=float(gamma_error),
        omitted_old_cold_tail_baryon_fraction_upper=tail_bound,
        base_luminosity=baseL,first_core_internal_face_luminosity=float(flux),
        old_reference_face_luminosity=float(original['luminosity'][model.i,0]),
        first_core_retained_half_cell_heat_erg_g_s=heat,
        frozen_density_logT_rate_s=float(temperature_rate),
        frozen_density_one_percent_timescale_s=float(.01/abs(temperature_rate)),
        thermal_interpretation='One luminosity is shared on the actual base face; the neighbouring core flux differs and gives retained-cell heat. This is a nonstationary initial background, not a steady core-envelope flux match or an evolved timescale.',
        native_calls=len(raws)+len(ids),seconds=seconds,total_compute_seconds=result['total_compute_seconds']+seconds,
        same_material_mechanical_background=True,steady_thermal_match=False,physical_zero_pressure_surface=False,
        final_dynamic_charge_solved=False,full_goal_complete=False)
    assert audit['native_calls']<=430 and audit['total_compute_seconds']<600
    task.write(out/'audit.json',audit);signal.alarm(0);print(json.dumps(audit),flush=True)


if __name__=='__main__':main()
