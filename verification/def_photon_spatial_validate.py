"""Bounded validation and stored-source sensitivity of the unsplit photon tangent."""
from pathlib import Path
import argparse
import json
import resource
import signal
import time
import numpy as np
import def_photon_spatial_coupling as m

OUT=m.OUT


def prepare():
    assert not (OUT/'coupling-plan.json').exists()
    paths=[Path(__file__),Path(m.__file__),OUT/'bank.npz',OUT/'bank.json',OUT/'reassessment.json',OUT/'pilot.json']
    m.ex.write(OUT/'coupling-plan.json',dict(classification='Counterexample candidate',
        claim='Evolve the same material plus spectral/angular photons with absorption, dipole angular gain, stimulated Kompaneets exchange and spatial streaming in one unsplit system.',
        duration='One local temperature scale-height light travel time for kH=0,1,10. The isolated Kompaneets control runs for 8/rate_C with coefficients frozen; it is not stellar evolution over that physical duration.',
        angular_modes=[16,32],steps=[16,32,64],
        choice_before_production='The kH=10 mode crosses 10 angular phase radians. Use 16/32 Legendre modes rather than the pilot six modes, which measured cost only. No result-driven angular enlargement.',
        gates=dict(time_order_min=1.8,time_difference_initial=.001,angle_difference_initial=.001,
            energy_balance=1e-9,energy_equation_residual=1e-9,entropy_growth=1e-10,tail_capacity_fraction=1e-17,
            Kompaneets_heating_relative=1e-5,Kompaneets_number=1e-10,equilibrium_projection=.001),
        controls=['k=0, Compton-off reproduces Phase60 arrowhead solver on the same cells',
            'Isolated Kompaneets: photon number, energy, LTE/Bose null modes, continuum 4*a*T^3*rate_C heating, 16/32/64 time order and asymptotic energy/number projection',
            'No-collision Fourier streaming: exact angular-average sinc(k*c*t) and dipole/quadrupole moments',
            'Split every fixed-opacity frequency cell in two; compare the coupled energy/temperature and isolated Kompaneets moments'],
        uncertainty='At actual matter state repeat kH=1 with each bracket-temperature opacity held at the same physical frequency, and with isotropic rather than dipole gain and Compton removed. These are model sensitivities, not certified envelopes or permission to replace an incomplete physical kernel.',
        budget=dict(measured_pilot_seconds=json.loads((OUT/'pilot.json').read_text())['seconds'],
            forecast_seconds=[15,80],forecast_basis='16/32 angular scaling and step count from the six-mode pilot; sparse factorization scaling is an assumption',
            hard_seconds=180,peak_memory_estimate_GB=[.3,1.5],CPU_workers=1,new_queries=0,new_EOS_calls=0,new_stellar_steps=0,automatic_expansion=False),
        bindings={p.relative_to(m.h.ROOT).as_posix():m.h.digest(p) for p in paths}))


def summary(r):return {k:v for k,v in r.items() if k not in ['T','E']}


def run():
    assert not (OUT/'result.json').exists();start=time.monotonic();signal.alarm(180)
    plan=json.loads((OUT/'coupling-plan.json').read_text())
    for p,sha in plan['bindings'].items():assert m.h.digest(m.h.ROOT/p)==sha,p
    b=dict(np.load(OUT/'bank.npz'));duration=float(b['H']/b['c']);rows=[];paths={}
    for kH in [0,1,10]:
        op=m.Operator(b,16,kH);runs=[op.evolve(duration,n) for n in plan['steps']]
        errors=[m.compare(a,z,op) for a,z in zip(runs[:-1],runs[1:])]
        order=float(np.log2(errors[0]/errors[1]));small=runs[-1]
        del op
        op=m.Operator(b,32,kH);fine=op.evolve(duration,64)
        # Orthogonal projection comparison counts the high angular modes too.
        common=small['E'].reshape(-1,16);large=fine['E'].reshape(-1,32)
        angular=float(np.sqrt(abs(small['T']-fine['T'])**2+np.linalg.norm(common-large[:,:16])**2+np.linalg.norm(large[:,16:])**2)/np.sqrt(op.Cm))
        row=dict(kH=kH,duration_seconds=duration,time_order=order,time_differences_initial=errors,
            angular_difference_initial=angular,tail_capacity_fraction=op.tail,
            final_temperature_K=float((fine['T']/np.sqrt(op.Cm)).real),
            energy_balance=max(r['balance'] for r in runs+[fine]),
            energy_equation_residual=max(r['energy_equation_residual'] for r in runs+[fine]),
            entropy_growth=max(r['entropy_growth'] for r in runs+[fine]))
        row['passed']=bool(order>=1.8 and errors[-1]<.001 and angular<.001 and op.tail<1e-17
            and row['energy_balance']<1e-9 and row['energy_equation_residual']<1e-9 and row['entropy_growth']<1e-10)
        np.savez_compressed(OUT/f'mode-{kH}.npz',T=fine['T'],E=fine['E'].reshape(-1,32),u=op.u,Ci=op.Ci,Cm=op.Cm,kc=op.kc,duration=duration)
        rows.append(row);paths[kH]=small
        print('MODE',row,flush=True);del op
    # Verify the actual evolved material force, not just a temperature readout.
    op=m.Operator(b,16,1);r=paths[1];T,E=r['T'],r['E'];Td,Ed=op.rhs(T,E)
    modes=E.reshape(-1,16);dmodes=Ed.reshape(-1,16);root=np.sqrt(op.Ci)
    E0=float(np.sum(root*modes[:,0]).real);F=float(b['c'])*(root@modes[:,1])/np.sqrt(3)
    pressure=root@(modes[:,0]/3+2*modes[:,2]/(3*np.sqrt(5)))
    force=(root@((op.rate_a+op.rate_s)*modes[:,1]))/(np.sqrt(3)*float(b['c']))
    momentum=(root@dmodes[:,1])/(np.sqrt(3)*float(b['c']))+1j*op.kc/float(b['c'])*pressure+force
    momentum_score=float(abs(momentum)/(1+abs(force)+abs(op.kc/float(b['c'])*pressure)))
    m.ex.write(OUT/'stress-source.json',dict(classification='Counterexample candidate',
        E0_erg_cm3=E0,flux_erg_cm2_s=[float(F.real),float(F.imag)],
        radial_pressure_erg_cm3=[float(pressure.real),float(pressure.imag)],
        force_on_matter_erg_cm4=[float(force.real),float(force.imag)],momentum_balance_score=momentum_score,
        matter_velocity_evolved=False,interpretation='Complex amplitudes of one local spatial Fourier mode. Force is the opposite photon collisional momentum rate; no global atmosphere or GR momentum trajectory is claimed.'))
    # Absorption-only homogeneous regression against the separate arrowhead.
    zero=m.Operator(b,2,0,compton=False);test=zero.evolve(duration,64)
    ref=m.previous.evolve(zero.Cm,zero.Ci,zero.rate_a,duration,64)
    e=test['E'].reshape(-1,2)[:,0]*np.sqrt(zero.Ci)
    arrow=float(np.real(m.previous.norm(test['T']/np.sqrt(zero.Cm)-ref[0],e-ref[1],zero.Cm,zero.Ci)/np.sqrt(zero.Cm)))
    # Free angular propagation has a continuum sinc energy response.
    from scipy.linalg import expm
    from scipy.special import spherical_jn
    free=[]
    for n in [16,32]:
        ell=np.arange(n-1);a=(ell+1)/np.sqrt((2*ell+1)*(2*ell+3));V=np.diag(a,1)+np.diag(a,-1)
        initial=np.eye(n)[:,0];actual=expm(-10j*V)@initial
        exact=np.array([np.sqrt(2*l+1)*(-1j)**l*spherical_jn(l,10) for l in range(3)])
        free.append(float(np.max(abs(actual[:3]-exact))))
    # Standalone Compton keeps the scalar opacity data out of the control.
    c=dict(b,rate_a=np.zeros_like(b['rate_a']),rate_s=np.zeros_like(b['rate_s']))
    op=m.Operator(c,1);cruns=[op.evolve(8/float(b['rate_C']),n) for n in [16,32,64]]
    ce=[m.compare(a,z,op) for a,z in zip(cruns[:-1],cruns[1:])];cp=float(np.log2(ce[0]/ce[1]))
    number=float(abs(op.number@cruns[-1]['E'])/(np.sqrt(op.Cm)*np.linalg.norm(op.number)))
    initial_T=np.sqrt(op.Cm);initial_E=np.zeros(op.size,complex);Td,Ed=op.rhs(initial_T,initial_E)
    heating=float(op.energy@Ed.real);exact_heat=4*float(b['arad'])*float(b['T'])**3*float(b['rate_C'])
    heating_error=abs(heating/exact_heat-1)
    v1=np.r_[np.sqrt(op.Cm),op.energy];v2=np.r_[0,op.number];Q=np.column_stack([v1,v2]);initial=np.r_[initial_T,initial_E]
    equilibrium=Q@np.linalg.solve(Q.T@Q,Q.T@initial)
    last=np.r_[cruns[-1]['T'],cruns[-1]['E']];equilibrium_error=float(np.linalg.norm(last-equilibrium)/initial_T)
    null=[]
    for vec in [v1,v2]:
        a,z=op.rhs(vec[0],vec[1:]);scale=op.aa*abs(vec[0])+np.linalg.norm(op.L@vec[1:])+np.linalg.norm(op.q)*np.linalg.norm(vec)
        null.append(float(np.sqrt(abs(a)**2+np.linalg.norm(z)**2)/scale))
    # Refine frequency transport on the same fixed opacity field, not a new fit.
    split_edges=np.sort(np.r_[b['edges_u'],b['u']]);newu=(split_edges[:-1]+split_edges[1:])/2
    newC=4*float(b['arad'])*float(b['T'])**3*m.previous.loss.weights(split_edges,1)[:,1]
    finebank=dict(b,edges_u=split_edges,u=newu,Ci=newC,rate_a=np.repeat(b['rate_a'],2),rate_s=np.repeat(b['rate_s'],2))
    fineop=m.Operator(finebank,16,1);fr=fineop.evolve(duration,64)
    # Compare integrated photon moments on the common parent cells, allowing
    # the small final partial cell at the declared u=60 truncation boundary.
    nr=min(len(paths[1]['E'])//16,len(fr['E'])//32)
    refE=paths[1]['E'].reshape(-1,16)[:nr]*np.sqrt(b['Ci'][:nr,None])
    fineE=(fr['E'].reshape(-1,16)[:2*nr]*np.sqrt(fineop.Ci[:2*nr,None])).reshape(nr,2,16).sum(1)
    freq=float(np.sqrt(abs(fr['T']-paths[1]['T'])**2+np.sum(abs(fineE-refE)**2/b['Ci'][:nr,None]))/np.sqrt(op.Cm))
    controls=dict(arrowhead_relative=arrow,free_angular_16_32=free,Kompaneets_time_order=cp,
        Kompaneets_time_differences_initial=ce,Kompaneets_number_relative=number,
        Kompaneets_energy_balance=max(r['balance'] for r in cruns),Kompaneets_heating_relative=float(heating_error),
        Kompaneets_equilibrium_projection=equilibrium_error,Kompaneets_null_residuals=null,momentum_balance_score=momentum_score,
        frequency_split_energy_norm_difference_initial=freq)
    controls['passed']=bool(arrow<1e-9 and max(free)<1e-9 and cp>=1.8 and ce[-1]<.001 and number<1e-10
        and controls['Kompaneets_energy_balance']<1e-9 and heating_error<1e-5 and equilibrium_error<.001
        and max(null)<1e-10 and momentum_score<1e-9 and freq<.001)
    np.savez_compressed(OUT/'compton-control.npz',T=cruns[-1]['T'],E=cruns[-1]['E'],equilibrium=equilibrium,Ci=op.Ci,Cm=op.Cm,u=op.u)
    m.ex.write(OUT/'controls.json',controls);print('CONTROLS',controls,flush=True);del op
    sensitivities=[]
    for name,alt,options in [
        ('lower-temperature-opacity',dict(b,rate_a=b['rate_a']*b['absorption_source_states'][0]/b['absorption'],rate_s=b['rate_s']*b['scattering_source_states'][0]/b['scattering']),{}),
        ('upper-temperature-opacity',dict(b,rate_a=b['rate_a']*b['absorption_source_states'][1]/b['absorption'],rate_s=b['rate_s']*b['scattering_source_states'][1]/b['scattering']),{}),
        ('isotropic-gain',b,dict(dipole=False)),('no-Compton',b,dict(compton=False))]:
        op=m.Operator(alt,16,1,**options);r=op.evolve(duration,64)
        sensitivities.append(dict(name=name,energy_norm_difference_initial=m.compare(r,paths[1],op),temperature_difference_K=float(((r['T']-paths[1]['T'])/np.sqrt(op.Cm)).real)))
        del op
    result=dict(classification='Counterexample candidate',passed=all(r['passed'] for r in rows) and controls['passed'],
        rows=rows,controls=controls,sensitivities=sensitivities,seconds=time.monotonic()-start,
        peak_memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
        unsplit_absorption_scattering_Compton_spatial_candidate_evolved=True,
        actual_temperature_native_request_passed=False,actual_opacity_interpolation_certified=False,
        native_angular_energy_kernel_certified=False,physical_atmosphere_closed=False,
        full_GR_photon_feedback_evolved=False,full_dynamic_charge_solved=False)
    m.ex.write(OUT/'result.json',result);print('RESULT',result,flush=True);signal.alarm(0)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
