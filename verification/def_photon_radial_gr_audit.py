"""Saved-state physical readout; no new evolution and no native EOS calls."""
from pathlib import Path
import hashlib
import json
import time
import sys
import signal
import numpy as np

OUT=Path('outputs/direct-eos-gr33/def-photon-radial-gr')
read=lambda p:json.loads(p.read_text())
write=lambda p,x:p.write_text(json.dumps(x,indent=2)+'\n')


def state_check():
    import def_photon_radial_gr as model
    start=time.monotonic();th=np.load(OUT.parent/'def-photon-population-coupled/thermo.npz');rows={}
    for label in ['fine-32','fine-64','fine-128','spatial-p2']:
        d=np.load(OUT/(label+'.npz'));change=d['y']-d['initial'];volume=d['volume']
        rows[label]=np.array([change[19*j]/np.sqrt(volume[j]*th['Cf']) for j in range(2)])
    fine=rows['fine-128'];norm=np.linalg.norm(fine)
    first=np.linalg.norm(rows['fine-32']-rows['fine-64'])/norm
    last=np.linalg.norm(rows['fine-64']-fine)/norm
    spatial=np.linalg.norm(rows['spatial-p2']-fine)/norm
    b,d,R,edges=model.geometry();s=np.load(model.patch.base.surface.OUT/'background.npz')
    rs=float(s['r'][-1]);r=s['r']/rs;gr=model.patch.base.h.gr;geo=gr.G*.1*R**2/gr.C**4
    p=np.interp(edges[1],r,s['p'])*rs**2/geo;e=np.interp(edges[1],r,s['e'])*rs**2/geo
    phi=np.interp(edges[1],r,s['phi']);native,_=model.patch.base.h.inputs();cx=float(native['CX'][0])
    expected=float(d['raw'][0,0]*(float(b['c'])**2*cx+d['raw'][0,2]))
    row=dict(classification='Counterexample candidate',seconds=time.monotonic()-start,
        material_temperature_time_relative=float(last),material_temperature_time_order=float(np.log2(first/last)),
        material_temperature_spatial_relative=float(spatial),GR_pressure_cgs=float(p),GR_energy_cgs=float(e),
        isotope_rest_CX=cx,native_energy_with_rest_CX_cgs=expected,
        pressure_match_relative=float(abs(p/d['raw'][0,1]-1)),energy_match_relative=float(abs(e/expected-1)),
        conformal_match_relative=float(abs(np.exp(-2*phi*phi)/b['local_A']-1)),
        interpretation='Use the saved isotope-rest CX multiplying c^2, distinct from the neutral-element gas energy normalization. The apparent0.7percent discrepancy from replacing CX by1 is a reference mismatch, not a measured GR/EOS energy defect.')
    assert max(row[k] for k in ['pressure_match_relative','energy_match_relative','conformal_match_relative'])<1e-12
    assert last<.02 and spatial<.02 and row['material_temperature_time_order']>1.5
    write(OUT/'state-check.json',row);print(json.dumps(row,indent=2))


def main():
    start=time.monotonic();plan=read(OUT/'plan.json')
    for name,h in plan['bindings'].items():assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==h,name
    th=np.load(OUT.parent/'def-photon-population-coupled/thermo.npz')
    bank=np.load(OUT.parent/'def-photon-hhe-coupled/bank.npz')
    info=read(OUT/'input.json');result=read(OUT/'result.json')
    n=info['photon_cells'];root=np.sqrt(bank['Ci'][:n]);weights=np.r_[np.sqrt(th['Cf']),th['g']]
    snapshots={};energy_rows={}
    for label in ['fine-128','one-way','spatial-p2']:
        if not (OUT/(label+'.npz')).exists():continue
        state=np.load(OUT/(label+'.npz'));delta=state['y']-state['initial'];V=state['volume']
        temperature=np.array([delta[19*j]/np.sqrt(V[j]*th['Cf']) for j in range(2)])
        material=np.array([np.sqrt(V[j])*(weights@delta[19*j:19*(j+1)]) for j in range(2)])
        photons=np.array([np.sqrt(V[j])*(root@delta[38+n*j:38+n*(j+1)]) for j in range(2)])
        ions=np.array([np.linalg.solve(th['U'],delta[19*j+1:19*(j+1)]/np.sqrt(V[j])) for j in range(2)])
        raw_energy=material+photons;canonical=V*state['ports'][:2]
        stored=read(OUT/(label+'.json'))['history'][-1]
        assert np.allclose(canonical,stored['cell_energy_erg'],rtol=1e-13,atol=0)
        assert np.linalg.norm(raw_energy-V*state['source_energy'])/np.linalg.norm(raw_energy)<1e-12
        energy_rows[label]=dict(material_energy_increment_erg=material.tolist(),
            photon_energy_increment_erg=photons.tolist(),nonadiabatic_energy_erg=canonical.tolist(),
            compression_subtraction_erg=(raw_energy-canonical).tolist(),
            material_temperature_increment_K=temperature.tolist(),
            ion_population_increment_norm_cm3=np.linalg.norm(ions,axis=1).tolist(),
            net_nonadiabatic_balance=float(abs(canonical.sum())/np.sum(abs(canonical))))
        snapshots[label]=(temperature,delta)
    fine=energy_rows['fine-128'];assert max(abs(np.array(fine['material_temperature_increment_K'])))>0
    row=dict(classification='Counterexample candidate',saved_input_bindings_passed=True,
        actual_material_energy_exchange_nonzero=True,endpoint=energy_rows,
        numerical_target_passed=result['passed'],seconds=time.monotonic()-start,
        scope='Endpoint component accounting and binding replay only; independent of the stage solver but not a continuum error certificate.')
    if 'one-way' in snapshots:
        dT=snapshots['fine-128'][0]-snapshots['one-way'][0]
        row['GR_return_temperature_increment_K']=dT.tolist()
        row['GR_return_thermal_state_norm_relative']=float(np.linalg.norm(snapshots['fine-128'][1]-snapshots['one-way'][1])/np.linalg.norm(snapshots['fine-128'][1]))
    if OUT.name.endswith('-resolved'):
        import def_photon_radial_gr_resolved as model
        signal.alarm(90);c=model.Coupled();m,p=c.m,c.p
        state=np.load(OUT/'fine-128.npz');rho=state['density'];delta=state['y']-state['initial'];V=state['volume']
        native=np.load(model.prior.ph.old.old.a.OUT/'eos-state.npz');st=np.load(model.prior.ph.old.old.RATES/'station.npz')
        convert=float(st['rho']*st['mass_scale'])*6.02214076e23
        def coords(i):return np.array([native['number_fractions'][i,e,q+1:].sum() for e,q,_ in th['edges']])*convert
        arho=(coords(6)*np.exp(1e-4)-coords(5)*np.exp(-1e-4))/(2e-4)
        U=th['U'];brho=arho-th['coordinates'];H=U.T@U/info['T_K'];latent=U.T@th['g']
        eos=native['eos'][0];pg=eos[1];pT=pg*eos[6]/info['T_K'];rhob=info['rho_g_cm3']
        material_T=snapshots['fine-128'][0];energy=[];pressure=[]
        for j in range(2):
            dn=np.linalg.solve(U,delta[19*j+1:19*(j+1)]/np.sqrt(V[j]))
            de_rad=root@delta[38+n*j:38+n*(j+1)]/np.sqrt(V[j])
            de_gas=float(th['Cf'])*material_T[j]+latent@dn+(rhob*eos[9]-latent@arho)*rho[j]
            energy.append(de_gas+de_rad-(pg+4*float(bank['arad'])*info['T_K']**4/3)*rho[j])
            dp_gas=pT*material_T[j]+pg*eos[5]*rho[j]-brho@H@(dn-th['a'][-1]*material_T[j]-arho*rho[j])
            pressure.append(dp_gas+de_rad/3-m.gammaP[j]*rho[j])
        manual=np.r_[energy,pressure];stored=state['ports'][:4]
        history=read(OUT/'fine-128.json')['history']
        scale=np.r_[np.max(abs(np.array([h['cell_energy_erg'] for h in history])/V),axis=0),
                    np.max(abs(np.array([h['pressure_increment_erg_cm3'] for h in history])),axis=0)]
        errors=abs(manual-stored)/scale
        assert max(errors)<1e-8,errors
        collision_state=state['y']-p.Eq@rho;x,E=collision_state[:38],collision_state[38:]
        material_rate=-(p.A@x+p.B.T@E);photon_rate=-(p.B@x+p.diag*E+p.off(E))
        material_power=np.array([np.sqrt(V[j])*weights@material_rate[19*j:19*(j+1)]/m.tc for j in range(2)])
        photon_power=np.array([np.sqrt(V[j])*root@photon_rate[n*j:n*(j+1)]/m.tc for j in range(2)])
        power_defect=float(np.linalg.norm(material_power+photon_power)/np.linalg.norm(material_power))
        assert power_defect<1e-8,power_defect
        row['independent_thermodynamic_ports_relative']=errors.tolist()
        row['material_collision_power_erg_s']=material_power.tolist()
        row['photon_collision_power_erg_s']=photon_power.tolist()
        row['collision_power_cancellation_relative']=power_defect
        row['seconds']=time.monotonic()-start;signal.alarm(0)
    write(OUT/'replay.json',row);print(json.dumps(row,indent=2))


if __name__=='__main__':
    if len(sys.argv)>1:OUT=Path(sys.argv[1])
    if '--state' in sys.argv:state_check()
    else:main()
