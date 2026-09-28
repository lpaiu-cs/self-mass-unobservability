"""Request 22: audit thermal and material closure of the frozen GR candidate.

Counterexample candidate: output-only MESA replay and explicitly frozen-coefficient
transport diagnostics. This module does not pretend to evolve a GR star.
"""
from pathlib import Path
import json, os, shutil, subprocess, sys, time
import numpy as np
from thermal_wd import mesa
from thermal_restart import sha
from thermal_robustness import gzcopy
import gr_mass as gr
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'outputs/thermal-closure22'
CACHE = Path('/home/lpaiu/work/thermal-closure22')
MESA = Path('/home/lpaiu/work/thermal-restart19/mesa-r7624')
SOURCE_RUN = Path('/home/lpaiu/work/thermal-robustness20/mesh')
OLD = ROOT/'outputs/gr-mass21'
LSUN = 3.8418e33  # exact value checked against the distributed const_def.f90
EXTRA = ['luminosity', 'eps_grav', 'eps_nuc', 'eps_nuc_neu_total', 'non_nuc_neu',
         'extra_heat', 'opacity', 'gradT', 'actual_gradT', 'gradr', 'grada',
         'mixing_type', 'lum_conv_div_L', 'cp', 'cv', 'entropy',
         'd_lnepsnuc_dlnT', 'd_lnepsnuc_dlnd']


def save(name, obj):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/name).write_text(json.dumps(obj, ensure_ascii=False, indent=2)+'\n')


def prepare():
    assert not CACHE.exists() and not OUT.exists()
    CACHE.mkdir(); OUT.mkdir(parents=True)
    save('plan.json', dict(classification='Counterexample candidate', checkpoint='cdab2c5',
        source_model=19041, restart_model=19000, stop_model=19042,
        replay='Same mesh trajectory and binary, output-only additional columns. Require exact equality of all shared profile values at model 19041 before using new columns.',
        tests=['Original discrete luminosity and nuclear/neutrino/gravothermal budget',
               'GR heat flux and redshifted luminosity balance with explicit positive-conductivity assumptions',
               'Proper baryon and isotope inventory after pressure-attached remapping',
               'Frozen original opacity and heat-source controls, never relabelled as re-evaluated microphysics'],
        invariant='Frozen Request21 mass match and primary optical failure remain unchanged.',
        boundaries='No new stellar evolution code, network or opacity re-evaluation at the GR state, observational inference, or transport error certificate. A failed steady-state closure does not rule out a time-dependent star.',
        previous_manifest_sha256=sha(OLD/'manifest.json')))
    run=CACHE/'replay'; run.mkdir(); (run/'photos1').mkdir()
    for p in SOURCE_RUN.iterdir():
        if p.is_file() and (p.name.startswith('inlist') or p.suffix in ('.list','.net') or p.name=='binary'):
            shutil.copy2(p,run/p.name)
    shutil.copy2(SOURCE_RUN/'photos1/19000',run/'photos1/19000')
    (run/'.restart').write_text('19000\n')
    p=run/'inlist1'; content=p.read_text()
    assert content.count('max_model_number = 40000')==1
    p.write_text(content.replace('max_model_number = 40000','max_model_number = 19042'))
    p=run/'profile_columns.list'
    p.write_text(p.read_text()+'\n! Request22: output only, unchanged physics.\n'+'\n'.join(EXTRA)+'\n')
    inputs={p.relative_to(run).as_posix():sha(p) for p in sorted(run.rglob('*')) if p.is_file()}
    dest=OUT/'inputs'; dest.mkdir()
    for rel in inputs:
        if rel in ('binary','photos1/19000'): continue
        shutil.copy2(run/rel,dest/('restart-input.txt' if rel=='.restart' else rel))
    save('inputs.json',dict(classification='Proven',sha256=inputs,
        run=str(run), binary_from=str(SOURCE_RUN/'binary'), restart_from=str(SOURCE_RUN/'photos1/19000')))
    print('Prepared output-only replay 19000 -> 19042',flush=True)


def run():
    folder=CACHE/'replay'; log=folder/'execution.log'; assert not log.exists()
    for rel,digest in json.loads((OUT/'inputs.json').read_text())['sha256'].items():
        assert sha(folder/rel)==digest,rel
    runtime=MESA.parent; env=os.environ.copy()
    env.update(MESA_DIR=str(MESA),LD_LIBRARY_PATH=str(runtime/'mesasdk/lib')+':'+str(runtime/'mesasdk/lib64'),OMP_NUM_THREADS='4')
    start=time.monotonic()
    with log.open('w') as f:
        result=subprocess.run(['./binary'],cwd=folder,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=1800)
    shutil.copy2(log,OUT/'replay.log')
    save('execution.json',dict(classification='Proven',returncode=result.returncode,
        elapsed_s=time.monotonic()-start,environment={k:env[k] for k in ['MESA_DIR','LD_LIBRARY_PATH','OMP_NUM_THREADS']}))
    assert result.returncode==0
    print('Replay process exited',result.returncode,flush=True)


def collect():
    logs=CACHE/'replay/LOGS1'
    index=np.loadtxt(logs/'profiles.index',skiprows=1,dtype=int)
    row=index[index[:,0]==19041]; assert len(row)==1
    source=logs/f'profile{row[0,2]}.data'
    h,d=mesa(source); oh,od=mesa(gr.PROFILE)
    assert int(h['model_number'])==19041
    checks={k:bool(np.array_equal(d[k],v)) for k,v in od.items()}
    assert all(checks.values()),[k for k,v in checks.items() if not v]
    assert float(h['star_age'])==float(oh['star_age'])
    assert set(EXTRA)<=set(d)
    gzcopy(source,OUT/'selected.data.gz');gzcopy(logs/'history.data',OUT/'history.data.gz')
    shutil.copy2(logs/'profiles.index',OUT/'profiles.index')
    save('replay-control.json',dict(classification='Proven',model=19041,
        zones=len(d['zone']),shared_column_equality=checks,age_equal=True,
        original_profile_sha256=sha(gr.PROFILE),full_profile_raw_sha256=sha(source),
        header_note='Headers include output locations and restarted elapsed diagnostics; all shared structure columns and age, not every header string, are required identical.'))
    print('PASS exact shared profile columns:',len(checks),'zones:',len(d['zone']))


def source_audit():
    names=['const/public/const_def.f90','star/defaults/profile_columns.list',
           'star/private/profile_getval.f90','star/private/hydro_eqns.f90',
           'star/private/eps_grav.f90','kap/private/kap_eval.f90','kap/private/condint.f90']
    bindings={}
    for name in names:
        src=MESA/name;dest=OUT/'sources'/name;dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(src,dest);bindings[name]=sha(dest)
    eq=(MESA/names[3]).read_text()
    assert 'non_nuc_neu = 0.5d0*(s% non_nuc_neu_start(k) + s% non_nuc_neu(k))' in eq
    assert 'eps_burn = eps_nuc - non_nuc_neu + s% extra_heat(k) + s% irradiation_heat(k)' in eq
    assert 'dLdm = s% eps_grav(k) + eps_burn' in eq
    assert 'lsol = 3.8418d33' in (MESA/names[0]).read_text()
    assert 'kap = 1d0 / (1d0/kap_rad + 1d0/kap_ec)' in (MESA/names[5]).read_text()
    save('source-audit.json',dict(classification='Imported from prior work',mesa_release=7624,
        sha256=bindings,Lsun_erg_s=LSUN,boltz_sigma_cgs=5.670400e-5,
        energy='eps_nuc already excludes reaction neutrinos. The discrete equation subtracts the mean of start/end non-nuclear neutrino loss, then adds eps_grav, extra and irradiation heating; viscosity is optional.',
        caution='The output contains end-state non_nuc_neu, not its remeshed start-state value. A budget using only end-state profile values is an instantaneous diagnostic, not the exact discrete solver residual.',
        opacity='Radiative and electron-conductive opacities combine harmonically. Convection is a separate channel.',
        reference=dict(title='Lander and Andersson (2018), Heat conduction in rotating relativistic stars',
            url='https://doi.org/10.1093/mnras/sty1725',
            use='Primary reference for the quasi-steady Fourier limit and its causal/dynamical limitations. Request22 spherical formulas are also derived and symbolically checked here.')))


def budget():
    h,d=mesa(OUT/'selected.data.gz');dm=d['dm']; L=d['luminosity']*LSUN
    integrals={k:float(np.dot(dm,d[k])/LSUN) for k in
               ['eps_nuc','eps_nuc_neu_total','non_nuc_neu','eps_grav','extra_heat']}
    # No second subtraction of reaction neutrinos from eps_nuc.
    q=d['eps_nuc']-d['non_nuc_neu']+d['eps_grav']+d['extra_heat']
    generated=dm*q/LSUN; emitted=(L-np.r_[L[1:],0.])/LSUN
    source_net=integrals['eps_nuc']-integrals['non_nuc_neu']+integrals['extra_heat']
    peak=int(np.argmax(d['eps_nuc']))
    _,hist=mesa(OUT/'history.data.gz'); ix=np.flatnonzero(hist['model_number']==19041)
    assert len(ix)==1; i=int(ix[0]); assert 0<i<len(hist['star_age'])-1
    dt=float(hist['star_age'][i+1]-hist['star_age'][i-1])
    slopes={k:float((hist[k][i+1]-hist[k][i-1])/dt) for k in ['log_Teff','log_R','log_L']}
    conv=np.abs(d['lum_conv_div_L'])>.01
    body=dict(classification='Counterexample candidate',model=19041,
        integrated_Lsun=integrals,surface_Lsun=float(L[0]/LSUN),
        photospheric_header_Lsun=float(h['photosphere_L']),
        net_nuclear_minus_thermal_neutrino_Lsun=source_net,
        steady_state_missing_sink_Lsun=source_net-float(L[0]/LSUN),
        steady_state_source_to_surface_ratio=source_net/float(L[0]/LSUN),
        gravothermal_to_surface_ratio=integrals['eps_grav']/float(L[0]/LSUN),
        instantaneous_budget_residual_Lsun=float(sum(generated-emitted)),
        instantaneous_budget_absolute_cell_residual_Lsun=float(sum(abs(generated-emitted))),
        discrete_solver_residual_measured=False,
        hydrogen_burning_peak=dict(zone=int(d['zone'][peak]),T_K=float(10**d['logT'][peak]),
            rho_cgs=float(d['rho'][peak]),enclosed_source_Msun=float(d['mass'][peak]),
            eps_nuc_erg_g_s=float(d['eps_nuc'][peak])),
        logarithmic_history_slopes_per_year=slopes,
        adjacent_history_span_year=dt,
        radius_e_folding_year=float(1/(np.log(10)*abs(slopes['log_R']))),
        orbit_days=1.6294,
        radius_fractional_linear_change_per_orbit=float(np.log(10)*slopes['log_R']*1.6294/365.25),
        time_scale_note='Adjacent-model secular slope, not a relaxation pole, formation probability or periodic response.',
        convective_flux_gt_1percent_zone_count=int(sum(conv)),
        convective_flux_gt_1percent_mass_fraction=float(np.dot(dm,conv)/sum(dm)),
        source_interpretation='Nonzero gravothermal storage/release is essential for this evolving MESA snapshot. These values are not re-evaluated rates of the GR candidate.')
    save('source-budget.json',body)
    np.savez_compressed(OUT/'source-budget.npz',zone=d['zone'],dm_g=dm,
        emitted_Lsun=emitted,instantaneous_generated_Lsun=generated)
    print(json.dumps(body,ensure_ascii=False,indent=2))


def load_star(label):
    calibrated=label=='calibrated'
    record=json.loads((OLD/('calibrated-final.json' if calibrated else 'gr-refined-4.json')).read_text())
    grid=np.load(OLD/('calibrated-eos-final.npz' if calibrated else 'eos-grid-4.npz'))
    path=gr.ThermalPath(temperature_scale=record.get('temperature_scale',1.))
    star=gr.Star(np.log(record['central_pressure_cgs']),gr.EOS(),path,
        rtol=1e-12,r0=10.,table=gr.table_interpolant(grid['logP'],grid['values']),max_step=.012)
    assert abs(star.mass/gr.TARGET-1)<1e-8
    return star


def remap(label, subdivision):
    h,d=mesa(OUT/'selected.data.gz');star=load_star(label);lp=star.path.lp
    # Quadrature uses exact cumulative proper masses at pressure-cell edges.
    # ponytail: midpoint quadrature with refinement is diagnostic, not interval certification.
    edge=np.unique(np.r_[np.concatenate([lp[:-1]+j/subdivision*np.diff(lp)
        for j in range(subdivision)]),lp[-1],star.sol.t[0]])
    edge=edge[edge<=star.sol.t[0]]
    state=star.sol.sol(edge); shell=-(np.diff(state[2]))*gr.C**2/gr.G*1000
    mid=(edge[1:]+edge[:-1])/2; states=star.sol.sol(mid)
    assert np.all(shell>0)
    r,m,mb,nu=states;nu=nu+star.nu_shift;f=1-2*m/r
    temp=np.exp(star.path.lt(np.clip(mid,lp[0],lp[-1])))*star.path.temperature_scale
    nabla=star.path.lt.derivative()(np.clip(mid,lp[0],lp[-1]));nabla[mid>lp[-1]]=0.
    cx=star.path.cx(np.clip(mid,lp[0],lp[-1]))
    eos=np.array([star.table(p) for p in mid]);rho_b=eos[:,0]/cx
    egeom=gr.G*eos[:,0]*1000/gr.C**2*(1+eos[:,2]*1e-4/gr.C**2)
    pgeom=gr.G*np.exp(mid)*.1/gr.C**4
    dlndr=-(egeom+pgeom)*(m+4*np.pi*r**3*pgeom)/(r*r*f*pgeom) # per metre
    tolman=pgeom/(egeom+pgeom)
    opacity=np.exp(np.interp(mid,lp,np.log(d['opacity'])))
    conductivity=16*5.670400e-5*temp**3/(3*opacity*rho_b)
    luminosity=-4*np.pi*(100*r)**2*conductivity*np.sqrt(f)*temp*(dlndr/100)*(nabla-tolman)
    source_L=np.interp(mid,lp,d['luminosity'])*LSUN
    source_conv=np.interp(mid,lp,d['lum_conv_div_L'])
    source_diff=source_L*(1-source_conv)
    q_source={k:np.interp(mid,lp,d[k]) for k in ['eps_nuc','non_nuc_neu','eps_grav','extra_heat']}
    power={k:float(np.dot(np.exp(2*nu)*shell,v)/LSUN) for k,v in q_source.items()}
    isos=json.loads((ROOT/'outputs/thermal-restart19/network-diagnosis.json').read_text())['profile_isotopes']
    keep=[k for k in isos if not k.startswith('f')]
    norm=sum(d[k] for k in keep)
    fractions=np.array([np.interp(mid,lp,d[k]/norm) for k in keep])
    inventories=fractions@shell; central_mass=state[2,-1]*gr.C**2/gr.G*1000
    # Include the tiny unresolved centre and matched-composition mathematical atmosphere.
    extra_mass=central_mass+(star.baryon-star.photosphere[2])*gr.C**2/gr.G*1000
    inventories+=central_mass*np.array([d[k][-1]/norm[-1] for k in keep])
    inventories+=(extra_mass-central_mass)*np.array([d[k][0]/norm[0] for k in keep])
    source_inventory=np.array([np.dot(d['dm'],d[k]/norm) for k in keep])
    source_baryon=float(sum(d['dm']));new_baryon=star.baryon*gr.C**2/gr.G*1000
    assert abs(sum(inventories)/new_baryon-1)<1e-10
    # Exclude convection, near-zero luminosity and the optically thin boundary in this diagnostic.
    domain=(abs(source_conv)<.01)&(abs(source_diff)>1e-3*LSUN)&(temp>1e5)&(mid<lp[-1])
    ratio=luminosity[domain]/source_diff[domain]
    steady=power['eps_nuc']-power['non_nuc_neu']+power['extra_heat']
    allq=steady+power['eps_grav'];target_red=float(h['photosphere_L'])*np.exp(2*(star.photosphere[3]+star.nu_shift))
    result=dict(classification='Counterexample candidate',label=label,subdivision=subdivision,
        quadrature_shell_count=len(shell),source_baryon_mass_g=source_baryon,proper_baryon_mass_g=float(new_baryon),
        baryon_inventory_relative_change=float(new_baryon/source_baryon-1),
        isotope_inventory={k:dict(original_renormalized_g=float(a),GR_g=float(b),relative_change=float(b/a-1),
            normalized_abundance_relative_change=float((b/new_baryon)/(a/source_baryon)-1))
            for k,a,b in zip(keep,source_inventory,inventories)},
        omitted_fluorine_original_inventory_g=float(sum(np.dot(d['dm'],d[k]) for k in isos if k.startswith('f'))),
        frozen_source_power_at_infinity_Lsun=power,
        frozen_steady_source_at_infinity_Lsun=steady,
        frozen_full_source_at_infinity_Lsun=allq,
        retained_photospheric_L_at_infinity_Lsun=target_red,
        frozen_full_budget_residual_Lsun=allq-target_red,
        frozen_steady_source_to_surface_ratio=steady/target_red,
        frozen_opacity_diffusion_to_original_flux_quantiles=dict(zip(['p01','p10','p50','p90','p99'],map(float,np.quantile(ratio,[.01,.1,.5,.9,.99])))),
        frozen_opacity_diffusion_sign_mismatch_samples=int(sum(ratio<0)),diagnostic_samples=int(sum(domain)),
        flux_domain='abs(source convective fraction)<0.01, abs(source diffusive L)>1e-3 Lsun, T>1e5 K, below central extension. Quantiles count pressure subcells, not mass weights.',
        temperature_scale=star.path.temperature_scale,
        material_entropy_evolution_solved=False,
        boundary='Frozen opacity and reaction/gravothermal rates from the MESA pressure path. They are not microphysics at the GR density/temperature. Failure diagnoses this transport transplant, not nonexistence of an evolving GR star.')
    save(f'{label}-{subdivision}.json',result)
    np.savez_compressed(OUT/f'{label}-{subdivision}.npz',logP=mid,r_m=r,dmB_g=shell,nu=nu,
        T_K=temp,rhoB_cgs=rho_b,nabla=nabla,tolman_nabla=tolman,
        frozen_kappa_cgs=opacity,local_diffusive_L_erg_s=luminosity,original_diffusive_L_erg_s=source_diff,
        diagnostic_domain=domain)
    print(label,subdivision,'baryon change',result['baryon_inventory_relative_change'],
          'H change',result['isotope_inventory']['h1']['relative_change'],
          'frozen budget',allq,target_red,flush=True)


def diagnose():
    for label in ['primary','calibrated']:
        for subdivision in [1,2,4]: remap(label,subdivision)


def independent_inventory(coordinate='sqrt_pressure'):
    h,d=mesa(OUT/'selected.data.gz');lp=d['logP']*np.log(10.)
    keep=[k for k in json.loads((ROOT/'outputs/thermal-restart19/network-diagnosis.json').read_text())['profile_isotopes'] if not k.startswith('f')]
    norm=sum(d[k] for k in keep);results={}
    for label in ['primary','calibrated']:
        star=load_star(label);values=[];powers=[]
        edges=np.r_[lp,star.sol.t[0]];assert np.all(np.diff(edges)>0)
        for order in [4,8]:
            nodes,weights=np.polynomial.legendre.leggauss(order)
            if coordinate=='pressure':
                x=(edges[:-1,None]+edges[1:,None])/2+np.diff(edges)[:,None]*nodes/2
                w=np.diff(edges)[:,None]*weights/2
            else:
                # At the regular centre r ~ sqrt(ln Pc-ln P). This coordinate
                # removes the square-root endpoint from the differential volume.
                a=np.sqrt(star.logpc-edges[1:]);b=np.sqrt(star.logpc-edges[:-1])
                u=(a[:,None]+b[:,None])/2+(b-a)[:,None]*nodes/2
                x=star.logpc-u*u;w=(b-a)[:,None]*weights*u
            x=x.ravel();w=w.ravel();r,m,mb,nu=star.sol.sol(x)
            eos=np.array([star.table(p) for p in x])
            cx=star.path.cx(np.clip(x,lp[0],lp[-1]));rhoB=eos[:,0]/cx*1000
            p=gr.G*np.exp(x)*.1/gr.C**4
            e=gr.G*eos[:,0]*1000/gr.C**2*(1+eos[:,2]*1e-4/gr.C**2)
            f=1-2*m/r
            dr=-p*r*r*f/((e+p)*(m+4*np.pi*r**3*p))
            # Independent differential proper-volume integral; no cumulative-MB differences.
            weight=-4*np.pi*r*r*rhoB/np.sqrt(f)*dr*w*1000
            tiny_centre=star.sol.y[2,0]*gr.C**2/gr.G*1000
            atmo=(star.baryon-star.photosphere[2])*gr.C**2/gr.G*1000
            masses={k:float(np.dot(weight,np.interp(x,lp,d[k]/norm))+
                tiny_centre*d[k][-1]/norm[-1]+atmo*d[k][0]/norm[0]) for k in ['h1','he4']}
            masses['baryon']=float(sum(weight)+tiny_centre+atmo);values.append(masses)
            q=np.interp(x,lp,d['eps_nuc']-d['non_nuc_neu']+d['eps_grav']+d['extra_heat'])
            powers.append(float(np.dot(weight*np.exp(2*(nu+star.nu_shift)),q)/LSUN))
        comparison=json.loads((OUT/f'{label}-4.json').read_text())
        errors={k:float(values[1][k]/comparison['isotope_inventory'][k]['GR_g']-1) for k in ['h1','he4']}
        errors['baryon']=float(values[1]['baryon']/comparison['proper_baryon_mass_g']-1)
        refinement={k:float(values[1][k]/values[0][k]-1) for k in values[0]}
        results[label]=dict(differential_Gauss_4_and_8_g=values,
            relative_difference_from_cumulative_midpoint=errors,Gauss_refinement=refinement,
            frozen_full_source_Gauss_4_and_8_Lsun=powers,
            frozen_full_source_difference_from_midpoint_Lsun=powers[1]-comparison['frozen_full_source_at_infinity_Lsun'])
    passed=all(max(abs(v) for v in row['relative_difference_from_cumulative_midpoint'].values())<1e-6 and
        max(abs(v) for v in row['Gauss_refinement'].values())<1e-7 for row in results.values())
    save('inventory-pressure-pilot.json' if coordinate=='pressure' else 'inventory-audit.json',dict(classification='Proven',
        coordinate=coordinate,passed=passed,
        result=('Independent inventory agreement passed.' if passed else
            'Direct lnP Gauss pilot failed near the square-root central endpoint; preserved before regularizing the quadrature coordinate.')+
            ' No interval or physical EOS error certificate is implied.',results=results))
    print(json.dumps(results,indent=2))
    if coordinate!='pressure': assert passed


def inventory_pressure_pilot():
    independent_inventory('pressure')


def source_flux_control():
    _,d=mesa(OUT/'selected.data.gz');lp=d['logP']*np.log(10.)
    mid=(lp[:-1]+lp[1:])/2
    # Matching the centre-sampled temperatures to volume-midpoint radii avoids
    # treating zone-face radius and zone-centre temperature as the same location.
    outer=d['radius']*gr.RSUN*100;inner=np.r_[outer[1:],0.]
    radius=np.cbrt((outer**3+inner**3)/2)
    rpath=PchipInterpolator(lp,np.log(radius));r=np.exp(rpath(mid))
    tpath=PchipInterpolator(lp,d['logT']*np.log(10.));T=np.exp(tpath(mid))
    rho=np.exp(np.interp(mid,lp,d['logRho']*np.log(10.)))
    kap=np.exp(np.interp(mid,lp,np.log(d['opacity'])))
    K=16*5.670400e-5*T**3/(3*kap*rho)
    gradient=T*tpath.derivative()(mid)/(r*rpath.derivative()(mid))
    calc=-4*np.pi*r*r*K*gradient
    L=np.interp(mid,lp,d['luminosity'])*LSUN;conv=np.interp(mid,lp,d['lum_conv_div_L'])
    actual=L*(1-conv);domain=(abs(conv)<.01)&(abs(actual)>1e-3*LSUN)&(T>1e5)
    ratio=calc[domain]/actual[domain]
    save('source-flux-control.json',dict(classification='Counterexample candidate',
        ratio_quantiles=dict(zip(['p01','p10','p50','p90','p99'],map(float,np.quantile(ratio,[.01,.1,.5,.9,.99])))),
        domain_samples=int(sum(domain)),sign_mismatch_samples=int(sum(ratio<0)),
        limitation='PCHIP centre gradients and volume-midpoint radii reconstruct a continuous diagnostic from finite MESA zones; this is not the exact face MLT stencil or a transport convergence certificate. The same finite-profile limitation applies to GR frozen-opacity flux comparisons.'))
    print('Original-profile flux control quantiles:',np.quantile(ratio,[.01,.1,.5,.9,.99]))


def entropy_control():
    a=np.load(OLD/'eos-grid-4.npz');b=np.load(OLD/'calibrated-eos-final.npz')
    assert np.array_equal(a['logP'],b['logP'])
    p=a['logP'];path=gr.ThermalPath();cx=path.cx(np.clip(p,path.lp[0],path.lp[-1]))
    delta=(b['values'][:,3]-a['values'][:,3])*cx
    assert np.all(delta>0)
    save('entropy-control.json',dict(classification='Counterexample candidate',
        same_pressure_and_composition_grid=True,states=len(p),
        positive_entropy_change_at_all_sampled_states=True,
        delta_specific_entropy_per_baryon_gram_cgs_min=float(min(delta)),
        delta_specific_entropy_per_baryon_gram_cgs_max=float(max(delta)),
        scope='FreeEOS entropy at fixed P,X before/after the uniform temperature calibration; sampled-state check, not a material-shell mapping or entropy-evolution solution.'))
    print('Positive fixed-P,X EOS entropy shift at',len(p),'states:',min(delta),max(delta))


def symbolic():
    import sympy as s
    r=s.symbols('r',positive=True);P,T,nu=s.symbols('P T nu',cls=s.Function)
    f,K,area,e=s.symbols('f K area e',positive=True)
    dlogP=s.diff(P(r),r)/P(r);nabla=(s.diff(T(r),r)/T(r))/dlogP
    flux=-K*s.sqrt(f)*(s.diff(T(r),r)+T(r)*s.diff(nu(r),r))
    formula=-K*s.sqrt(f)*T(r)*dlogP*(nabla-P(r)/(e+P(r)))
    assert s.simplify((flux-formula).subs(s.diff(nu(r),r),-s.diff(P(r),r)/(e+P(r))))==0
    Tinf=s.symbols('Tinf',positive=True)
    assert s.simplify(flux.subs(T(r),Tinf*s.exp(-nu(r))).doit())==0
    q,mbp,rhob=s.symbols('q mbp rhob',positive=True);L=s.Function('L')
    mbprime=area*rhob/s.sqrt(f)
    luminosity_prime=s.exp(2*nu(r))*q*mbprime
    assert s.simplify(luminosity_prime/mbprime/s.exp(2*nu(r))-q)==0
    source,sdot=s.symbols('source sdot')
    residual=source-T(r)*sdot-q
    assert s.simplify(residual.subs(sdot,(source-q)/T(r)))==0
    # Uniform-source flat-space sphere is an independent sign/normalization control.
    rho0,k0,q0,t0=s.symbols('rho0 k0 q0 t0',positive=True)
    trial=t0-rho0*q0*r*r/(6*k0)
    trial_L=s.simplify(-4*s.pi*r*r*k0*s.diff(trial,r))
    assert s.simplify(trial_L-4*s.pi*rho0*q0*r**3/3)==0
    # Exact static diagonal metric: R_tr=G_tr=0. A nonzero net radial
    # luminosity therefore requires the quasi-static approximation or evolution.
    t,theta,phi=s.symbols('t theta phi');coords=[t,r,theta,phi]
    ff=s.Function('f')(r)
    metric=s.diag(-s.exp(2*nu(r)),1/ff,r*r,r*r*s.sin(theta)**2)
    inv=metric.inv()
    def christ(a,b,c):
        return s.simplify(sum(inv[a,j]*(s.diff(metric[j,b],coords[c])+
            s.diff(metric[j,c],coords[b])-s.diff(metric[b,c],coords[j]))/2 for j in range(4)))
    gam=[[[christ(a,b,c) for c in range(4)] for b in range(4)] for a in range(4)]
    Rtr=sum(s.diff(gam[a][0][1],coords[a])-s.diff(gam[a][0][a],r)+
        sum(gam[a][0][1]*gam[b][a][b]-gam[b][0][a]*gam[a][1][b] for b in range(4)) for a in range(4))
    assert s.simplify(Rtr)==0 and metric[0,1]==0
    F=s.symbols('F');Ttr=-s.exp(nu(r))*F/s.sqrt(ff)
    assert s.solve(Ttr,F)==[0]
    save('symbolic-checks.json',dict(classification='Proven',
        pressure_gradient_flux=True,Tolman_zero_flux=True,proper_volume_redshifted_energy_balance=True,
        undetermined_entropy_rate_identity=True,uniform_source_sphere_control=True,
        static_metric_Gtr_zero=True,static_comoving_net_heat_flux_zero=True,
        conditions='Spherical quasi-static metric, local positive scalar conductivity and negligible heat-flux backreaction. The steady Fourier law is not a proof of a causal relaxation time or dynamical stability.'))
    print('PASS: seven symbolic transport, conservation and exact-static-boundary controls')


def finalize():
    budget=json.loads((OUT/'source-budget.json').read_text())
    audit=json.loads((OUT/'inventory-audit.json').read_text())
    pilot=json.loads((OUT/'inventory-pressure-pilot.json').read_text())
    primary=json.loads((OUT/'primary-4.json').read_text())
    calibrated=json.loads((OUT/'calibrated-4.json').read_text())
    assert audit['passed'] and not pilot['passed']
    assert budget['steady_state_source_to_surface_ratio']>100
    assert abs(primary['isotope_inventory']['h1']['relative_change'])>.1
    assert abs(calibrated['isotope_inventory']['h1']['relative_change'])>.05
    assert all(abs(v['frozen_full_source_difference_from_midpoint_Lsun'])<1e-5 for v in audit['results'].values())
    save('gates.json',dict(classification='Counterexample candidate',
        original_shared_structure_exactly_reproduced=True,independent_inventory_audit_passed=True,
        prior_GR_mass_match_preserved=True,original_snapshot_thermal_steady_state=False,
        pressure_attached_same_material_interpretation_supported=False,
        unchanged_opacity_heat_source_transplant_closes=False,
        exact_static_metric_supports_nonzero_net_comoving_luminosity=False,
        seven_symbolic_checks_passed=True,original_discrete_MESA_energy_residual_certified=False,
        actual_GR_thermal_composition_evolution_solved=False,full_22_isotope_microphysics=False,
        GR_fluid_scalar_dynamics_or_stability_certified=False,
        complete_nonlinear_observational_inference=False,final_PDF_or_ZIP_updated=False,
        classification_detail='Theorem progress: exact-static heat-flux boundary and redshifted conservation identities. Loophole progress: quantified failure of same-material and frozen thermal-coefficient interpretations; mathematical GR mass match remains.'))


def maintain():
    old=json.loads((OLD/'manifest.json').read_text())['sha256']
    additions={
        'model-definition': '분류: Counterexample candidate. Request21의 조정 후보는 원래 재규격화 프로필보다 총 바리온 질량이 0.069465%, 수소 재고가 5.42011% 작다. 같은 X(P)는 같은 물질을 보존하는 변환이 아니다. 별도 모형족의 조건부 GR 질량 해는 유지하지만, 원래 항성의 보존적 GR 변환이라는 해석은 통과하지 않는다.',
        'observable-targets': '분류: Counterexample candidate. 원래 프로필을 정확히 재현해 추가 열수송 자료를 복원했다. 조정 GR 구조에 원래 불투명도·열원을 동결 이식하면 확산 광도비 중앙값은 1.08982이고, 중력열을 포함한 적색편이 순 광도는 0.274013 L_sun으로 유지한 0.552221과 맞지 않는다. 이는 동결 계수 실험이며 새 상태의 실제 반응·수송 계산이나 독립 광학 예측이 아니다.',
        'adiabatic-limit': '분류: Proven. 양의 국소 전도계수와 준정적 구대칭 계량에서 열유속은 -K sqrt(f) (dT/dr + T dnu/dr)이며 영유속 조건은 T exp(nu)=상수다. 정확한 정적 대각 계량과 정지 물질, 다른 상쇄 흐름 부재에서는 G_tr=0이 순 열유속 0을 요구한다. 유한 광도를 쓰는 항성에는 준정적 근사 또는 시간 진화를 명시해야 한다. 이 경계는 궤도 완화시간이나 보편 동적 no-go를 증명하지 않는다.',
        'nonadiabatic-regime': '분류: Counterexample candidate. 원래 스냅숏의 순 핵반응 열원은 표면 광도의 약 284배이고 큰 음의 중력열 항을 가진다. 끝 상태 출력만의 에너지 잔차도 보존했다. 주변 모델에서 추정한 반지름 변화 시간 약 11617년은 진화 기울기이며 동적 응답의 pole이 아니다.\n\n분류: Conjectural. 다음 경계는 바리온 좌표의 물질·조성 보존, 새 상태의 실제 반응·수송 계수, GR 열·조성 진화와 안정성·유체·metric·scalar 전달함수를 순서대로 연결하는 것이다.',
        'failure-ledger-dynamic-chi': '분류: Counterexample candidate. 실패 단계는 질량 적분 이후에 원래 X(P)·광도·열 저장 항을 그대로 옮겨 동일한 실제 항성으로 해석하는 부분이다. 총 질량을 맞춰도 수소 재고와 열수송 수지가 남는다. 미정 엔트로피 변화율이나 온도 배율을 조정하는 것만으로 실제 진화 해결로 세지 않는다.\n\n분류: Proven. 직접 압력 좌표의 독립 Gauss 재고 검산은 중심 끝점 때문에 실패했다. 실패 수치를 보존하고 sqrt(ln Pc-ln P) 좌표로 정칙화한 독립 적분이 수소 약 2.03e-8, 총 바리온 약 1.22e-12의 상대 차이로 통과했다.\n\n분류: Conjectural. 최소 추가 조건은 물질 좌표와 각 반응·손실·수송 계수를 닫은 준정적 GR 진화다. Request21 질량 매칭과 고정 온도 광학 실패는 그대로 보존하며, 전체 항성·미분 인증·관측 추론 완료를 선언하지 않는다.'}
    dest=OUT/'request21-notes';dest.mkdir(exist_ok=False);bindings={}
    for rel in ['docs/'+k+'.md' for k in additions]+['paper/revision-manifest.json']:
        assert sha(ROOT/rel)==old[rel],rel
        snap=dest/Path(rel).name;shutil.copy2(ROOT/rel,snap)
        bindings[rel]=dict(snapshot=snap.relative_to(ROOT).as_posix(),sha256=old[rel],historical_manifest=rel.startswith('paper/'))
    save('historical-note-bindings.json',bindings)
    revision=json.loads((ROOT/'paper/revision-manifest.json').read_text())
    for stem,body in additions.items():
        rel='docs/'+stem+'.md'
        with (ROOT/rel).open('ab') as f:
            f.write(('\n\n## Request 22 GR 후보의 열수송과 물질 보존\n\n'+body+
                '\n\n세부 근거: [한글 보고서](../notes/REQUEST22_THERMAL_CLOSURE_KO.md).\n').encode())
        revision['sha256'][rel]=sha(ROOT/rel)
    revision['request22_supporting_note_update']=dict(evidence_manifest='outputs/thermal-closure22/manifest.json',
        historical_notes='outputs/thermal-closure22/historical-note-bindings.json',
        status='GR 질량 해 보존; 동일 물질·동결 열수송 이식의 실패 수치와 정확한 정적 열유속 no-go 확인',
        artifact_status='원고 PDF와 ZIP은 Request 12 동결본이다.')
    (ROOT/'paper/revision-manifest.json').write_text(json.dumps(revision,ensure_ascii=False,indent=2)+'\n')


def seal():
    import scipy,sympy
    save('provenance.json',dict(classification='Proven',before_task_checkpoint='cdab2c5',
        previous_manifest_sha256=sha(OLD/'manifest.json'),interpreter=sys.executable,
        versions={m.__name__:m.__version__ for m in [np,scipy,sympy]},
        runtime_provenance='outputs/gr-mass21/provenance.json',
        replay_runtime_and_photo_bindings='outputs/thermal-closure22/inputs.json',
        producer='verification/thermal_closure.py'))
    files=[p for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json']
    files += [ROOT/'verification/thermal_closure.py',ROOT/'notes/REQUEST22_THERMAL_CLOSURE_KO.md',ROOT/'.gitattributes']
    files += [ROOT/k for k in json.loads((OUT/'historical-note-bindings.json').read_text())]
    save('manifest.json',dict(classification='Proven',scope='Thermal and material closure audit of the frozen matched GR candidate',
        sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(files)}))


def verify():
    stages=['validated-variational','remaining-levers15','nbody-readout16','nonzero-drive17',
            'thermal-wd18','thermal-restart19','thermal-robustness20','gr-mass21','thermal-closure22']
    histories={f'outputs/{a}/manifest.json':f'outputs/{b}/historical-note-bindings.json' for a,b in zip(stages,stages[1:])}
    count=0
    for label in ['outputs/research-remediation/manifest.json',*histories,'paper/revision-manifest.json','outputs/thermal-closure22/manifest.json']:
        old=json.loads((ROOT/histories[label]).read_text()) if label in histories else {}
        for name,expected in json.loads((ROOT/label).read_text())['sha256'].items():
            path=ROOT/name
            if name in old:
                bind=old[name];path=ROOT/bind['snapshot'];assert expected==bind['sha256']
                if bind.get('historical_manifest'):
                    before=json.loads(path.read_text());after=json.loads((ROOT/name).read_text())
                    for k,value in before.items():
                        if k!='sha256': assert after[k]==value,k
                    for k,value in before['sha256'].items():
                        if k not in old: assert after['sha256'][k]==value,k
                else: assert (ROOT/name).read_bytes().startswith(path.read_bytes())
            assert sha(path)==expected,(label,name);count+=1
    _,actual=mesa(OUT/'selected.data.gz');_,previous=mesa(gr.PROFILE)
    assert all(np.array_equal(actual[k],v) for k,v in previous.items())
    assert len(actual['zone'])==5735 and len(previous)==37
    inputs=json.loads((OUT/'inputs.json').read_text())['sha256']
    old_inputs=json.loads((ROOT/'outputs/thermal-robustness20/mesh-inputs.json').read_text())['sha256']
    changed={'.restart','inlist1','profile_columns.list','photos1/18000'}
    assert all(inputs[k]==v for k,v in old_inputs.items() if k not in changed)
    before=(ROOT/'outputs/thermal-robustness20/runs/mesh/inlist1').read_text()
    assert (OUT/'inputs/inlist1').read_text()==before.replace('max_model_number = 40000','max_model_number = 19042')
    original_columns=(ROOT/'outputs/thermal-robustness20/runs/mesh/profile_columns.list').read_text()
    assert (OUT/'inputs/profile_columns.list').read_text().startswith(original_columns)
    assert json.loads((OUT/'inventory-audit.json').read_text())['passed']
    assert not json.loads((OUT/'inventory-pressure-pilot.json').read_text())['passed']
    previous_gates=json.loads((OLD/'gates.json').read_text())
    assert previous_gates['fixed_temperature_mass_match'] and not previous_gates['fixed_temperature_optical_pass']
    gates=json.loads((OUT/'gates.json').read_text())
    for k in ['original_snapshot_thermal_steady_state','pressure_attached_same_material_interpretation_supported',
              'actual_GR_thermal_composition_evolution_solved','complete_nonlinear_observational_inference','final_PDF_or_ZIP_updated']:
        assert not gates[k],k
    # Imported binary/module integrity is checked separately from source/artifact integrity.
    provenance=json.loads((OLD/'provenance.json').read_text())
    for key in ['runtime_sha256','module_sha256']:
        for path,digest in provenance[key].items(): assert sha(path)==digest,path
    print('PASS:',count,'현재·역사 SHA, 37열 전체 재현, 원래 GR 판정 및 열·물질 실패 경계')


if __name__=='__main__':
    globals()[sys.argv[1]]()
