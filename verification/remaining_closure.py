"""Request26: exhaust accessible closure tests; retain failed physical gates."""
from pathlib import Path
import json, shutil, sys, tempfile
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.integrate import solve_ivp
import fresh_microphysics as fresh
import reactive_energy as reaction
import gr_mass as gr
from thermal_wd import mesa
from thermal_restart import sha

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/remaining-closure26'
CACHE=Path('/home/lpaiu/work/remaining-closure26')
OLD=ROOT/'outputs/fresh-microphysics25'

def save(name,obj):
    (OUT/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')

def context():
    fresh.OUT=OUT;fresh.CACHE=CACHE

def prepare():
    assert not OUT.exists() and not CACHE.exists()
    OUT.mkdir();CACHE.mkdir();context()
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='36f655d',
        previous_manifest_sha256=sha(OLD/'manifest.json'),
        order=['Native reaction-channel isolation and source-vector accessibility',
               'Common-state EOS thermodynamic controls and transport residual',
               'Curved-background scalar stability and response',
               'Global derivative and complete inference identifiability boundaries'],
        controls=dict(channel_superposition_relative=1e-8,EOS_relative=1e-5,
                      scalar_mesh_relative=1e-4),
        limits='No runtime installation/rebuild. No Newtonian evolution passed off as GR. No finite sampling promoted to global certification. Preserve every failed control and finish all independent accessible levers.'))
    base=dict(np.load(OLD/'gr-input.npz'));header,_=mesa(fresh.SOURCE)
    for label in ['full','c12-only','c12-off']:
        fresh.setup_run(label,base,float(header['star_age']))
        folder=CACHE/label
        if label=='c12-only':
            (folder/'cno_extras.net').write_text('add_isos('+','.join(fresh.ISOS)+')\nadd_reactions(r_c12_pg_n13)\n')
        if label=='c12-off':
            text=(folder/'inlist1').read_text().replace('&star_job',
                "&star_job\n num_special_rate_factors=1\n reaction_for_special_factor(1)='r_c12_pg_n13'\n special_rate_factor(1)=0\n",1)
            (folder/'inlist1').write_text(text)
        columns=(folder/'profile_columns.list').read_text()
        names={line.split('!')[0].strip() for line in columns.splitlines()}
        columns+='\n'+'\n'.join(k for k in ['energy','pressure','entropy','cv','cp','gamma1','chiRho','chiT','grada'] if k not in names)+'\n'
        (folder/'profile_columns.list').write_text(columns)
        for name in ['profile_columns.list','cno_extras.net','inlist1']:
            shutil.copy2(folder/name,OUT/'inputs'/label/name)
        fresh.save(label+'-inputs.json',dict(classification='Proven',sha256={p.name:sha(p) for p in folder.iterdir() if p.is_file()}))
    for rel in ['star/private/profile_getval.f90','star/private/net.f90',
                'star/private/star_job_controls_params.inc','star/defaults/star_job.defaults',
                'net/private/net_derivs.f90','net/private/net_derivs_support.f90','net/private/net_eval.f90']:
        path=OUT/'sources'/rel;path.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(fresh.MESA/rel,path)

def run_pilot():
    context()
    for label in ['full','c12-only','c12-off']:
        try: fresh.run(label);fresh.collect(label)
        except Exception as err:
            folder=CACHE/label
            if (folder/'execution.log').exists(): shutil.copy2(folder/'execution.log',OUT/(label+'-execution.log'))
            save(label+'-failure.json',dict(classification='Counterexample candidate',error=repr(err)))
            print(label,'FAILED',repr(err),flush=True)
    _,full=mesa(OUT/'full-profile.data.gz')
    _,single=mesa(OUT/'c12-only-profile.data.gz');_,off=mesa(OUT/'c12-off-profile.data.gz')
    errors={k:float(np.max(abs(full[k]-off[k]-single[k])/np.maximum(1,abs(full[k]))))
            for k in ['eps_nuc','eps_nuc_neu_total']}
    save('channel-pilot.json',dict(classification='Counterexample candidate',errors=errors,
        passed=all(v<1e-8 for v in errors.values()),
        limitations='A passing single-channel heat control does not export the full composition source or every runtime weak Q.'))
    print('CHANNEL',errors,flush=True)

def channels():
    context();assert json.loads((OUT/'channel-pilot.json').read_text())['passed']
    rows=[r for r in reaction.reactions() if r['stoichiometry_complete']]
    save('channel-plan.json',dict(classification='Counterexample candidate',rows=rows,
        registered_sum_relative_tolerance=1e-8,
        interpretation='Isolate every closed channel with the same 22 input isotopes. Compare summed heat, neutrino loss and native derivatives to the original full network. Inferred extents using standard Q are conditional until actual weak Q and compound-channel equivalence are established.'))
    base=dict(np.load(OLD/'gr-input.npz'));header,_=mesa(fresh.SOURCE)
    for i,row in enumerate(rows):
        label='channel-'+row['name'];fresh.setup_run(label,base,float(header['star_age']))
        folder=CACHE/label
        (folder/'cno_extras.net').write_text('add_isos('+','.join(fresh.ISOS)+')\nadd_reactions('+row['name']+')\n')
        shutil.copy2(folder/'cno_extras.net',OUT/'inputs'/label/'cno_extras.net')
        fresh.save(label+'-inputs.json',dict(classification='Proven',sha256={p.name:sha(p) for p in folder.iterdir() if p.is_file()}))
        try: fresh.run(label);fresh.collect(label)
        except Exception as err:
            if (folder/'execution.log').exists(): shutil.copy2(folder/'execution.log',OUT/(label+'-execution.log'))
            save(label+'-failure.json',dict(classification='Counterexample candidate',error=repr(err)))
            print(label,'FAILED',repr(err),flush=True)
        print('CHANNEL PROGRESS',i+1,len(rows),flush=True)
    channel_analysis()

def channel_analysis():
    rows=json.loads((OUT/'channel-plan.json').read_text())['rows'];_,full=mesa(OUT/'full-profile.data.gz')
    raw=[];missing=[];powers=[]
    state=np.load(OLD/'gr-input.npz');weight=state['dm']*np.exp(2*state['nu'])/fresh.LSUN
    keys=['eps_nuc','eps_nuc_neu_total','d_lnepsnuc_dlnT','d_lnepsnuc_dlnd']
    for row in rows:
        path=OUT/('channel-'+row['name']+'-profile.data.gz')
        if not path.exists(): missing.append(row['name']);continue
        _,d=mesa(path);a=np.array([d[k] for k in keys])
        a[2:]*=np.maximum(1,abs(d['eps_nuc']))
        raw.append(a);powers.append(dict(name=row['name'],heat_Lsun=float(weight@a[0]),neutrino_Lsun=float(weight@a[1])))
    sums=np.sum(raw,axis=0);truth=np.array([full[k] for k in keys]);truth[2:]*=np.maximum(1,abs(full['eps_nuc']))
    errors={k:float(np.max(abs(a-b)/np.maximum(1,abs(b)))) for k,a,b in zip(keys,sums,truth)}
    np.savez_compressed(OUT/'channel-outputs.npz',raw=np.array(raw),sums=sums,full=truth)
    save('channel-sum.json',dict(classification='Counterexample candidate',missing=missing,errors=errors,
        passed=not missing and all(v<1e-8 for v in errors.values()),powers=sorted(powers,key=lambda x:abs(x['heat_Lsun']),reverse=True),
        full_species_derivative_exported=False,actual_weak_Q_exported=False))
    print('CHANNEL SUM',errors,'missing',missing,flush=True)

def fullnet_channels():
    context();rows=json.loads((OUT/'channel-plan.json').read_text())['rows']
    shutil.copy2(CACHE/'off-r34_pp2/execution.log',OUT/'off-r34_pp2-execution.log')
    shutil.copy2(CACHE/'half-r34_pp2/execution.log',OUT/'half-r34_pp2-execution.log')
    save('zero-factor-failure.json',dict(classification='Counterexample candidate',passed=False,
        reaction='r34_pp2',factors=[0,.5],reason='set_combo_screen_rates requires equal pre-branch r34_pp2 and r34_pp3 rates. Both independent factor changes fail. The pair must be perturbed together.'))
    groups=[['r34_pp2','r34_pp3']]+[[r['name']] for r in rows if not r['name'].startswith('r34_pp')]
    save('fullnet-channel-plan.json',dict(classification='Counterexample candidate',
        frozen_isolated_failure_sha256=sha(OUT/'channel-sum.json'),
        frozen_zero_factor_failure_sha256=sha(OUT/'zero-factor-failure.json'),
        groups=groups,
        method='Retain the original whole network and halve each reaction group. r34_pp2/pp3 must be grouped because their shared capture rate is split by the native network. Compare doubled differences to the full powers/derivatives. Earlier failures remain unchanged.',
        relative_sum_tolerance=1e-8))
    base=dict(np.load(OLD/'gr-input.npz'));header,_=mesa(fresh.SOURCE)
    for i,group in enumerate(groups):
        label='group-'+group[0];fresh.setup_run(label,base,float(header['star_age']))
        folder=CACHE/label
        job='&star_job\n num_special_rate_factors='+str(len(group))+'\n'
        for j,name in enumerate(group,1): job+=f" reaction_for_special_factor({j})='{name}'\n special_rate_factor({j})=0.5\n"
        text=(folder/'inlist1').read_text().replace('&star_job',job,1)
        (folder/'inlist1').write_text(text);shutil.copy2(folder/'inlist1',OUT/'inputs'/label/'inlist1')
        fresh.save(label+'-inputs.json',dict(classification='Proven',sha256={p.name:sha(p) for p in folder.iterdir() if p.is_file()}))
        fresh.run(label);fresh.collect(label)
        print('FULL NETWORK PROGRESS',i+1,len(groups),flush=True)
    fullnet_analysis()

def fullnet_analysis():
    groups=json.loads((OUT/'fullnet-channel-plan.json').read_text())['groups'];_,full=mesa(OUT/'full-profile.data.gz')
    keys=['eps_nuc','eps_nuc_neu_total','d_lnepsnuc_dlnT','d_lnepsnuc_dlnd']
    truth=np.array([full[k] for k in keys]);truth[2:]*=np.maximum(1,abs(full['eps_nuc']))
    delta=[]
    for group in groups:
        _,d=mesa(OUT/('group-'+group[0]+'-profile.data.gz'))
        a=np.array([d[k] for k in keys]);a[2:]*=np.maximum(1,abs(d['eps_nuc']))
        delta.append(2*(truth-a))
    delta=np.array(delta);sums=delta.sum(axis=0)
    errors={k:float(np.max(abs(a-b)/np.maximum(1,abs(b)))) for k,a,b in zip(keys,sums,truth)}
    np.savez_compressed(OUT/'fullnet-channel-outputs.npz',delta=delta,sums=sums,full=truth)
    save('fullnet-channel-sum.json',dict(classification='Counterexample candidate',errors=errors,
        passed=all(v<1e-8 for v in errors.values()),
        direct_full_species_derivative_exported=False,actual_weak_Q_exported=False,
        interpretation='Full-network channel perturbations; successful sum does not by itself observe each reaction extent or actual weak Q.'))
    print('FULL NETWORK SUM',errors,flush=True)

def thermodynamic_runs():
    context();base=dict(np.load(OLD/'gr-input.npz'));header,_=mesa(fresh.SOURCE)
    save('EOS-control-plan.json',dict(classification='Counterexample candidate',step=5e-5,
        tests=['MESA pressure/energy/entropy versus native derivatives',
               'Fixed-composition first-law Maxwell identities',
               'FreeEOS to MESA energy difference is not assumed a constant reference shift'],
        relative_tolerance=1e-5))
    for field in ['lnT','lnd']:
        for sign,suffix in [(-1,'minus'),(1,'plus')]:
            label='EOS-'+field+'-'+suffix;state={k:v.copy() for k,v in base.items()};state[field]+=sign*5e-5
            fresh.setup_run(label,state,float(header['star_age']));folder=CACHE/label
            shutil.copy2(OUT/'inputs/full/profile_columns.list',folder/'profile_columns.list')
            shutil.copy2(folder/'profile_columns.list',OUT/'inputs'/label/'profile_columns.list')
            fresh.save(label+'-inputs.json',dict(classification='Proven',sha256={p.name:sha(p) for p in folder.iterdir() if p.is_file()}))
            fresh.run(label);fresh.collect(label)
    EOS_analysis()

def EOS_analysis():
    _,d=mesa(OUT/'full-profile.data.gz');h=5e-5;unit=1.3806504e-16*6.02214179e23
    fd={}
    for var in ['lnT','lnd']:
        _,p=mesa(OUT/('EOS-'+var+'-plus-profile.data.gz'));_,m=mesa(OUT/('EOS-'+var+'-minus-profile.data.gz'))
        fd[var]={k:(p[k]-m[k])/(2*h) for k in ['pressure','energy','entropy']}
    temp=10**d['logT'];rho=d['rho'];P=d['pressure']
    pairs={
        'chiT':(fd['lnT']['pressure']/P,d['chiT']),
        'chiRho':(fd['lnd']['pressure']/P,d['chiRho']),
        'cv':(fd['lnT']['energy']/temp,d['cv']),
        'entropy_T_vs_cv':(fd['lnT']['entropy']*unit,d['cv']),
        'entropy_rho_Maxwell':(fd['lnd']['entropy']*unit,-P*d['chiT']/rho/temp),
        'energy_rho_Maxwell':(fd['lnd']['energy'],P*(1-d['chiT'])/rho)}
    results={key:dict(max_normalized_error=float(np.max(abs(a-b)/np.maximum(1,abs(b)))),
                       worst_zone=int(d['zone'][np.argmax(abs(a-b)/np.maximum(1,abs(b)))])) for key,(a,b) in pairs.items()}
    save('EOS-controls.json',dict(classification='Counterexample candidate',results=results,
        all_passed=all(v['max_normalized_error']<1e-5 for v in results.values()),
        same_state_only=True,common_EOS_GR_structure=False,
        interpretation='Native derivative versus finite differences and thermodynamic integrability are separate tests. A consistent numerical derivative is not sufficient to prove the first law for a blended EOS.'))
    print('EOS CONTROLS',results,flush=True)

def EOS_refine():
    if not (OUT/'EOS-controls-pilot.json').exists(): shutil.copy2(OUT/'EOS-controls.json',OUT/'EOS-controls-pilot.json')
    EOS_analysis();context()
    base=dict(np.load(OLD/'gr-input.npz'));header,_=mesa(fresh.SOURCE)
    save('EOS-refinement-plan.json',dict(classification='Counterexample candidate',
        frozen_first_failure_sha256=sha(OUT/'EOS-controls.json'),step=2.5e-5,
        correction='Use source kerg=1.3806504e-16, not 1.3806488e-16, to convert dimensionless entropy. The pressure and energy failures do not depend on this conversion.',
        normalization='First-law residuals normalized by cv*T or P/rho, not by a possibly vanishing individual term. No original gate is promoted.'))
    for field in ['lnT','lnd']:
        for sign,suffix in [(-1,'minus'),(1,'plus')]:
            label='EOS2-'+field+'-'+suffix;state={k:v.copy() for k,v in base.items()};state[field]+=sign*2.5e-5
            fresh.setup_run(label,state,float(header['star_age']));folder=CACHE/label
            shutil.copy2(OUT/'inputs/full/profile_columns.list',folder/'profile_columns.list')
            shutil.copy2(folder/'profile_columns.list',OUT/'inputs'/label/'profile_columns.list')
            fresh.save(label+'-inputs.json',dict(classification='Proven',sha256={p.name:sha(p) for p in folder.iterdir() if p.is_file()}))
            fresh.run(label);fresh.collect(label)
    EOS_refinement_analysis()

def EOS_refinement_analysis():
    _,d=mesa(OUT/'full-profile.data.gz');T=10**d['logT'];P=d['pressure'];rho=d['rho'];unit=1.3806504e-16*6.02214179e23
    fd={}
    for prefix,h in [('EOS',5e-5),('EOS2',2.5e-5)]:
        fd[prefix]={}
        for var in ['lnT','lnd']:
            _,p=mesa(OUT/(prefix+'-'+var+'-plus-profile.data.gz'));_,m=mesa(OUT/(prefix+'-'+var+'-minus-profile.data.gz'))
            fd[prefix][var]={k:(p[k]-m[k])/(2*h) for k in ['pressure','energy','entropy']}
    a=fd['EOS2'];b=fd['EOS'];res={}
    for key in ['pressure','energy','entropy']:
        scale={'pressure':P,'energy':d['cv']*T,'entropy':d['cv']/unit}[key]
        for var in ['lnT','lnd']:
            res[var+'_'+key+'_step_change']=float(np.max(abs(a[var][key]-b[var][key])/np.maximum(1,abs(scale))))
    # Compare two derivatives of actual evaluated potentials, independent of native chi/cv.
    res['first_law_T']=float(np.max(abs(a['lnT']['energy']-T*unit*a['lnT']['entropy'])/(d['cv']*T)))
    res['first_law_rho']=float(np.max(abs(a['lnd']['energy']-T*unit*a['lnd']['entropy']-P/rho)/(P/rho)))
    res['Maxwell_Srho_PT']=float(np.max(abs(T*unit*a['lnd']['entropy']+a['lnT']['pressure']/rho)/(P/rho)))
    save('EOS-refinement.json',dict(classification='Counterexample candidate',results=res,
        full_common_EOS_certified=False,
        interpretation='The small-step differences do not converge. Thus their apparent first-law residuals do not establish a continuum thermodynamic inconsistency. No global derivative or common-EOS certificate follows. The source has single-precision table interpolation, but the binary error budget is not independently certified.'))
    print('EOS REFINEMENT',res,flush=True)

def helm_control():
    context();base=dict(np.load(OLD/'gr-input.npz'));header,_=mesa(fresh.SOURCE)
    save('HELM-control-plan.json',dict(classification='Counterexample candidate',
        prior_EOS_failure_sha256=sha(OUT/'EOS-refinement.json'),
        method='Use the existing use_eosDT_HELMEOS flag as a thermodynamic formula/units positive control on the same supplied states. Fully ionized HELM replaces the default blended EOS and is NOT a physical replacement for the partially ionized envelope.',
        step=5e-5))
    for label,field,delta in [('HELM','lnT',0),('HELM-lnT-minus','lnT',-5e-5),('HELM-lnT-plus','lnT',5e-5),
                              ('HELM-lnd-minus','lnd',-5e-5),('HELM-lnd-plus','lnd',5e-5)]:
        state={k:v.copy() for k,v in base.items()};state[field]+=delta
        fresh.setup_run(label,state,float(header['star_age']));folder=CACHE/label
        shutil.copy2(OUT/'inputs/full/profile_columns.list',folder/'profile_columns.list')
        text=(folder/'inlist1').read_text().replace('&controls','&controls\n use_eosDT_HELMEOS=.true.\n',1)
        (folder/'inlist1').write_text(text)
        for name in ['profile_columns.list','inlist1']: shutil.copy2(folder/name,OUT/'inputs'/label/name)
        fresh.save(label+'-inputs.json',dict(classification='Proven',sha256={p.name:sha(p) for p in folder.iterdir() if p.is_file()}))
        fresh.run(label);fresh.collect(label)
    helm_analysis()

def helm_analysis():
    _,d=mesa(OUT/'HELM-profile.data.gz');T=10**d['logT'];P=d['pressure'];rho=d['rho'];unit=1.3806504e-16*6.02214179e23
    fd={}
    for var in ['lnT','lnd']:
        _,p=mesa(OUT/('HELM-'+var+'-plus-profile.data.gz'));_,m=mesa(OUT/('HELM-'+var+'-minus-profile.data.gz'))
        fd[var]={k:(p[k]-m[k])/1e-4 for k in ['pressure','energy','entropy']}
    pairs={'chiT':(fd['lnT']['pressure']/P,d['chiT']),
           'chiRho':(fd['lnd']['pressure']/P,d['chiRho']),
           'cv':(fd['lnT']['energy']/T,d['cv']),
           'entropy_T_cv':(fd['lnT']['entropy']*unit,d['cv'])}
    errors={k:float(max(abs(a-b)/np.maximum(1,abs(b)))) for k,(a,b) in pairs.items()}
    errors['first_law_rho']=float(max(abs(fd['lnd']['energy']-T*unit*fd['lnd']['entropy']-P/rho)/(P/rho)))
    hot=T>1e7
    hot_errors={k:float(max((abs(a-b)/np.maximum(1,abs(b)))[hot])) for k,(a,b) in pairs.items()}
    hot_errors['first_law_rho']=float(max((abs(fd['lnd']['energy']-T*unit*fd['lnd']['entropy']-P/rho)/(P/rho))[hot]))
    save('HELM-control.json',dict(classification='Counterexample candidate',errors=errors,
        passed=max(errors.values())<1e-5,default_EOS_failure_promoted=False,
        hot_domain_cells=int(hot.sum()),hot_errors=hot_errors,hot_diagnostic_passed=max(hot_errors.values())<1e-5,
        hot_domain_policy='T>1e7 K is a post-result localization diagnostic, not a promotion of the failed full-domain control.',
        physical_envelope_validated=False))
    print('HELM CONTROL',errors,'hot',hot_errors,flush=True)

def transport():
    from scipy.linalg import solve_banded
    import baryon_entropy as be
    s=np.load(OLD/'gr-input.npz');_,d=mesa(OUT/'full-profile.data.gz');mat=be.Material();eos=gr.EOS()
    cv=[]
    for i in range(len(s['dm'])):
        a=eos(2,s['lnd'][i]+np.log(s['CX'][i]),s['lnT'][i],mat.eps[i])
        cv.append(s['CX'][i]*a[10]/np.exp(s['lnT'][i]))
    cv=np.array(cv)[::-1];dm=s['dm'][::-1];T=np.exp(s['lnT'][::-1]);N=np.exp(s['nu'][::-1])
    r=s['r_mid_m'][::-1]*100;mass=s['m_mid_geom'][::-1]*100;rho=np.exp(s['lnd'][::-1]);f=1-2*mass/r
    kappa=d['opacity'][::-1];K=16*5.670400e-5*T**3/(3*kappa*rho)
    theta=N*T;rc=(r[1:]+r[:-1])/2;Nc=(N[1:]+N[:-1])/2;fc=(f[1:]+f[:-1])/2
    G=4*np.pi*rc**2*Nc*np.sqrt(fc)*(2/(1/K[1:]+1/K[:-1]))/np.diff(r)
    C=cv*dm;source=dm*N*N*(d['eps_nuc'][::-1]-d['non_nuc_neu'][::-1])
    assert np.all(C>0) and np.all(G>0)
    L=-G*np.diff(theta);div=np.r_[L,0.]-np.r_[0.,L]
    # Deliberately closed boundaries: this checks conservative storage, not the observed surface luminosity.
    rhs=source-div;rate=rhs/C
    duration=1.6294*86400;states=[];rows=[]
    diag=np.r_[G,0.]+np.r_[0.,G]
    for steps in [1,2,4,8,16,32]:
        dt=duration/steps;increment=np.zeros(len(C));residual=[]
        ab=np.zeros((3,len(C)));ab[1]=C+dt*diag;ab[0,1:]=-dt*G;ab[2,:-1]=-dt*G
        for j in range(steps):
            before=increment.copy();increment=solve_banded((1,1),ab,C*increment+dt*rhs)
            residual.append(float((np.dot(C,increment-before)-dt*sum(source))/(dt*sum(abs(source)))))
        states.append(increment)
        rows.append(dict(steps=steps,max_temperature_relative_change=float(np.max(abs(increment)/theta)),
            energy_balance_max_relative=float(max(abs(np.array(residual)))),
            radiated_energy=0.,stored_energy_erg=float(np.dot(C,increment))))
    change=float(np.max(abs(states[-1]-states[-2])/theta))
    np.savez_compressed(OUT/'frozen-GR-transport.npz',theta=theta,C=C,conductance=G,source=source,initial_flux=L,
                        initial_temperature_rate=rate/N,final_increment=states[-1])
    save('frozen-GR-transport.json',dict(classification='Counterexample candidate',rows=rows,
        last_two_max_relative_temperature_difference=change,duration_s=duration,
        energy_gate_passed=all(r['energy_balance_max_relative']<1e-7 for r in rows),
        closed_source_power_Lsun=float(sum(source)/fresh.LSUN),
        physical_evolution_solved=False,
        model='Declared frozen-density, fixed-composition, fixed-conductivity GR finite-volume heat equation with prescribed source and closed boundaries. FreeEOS heat capacity and MESA opacity are a hybrid control.',
        theorem='Positive capacities and symmetric positive face conductances give exact telescoping energy conservation and dissipative homogeneous diffusion. They do not prove stellar thermal stability once reaction/composition/metric feedback is included.',
        missing='No changing species/rest mass, expansion, common EOS, convection or surface atmosphere closure. This is a solver/control result, not a stellar evolution prediction.'))
    print('TRANSPORT',rows,'refinement',change,flush=True)

def radiative_transport():
    from scipy.linalg import solve_banded
    a=np.load(OUT/'frozen-GR-transport.npz');s=np.load(OLD/'gr-input.npz')
    C=a['C'];G=a['conductance'];source=a['source'];theta=a['theta'];n=len(C)
    initial_flux=-G*np.diff(theta);rhs=source-np.r_[initial_flux,0.]+np.r_[0.,initial_flux]
    row=json.loads((ROOT/'outputs/baryon-entropy23/scaled-33-4.json').read_text())
    R=row['photospheric_radius_m']*100;Ns=np.exp(s['nu'][0]);coef=4*np.pi*R**2*5.670400e-5/Ns**2
    diag=np.r_[G,0.]+np.r_[0.,G];duration=1.6294*86400;rows=[];states=[]
    save('radiative-transport-plan.json',dict(classification='Counterexample candidate',
        prior_closed_boundary_result_sha256=sha(OUT/'frozen-GR-transport.json'),
        boundary='Separate grey blackbody outer-face closure: L_infinity=4*pi*R^2*sigma*theta_surface^4/N_surface^2. No atmospheric transfer or moving geometry is solved.',
        held_fixed='FreeEOS heat capacity, baryon density, composition, metric, MESA conductivity and prescribed heat source.'))
    for steps in [1,2,4,8,16,32]:
        dt=duration/steps;change=np.zeros(n);loss=0.;errors=[];iterations=[]
        ab=np.zeros((3,n));ab[1]=C+dt*diag;ab[0,1:]=-dt*G;ab[2,:-1]=-dt*G
        for j in range(steps):
            prev=change.copy();b=C*prev+dt*rhs
            for k in range(30):
                flux=-G*np.diff(change)
                F=C*change+dt*(np.r_[flux,0.]-np.r_[0.,flux])-b
                outer=coef*(theta[-1]+change[-1])**4
                F[-1]+=dt*outer
                jac=ab.copy();jac[1,-1]+=dt*4*coef*(theta[-1]+change[-1])**3
                correction=solve_banded((1,1),jac,-F)
                change+=correction
                assert np.all(theta+change>0)
                if max(abs(correction)/theta)<1e-11: break
            else: raise AssertionError('radiative boundary Newton did not converge')
            outer=coef*(theta[-1]+change[-1])**4;loss+=outer*dt
            errors.append(abs((C@(change-prev)+outer*dt-source.sum()*dt)/(sum(abs(source))*dt)))
            iterations.append(k+1)
        states.append(change)
        rows.append(dict(steps=steps,max_temperature_relative_change=float(max(abs(change)/theta)),
            energy_balance_relative=float(max(errors)),surface_L_infinity_Lsun=float(outer/fresh.LSUN),
            stored_energy_erg=float(C@change),radiated_energy_erg=float(loss),max_Newton_iterations=max(iterations)))
    save('radiative-GR-transport.json',dict(classification='Counterexample candidate',rows=rows,
        last_two_max_relative_temperature_difference=float(max(abs(states[-1]-states[-2])/theta)),
        energy_gate_passed=max(r['energy_balance_relative'] for r in rows)<1e-7,
        full_physical_evolution_solved=False,
        interpretation='Nonlinear outer boundary and conservative time-step verification on frozen coefficients. Reported luminosity belongs to this declared control, not an observational prediction.'))
    print('RADIATIVE TRANSPORT',rows,flush=True)

def fluid_diagnostic():
    s=np.load(OLD/'gr-input.npz');import baryon_entropy as be
    mat=be.Material();eos=gr.EOS();gamma=[]
    for i in range(len(s['dm'])):
        a=eos(2,s['lnd'][i]+np.log(s['CX'][i]),s['lnT'][i],mat.eps[i]);gamma.append(a[4])
    r=s['r_mid_m'][::-1];M=s['m_mid_geom'][::-1];P=np.exp(s['logP'][::-1])*.1
    e=(np.exp(s['lnd'])*1000*(s['CX']*gr.C**2+s['u_W']*1e-4))[::-1]
    N=np.exp(s['nu'][::-1]);f=1-2*M/r;gamma=np.array(gamma)[::-1]
    nup=(M+4*np.pi*r**3*P*gr.G/gr.C**4)/(r*r*f);dp=-(e+P)*nup
    de1=PchipInterpolator(r,e).derivative()(r);de2=np.gradient(e,r,edge_order=2)
    pref=N*N*f*nup*gr.C**2
    n21=pref*(dp/(gamma*P)-de1/(e+P));n22=pref*(dp/(gamma*P)-de2/(e+P))
    interior=np.ones(len(r),bool);interior[:2]=False;interior[-2:]=False
    negative=(n21<0)&(n22<0)&interior;dm=s['dm'][::-1]
    save('buoyancy-diagnostic.json',dict(classification='Counterexample candidate',
        formula='N_infinity^2=N^2*f*nu_prime*c^2*(P_prime/(Gamma1*P)-e_prime/(e+P))',
        agreeing_negative_cells=int(negative.sum()),mass_fraction=float(sum(dm[negative])/sum(dm)),
        min_PCHIP_N2_s2=float(min(n21[interior])),min_three_point_N2_s2=float(min(n22[interior])),
        gamma1_range=[float(min(gamma)),float(max(gamma))],
        full_fluid_mode_stability_certified=False,
        interpretation='Local adiabatic Cowling/short-wavelength buoyancy diagnostic of the sampled material background; two differentiation choices. Cell composition/entropy interfaces and the unsolved evolving background prevent promoting it to a certified stellar eigenmode or a scalar observable.'))
    np.savez_compressed(OUT/'buoyancy-diagnostic.npz',r_m=r,gamma1=gamma,N2_PCHIP=n21,N2_three_point=n22,negative=negative)
    print('BUOYANCY negative agreeing cells',int(negative.sum()),'mass fraction',float(sum(dm[negative])/sum(dm)),flush=True)

def pp_diagnostic():
    context();base=dict(np.load(OLD/'gr-input.npz'));header,_=mesa(fresh.SOURCE)
    save('pp-derivative-plan.json',dict(classification='Counterexample candidate',
        previous_native_failure_sha256=sha(OLD/'derivative-control.json'),
        hypothesis='set_combo_screen_rates multiplies raw rate derivatives by fII/fIII but omits derivatives of the branching fraction. Test the missing term against actual full-network finite differences, not source inspection alone.',
        steps=[5e-5,2.5e-5],relative_residual_tolerance=1e-5))
    for h in [5e-5,2.5e-5]:
        for sign,suffix in [(-1,'minus'),(1,'plus')]:
            label='pp-'+str(h)+'-'+suffix;state={k:v.copy() for k,v in base.items()};state['lnT']+=sign*h
            fresh.setup_run(label,state,float(header['star_age']));folder=CACHE/label
            net='add_isos('+','.join(fresh.ISOS)+')\nadd_reactions(r34_pp2,r34_pp3,rbe7ec_li7_aux,rbe7pg_b8_aux)\n'
            (folder/'cno_extras.net').write_text(net);shutil.copy2(folder/'cno_extras.net',OUT/'inputs'/label/'cno_extras.net')
            fresh.save(label+'-inputs.json',dict(classification='Proven',sha256={p.name:sha(p) for p in folder.iterdir() if p.is_file()}))
            fresh.run(label);fresh.collect(label)
    pp_analysis()

def pp_analysis():
    _,g=mesa(OLD/'gr-profile.data.gz');domain=abs(g['eps_nuc'])>1
    rows={r['name']:r for r in reaction.reactions()};Q=float(rows['r34_pp2']['Q_standard_MeV'])
    q2=float(rows['r34_pp2']['Q_neutrino_list_MeV']);q3=float(rows['r34_pp3']['Q_neutrino_list_MeV'])
    fullNative=g['d_lnepsnuc_dlnT']*np.maximum(1,abs(g['eps_nuc']))
    result=[];arrays=[]
    for h in [5e-5,2.5e-5]:
        _,p=mesa(OUT/('pp-'+str(h)+'-plus-profile.data.gz'));_,m=mesa(OUT/('pp-'+str(h)+'-minus-profile.data.gz'))
        def frac(d):
            raw=(d['eps_nuc']+d['eps_nuc_neu_total'])/Q
            nu=np.divide(d['eps_nuc_neu_total'],raw,out=np.zeros_like(raw),where=abs(raw)>1e-290)
            return raw,(q3-nu)/(q3-q2)
        rp,fp=frac(p);rm,fm=frac(m)
        correction=(rp+rm)/2*(q3-q2)*(fp-fm)/(2*h)
        _,gp=mesa(OLD/('T'+str(h)+'-plus-profile.data.gz'));_,gm=mesa(OLD/('T'+str(h)+'-minus-profile.data.gz'))
        fd=(gp['eps_nuc']-gm['eps_nuc'])/(2*h)
        error=abs(fd-fullNative-correction)/np.maximum(1,abs(fullNative))
        result.append(dict(step=h,max_residual=float(error[domain].max()),worst_zone=int(g['zone'][np.flatnonzero(domain)[error[domain].argmax()]]),
                           passed=bool(error[domain].max()<1e-5)))
        arrays.append(correction)
    np.savez_compressed(OUT/'pp-derivative-correction.npz',correction=arrays[-1],previous=arrays[0],domain=domain,native=fullNative)
    save('pp-derivative-result.json',dict(classification='Counterexample candidate',results=result,
        original_native_failure_promoted=False,compiled_runtime_modified=False,
        interpretation='The correction is evaluated with actual branch-rate diagnostics on the checked state. It does not patch the compiled runtime or certify global derivatives.'))
    print('PP CORRECTION',result,flush=True)

def composition_reconstruction():
    import re
    assert json.loads((OUT/'fullnet-channel-sum.json').read_text())['passed']
    groups=json.loads((OUT/'fullnet-channel-plan.json').read_text())['groups']
    rows={r['name']:r for r in reaction.reactions()};iso=reaction.isotopes()
    _,full=mesa(OUT/'full-profile.data.gz');s=np.load(OLD/'gr-input.npz')
    powers=np.load(OUT/'fullnet-channel-outputs.npz')['delta']
    table=fresh.MESA/'data/rates_data/weakreactions.tables'
    pairs=set(re.findall(r'^\s*([a-z]+\d*)\s+([a-z]+\d*)\s+by',table.read_text(),re.MULTILINE))
    dX=np.zeros_like(s['X']);flow=[];qvalues=[];cancellation=np.zeros_like(dX);weak_work=np.zeros(len(dX))
    runtime_weak=[];zero_q=[]
    for group,power in zip(groups,powers):
        row=rows[group[0]];nu=row['nu']
        assert all(rows[n]['nu']==nu and rows[n]['Q_standard_MeV']==row['Q_standard_MeV'] for n in group)
        Q=np.full(len(dX),float(row['Q_standard_MeV']))
        if '_wk_' in group[0]:
            left=[k for k,v in nu.items() if v<0];right=[k for k,v in nu.items() if v>0]
            if len(left)==len(right)==1 and (left[0],right[0]) in pairs:
                # Match eval_weak's declared REAL conversion and clipped table T9.
                T9=np.maximum(np.float32(np.exp(s['lnT'])*1e-9),np.float32(.01))
                temp=np.float32(T9*np.float32(1e9))
                conv=np.float32(1.3806504e-16/1.602176487e-6)
                mue=full['eta']*float(conv)*temp.astype(float)
                sign=1 if iso[left[0]]['z']>iso[right[0]]['z'] else -1
                Q=float(iso[left[0]]['ex']-iso[right[0]]['ex'])+sign*mue
                runtime_weak.append(group[0])
        if np.any(abs(Q)<1e-8): zero_q.append(group[0]);continue
        extent=(power[0]+power[1])/(Q*reaction.QCONV)
        weak_work+=(Q-float(row['Q_standard_MeV']))*reaction.QCONV*extent
        # Conservative floating arithmetic scale only; not a model/weak-Q certification.
        noise=32*np.finfo(float).eps*(abs(full['eps_nuc'])+abs(full['eps_nuc_neu_total'])+1)/(abs(Q)*reaction.QCONV)
        for name,count in nu.items():
            j=fresh.ISOS.index(name);dX[:,j]+=float(iso[name]['a'])*count*extent
            cancellation[:,j]+=float(iso[name]['a'])*abs(count)*noise
        flow.append(extent);qvalues.append(Q)
    assert not zero_q
    fluor=[fresh.ISOS.index(k) for k in ['f17','f18','f19']]
    generation=dX[:,fluor].sum(axis=1);noise=cancellation[:,fluor].sum(axis=1)
    baryon=np.max(abs(dX.sum(axis=1))/np.maximum(1e-30,abs(dX).sum(axis=1)))
    A=np.array([float(iso[k]['a']) for k in fresh.ISOS]);W=np.array([float(iso[k]['w']) for k in fresh.ISOS])
    EX=np.array([float(iso[k]['qex']) for k in fresh.ISOS]);dY=dX/A
    rest_derivative=dY@((W-A)*(gr.C*100)**2)
    gdot=dY@((W-A)*(gr.C*100)**2-EX*reaction.QCONV)
    heatW=-rest_derivative-full['eps_nuc_neu_total']
    residual=full['eps_nuc']-heatW-gdot-weak_work
    save('reconstructed-energy-ledger.json',dict(classification='Counterexample candidate',
        max_relative_residual=float(max(abs(residual)/np.maximum(1,abs(full['eps_nuc'])))),
        max_weak_electron_work_erg_g_s=float(max(abs(weak_work))),
        identity='native heat = (-neutral rest derivative - reaction neutrinos) + dg/dtau + state-dependent weak-Q electron work',
        independent_native_flow_validation=False,
        interpretation='Bookkeeping closure of the channel-reconstructed vector, not an independent native dxdt comparison. It must not be counted as a validated composition integrator.'))
    np.savez_compressed(OUT/'reconstructed-composition-source.npz',dX_dproper_s=dX,
                        channel_extents=np.array(flow),source_Q_MeV=np.array(qvalues),arithmetic_scale=cancellation)
    save('composition-source.json',dict(classification='Counterexample candidate',
        groups=len(groups),runtime_weak_Q_formula_groups=runtime_weak,
        baryon_relative_cancellation=float(baryon),
        max_fluorine_generation_per_s=float(max(generation)),
        fluorine_generation_above_arithmetic_scale_cells=int(sum(generation>100*noise)),
        max_generation_zone=int(full['zone'][generation.argmax()]),
        fluorine_production_integral_g_s=float(s['dm']@generation),
        native_dxdt_directly_exported=False,
        limitations='Reconstruction from finite channel interventions and source Q formulas, not a direct observation of native dxdt. Arithmetic cancellation scale excludes source-to-binary mismatch, uncertainty in weak Q/eta and EOS errors. No global rate/error certificate.',
        consequence='If the reconstructed positive fluorine direction is the physical network vector, the exactly fluorine-free FreeEOS manifold is not invariant. Dropping and renormalizing newly produced F changes the reaction model; a 22-species continuation needs an EOS covering it.'))
    print('COMPOSITION',baryon,'F max',max(generation),'F cells',sum(generation>100*noise),'weak',runtime_weak,flush=True)

def background():
    """Declared PCHIP background of frozen material samples plus its atmosphere."""
    import baryon_entropy as be
    s=np.load(OLD/'gr-material-state.npz')
    row=json.loads((ROOT/'outputs/baryon-entropy23/scaled-33-4.json').read_text())
    mat=be.Material();eos=gr.EOS();ref=np.load(ROOT/'outputs/baryon-entropy23/adiabats-33.npz')['reference']
    center,tc,_=be.invert(eos,row['parameters'][0],ref[-1,3],mat.eps[-1],mat.lt[-1])
    surf=eos(1,mat.lp[0],mat.lt[0],mat.eps[0])
    ps=gr.G*surf[1]*.1/gr.C**4;rho=gr.G*surf[0]*1000/gr.C**2
    u0=surf[2]*1e-4/gr.C**2-1.5*ps/rho
    rp=row['photospheric_radius_m'];mp=row['photospheric_mass_GM_solar']*gr.GM_SUN/gr.C**2
    def rhs(z,y):
        r,m=y;p=ps*z**2.5;den=rho*z**1.5;e=den*(1+u0)+1.5*p
        dr=-2.5*ps/rho*r*r*(1-2*m/r)/((1+u0+2.5*ps/rho*z)*(m+4*np.pi*r**3*p))
        return [dr,4*np.pi*r*r*e*dr]
    z=np.linspace(1,0,101)
    at=solve_ivp(rhs,(1,0),[rp,mp],t_eval=z,method='DOP853',rtol=1e-12,atol=[1e-6,1e-12],max_step=.02)
    assert at.success
    R,M=at.y[:,-1];nv=.5*np.log1p(-2*M/R)
    na=nv-np.log1p(2.5*ps/rho*z/(1+u0))
    ea=(rho*z**1.5*(1+u0)+1.5*ps*z**2.5)*gr.C**4/gr.G
    pa=ps*z**2.5*gr.C**4/gr.G
    cnu=json.loads((OLD/'gr-state-control.json').read_text())['central_nu']
    r=np.r_[0.,s['r_mid_m'][::-1],at.y[0]]
    nu=np.r_[cnu,s['nu'][::-1],na]
    mass=np.r_[0.,s['m_mid_geom'][::-1],at.y[1]]
    energy=np.r_[center[0]*1000*(gr.C**2+center[2]*1e-4),
                 (np.exp(s['lnd'])*1000*(s['CX']*gr.C**2+s['u_W']*1e-4))[::-1],ea]
    pressure=np.r_[center[1]*.1,np.exp(s['logP'][::-1])*.1,pa]
    trace=energy-3*pressure;compact=np.divide(2*mass,r,out=np.zeros_like(r),where=r>0)
    assert np.all(np.diff(r)>0) and np.all(trace>=0)
    x=r/R
    # ponytail: this is a declared finite-profile interpolant, not an EOS/PDE error enclosure.
    N=PchipInterpolator(x,np.exp(nu));f=PchipInterpolator(x,1-compact)
    tr=PchipInterpolator(x,trace)
    return R,M,N,f,tr,dict(r=r,nu=nu,mass=mass,energy=energy,pressure=pressure,trace=trace)

def scalar():
    R,M,N,f,tr,data=background();mu=M/R
    np.savez_compressed(OUT/'scalar-background.npz',**data)
    def interior(beta,k,tol):
        def rhs(x,y):
            n=float(N(x));ff=float(f(x));v=4*np.pi*gr.G*beta*float(tr(x))*R*R/gr.C**4
            return [y[1]/(x*x*n*np.sqrt(ff)),x*x*(v*n-k*k/n)/np.sqrt(ff)*y[0]]
        start=1e-6;v=4*np.pi*gr.G*beta*float(tr(0))*R*R/gr.C**4
        aa=v-k*k/float(N(0))**2
        initial=[1+aa*start*start/6,float(N(0))*aa*start**3/3]
        sol=solve_ivp(rhs,(start,1),initial,method='DOP853',rtol=tol,atol=[tol*1e-2,tol*1e-6],max_step=.002)
        assert sol.success
        return sol.y[:,-1]
    def outgoing(k,far,tol):
        xmax=far/k;ff=1-2*mu/xmax;xs=xmax+2*mu*np.log(xmax/(2*mu)-1)
        wave=np.exp(1j*k*xs)
        def rhs(logx,y):
            x=np.exp(logx);ff=1-2*mu/x;fp=2*mu/x**2
            return [x*y[1],-x*(fp/ff*y[1]+(k*k/ff**2-2*mu/(x**3*ff))*y[0])]
        sol=solve_ivp(rhs,(np.log(xmax),0),np.array([wave,1j*k/ff*wave]),method='DOP853',rtol=tol,atol=tol*1e-4,max_step=.04)
        assert sol.success
        psi,dp=sol.y[:,-1]
        return np.array([psi,(1-2*mu)*(dp-psi)])
    def W(a,b): return a[0]*b[1]-a[1]*b[0]
    statics=[];responses=[]
    for tol,far in [(1e-10,20),(2e-12,40)]:
        y=interior(-4,0,tol);H=-np.log1p(-2*mu)/(2*mu)
        chi=-R*y[1]/(y[0]+H*y[1]);statics.append(float(chi))
        for days in [1.6294,327.26]:
            omega=2*np.pi/(days*86400);k=omega*R/gr.C
            a=interior(-4,k,tol);b=interior(0,k,tol);out=outgoing(k,far,tol);inc=out.conjugate()
            # Exact Wronskian identity avoids subtracting two nearly equal S matrices.
            response=-R*W(a,b)/(W(out,a)*W(b,inc))
            responses.append(dict(tol=tol,outer_phase=far,period_days=days,kR=k,
                real_m=float(response.real),imag_m=float(response.imag),static_m=chi,
                phase_rad=float(np.angle(response)),phase_delay_s=float(np.angle(response)/omega)))
    # Conditional analytic envelope uses rational bounds and pi < 22/7.
    import fractions
    F=fractions.Fraction
    save('scalar-envelope-pilot-failure.json',dict(classification='Counterexample candidate',
        original_density_upper_kg_m3=1e8,actual_max=float(max(data['energy']/gr.C**2)),passed=False,
        reason='The first proposed value envelope did not contain the central density. A separately stated analytic envelope uses 2e8; no original containment pass is asserted.'))
    kappa=2*F(22,7)*4*F('6.67428e-11')*F('2e8')*F('7e7')**2/(F(299792458)**2*F('0.99')**2)
    assert kappa<1
    enclosure=dict(R_upper_m=7e7,energy_over_c2_upper_kg_m3=2e8,N_lower=.99,N_upper=1.,f_lower=.99,
                   beta_abs_upper=4.,kappa_upper=float(kappa),kappa_exact=str(kappa))
    assert R<7e7 and max(data['energy']/gr.C**2)<2e8 and min(np.exp(data['nu']))>.99
    assert min(1-2*np.divide(data['mass'],data['r'],out=np.zeros_like(data['r']),where=data['r']>0))>.99
    errors=[abs(responses[i+2]['real_m']/responses[i]['real_m']-1) for i in range(2)]
    save('scalar-response.json',dict(classification='Counterexample candidate',R_m=R,M_geom_m=M,
        static_susceptibility_m=statics,static_relative_change=abs(statics[1]/statics[0]-1),
        response=responses,real_response_refinement_relative=errors,
        passed=max(errors)<1e-4,
        convention='l=0 scattering response R*(S_beta/S_0-1)/(2 i kR) on the declared GR-background interpolant, e^(-i omega t). S_0 removes the beta=0 curvature scattering. Not a binary force/TOA transfer function.',
        limits='Finite-profile and exterior-asymptotic numerical controls only. No full fluid/thermal stability or finite-background nonlinear inference.'))
    save('scalar-stability-bound.json',dict(classification='Proven',conditional_envelope=enclosure,
        theorem='For regular asymptotically vanishing perturbations, integral V phi^2 <= kappa integral p (phi prime)^2. kappa<1 excludes a negative omega^2 scalar bound state; angular gradients are nonnegative.',
        geometric_conditions='Static spherical metric N in [0.99,1], f>=0.99, compact support radius<=7e7 m, 0<=e-3P<=e<=2e8*c^2 SI, |beta|<=4. Schwarzschild exterior included.',
        membership='The declared shape-preserving finite-profile interpolation satisfies these value envelopes; no interval enclosure of the original EOS/TOV continuum is claimed.',
        exclusions='Not fluid, convection or thermal stability. No theorem that every stable response has a single relaxation pole.'))
    print('SCALAR',statics,'errors',errors,'conditional kappa',float(kappa),flush=True)

def symbolic():
    import sympy as sp
    x,b,eps,z,lam,ds,B,C=sp.symbols('x b eps z lam ds B C')
    # Finite values and first derivatives do not determine an interval derivative bound.
    hidden=sp.prod((x-i)**2 for i in [-2,-1,0,1,2])
    assert all(hidden.subs(x,i)==0 and sp.diff(hidden,x).subs(x,i)==0 for i in [-2,-1,0,1,2])
    assert sp.diff(hidden,x).subs(x,sp.Rational(1,3))!=0
    # At the even-coupling GR branch, the matter source has no linear thermal cross term.
    phi0,delta,t0,thermal=sp.symbols('phi0 delta t0 thermal')
    source=b*(phi0+eps*delta)*(t0+eps*thermal)
    assert sp.diff(source,eps).subs({eps:0,phi0:0})==b*delta*t0
    matrix=sp.Matrix([[ds,-eps*B],[-eps*C,z+lam]])
    transfer=sp.factor(matrix.inv()[0,0])
    assert sp.simplify(transfer-(z+lam)/(ds*(z+lam)-eps**2*B*C))==0
    assert sp.simplify(transfer.subs(eps,0)-1/ds)==0
    # Stable Wronskian scattering identity avoids an O(k*chi/R) subtraction.
    v,w,out,inc=[sp.Matrix(sp.symbols(prefix+'0 '+prefix+'1')) for prefix in ['v','w','o','i']]
    def W(a,b): return a[0]*b[1]-a[1]*b[0]
    assert sp.expand(W(v,inc)*W(out,w)-W(w,inc)*W(out,v)-W(v,w)*W(out,inc))==0
    # Fourier-law finite-volume energy and quadratic dissipation, arbitrary positive conductances.
    t=sp.symbols('t0:4');g=sp.symbols('g0:3',positive=True)
    flux=[-g[i]*(t[i+1]-t[i]) for i in range(3)]
    rates=[-flux[0],flux[0]-flux[1],flux[1]-flux[2],flux[2]]
    assert sum(rates)==0
    assert sp.expand(sum(t[i]*rates[i] for i in range(4))+sum(g[i]*(t[i+1]-t[i])**2 for i in range(3)))==0
    # Analytic uniform sphere independently checks scalar source sign and matching.
    k=.04;start=1e-6
    def rhs(r,y): return [y[1]/r**2,-k*k*r*r*y[0]]
    sol=solve_ivp(rhs,(start,1),[1-k*k*start**2/6,-k*k*start**3/3],method='DOP853',rtol=1e-12,atol=1e-14,max_step=.01)
    a,j=sol.y[:,-1];got=-j/(a+j);expected=np.tan(k)/k-1
    assert abs(got/expected-1)<1e-10
    cert=json.loads((OUT/'scalar-stability-bound.json').read_text())['conditional_envelope']
    kap=sp.Rational(cert['kappa_exact']);a=kap/4;R=sp.Integer(70000000)
    derivative_bounds=[sp.factor(R*sp.factorial(n)*a**n/(1-kap)**(n+1)) for n in [1,2]]
    save('symbolic.json',dict(classification='Proven',
        finite_sampling_counterexample=str(hidden),even_GR_branch_thermal_decoupling=True,
        coupled_transfer=str(transfer),Wronskian_identity=True,conservative_diffusion_dissipation=True,
        uniform_sphere_relative_error=float(abs(got/expected-1)),
        static_beta_derivative_bounds_m=[str(v) for v in derivative_bounds],
        static_beta_derivative_bounds_m_decimal=[float(v) for v in derivative_bounds],
        derivative_domain='b=abs(beta) in [0,4], fixed metric/matter satisfying the separately declared scalar value envelope; not a derivative certificate for EOS-dependent stellar reconstruction or timing.'))
    print('PASS symbolic closure, decoupling, Wronskian, finite-sampling no-go and uniform sphere',flush=True)

def source_bindings():
    names=['star/private/micro.f90','star/private/report.f90','eos/private/eosdt_eval.f90',
           'net/private/net_screen.f90','net/private/net_initialize.f90','rates/private/eval_weak.f90',
           'data/rates_data/weakreactions.tables','const/public/const_def.f90']
    for rel in names:
        path=OUT/'sources'/rel;path.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(fresh.MESA/rel,path)
    save('source-bindings.json',dict(classification='Imported from prior work',
        sha256={p.relative_to(OUT/'sources').as_posix():sha(p) for p in sorted((OUT/'sources').rglob('*')) if p.is_file()},
        literature=[dict(title='Damour and Esposito-Farese (1996), tensor-scalar field equations',url='https://arxiv.org/abs/gr-qc/9602056'),
                    dict(title='Potekhin et al. (2020), redshifted thermal evolution',url='https://doi.org/10.1093/mnras/staa1871')],
        note='External papers supply framework context; all new numerical claims are derived from this repository. Sources and old runtime build identity remain distinct.'))

def ledger():
    rows=[
        ('native_temperature_derivative','conditional_pass','PP branching derivative omission matches the frozen-state error; correction residual <3e-8. The original binary and failure are preserved.',[]),
        ('reaction_channels','conditional_pass','63 whole-network grouped interventions reproduce heat/nu/native derivatives; isolated-channel and individual-PP failures retained.',[]),
        ('composition_and_energy','partial','22-component vector reconstructed from channel outputs/source Q. Direct native vector and independent time-step comparison absent; F-free manifold not invariant.',['native_species_export','EOS_including_generated_F']),
        ('common_EOS_and_global_derivatives','blocked','Default small-step EOS differences do not converge; HELM full-domain control also fails. Local numerical agreement cannot certify global errors.',['verified_common_EOS_evaluator','controlled_numerical_error_and_regularity_bounds']),
        ('GR_heat_and_composition_evolution','blocked','Conservative frozen-coefficient controls pass energy tests. Closed boundary exits the local regime; radiative boundary still changes local T by about 15%. No species/metric/density/opacity updates.',['common_EOS_and_global_derivatives','composition_and_energy','moving_structure_convection_atmosphere']),
        ('GR_scalar_stability_response','conditional_pass','Curved-background monopole scattering and an explicit conditional scalar stability bound; no full stellar mode or binary transfer claim.',[]),
        ('fluid_thermal_scalar_coupling','blocked','Buoyancy diagnostic has agreeing negative regions. At phi0=0 thermal modes decouple linearly; finite-background Schur residue needs the missing cross derivatives.',['GR_heat_and_composition_evolution','fluid_mode_boundary_conditions','nonzero_background_couplings']),
        ('complete_nonlinear_observation','blocked','No certified time-dependent three-body force and photon/noise/pulse forward model follows from a static mass or a monopole scattering coefficient. Existing raw observational failures remain frozen.',['fluid_thermal_scalar_coupling','full_initialization_flow_readout_error_bounds']),
        ('final_submission','pending_integration','Request12 PDF/ZIP are historical. A conditional methodology paper needs manuscript integration and claim-scoped review; full observational closure is required only for the stronger physical claims.',['manuscript_integration','claim_scoped_submission_review'])]
    save('lever-ledger.json',dict(classification='Proven',
        scope='Every currently named lever has a concrete evaluation, theorem boundary, or unresolved dependency. Exhausting accessible tests is not proving every scientific objective.',
        levers=[dict(name=n,status=s,evidence=e,missing=m) for n,s,e,m in rows],
        permission_needed_for_completed_work=False,
        repo_boundary='AGENTS.md excludes runtime/build-environment and additional pulsar/LLR empirical work. No such work was reopened. A mathematical missing specification is not solved by requesting the same approval again.'))

def maintain():
    docs=['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi']
    prior=json.loads((OLD/'manifest.json').read_text())['sha256'];bindings={};dest=OUT/'request25-notes';dest.mkdir()
    for rel in ['docs/'+k+'.md' for k in docs]+['paper/revision-manifest.json']:
        assert sha(ROOT/rel)==prior[rel]
        snap=dest/Path(rel).name;shutil.copy2(ROOT/rel,snap)
        bindings[rel]=dict(snapshot=snap.relative_to(ROOT).as_posix(),sha256=prior[rel],historical_manifest=rel.startswith('paper/'))
    save('historical-note-bindings.json',bindings)
    bodies=[
        '분류: Counterexample candidate. Request26은 실제 네트워크를 63개 반응 그룹으로 조정하여 가열·중성미자·미분을 재구성하고, 소스 Q 정의로 22종 조성 변화 방향과 에너지 장부를 복원했다. 이는 native 전체 변화율의 독립 출력 검증과 구분한다. 초기 불소가 없어도 양의 생성 방향이 있어 19종 불소 생략 모형의 정확한 시간 불변성을 주장할 수 없다.',
        '분류: Counterexample candidate. 물질 표본·수학적 대기를 연결한 GR 계량/trace 보간 배경에서 scalar 정적 감수율 약 1166.852694 m를 얻었다. 두 지정 일 단위 주파수의 곡률 제거 monopole 산란 위상/ω는 약 3.8922 마이크로초다. 다체 힘·광자·TOA 전달함수나 단일 완화 pole의 증명은 아니다.',
        '분류: Proven. 명시한 계량·trace 값 영역에서 scalar 에너지의 음의 퍼텐셜 항은 양의 구배 항의 0.018670 미만으로 상계되어 음의 주파수 제곱 모드를 배제한다. 같은 고정 배경의 정적 beta 도함수 상계도 얻었다. 전체 EOS/TOV 연속 해의 영역 소속이나 유체·열 안정성 인증은 아니다. 영 scalar 가지에서는 열 교란이 선형 scalar 방정식에서 분리된다.',
        '분류: Counterexample candidate. PP 분기비 미분 누락을 실제 채널 출력으로 검산하여 지정 상태의 기존 온도 미분 오차를 설명했다. 보정 후 차분 잔차는 최대 약 2.74e-8이다. 비선형 복사 경계를 가진 GR 동결 계수 열수송 대조는 에너지 수지를 통과했지만 국소 온도 변화가 약 14.86%여서 물리적 시간 진화 인증으로 세지 않는다.\n\n분류: Proven. 열 변수를 제거한 scalar 응답에는 교차 결합의 곱이 들어간다. 긴 열 시간만으로 관측 가능한 느린 scalar pole을 주장할 수 없다.',
        '분류: Counterexample candidate. 단독 반응 합산과 개별 PP 비율 조정의 실패를 보존하고, 공유 PP 분기를 묶은 전체 네트워크 대조를 별도로 통과했다. EOS의 작은 간격 차분은 수렴하지 않았고 HELM 전체 영역 대조도 실패했다. 이를 연속 EOS의 제1법칙 위반으로 확정하지 않으며 전역 미분 보증으로도 사용하지 않는다. 두 미분법의 국소 부력 진단에서 음의 영역을 확인했지만 고유모드 인증은 아니다.\n\n분류: Conjectural. 공통 EOS·독립 native 조성 벡터·열/유체/metric 진화·비영 배경의 교차 결합과 완전한 관측 전방 모형이 남는다. 모든 현재 레버를 실행 결과 또는 정확한 누락 조건으로 원장에 기록했으며, 접근 가능한 검사의 소진을 연구 전체의 완료로 바꾸지 않는다.'
    ]
    revision=json.loads((ROOT/'paper/revision-manifest.json').read_text())
    for stem,body in zip(docs,bodies):
        rel='docs/'+stem+'.md'
        with (ROOT/rel).open('ab') as f:
            f.write(('\n\n## Request 26 남은 폐쇄 조건의 실행 검산\n\n'+body+'\n\n세부 근거: [한글 보고서](../notes/REQUEST26_REMAINING_CLOSURE_KO.md).\n').encode())
        revision['sha256'][rel]=sha(ROOT/rel)
    revision['request26_supporting_note_update']=dict(evidence_manifest='outputs/remaining-closure26/manifest.json',
        historical_notes='outputs/remaining-closure26/historical-note-bindings.json',
        status='PP 미분 원인·반응 채널·조성 방향·GR scalar 응답 및 보존형 열 대조 검산; 물리적 미완료 조건 명시',
        artifact_status='원고 PDF와 ZIP은 Request12 역사 산출물이다.')
    (ROOT/'paper/revision-manifest.json').write_text(json.dumps(revision,ensure_ascii=False,indent=2)+'\n')

def seal():
    imports=list(OUT.glob('*-import-control.json'));assert len(imports)==146
    assert all(json.loads(p.read_text())['passed'] for p in imports)
    save('gates.json',dict(classification='Proven',valid_initial_models=146,exported_zone_states=146*5735,
        failed_initial_configurations=3,grouped_channel_reconstruction_passed=True,
        PP_local_derivative_correction_passed=True,reconstructed_composition_energy_checked=True,
        conditional_GR_scalar_bound_and_response_checked=True,conservative_thermal_solver_controls_passed=True,
        native_dxdt_independent_validation=False,common_EOS_certified=False,global_derivative_certificate=False,
        complete_GR_thermal_composition_evolution=False,full_fluid_thermal_mode_stability=False,
        nonlinear_binary_force_and_readout_certified=False,complete_nonlinear_observational_inference=False,
        final_submission_package_updated=False,
        classification_detail='Theorem progress: conditional GR scalar energy/derivative bounds, thermal decoupling, conservative diffusion and finite-sampling boundary. Loophole progress: fresh channel/PP derivative evidence, reconstructed composition direction and curved-background response. Physical closure remains incomplete.'))
    save('provenance.json',dict(classification='Proven',checkpoint='36f655d',
        previous_manifest_sha256=sha(OLD/'manifest.json'),
        runtime_binding='outputs/fresh-microphysics25/runtime.json',data_binding='outputs/fresh-microphysics25/runtime-data-bindings.json',
        policies='Keep all original failures and later diagnostic domains separate. No runtime rebuild or new empirical pulsar/LLR work. No approval is missing for the completed tests.'))
    files=[p for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json']
    files += [ROOT/'verification/remaining_closure.py',ROOT/'notes/REQUEST26_REMAINING_CLOSURE_KO.md',ROOT/'.gitattributes']
    files += [ROOT/k for k in json.loads((OUT/'historical-note-bindings.json').read_text())]
    save('manifest.json',dict(classification='Proven',sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(files)}))

def verify():
    stages=['validated-variational','remaining-levers15','nbody-readout16','nonzero-drive17','thermal-wd18',
        'thermal-restart19','thermal-robustness20','gr-mass21','thermal-closure22','baryon-entropy23','reactive-energy24',
        'fresh-microphysics25','remaining-closure26']
    histories={f'outputs/{a}/manifest.json':f'outputs/{b}/historical-note-bindings.json' for a,b in zip(stages,stages[1:])};count=0
    for label in ['outputs/research-remediation/manifest.json',*histories,'paper/revision-manifest.json','outputs/remaining-closure26/manifest.json']:
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
    assert sha(fresh.BINARY)==json.loads((OLD/'runtime.json').read_text())['sha256']
    data=json.loads((OLD/'runtime-data-bindings.json').read_text())
    for path,digest in data['sha256'].items(): assert sha(Path(data['root'])/path)==digest,path
    for key in ['runtime_sha256','module_sha256']:
        for path,digest in json.loads((ROOT/'outputs/gr-mass21/provenance.json').read_text())[key].items(): assert sha(path)==digest,path
    assert not json.loads((OUT/'channel-sum.json').read_text())['passed']
    assert json.loads((OUT/'fullnet-channel-sum.json').read_text())['passed']
    assert not json.loads((OUT/'EOS-controls.json').read_text())['all_passed']
    assert not json.loads((OUT/'HELM-control.json').read_text())['passed']
    gates=json.loads((OUT/'gates.json').read_text())
    for key in ['native_dxdt_independent_validation','common_EOS_certified','global_derivative_certificate',
        'complete_GR_thermal_composition_evolution','full_fluid_thermal_mode_stability',
        'nonlinear_binary_force_and_readout_certified','complete_nonlinear_observational_inference','final_submission_package_updated']:
        assert gates[key] is False,key
    plans=[('fullnet-channel-plan.json','frozen_isolated_failure_sha256','channel-sum.json'),
           ('EOS-refinement-plan.json','frozen_first_failure_sha256','EOS-controls.json')]
    for plan,key,target in plans: assert json.loads((OUT/plan).read_text())[key]==sha(OUT/target)
    print('PASS',count,'artifact/history SHA;',len(data['sha256']),'MESA data SHA; scientific failure gates',flush=True)

def recheck():
    import contextlib,io
    global OUT,CACHE
    original,cache=OUT,CACHE;oldfresh=(fresh.OUT,fresh.CACHE);worst=0.;valid=failed=0
    expected_failures={'off-r34_pp2','half-r34_pp2','channel-rne20ap_to_mg24'}
    with tempfile.TemporaryDirectory(prefix='remaining-closure26-') as folder:
        OUT=Path(folder)/'out';CACHE=Path(folder)/'runs';CACHE.mkdir();shutil.copytree(original,OUT);context()
        try:
            for inputs in sorted(original.glob('*-inputs.json')):
                label=inputs.name[:-len('-inputs.json')];run_dir=CACHE/label;run_dir.mkdir();shutil.copy2(fresh.BINARY,run_dir/'binary')
                for p in (original/'inputs'/label).iterdir():
                    if p.name=='input.mod.gz':
                        with fresh.gzip.open(p,'rb') as src,(run_dir/'input.mod').open('wb') as dst: shutil.copyfileobj(src,dst)
                    else: shutil.copy2(p,run_dir/p.name)
                try:
                    with contextlib.redirect_stdout(io.StringIO()): fresh.run(label);fresh.collect(label)
                except AssertionError as err:
                    if label not in expected_failures: raise
                    assert isinstance(err.args[0],tuple) and err.args[0][0]=='No complete starting state',err
                    marker='rate_for_pg_pa_branches' if label.startswith('channel-') else 'set_combo_screen_rates'
                    assert marker in (run_dir/'execution.log').read_text(),label
                    failed+=1;continue
                assert label not in expected_failures,label
                _,a=mesa(OUT/(label+'-profile.data.gz'));_,b=mesa(original/(label+'-profile.data.gz'))
                for key in ['rho','logT',*fresh.ISOS,'eps_nuc','eps_nuc_neu_total','non_nuc_neu','opacity','d_lnepsnuc_dlnT','d_lnepsnuc_dlnd']:
                    error=float(max(abs(a[key]-b[key])/np.maximum(1e-30,abs(b[key]))));worst=max(worst,error)
                    assert error<1e-8,(label,key,error)
                valid+=1
                if valid%10==0: print('RECHECK',valid,'valid profiles;',failed,'expected failed configurations',flush=True)
            assert valid==146 and failed==3,(valid,failed)
            with contextlib.redirect_stdout(io.StringIO()):
                channel_analysis();fullnet_analysis();EOS_analysis();EOS_refinement_analysis();helm_analysis()
                pp_analysis();composition_reconstruction();scalar();transport();radiative_transport();fluid_diagnostic();symbolic()
            assert json.loads((OUT/'fullnet-channel-sum.json').read_text())['passed']
            assert all(r['passed'] for r in json.loads((OUT/'pp-derivative-result.json').read_text())['results'])
            assert not json.loads((OUT/'HELM-control.json').read_text())['passed']
            assert json.loads((OUT/'radiative-GR-transport.json').read_text())['energy_gate_passed']
            assert json.loads((OUT/'scalar-response.json').read_text())['passed']
        finally: OUT,CACHE=original,cache;fresh.OUT,fresh.CACHE=oldfresh
    print('PASS 146 actual profiles, 3 original failed configurations, all numerical/symbolic checks; max repeat difference',worst,flush=True)

if __name__=='__main__': globals()[sys.argv[1]]()
