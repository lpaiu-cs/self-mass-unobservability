"""Request25: evaluate existing MESA microphysics at supplied thermodynamic states.

Counterexample candidate: imported state evaluation, not a GR evolution run.
"""
from pathlib import Path
import json, os, shutil, subprocess, sys, time
import gzip, tempfile
import numpy as np
from thermal_wd import mesa
from thermal_restart import sha
from thermal_robustness import gzcopy
import gr_mass as gr
import baryon_entropy as be

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/fresh-microphysics25'
CACHE=Path('/home/lpaiu/work/fresh-microphysics25')
MESA=Path('/home/lpaiu/work/thermal-restart19/mesa-r7624')
OLD=ROOT/'outputs/reactive-energy24'
SOURCE=ROOT/'outputs/thermal-closure22/selected.data.gz'
BINARY=Path('/home/lpaiu/work/thermal-closure22/replay/binary')
LSUN=3.8418e33
ISOS=json.loads((ROOT/'outputs/thermal-restart19/network-diagnosis.json').read_text())['profile_isotopes']


def save(name,obj):
    path=OUT/name;path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists() and not CACHE.exists()
    OUT.mkdir(parents=True);CACHE.mkdir()
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='367bc1e',
        previous_manifest_sha256=sha(OLD/'manifest.json'),
        method='Reuse the unchanged MESA binary and its saved-model loader. finish_load_model calls set_vars with net/neu/kap enabled; profile_starting_model exports the pre-evolution state. Require the exported rho,T,X to match the supplied state before using powers or opacity.',
        stages=['Original source-state import control','Fresh evaluation at the scaled Request23 GR material state',
                'Finite temperature/density perturbations versus reported nuclear derivatives',
                'Redshifted power and frozen-rate comparison on the same material template'],
        tolerances=dict(import_logrho=1e-10,import_logT=1e-10,import_X_absolute=1e-12,
            source_power_relative=1e-6,nuclear_derivative_relative=1e-3),
        gates='Record any failed source or state control without replacing it with a later fit. Imported models are a carrier for state evaluation, not Newtonian or GR evolution solutions.',
        limitations='MESA supplies its own EOS electron/ion auxiliary quantities at the requested baryon rho,T,X. This does not certify equality to the FreeEOS auxiliaries, expose every runtime reaction flow, close the composition derivative, or solve transport and GR time evolution.'))
    for name in ['star/private/read_model.f90','star/private/hydro_vars.f90','star/private/net.f90',
        'star/job/run_star_support.f90','star/private/profile_getval.f90','star/private/write_model.f90',
        'star/defaults/star_job.defaults','star/defaults/controls.defaults']:
        path=OUT/'sources'/name;path.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(MESA/name,path)
    save('runtime.json',dict(classification='Imported from prior work',binary=str(BINARY),sha256=sha(BINARY),
        original_binding='outputs/thermal-closure22/inputs.json',toolchain_changed=False))


def write_model(folder,data,mass_g,age):
    fields=['lnd','lnT','lnR','L','dq']+ISOS
    rows=np.c_[data['lnd'],data['lnT'],data['lnR'],data['L'],data['dm']/mass_g,data['X']]
    assert np.max(abs(data['X'].sum(axis=1)-1))<1e-12
    text='! Request25 supplied state; initialization only.\n!\n0 -- model for mesa/star\n\n'
    for key,value in [('M/Msun',mass_g/1.9892e33),('model_number',1),('star_age',age),
        ('initial_z',.02),('n_shells',len(rows)),('net_name',"'cno_extras.net'"),('species',22)]:
        text+=f'{key:>32}  {value if isinstance(value,(str,int)) else format(value,".17e")}\n'
    text+='\n'+' '.join(fields)+'\n'
    with (folder/'input.mod').open('w') as f:
        f.write(text)
        for i,row in enumerate(rows): f.write(str(i+1)+' '+' '.join(f'{v:.17e}' for v in row)+'\n')


def setup_run(label,data,age,screening=None,weak_factor=None):
    folder=CACHE/label;folder.mkdir(exist_ok=False)
    shutil.copy2(BINARY,folder/'binary')
    shutil.copy2(ROOT/'outputs/thermal-closure22/inputs/cno_extras.net',folder/'cno_extras.net')
    shutil.copy2(ROOT/'outputs/thermal-closure22/inputs/profile_columns.list',folder/'profile_columns.list')
    with (folder/'profile_columns.list').open('a') as f: f.write('\neta\ntheta_e\nabar\nzbar\nz2bar\nye\n')
    # Keep the source microphysics and atmosphere controls, disable initial
    # physical changes; the starting profile is selected before any evolution.
    text=(ROOT/'outputs/thermal-closure22/inputs/inlist1').read_text()
    start=text.index('&controls');text=text[start:]
    text=text.replace('max_model_number = 19042','max_model_number = 1')
    if screening is not None: text=text.replace('&controls',"&controls\n screening_mode = '"+screening+"'\n",1)
    if weak_factor is not None: text=text.replace('&controls',f'&controls\n weak_rate_factor = {weak_factor}\n',1)
    job="""&star_job
 load_saved_model = .true.
 saved_model_name = 'input.mod'
 profile_starting_model = .true.
 show_log_description_at_start = .false.
 pgstar_flag = .false.
 save_model_when_terminate = .true.
 save_model_filename = 'final.mod'
 write_profile_when_terminate = .true.
 filename_for_profile_when_terminate = 'final_profile.data'
/
"""
    (folder/'inlist1').write_text(job+text)
    (folder/'inlist').write_text("""&binary_job
 inlist_names(1) = 'inlist1'
 inlist_names(2) = 'inlist2'
 evolve_both_stars = .false.
/
&binary_controls
 m1 = 0.1975
 m2 = 1.4
 initial_period_in_days = 3.4
 do_tidal_sync = .false.
 do_initial_orbit_sync_1 = .false.
 do_jdot_mb = .false.
/
&binary_pgstar
/
""")
    (folder/'inlist2').write_text('&star_job\n/\n&controls\n/\n&pgstar\n/\n')
    write_model(folder,data,float(sum(data['dm'])),age)
    np.savez_compressed(OUT/(label+'-input.npz'),**data)
    dest=OUT/'inputs'/label;dest.mkdir(parents=True)
    for path in folder.iterdir():
        if path.name=='binary': continue
        if path.name=='input.mod': gzcopy(path,dest/'input.mod.gz')
        else: shutil.copy2(path,dest/path.name)
    save(label+'-inputs.json',dict(classification='Proven',sha256={p.name:sha(p) for p in folder.iterdir() if p.is_file()}))


def source_input():
    h,d=mesa(SOURCE)
    data=dict(lnd=np.log(d['rho']),lnT=d['logT']*np.log(10),lnR=np.log(d['radius']*gr.RSUN*100),
        L=d['luminosity']*LSUN,dm=d['dm'],X=np.array([d[k] for k in ISOS]).T)
    setup_run('source',data,float(h['star_age']))


def run(label):
    folder=CACHE/label
    for name,digest in json.loads((OUT/(label+'-inputs.json')).read_text())['sha256'].items(): assert sha(folder/name)==digest,name
    env=os.environ.copy();runtime=MESA.parent
    env.update(MESA_DIR=str(MESA),LD_LIBRARY_PATH=str(runtime/'mesasdk/lib')+':'+str(runtime/'mesasdk/lib64'),OMP_NUM_THREADS='4')
    assert not (folder/'execution.log').exists()
    begin=time.monotonic()
    with np.load(OUT/(label+'-input.npz')) as data: expected_rows=len(data['dm'])
    with (folder/'execution.log').open('w') as f:
        process=subprocess.Popen(['./binary'],cwd=folder,env=env,stdout=f,stderr=subprocess.STDOUT)
        complete=False
        try:
            while time.monotonic()-begin<180:
                path=folder/'LOGS1/profile1.data'
                if path.exists():
                    try:
                        h,d=mesa(path)
                        complete=int(h['model_number'])==1 and len(d['zone'])==expected_rows
                        complete=complete and d['zone'][-1]==len(d['zone']) and np.all(np.isfinite(d['eps_nuc']))
                    except (ValueError,IndexError,KeyError,OSError): complete=False
                    if complete: break
                if process.poll() is not None: break
                time.sleep(.1)
        finally:
            if process.poll() is None: process.terminate()
            process.wait(timeout=10)
        assert complete,('No complete starting state',label,process.returncode)
    shutil.copy2(folder/'execution.log',OUT/(label+'-execution.log'))
    save(label+'-execution.json',dict(classification='Proven',returncode=process.returncode,
        elapsed_s=time.monotonic()-begin,stopped_after_complete_starting_profile=True,
        selection='Only model 1 is used, before the first evolution step. Terminate the disposable process after the complete profile appears.'))
    print(label,'starting profile captured; returncode',process.returncode,flush=True)


def run_source(): run('source')


def collect(label):
    folder=CACHE/label;h,d=mesa(folder/'LOGS1/profile1.data');inp=np.load(OUT/(label+'-input.npz'))
    assert int(h['model_number'])==1 and len(d['zone'])==len(inp['dm'])
    actual=np.array([d[k] for k in ISOS]).T
    errors=dict(logrho=float(np.max(abs(np.log(d['rho'])-inp['lnd']))),
        logT=float(np.max(abs(d['logT']*np.log(10)-inp['lnT']))),
        X_absolute=float(np.max(abs(actual-inp['X']))),
        dm_relative=float(np.max(abs(d['dm']/inp['dm']-1))))
    passed=errors['logrho']<1e-10 and errors['logT']<1e-10 and errors['X_absolute']<1e-12 and errors['dm_relative']<1e-10
    if label=='source':
        _,old=mesa(SOURCE)
        controls={k:float(np.max(abs(d[k]-old[k])/np.maximum(abs(old[k]),1e-30))) for k in
            ['eps_nuc','eps_nuc_neu_total','non_nuc_neu','opacity']}
        passed=passed and all(v<1e-6 for v in controls.values())
        errors['source_microphysics_max_relative']=controls
    save(label+'-import-control.json',dict(classification='Counterexample candidate',passed=bool(passed),errors=errors,
        evolved_profile_used=False,version_note='Reused data/version_number contains 7623; source archive is named r7624. This metadata is not an independent compiled-build certificate. Exact executable SHA and numerical import controls define the runtime boundary.'))
    gzcopy(folder/'LOGS1/profile1.data',OUT/(label+'-profile.data.gz'))
    assert passed,(label,errors)
    print(label,'PASS import',errors,flush=True)


def collect_source(): collect('source')


def gr_states(emit_run=True):
    solver=be.Solver(33,4);m=solver.mat;n=len(m.lp)
    row=json.loads((ROOT/'outputs/baryon-entropy23/scaled-33-4.json').read_text())
    scale=row['baryon_scale'];B=m.B*scale
    error,inner,outer=solver.branches([*row['parameters'][:2],0.],scale,record=True)
    assert max(abs(error))<1e-8
    faces=np.zeros((n+1,3));faces[0]=[row['photospheric_radius_m']/m.R,gr.TARGET/B,m.lp[0]]
    faces[1:m.split+1]=outer[1:,1:];faces[m.split:n]=inner[1:,1:][::-1]
    faces[n]=[0.,0.,row['parameters'][0]]
    states=[];lp=[];lt=[];rho=[];energy=[];entropy=[];nu=[]
    atm=be.atmosphere(row);vacM=atm['ADM_mass_GM_solar']*gr.GM_SUN/gr.C**2
    nvac=.5*np.log1p(-2*vacM/atm['vacuum_radius_m'])
    surface=solver.ref[0];pv=surface[1]/surface[0];u0=surface[2]-1.5*pv
    nuface=nvac-np.log1p(2.5*pv/((gr.C*100)**2+u0))
    nuphoto=nuface
    for i in range(n):
        isouter=i<m.split
        if isouter:
            begin=outer[i,0];middle=(m.outer[i]+m.outer[i+1])/2;start=outer[i,1:]
        else:
            j=n-1-i;begin=inner[j,0];middle=(m.inner[i]+m.inner[i+1])/2;start=inner[j,1:]
        y=solver.step(np.log(begin),np.log(middle),start,i,B,isouter)
        vals=[]
        for p in [faces[i,2],y[2],faces[i+1,2]]:
            a,t,_=be.invert(solver.eos,p,solver.ref[i,3],m.eps[i],m.lt[i]);vals.append((a,t))
        a,t=vals[1];ab,ae=vals[0][0],vals[2][0]
        # At fixed X,s per cell, dnu=-d ln h. Evaluate thermal differences
        # before division by c^2, avoiding subtraction of rest energies.
        ht=ab[2]+ab[1]/ab[0];hm=a[2]+a[1]/a[0];hb=ae[2]+ae[1]/ae[0]
        nu.append(nuface-np.log1p((hm-ht)/((gr.C*100)**2+ht)))
        nuface-=np.log1p((hb-ht)/((gr.C*100)**2+ht))
        states.append(y);lp.append(y[2]);lt.append(t);rho.append(a[0]/m.cx[i]);energy.append(a[2]*m.cx[i]);entropy.append(a[3]*m.cx[i])
    states=np.array(states)
    X=np.array([m.d[k] if not k.startswith('f') else np.zeros(n) for k in ISOS]).T
    X/=X.sum(axis=1)[:,None]
    data=dict(lnd=np.log(rho),lnT=np.array(lt),lnR=np.log(faces[:-1,0]*m.R*100),
        L=m.d['luminosity']*LSUN,dm=m.d['dm']*scale,X=X,
        logP=np.array(lp),r_mid_m=states[:,0]*m.R,m_mid_geom=states[:,1]*B,
        nu=np.array(nu),u_W=np.array(energy),s_B=np.array(entropy),CX=m.cx)
    np.savez_compressed(OUT/'gr-material-state.npz',**data)
    save('gr-state-control.json',dict(classification='Counterexample candidate',zones=n,
        midpoint='Integrate the frozen cell TOV branch to each baryon midpoint; direct FreeEOS entropy inversion there and at both faces.',
        interface_max=float(max(abs(error))),surface_nu=float(nuphoto),central_nu=float(nuface),
        entropy_max_relative=float(np.max(abs(data['s_B']/(solver.ref[:,3]*m.cx)-1))),
        normalization='nu matched through the same mathematical isentropic atmosphere; within each material cell use dnu=-d ln h. No heat transport in atmosphere.',
        EOS='19-isotope non-F composition, grouped by element in FreeEOS. The 22-isotope reaction network receives zero initial F fractions.',
        closure='Numerical material-state sampling, no new mass fit or time evolution; original 5735-cell template retained.'))
    if emit_run: setup_run('gr',data,float(m.h['star_age']))
    print('GR material states prepared',n,flush=True)


def run_gr(): run('gr');collect('gr')


def perturb_inputs():
    h,d=mesa(SOURCE);base=dict(np.load(OUT/'gr-input.npz'))
    _,source=mesa(SOURCE)
    x=np.array([source[k] if not k.startswith('f') else np.zeros(len(source['dm'])) for k in ISOS]).T
    x/=x.sum(axis=1)[:,None]
    source19=dict(np.load(OUT/'source-input.npz'));source19['X']=x
    setup_run('source19',source19,float(h['star_age']))
    labels=['source19']
    for field,prefix in [('lnT','T'),('lnd','rho')]:
        for step in [1e-3,5e-4]:
            for sign,suffix in [(-1,'minus'),(1,'plus')]:
                label=prefix+str(step)+'-'+suffix;data={k:v.copy() for k,v in base.items()}
                data[field]+=sign*step
                setup_run(label,data,float(h['star_age']));labels.append(label)
    save('perturbation-plan.json',dict(classification='Counterexample candidate',labels=labels,
        derivative_domain='abs(eps_nuc)>1 erg/g/s; report covered integrated absolute nuclear power and max derivative discrepancy with denominator max(1,abs(reported normalized derivative)).',
        initial_steps=[1e-3,5e-4],tolerance=1e-3))


def run_perturbations():
    for label in json.loads((OUT/'perturbation-plan.json').read_text())['labels']:
        run(label);collect(label)


def analyze():
    state=np.load(OUT/'gr-material-state.npz');h,g=mesa(OUT/'gr-profile.data.gz')
    _,source=mesa(OUT/'source-profile.data.gz');_,s19=mesa(OUT/'source19-profile.data.gz')
    dm=state['dm'];weight=dm*np.exp(2*state['nu'])/LSUN
    powers={}
    for label,d in [('frozen_source22',source),('fresh_source19',s19),('fresh_GR',g)]:
        values={k:float(np.dot(weight,d[k])) for k in ['eps_nuc','eps_nuc_neu_total','non_nuc_neu']}
        values['nuclear_minus_thermal_neutrino']=values['eps_nuc']-values['non_nuc_neu']
        powers[label]=values
    boundary=json.loads((OUT/'gr-state-control.json').read_text());oh,_=mesa(SOURCE)
    Linf=float(oh['photosphere_L'])*np.exp(2*boundary['surface_nu'])
    peak=int(np.argmax(g['eps_nuc']));dense=(np.exp(state['lnT'])>1e6)&(g['eps_nuc']>1.)
    ratio=g['opacity']/source['opacity']
    result=dict(classification='Counterexample candidate',power_infinity_Lsun=powers,
        retained_photospheric_L_infinity_Lsun=float(Linf),
        fresh_source_to_retained_luminosity=powers['fresh_GR']['nuclear_minus_thermal_neutrino']/Linf,
        fresh_minus_retained_Lsun=powers['fresh_GR']['nuclear_minus_thermal_neutrino']-Linf,
        fresh_vs_frozen_nuclear_relative=powers['fresh_GR']['eps_nuc']/powers['frozen_source22']['eps_nuc']-1,
        fluorine_omission_at_source_nuclear_relative=powers['fresh_source19']['eps_nuc']/powers['frozen_source22']['eps_nuc']-1,
        peak=dict(zone=int(g['zone'][peak]),rho_B_g_cm3=float(g['rho'][peak]),T_K=float(10**g['logT'][peak]),
            eps_nuc_erg_g_s=float(g['eps_nuc'][peak]),reaction_neutrino_erg_g_s=float(g['eps_nuc_neu_total'][peak])),
        opacity_fresh_to_original=dict(all_cells_min=float(ratio.min()),all_cells_max=float(ratio.max()),
            burning_cells_unweighted_median=float(np.median(ratio[dense])),burning_cell_count=int(dense.sum())),
        interpretation='Fresh MESA net heat and neutrino/opacity evaluation on the declared GR material states. Luminosity and gravothermal/composition derivatives are not solved; this is an instantaneous source diagnostic. Do not call the mismatch an energy violation of an evolving star. MESA EOS auxiliaries are not yet certified against FreeEOS.',
        source_weights='All three comparisons use the same GR proper baryon cell masses and metric redshifts. Historical source powers keep their original verdict elsewhere.')
    save('power-comparison.json',result)
    derivatives=[];domain=abs(g['eps_nuc'])>1.;denom=np.maximum(1.,abs(g['eps_nuc']))
    coverage=float(np.dot(dm[domain],abs(g['eps_nuc'][domain]))/np.dot(dm,abs(g['eps_nuc'])))
    for prefix,field in [('T','d_lnepsnuc_dlnT'),('rho','d_lnepsnuc_dlnd')]:
        for step in [1e-3,5e-4]:
            _,p=mesa(OUT/(prefix+str(step)+'-plus-profile.data.gz'))
            _,n=mesa(OUT/(prefix+str(step)+'-minus-profile.data.gz'))
            fd=(p['eps_nuc']-n['eps_nuc'])/(2*step)/denom
            rel=abs(fd-g[field])/np.maximum(1.,abs(g[field]))
            index=np.flatnonzero(domain)[np.argmax(rel[domain])]
            derivatives.append(dict(variable=prefix,step=step,domain_cells=int(domain.sum()),
                covered_absolute_nuclear_power_fraction=coverage,max_normalized_discrepancy=float(rel[domain].max()),
                worst_zone=int(g['zone'][index]),analytic=float(g[field][index]),finite_difference=float(fd[index]),
                passed=bool(rel[domain].max()<1e-3)))
    save('derivative-control.json',dict(classification='Counterexample candidate',results=derivatives,
        all_passed=all(r['passed'] for r in derivatives),
        boundary='Sampled native MESA net derivatives including its EOS auxiliary response; not a global interval certificate or full joint GR-fluid derivative.'))
    print(json.dumps(result,indent=2),flush=True);print('Derivatives',json.dumps(derivatives,indent=2),flush=True)


def symbolic():
    import sympy as s
    P,rho,u,cx,c=s.symbols('P rho u cx c',positive=True)
    dP,drho,du=s.symbols('dP drho du')
    h=cx*c*c+u+P/rho;dh=du+dP/rho-P*drho/rho**2
    assert s.simplify(dh.subs(du,P*drho/rho**2)/h-dP/(rho*h))==0
    N,q,dB,dTau,dt=s.symbols('N q dB dTau dt',positive=True)
    assert s.simplify((N*q*dB*dTau).subs(dTau,N*dt)/dt-N**2*q*dB)==0
    save('symbolic.json',dict(classification='Proven',isentropic_cell_metric_enthalpy_identity=True,
        two_lapse_luminosity_weight=True,
        assumptions='Fixed X,s within each declared material cell for metric reconstruction; continuous pressure and metric across material interfaces. Frozen static background only for instantaneous redshifted source diagnostics.'))
    print('PASS symbolic metric enthalpy and luminosity redshifts')


def derivative_diagnostic_inputs():
    assert not (OUT/'derivative-diagnostic-plan.json').exists()
    failed=json.loads((OUT/'derivative-control.json').read_text());assert not failed['all_passed']
    base=dict(np.load(OUT/'gr-input.npz'));h,_=mesa(SOURCE);labels=[]
    for field,prefix in [('lnT','T'),('lnd','rho')]:
        for step in [1e-4,5e-5,2.5e-5]:
            for sign,suffix in [(-1,'minus'),(1,'plus')]:
                label=prefix+str(step)+'-'+suffix;data={k:v.copy() for k,v in base.items()};data[field]+=sign*step
                setup_run(label,data,float(h['star_age']));labels.append(label)
    setup_run('noscreen',base,float(h['star_age']),screening='');labels.append('noscreen')
    for field,prefix in [('lnT','T'),('lnd','rho')]:
        for sign,suffix in [(-1,'minus'),(1,'plus')]:
            label='noscreen-'+prefix+'-'+suffix;data={k:v.copy() for k,v in base.items()};data[field]+=sign*5e-5
            setup_run(label,data,float(h['star_age']),screening='');labels.append(label)
    save('derivative-diagnostic-plan.json',dict(classification='Counterexample candidate',labels=labels,
        frozen_primary_failure_sha256=sha(OUT/'derivative-control.json'),
        goals=['Test local finite-difference convergence at smaller steps without changing original failures',
               'Isolate screening contribution by comparing otherwise identical no-screening states'],
        local_FD_convergence_tolerance=1e-5,screening_diagnostic='Changing screening is a diagnostic model only, never a replacement physical fit.'))


def run_diagnostics():
    for label in json.loads((OUT/'derivative-diagnostic-plan.json').read_text())['labels']:
        run(label);collect(label)


def diagnostic_analysis():
    _,g=mesa(OUT/'gr-profile.data.gz');domain=abs(g['eps_nuc'])>1.;denom=np.maximum(1,abs(g['eps_nuc']))
    _,ng=mesa(OUT/'noscreen-profile.data.gz');result={};arrays={}
    for prefix,field in [('T','d_lnepsnuc_dlnT'),('rho','d_lnepsnuc_dlnd')]:
        values=[]
        for step in [1e-4,5e-5,2.5e-5]:
            _,p=mesa(OUT/(prefix+str(step)+'-plus-profile.data.gz'));_,n=mesa(OUT/(prefix+str(step)+'-minus-profile.data.gz'))
            values.append((p['eps_nuc']-n['eps_nuc'])/(2*step)/denom)
        last=values[-1];prev=values[-2];norm=np.maximum(1.,abs(last));conv=abs(last-prev)/norm
        error=abs(last-g[field])/np.maximum(1.,abs(g[field]))
        _,p=mesa(OUT/('noscreen-'+prefix+'-plus-profile.data.gz'));_,n=mesa(OUT/('noscreen-'+prefix+'-minus-profile.data.gz'))
        fd=(p['eps_nuc']-n['eps_nuc'])/1e-4/np.maximum(1,abs(ng['eps_nuc']))
        noerror=abs(fd-ng[field])/np.maximum(1,abs(ng[field]))
        result[prefix]=dict(local_convergence_max=float(conv[domain].max()),
            local_convergence_passed=bool(conv[domain].max()<1e-5),
            native_derivative_discrepancy_max=float(error[domain].max()),
            native_worst_zone=int(g['zone'][np.flatnonzero(domain)[np.argmax(error[domain])]]),
            noscreen_derivative_discrepancy_max=float(noerror[domain].max()),
            noscreen_worst_zone=int(g['zone'][np.flatnonzero(domain)[np.argmax(noerror[domain])]]))
        arrays[prefix]=last;arrays[prefix+'_previous']=prev;arrays[prefix+'_native']=g[field]
    np.savez_compressed(OUT/'local-nuclear-derivatives.npz',zone=g['zone'],domain=domain,
        d_eps_dlnT=arrays['T']*denom,d_eps_dlnrho=arrays['rho']*denom,**arrays)
    save('derivative-diagnostic.json',dict(classification='Counterexample candidate',results=result,
        original_native_derivative_gate_promoted=False,reference_state_sha256=sha(OUT/'gr-input.npz'),
        derivative_units='T and rho arrays are d eps_nuc / d ln(variable) divided by max(1,abs(eps_nuc)); d_eps_* arrays are erg/g/s per logarithmic variable.',
        interpretation='The local derivative arrays differentiate the actual screened evaluator at specified states. They are not a global error certificate; original native-derivative failures remain frozen.'))
    print('Derivative diagnostic',json.dumps(result,indent=2))


def weak_diagnostic():
    base=dict(np.load(OUT/'gr-input.npz'));h,_=mesa(SOURCE)
    save('weak-diagnostic-plan.json',dict(classification='Counterexample candidate',
        primary_failure_sha256=sha(OUT/'derivative-control.json'),
        diagnostic='Disable weak rates and screening together; compare T finite difference with native derivative to isolate raw strong-network contribution. This is not a replacement physical model.'))
    for label,step in [('noweak',0.),('noweak-T-minus',-5e-5),('noweak-T-plus',5e-5)]:
        data={k:v.copy() for k,v in base.items()};data['lnT']+=step
        setup_run(label,data,float(h['star_age']),screening='',weak_factor=0.)
        run(label);collect(label)
    weak_analysis()


def weak_analysis():
    _,g=mesa(OUT/'noweak-profile.data.gz');_,p=mesa(OUT/'noweak-T-plus-profile.data.gz');_,n=mesa(OUT/'noweak-T-minus-profile.data.gz')
    _,full=mesa(OUT/'gr-profile.data.gz');domain=abs(full['eps_nuc'])>1
    fd=(p['eps_nuc']-n['eps_nuc'])/1e-4/np.maximum(1,abs(g['eps_nuc']))
    error=abs(fd-g['d_lnepsnuc_dlnT'])/np.maximum(1,abs(g['d_lnepsnuc_dlnT']))
    _,noscreen=mesa(OUT/'noscreen-profile.data.gz')
    result=dict(classification='Counterexample candidate',max_discrepancy=float(error[domain].max()),
        worst_zone=int(g['zone'][np.flatnonzero(domain)[np.argmax(error[domain])]]),
        remaining_reaction_neutrino_max_erg_g_s=float(max(g['eps_nuc_neu_total'])),
        heating_change_max_erg_g_s=float(max(abs(g['eps_nuc']-noscreen['eps_nuc']))),
        raw_strong_isolation_valid=False,
        interpretation='The configured weak_rate_factor=0 did not remove all reaction neutrino power. This run is not a successful removal of every weak/compound channel and cannot isolate a raw strong-only network. No root-cause claim follows from the unchanged discrepancy.')
    save('weak-diagnostic.json',result);print(result,flush=True)


def local_derivatives(state_path):
    info=json.loads((OUT/'derivative-diagnostic.json').read_text())
    if sha(state_path)!=info['reference_state_sha256']:
        raise ValueError('Derivative arrays are bound to the frozen GR thermodynamic state')
    assert all(v['local_convergence_passed'] for v in info['results'].values())
    data=np.load(OUT/'local-nuclear-derivatives.npz');domain=data['domain']
    # Only expose the checked 896-cell domain, never unvalidated extrapolation.
    return {k:data[k][domain] for k in ['zone','d_eps_dlnT','d_eps_dlnrho']}


def eos_boundary():
    state=np.load(OUT/'gr-material-state.npz');_,g=mesa(OUT/'gr-profile.data.gz')
    ratio=g['pressure']/np.exp(state['logP']);domain=abs(g['eps_nuc'])>1
    result=dict(classification='Counterexample candidate',
        MESA_to_FreeEOS_pressure=dict(all_cells_min=float(min(ratio)),all_cells_max=float(max(ratio)),
            burning_cells_min=float(min(ratio[domain])),burning_cells_max=float(max(ratio[domain])),
            burning_cells_unweighted_median=float(np.median(ratio[domain]))),
        interpretation='MESA recomputes its EOS auxiliaries while FreeEOS defines the frozen TOV pressure. Equal rho_B,T,X does not prove equal EOS. Fresh rates/opacity are a declared hybrid state evaluation, not a fully common-EOS GR evolutionary solution.')
    save('EOS-boundary.json',result);print(result)


def source_bindings():
    names=['rates/private/rates_support.f90','rates/private/screen.f90','rates/private/screen5.f90',
        'rates/private/eval_weak.f90','net/private/net_screen.f90','net/private/net_derivs.f90',
        'net/private/net_derivs_support.f90','net/private/net_eval.f90','star/private/star_private_def.f90',
        'data/version_number','const/public/const_def.f90']
    for name in names:
        dest=OUT/'sources'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(MESA/name,dest)
    save('source-bindings.json',dict(classification='Imported from prior work',
        sha256={p.relative_to(OUT/'sources').as_posix():sha(p) for p in sorted((OUT/'sources').rglob('*')) if p.is_file()},
        runtime_data_version=(MESA/'data/version_number').read_text().strip(),
        source_archive='outputs/thermal-restart19/source-bindings.json',
        caution='Source inspection exposes possible blend-derivative omissions, but no specific source line has been shown to explain the measured binary temperature discrepancy. Native derivative failure is established by evaluator differences; root cause is not claimed.'))


def runtime_bindings():
    # Bind the actual microphysics data, including pre-existing rate caches.
    paths=[p for p in (MESA/'data').rglob('*') if p.is_file()]
    save('runtime-data-bindings.json',dict(classification='Proven',root=str(MESA/'data'),
        sha256={p.relative_to(MESA/'data').as_posix():sha(p) for p in sorted(paths)},
        files=len(paths),bytes=sum(p.stat().st_size for p in paths),
        note='Identity of the reused local data tree; not a claim every file was read or an independent build certificate.'))
    print('Runtime data bindings',len(paths),flush=True)


def maintain():
    prior=json.loads((OLD/'manifest.json').read_text())['sha256'];bindings={}
    docs=['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi']
    dest=OUT/'request24-notes';dest.mkdir(exist_ok=False)
    for rel in ['docs/'+k+'.md' for k in docs]+['paper/revision-manifest.json']:
        assert sha(ROOT/rel)==prior[rel],rel
        snap=dest/Path(rel).name;shutil.copy2(ROOT/rel,snap)
        bindings[rel]=dict(snapshot=snap.relative_to(ROOT).as_posix(),sha256=prior[rel],historical_manifest=rel.startswith('paper/'))
    save('historical-note-bindings.json',bindings)
    bodies=[
        '분류: Counterexample candidate. Request25는 기존 MESA 실행 파일의 모형 입력·초기 출력 경로로 새 GR 물질 상태의 핵 가열·반응 및 열적 중성미자·불투명도를 재평가했다. 출력의 밀도·온도·조성·구역 질량이 입력과 일치하는지 먼저 확인했다. FreeEOS가 정한 TOV 상태를 MESA 고유 EOS 보조량으로 평가한 혼합 모형이며 같은 EOS의 GR 진화 해는 아니다.',
        '분류: Counterexample candidate. 같은 GR 바리온 질량과 적색편이 가중치에서 재평가 핵 가열은 219.464816 Lsun으로 기존 값을 옮긴 156.792765 Lsun보다 약 39.9713% 높다. 열적 중성미자를 뺀 원천은 유지한 광도의 약 397.42배다. 이는 저장·팽창·조성 진화를 풀어야 한다는 순간 원천 진단이며 관측 광도 예측이나 에너지 보존 위반 판정이 아니다.',
        '분류: Proven. 구역 내부에서 조성·엔트로피가 고정되면 TOV의 계량 퍼텐셜은 dν=−d ln h를 만족한다. 이를 같은 수학적 대기와 연결해 적색편이를 계산했다. 국소 에너지와 고유 시간의 변환은 무한원 광도 원천에 e^(2ν)를 준다. 반응이 켜진 후의 엔트로피 보존이나 정확한 정적 열유속의 허용을 뜻하지 않는다.',
        '분류: Counterexample candidate. 기존 핵 가열 미분의 직접 차분 대조는 사전 문턱을 실패했다. 작은 간격의 실제 평가 함수 차분으로 896개 구역의 국소 미분 자료를 별도 구성했다. 이 영역은 절대 핵 가열의 약 99.999808%를 포함하며 특정 상태 SHA에 결속되어 다른 상태로 사용할 수 없다. 유한 차분 수렴은 전역 미분 오차 보장이나 미래 GR 동역학의 검증이 아니다.',
        '분류: Counterexample candidate. 기존 온도 미분과 실제 평가 함수 미분의 최대 차이는 작은 간격에서도 약 0.780312%로 남았다. 스크리닝 제거만으로 사라지지 않았고 weak_rate_factor=0도 반응 중성미자를 모두 제거하지 않아 순수 강반응 분리 대조로 인정하지 않는다. 근본 원인은 아직 특정하지 않았으며 기존 native 미분 통과로 승격하지 않는다.\n\n분류: Conjectural. 같은 입력 상태에서도 MESA/FreeEOS 압력 비는 전체 구역에서 약 0.99172–1.02410이다. 공통 EOS 보조량·전체 조성 변화율·미분 오차 인증과 GR 열수송/시간 적분은 미완료다. 정적 질량 매칭과 실제 반응 원천 평가를 완전한 진화·동적 관측량으로 세지 않는다.'
    ]
    revision=json.loads((ROOT/'paper/revision-manifest.json').read_text())
    for stem,body in zip(docs,bodies):
        rel='docs/'+stem+'.md'
        with (ROOT/rel).open('ab') as f:
            f.write(('\n\n## Request 25 새 GR 상태의 미세물리 평가\n\n'+body+'\n\n세부 근거: [한글 보고서](../notes/REQUEST25_FRESH_MICROPHYSICS_KO.md).\n').encode())
        revision['sha256'][rel]=sha(ROOT/rel)
    revision['request25_supporting_note_update']=dict(evidence_manifest='outputs/fresh-microphysics25/manifest.json',
        historical_notes='outputs/fresh-microphysics25/historical-note-bindings.json',
        status='새 GR 상태의 핵 가열·중성미자·불투명도 재평가; native 미분 실패 보존 및 별도 국소 미분 자료',
        artifact_status='원고 PDF와 ZIP은 Request 12 동결본이다.')
    (ROOT/'paper/revision-manifest.json').write_text(json.dumps(revision,ensure_ascii=False,indent=2)+'\n')


def seal():
    imports=sorted(OUT.glob('*-import-control.json'))
    assert len(imports)==31 and all(json.loads(p.read_text())['passed'] for p in imports)
    assert not json.loads((OUT/'derivative-control.json').read_text())['all_passed']
    assert all(r['local_convergence_passed'] for r in json.loads((OUT/'derivative-diagnostic.json').read_text())['results'].values())
    save('gates.json',dict(classification='Proven',imported_state_controls_passed=31,
        fresh_native_net_neutrino_opacity_evaluated=True,native_derivative_primary_gate_passed=False,
        separate_state_bound_local_derivative_data_passed=True,raw_strong_only_diagnostic_valid=False,
        common_EOS_auxiliaries_certified=False,full_composition_source_vector_exported=False,
        runtime_Q_to_rest_mass_energy_closure_certified=False,GR_heat_transport_solved=False,
        GR_thermal_composition_time_evolution_solved=False,global_derivative_error_certificate=False,
        full_nonlinear_observational_inference=False,final_PDF_or_ZIP_updated=False,
        classification_detail='Theorem progress: metric enthalpy/redshift identities. Loophole progress: fresh microphysics source evaluation and state-bound local derivative data; original native derivative failure retained.'))
    save('provenance.json',dict(classification='Proven',checkpoint='367bc1e',previous_manifest_sha256=sha(OLD/'manifest.json'),
        population='31 initial model profiles x 5735 zones = 177785 exported zone-states, not a count of internal kernel calls.',
        first_pilot='max_model_number=1 still allowed one evolution step after profile1; only the pre-evolution model-1 profile was selected. Later runs explicitly stop their own process after capturing a complete model-1 profile.',
        derivative_policy='Primary failure is bound before smaller-step/no-screening/weak-factor diagnostics. No threshold was changed to promote native derivatives. Local data are usable only on the checked state and domain.',
        dependencies='Reuse the frozen FreeEOS runtime and the exact MESA executable/data in runtime*.json. No toolchain installation or rebuild.',
        limitations='No proof of common MESA/FreeEOS EOS auxiliaries or of every actual reaction Q and composition flow. Runtime metadata says 7623 and is separate from the r7624-named downloaded package.'))
    files=[p for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json']
    files += [ROOT/'verification/fresh_microphysics.py',ROOT/'notes/REQUEST25_FRESH_MICROPHYSICS_KO.md',ROOT/'.gitattributes']
    files += [ROOT/k for k in json.loads((OUT/'historical-note-bindings.json').read_text())]
    save('manifest.json',dict(classification='Proven',sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(files)}))


def verify():
    stages=['validated-variational','remaining-levers15','nbody-readout16','nonzero-drive17','thermal-wd18',
        'thermal-restart19','thermal-robustness20','gr-mass21','thermal-closure22','baryon-entropy23','reactive-energy24','fresh-microphysics25']
    histories={f'outputs/{a}/manifest.json':f'outputs/{b}/historical-note-bindings.json' for a,b in zip(stages,stages[1:])}
    count=0
    for label in ['outputs/research-remediation/manifest.json',*histories,'paper/revision-manifest.json','outputs/fresh-microphysics25/manifest.json']:
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
    plan=json.loads((OUT/'derivative-diagnostic-plan.json').read_text())
    assert plan['frozen_primary_failure_sha256']==sha(OUT/'derivative-control.json')
    gates=json.loads((OUT/'gates.json').read_text())
    for k in ['native_derivative_primary_gate_passed','raw_strong_only_diagnostic_valid','common_EOS_auxiliaries_certified',
        'full_composition_source_vector_exported','runtime_Q_to_rest_mass_energy_closure_certified','GR_heat_transport_solved',
        'GR_thermal_composition_time_evolution_solved','global_derivative_error_certificate',
        'full_nonlinear_observational_inference','final_PDF_or_ZIP_updated']: assert not gates[k],k
    assert len(local_derivatives(OUT/'gr-input.npz')['zone'])==896
    try: local_derivatives(OUT/'source-input.npz')
    except ValueError: pass
    else: raise AssertionError('Wrong derivative state was accepted')
    assert sha(BINARY)==json.loads((OUT/'runtime.json').read_text())['sha256']
    data=json.loads((OUT/'runtime-data-bindings.json').read_text())
    for name,expected in data['sha256'].items(): assert sha(Path(data['root'])/name)==expected,name
    for key in ['runtime_sha256','module_sha256']:
        for path,digest in json.loads((ROOT/'outputs/gr-mass21/provenance.json').read_text())[key].items(): assert sha(path)==digest,path
    print('PASS',count,'artifact/history SHA;',len(data['sha256']),'MESA data SHA; derivative and evolution boundaries',flush=True)


def recheck():
    global OUT,CACHE
    original,cache=OUT,CACHE;worst=0.
    with tempfile.TemporaryDirectory(prefix='fresh-microphysics25-') as folder:
        OUT=Path(folder)/'out';CACHE=Path(folder)/'runs';CACHE.mkdir();shutil.copytree(original,OUT)
        try:
            for input_file in sorted(original.glob('*-inputs.json')):
                label=input_file.name[:-len('-inputs.json')];run_dir=CACHE/label;run_dir.mkdir()
                shutil.copy2(BINARY,run_dir/'binary')
                for p in (original/'inputs'/label).iterdir():
                    if p.name=='input.mod.gz':
                        with gzip.open(p,'rb') as src,(run_dir/'input.mod').open('wb') as dst: shutil.copyfileobj(src,dst)
                    else: shutil.copy2(p,run_dir/p.name)
                run(label);collect(label)
                _,new=mesa(OUT/(label+'-profile.data.gz'));_,old=mesa(original/(label+'-profile.data.gz'))
                for key in ['rho','logT',*ISOS,'eps_nuc','eps_nuc_neu_total','non_nuc_neu','opacity','d_lnepsnuc_dlnT','d_lnepsnuc_dlnd']:
                    error=float(np.max(abs(new[key]-old[key])/np.maximum(1e-30,abs(old[key]))));worst=max(worst,error)
                    assert error<1e-8,(label,key,error)
            # Recompute diagnostics from the fresh profiles, without promoting failures.
            analyze();diagnostic_analysis();weak_analysis();eos_boundary();symbolic()
            assert not json.loads((OUT/'derivative-control.json').read_text())['all_passed']
            assert all(r['local_convergence_passed'] for r in json.loads((OUT/'derivative-diagnostic.json').read_text())['results'].values())
        finally: OUT,CACHE=original,cache
    print('PASS 31 fresh runtime profiles; max relative repeat difference',worst,flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
