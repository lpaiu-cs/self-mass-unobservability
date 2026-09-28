"""Request 21: source-normalized thermal EOS and spherical GR mass matching.

Counterexample candidate: prescribed T(P), X(P) from an evolved thermal WD.
The EOS is evaluated afresh at each pressure. This is not a GR evolution run.
"""
from pathlib import Path
import ctypes, hashlib, json, shutil, subprocess, sys, time
import numpy as np
from scipy.interpolate import PchipInterpolator, interp1d
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, root
from thermal_wd import mesa
from thermal_restart import sha

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'outputs/gr-mass21'
CACHE = Path('/home/lpaiu/work/gr-mass21')
SOURCE = CACHE/'free_eos-3.0.0'
PROFILE = ROOT/'outputs/thermal-robustness20/runs/mesh/selected.data.gz'
G, C, GM_SUN, RSUN = 6.67428e-11, 299792458., 1.3271244e20, 695980000.
TARGET = .197536385307*GM_SUN/C**2
ELEMENTS = ['h','he','c','n','o','ne','na','mg','al','si','p','s','cl','ar','ca','ti','cr','mn','fe','ni']


def save(name, obj):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/name).write_text(json.dumps(obj, ensure_ascii=False, indent=2)+'\n')


def prepare():
    assert not (OUT/'plan.json').exists()
    save('plan.json', dict(classification='Counterexample candidate', checkpoint='b5e9028',
        target_GM_solar=.197536385307, mass_relative_tolerance=1e-8,
        source_profile=str(PROFILE.relative_to(ROOT)), profile_sha256=sha(PROFILE),
        eos='Unmodified FreeEOS 3.0.0 EOS1 (3,1,-2); all supported ionization stages.',
        thermal_closure='Prescribed log T(log P) and isotope number abundances(log P), interpolated from the mesh-refined evolved profile. Above the last central sample hold T and composition fixed. Re-evaluate rho(P,T,X) and u(P,T,X) with the EOS; solve the full TOV equations and proper baryon integral.',
        surface='Original first-zone pressure; retain original photospheric luminosity only for a conditional Teff diagnostic. No luminosity transport or GR evolution is solved.',
        composition='FreeEOS supports 20 elements but excludes fluorine. Primary model omits f17,f18,f19 and renormalizes the remaining baryon fractions; record exact removed abundance. Trace substitution controls are separate model sensitivity tests, not certified error bounds.',
        optical=dict(Teff_K=[15500,16100], logg_cgs=[5.67,5.97]),
        controls=['neutral dilute gas energy zero','isothermal first law','EOS density-pressure inversion',
                  'Newtonian same-EOS structure','GR pressure-coordinate versus radius-coordinate integration',
                  'table resolution and solver tolerance refinement'],
        claim_boundary='Numerical GR mass match for the declared EOS and thermal/composition closure only. This does not complete the original 22-isotope evolving rotating observed star or observational inference.'))
    (OUT/'sources').mkdir(exist_ok=True)
    shutil.copy2(CACHE/'free_eos-3.0.0.tar.gz',OUT/'sources/free_eos-3.0.0.tar.gz')
    save('source-bindings.json', dict(classification='Imported from prior work',
        url='https://sourceforge.net/projects/freeeos/files/freeeos/3.0.0%20Source/free_eos-3.0.0.tar.gz/download',
        archive_sha256=sha(OUT/'sources/free_eos-3.0.0.tar.gz'),
        downloaded_archive_signature_verified=False,
        upstream_source_unchanged=True))


def build():
    cmd=['gfortran','-O2','-fPIC','-shared','-I'+str(CACHE/'build/modules'),
         str(ROOT/'verification/gr_eos_bridge.f90'),'-L'+str(CACHE/'build/src'),
         '-Wl,-rpath,'+str(CACHE/'build/src'),'-lfree_eos','-o',str(CACHE/'gr_eos_bridge.so')]
    # Actual CMake module directory is discovered without changing the toolchain.
    module=next((CACHE/'build').rglob('mod_free_eos.mod')).parent
    cmd[4]='-I'+str(module)
    result=subprocess.run(cmd,capture_output=True,text=True)
    save('bridge-build.json',dict(command=cmd,returncode=result.returncode,stdout=result.stdout,stderr=result.stderr))
    assert result.returncode==0,result.stderr


class EOS:
    def __init__(self):
        self.lib=ctypes.CDLL(str(CACHE/'gr_eos_bridge.so'))
        self.call=self.lib.gr_eos
        self.call.argtypes=[ctypes.c_int,ctypes.c_double,ctypes.c_double,
            np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS'),ctypes.POINTER(ctypes.c_int)]
        self.call.restype=None

    def __call__(self, mode, logvalue, logT, eps):
        result=np.full(12,np.nan); info=ctypes.c_int(-999)
        self.call(mode,float(logvalue),float(logT),np.ascontiguousarray(eps),result,ctypes.byref(info))
        if info.value: raise ValueError((info.value,mode,logvalue,logT,eps.tolist()))
        assert np.all(np.isfinite(result)) and result[0]>0 and result[1]>0
        return result


class ThermalPath:
    def __init__(self, profile=PROFILE, fluorine='omit', temperature_scale=1.):
        self.header,d=mesa(profile)
        self.lp=d['logP']*np.log(10.)
        assert np.all(np.diff(self.lp)>0)
        self.lt=PchipInterpolator(self.lp,d['logT']*np.log(10.),extrapolate=False)
        isos=json.loads((ROOT/'outputs/thermal-restart19/network-diagnosis.json').read_text())['profile_isotopes']
        lines=(ROOT/'outputs/thermal-robustness20/sources/data/chem_data/isotopes.data').read_text().splitlines()[1:]
        weights={}
        for line in lines[::4]:
            w=line.split()
            if w: weights[w[0]]=(float(w[1]),int(w[2])+int(w[3]))
        import re
        abundance=np.array([d[k] for k in isos])
        self.removed=sum(d[k] for k in isos if k.startswith('f'))
        keep=np.array([not k.startswith('f') for k in isos])
        norm=abundance[keep].sum(axis=0)
        self.removed_mean=float(np.dot(self.removed,d['dm'])/sum(d['dm']))
        eps=np.zeros((len(self.lp),20)); cx=np.zeros(len(self.lp))
        for k,arr in zip(isos,abundance):
            element=re.sub('[0-9]','',k)
            if element=='f': continue
            w,a=weights[k]; x=arr/norm
            eps[:,ELEMENTS.index(element)]+=x/a
            cx+=x*w/a
        # Linear abundance interpolation preserves normalization and atom counts.
        self.cx=interp1d(self.lp,cx,bounds_error=True)
        self.ep=interp1d(self.lp,eps,axis=0,bounds_error=True)
        self.fluorine=fluorine
        self.temperature_scale=temperature_scale
        assert fluorine=='omit'

    def __call__(self, lp):
        x=np.clip(lp,self.lp[0],self.lp[-1])
        cx=float(self.cx(x)); ep=np.maximum(self.ep(x),0)/cx
        return float(self.lt(x))+np.log(self.temperature_scale),cx,np.ascontiguousarray(ep)


def probe():
    eos=EOS(); path=ThermalPath(); result=[]
    for i in np.linspace(0,len(path.lp)-1,16).astype(int):
        lp=path.lp[i]; lt,cx,eps=path(lp); t=time.monotonic()
        try:
            a=eos(1,lp,lt,eps)
            result.append(dict(index=int(i),logP=lp/np.log(10),logT=lt/np.log(10),
                logrho=np.log10(a[0]),energy=a[2],runtime=time.monotonic()-t))
        except ValueError as exc: result.append(dict(index=int(i),error=str(exc)))
    save('probe.json',dict(rows=result,omitted_F_baryon_fraction=path.removed_mean))
    print(json.dumps(result,indent=2))


def atmosphere_plan():
    assert not (OUT/'atmosphere-plan.json').exists()
    save('atmosphere-plan.json',dict(classification='Counterexample candidate',
        reason='A nonzero-pressure photosphere cannot be identified with an exact vacuum boundary.',
        atmosphere='Below the first-zone pressure Ps attach a gamma=5/3 polytrope, matching rho and total energy at Ps. P=Ps*z^(5/2), rho=rhos*z^(3/2), u=u0+3P/(2rho), u0=us-3Ps/(2rhos). Integrate z from 1 to 0. This explicitly specified atmosphere reaches P=rho=0 and permits a Schwarzschild exterior.',
        observables='Optical radius is the base of this mathematical atmosphere. Its mass and thickness are reported separately. Atmospheric transport and spectroscopy are not solved.',
        checks='Continuity at Ps, first-law identity, and atmosphere mass relative to the numerical matching tolerance.'))


class Star:
    def __init__(self, logpc, eos, path, relativistic=True, rtol=2e-8, r0=100., table=None, max_step=.08):
        self.eos,self.path,self.rtol,self.relativistic=eos,path,rtol,relativistic
        self.pc=np.exp(logpc); self.r0=r0; self.table=table
        p0,e0,b0,_=self.state(logpc)
        m0=4*np.pi*e0*r0**3/3
        pressure_drop=2*np.pi/3*(e0+p0)*(e0+3*p0)*r0*r0 if relativistic else 2*np.pi/3*e0*e0*r0*r0
        lp0=logpc+np.log1p(-pressure_drop/p0)
        self.logpc=logpc
        def rhs(lp,y):
            r,m,mb,nu=y
            p,e,b,_=self.state(lp)
            f=1-2*m/r if relativistic else 1.
            assert f>0
            grad=(e+p)*(m+4*np.pi*r**3*p)/(r*r*f) if relativistic else e*m/(r*r)
            dr=-p/grad
            return [dr,4*np.pi*r*r*e*dr,4*np.pi*r*r*b/np.sqrt(f)*dr,-p/(e+p) if relativistic else 0.]
        # Geometric units; nu is normalized at the vacuum boundary below.
        self.sol=solve_ivp(rhs,(lp0,path.lp[0]),[r0,m0,4*np.pi*b0*r0**3/3,0.],
            method='DOP853',rtol=rtol,atol=[1e-5,1e-12,1e-12,1e-16],max_step=max_step,dense_output=True)
        assert self.sol.success,self.sol.message
        self.photosphere=self.sol.y[:,-1].copy()
        ps,es,bs,us=self.state(path.lp[0])
        # Neutral-atom rest density in geometric units; finite-polytrope atmosphere.
        lt,cx,ep=path(path.lp[0]); a=eos(1,path.lp[0],lt,ep)
        rho_s=G*(a[0]*1000)/C**2
        u0=us/C**2-1.5*ps/rho_s
        self.atmosphere_eos=(ps,rho_s,bs,u0)
        def atmosphere(z,y):
            r,m,mb,nu=y
            p=ps*z**2.5; rho=rho_s*z**1.5
            e=rho*(1+u0)+1.5*p if relativistic else rho
            b=bs*z**1.5; f=1-2*m/r if relativistic else 1.
            dr=-2.5*(ps/rho_s)*r*r*f/((1+u0+2.5*ps/rho_s*z)*(m+4*np.pi*r**3*p)) if relativistic else -2.5*(ps/rho_s)*r*r/m
            return [dr,4*np.pi*r*r*e*dr,4*np.pi*r*r*b/np.sqrt(f)*dr,
                    -2.5*(ps/rho_s)/(1+u0+2.5*ps/rho_s*z) if relativistic else 0.]
        self.tail=solve_ivp(atmosphere,(1.,0.),self.photosphere,method='DOP853',
            rtol=rtol,atol=[1e-4,1e-12,1e-12,1e-15],max_step=.05,dense_output=True)
        assert self.tail.success
        self.radius,self.mass,self.baryon,self.nu=self.tail.y[:,-1]
        self.nu_shift=.5*np.log1p(-2*self.mass/self.radius)-self.nu

    def state(self,lp):
        lt,cx,ep=self.path(lp)
        a=self.eos(1,lp,lt,ep) if self.table is None else self.table(lp)
        rho=a[0]*1000; u=a[2]*1e-4
        p=G*np.exp(lp)*.1/C**4
        e=G*rho/C**2*(1+u/C**2) if self.relativistic else G*rho/C**2
        return p,e,G*rho/cx/C**2,u

    def record(self):
        r=self.photosphere[0]
        # Luminosity fixed at the photosphere, used solely as an optical diagnostic.
        original_r=float(self.path.header['photosphere_r'])*RSUN
        teff=float(self.path.header['Teff'])*np.sqrt(original_r/r)
        grav=self.mass*C**2/r**2/np.sqrt(1-2*self.mass/r) if self.relativistic else self.mass*C**2/r**2
        return dict(classification='Counterexample candidate',
            GR=self.relativistic,central_pressure_cgs=self.pc,
            mass_GM_solar=float(self.mass*C**2/GM_SUN),
            mass_relative_residual=float(self.mass/TARGET-1),
            proper_baryon_mass_GM_solar=float(self.baryon*C**2/GM_SUN),
            photospheric_radius_m=float(r),photospheric_radius_source_Rsun=float(r/RSUN),
            vacuum_radius_m=float(self.radius),
            atmosphere_mass_fraction=float((self.mass-self.photosphere[1])/self.mass),
            atmosphere_thickness_m=float(self.radius-r),
            conditional_Teff_K=float(teff),surface_logg_cgs=float(np.log10(grav*100)),
            optical_cut_passed=bool(15500<=teff<=16100 and 5.67<=np.log10(grav*100)<=5.97),
            rtol=self.rtol,r0_m=self.r0,nfev=int(self.sol.nfev+self.tail.nfev),
            closure='Prescribed T(P), X(P); independently evaluated FreeEOS; finite polytropic atmosphere; no GR evolution or transport.')


def match():
    assert (OUT/'atmosphere-plan.json').exists()
    eos=EOS(); path=ThermalPath(); trials=[]
    for gr in [False,True]:
        def objective(lp):
            star=Star(lp,eos,path,relativistic=gr)
            r=star.record();trials.append(r)
            print('trial',gr,lp,r['mass_GM_solar'],flush=True)
            return star.mass/TARGET-1
        fit=brentq(objective,path.lp[-1]-.15,path.lp[-1]+.15,xtol=2e-10)
        star=Star(fit,eos,path,relativistic=gr,rtol=2e-9,r0=50.)
        save(('gr' if gr else 'newtonian')+'-match.json',star.record())
        np.savez_compressed(OUT/(('gr' if gr else 'newtonian')+'-structure.npz'),
            logP=star.sol.t,r_m=star.sol.y[0],mass_geom_m=star.sol.y[1],
            baryon_geom_m=star.sol.y[2],nu=star.sol.y[3]+star.nu_shift,
            atmosphere_z=star.tail.t,atmosphere_state=star.tail.y)
    save('matching-trials.json',trials)
    print(json.dumps({k:json.loads((OUT/(k+'-match.json')).read_text()) for k in ['newtonian','gr']},indent=2))


def refinement_plan():
    save('refinement-plan.json',dict(classification='Counterexample candidate',
        pilot='Direct-call pilot results retained in gr-match.json and newtonian-match.json. Pilot optical cut failed and mass residual did not reach 1e-8 after tighter integration.',
        repair='Use linear abundance interpolation to preserve normalization; use fixed EOS samples and refine both EOS sampling and TOV integration. Retain pilot results, do not infer an ODE error certificate from solver tolerance.',
        original_temperature_scale=1., temperature_scale_fit='Not part of the primary frozen-temperature test.',
        grid='All original pressure knots, subdivisions of each interval, and a central extension of 0.4 in ln P with frozen central T and composition.',
        model_source_unchanged=True))


def make_table(path,eos,subdivision=1):
    lp=path.lp
    grid=np.unique(np.concatenate([lp[:-1]+j/subdivision*np.diff(lp) for j in range(subdivision)]+
                          [np.linspace(lp[-1],lp[-1]+.4,101*subdivision)]))
    values=[]
    for p in grid:
        lt,cx,ep=path(p); values.append(eos(1,p,lt,ep))
    values=np.array(values)
    assert np.all(values[:,4]>0) and np.all(values[:,7]>0)
    # Pressure and density are positive; logarithmic interpolation preserves it.
    return grid,values,table_interpolant(grid,values)


def table_interpolant(grid,values):
    log_density=PchipInterpolator(grid,np.log(values[:,0]),extrapolate=False)
    rest=PchipInterpolator(grid,values[:,1:],axis=0,extrapolate=False)
    def table(p): return np.r_[np.exp(log_density(p)),rest(p)]
    return table


def refined():
    assert (OUT/'refinement-plan.json').exists()
    eos=EOS(); path=ThermalPath(); rows=[]
    for sub in [1,2,4]:
        t=time.monotonic(); grid,values,table=make_table(path,eos,sub)
        np.savez_compressed(OUT/f'eos-grid-{sub}.npz',logP=grid,values=values)
        for gr in [False,True]:
            def objective(lp):
                return Star(lp,eos,path,gr,rtol=2e-11,r0=20.,table=table,max_step=.03).mass/TARGET-1
            lp=brentq(objective,path.lp[-1]-.15,path.lp[-1]+.15,xtol=3e-12)
            star=Star(lp,eos,path,gr,rtol=2e-12,r0=10.,table=table,max_step=.015)
            row=star.record();row['subdivision']=sub;row['table_rows']=len(grid)
            rows.append(row)
            name=('gr' if gr else 'newtonian')+f'-refined-{sub}'
            save(name+'.json',row)
            np.savez_compressed(OUT/(name+'-structure.npz'),logP=star.sol.t,state=star.sol.y,
                atmosphere_z=star.tail.t,atmosphere_state=star.tail.y,nu_shift=star.nu_shift)
            print(name,row['mass_relative_residual'],row['photospheric_radius_m'],flush=True)
        print('elapsed',sub,time.monotonic()-t,flush=True)
    save('refined-matching.json',rows)


def calibration_plan():
    assert not (OUT/'calibration-plan.json').exists()
    save('calibration-plan.json',dict(classification='Counterexample candidate',
        registered_after='Frozen-temperature optical failure in the pilot and subdivision-1 refinement; the original failure remains.',
        parameters='Central pressure and a uniform multiplier of the prescribed T(P); composition pressure path unchanged.',
        bounds=dict(temperature_multiplier=[.9,1.1],central_logP_offset=[-.4,.4]),
        targets='ADM GM mass and Teff=15800 K conditional on the retained MESA photospheric luminosity.',
        meaning='A constructed two-parameter equilibrium candidate. The fitted temperature and radius are not independent predictions, evidence for formation, or an observational posterior.',
        gates='Both residuals below 1e-8 within the sampled EOS model; refine EOS samples and use independent radial TOV integration. Preserve any failed calibration.'))


def calibrate():
    assert (OUT/'calibration-plan.json').exists()
    eos=EOS(); base=ThermalPath(); trials=[]
    radius_target=float(base.header['photosphere_r'])*RSUN*(float(base.header['Teff'])/15800.)**2
    x0=[np.log(json.loads((OUT/'gr-refined-4.json').read_text())['central_pressure_cgs']),0.]
    for sub in [2,4]:
        def objective(x):
            pc,ts=x
            assert abs(pc-base.lp[-1])<.4 and .9<np.exp(ts)<1.1,x
            path=ThermalPath(temperature_scale=np.exp(ts))
            grid,values,table=make_table(path,eos,sub)
            star=Star(pc,eos,path,rtol=1e-11,r0=10.,table=table,max_step=.02)
            error=[np.log(star.mass/TARGET),np.log(star.photosphere[0]/radius_target)]
            trials.append(dict(subdivision=sub,logpc=float(pc),temperature_scale=float(np.exp(ts)),residual=error))
            save('calibration-trials.json',trials)
            print('calibrate',sub,x,error,flush=True)
            return error
        fit=root(objective,x0,method='hybr',options={'eps':1e-8,'xtol':3e-9})
        residual=objective(fit.x)
        assert max(abs(np.array(residual)))<1e-8,(fit.message,residual)
        path=ThermalPath(temperature_scale=np.exp(fit.x[1]))
        grid,values,table=make_table(path,eos,sub)
        star=Star(fit.x[0],eos,path,rtol=2e-12,r0=10.,table=table,max_step=.015)
        row=star.record(); row.update(temperature_scale=float(np.exp(fit.x[1])),
            target_radius_m=radius_target,subdivision=sub,root_success=bool(fit.success),
            fitted_Teff_is_independent_validation=False)
        save(f'calibrated-{sub}.json',row)
        np.savez_compressed(OUT/f'calibrated-eos-{sub}.npz',logP=grid,values=values)
        np.savez_compressed(OUT/f'calibrated-structure-{sub}.npz',logP=star.sol.t,state=star.sol.y,
            atmosphere_z=star.tail.t,atmosphere_state=star.tail.y,nu_shift=star.nu_shift)
        x0=fit.x


def source_audit():
    names=['src/mod_free_eos.f90','src/free_eos_detailed.f90','src/ionize.f90',
           'src/mod_free_eos_constants.f90','src/mod_ionization_data.f90',
           'src/mod_isotopic_mass_data.f90','COPYING','INSTALL','README.release']
    files={}
    import tarfile
    with tarfile.open(OUT/'sources/free_eos-3.0.0.tar.gz') as archive:
        for name in names:
            content=archive.extractfile('free_eos-3.0.0/'+name).read()
            assert content==(SOURCE/name).read_bytes(),name
            dest=OUT/'sources'/name;dest.parent.mkdir(exist_ok=True,parents=True);dest.write_bytes(content)
            files[dest.relative_to(ROOT).as_posix()]=sha(dest)
    for name in ['configure.log','configure-test-on.log','build.log','match.log','refined.log']:
        if (CACHE/name).exists(): shutil.copy2(CACHE/name,OUT/name)
    save('source-audit.json',dict(classification='Proven',exact_archive_member_sha256=files,
        energy_reference='Non-H neutral electronic ground states; H2 ground state for hydrogen. Wrapper subtracts exactly 0.5*c2*cr*h2diss*eps_H to use neutral atoms for all species.',
        total_energy_cgs='epsilon = rho_atomic * (c^2 + u_neutral)',
        proper_baryon_density_cgs='rho_B = rho_atomic / C_X; C_X=sum X_baryon_i*W_i/A_i',
        no_additional_atomic_mass_factor_on_rho_atomic=True,
        electron_rest_mass_not_added_again=True,
        provenance='EOS itself remains unmodified; only its documented API is wrapped.'))


def symbolic():
    import sympy as s
    rho,c,u,d,p,r,m,z,ps,rs,u0=s.symbols('rho c u d p r m z ps rs u0',positive=True)
    assert s.expand(rho*(c*c-d+u+d)-rho*(c*c+u))==0
    pp=ps*z**s.Rational(5,2); rr=rs*z**s.Rational(3,2)
    uu=u0+s.Rational(3,2)*pp/rr
    assert s.simplify(s.diff(uu,z)-pp/rr**2*s.diff(rr,z))==0
    f=1-2*m/r; e=s.symbols('e',positive=True)
    pr=-(e+p)*(m+4*s.pi*r**3*p)/(r*r*f)
    rp=-p*r*r*f/((e+p)*(m+4*s.pi*r**3*p))
    assert s.simplify(pr/p*rp-1)==0
    assert s.simplify((s.symbols('nu_r')+0).subs('nu_r',(m+4*s.pi*r**3*p)/(r*r*f))*rp+p/(e+p))==0
    save('symbolic-audit.json',dict(classification='Proven',energy_reference_shift_invariance=True,
        polytropic_atmosphere_first_law=True,TOV_log_pressure_coordinate=True,
        metric_time_potential_chain_rule=True))


def eos_checks():
    eos=EOS(); path=ThermalPath(); rows=[]
    for i in np.linspace(0,len(path.lp)-1,31).astype(int):
        lp=path.lp[i];lt,cx,ep=path(lp)
        a=eos(1,lp,lt,ep); b=eos(2,np.log(a[0]),lt,ep)
        # Maxwell relation at fixed composition; reference offset has no rho derivative.
        expected=b[1]/b[0]*(1-b[6])
        maxwell=abs(b[9]-expected)/max(abs(b[9]),abs(expected),1.)
        rows.append(dict(index=int(i),pressure_inverse_error=float(abs(b[1]/np.exp(lp)-1)),
                         density_inverse_error=float(abs(a[0]/b[0]-1)),maxwell_relative_error=float(maxwell)))
    ep=np.zeros(20);ep[1]=1/4.00260325413
    a=eos(2,np.log(1e-8),np.log(1000.),ep)
    k,h,c,na=1.380649e-16,6.62607015e-27,C*100,6.02214076e23
    arad=8*np.pi**5*k**4/(15*h**3*c**3)
    ideal=1.5*k*na*1000*ep[1]+arad*1000**4/1e-8
    neutral_error=abs(a[2]/ideal-1)
    assert max(x['pressure_inverse_error'] for x in rows)<1e-7
    assert max(x['maxwell_relative_error'] for x in rows)<1e-6
    assert neutral_error<1e-5
    save('eos-checks.json',dict(classification='Proven',scope='specified source implementation and sampled states',
        rows=rows,neutral_helium_energy_relative_error=float(neutral_error),
        maximum_omitted_fluorine_baryon_fraction=float(path.removed.max()),
        mass_weighted_omitted_fluorine_baryon_fraction=path.removed_mean,
        thermodynamic_error_certificate=False))
    print('PASS: source energy reference, EOS inversion and sampled Maxwell controls')


def independent_radius(star, direct=False):
    """Independent r-coordinate TOV solve of the same declared model."""
    path=star.path; eos=star.eos
    def state(lp):
        lp=max(lp,path.lp[0])
        lt,cx,ep=path(lp)
        a=eos(1,lp,lt,ep) if direct else star.table(lp)
        p=G*np.exp(lp)*.1/C**4
        rho=G*a[0]*1000/C**2
        return p,rho*(1+a[2]*1e-4/C**2),rho/cx
    r0=15.;pc,ec,bc=state(star.logpc)
    p0=pc-2*np.pi/3*(ec+pc)*(ec+3*pc)*r0*r0
    y0=[4*np.pi*ec*r0**3/3,np.log(p0*C**4/G*10),4*np.pi*bc*r0**3/3]
    def rhs(r,y):
        m,lp,mb=y;p,e,b=state(lp);f=1-2*m/r
        return [4*np.pi*r*r*e,-(e+p)/p*(m+4*np.pi*r**3*p)/(r*r*f),4*np.pi*r*r*b/np.sqrt(f)]
    def surface(r,y): return y[1]-path.lp[0]
    surface.terminal=True;surface.direction=-1
    sol=solve_ivp(rhs,(r0,1.1*star.photosphere[0]),y0,method='DOP853',
        events=surface,rtol=2e-11,atol=[1e-12,1e-13,1e-12],max_step=star.photosphere[0]/3000)
    assert sol.success and len(sol.t_events[0])==1
    r=sol.t[-1];m,lp,b=sol.y[:,-1]
    return dict(radius_m=float(r),mass_geom_m=float(m),baryon_geom_m=float(b),nfev=sol.nfev,
        relative_difference_radius=float(r/star.photosphere[0]-1),
        relative_difference_mass=float(m/star.photosphere[1]-1),
        relative_difference_baryon=float(b/star.photosphere[2]-1),direct_EOS=direct)


def uniform_control():
    # Independent closed-form Schwarzschild interior positive control.
    mass=.2; energy=3*mass/(4*np.pi); r0=1e-5
    def pressure(r):
        v=np.sqrt(1-2*mass*r*r);s=np.sqrt(1-2*mass)
        return energy*(v-s)/(3*s-v)
    def rhs(r,y):
        m,p=y
        return [4*np.pi*r*r*energy,-(energy+p)*(m+4*np.pi*r**3*p)/(r*r*(1-2*m/r))]
    sol=solve_ivp(rhs,(r0,1.),[mass*r0**3,pressure(r0)],method='DOP853',
        rtol=2e-12,atol=1e-14,dense_output=True,max_step=.002)
    points=np.linspace(r0,1.,1001);vals=sol.sol(points)
    error=float(np.max(abs(vals[1]-pressure(points)))/pressure(0.))
    assert error<1e-9 and abs(vals[0,-1]/mass-1)<1e-11
    return dict(classification='Proven',uniform_density_TOV_relative_pressure_error=error)


def audit():
    calibrated=len(sys.argv)>2 and sys.argv[2]=='calibrated'
    name='calibrated' if calibrated else 'primary'
    record=json.loads((OUT/('calibrated-final.json' if calibrated else 'gr-refined-4.json')).read_text())
    eos=EOS();path=ThermalPath(temperature_scale=record.get('temperature_scale',1.))
    data=np.load(OUT/('calibrated-eos-final.npz' if calibrated else 'eos-grid-4.npz'))
    table=table_interpolant(data['logP'],data['values'])
    star=Star(np.log(record['central_pressure_cgs']),eos,path,rtol=1e-12,r0=5.,table=table,max_step=.01)
    radial=independent_radius(star)
    direct=independent_radius(star,direct=True)
    control=uniform_control()
    sample=np.linspace(path.lp[0],star.logpc,67)
    density_errors=[]
    for p in sample:
        lt,cx,ep=path(p);real=eos(1,p,lt,ep);cached=table(p)
        density_errors.append(abs(real[0]/cached[0]-1))
    result=dict(classification='Proven',scope='Numerical cross-checks of the declared closure, not a continuum certificate.',
        mass_residual_refined=star.mass/TARGET-1,
        pressure_coordinate_record=star.record(),independent_radial=radial,
        independent_radial_direct_EOS=direct,uniform_sphere=control,
        maximum_sampled_table_density_difference=float(max(density_errors)),
        numerical_mass_tolerance_passed=bool(abs(star.mass/TARGET-1)<1e-8),
        direct_EOS_mass_tolerance_passed=bool(abs(direct['mass_geom_m']/TARGET-1)<1e-8))
    save(name+'-audit.json',result)
    assert max(abs(radial[k]) for k in ['relative_difference_radius','relative_difference_mass'])<1e-7
    assert abs(direct['relative_difference_mass'])<1e-7
    print(json.dumps(result,indent=2))


def polish():
    """Reuse the already measured calibration Jacobian at tighter TOV tolerance."""
    prior=[x for x in json.loads((OUT/'calibration-trials.json').read_text()) if x['subdivision']==4]
    def vector(row): return np.array([row['logpc'],np.log(row['temperature_scale'])])
    base=prior[0]; xb=vector(base); eb=np.array(base['residual']); columns=[]
    for j in range(2):
        row=next(x for x in prior if abs(vector(x)[j]-xb[j])>1e-6 and abs(vector(x)[1-j]-xb[1-j])<1e-14)
        columns.append((np.array(row['residual'])-eb)/(vector(row)[j]-xb[j]))
    jac=np.array(columns).T
    record=json.loads((OUT/'calibrated-4.json').read_text())
    x=np.array([np.log(record['central_pressure_cgs']),np.log(record['temperature_scale'])])
    eos=EOS(); trials=[]
    for _ in range(5):
        path=ThermalPath(temperature_scale=np.exp(x[1]));grid,values,table=make_table(path,eos,4)
        star=Star(x[0],eos,path,rtol=2e-13,r0=5.,table=table,max_step=.008)
        error=np.array([np.log(star.mass/TARGET),np.log(star.photosphere[0]/record['target_radius_m'])])
        trials.append(dict(logpc=float(x[0]),temperature_scale=float(np.exp(x[1])),residual=error.tolist()))
        print('polish',trials[-1],flush=True)
        if max(abs(error))<2e-9: break
        x-=np.linalg.solve(jac,error)
    save('polish-trials.json',dict(jacobian=jac.tolist(),trials=trials))
    assert max(abs(error))<1e-8,error
    result=star.record();result.update(temperature_scale=float(np.exp(x[1])),target_radius_m=record['target_radius_m'],
        radius_relative_residual=float(error[1]),fitted_Teff_is_independent_validation=False)
    save('calibrated-final.json',result)
    np.savez_compressed(OUT/'calibrated-eos-final.npz',logP=grid,values=values)
    np.savez_compressed(OUT/'calibrated-structure-final.npz',logP=star.sol.t,state=star.sol.y,
        atmosphere_z=star.tail.t,atmosphere_state=star.tail.y,nu_shift=star.nu_shift)


def finalize():
    primary=json.loads((OUT/'primary-audit.json').read_text())
    calibrated=json.loads((OUT/'calibrated-audit.json').read_text())
    fitted=json.loads((OUT/'calibrated-final.json').read_text())
    newton=json.loads((OUT/'newtonian-refined-4.json').read_text())
    gr=json.loads((OUT/'gr-refined-4.json').read_text())
    gates=dict(classification='Counterexample candidate',
        source_normalized_thermal_EOS=True,full_spherical_TOV_and_proper_baryon_integral=True,
        zero_pressure_vacuum_boundary=True,
        fixed_temperature_mass_match=primary['direct_EOS_mass_tolerance_passed'],
        fixed_temperature_optical_pass=gr['optical_cut_passed'],
        separately_calibrated_mass_match=calibrated['direct_EOS_mass_tolerance_passed'],
        separately_calibrated_optical_candidate=fitted['optical_cut_passed'],
        original_22_isotope_GR_evolution=False,rotation_and_thermal_transport_solved=False,
        stellar_stability_certified=False,full_EOS_or_TOV_error_certificate=False,
        complete_force_derivative_certificate=False,complete_nonlinear_observational_inference=False,
        final_manuscript_PDF_and_ZIP_updated=False,
        result_class='Theorem progress: energy-reference and proper-volume closure. Loophole progress: constructed thermal GR equilibrium candidate, not an observed dynamic loophole.')
    assert gates['fixed_temperature_mass_match'] and not gates['fixed_temperature_optical_pass']
    assert gates['separately_calibrated_mass_match'] and gates['separately_calibrated_optical_candidate']
    save('gates.json',gates)
    save('newtonian-GR-comparison.json',dict(classification='Counterexample candidate',same_EOS_and_thermal_path=True,
        photospheric_radius_difference_m=gr['photospheric_radius_m']-newton['photospheric_radius_m'],
        conditional_Teff_difference_K=gr['conditional_Teff_K']-newton['conditional_Teff_K'],
        warning='The much larger difference from the original evolved MESA candidate includes EOS replacement and thermal/composition closure changes. It is not attributable to GR alone.'))
    for name in ['calibrate.log','polish.log','primary-audit.log','calibrated-audit.log']:
        shutil.copy2(CACHE/name,OUT/name)


def maintain():
    old=json.loads((ROOT/'outputs/thermal-robustness20/manifest.json').read_text())['sha256']
    final=json.loads((OUT/'calibrated-final.json').read_text())
    additions={
        'model-definition':f"분류: Counterexample candidate. FreeEOS 3.0.0 EOS1을 실제로 재평가하는 구대칭 TOV 모형을 구성했다. 중성 원자 에너지 기준, 고유 바리온 부피, 영압 외곽을 명시했다. 고정 T(P)·조성 조건에서는 GR 질량 일치만 통과하고 광학 조건은 실패했다. 별도 온도 배율 {final['temperature_scale']:.12f}의 후보는 목표 질량과 조건부 Teff를 맞춘다. 불소를 제외해 재규격화한 조성, 정해 둔 T(P), 수학적 외곽 대기는 모형 가정이다.",
        'observable-targets':f"분류: Counterexample candidate. 별도 조정 후보의 GR 질량은 {final['mass_GM_solar']:.12f} GM_sun 단위, 반지름은 {final['photospheric_radius_source_Rsun']:.10f} 배포 R_sun, 조건부 Teff는 {final['conditional_Teff_K']:.6f} K, logg는 {final['surface_logg_cgs']:.8f}다. Teff와 반지름은 유지한 모형 광도를 조건으로 조정에 사용한 값이며 독립 관측 예측이 아니다. 두 적분 좌표와 직접 EOS 재평가가 질량 수치 허용오차를 통과했다.",
        'adiabatic-limit':"분류: Proven. 정지 에너지와 내부에너지의 반대 방향 기준 이동은 총 에너지 밀도를 보존한다. FreeEOS의 H2 영점 이동을 소스의 정확한 상수로 제거하고 중성 원자 질량을 더했다. 압력 좌표 TOV 변환과 외곽 폴리트로프의 제1법칙도 기호적으로 검증했다. 기존 평탄·고정 구각 scalar 오차 상계는 이 새 GR 유체 모형의 오차 보장이 아니다.",
        'nonadiabatic-regime':"분류: Conjectural. 이번 계산은 열 EOS를 가진 GR 평형의 질량·부피 경계를 진전시켰다. 이 구조의 유체·metric·scalar 결합 동역학과 궤도 주파수의 관측 전달함수는 아직 계산하지 않았다. 정적 GR 질량 일치를 새로운 위상 지연이나 동적 관측량의 증명으로 세지 않는다.",
        'failure-ledger-dynamic-chi':"분류: Counterexample candidate. 고정 온도 경로의 실제 EOS/GR 재구성은 약 16566 K로 원래 광학 컷을 실패했다. 별도 온도 배율 조정으로 얻은 후보는 그 실패를 대체하지 않는다. 같은 EOS의 Newtonian–GR 비교에서 온도 차이는 약 3.24 K이므로 원래 MESA 후보와의 전체 차이를 GR 효과로 해석하지 않는다.\n\n분류: Conjectural. 남은 경계는 원래 22개 동위원소 전체와 진화·열수송·회전을 일관되게 연결하는 실제 항성 모형, 안정성, 유체·metric·scalar 동역학, 전구간 미분 오차와 전체 관측 비선형 추론이다. 이번 수치 질량 허용오차는 EOS 정확도나 관측 질량 오차 인증이 아니다."}
    dest=OUT/'request20-notes';dest.mkdir(exist_ok=False);bindings={}
    for rel in ['docs/'+k+'.md' for k in additions]+['paper/revision-manifest.json']:
        assert sha(ROOT/rel)==old[rel],rel
        snap=dest/Path(rel).name;shutil.copy2(ROOT/rel,snap)
        bindings[rel]=dict(snapshot=snap.relative_to(ROOT).as_posix(),sha256=old[rel],historical_manifest=rel.startswith('paper/'))
    save('historical-note-bindings.json',bindings)
    revision=json.loads((ROOT/'paper/revision-manifest.json').read_text())
    for stem,body in additions.items():
        rel='docs/'+stem+'.md'
        with (ROOT/rel).open('ab') as f:
            f.write(('\n\n## Request 21 열 EOS와 GR 질량 매칭\n\n'+body+'\n\n세부 근거: [한글 보고서](../notes/REQUEST21_GR_MASS_MATCHING_KO.md).\n').encode())
        revision['sha256'][rel]=sha(ROOT/rel)
    revision['request21_supporting_note_update']=dict(evidence_manifest='outputs/gr-mass21/manifest.json',
        historical_notes='outputs/gr-mass21/historical-note-bindings.json',
        status='명시한 열 EOS·열 경로의 GR 질량 매칭 및 별도 광학 조정 후보; 실제 GR 진화·전체 인증·추론 미완료',
        artifact_status='원고 PDF와 ZIP은 Request 12 동결본이다.')
    (ROOT/'paper/revision-manifest.json').write_text(json.dumps(revision,ensure_ascii=False,indent=2)+'\n')


def seal():
    import scipy, mpmath, sympy
    save('provenance.json',dict(before_task_checkpoint='b5e9028',
        interpreter=sys.executable,versions={m.__name__:m.__version__ for m in [np,scipy,mpmath,sympy]},
        producer='verification/gr_mass.py',bridge='verification/gr_eos_bridge.f90',
        runtime_sha256={str(p):sha(p) for p in [CACHE/'gr_eos_bridge.so',CACHE/'build/src/libfree_eos.so.1.0.0']},
        module_sha256={str(Path(m.__file__)):sha(m.__file__) for m in [np,scipy,mpmath,sympy]},
        source_archive_sha256=sha(OUT/'sources/free_eos-3.0.0.tar.gz'),
        previous_manifest_sha256=sha(ROOT/'outputs/thermal-robustness20/manifest.json')))
    files=[p for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json']
    files += [ROOT/'verification/gr_mass.py',ROOT/'verification/gr_eos_bridge.f90',
              ROOT/'notes/REQUEST21_GR_MASS_MATCHING_KO.md',ROOT/'.gitattributes']
    files += [ROOT/k for k in json.loads((OUT/'historical-note-bindings.json').read_text())]
    save('manifest.json',dict(classification='Proven',scope='Source-normalized thermal GR mass matching and preserved failure/calibration boundaries',
        sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(files)}))


def verify():
    histories={
        'outputs/validated-variational/manifest.json':'outputs/remaining-levers15/historical-note-bindings.json',
        'outputs/remaining-levers15/manifest.json':'outputs/nbody-readout16/historical-note-bindings.json',
        'outputs/nbody-readout16/manifest.json':'outputs/nonzero-drive17/historical-note-bindings.json',
        'outputs/nonzero-drive17/manifest.json':'outputs/thermal-wd18/historical-note-bindings.json',
        'outputs/thermal-wd18/manifest.json':'outputs/thermal-restart19/historical-note-bindings.json',
        'outputs/thermal-restart19/manifest.json':'outputs/thermal-robustness20/historical-note-bindings.json',
        'outputs/thermal-robustness20/manifest.json':'outputs/gr-mass21/historical-note-bindings.json'}
    count=0
    for label in ['outputs/research-remediation/manifest.json',*histories,'paper/revision-manifest.json','outputs/gr-mass21/manifest.json']:
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
            assert sha(path)==expected,(label,name)
            count+=1
    provenance=json.loads((OUT/'provenance.json').read_text())
    for key in ['runtime_sha256','module_sha256']:
        for path,digest in provenance[key].items(): assert sha(path)==digest,path
    gates=json.loads((OUT/'gates.json').read_text())
    assert not gates['fixed_temperature_optical_pass'] and not gates['original_22_isotope_GR_evolution']
    assert gates['fixed_temperature_mass_match'] and gates['separately_calibrated_mass_match']
    assert not gates['complete_nonlinear_observational_inference']
    print('PASS:',count,'현재·역사 SHA, 동결 원고 및 실패/조정/미완료 경계')


if __name__=='__main__':
    globals()[sys.argv[1]]()
