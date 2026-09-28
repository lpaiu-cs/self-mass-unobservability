"""Use native EOS1 without LTE photons before evolving radiation independently."""
import ctypes, json, shutil, subprocess, sys
import numpy as np
import sympy as sp
import gr_molecular_reference_runner as reference

g=reference.g;model=reference.original.model;OUT=g.OUT/'gr-radiation-eos-split'
CACHE=g.CACHE/'radiation-eos-split';BRIDGE=CACHE/'gas.so'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists() and not CACHE.exists();OUT.mkdir();CACHE.mkdir();reference.verify()
    before=(model.OUT/'direct_ion_bridge.f90').read_text()
    old='call free_eos(0,3,1,-2,kif';new='call free_eos(0,3,11,-2,kif'
    assert before.count(old)==1;after=before.replace(old,new);assert after.replace(new,old)==before
    after+='''
subroutine photon_constants(values) bind(C)
  use iso_c_binding
  use mod_free_eos_constants, only: prad_const,clight
  implicit none
  real(c_double),intent(out) :: values(2)
  values=[3._c_double*prad_const,clight]
end subroutine photon_constants
'''
    (OUT/'gas-bridge.f90').write_text(after)
    for name in ['mod_free_eos.f90','mod_free_eos_constants.f90']:
        shutil.copy2(model.CACHE/'source/src'/name,OUT/name)
    paths=[g.ROOT/'verification/gr_radiation_eos_split.py',model.OUT/'manifest.json',reference.OUT/'manifest.json',
        OUT/'gas-bridge.f90',OUT/'mod_free_eos.f90',OUT/'mod_free_eos_constants.f90',g.OUT/'gr-transport/plan.json']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='d79466b',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        cells=[0,1175,1176,2972,5734],temperature_offsets=[-.001,0.,.001],
        decomposition_relative_tolerance=1e-10,inverse_log_density_tolerance=1e-10,
        gamma_identity_relative_tolerance=1e-4,
        intervention='Call the same frozen two-spectrum library through its existing EOS1 ifmodified=11 option, which selects ifrad=0. The original EOS1 option is ifmodified=1. The library/source, molecular flag and physical constants are not rebuilt or changed. Read the compiled radiation constant instead of substituting another constants convention.',
        identities='f_total=f_gas-a_rad*T^4/(3*rho_B); P_total=P_gas+a_rad*T^4/3; u_total=u_gas+a_rad*T^4/rho_B; s_total=s_gas+4*a_rad*T^3/(3*rho_B). Subtract their derivatives consistently.',
        controls='At five real states and three temperature offsets compare native gas output with analytic photon subtraction, pressure inverse and alternating gas/total call history. Retain the existing 2 erg/g or 32 ulp energy/entropy gate. Retain the previous GR-transport 1e-4 native gamma-identity criterion separately.',
        full_reference='Apply the algebraic decomposition to all 5735 saved new-EOS density/temperature evaluations; report gas pressure/capacity positivity. This full-grid algebraic pass is distinct from the 15 native gas comparisons.',
        boundary='An LTE matter/photon split, not a physical opacity or radiation-transport solution. Do not assign the total effective conductive opacity to photon absorption/scattering, infer a Planck mean from a Rosseland mean, or count LTE photons twice. New GR geometry is still being reconstructed.'))


def build():
    module=next((model.CACHE/'build').rglob('mod_free_eos.mod')).parent
    cmd=['gfortran','-O2','-fPIC','-shared','-I'+str(module),str(OUT/'gas-bridge.f90'),
        '-L'+str(model.LIB.parent),'-Wl,-rpath,'+str(model.LIB.parent),'-l'+model.NAME,'-o',str(BRIDGE)]
    before=g.c.sha(model.LIB);r=subprocess.run(cmd,capture_output=True,text=True)
    (OUT/'build.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr
    assert g.c.sha(model.LIB)==before
    save('build.json',dict(classification='Counterexample candidate',command=cmd,completed=True,
        runtime={str(p):g.c.sha(p) for p in [BRIDGE,model.LIB]}))


class GasEOS(g.EOS):
    def __init__(self):
        super().__init__();self.gas_lib=ctypes.CDLL(str(BRIDGE));native=self.gas_lib.ionization_inventory
        array=np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS')
        native.argtypes=[ctypes.c_int,ctypes.c_double,ctypes.c_double,array,array,ctypes.POINTER(ctypes.c_int)]
        native.restype=None
        ctypes.c_int.in_dll(self.gas_lib,'__mod_free_eos_MOD_retain_molecules').value=1
        def capture(mode,value,t,eps,out,info):
            raw=np.full(24,np.nan);native(mode,value,t,eps,raw,info);out[:]=raw[:22];out[20]=0.
        self.call=capture
        values=np.zeros(2);fn=self.gas_lib.photon_constants;fn.argtypes=[array];fn.restype=None;fn(values)
        self.a_rad,self.c_light=map(float,values)


def split(a,T,arad):
    rho=a[...,0];prad=arad*T**4/3;urad=3*prad/rho
    P=a[...,1]-prad;cvT=a[...,10]-4*urad
    cr=a[...,1]*a[...,5]/P;ct=(a[...,1]*a[...,6]-4*prad)/P
    return dict(P_gas=P,u_gas=a[...,2]-urad,s_gas=a[...,3]-4*urad/(3*T),
        cvT_gas=cvT,du_dlnrho_gas=a[...,9]+urad,chi_rho_gas=cr,chi_T_gas=ct,
        gamma1_gas=cr+P/rho*ct**2/cvT,energy_radiation=3*prad,pressure_radiation=prad)


def symbolic():
    rho,T,a=sp.symbols('rho T a',positive=True);F=sp.Function('F')(rho,T);radiation=-a*T**4/(3*rho)
    pressure=rho*rho*sp.diff(radiation,rho);entropy=-sp.diff(radiation,T);energy=radiation+T*entropy
    assert sp.simplify(pressure-a*T**4/3)==0
    assert sp.simplify(energy-a*T**4/rho)==0 and sp.simplify(entropy-4*a*T**3/(3*rho))==0
    assert sp.simplify(T*sp.diff(energy,T)-4*energy)==0
    assert sp.simplify(rho*sp.diff(energy,rho)+energy)==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        identity='Subtracting the photon free energy and all its derivatives gives the stated matter-only P,u,s,cvT and density derivative. Reintroducing an independent radiation stress tensor must use these matter-only quantities; retaining the original total EOS would double count photons.',
        limitation='Conditional thermodynamic algebra. Native consistency and physical radiation closure are separate.'))


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    symbolic();gas=GasEOS();total=model.EOS();data=dict(np.load(reference.OUT/'reference-state.npz'));rows=[]
    for i in plan['cells']:
        for offset in plan['temperature_offsets']:
            r=float(data['lnd'][i]);t=float(data['lnT'][i]+offset);X=data['X'][i];T=np.exp(t)
            a=total(2,r,t,X);b=gas(2,r,t,X);again=total(2,r,t,X);expected=split(a,T,gas.a_rad)
            mapping=dict(P_gas=1,u_gas=2,s_gas=3,chi_rho_gas=5,chi_T_gas=6,du_dlnrho_gas=9,cvT_gas=10)
            errors={k:float(abs(b[j]-expected[k])/max(abs(expected[k]),1.)) for k,j in mapping.items()}
            budget=max(2.,32*np.spacing(abs(b[2]+b[1]/b[0])))
            energy_score=float(abs(b[2]-expected['u_gas'])/budget)
            entropy_score=float(T*abs(b[3]-expected['s_gas'])/budget)
            back=gas(1,float(np.log(b[1])),t,X);density_error=float(abs(np.log(back[0])-r))
            gamma_error=float(abs(b[4]/expected['gamma1_gas']-1))
            history=bool(np.array_equal(a,again));passed=bool(max(errors.values())<plan['decomposition_relative_tolerance']
                and energy_score<=1 and entropy_score<=1 and history
                and density_error<plan['inverse_log_density_tolerance'] and gamma_error<plan['gamma_identity_relative_tolerance'])
            rows.append(dict(cell=i,temperature_offset=offset,relative_errors=errors,energy_score=energy_score,
                entropy_score=entropy_score,pressure_inverse_log_density_error=density_error,
                gamma_identity_error=gamma_error,call_history_bitwise=history,passed=passed))
            np.savez_compressed(OUT/f'control-{i}-{offset}.npz',total=a,gas=b,gas_inverse=back,**expected)
    alltotal=np.concatenate([np.load(reference.OUT/f'block-{i:04}.npz')['new'] for i in range(0,len(data['X']),128)])
    values=split(alltotal,np.exp(data['lnT']),gas.a_rad)
    positive=bool(np.all(values['P_gas']>0)&np.all(values['cvT_gas']>0)&np.all(values['chi_rho_gas']>0))
    np.savez_compressed(OUT/'decomposed-reference.npz',**values,dm=data['dm'],X=data['X'],lnT=data['lnT'],lnd=data['lnd'])
    save('result.json',dict(classification='Counterexample candidate',completed=True,rows=rows,
        native_controls_passed=all(r['passed'] for r in rows),algebraic_reference_cells=len(alltotal),
        algebraic_gas_positivity_passed=positive,compiled_a_rad_cgs=gas.a_rad,compiled_light_speed_cgs=gas.c_light,
        maximum_pressure_radiation_fraction=float(np.max(values['pressure_radiation']/alltotal[:,1])),
        maximum_heat_capacity_radiation_fraction=float(np.max(4*values['energy_radiation']/alltotal[:,0]/alltotal[:,10])),
        full_grid_native_gas_comparison=False,physical_EOS_certified=False,photon_opacity_identified=False,
        radiation_transport_evolved=False,full_GR_evolution=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify()


def verify():
    for name,key in [('plan.json','bindings'),('manifest.json','sha256')]:
        for rel,digest in json.loads((OUT/name).read_text())[key].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for path,digest in json.loads((OUT/'build.json').read_text())['runtime'].items():assert g.c.sha(path)==digest,path
    assert json.loads((OUT/'result.json').read_text())['completed']
    print('PASS matter/photon EOS split bindings; inspect native comparison and positivity gates',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
