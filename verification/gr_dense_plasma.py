"""Fixed-inventory dense-plasma free energies and an independent provider audit.

Counterexample candidate. A model comparison, not an EOS replacement or
physical-error certificate. Keep the distributed Fortran source unchanged.
"""
import ctypes, json, shutil, subprocess, sys
import mpmath as mp
import numpy as np
import sympy as sp
import gr_fermi_plasma as plasma
import eos_species_inventory as inventory

g=plasma.g;model=plasma.split.model;OUT=g.OUT/'gr-dense-plasma'
CACHE=g.CACHE/'dense-plasma';FIELDS=['F','U','P','S','CV','PDT','PDR']


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists() and not CACHE.exists();OUT.mkdir();CACHE.mkdir();plasma.verify()
    names=['potekhin-chabrier-eos22.f','potekhin-chabrier-eipsubr.html',
        'potekhin-chabrier2010.pdf','potekhin-chabrier2010.txt','baiko-chugunov2022.pdf','baiko-chugunov2022.txt']
    for name in names:shutil.copy2(g.ROOT/'outputs'/name,OUT/name)
    for name in ['mod_coulomb.f90','mod_exchange.f90','master_exchange.f90']:
        shutil.copy2(model.CACHE/'source/src'/name,OUT/name)
    (OUT/'cached-free.f90').write_text('''subroutine cached_free(res) bind(C)
  use iso_c_binding
  use mod_free_eos_constants, only: boltzmann,cpe
  use mod_master_coulomb_data, only: fcoulomb
  use mod_master_exchange_data, only: iforder,fexprime,n_e,t,p_e,pstarprime,psiprime,psi
  implicit none
  real(c_double),intent(out) :: res(5)
  real(c_double) :: fex2
  fex2=boltzmann*t*(psiprime-psi)*n_e-(cpe*pstarprime(1)-p_e)
  res=[fcoulomb,fexprime+fex2,n_e,t,real(iforder,c_double)]
end subroutine cached_free
''')
    paths=[g.ROOT/'verification/gr_dense_plasma.py',g.ROOT/'verification/eos_species_inventory.py',
        plasma.OUT/'manifest.json',plasma.split.reference.OUT/'manifest.json',model.OUT/'manifest.json']+list(OUT.iterdir())
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='e9500e5',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        native_library={str(model.LIB):g.c.sha(model.LIB)},
        sources=['https://www.ioffe.ru/astro/EIP/eos22.f','https://www.ioffe.ru/astro/EIP/eipsubr.html',
            'https://arxiv.org/abs/1001.0690','https://arxiv.org/abs/2112.04822'],
        provider_revision='09.06.22; nominal ionic mass A_i*m_u, internal AUM=1822.88848 retained.',
        target='Call classical EOSFI22 excess outputs and CORMIX directly at fixed baryon/nuclear inventory. Add the quantum free energy of BC2022 eq34 once, differentiating that same free energy. Audit the distributed EOSFI22 and LIQUBC outputs independently before using them.',
        domain='Use all positive isotope number fractions X_i/A_i, with no trace cutoff; neutrality fixes n_e=N_A*rho_B*sum(Z_i X_i/A_i). A_i is the stored baryon mass number, not a newly measured nuclear mass. Native constants convert physical density and temperature to r_s and Gamma_e. No MELANGE density inversion, ideal-electron addition, photons, spin or ideal-ion free energy.',
        coverage='Audit native 5735 states first. For the comparison require abs(1-ne_native/ne_fully_ionized)<=1e-8. This is a charge-deficit selection, not a proof that every trace ion is fully stripped. Do not apply this fully ionized EOS to the excluded partially ionized states.',
        ionization_selection_tolerance=1e-8,inventory_tolerance=1e-10,cache_density_tolerance=1e-10,
        controls=[0,1175,2972,3043,4352,5734],
        quantum_controls_R=[90.,500.,1900.,120000.],quantum_controls_theta=[1e-5,.002,.1,1.,10.,100.],
        quantum_mp_score_tolerance=1e-8,quantum_mp_score_floor=1e-10,
        finite_log_steps=[2e-4,1e-4],finite_derivative_absolute_normalized_tolerance=1e-5,
        original_provider_failures='Retain raw first/second EOSFI22 output groups, LIQUBC and all failed finite derivative gates. Source inspection indicates first excess group omits quantum additions and liquid PDR2 uses PDTQL. No retrospective threshold change.',
        scope='Finite implementation and source consistency. Fitted plasma correlations, mixture quantum physics, partial ionization, physical and continuous error, EOS/GR replacement, transport, scalar and observation closure are not certified.'))


def bindings():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for path,digest in plan['native_library'].items():assert g.c.sha(path)==digest,path
    return plan


def build():
    bindings();commands=[['gfortran','-O2','-fPIC','-shared','-std=legacy',
        str(OUT/'potekhin-chabrier-eos22.f'),'-o',str(CACHE/'pc.so')],
        ['gfortran','-O2','-fPIC','-shared','-I'+str(model.LIB.parent),str(OUT/'cached-free.f90'),
         '-L'+str(model.LIB.parent),'-Wl,-rpath,'+str(model.LIB.parent),'-l'+model.NAME,'-o',str(CACHE/'cached.so')]]
    logs=[]
    for command in commands:
        result=subprocess.run(command,capture_output=True,text=True);logs.append(dict(command=command,returncode=result.returncode,stdout=result.stdout,stderr=result.stderr))
        save('build-log.json',logs);assert result.returncode==0,result.stderr
    paths=[CACHE/'pc.so',CACHE/'cached.so',model.LIB];linked=[]
    for path in paths:
        result=subprocess.run(['ldd',str(path)],capture_output=True,text=True);assert result.returncode==0
        linked.append(dict(path=str(path),output=result.stdout))
    dependencies={}
    for entry in linked:
        for line in entry['output'].splitlines():
            for token in line.split():
                if token.startswith('/') and g.Path(token).is_file():dependencies[token]=g.c.sha(token)
    save('runtime.json',dict(classification='Imported from prior work',compiler=subprocess.check_output(['gfortran','--version'],text=True),
        sha256={**{str(p):g.c.sha(p) for p in paths},**dependencies},ldd=linked))


class Provider:
    def __init__(self):self.lib=ctypes.CDLL(str(CACHE/'pc.so'))

    def call(self,name,inputs,nout):
        args=[ctypes.c_int(v) if type(v) is int else ctypes.c_double(v) for v in inputs]
        outputs=[ctypes.c_double(float('nan')) for _ in range(nout)];args+=outputs
        fn=getattr(self.lib,name+'_');fn.restype=None;fn.argtypes=[ctypes.POINTER(type(x)) for x in args]
        fn(*[ctypes.byref(x) for x in args]);values=np.array([x.value for x in outputs])
        assert np.all(np.isfinite(values)),(name,inputs,values)
        return values


def seven(a):return np.r_[a[:3],a[1]-a[0],a[3:]]


def quantum(R,theta):
    """BC2022 eq34, with all thermodynamics derived from that free energy."""
    assert R>0 and theta>0
    c1=.351*R/(90+R);c=np.array([c1,.294,np.sqrt(1-c1*c1-.294**2)])
    d1=90/(90+R);d=np.array([d1,0.,-(c1/c[2])**2*d1])
    d1p=-R/(90+R)*d1;dp=np.array([d1p,0.,-(c1/c[2])**2*(2*(d1-d[2])*d1+d1p)])
    y=c*theta;f=np.empty(3);u=np.empty(3);cv=np.empty(3);small=y<.1;z=y[small]**2
    # Same analytic series for f, y*f_y and u-y*u_y, through y^12.
    coefficients=np.array([1/24,-1/2880,1/181440,-1/9676800,1/479001600,-691/15692092416000])
    powers=2*np.arange(1,7)
    f[small]=z*np.polynomial.polynomial.polyval(z,coefficients)
    u[small]=z*np.polynomial.polynomial.polyval(z,coefficients*powers)
    cv[small]=z*np.polynomial.polynomial.polyval(z,coefficients*powers*(1-powers))
    yy=y[~small];ex=np.exp(-yy);den=-np.expm1(-yy)
    f[~small]=np.log(den/yy)+yy/2
    u[~small]=yy/den-yy/2-1
    cv[~small]=yy*yy*ex/(den*den)-1
    a=.5-d/3;p=a*u;pdt=a*cv;pdr=(a+a*a+dp/9)*u-a*a*cv
    return np.array([f.sum(),u.sum(),p.sum(),(u-f).sum(),cv.sum(),pdt.sum(),pdr.sum()])


def independent_quantum(R,theta):
    mp.mp.dps=60;R=mp.mpf(str(R));theta=mp.mpf(str(theta))
    def f(s,t):
        r=R*mp.exp(-s/3);th=theta*mp.exp(s/2-t)
        c1=mp.mpf('.351')*r/(90+r);c2=mp.mpf('.294');c3=mp.sqrt(1-c1*c1-c2*c2)
        return sum(mp.log(2*mp.sinh(c*th/2)/(c*th)) for c in [c1,c2,c3])
    F=f(0,0);U=-mp.diff(f,(0,0),(0,1));P=mp.diff(f,(0,0),(1,0))
    return np.array(list(map(float,[F,U,P,U-F,U-mp.diff(f,(0,0),(0,2)),
        P+mp.diff(f,(0,0),(1,1)),P+mp.diff(f,(0,0),(2,0))])))


def symbolic():
    r=sp.symbols('r',positive=True);c1=sp.Rational(351,1000)*r/(90+r);c2=sp.Rational(294,1000);c3=sp.sqrt(1-c1*c1-c2*c2)
    D=[r*sp.diff(sp.log(c),r) for c in [c1,c2,c3]];d1=90/(90+r)
    assert sp.simplify(D[0]-d1)==0 and D[1]==0
    assert sp.simplify(D[2]+(c1/c3)**2*d1)==0
    assert sp.simplify(r*sp.diff(D[2],r)+(c1/c3)**2*(2*(d1-D[2])*d1-r/(90+r)*d1))==0
    a,d,dp,u,cv=sp.symbols('a d dp u cv');a=sp.Rational(1,2)-d/3
    derived=(a+a*a+dp/9)*u-a*a*cv
    printed=sp.Rational(3,4)*u-cv/4-d*(u-cv/3)/2+dp*u/9
    assert sp.factor(derived-printed)==sp.factor(d*(2*d-3)*(u-cv)/18)
    y=sp.symbols('y');series=sp.series(sp.log(sp.sinh(y/2)/(y/2)),y,0,14).removeO()
    expected=sum(v*y**(2*j) for j,v in enumerate([sp.Rational(1,24),-sp.Rational(1,2880),sp.Rational(1,181440),
        -sp.Rational(1,9676800),sp.Rational(1,479001600),-sp.Rational(691,15692092416000)],1))
    assert sp.expand(series-expected)==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        free='f_q=sum log(2*sinh(C_i(R)*theta/2)/(C_i(R)*theta)); R proportional n^(-1/3), theta proportional n^(1/2)/T at fixed species.',
        derivatives='D_i=d ln C_i/d ln R; a_i=1/2-D_i/3; u_i=y_i*f_i_prime, cv_i=u_i-y_i*u_i_prime; p_i=a_i*u_i; PDT_i=a_i*cv_i; PDR_i=(a_i+a_i^2+D_i_prime/9)*u_i-a_i^2*cv_i, D_i_prime=dD_i/dlnR.',
        printed_equation39='The derivative of eq34/36 differs term by term from printed eq39 by D_i*(2D_i-3)*(u_i-cv_i)/18 even if the exact D_i_prime is used. This need not sum to zero. The distributed source also substitutes -R/(90+R)*D_i for every D_i_prime, which is not the exact D_3_prime.',
        scope='Conditional algebra for the declared BC2022 free-energy fit. No new many-body calculation or physical uncertainty bound.'))


def controls():
    plan=bindings();p=Provider();symbolic();rows=[]
    for R in plan['quantum_controls_R']:
        for theta in plan['quantum_controls_theta']:
            ref=independent_quantum(R,theta);new=quantum(R,theta);raw=seven(p.call('liqubc',[R,theta],6))
            scale=np.maximum(abs(ref),plan['quantum_mp_score_floor'])
            score=float(np.max(abs(new-ref)/scale));rawscore=float(np.max(abs(raw-ref)/scale))
            rows.append(dict(R=R,theta=theta,independent=ref.tolist(),derived=new.tolist(),provider=raw.tolist(),
                derived_score=score,provider_score=rawscore,derived_passed=score<plan['quantum_mp_score_tolerance'],
                provider_passed=rawscore<plan['quantum_mp_score_tolerance']))
    save('quantum-controls.json',dict(classification='Counterexample candidate',rows=rows,
        derived_passed=all(r['derived_passed'] for r in rows),provider_passed=all(r['provider_passed'] for r in rows)))
    assert all(r['derived_passed'] for r in rows)
    print('QUANTUM',len(rows),'independent controls; derived max',max(r['derived_score'] for r in rows),
        'provider max',max(r['provider_score'] for r in rows),flush=True)


def state_data():
    state=dict(np.load(plasma.split.reference.OUT/'reference-state.npz'))
    native=np.concatenate([np.load(plasma.split.reference.OUT/f'block-{i:04}.npz')['new'] for i in range(0,5735,128)])
    return state,native


def mixture(p,r,t,X):
    c=plasma.constants();count=X/g.c.A;number=count.sum();x=count/number;z=g.c.Z
    ne=c['N_A']*np.exp(r)*float(count@z);ae=(3/(4*np.pi*ne))**(1/3)
    a0=c['hbar']**2/(c['m_e_g']*c['e_esu']**2);rs=ae/a0;ge=c['e_esu']**2/(ae*c['k_B']*np.exp(t))
    moments=[float(x@z),float(x@(z*z)),float(x@(z**2.5)),float(x@(z**(5/3))),float(x@(z*(z+1)**1.5))]
    assert ge*moments[3]<175,'Liquid phase is required'
    original=np.zeros((2,7));ql=np.zeros(7);rawql=np.zeros(7);component_parameters=[]
    for j in np.flatnonzero(x>0):
        gamma=ge*z[j]**(5/3);theta=gamma/np.sqrt(rs)*np.sqrt(3/(1822.88848*g.c.A[j]))/z[j]**(7/6)
        R=3*(gamma/theta)**2
        original+=x[j]*p.call('eosfi22',[0,float(g.c.A[j]),float(z[j]),rs,gamma],14).reshape(2,7)
        ql+=x[j]*quantum(R,theta);rawql+=x[j]*seven(p.call('liqubc',[R,theta],6))
        component_parameters.append([int(j),float(x[j]),float(R),float(theta),float(gamma)])
    mix=seven(p.call('cormix',[rs,ge,*moments],6));classical=original[0]+mix
    return dict(corrected=classical+ql,classical=classical,quantum=ql,raw_quantum=rawql,raw_groups=original,
        mixing=mix,rs=rs,Gamma_e=ge,Gamma_mean=ge*moments[3],ni_kT=c['N_A']*np.exp(r)*number*c['k_B']*np.exp(t),
        ne_full=ne,components=component_parameters)


def run():
    plan=bindings();assert json.loads((OUT/'quantum-controls.json').read_text())['derived_passed']
    state,native=state_data();c=plasma.constants();ne=native[:,13]*c['N_A']
    full=c['N_A']*np.exp(state['lnd'])*((state['X']/g.c.A)@g.c.Z);deficit=1-ne/full
    mask=abs(deficit)<=plan['ionization_selection_tolerance'];indices=np.flatnonzero(mask)
    save('selection.json',dict(classification='Counterexample candidate',total_cells=len(mask),selected_cells=len(indices),
        selected_baryon_mass_fraction=float(state['dm'][mask].sum()/state['dm'].sum()),
        deficit_range=[float(deficit.min()),float(deficit.max())],
        rule='Absolute aggregate charge deficit <=1e-8; near complete ionization, not an exact trace-species statement.'))
    p=Provider();eos=model.EOS();lib=ctypes.CDLL(str(CACHE/'cached.so'));fn=lib.cached_free
    fn.argtypes=[np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS')];fn.restype=None
    # Reuse the existing deterministic native snapshot and inventory checker.
    controls=[]
    for i in plan['controls']:
        snap=eos.snapshot(state['lnd'][i],state['lnT'][i],state['X'][i]);cached=np.empty(5);fn(cached)
        check=inventory.check(snap,state['X'][i],state['lnd'][i],eos)
        same=np.array_equal(snap['eos'],native[i]);density=abs(cached[2]/ne[i]-1)
        assert not check['missing_nonzero_elements'] and check['inventory_error']<plan['inventory_tolerance'] and check['charge_error']<plan['inventory_tolerance']
        assert same and density<plan['cache_density_tolerance'] and cached[4]==2
        controls.append(dict(cell=i,**check,native_outputs_bitwise=same,cache_ne_relative_error=float(density),cached=cached.tolist()))
    save('native-controls.json',dict(classification='Counterexample candidate',passed=True,rows=controls))
    vals=[];caches=[];quantums=[];mixings=[];params=[];rows=[];originals=[]
    for i in indices:
        r,t,X=state['lnd'][i],state['lnT'][i],state['X'][i];a=eos(2,r,t,X);cached=np.empty(5);fn(cached)
        assert np.array_equal(a,native[i]),int(i)
        assert cached[4]==2 and abs(cached[2]/ne[i]-1)<plan['cache_density_tolerance']
        value=mixture(p,r,t,X);vals.append(value['corrected']);caches.append(cached)
        quantums.append(value['quantum']);mixings.append(value['mixing']);originals.append(value['raw_groups'])
        params.append([value['rs'],value['Gamma_e'],value['Gamma_mean'],value['ni_kT'],value['ne_full']])
        row=dict(cell=int(i),raw_group_F_difference=float(value['raw_groups'][1,0]-value['raw_groups'][0,0]),
            quantum_free=float(value['quantum'][0]))
        if i in plan['controls']:
            finite=[]
            for h in plan['finite_log_steps']:
                rminus=mixture(p,r-h,t,X)['corrected'];rplus=mixture(p,r+h,t,X)['corrected']
                tminus=mixture(p,r,t-h,X)['corrected'];tplus=mixture(p,r,t+h,X)['corrected'];b=value['corrected']
                estim=np.array([-(tplus[0]-tminus[0])/(2*h),(rplus[0]-rminus[0])/(2*h),
                    b[1]+(tplus[1]-tminus[1])/(2*h),b[2]+(tplus[2]-tminus[2])/(2*h),b[2]+(rplus[2]-rminus[2])/(2*h)])
                expected=b[[1,2,4,5,6]];score=float(np.max(abs(estim-expected)/np.maximum(1,abs(expected))))
                finite.append(dict(step=h,estimated=estim.tolist(),expected=expected.tolist(),score=score,
                    passed=score<plan['finite_derivative_absolute_normalized_tolerance']))
            row['finite_derivatives']=finite
        rows.append(row)
        if len(rows)%128==0: print('DENSE PLASMA',len(rows),'/',len(indices),flush=True)
    vals=np.array(vals);caches=np.array(caches);params=np.array(params);native_free=caches[:,:2].sum(1)/params[:,3]
    residual=vals[:,0]-native_free
    np.savez_compressed(OUT/'stellar-comparison.npz',cells=indices,deficit_all=deficit,selected_mask=mask,
        corrected=vals,native_cached=caches,quantum=np.array(quantums),mixing=np.array(mixings),raw_groups=np.array(originals),
        parameters=params,native_reduced_free=native_free,residual_reduced_free=residual)
    finite=[a for r in rows for a in r.get('finite_derivatives',[])];save('rows.json',rows)
    save('result.json',dict(classification='Counterexample candidate',completed=True,cells=len(indices),
        all_native_outputs_bitwise=True,finite_derivatives_passed=all(a['passed'] for a in finite),
        maximum_finite_derivative_score=max(a['score'] for a in finite),
        residual_free_range=[float(residual.min()),float(residual.max())],
        maximum_residual_over_native=float(np.max(abs(residual)/abs(native_free))),
        quantum_free_range=[float(x) for x in [np.min(quantums,axis=0)[0],np.max(quantums,axis=0)[0]]],
        gamma_mean_range=[float(params[:,2].min()),float(params[:,2].max())],
        interpretations='PC full-ion free energy minus cached native Coulomb+exchange at the saved near-fully-ionized states. Different equilibrium/ionization and provider fit conventions remain; this difference is not a physical error bound or a time-dependent released energy.',
        physical_EOS_certified=False,native_EOS_replaced=False,full_GR_evolution=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for path,digest in json.loads((OUT/'runtime.json').read_text())['sha256'].items():assert g.c.sha(path)==digest,path
    assert json.loads((OUT/'result.json').read_text())['completed']
    print('PASS dense plasma provenance; raw provider failures and separate finite/physical verdicts retained',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
