"""Actual MDH pressure-ionization moments and constrained component curvature."""
from concurrent.futures import ProcessPoolExecutor, as_completed
import ctypes, json, shutil, subprocess, sys
import numpy as np
import sympy as sp
import exchange_curvature as ex

c=ex.c;g=c.g;OUT=g.OUT/'pressure-ionization';CACHE=g.CACHE/'pressure-ionization'
NAME='free_eos_direct24_pi_trace';LIB=CACHE/'build/src'/('lib'+NAME+'.so.1.0.0')
L=1e-8


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')
def read(path): return json.loads(path.read_text())


def polynomials():
    a=sp.symbols('a0:7');b,y=sp.symbols('b y');u=(*a,b,y)
    p2=2*(3*a[1]*a[2]+a[0]*a[3])+b*y
    p3=3*a[0]*a[3]**2+30*a[1]*a[2]*a[3]+9*a[2]**3+6*a[0]*a[2]*a[4]+9*a[1]**2*a[4]+6*a[0]*a[1]*a[5]+a[0]**2*a[6]
    return u,p2,p3


def functions():
    u,p2,p3=polynomials()
    return [sp.lambdify(u,v,'numpy') for v in [p2,p3,sp.Matrix([sp.diff(p2,x) for x in u]),
        sp.Matrix([sp.diff(p3,x) for x in u]),sp.hessian(p2,u),sp.hessian(p3,u)]]


def prepare():
    assert not OUT.exists() and not CACHE.exists();OUT.mkdir();CACHE.mkdir()
    original=g.d.CACHE/'full-integral-source';source=CACHE/'source';shutil.copytree(original,source);changed={}
    patches={
        'mod_pi.f90':[('  public mdh_pi_called, quad','''  public mdh_pi_called, quad, pi_trace_state, pi_trace_extra
  real(fp_kind), save :: pi_trace_state(12)=0._fp_kind, pi_trace_extra(9,3)=0._fp_kind''')],
        'mdh_pi.f90':[('end subroutine mdh_pi_end','''  pi_trace_extra(:,1)=extrasum
  pi_trace_extra(:,2)=extrasumf
  pi_trace_extra(:,3)=extrasumt
  pi_trace_state=[t,rho,rf,rt,ppi,ppif,ppit,spi,spif,spit,upi,quad]
end subroutine mdh_pi_end''')],
        'mod_pi_fit.f90':[('  implicit none','''  implicit none
  real(fp_kind), save, public :: pi_trace_neutral(26)=0._fp_kind, pi_trace_ion3(318)=0._fp_kind''')],
        'effective_radius.f90':[('end subroutine effective_radius','''  pi_trace_neutral=r_neutral
  pi_trace_ion3=r_ion3
end subroutine effective_radius''')],
        'CMakeLists.txt':[('OUTPUT_NAME free_eos_direct24_integral_full','OUTPUT_NAME '+NAME)]}
    # mdh_pi_end has its own ONLY import; change only that subroutine.
    path=source/'src/mdh_pi.f90';raw=path.read_text();head,tail=raw.split('subroutine mdh_pi_end(',1)
    tail=tail.replace('  use mod_mdh_pi_data, only: mdh_pi_called, quad',
        '  use mod_mdh_pi_data, only: mdh_pi_called, quad, pi_trace_state, pi_trace_extra',1)
    path.write_text(head+'subroutine mdh_pi_end('+tail)
    for name,pairs in patches.items():
        path=source/'src'/name;text=path.read_text()
        for old,new in pairs:
            assert text.count(old)==1,(name,old,text.count(old));text=text.replace(old,new)
        path.write_text(text);shutil.copy2(original/'src'/name,OUT/('before-'+name));shutil.copy2(path,OUT/name)
        changed[name]=dict(before=g.c.sha(original/'src'/name),after=g.c.sha(path))
    for name in ['mod_free_eos_constants.f90','free_eos_detailed.f90']:
        shutil.copy2(original/'src'/name,OUT/name)
    shutil.copy2(c.s.OUT/'inventory_bridge.f90',OUT/'direct_ion_bridge.f90')
    paths=[g.ROOT/'verification/pressure_ionization_curvature.py',c.OUT/'manifest.json',ex.OUT/'manifest.json',
        g.OUT/'reference-state.npz',OUT/'direct_ion_bridge.f90']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='0a5ffca',
        cells=5735,processes=4,block_cells=128,control_cells=[0,2972,5734],source_changes=changed,
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},original_library_sha256=g.c.sha(c.s.LIB),
        original_EOS_Coulomb_bitwise_required=True,moment_relative_tolerance=1e-10,
        thermodynamic_normalized_tolerance=1e-9,projected_Coulomb_Gram_tolerance=1e-12,
        dimensions='Three electron/Coulomb coordinates followed by seven neutral-radius moments alpha0..alpha6, ionic charge^(3/2) beta and non-bare radius-cube gamma. H2+ contributes to the hard-sphere neutral moments as in the actual MDH convention.',
        moment_units='Capture native extrasum in nu=n/(rho_native*N_A) form and actual radii. Scale radii by L=1e-8 cm for matrix conditioning only. Nuclear inventory constraints count two H nuclei per molecule.',
        method='An isolated read-only trace captures final mdh_pi_end inputs/outputs and effective radii. Reconstruct moments from saved species, replay pressure and entropy including their supplied chain-rule derivatives, then differentiate the same explicit quadratic/cubic free-energy polynomial.',
        boundary='Positive saved support and declared finite matrices only. Excitation/partition terms, omitted populations, coupled global roots, native continuous evaluation and physical EOS errors are not certified.',
        physical_EOS_certified=False,full_GR_evolution=False))


def build():
    old=g.d.OUT
    try:
        g.d.OUT=OUT;g.d.build_at(CACHE/'source',CACHE/'build',NAME,CACHE/'unused-inventory.so')
    finally: g.d.OUT=old
    command=['gfortran','-O2','-fPIC','-shared','-I'+str(LIB.parent),str(c.OUT/'coulomb_bridge.f90'),
        '-L'+str(LIB.parent),'-Wl,-rpath,'+str(LIB.parent),'-l'+NAME,'-o',str(CACHE/'pi.so')]
    result=subprocess.run(command,capture_output=True,text=True);(OUT/'pi-bridge-build.log').write_text(result.stdout+result.stderr)
    assert result.returncode==0,result.stderr
    save('bridge.json',dict(classification='Counterexample candidate',command=command,bridge_sha256=g.c.sha(CACHE/'pi.so'),library_sha256=g.c.sha(LIB)))
    assert g.c.sha(c.s.LIB)==read(OUT/'plan.json')['original_library_sha256']


class EOS(c.EOS):
    def __init__(self):
        super().__init__();self.inventory_lib=ctypes.CDLL(str(CACHE/'pi.so'))
        array=np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS');call=self.inventory_lib.coulomb_inventory
        call.argtypes=[ctypes.c_int,ctypes.c_double,ctypes.c_double,array,array,ctypes.POINTER(ctypes.c_int)];call.restype=None
        def capture(mode,value,t,eps,out,info):
            raw=np.full(25,np.nan);call(mode,value,t,eps,raw,info)
            out[:]=raw[:22];out[20]=0.;self.molecules=raw[22:24].copy();self.fl=float(raw[24]);self.rho_native=float(raw[0])
        self.call=capture;self.probe_call=self.inventory_lib.coulomb_probe
        self.probe_call.argtypes=[ctypes.c_double]*4+[array];self.probe_call.restype=None
        self.poly=functions()

    def array(self,module,name,size):
        return np.ctypeslib.as_array((ctypes.c_double*size).in_dll(self.inventory_lib,'__'+module+'_MOD_'+name)).copy()

    def full(self,r,t,x,electron):
        snap,a,_,arg,pops,mols,row=self.actual(r,t,x)
        state=self.array('mod_mdh_pi_data','pi_trace_state',12)
        extra=self.array('mod_mdh_pi_data','pi_trace_extra',27).reshape((3,9))
        neutral=self.array('mod_pi_fit','pi_trace_neutral',26)
        ion3=self.array('mod_pi_fit','pi_trace_ion3',318)
        assert state[0]>0 and state[1]>0 and state[11]==10 and state[10]==0
        assert np.all(neutral>0) and np.all(ion3>0) and np.all(np.isfinite(extra))
        ym=(x/g.c.A)@self.mapping;eps=ym/float(ym@self.weights)
        molecular=eps[0]*snap['molecular_H_fractions']/2
        n,B,atoms,elements=species_coordinates(snap['number_fractions'],molecular,neutral,ion3)
        units=np.r_[L**np.arange(7),1.,L**3]
        u,uf,ut=extra/units;reconstructed=B[3:]@n
        moment_error=float(np.max(abs(reconstructed-u)/np.maximum(abs(u),1e-300)))
        S=state[1]*self.constants[0];kT=self.constants[1]*state[0];kappa=S*(4*np.pi/3)*L**3;quad=state[11]
        p2,p3,grad2,grad3,h2,h3=[np.asarray(fun(*u),dtype=float) for fun in self.poly]
        grad2=grad2.ravel();grad3=grad3.ravel()
        f=kappa*p2+quad*kappa**2*p3;pressure_reduced=kappa*p2+2*quad*kappa**2*p3
        free=kT*S*f;pressure=kT*S*pressure_reduced;entropy=-self.constants[1]*self.constants[0]*f
        grad_f=kappa*grad2+quad*kappa**2*grad3;grad_p=kappa*grad2+2*quad*kappa**2*grad3
        density_p=2*kappa*p2+6*quad*kappa**2*p3
        sf=-self.constants[1]*self.constants[0]*(grad_f@uf+pressure_reduced*state[2])
        st=-self.constants[1]*self.constants[0]*(grad_f@ut+pressure_reduced*state[3])
        pf=kT*S*(grad_p@uf+density_p*state[2]);pt=pressure+kT*S*(grad_p@ut+density_p*state[3])
        native=np.array([-state[7]*state[0]*state[1],state[4],state[7],state[8],state[9],state[5],state[6]])
        computed=np.array([free,pressure,entropy,sf,st,pf,pt])
        normal=np.array([abs(native[0]),abs(native[1]),abs(native[2]),max(abs(native[2]),abs(native[3])),
            max(abs(native[2]),abs(native[4])),max(abs(native[1]),abs(native[5])),max(abs(native[1]),abs(native[6]))])
        thermo=float(np.max(abs(native-computed)/np.maximum(normal,1e-300)))
        G=gram(n,B,atoms,elements);Gprevious=c.tangent_gram(pops,mols)
        gram_error=float(abs(S*G[:3,:3]-Gprevious).max()/max(1.,float(np.sum(pops)+np.sum(mols))))
        H=np.zeros((12,12));H[:3,:3]=S*c.matrix(a)/a[15]
        H[0,0]+=S*(a[13]+electron[2])/(a[10]*a[12])
        H[3:,3:]=kappa*h2+quad*kappa**2*h3
        eigen,Q=np.linalg.eigh(G);assert eigen.min()>=-1e-12*max(1.,eigen.max())
        root=(Q*np.sqrt(np.maximum(eigen,0)))@Q.T;C=root@H@root;C=(C+C.T)/2
        margin=float(1+min(0.,np.linalg.eigvalsh(C)[0]))
        row.update(moment_relative_error=moment_error,thermodynamic_normalized_error=thermo,
            projected_Coulomb_Gram_error=gram_error,combined_component_margin=margin)
        return snap,a,state,extra,neutral,ion3,G,H,computed,native,row


def species_coordinates(populations,molecules,neutral,ion3):
    n=[];columns=[];atoms=[];elements=[];offset=0
    for element,Z in enumerate(g.d.CHARGES):
        for j in range(Z+1):
            if populations[element,j]==0: continue
            alpha=(neutral[element]/L)**np.arange(7) if j==0 else np.zeros(7)
            b=[j,int(j>0),j*j,*alpha,j**1.5,ion3[offset+j]/L**3 if j<Z else 0.]
            n.append(populations[element,j]);columns.append(b);atoms.append(1);elements.append(element)
        offset+=Z
    assert offset==316
    for j,v in enumerate(molecules):
        if v==0: continue
        n.append(v);columns.append([j,j,j,*((neutral[24+j]/L)**np.arange(7)),j,ion3[316+j]/L**3]);atoms.append(2);elements.append(0)
    return np.array(n),np.array(columns).T,np.array(atoms),np.array(elements)


def gram(n,B,atoms,elements):
    n=n.astype(np.longdouble);B=B.astype(np.longdouble);atoms=atoms.astype(np.longdouble);G=np.zeros((12,12),dtype=np.longdouble)
    for element in np.unique(elements):
        selected=elements==element;v=n[selected];A=atoms[selected];b=B[:,selected]
        center=np.sum(b*(v*A),axis=1)/np.sum(v*A*A);residual=b-center[:,None]*A
        G+=(residual*v)@residual.T
    return np.asarray(G,dtype=float)


def symbolic():
    u,p2,p3=polynomials()
    for degree,p in [(2,p2),(3,p3)]:
        grad=sp.Matrix([sp.diff(p,x) for x in u]);H=sp.hessian(p,u)
        assert sp.expand(sum(x*sp.diff(p,x) for x in u)-degree*p)==0
        assert (H*sp.Matrix(u)-(degree-1)*grad).applyfunc(sp.expand)==sp.zeros(9,1)
    save('symbolic.json',dict(classification='Proven',passed=True,
        quadratic=str(p2),cubic=str(p3),
        scope='Explicit MDH ground-state pressure-ionization polynomial. All moments are linear in species densities at fixed radii; Hessians pull back by B^T H B. The quadratic and cubic homogeneity identities give pressure F2+2F3 and zero internal energy for terms linear in T at fixed species density. Excitation/partition free energies are separate.',
        entropy='The entropy and pressure derivatives include supplied moment derivatives and the native logarithmic density derivatives rf,rt. A fixed-moment derivative is not substituted for the full EOS chain rule.'))


def gates(row,plan):
    return row['moment_relative_error']<plan['moment_relative_tolerance'] and row['thermodynamic_normalized_error']<plan['thermodynamic_normalized_tolerance'] and row['projected_Coulomb_Gram_error']<plan['projected_Coulomb_Gram_tolerance']


def control():
    plan=read(OUT/'plan.json');state=dict(np.load(g.OUT/'reference-state.npz'));eos=EOS();rows=[];symbolic()
    for i in plan['control_cells']:
        electron=dict(np.load(ex.OUT/f'control-{i}.npz'))['values'];data=eos.full(state['lnd'][i],state['lnT'][i],state['X'][i],electron)
        snap,a,meta,extra,rn,ri,G,H,computed,native,row=data
        old=dict(np.load(c.OUT/f'control-{i}.npz'));assert np.array_equal(a,old['values'])
        species=dict(np.load(c.s.OUT/f'control-{i}.npz'));assert all(np.array_equal(snap[k],v) for k,v in species.items())
        row.update(cell=i,passed=gates(row,plan));rows.append(row)
        np.savez_compressed(OUT/f'control-{i}.npz',state=meta,extra=extra,neutral=rn,ion3=ri,Gram=G,Hessian=H,computed=computed,native=native)
        print('PRESSURE IONIZATION CONTROL',row,flush=True)
    save('control.json',dict(classification='Counterexample candidate',passed=all(r['passed'] for r in rows),rows=rows));assert all(r['passed'] for r in rows)


def block(start):
    plan=read(OUT/'plan.json');state=dict(np.load(g.OUT/'reference-state.npz'));eos=EOS();rows=[];arrays=[]
    electron=dict(np.load(ex.OUT/f'block-{start}.npz'));old=dict(np.load(c.OUT/f'block-{start}.npz'));stop=min(start+128,plan['cells'])
    for i in range(start,stop):
        snap,a,meta,extra,rn,ri,G,H,computed,native,row=eos.full(state['lnd'][i],state['lnT'][i],state['X'][i],electron['values'][i-start])
        assert np.array_equal(a,old['values'][i-start]) and np.array_equal(snap['eos'],old['eos'][i-start])
        row.update(cell=i,passed=gates(row,plan));rows.append(row)
        arrays.append(dict(state=meta,extra=extra,neutral=rn,ion3=ri,Gram=G,Hessian=H,computed=computed,native=native))
    path=OUT/f'block-{start}.npz';np.savez_compressed(path,**{k:np.array([a[k] for a in arrays]) for k in arrays[0]})
    record=dict(classification='Counterexample candidate',start=start,stop=stop,rows=rows,passed=all(r['passed'] for r in rows),
        output_sha256=g.c.sha(path),plan_sha256=g.c.sha(OUT/'plan.json'))
    save(f'block-{start}.json',record);assert record['passed'],start;return record


def run():
    plan=read(OUT/'plan.json');assert read(OUT/'control.json')['passed']
    for rel,digest in plan['bindings'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    bridge=read(OUT/'bridge.json');assert g.c.sha(LIB)==bridge['library_sha256'] and g.c.sha(CACHE/'pi.so')==bridge['bridge_sha256']
    assert g.c.sha(c.s.LIB)==plan['original_library_sha256'];rows=[]
    with ProcessPoolExecutor(max_workers=plan['processes']) as pool:
        for done in as_completed([pool.submit(block,i) for i in range(0,plan['cells'],128)]):
            rows+=done.result()['rows'];print('PRESSURE IONIZATION',len(rows),'/',plan['cells'],flush=True)
    save('result.json',dict(classification='Counterexample candidate',passed=True,cells=len(rows),
        original_EOS_and_Coulomb_bitwise=True,maximum_moment_relative_error=max(r['moment_relative_error'] for r in rows),
        maximum_thermodynamic_normalized_error=max(r['thermodynamic_normalized_error'] for r in rows),
        maximum_projected_Coulomb_Gram_error=max(r['projected_Coulomb_Gram_error'] for r in rows),
        minimum_combined_component_margin=min(r['combined_component_margin'] for r in rows),
        nonpositive_margin_cells=sum(r['combined_component_margin']<=0 for r in rows),
        excitation_included=False,full_EOS_Hessian_certified=False,physical_EOS_certified=False,full_GR_evolution=False))
    print('PRESSURE IONIZATION COMPLETE',len(rows),'minimum margin',min(r['combined_component_margin'] for r in rows),flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
