"""Trace actual density-dependent excitation contributions and their curvature."""
from concurrent.futures import ProcessPoolExecutor, as_completed
import ctypes, json, shutil, subprocess, sys
import numpy as np
import sympy as sp
import pressure_ionization_curvature as p

c=p.c;g=p.g;OUT=g.OUT/'excitation-curvature';CACHE=g.CACHE/'excitation-curvature'
NAME='free_eos_direct24_excitation_trace';LIB=CACHE/'build/src'/('lib'+NAME+'.so.1.0.0')
TAGS=[(1,0),(2,0),(2,1),(25,0),(26,1)]
SELECT=np.array([3,4,5,10])


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')
def read(path): return json.loads(path.read_text())


def prepare():
    assert not OUT.exists() and not CACHE.exists();OUT.mkdir();CACHE.mkdir()
    source=CACHE/'source';shutil.copytree(p.CACHE/'source',source);changes={}
    path=source/'src/mod_excitation_block.f90';text=path.read_text()
    text=text.replace('  implicit none','''  implicit none
  integer, save, public :: extrace_count=0, extrace_tag(2)=0, extrace_ids(2,318)=0
  integer, save, public :: extrace_mode=0, extrace_nr=-1
  real(fp_kind), save, public :: extrace_value(6,318)=0._fp_kind
  real(fp_kind), save, public :: extrace_grad(4,318)=0._fp_kind, extrace_hess(4,4,318)=0._fp_kind
  real(fp_kind), save, public :: extrace_hf(4,318)=0._fp_kind, extrace_ht(4,318)=0._fp_kind
  real(fp_kind), save, public :: extrace_moments(4,3)=0._fp_kind''',1);path.write_text(text)
    path=source/'src/excitation_sum.f90';text=path.read_text()
    text=text.replace('  use mod_free_eos_constants, only: pi, c2, electron_mass, h_mass','''  use mod_excitation_block, only: extrace_count, extrace_tag, extrace_ids, &
       extrace_value, extrace_grad, extrace_hess, extrace_hf, extrace_ht, &
       extrace_mode, extrace_nr, extrace_moments
  use mod_free_eos_constants, only: pi, c2, electron_mass, h_mass''',1)
    text=text.replace('  free_sum = 0._fp_kind','''  free_sum = 0._fp_kind
  extrace_count=0;extrace_ids=0;extrace_value=0._fp_kind
  extrace_grad=0._fp_kind;extrace_hess=0._fp_kind
  extrace_hf=0._fp_kind;extrace_ht=0._fp_kind
  extrace_mode=ifexcited;extrace_nr=ifnr
  extrace_moments(:,1)=extrasum([1,2,3,nextrasum-1])
  extrace_moments(:,2)=extrasumf([1,2,3,nextrasum-1])
  extrace_moments(:,3)=extrasumt([1,2,3,nextrasum-1])''',1)
    call='call exsum_component_add(&';parts=text.split(call);assert len(parts)==5
    tags=['[ielement,iz-1]','[1,0]','[25,0]','[26,1]']
    text=parts[0]+''.join('extrace_tag='+tag+'\n              '+call+part for tag,part in zip(tags,parts[1:]))
    head,tail=text.split('subroutine exsum_component_add(',1)
    tail=tail.replace('  ! Arguments','''  use mod_excitation_block, only: extrace_count, extrace_tag, extrace_ids, &
       extrace_value, extrace_grad, extrace_hess, extrace_hf, extrace_ht
  ! Arguments''',1)
    tail=tail.replace('  ! Local variables:','  integer :: trace_i,trace_j,trace_slot\n  ! Local variables:',1)
    tail=tail.replace('  nuvarsum_scale = nuvar/qratio_scale','''  extrace_count=extrace_count+1;trace_slot=extrace_count
  if(trace_slot.gt.318) error stop 'excitation trace capacity exceeded'
  extrace_ids(:,trace_slot)=extrace_tag
  extrace_value(:,trace_slot)=[nuvar,nuvarf,nuvart,qratio/qratio_scale, &
       qratiot/qratio_scale,qratiov/qratio_scale]
  if(ifmhd_logical) then
    do trace_i=1,4
      if(iz_logical.or.trace_i.eq.4) then
        extrace_grad(trace_i,trace_slot)=qratio_dx(trace_i)/qratio_scale
        if(ifnr03) then
          extrace_hf(trace_i,trace_slot)=dqratio_dxf(trace_i)/qratio_scale
          extrace_ht(trace_i,trace_slot)=dqratio_dxt(trace_i)/qratio_scale
        endif
        do trace_j=1,trace_i
          if(iz_logical.or.trace_j.eq.4) then
            extrace_hess(trace_i,trace_j,trace_slot)=qratio_dx2(trace_i,trace_j)/qratio_scale
            extrace_hess(trace_j,trace_i,trace_slot)=extrace_hess(trace_i,trace_j,trace_slot)
          endif
        enddo
      endif
    enddo
  endif
  nuvarsum_scale = nuvar/qratio_scale''',1)
    path.write_text(head+'subroutine exsum_component_add('+tail)
    path=source/'src/CMakeLists.txt';text=path.read_text();assert 'OUTPUT_NAME '+p.NAME in text
    path.write_text(text.replace('OUTPUT_NAME '+p.NAME,'OUTPUT_NAME '+NAME))
    for name in ['mod_excitation_block.f90','excitation_sum.f90','CMakeLists.txt']:
        shutil.copy2(p.CACHE/'source/src'/name,OUT/('before-'+name));shutil.copy2(source/'src'/name,OUT/name)
        changes[name]=dict(before=g.c.sha(p.CACHE/'source/src'/name),after=g.c.sha(source/'src'/name))
    shutil.copy2(c.s.OUT/'inventory_bridge.f90',OUT/'direct_ion_bridge.f90')
    paths=[g.ROOT/'verification/excitation_curvature.py',p.OUT/'manifest.json',g.OUT/'reference-state.npz',OUT/'direct_ion_bridge.f90']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='96f08fa',cells=5735,processes=4,block_cells=128,
        control_cells=[0,2972,5734],source_changes=changes,
        bindings={a.relative_to(g.ROOT).as_posix():g.c.sha(a) for a in paths},original_library_sha256=g.c.sha(c.s.LIB),
        expected_excitation_modes=[2,12],allowed_component_tags=TAGS,
        original_EOS_and_MDH_bitwise_required=True,moment_relative_tolerance=1e-10,
        component_population_relative_tolerance=1e-10,aggregate_normalized_tolerance=1e-10,
        moment_Hessian_direction_tolerance=1e-10,
        method='Capture the final excitation_sum contributions after their log(1+ratio) transformation, divided by their explicit exp(600) storage scale. Read only source-defined lower-triangular Hessian entries and only the fourth coordinate for non-neutral atomic ions. Include H2+ neutral-radius convention.',
        coordinates='Four physical moments alpha0,alpha1,alpha2,beta; at fixed native density scale S and radius unit L, use S*diag(1,L,L^2,1) to transform derivatives to existing nu-moment coordinates.',
        scope='Saved positive species support and fixed branch of the declared native model. Partition cutoffs, infinite-state tails, rounding/continuous errors and physical calibration are not certified.',
        full_EOS_Hessian_certified=False,physical_EOS_certified=False,full_GR_evolution=False))


def build():
    old=g.d.OUT
    try:
        g.d.OUT=OUT;g.d.build_at(CACHE/'source',CACHE/'build',NAME,CACHE/'unused-inventory.so')
    finally: g.d.OUT=old
    command=['gfortran','-O2','-fPIC','-shared','-I'+str(LIB.parent),str(c.OUT/'coulomb_bridge.f90'),
        '-L'+str(LIB.parent),'-Wl,-rpath,'+str(LIB.parent),'-l'+NAME,'-o',str(CACHE/'excitation.so')]
    result=subprocess.run(command,capture_output=True,text=True);(OUT/'excitation-bridge-build.log').write_text(result.stdout+result.stderr)
    assert result.returncode==0,result.stderr
    save('bridge.json',dict(classification='Counterexample candidate',command=command,bridge_sha256=g.c.sha(CACHE/'excitation.so'),library_sha256=g.c.sha(LIB)))


def gram(n,B,atoms,elements):
    n=n.astype(np.longdouble);B=B.astype(np.longdouble);A=atoms.astype(np.longdouble);G=np.zeros((len(B),len(B)),dtype=np.longdouble)
    for element in np.unique(elements):
        selected=elements==element;v=n[selected];a=A[selected];b=B[:,selected]
        center=np.sum(b*(v*a),axis=1)/np.sum(v*a*a);r=b-center[:,None]*a
        G+=(r*v)@r.T
    return np.asarray(G,float)


class EOS(p.EOS):
    def __init__(self):
        super().__init__();self.inventory_lib=ctypes.CDLL(str(CACHE/'excitation.so'))
        array=np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS');call=self.inventory_lib.coulomb_inventory
        call.argtypes=[ctypes.c_int,ctypes.c_double,ctypes.c_double,array,array,ctypes.POINTER(ctypes.c_int)];call.restype=None
        def capture(mode,value,t,eps,out,info):
            raw=np.full(25,np.nan);call(mode,value,t,eps,raw,info)
            out[:]=raw[:22];out[20]=0.;self.molecules=raw[22:24].copy();self.fl=float(raw[24]);self.rho_native=float(raw[0])
        self.call=capture;self.probe_call=self.inventory_lib.coulomb_probe
        self.probe_call.argtypes=[ctypes.c_double]*4+[array];self.probe_call.restype=None

    def trace(self):
        module='mod_excitation_block';prefix='__'+module+'_MOD_'
        def integer(name): return ctypes.c_int.in_dll(self.inventory_lib,prefix+name).value
        count=integer('extrace_count');mode=integer('extrace_mode');nr=integer('extrace_nr')
        assert mode in [2,12] and nr==0 and 0<=count<=len(TAGS)
        ids=np.ctypeslib.as_array((ctypes.c_int*(318*2)).in_dll(self.inventory_lib,prefix+'extrace_ids')).reshape(318,2)[:count].copy()
        assert len({tuple(v) for v in ids})==count and all(tuple(v) in TAGS for v in ids)
        result=dict(count=count,mode=mode,nr=nr,ids=np.zeros((5,2),int))
        result['ids'][:count]=ids
        for name,size,shape in [('value',6,(5,6)),('grad',4,(5,4)),('hess',16,(5,4,4)),('hf',4,(5,4)),('ht',4,(5,4))]:
            source=self.array(module,'extrace_'+name,318*size).reshape((318,*shape[1:]))[:count]
            assert np.all(np.isfinite(source));a=np.zeros(shape);a[:count]=source;result[name]=a
        result['moments']=self.array(module,'extrace_moments',12).reshape(3,4)
        result['sums']=np.array([ctypes.c_double.in_dll(self.inventory_lib,prefix+k).value for k in ['free_sum','usum','ssum','psum']])
        return result

    def full_excitation(self,r,t,x,electron):
        base=super().full(r,t,x,electron);snap,a,meta,extra,rn,ri,G0,H0,computed,native,row=base;trace=self.trace()
        S=meta[1]*self.constants[0];M=trace['moments'];expected=S*extra[:,[0,1,2,7]]
        moment_error=float(np.max(abs(M[0]-expected[0])/np.maximum(abs(expected[0]),1e-300)))
        ym=(x/g.c.A)@self.mapping;eps=ym/float(ym@self.weights);mol=eps[0]*snap['molecular_H_fractions']/2
        n,B,atoms,elements=p.species_coordinates(snap['number_fractions'],mol,rn,ri)
        keys=[(el+1,j) for el,Z in enumerate(g.d.CHARGES) for j in range(Z+1) if snap['number_fractions'][el,j]>0]
        keys += [(25+j,j) for j,v in enumerate(mol) if v>0];assert len(keys)==len(n)
        factor=S*np.array([1.,p.L,p.L*p.L,1.]);D=np.zeros((4,len(n)));W=np.zeros((4,4));aggregate=np.zeros(4)
        population_error=0.;direction_error=0.;pressure_error=0.
        for j in range(trace['count']):
            tag=tuple(trace['ids'][j]);value=trace['value'][j];gradient=trace['grad'][j];hess=trace['hess'][j]
            nu=value[0];expected_nu=n[keys.index(tag)] if tag in keys else 0.
            population_error=max(population_error,float(abs(nu-expected_nu)/max(abs(expected_nu),1e-300)))
            if tag in keys: D[:,keys.index(tag)]=factor*gradient
            W+=nu*(factor[:,None]*hess*factor[None,:])
            aggregate+=nu*np.array([value[3],value[4],value[3]+value[4],value[5]])
            euler=-gradient@M[0];normal=max(abs(value[5]),float(abs(gradient)@abs(M[0])),1e-300)
            pressure_error=max(pressure_error,float(abs(euler-value[5])/normal))
            for direction,target in [(M[1],trace['hf'][j]),(M[2],trace['ht'][j])]:
                predicted=hess@direction;normal=np.maximum(abs(hess)@abs(direction),1e-300)
                direction_error=max(direction_error,float(np.max(abs(predicted-target)/normal)))
        aggregate_error=float(np.max(abs(aggregate-trace['sums']))/max(float(abs(trace['sums']).max()),1e-300))
        B=np.vstack([B,D]);G=gram(n,B,atoms,elements)
        Hnon=np.zeros((16,16));Hnon[:12,:12]=H0;Hnon[0,0]=S*a[4]/a[15]
        Hnon[np.ix_(SELECT,SELECT)]-=W
        for j,k in enumerate(SELECT): Hnon[k,12+j]=Hnon[12+j,k]=-1.
        H=Hnon.copy();H[0,0]+=S*(a[13]+electron[2])/(a[10]*a[12])
        eig,Q=np.linalg.eigh(G);assert eig.min()>=-1e-12*max(1.,eig.max())
        root=(Q*np.sqrt(np.maximum(eig,0)))@Q.T;C=root@H@root;C=(C+C.T)/2
        margin=float(1+min(0.,np.linalg.eigvalsh(C)[0]))
        row.update(excitation_moment_relative_error=moment_error,component_population_relative_error=population_error,
            aggregate_normalized_error=aggregate_error,moment_Hessian_direction_error=direction_error,
            pressure_Euler_normalized_error=pressure_error,combined_excitation_margin=margin,excitation_components=trace['count'])
        return base,trace,G,H,Hnon,row


def symbolic():
    n0,n1=sp.symbols('n0 n1',positive=True);m=sp.symbols('m',real=True);n=sp.Matrix([n0,n1]);C=sp.Matrix([[1,2]])
    L=[sp.log(1+sp.exp(-m)),m*m/3];F=-sum(n[i]*L[i].subs(m,(C*n)[0]) for i in range(2))
    D=sp.Matrix([[sp.diff(v,m).subs(m,(C*n)[0]) for v in L]])
    W=sp.Matrix([[sum(n[i]*sp.diff(L[i],m,2).subs(m,(C*n)[0]) for i in range(2))]])
    assert (sp.hessian(F,[n0,n1])+D.T*C+C.T*D+C.T*W*C).applyfunc(sp.simplify)==sp.zeros(2,2)
    save('symbolic.json',dict(classification='Proven',passed=True,
        free_energy='F_exc/(kT)=-sum_j n_j L_j(M), M=C*n.',
        Hessian='H_exc=-(D^T C+C^T D+C^T W C), D_aj=dL_j/dM_a, W_ab=sum_j n_j d2L_j/dM_a dM_b.',
        scope='Fixed T, linear physical moments and differentiable declared branch. A two-species nonlinear manufactured example was differentiated independently. Hard cutoffs, sum tails and numerical derivative correctness are separate.'))


def gates(row,plan):
    return row['excitation_moment_relative_error']<plan['moment_relative_tolerance'] and row['component_population_relative_error']<plan['component_population_relative_tolerance'] and max(row['aggregate_normalized_error'],row['pressure_Euler_normalized_error'])<plan['aggregate_normalized_tolerance'] and row['moment_Hessian_direction_error']<plan['moment_Hessian_direction_tolerance']


def control():
    plan=read(OUT/'plan.json');state=dict(np.load(g.OUT/'reference-state.npz'));eos=EOS();rows=[];symbolic()
    for i in plan['control_cells']:
        electron=dict(np.load(p.ex.OUT/f'control-{i}.npz'))['values']
        base,trace,G,H,Hnon,row=eos.full_excitation(state['lnd'][i],state['lnT'][i],state['X'][i],electron)
        snap,a,meta,extra,rn,ri,G0,H0,computed,native,_=base;old=dict(np.load(p.OUT/f'control-{i}.npz'))
        for k,v in dict(state=meta,extra=extra,neutral=rn,ion3=ri,Gram=G0,Hessian=H0,computed=computed,native=native).items(): assert np.array_equal(v,old[k]),(i,k)
        assert np.array_equal(snap['eos'],dict(np.load(c.s.OUT/f'control-{i}.npz'))['eos'])
        row.update(cell=i,passed=gates(row,plan));rows.append(row)
        np.savez_compressed(OUT/f'control-{i}.npz',**trace,Gram=G,Hessian=H,nonideal_Hessian=Hnon)
        print('EXCITATION CONTROL',row,flush=True)
    save('control.json',dict(classification='Counterexample candidate',passed=all(r['passed'] for r in rows),rows=rows));assert all(r['passed'] for r in rows)


def block(start):
    plan=read(OUT/'plan.json');state=dict(np.load(g.OUT/'reference-state.npz'));eos=EOS();rows=[];arrays=[]
    electron=dict(np.load(p.ex.OUT/f'block-{start}.npz'));old=dict(np.load(p.OUT/f'block-{start}.npz'));stop=min(start+128,plan['cells'])
    for i in range(start,stop):
        base,trace,G,H,Hnon,row=eos.full_excitation(state['lnd'][i],state['lnT'][i],state['X'][i],electron['values'][i-start])
        snap,a,meta,extra,rn,ri,G0,H0,computed,native,_=base
        for k,v in dict(state=meta,extra=extra,neutral=rn,ion3=ri,Gram=G0,Hessian=H0,computed=computed,native=native).items(): assert np.array_equal(v,old[k][i-start]),(i,k)
        row.update(cell=i,passed=gates(row,plan));rows.append(row);arrays.append(dict(**trace,Gram=G,Hessian=H,nonideal_Hessian=Hnon))
    path=OUT/f'block-{start}.npz';np.savez_compressed(path,**{k:np.array([a[k] for a in arrays]) for k in arrays[0]})
    record=dict(classification='Counterexample candidate',start=start,stop=stop,rows=rows,passed=all(r['passed'] for r in rows),
        output_sha256=g.c.sha(path),plan_sha256=g.c.sha(OUT/'plan.json'))
    save(f'block-{start}.json',record);assert record['passed'],start;return record


def run():
    plan=read(OUT/'plan.json');assert read(OUT/'control.json')['passed']
    for rel,digest in plan['bindings'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    bridge=read(OUT/'bridge.json');assert g.c.sha(LIB)==bridge['library_sha256'] and g.c.sha(CACHE/'excitation.so')==bridge['bridge_sha256']
    assert g.c.sha(c.s.LIB)==plan['original_library_sha256'];rows=[]
    with ProcessPoolExecutor(max_workers=plan['processes']) as pool:
        for done in as_completed([pool.submit(block,i) for i in range(0,plan['cells'],128)]):
            rows+=done.result()['rows'];print('EXCITATION',len(rows),'/',plan['cells'],flush=True)
    save('result.json',dict(classification='Counterexample candidate',passed=True,cells=len(rows),
        original_MDH_outputs_and_matrices_bitwise=True,
        maximum_component_population_relative_error=max(r['component_population_relative_error'] for r in rows),
        maximum_aggregate_normalized_error=max(r['aggregate_normalized_error'] for r in rows),
        maximum_moment_Hessian_direction_error=max(r['moment_Hessian_direction_error'] for r in rows),
        maximum_pressure_Euler_normalized_error=max(r['pressure_Euler_normalized_error'] for r in rows),
        minimum_combined_excitation_margin=min(r['combined_excitation_margin'] for r in rows),
        nonpositive_margin_cells=sum(r['combined_excitation_margin']<=0 for r in rows),
        full_EOS_Hessian_certified=False,physical_EOS_certified=False,full_GR_evolution=False))
    print('EXCITATION COMPLETE',len(rows),'minimum margin',min(r['combined_excitation_margin'] for r in rows),flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
