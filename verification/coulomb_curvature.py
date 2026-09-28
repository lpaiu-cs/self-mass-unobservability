"""Actual EOS1 Coulomb value, tangent and constrained curvature controls.

Saved positive populations only; exchange, pressure ionization and the
global nonideal feasible domain are not certified by this component audit.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import ctypes, json, shutil, subprocess, sys
import numpy as np
import sympy as sp
import eos_species_inventory as s

g=s.g;OUT=g.OUT/'coulomb-curvature';CACHE=g.CACHE/'coulomb-curvature'
FIELDS=['free','mu_e','mu_0','mu_2','H_ee','H_e0','H_e2','H_00','H_02','H_22',
    'ne','pe','dlnne_dlnf','deta_dlnf','lambda','kT','H_0e','H_2e','theta_e']


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();CACHE.mkdir(exist_ok=True)
    paths=[g.ROOT/'verification/coulomb_curvature.py',s.OUT/'manifest.json',s.OUT/'inventory_bridge.f90',
        g.ROOT/'verification/eos_species_inventory.py',g.OUT/'reference-state.npz']
    for name in ['mod_free_eos.f90','mod_coulomb.f90','master_coulomb.f90','coulomb.f90','mod_free_eos_constants.f90']:
        target=OUT/name;shutil.copy2(g.d.CACHE/'full-integral-source/src'/name,target);paths.append(target)
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='8f878bb',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},library_sha256=g.c.sha(s.LIB),
        state='outputs/direct-eos-gr33/reference-state.npz',cells=5735,processes=4,block_cells=128,
        control_cells=[0,2972,5734],same_EOS_free_value_relative_tolerance=1e-9,
        same_electron_density_relative_tolerance=1e-10,normalized_Hessian_symmetry_tolerance=1e-12,
        finite_log_steps=[1e-4,5e-5],normalized_Hessian_derivative_tolerance=1e-5,
        settings=dict(fermi_morder=21,ifcoulomb=5,if_dc=0,if_pteh=0,master_ifnr=3),
        hypothesis='EOS1 (3,1,-2) uses non-PTEH Coulomb sums with no metal-sum approximation, ifcoulomb=5 and no diffraction. Reconstruct sum0,sum2 from the actual saved ion and H2+ populations, and verify the Coulomb free energy before interpreting curvature.',
        derivative_coordinates='Independent ln f, ln sum0, ln sum2 at fixed temperature; convert the ln f column to ln ne using dlnne/dlnf. Normalize the physical Hessian as diag(ne,sum0,sum2) H diag(ne,sum0,sum2)/(kT ne).',
        curvature='Use ideal-ion Fisher metric diag(1/n_species), ideal electron compressibility, and actual Coulomb Hessian; project onto each element inventory, counting two H nuclei in H2/H2+. Report signed eigenvalues on saved positive population support, with no pass criterion implying whole-EOS convexity.',
        full_EOS_Hessian_certified=False,physical_EOS_certified=False,full_GR_evolution=False))


def build():
    text=(s.OUT/'inventory_bridge.f90').read_text().replace('ionization_inventory','coulomb_inventory').replace('res(24)','res(25)')
    text=text.replace('end subroutine coulomb_inventory','  res(25) = fl\nend subroutine coulomb_inventory')
    text+='''
subroutine native_constants(res) bind(C)
  use iso_c_binding
  use mod_free_eos_constants, only: avogadro,boltzmann,ct,c_e,cpe
  real(c_double), intent(out) :: res(5)
  res=[avogadro,boltzmann,ct,c_e,cpe]
end subroutine native_constants

subroutine coulomb_probe(fl,tl,sum0,sum2,res) bind(C)
  use iso_c_binding
  use mod_free_eos_constants, only: boltzmann,ct,c_e,cpe
  use mod_fermi_dirac, only: fermi_dirac
  use mod_coulomb, only: master_coulomb
  use mod_master_coulomb_data, only: fcoulomb
  implicit none
  real(c_double), value :: fl,tl,sum0,sum2
  real(c_double), intent(out) :: res(19)
  real(c_double) :: r(9),p(9),ss(3),uu(3),f,t,ne,pe,kt,lam,gamma_e,den, &
       dve,dvef,dvet,dve0,dve2,dv0,dv0f,dv0t,dv2,dv2f,dv2t,dv00,dv02,dv22
  call fermi_dirac(0,fl,tl+log(ct),r,p,ss,uu,21)
  f=exp(fl);t=exp(tl);ne=c_e*r(1);pe=cpe*p(1);kt=boltzmann*t
  call master_coulomb(3,r,f,sum0,0.d0,0.d0,0.d0,sum2,0.d0,0.d0,0.d0, &
       ne,t,pe,p,lam,gamma_e,5,0,0,dve,dvef,dvet,dve0,dve2, &
       dv0,dv0f,dv0t,dv2,dv2f,dv2t,dv00,dv02,dv22)
  den=ne*r(2)
  res=[fcoulomb,-kt*dve,-kt*dv0,-kt*dv2,-kt*dvef/den, &
       -kt*dve0,-kt*dve2,-kt*dv00,-kt*dv02,-kt*dv22, &
       ne,pe,r(2),sqrt(1.d0+f),lam,kt,-kt*dv0f/den,-kt*dv2f/den,r(2)/sqrt(1.d0+f)]
end subroutine coulomb_probe
'''
    source=OUT/'coulomb_bridge.f90';source.write_text(text)
    command=['gfortran','-O2','-fPIC','-shared','-I'+str(s.LIB.parent),str(source),
        '-L'+str(s.LIB.parent),'-Wl,-rpath,'+str(s.LIB.parent),'-lfree_eos_direct24_integral_full','-o',str(CACHE/'coulomb.so')]
    result=subprocess.run(command,capture_output=True,text=True);(OUT/'build.log').write_text(result.stdout+result.stderr)
    assert result.returncode==0,result.stderr
    plan=json.loads((OUT/'plan.json').read_text());assert g.c.sha(s.LIB)==plan['library_sha256']
    save('bridge.json',dict(classification='Counterexample candidate',command=command,source_sha256=g.c.sha(source),
        bridge_sha256=g.c.sha(CACHE/'coulomb.so'),original_library_unchanged=True))


class EOS(s.InventoryEOS):
    def __init__(self):
        super().__init__();self.inventory_lib=ctypes.CDLL(str(CACHE/'coulomb.so'));call=self.inventory_lib.coulomb_inventory
        array=np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS')
        call.argtypes=[ctypes.c_int,ctypes.c_double,ctypes.c_double,array,array,ctypes.POINTER(ctypes.c_int)];call.restype=None
        def capture(mode,value,t,eps,out,info):
            raw=np.full(25,np.nan);call(mode,value,t,eps,raw,info)
            out[:]=raw[:22];out[20]=0.;self.molecules=raw[22:24].copy();self.fl=float(raw[24]);self.rho_native=float(raw[0])
        self.call=capture;self.constants=np.empty(5)
        self.inventory_lib.native_constants.argtypes=[array];self.inventory_lib.native_constants.restype=None
        self.inventory_lib.native_constants(self.constants);assert np.all(self.constants>0)
        self.probe_call=self.inventory_lib.coulomb_probe;self.probe_call.argtypes=[ctypes.c_double]*4+[array];self.probe_call.restype=None

    def probe(self,fl,tl,sum0,sum2):
        out=np.full(len(FIELDS),np.nan);self.probe_call(float(fl),float(tl),float(sum0),float(sum2),out)
        assert np.all(np.isfinite(out));return out

    def actual(self,r,t,x):
        snapshot=self.snapshot(r,t,x)
        old_free=ctypes.c_double.in_dll(self.inventory_lib,'__mod_master_coulomb_data_MOD_fcoulomb').value
        ym=(x/g.c.A)@self.mapping;cx=float(ym@self.weights);eps=ym/cx
        scale=self.rho_native*self.constants[0]
        populations=snapshot['number_fractions']*scale
        molecules=eps[0]*snapshot['molecular_H_fractions']*scale/2
        sum0=float(populations[:,1:].sum()+molecules[1]);sum2=float((populations*np.arange(29)**2).sum()+molecules[1])
        fl=self.fl;values=self.probe(fl,t,sum0,sum2)
        expected_ne=snapshot['eos'][13]*self.constants[0]
        q=np.array([values[10],sum0,sum2]);H=matrix(values)
        normalized=H*q[:,None]*q[None,:]/(values[15]*q[0])
        symmetry=np.array([values[5]-values[16],values[6]-values[17]])*q[1:]/values[15]
        report=dict(original_free=float(old_free),reconstructed_free=float(values[0]),
            free_relative_error=float(abs(values[0]/old_free-1)),electron_density_relative_error=float(abs(values[10]/expected_ne-1)),
            normalized_symmetry_error=float(abs(symmetry).max()))
        return snapshot,values,normalized,np.array([fl,t,sum0,sum2]),populations,molecules,report


def matrix(a): return np.array([[a[4],a[5],a[6]],[a[5],a[7],a[8]],[a[6],a[8],a[9]]])


def tangent_gram(populations,molecules):
    """Fisher-metric projection, evaluated as a positive sum of outer products."""
    G=np.zeros((3,3),dtype=np.longdouble)
    for i,row in enumerate(populations):
        n=row.astype(np.longdouble);z=np.arange(29,dtype=np.longdouble);atoms=np.ones(29,dtype=np.longdouble)
        if i==0:
            n=np.r_[n,molecules.astype(np.longdouble)];z=np.r_[z,0.,1.];atoms=np.r_[atoms,2.,2.]
        if not n.any(): continue
        b=np.array([z,(z>0).astype(np.longdouble),z*z]);den=np.sum(n*atoms*atoms)
        center=np.sum(b*(n*atoms),axis=1)/den;v=b-center[:,None]*atoms
        G+=(v*n)@v.T
    return np.asarray(G,dtype=float)


def curvature(populations,molecules,a):
    G=tangent_gram(populations,molecules);ev,Q=np.linalg.eigh(G)
    assert ev.min()>=-1e-12*max(1.,ev.max())
    root=(Q*np.sqrt(np.maximum(ev,0)))@Q.T;H=matrix(a)/a[15]
    C=root@H@root;C=(C+C.T)/2
    H[0,0]+=a[13]/(a[10]*a[12]);total=root@H@root;total=(total+total.T)/2
    cmin=float(np.linalg.eigvalsh(C)[0]);emin=float(np.linalg.eigvalsh(total)[0])
    return np.array([1+min(0.,cmin),1+min(0.,emin)]),G


def gates(row,plan):
    return row['free_relative_error']<plan['same_EOS_free_value_relative_tolerance'] and row['electron_density_relative_error']<plan['same_electron_density_relative_tolerance'] and row['normalized_symmetry_error']<plan['normalized_Hessian_symmetry_tolerance']


def symbolic():
    # Include H2 and H2+ in the two-element inventory constraint.
    n=sp.diag(1,4,9,16,25,36,49);A=sp.Matrix([[1,1,2,2,0,0,0],[0,0,0,0,1,1,1]])
    B=sp.Matrix([[0,1,0,1,0,1,2],[0,1,0,1,0,1,1],[0,1,0,1,0,1,4]])
    root=sp.diag(1,2,3,4,5,6,7);P=sp.eye(7)-root*A.T*(A*n*A.T).inv()*A*root
    assert P*P==P and A*root*P==sp.zeros(2,7)
    gram=B*n*B.T-B*n*A.T*(A*n*A.T).inv()*A*n*B.T
    assert gram==B*root*P*root*B.T
    residue=B-B*n*A.T*(A*n*A.T).inv()*A
    assert residue*n*residue.T==gram
    save('symbolic.json',dict(classification='Proven',passed=True,
        constrained_Gram='G=B diag(n) B^T-B diag(n) A^T [A diag(n) A^T]^-1 A diag(n) B^T. The centered positive outer-product formula is algebraically identical, including two H nuclei per molecule.',
        spectrum='On the saved positive population support, congruence by sqrt(diag(n)) turns ideal-ion translational entropy into identity. Nonzero eigenvalues of the projected rank-three Coulomb/electron correction equal those of sqrt(G) H sqrt(G)/(kT). Include zero eigenvalues when the tangent dimension exceeds this rank.',
        scope='Algebraic statement for specified populations and component Hessians; it neither validates their floating evaluation nor certifies the omitted EOS components or domain.'))


def control():
    plan=json.loads((OUT/'plan.json').read_text());state=dict(np.load(g.ROOT/plan['state']));eos=EOS();rows=[];symbolic()
    for i in plan['control_cells']:
        snap,a,H,arg,pops,mols,row=eos.actual(state['lnd'][i],state['lnT'][i],state['X'][i])
        old=dict(np.load(s.OUT/f'control-{i}.npz'));assert all(np.array_equal(snap[k],old[k]) for k in old)
        q=np.array([a[10],arg[2],arg[3]]);finite=[]
        for h in plan['finite_log_steps']:
            columns=[]
            for j in range(3):
                sides=[]
                for sign in [-1,1]:
                    v=arg.copy()
                    if j==0: v[0]+=sign*h
                    else: v[j+1]*=np.exp(sign*h)
                    b=eos.probe(*v);sides.append(b[1:4]/b[15])
                column=(sides[1]-sides[0])/(2*h)
                if j==0: column/=a[12]
                columns.append(column*q/q[0])
            finite.append(np.array(columns).T)
        finite=np.array(finite);score=float(np.max(abs(finite-H)/np.maximum(1,abs(H))))
        margins,G=curvature(pops,mols,a)
        row.update(cell=i,normalized_finite_Hessian_error=score,component_curvature_margins=margins.tolist())
        row['passed']=gates(row,plan) and score<plan['normalized_Hessian_derivative_tolerance'];rows.append(row)
        np.savez_compressed(OUT/f'control-{i}.npz',values=a,normalized_Hessian=H,finite_Hessians=finite,arguments=arg,Gram=G,margins=margins)
        print('COULOMB CONTROL',row,flush=True)
    save('control.json',dict(classification='Counterexample candidate',passed=all(r['passed'] for r in rows),rows=rows,
        original_species_and_EOS_outputs_bitwise=True));assert all(r['passed'] for r in rows)


def block(start):
    plan=json.loads((OUT/'plan.json').read_text());state=dict(np.load(g.ROOT/plan['state']));stop=min(start+plan['block_cells'],len(state['X']));eos=EOS()
    rows=[];values=[];margins=[];grams=[];outputs=[]
    for i in range(start,stop):
        snap,a,H,arg,pops,mols,row=eos.actual(state['lnd'][i],state['lnT'][i],state['X'][i]);row.update(cell=i,passed=gates(row,plan))
        margin,gram=curvature(pops,mols,a);rows.append(row);values.append(a);margins.append(margin);grams.append(gram);outputs.append(snap['eos'])
    target=OUT/f'block-{start}.npz';np.savez_compressed(target,values=np.array(values),margins=np.array(margins),Gram=np.array(grams),eos=np.array(outputs))
    result=dict(classification='Counterexample candidate',start=start,stop=stop,rows=rows,passed=all(r['passed'] for r in rows),
        output_sha256=g.c.sha(target),state_sha256=g.c.sha(g.ROOT/plan['state']),plan_sha256=g.c.sha(OUT/'plan.json'))
    save(f'block-{start}.json',result);assert result['passed'],start;return result


def run():
    plan=json.loads((OUT/'plan.json').read_text());assert json.loads((OUT/'control.json').read_text())['passed']
    for rel,digest in plan['bindings'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    assert g.c.sha(s.LIB)==plan['library_sha256']
    assert g.c.sha(CACHE/'coulomb.so')==json.loads((OUT/'bridge.json').read_text())['bridge_sha256']
    records=[]
    with ProcessPoolExecutor(max_workers=plan['processes']) as pool:
        for done in as_completed([pool.submit(block,i) for i in range(0,plan['cells'],plan['block_cells'])]):
            records.append(done.result());save('progress.json',dict(completed_cells=sum(r['stop']-r['start'] for r in records),total_cells=plan['cells']))
            print('COULOMB CURVATURE',len(records),'/',45,flush=True)
    records.sort(key=lambda r:r['start']);parts=[dict(np.load(OUT/f"block-{r['start']}.npz")) for r in records]
    arrays={k:np.concatenate([p[k] for p in parts]) for k in parts[0]};rows=[r for rec in records for r in rec['rows']]
    table=dict(np.load(g.OUT/'initial-adiabats-17.npz'));assert np.array_equal(arrays['eos'],table['reference'])
    save('result.json',dict(classification='Counterexample candidate',passed=True,cells=len(rows),all_21_EOS_outputs_bitwise_equal_independent_table=True,
        maximum_replay_free_relative_error=max(r['free_relative_error'] for r in rows),maximum_electron_density_relative_error=max(r['electron_density_relative_error'] for r in rows),
        maximum_normalized_symmetry_error=max(r['normalized_symmetry_error'] for r in rows),
        margin_labels=['ideal_ions_plus_Coulomb','ideal_ions_plus_ideal_electrons_plus_Coulomb'],
        minimum_component_margins=arrays['margins'].min(0).tolist(),nonpositive_component_margin_cells=(arrays['margins']<=0).sum(0).tolist(),
        omitted=['Exchange','Pressure ionization','Density-dependent partition/excitation','Zero-mask uncertainty','Global feasible-domain and physical-model error'],
        full_EOS_Hessian_certified=False,physical_EOS_certified=False,full_GR_evolution=False))
    print('COULOMB CURVATURE COMPLETE',len(rows),arrays['margins'].min(0),flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
