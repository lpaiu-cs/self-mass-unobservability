"""Test the original GR EOS with upstream extended precision and split inputs."""
from types import FunctionType
import ctypes, json, shutil, sys
import mpmath as mp
import numpy as np
import gr_subcell_root_audit as original

g=original.g;OUT=g.OUT/'gr-eos-extended-precision';CACHE=g.CACHE/'extended-precision'
NAME='free_eos_direct24_extended';LIB=CACHE/'build/src'/('lib'+NAME+'.so.1.0.0')
BRIDGE=CACHE/'split-input.so'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists() and not CACHE.exists();OUT.mkdir();CACHE.mkdir();original.verify()
    source=CACHE/'source';parent=g.d.CACHE/'full-integral-source';shutil.copytree(parent,source)
    path=source/'src/CMakeLists.txt';before=path.read_text();after=before
    replacements=[('set(supported_fp_precision 15 CACHE','set(supported_fp_precision 33 CACHE'),
        ('set(supported_fp_exponent 300 CACHE','set(supported_fp_exponent 4000 CACHE'),
        ('OUTPUT_NAME free_eos_direct24_integral_full','OUTPUT_NAME '+NAME)]
    for old,new in replacements:assert after.count(old)==1;after=after.replace(old,new)
    reverse=after
    for old,new in reversed(replacements):assert reverse.count(new)==1;reverse=reverse.replace(new,old)
    assert reverse==before;path.write_text(after)
    shutil.copy2(parent/'src/CMakeLists.txt',OUT/'before-CMakeLists.txt');shutil.copy2(path,OUT/path.name)
    save('source-tree.json',dict(before={p.relative_to(parent).as_posix():g.c.sha(p) for p in parent.rglob('*') if p.is_file()},
        after={p.relative_to(source).as_posix():g.c.sha(p) for p in source.rglob('*') if p.is_file()},substitutions=replacements))
    bridge='''! Extended arithmetic candidate. No physical EOS or native error certificate.
subroutine direct_ion_eos(kif, vh, vl, th, tlow, cx, eps, hi, lo, info) bind(C)
  use iso_c_binding
  use mod_free_eos_types, only: fp_kind
  use mod_free_eos, only: free_eos
  use mod_free_eos_constants, only: c2, cr
  use mod_ionization_data, only: h2diss
  implicit none
  integer(c_int), value :: kif
  real(c_double), value :: vh, vl, th, tlow, cx
  real(c_double), intent(in) :: eps(24)
  real(c_double), intent(out) :: hi(22),lo(22)
  integer(c_int), intent(out) :: info
  integer :: iterations
  real(fp_kind) :: match_value,tl,epsq(24),res(22),cxq
  real(fp_kind) :: fl,t,rho,rl,p,pl,cf,cp,s,sf,st,grada,rtp,qe,qv,rmue, &
       fh2,fhe2,fhe3,xmu1,xmu3,eta,gamma1,gamma2,gamma3,h2rat,h2plusrat, &
       lambda,gamma_e,sound2,pressure(3),density(3),energy(3),entropy(3)
  match_value=real(vh,fp_kind)+real(vl,fp_kind)
  tl=real(th,fp_kind)+real(tlow,fp_kind)
  epsq=real(eps,fp_kind);cxq=real(cx,fp_kind)
  ! The disabled diffraction output is padding, discarded by the Python caller.
  gamma_e=0._fp_kind
  call free_eos(0,3,1,-2,kif,epsq,match_value,tl,fl,t,rho,rl,p,pl, &
       cf,cp,s,sf,st,grada,rtp,qe,qv,rmue,fh2,fhe2,fhe3,xmu1,xmu3,eta, &
       gamma1,gamma2,gamma3,h2rat,h2plusrat,lambda,gamma_e,sound2, &
       iterations,info,pressure=pressure,density=density,energy=energy,entropy=entropy)
  if (info /= 0) return
  res(:12) = [rho,p,qe-0.5_fp_kind*c2*cr*h2diss*epsq(1),s,gamma1, &
       pressure(2),pressure(3),density(2),density(3),energy(2),energy(3), &
       0.5_fp_kind*c2*cr*h2diss*epsq(1)]
  res(13:) = [eta,rmue,fh2,fhe2,fhe3,xmu1,xmu3,lambda,gamma_e,sound2]
  res(1)=res(1)/cxq;res(3:4)=res(3:4)*cxq;res(10:11)=res(10:11)*cxq
  hi=real(res,c_double);lo=real(res-real(hi,fp_kind),c_double)
end subroutine direct_ion_eos

subroutine precision_info(values) bind(C)
  use iso_c_binding
  use mod_free_eos_types, only: fp_kind,lapack_fp_kind
  implicit none
  integer(c_int),intent(out) :: values(5)
  values=[precision(1._fp_kind),digits(1._fp_kind),range(1._fp_kind), &
       storage_size(1._fp_kind),digits(1._lapack_fp_kind)]
end subroutine precision_info
'''
    (OUT/'direct_ion_bridge.f90').write_text(bridge)
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='1252b97',
        intervention='Use the original full-integral GR EOS source, not the newer molecular model. Select the upstream 33-digit/4000-exponent real kind. Keep coefficients, physical options, solver and Fermi tolerances unchanged. Split binary64 high/low inputs and outputs carry extended values across the C ABI; composition normalization follows the old binary64 map and its exact lifted values.',
        root_iterations=15,energy_gate='Unchanged max(2 erg/g,32 ulp(float(abs(u+P/rho)))); require score <=1, aim <0.25.',
        known_target_cells=[0,1175,2972,5734],known_root_logT_tolerance=1e-10,
        comparison='Solve the preserved failed entropy target with this new evaluator. Record its difference from the old evaluator at the same high input, and the residual when its refined root input is rounded back to one binary64. Never relabel the original evaluator failure as passed.',
        boundary='Upstream LAPACK remains binary64. Higher floating precision does not raise source-constant accuracy, eliminate physical-model error, or certify Fermi quadrature/native/continuous-root error. No original GR library, state or running path is replaced.',
        runtime_original=json.loads((original.OUT/'plan.json').read_text())['runtime'],
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in [g.ROOT/'verification/gr_eos_extended_precision.py',
            original.OUT/'manifest.json',OUT/'source-tree.json',OUT/'direct_ion_bridge.f90',g.OUT/'initial-state-17-4.npz']}))


def build():
    fn=g.d.build_at
    FunctionType(fn.__code__,dict(fn.__globals__,OUT=OUT,save=save))(CACHE/'source',CACHE/'build',NAME,BRIDGE)
    lib=ctypes.CDLL(str(BRIDGE));a=np.zeros(5,dtype=np.int32)
    lib.precision_info.argtypes=[np.ctypeslib.ndpointer(np.int32,flags='C_CONTIGUOUS')];lib.precision_info(a)
    assert list(a)==[33,113,4931,128,53],a
    save('precision.json',dict(classification='Counterexample candidate',decimal_digits=int(a[0]),binary_significand_bits=int(a[1]),
        decimal_exponent_range=int(a[2]),storage_bits=int(a[3]),LAPACK_significand_bits=int(a[4])))


class EOS:
    def __init__(self):
        self.lib=ctypes.CDLL(str(BRIDGE));self.call=self.lib.direct_ion_eos
        array=np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS')
        self.call.argtypes=[ctypes.c_int]+[ctypes.c_double]*5+[array]*3+[ctypes.POINTER(ctypes.c_int)]
        self.call.restype=None;old=g.EOS();self.mapping=old.mapping;self.weights=old.weights

    def raw(self,mode,value,t,eps,cx):
        hi=np.full(22,np.nan);lo=hi.copy();info=ctypes.c_int(-999)
        vh,th=float(value),float(t);vl,tl=float(mp.mpf(value)-vh),float(mp.mpf(t)-th)
        self.call(mode,vh,vl,th,tl,cx,np.ascontiguousarray(eps),hi,lo,ctypes.byref(info))
        assert info.value==0 and np.all(np.isfinite(hi)) and np.all(np.isfinite(lo)),(info.value,value,t)
        return [mp.mpf(float(a))+mp.mpf(float(b)) for i,(a,b) in enumerate(zip(hi,lo)) if i!=20]

    def __call__(self,lp,t,X):
        ym=(X/g.c.A)@self.mapping;cx=float(ym@self.weights);eps=ym/cx
        seed=np.zeros(24);seed[2 if X[g.c.NAMES.index('c12')]<.5 else 0]=1;seed/=seed@self.weights
        self.raw(0,mp.mpf(-20),mp.mpf(float(np.log(1e6))),seed,1.)
        return self.raw(1,lp,t,eps,cx)


def solve(eos,lp,target,X,guess,iterations):
    t=mp.mpf(guess);target=mp.mpf(target);rows=[];best=None
    for _ in range(iterations):
        a=eos(lp,t,X);T=mp.exp(t);H=a[2]+a[1]/a[0];C=a[10]-a[1]/a[0]*a[8]
        assert C>0;budget=max(2.,32*np.spacing(abs(float(H))));error=a[3]-target
        score=abs(T*error)/budget;rows.append(dict(lnT=mp.nstr(t,80),entropy_error=mp.nstr(error,80),score=float(score),budget_erg_g=float(budget)))
        if best is None or score<best[0]:best=score,t,a
        if score<mp.mpf('.25'):break
        t-=max(mp.mpf('-.15'),min(mp.mpf('.15'),error*T/C))
    return best,rows


def run():
    mp.mp.dps=80;plan=json.loads((OUT/'plan.json').read_text());failure=json.loads((original.OUT/'original-failure.json').read_text())
    eos=EOS();X=np.array(failure['X']);lp=mp.mpf(failure['logP']);t0=mp.mpf(failure['best_lnT'])
    held=eos(lp,t0,X);old=g.EOS()(1,float(lp),float(t0),X)
    best,trace=solve(eos,lp,mp.mpf(failure['entropy']),X,failure['guess'],plan['root_iterations'])
    rounded=eos(lp,mp.mpf(float(best[1])),X)
    rounded_score=abs(mp.exp(float(best[1]))*(rounded[3]-failure['entropy']))/trace[-1]['budget_erg_g']
    result=dict(classification='Counterexample candidate',failed_target_trace=trace,failed_target_new_evaluator_passed=bool(best[0]<=1),
        refined_lnT=mp.nstr(best[1],80),refined_score=float(best[0]),rounded_input_new_evaluator_score=float(rounded_score),
        original_evaluator_same_input_differences=[mp.nstr(a-mp.mpf(float(b)),40) for a,b in zip(held,old)],
        old_target_and_old_evaluator_failure_preserved=True,physical_or_continuous_root_certified=False)
    save('failed-target-result.json',result);print('EXTENDED ROOT',float(best[0]),float(rounded_score),flush=True)
    state=dict(np.load(g.OUT/'initial-state-17-4.npz'));controls=[]
    for i in plan['known_target_cells']:
        X=state['X'][i];lt=mp.mpf(float(state['lnT'][i]));lp=mp.mpf(float(np.log(state['pressure'][i]))) if 'pressure' in state else mp.mpf(float(state['logP'][i]))
        known=eos(lp,lt,X);best,trace=solve(eos,lp,known[3],X,lt+mp.mpf('.01'),plan['root_iterations'])
        error=abs(best[1]-lt);passed=best[0]<=1 and error<plan['known_root_logT_tolerance']
        controls.append(dict(cell=i,trace=trace,score=float(best[0]),lnT_error=mp.nstr(error,40),passed=bool(passed)))
        print('EXTENDED KNOWN ROOT',i,float(best[0]),float(error),flush=True)
    save('controls.json',dict(classification='Counterexample candidate',rows=controls,all_passed=all(r['passed'] for r in controls)))
    save('result.json',dict(classification='Counterexample candidate',completed=True,
        finite_passed=result['failed_target_new_evaluator_passed'] and all(r['passed'] for r in controls),
        physical_EOS_certified=False,continuous_native_or_root_error_certified=False,original_GR_path_replaced=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify()


def verify():
    for name,key in [('plan.json','bindings'),('manifest.json','sha256')]:
        for rel,digest in json.loads((OUT/name).read_text())[key].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for name,key in [('plan.json','runtime_original'),('build.json','sha256')]:
        for path,digest in json.loads((OUT/name).read_text())[key].items():assert g.c.sha(path)==digest,path
    trees=json.loads((OUT/'source-tree.json').read_text())
    for name,folder in [('before',g.d.CACHE/'full-integral-source'),('after',CACHE/'source')]:
        for rel,digest in trees[name].items():assert g.c.sha(folder/rel)==digest,rel
    assert json.loads((OUT/'result.json').read_text())['completed']
    print('PASS extended-precision candidate bindings; consult finite_passed, no physical or continuum certificate',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
