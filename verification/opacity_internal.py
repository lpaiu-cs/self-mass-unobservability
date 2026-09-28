"""Isolate the live opacity branch and its single-precision boundary."""
import json, shutil, subprocess, sys
import numpy as np
import sympy as sp
import native_opacity as o

g=o.g;OUT=o.OUT/'internal'


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    old=json.loads((o.OUT/'derivatives/plan.json').read_text())
    assert g.c.sha(o.OUT/'before-internal-native_opacity.py')==old['inputs_sha256']['verification/native_opacity.py']
    assert not json.loads((o.OUT/'derivatives/result.json').read_text())['passed']
    sources=['kap/public/kap_def.f90','kap/private/kap_eval_fixed.f90',
        'kap/private/kap_eval.f90','kap/public/kap_lib.f90','kap/private/kap_eval_support.f90']
    bindings={}
    for rel in sources:
        path=g.c.fresh.MESA/rel;target=OUT/'sources'/rel;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(path,target);bindings[rel]=g.c.sha(target)
    binary=g.c.fresh.BINARY
    disassembly=subprocess.check_output(['objdump','-d','--start-address=0x7d42c0','--stop-address=0x7d43d0',str(binary)],text=True)
    assert 'cvtsd2ss' in disassembly
    (OUT/'native-input-conversions.txt').write_text(disassembly)
    (OUT/'native-table-layout.txt').write_text(subprocess.check_output([
        'objdump','-d','--start-address=0x7e6cb0','--stop-address=0x7e6cf2',str(binary)],text=True))
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='90a5106',
        input_sha256={str(p.relative_to(g.ROOT)):g.c.sha(p) for p in [
            o.OUT/'baseline-type1-captured.npz',o.OUT/'derivatives/result.json',
            g.ROOT/'verification/native_opacity.py',g.ROOT/'verification/opacity_internal.py']},
        source_sha256=bindings,executable_sha256=g.c.sha(binary),
        purpose='Read actual inputs to the radiative/conductive combiner, its radiative component and runtime blend boundaries. Reproduce the original outputs bitwise before inferring which branch was active.',
        internal_columns=['rho','log10rho','T','log10T','zbar','lnfree_e','free_e_rho','free_e_T','kap_rad','kap_rad_rho','kap_rad_T'],
        physical_opacity_certified=False))
    blend_identity()


def blend_identity():
    t,lo,hi=sp.symbols('t lo hi',real=True);a,b=sp.Function('a')(t),sp.Function('b')(t)
    w=(t-lo)/(hi-lo);k=w*a+(1-w)*b
    returned=(w*sp.diff(a,t)+(1-w)*sp.diff(b,t))/k
    missing=(a-b)/((hi-lo)*k)
    assert sp.simplify(sp.diff(sp.log(k),t)-returned-missing)==0
    # Positive constant components expose the omitted weight derivative.
    assert sp.simplify(missing.subs({a:2,b:1,lo:0,hi:1}))==1/(t+1)
    save('blend-identity.json',dict(classification='Proven',passed=True,
        assumptions='Positive differentiable component opacities, interior of a blend interval, t=log10 T, and fixed blend endpoints lo<hi.',
        identity='For k=w*a+(1-w)*b and w=(t-lo)/(hi-lo), d ln k/d ln T equals the component-weighted logarithmic derivatives plus (a-b)/(ln(10)*(hi-lo)*k).',
        source_boundary='The saved low/high-temperature and Compton blend formulas omit this weight-derivative term. This is an algebraic statement about the displayed smooth formulas, not a claim that a current cell entered a blend region.',
        positive_control='Constant a=2,b=1,lo=0,hi=1 gives nonzero missing derivative 1/(ln(10)*(1+t)); equal constants give zero.'))


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['input_sha256'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    label='internal-baseline';o.setup(label);o.trace(label,inspect_internal=True)
    original=dict(np.load(o.OUT/'baseline-type1-captured.npz'));actual=dict(np.load(o.OUT/(label+'-captured.npz')))
    equal={key:bool(np.array_equal(original[key],actual[key])) for key in original};assert all(equal.values()),equal
    row=json.loads((o.OUT/(label+'-trace.json')).read_text());bounds=row['global_controls'];a=actual['inner']
    low=max(bounds['kap_blend_logt_lower_bdy'],bounds['kap_z_tables']['logT_min'])
    high=min(bounds['kap_blend_logt_upper_bdy'],bounds['kap_lowt_z_tables']['logT_max'])
    compton_low=bounds['compton_blend_hi']-.5
    lower_blend=(a[:,3]>low)&(a[:,3]<high);compton=a[:,3]>compton_low
    rho_changes=a[:,1]-original['parameters'][:,3];temp_changes=a[:,3]-original['parameters'][:,4]
    rounded=np.float32(original['parameters'][:,3:5]).astype(float)
    save('result.json',dict(classification='Counterexample candidate',passed=True,bitwise_equal=equal,
        actual_logT_range=[float(a[:,3].min()),float(a[:,3].max())],
        lowT_blend=[low,high],lowT_blend_cells=np.flatnonzero(lower_blend).tolist(),
        compton_blend_low=compton_low,compton_active_cells=np.flatnonzero(compton).tolist(),
        input_rounding_only_rows=int(np.count_nonzero(np.all(a[:,[1,3]]==rounded,axis=1))),
        altered_log_rho_maximum=float(abs(rho_changes).max()),altered_log_T_maximum=float(abs(temp_changes).max()),
        source_of_all_derivative_failures_identified=False,physical_or_continuous_certificate=False))
    print('OPACITY INTERNAL',len(a),'identical outputs;',int(lower_blend.sum()),'lowT blend;',int(compton.sum()),'Compton',flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
