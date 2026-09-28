"""Differentiate the selected monotone cubic, not a new cubic of derivatives.

Proven: on an open fixed-selector region this interpolator is linear in its
four ordinates. Counterexample candidate: native calls and finite controls.
"""
import json, shutil, sys
import numpy as np
import native_opacity as o

g=o.g;OUT=o.OUT/'cubic'


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def weights(x,y,at):
    x=np.asarray(x);y=np.asarray(y);h=np.diff(x);identity=np.eye(4)
    assert x.shape==y.shape==(4,) and np.all(h>0) and x[0]<=at<=x[-1]
    slope=(identity[1:]-identity[:-1])/h[:,None];s=slope@y
    m=np.zeros((4,4));stable=True
    for i in [1,2]:
        p=(slope[i-1]*h[i]+slope[i]*h[i-1])/(h[i-1]+h[i])
        candidates=np.array([slope[i-1],slope[i],.5*p]);v=candidates@y
        factor=np.copysign(1,s[i-1])+np.copysign(1,s[i])
        stable &= bool(s[i-1]!=0 and s[i]!=0)
        if factor!=0:
            selected=int(np.argmin(abs(v)));ordered=np.sort(abs(v))
            stable &= bool(ordered[0]!=ordered[1])
            m[i]=factor*np.copysign(1,v[selected])*candidates[selected]
    for i,j,k in [(0,0,1),(3,2,1)]:
        p=slope[j]*(1+h[j]/(h[j]+h[k]))-slope[k]*h[j]/(h[j]+h[k]);value=p@y
        stable &= bool(value*s[j]!=0 and abs(value)!=2*abs(s[j]))
        if value*s[j]<=0: continue
        m[i]=2*slope[j] if abs(value)>2*abs(s[j]) else p
    j=min(int(np.searchsorted(x,at,side='right'))-1,2);dx=at-x[j]
    quadratic=(3*slope[j]-2*m[j]-m[j+1])/h[j]
    cubic=(m[j]+m[j+1]-2*slope[j])/(h[j]*h[j])
    return identity[j]+dx*(m[j]+dx*(quadratic+dx*cubic)),stable


def prepare():
    assert not OUT.exists();OUT.mkdir()
    assert not json.loads((o.OUT/'derivatives/result.json').read_text())['passed']
    older=json.loads((o.OUT/'internal/plan.json').read_text())
    assert g.c.sha(o.OUT/'before-cubic-native_opacity.py')==older['input_sha256']['verification/native_opacity.py']
    sources=['interp_1d/private/interp_1d_pm.f90','interp_1d/public/interp_1d_lib.f90']
    for rel in sources:
        target=OUT/'sources'/rel;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(g.c.fresh.MESA/rel,target)
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='27a404e',
        inputs_sha256={str(p.relative_to(g.ROOT)):g.c.sha(p) for p in [
            o.OUT/'baseline-type1-captured.npz',o.OUT/'internal/result.json',
            g.ROOT/'verification/native_opacity.py',g.ROOT/'verification/opacity_cubic.py']},
        value_relative_tolerance=2e-13,finite_log_steps=[1e-5,5e-6],finite_derivative_relative_tolerance=1e-7,
        source='The actual native interpolate_vector calls with n_old=4 and interp_pm in the Type1 opacity call. Read both inputs and outputs; preserve every result.',
        formula='Within an open selector region, I(y)=w(y) dot y and D I(y)[v]=w(y) dot v. I(v) generally selects a different monotone interpolant and is not that derivative.',
        boundary='Selector ties are explicitly marked; a continuous full opacity/EOS/trajectory certificate is not inferred.',
        physical_or_continuous_certificate=False))


def capture():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['inputs_sha256'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    label='cubic-baseline';o.setup(label);o.trace(label,inspect_internal=True,capture_cubic=True)
    baseline=dict(np.load(o.OUT/'baseline-type1-captured.npz'));actual=dict(np.load(o.OUT/(label+'-captured.npz')))
    equal={k:bool(np.array_equal(v,actual[k])) for k,v in baseline.items()};assert all(equal.values()),equal
    save('capture-control.json',dict(classification='Counterexample candidate',passed=True,bitwise_equal=equal,
        captured_sha256=g.c.sha(o.OUT/(label+'-cubic.npz'))))


def audit():
    assert json.loads((OUT/'capture-control.json').read_text())['passed']
    plan=json.loads((OUT/'plan.json').read_text());data=dict(np.load(o.OUT/'cubic-baseline-cubic.npz'))
    n=len(data['cell']);assert n%3==0;rows=[];max_value=0.;max_finite=0.;stable_count=0
    for i in range(n):
        w,stable=weights(data['x'][i],data['values'][i],data['at'][i]);value=w@data['values'][i]
        score=abs(value-data['result'][i])/max(1,abs(data['result'][i]));max_value=max(max_value,score)
        assert score<plan['value_relative_tolerance'],(i,score)
    for i in range(0,n,3):
        assert np.all(data['cell'][i:i+3]==data['cell'][i]) and np.all(data['axis'][i:i+3]==data['axis'][i])
        assert np.all(data['x'][i:i+3]==data['x'][i]) and np.all(data['at'][i:i+3]==data['at'][i])
        x=data['x'][i];y=data['values'][i];at=data['at'][i];w,stable=weights(x,y,at)
        corrected=[];original=[];finite_errors=[];path_stable=stable
        for direction in [1,2]:
            v=data['values'][i+direction];corrected.append(float(w@v));original.append(float(data['result'][i+direction]))
            for h in plan['finite_log_steps']:
                endpoints=[]
                for sign in [-1,1]:
                    shifted=y+sign*h*v;wj,sj=weights(x,shifted,at);endpoints.append(wj@shifted)
                    path_stable &= bool(sj and np.allclose(wj,w,rtol=0,atol=2e-14))
                finite=(endpoints[1]-endpoints[0])/(2*h)
                finite_errors.append(abs(finite-corrected[-1])/max(1,abs(corrected[-1])))
        if path_stable:
            stable_count+=1;max_finite=max(max_finite,max(finite_errors))
            assert max(finite_errors)<plan['finite_derivative_relative_tolerance'],(i,finite_errors)
        rows.append(dict(cell=int(data['cell'][i]),axis=str(data['axis'][i]),corrected=corrected,original=original,
            derivative_difference=(np.array(corrected)-original).tolist(),stable_tested_path=bool(path_stable)))
    assert stable_count>0
    worst=max(rows,key=lambda r:max(abs(np.array(r['derivative_difference']))))
    save('result.json',dict(classification='Counterexample candidate',passed=True,captured_calls=n,
        interpolation_triples=len(rows),maximum_native_value_score=max_value,
        stable_finite_paths=stable_count,maximum_stable_path_finite_score=max_finite,worst=worst,rows=rows,
        scope='The selected four-ordinate interpolator derivative is repaired and checked. This does not yet replace all native opacity table/rounding paths or certify physical opacity.'))
    print('CUBIC DERIVATIVE',n,'native values;',stable_count,'stable paths; worst',worst,flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
