"""Exact rational selector proof on saved affine opacity-ordinate paths.

Proven, conditional: this encloses the four-ordinate interpolation map, not
the unknown curvature or error of the underlying physical opacity tables.
"""
from fractions import Fraction as F
import json, sys
import numpy as np
import opacity_cubic as c

OUT=c.OUT/'certificate'


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def affine_interpolant(x,y,v,at,radius):
    x=np.array([F(float(a)) for a in x],dtype=object);at=F(float(at));h=np.diff(x)
    values=np.array([[F(float(a)),F(float(b))] for a,b in zip(y,v)],dtype=object)
    assert len(x)==4 and all(a>0 for a in h) and x[0]<=at<=x[-1] and radius>0
    def bounds(a): return min(a[0]-radius*a[1],a[0]+radius*a[1]),max(a[0]-radius*a[1],a[0]+radius*a[1])
    def sign(a):
        lo,hi=bounds(a)
        if lo>0: return 1
        if hi<0: return -1
        if lo==hi==0: return 0
        raise ValueError('A sign may change on the declared affine path')
    def minimum(candidates):
        j=min(range(len(candidates)),key=lambda j:candidates[j][0])
        assert all(bounds(a-candidates[j])[0]>=0 for a in candidates),'A minimum selector may change'
        return candidates[j]
    s=(values[1:]-values[:-1])/h[:,None];m=np.array([[F(0),F(0)] for _ in x],dtype=object)
    for i in [1,2]:
        p=(s[i-1]*h[i]+s[i]*h[i-1])/(h[i-1]+h[i])
        a,b,d=sign(s[i-1]),sign(s[i]),sign(p)
        if a+b==0: continue
        m[i]=(a+b)*minimum([a*s[i-1],b*s[i],F(1,2)*d*p])
    for i,j,k in [(0,0,1),(3,2,1)]:
        p=s[j]*(1+h[j]/(h[j]+h[k]))-s[k]*h[j]/(h[j]+h[k]);a,b=sign(p),sign(s[j])
        if a*b<=0: continue
        difference=a*p-2*b*s[j];lo,hi=bounds(difference)
        if lo>0: m[i]=2*s[j]
        elif hi<=0: m[i]=p
        else: raise ValueError('An endpoint limiter may change on the declared path')
    j=min(sum(a<=at for a in x)-1,2);dx=at-x[j]
    q=(3*s[j]-2*m[j]-m[j+1])/h[j];r=(m[j]+m[j+1]-2*s[j])/(h[j]*h[j])
    return values[j]+dx*(m[j]+dx*(q+dx*r))


def prepare():
    assert not OUT.exists();OUT.mkdir()
    assert json.loads((c.OUT/'result.json').read_text())['passed']
    save('plan.json',dict(classification='Proven',
        source_sha256=c.g.c.sha(c.g.ROOT/'verification/opacity_cubic_certificate.py'),
        input_sha256=c.g.c.sha(c.o.OUT/'cubic-baseline-cubic.npz'),
        affine_parameter_radius='0.00001',
        interpretation='The saved binary64 x, at, value and derivative ordinates are exact rational inputs. At fixed composition coordinate, use y(t)=y0+t*v for each saved density and temperature direction and |t|<=1e-5. Verify every sign/minimum/endpoint limiter over the whole interval by exact endpoint inequalities of affine expressions.',
        claim='For each accepted path the interpolation value is an exact affine function a+b*t throughout the interval, with exact derivative b on its interior. Compare b with the saved native interpolation of derivative ordinates.',
        exclusions='This certificate does not establish that physical or tabulated opacity ordinates follow these affine paths. It does not certify table data, native floating-point evaluation, full opacity derivatives, EOS, GR evolution or observations.'))


def run():
    plan=json.loads((OUT/'plan.json').read_text());path=c.o.OUT/'cubic-baseline-cubic.npz'
    assert c.g.c.sha(path)==plan['input_sha256']
    assert c.g.c.sha(c.g.ROOT/'verification/opacity_cubic_certificate.py')==plan['source_sha256']
    data=dict(np.load(path));radius=F(plan['affine_parameter_radius']);rows=[];unresolved=[];worst=None
    for i in range(0,len(data['cell']),3):
        for direction in [1,2]:
            try:
                affine=affine_interpolant(data['x'][i],data['values'][i],data['values'][i+direction],data['at'][i],radius)
            except (AssertionError,ValueError) as exc:
                unresolved.append(dict(call=i,direction=direction,reason=str(exc)));continue
            wrong=F(float(data['result'][i+direction]));error=affine[1]-wrong
            row=dict(call=i,cell=int(data['cell'][i]),direction=direction,
                exact_value_at_zero=str(affine[0]),exact_derivative=str(affine[1]),
                exact_native_derivative_difference=str(error))
            rows.append(row)
            if worst is None or abs(error)>abs(F(worst['exact_native_derivative_difference'])): worst=row
        if i%3000==0: print('EXACT CUBIC PATHS',i,len(rows),len(unresolved),flush=True)
    assert rows and worst and F(worst['exact_native_derivative_difference'])!=0
    # Independent affine reproduction and a nonlinear-interpolator counterexample.
    x=[0.,1.,2.,3.];y=[1.,2.,3.,4.];v=[2.,2.,2.,2.]
    a=affine_interpolant(x,y,v,1.25,radius);assert a[0]==F(9,4) and a[1]==2
    save('result.json',dict(classification='Proven',completed=True,certified_paths=len(rows),
        unresolved_paths=unresolved,all_paths_certified=len(unresolved)==0,
        worst=worst,rows=rows,physical_or_full_derivative_certificate=False))
    print('EXACT CUBIC CERTIFICATES',len(rows),'unresolved',len(unresolved),'worst difference',float(F(worst['exact_native_derivative_difference'])),flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
