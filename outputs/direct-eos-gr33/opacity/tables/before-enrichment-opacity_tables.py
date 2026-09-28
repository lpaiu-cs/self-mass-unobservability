"""Double-precision evaluation of the actual saved opacity spline coefficients.

Counterexample candidate. Original data/composition assumptions remain. This
repairs arithmetic and the derivative of the declared clipped/interpolated map;
it is not physical opacity calibration or an everywhere-C1 model.
"""
import json, shutil, sys
import numpy as np
import sympy as sp
import opacity_cubic as c

o=c.o;g=c.g;OUT=o.OUT/'tables'


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    before=json.loads((c.OUT/'propagation/plan.json').read_text())
    assert g.c.sha(o.OUT/'before-tables-native_opacity.py')==before['inputs_sha256']['verification/native_opacity.py']
    for rel in ['kap/private/kap_eval_support.f90','kap/private/kap_eval_fixed.f90','kap/private/condint.f90','kap/public/kap_def.f90']:
        target=OUT/'sources'/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(g.c.fresh.MESA/rel,target)
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='2438b36',
        inputs_sha256={str(p.relative_to(g.ROOT)):g.c.sha(p) for p in [
            o.OUT/'baseline-type1-captured.npz',g.ROOT/'verification/native_opacity.py',g.ROOT/'verification/opacity_tables.py']},
        method='Read the loaded native radiative spline coefficients and conduction coefficients from a disposable child, with unchanged stellar outputs. Retain the saved binary32 data exactly; evaluate the nominal bicubic polynomials in binary64, differentiate the selected X cubic, retain linear Z interpolation and the harmonic radiative/conductive combination.',
        boundary='Restricted to the actually loaded Type1 high-temperature tables below Compton blending, logR>-19.5 and logT>3.5. Reject missing tables. Differentiate coordinate clipping explicitly; do not claim differentiability on a clipping or limiter boundary.',
        value_relative_tolerance=1e-4,finite_log_steps=[5e-5,2.5e-5],finite_derivative_relative_tolerance=1e-3,
        model_change='Remove intermediate binary32 rounding and the artificial 1e-7 upper-bound snap; use nominal exact 1/6 and 1/36 spline coefficients. Preserve all original failed derivative verdicts.',
        physical_opacity_certified=False,continuous_full_derivatives_certified=False))
    polynomial_control()


def capture():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['inputs_sha256'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    label='tables-baseline';o.setup(label);o.trace(label,inspect_internal=True,dump_tables=True)
    base=dict(np.load(o.OUT/'internal-baseline-captured.npz'));actual=dict(np.load(o.OUT/(label+'-captured.npz')))
    equal={key:bool(np.array_equal(value,actual[key])) for key,value in base.items()};assert all(equal.values()),equal
    save('capture-control.json',dict(classification='Counterexample candidate',passed=True,bitwise_equal=equal,
        coefficients_sha256=g.c.sha(o.OUT/(label+'-tables.npz'))))


def spline(gridx,gridy,coeff,x,y):
    """Nominal PSPLINE bicubic, with the derivative of coordinate clipping."""
    clipped_x=np.clip(x,gridx[0],gridx[-1]);clipped_y=np.clip(y,gridy[0],gridy[-1])
    i=np.clip(np.searchsorted(gridx,clipped_x,side='right')-1,0,len(gridx)-2)
    j=np.clip(np.searchsorted(gridy,clipped_y,side='right')-1,0,len(gridy)-2)
    def basis(grid,k,z):
        h=float(grid[k+1]-grid[k]);t=(z-grid[k])/h;s=1-t
        return np.array([[s,t],[h*h*(s*s*s-s)/6,h*h*(t*t*t-t)/6]]),np.array([[-1/h,1/h],[h*(1-3*s*s)/6,h*(3*t*t-1)/6]])
    bx,dx=basis(gridx,i,clipped_x);by,dy=basis(gridy,j,clipped_y);answer=np.zeros(3)
    for a,b,k in [(0,0,0),(1,0,1),(0,1,2),(1,1,3)]:
        nodes=coeff[k,i:i+2,j:j+2]
        answer += [bx[a]@nodes@by[b],dx[a]@nodes@by[b],bx[a]@nodes@dy[b]]
    answer[1]*=float(x==clipped_x);answer[2]*=float(y==clipped_y)
    return answer


def polynomial_control():
    x,y=sp.symbols('x y');f=x**3*y**2+2*x**2*y+4*x-3*y+7
    fields=[f,sp.diff(f,x,2),sp.diff(f,y,2),sp.diff(f,x,2,y,2)]
    grid=np.array([0.,1.]);coeff=np.array([[[float(a.subs({x:float(i),y:float(j)})) for j in grid] for i in grid] for a in fields])
    expected=sp.lambdify((x,y),[f,sp.diff(f,x),sp.diff(f,y)],'numpy');errors=[]
    for a in np.linspace(0,1,9):
        for b in np.linspace(0,1,9): errors.append(float(abs(spline(grid,grid,coeff,a,b)-expected(a,b)).max()))
    assert max(errors)<2e-13,max(errors)
    clipped=spline(grid,grid,coeff,-.1,.5);assert clipped[1]==0
    save('polynomial-control.json',dict(classification='Proven',passed=True,
        declared_polynomial=str(f),points=81,maximum_binary64_error=max(errors),
        scope='Manufactured bicubic reproduction and a constant clipped-coordinate control; not a physical data certificate.'))


class Opacity:
    def __init__(self):
        control=json.loads((OUT/'capture-control.json').read_text());assert control['passed']
        path=o.OUT/'tables-baseline-tables.npz';assert g.c.sha(path)==control['coefficients_sha256']
        self.data={k:a.astype(float) for k,a in dict(np.load(path)).items()}
        self.controls=json.loads((o.OUT/'tables-baseline-trace.json').read_text())['global_controls']
        assert self.controls['choices']==dict(cubic_interpolation_in_X=1,cubic_interpolation_in_Z=0,include_electron_conduction=1)
        assert self.data['clip_to_kap_table_boundaries']==1
        self.zgrid=np.array(self.controls['kap_z_tables']['Z'])

    def table(self,zi,xi,r,t):
        stem=f'kap_z_tables-{zi}-{xi}'
        assert stem+'-coeff' in self.data,('Unloaded opacity table',zi,xi)
        a=spline(self.data[stem+'-R'],self.data[stem+'-T'],self.data[stem+'-coeff'],r-3*t+18,t)
        return np.array([a[0],a[1],a[2]-3*a[1]])

    def hydrogen(self,zi,X,r,t):
        grid=self.data[f'kap_z_tables-{zi}-Xgrid'];n=len(grid)
        if X<=grid[0]: return self.table(zi,0,r,t)
        if X>=grid[-1]: return self.table(zi,n-1,r,t)
        j=int(np.searchsorted(grid,X,side='right'))-1
        if X==grid[j]: return self.table(zi,j,r,t)
        if n<4:
            u=(X-grid[j])/(grid[j+1]-grid[j]);return (1-u)*self.table(zi,j,r,t)+u*self.table(zi,j+1,r,t)
        start=np.clip(j-1,0,n-4);indices=range(start,start+4)
        values=np.array([self.table(zi,k,r,t) for k in indices]);w,_=c.weights(grid[start:start+4],values[:,0],X)
        return w@values

    def __call__(self,p):
        zbar,X,Z,r,t=p[:5]
        assert t>max(3.5,self.controls['kap_blend_logt_upper_bdy'])
        assert t<self.controls['compton_blend_hi']-.5 and r-3*t+18>-19.5
        grid=self.zgrid
        if Z<=grid[0]: rad=self.hydrogen(0,X,r,t)
        elif Z>=grid[-1]: rad=self.hydrogen(len(grid)-1,X,r,t)
        else:
            j=int(np.searchsorted(grid,Z,side='right'))-1;u=(Z-grid[j])/(grid[j+1]-grid[j])
            rad=(1-u)*self.hydrogen(j,X,r,t)+u*self.hydrogen(j+1,X,r,t)
        cg=self.data['conduction-logzs'];z=np.clip(np.log10(zbar),cg[0],cg[-1])
        j=np.clip(np.searchsorted(cg,z,side='right')-1,0,len(cg)-2);u=(z-cg[j])/(cg[j+1]-cg[j])
        cr=self.data['conduction-logrhos'];ct=self.data['conduction-logts'];coeff=self.data['conduction-f_ary']
        cond=(1-u)*spline(cr,ct,coeff[:,:,:,j],r,t)+u*spline(cr,ct,coeff[:,:,:,j+1],r,t)
        rc=np.clip(r,cr[0],cr[-1]);tc=np.clip(t,ct[0],ct[-1])
        cond=np.array([3*tc-rc-cond[0]-3.51937938116756,(-1-cond[1])*float(rc==r),(3-cond[2])*float(tc==t)])
        kr=10**rad[0];kc=10**cond[0];k=kr*kc/(kr+kc)
        result=np.r_[k,(k/kr)*rad[1:]+(k/kc)*cond[1:]];assert np.all(np.isfinite(result)) and k>0
        return result


def run():
    plan=json.loads((OUT/'plan.json').read_text());source=dict(np.load(o.OUT/'baseline-type1-captured.npz'))
    model=Opacity();values=[];errors=[];mutual=[]
    for i,p in enumerate(source['parameters']):
        a=model(p);values.append(a);directions=[]
        for j in [0,1]:
            estimates=[]
            for h in plan['finite_log_steps']:
                pair=[]
                for sign in [-1,1]:
                    shifted=p.copy();shifted[3+j]+=sign*h/np.log(10);pair.append(model(shifted)[0])
                estimates.append(np.log(pair[1]/pair[0])/(2*h))
            directions.append(estimates)
        directions=np.array(directions);errors.append(abs(directions-a[1:,None])/np.maximum(1,abs(a[1:,None])))
        mutual.append(abs(directions[:,0]-directions[:,1])/np.maximum(1,abs(directions[:,1])))
        if i%1000==0: print('DOUBLE OPACITY',i,flush=True)
    values=np.array(values);errors=np.array(errors);mutual=np.array(mutual)
    relative=abs(values[:,0]/source['outputs'][:,0]-1)
    passed=bool(relative.max()<plan['value_relative_tolerance'] and errors.max()<plan['finite_derivative_relative_tolerance'] and mutual.max()<plan['finite_derivative_relative_tolerance'])
    np.savez_compressed(OUT/'evaluation.npz',values=values,derivative_scores=errors,mutual_scores=mutual,value_relative_errors=relative)
    save('result.json',dict(classification='Counterexample candidate',passed=passed,cells=len(values),
        maximum_value_relative_change=float(relative.max()),maximum_derivative_score=float(errors.max()),
        derivative_worst_index=[int(a) for a in np.unravel_index(np.argmax(errors),errors.shape)],
        maximum_two_step_difference=float(mutual.max()),
        physical_opacity_certified=False,continuous_full_derivatives_certified=False))
    print('DOUBLE OPACITY COMPLETE',passed,float(relative.max()),float(errors.max()),float(mutual.max()),flush=True)
    assert passed,'Preserved failed double opacity gate'


def scalar_compatibility():
    plan=json.loads((OUT/'plan.json').read_text())
    assert g.c.sha(OUT/'before-scalar-compat-opacity_tables.py')==plan['inputs_sha256']['verification/opacity_tables.py']
    plan['inputs_sha256']['verification/opacity_tables.py']=g.c.sha(g.ROOT/'verification/opacity_tables.py')
    save('plan.json',plan)
    save('scalar-compatibility.json',dict(classification='Proven',
        preflight_failure="SymPy/mpmath could not convert the NumPy scalar string 'np.float64(1.0)'",
        repair='Convert manufactured-control substitution values to Python float and JSON indices to int. No scientific tolerance or model formula changed.',
        before_plan_sha256=g.c.sha(OUT/'before-scalar-compat-plan.json')))
    polynomial_control()


if __name__=='__main__': globals()[sys.argv[1]]()
