"""Separate the native zero mixed term from the restored spline polynomial."""
from fractions import Fraction as F
import json, shutil, sys
import numpy as np
import opacity_tables as t

OUT=t.OUT/'mixed-term'


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    prior=json.loads((t.OUT/'enrichment-plan.json').read_text())
    assert t.g.c.sha(t.OUT/'before-mixed-term-opacity_tables.py')==prior['inputs_sha256']['verification/opacity_tables.py']
    source=t.g.c.fresh.MESA/'interp_2d/private/bicub_sg.f';text=source.read_text()
    start=text.index('      subroutine fvbicub(');body=text[start:]
    assert 'integer v,z36th,iadr,i,j' in body and 'z36th=sixth*sixth' in body
    shutil.copy2(source,OUT/'bicub_sg.f')
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='3f28e87',
        inputs_sha256={str(p.relative_to(t.g.ROOT)):t.g.c.sha(p) for p in [
            t.OUT/'result.json',t.OUT/'evaluation.npz',t.o.OUT/'tables-enriched-tables.npz',
            t.g.ROOT/'verification/opacity_tables.py',t.g.ROOT/'verification/opacity_mixed_term.py']},
        native_defect='The saved fvbicub declares z36th as INTEGER and assigns sixth*sixth. Its value is zero, so the mixed-curvature term is omitted from value and derivative outputs.',
        controls='Use the same double-precision evaluator with only the conduction mixed field set to zero for native-implementation matching. Replay the restored result exactly. Keep the original 1e-4 value and 1e-3 finite-derivative gates. Preserve the failed native-value match of the restored model.',
        physical_opacity_certified=False))
    # f=x^2*y^2 has a nonzero mixed term; native omission is visible at the midpoint.
    grid=np.array([0.,1.]);coeff=np.zeros((4,2,2))
    for i,x in enumerate(grid):
        for j,y in enumerate(grid): coeff[:,i,j]=[x*x*y*y,2*y*y,2*x*x,4]
    full=t.spline(grid,grid,coeff,.5,.5);coeff[3]=0;omitted=t.spline(grid,grid,coeff,.5,.5)
    assert abs(full[0]-1/16)<1e-15 and abs(omitted[0])<1e-15
    save('positive-control.json',dict(classification='Proven',passed=True,
        polynomial='x^2*y^2',point=[.5,.5],full_value=float(full[0]),omitted_value=float(omitted[0]),
        integer_assignment_result=int(np.float32(1/6)*np.float32(1/6))))


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['inputs_sha256'].items(): assert t.g.c.sha(t.g.ROOT/rel)==digest,rel
    original=dict(np.load(t.o.OUT/'baseline-type1-captured.npz'));restored=dict(np.load(t.OUT/'evaluation.npz'))
    legacy=t.Opacity(include_mixed=False);full=t.Opacity();values=[];scores=[];mutual=[]
    for i,p in enumerate(original['parameters']):
        assert np.array_equal(full(p),restored['values'][i]),i
        a=legacy(p);values.append(a);derivatives=[]
        for j in [0,1]:
            direction=[]
            for h in [5e-5,2.5e-5]:
                pair=[]
                for sign in [-1,1]:
                    shifted=p.copy();shifted[3+j]+=sign*h/np.log(10);pair.append(legacy(shifted)[0])
                direction.append(np.log(pair[1]/pair[0])/(2*h))
            derivatives.append(direction)
        derivatives=np.array(derivatives)
        scores.append(abs(derivatives-a[1:,None])/np.maximum(1,abs(a[1:,None])))
        mutual.append(abs(derivatives[:,0]-derivatives[:,1])/np.maximum(1,abs(derivatives[:,1])))
        if i%1000==0: print('NATIVE MIXED-OMISSION MODEL',i,flush=True)
    values=np.array(values);scores=np.array(scores);mutual=np.array(mutual)
    change=abs(values[:,0]/original['outputs'][:,0]-1);passed=bool(change.max()<1e-4 and scores.max()<1e-3 and mutual.max()<1e-3)
    np.savez_compressed(OUT/'native-compatible.npz',values=values,derivative_scores=scores,mutual_scores=mutual,native_value_relative_change=change)
    save('result.json',dict(classification='Counterexample candidate',passed=passed,cells=len(values),
        restored_model_replay_bitwise=True,maximum_native_value_relative_change=float(change.max()),
        maximum_derivative_score=float(scores.max()),maximum_two_step_difference=float(mutual.max()),
        restored_spline_maximum_relative_correction=float(abs(restored['values'][:,0]/values[:,0]-1).max()),
        original_restored_vs_native_value_gate_still_failed=True,
        physical_opacity_or_full_continuous_certificate=False))
    print('NATIVE OMISSION MATCH',passed,float(change.max()),float(scores.max()),float(mutual.max()),flush=True)
    assert passed,'Preserved failed native omission match'


def certificate():
    assert not (OUT/'continuous-correction-bound.json').exists()
    model=t.Opacity();data=model.data;rx=data['conduction-logrhos'];ty=data['conduction-logts'];coeff=data['conduction-f_ary'][3]
    maxima=[F(0),F(0),F(0)];locations=[None,None,None];count=0
    for i in range(len(rx)-1):
        hx=F(float(rx[i+1]))-F(float(rx[i]))
        for j in range(len(ty)-1):
            hy=F(float(ty[j+1]))-F(float(ty[j]))
            for z in range(coeff.shape[2]):
                total=sum((abs(F(float(v))) for v in coeff[i:i+2,j:j+2,z].flat),F(0))
                bounds=[hx*hx*hy*hy*total/225,hx*hy*hy*total/45,hx*hx*hy*total/45]
                for k,bound in enumerate(bounds):
                    if bound>maxima[k]: maxima[k]=bound;locations[k]=[i,j,z]
                count+=1
    assert count==29400
    save('continuous-correction-bound.json',dict(classification='Proven',completed=True,
        coefficient_sha256=t.g.c.sha(t.o.OUT/'tables-enriched-tables.npz'),rectangles=count,
        proof='For u in [0,1], abs(u^3-u)<=2/5 and abs(3*u^2-1)<=2. The omitted mixed term is hx^2*hy^2/36 times the sum of four mixed-curvature coefficients times qx*qy. The triangle inequality gives value <=hx^2*hy^2*sum_abs/225, x derivative <=hx*hy^2*sum_abs/45, and y derivative <=hx^2*hy*sum_abs/45. The extrema of u-u^3 are endpoints and u=1/sqrt(3); 2/(3*sqrt(3))<2/5 follows from 25<27. All coefficient and grid arithmetic here is exact rational arithmetic.',
        quantities=['delta log10 conductivity','d delta log10 conductivity/d log10 rho','d delta log10 conductivity/d log10 T'],
        exact_global_bounds=[str(a) for a in maxima],approximate_global_bounds=[float(a) for a in maxima],locations=locations,
        scope='Bounds over every closed saved conduction interpolation rectangle and every convex Z interpolation. They bound the declared missing polynomial term, not table-data uncertainty, total opacity error, EOS, or GR evolution.'))
    print('CONTINUOUS MIXED-TERM BOUND',count,[float(a) for a in maxima],flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
