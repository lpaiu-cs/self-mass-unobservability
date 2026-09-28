"""Isolate the retained density-derivative failure by interaction component."""
import json,sys
import numpy as np
import gr_dense_plasma as d

g=d.g;OUT=g.OUT/'gr-dense-plasma-derivative-probe'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def parts(p,rs,ge,X):
    x=X/g.c.A;x=x/x.sum();z=g.c.Z
    electron=p.call('excor7',[rs,ge],7)[[0,1,2,4,3,5,6]]*(x@z)
    ii=np.zeros(7);ie=np.zeros(7)
    for j in np.flatnonzero(x>0):
        ii+=x[j]*d.seven(p.call('fition9',[float(ge*z[j]**(5/3))],6))
        ie+=x[j]*d.seven(p.call('fscrliq8',[rs,ge,float(z[j])],6))
    moments=[float(x@z),float(x@(z*z)),float(x@(z**2.5)),float(x@(z**(5/3))),float(x@(z*(z+1)**1.5))]
    mixing=d.seven(p.call('cormix',[rs,ge,*moments],6))
    return np.array([ii,ie,electron,mixing])


def run():
    assert not OUT.exists();OUT.mkdir();d.verify()
    paths=[g.ROOT/'verification/gr_dense_plasma_derivative_probe.py',d.OUT/'manifest.json']
    save('plan.json',dict(classification='Counterexample candidate',bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        cells=[2972,3043,4352,5734],steps=[2e-4,1e-4],
        purpose='Identify interaction component of the original failed PDR finite derivative; retain all signed residuals without setting a new pass threshold.'))
    state,_=d.state_data();a=np.load(d.OUT/'stellar-comparison.npz');p=d.Provider();rows=[]
    for i in [2972,3043,4352,5734]:
        k=int(np.flatnonzero(a['cells']==i)[0]);rs,ge=a['parameters'][k,:2];X=state['X'][i];base=parts(p,rs,ge,X)
        for h in [2e-4,1e-4]:
            minus=parts(p,rs*np.exp(h/3),ge*np.exp(-h/3),X)
            plus=parts(p,rs*np.exp(-h/3),ge*np.exp(h/3),X)
            estimate=base[:,2]+(plus[:,2]-minus[:,2])/(2*h);error=estimate-base[:,6]
            rows.append(dict(cell=i,step=h,parts=['ii','ie','ee','mix'],expected=base[:,6].tolist(),
                finite=estimate.tolist(),difference=error.tolist()))
    save('result.json',dict(classification='Counterexample candidate',rows=rows,physical_EOS_certified=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    print(json.dumps(rows,indent=2),flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
