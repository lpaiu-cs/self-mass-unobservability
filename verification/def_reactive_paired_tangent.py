"""Paired native baselines and the existing strict representable inverse."""
from pathlib import Path
from types import FunctionType
import inspect
import numpy as np
import def_reactive_thermal_tangent as old

OUT=old.OUT/'paired-native'
h=old.h


def initialize():
    global data,sources,eos,rest,background
    old.initialize();data,sources,eos,rest=old.data,old.sources,old.eos,old.rest
    background=np.load(h.OUT/'absolute-shoot/background-0.001.npz')['states']


def cell(i):
    X=data['X'][i];lapse=data['A'][i]*data['N'][i];T=data['lnT'][i]
    lp=float(background[i,2]);calls=0
    def query(mode,p,t,x):
        nonlocal calls
        calls+=1;return eos(mode,p,t,x)
    # Pair at the exact saved pressure coordinate: log(exp(lp)) and independent
    # native call rounding must not become a spurious reaction heat source.
    raw=query(1,lp,T,X);rho=raw[0];p=raw[1]
    cp=raw[10]-p/rho*raw[8]
    h0=np.longdouble(raw[2])+np.longdouble(p)/rho
    total0=np.longdouble((X/old.thermal.g.c.A)@old.thermal.g.c.W)*(h.gr.C*100)**2+np.longdouble(raw[2])
    rows=[]
    for dt in [8.,4.]:
        x=np.asarray(X.astype(np.longdouble)+np.longdouble(dt*lapse)*sources['dxdt'][i],float)
        assert x.min()>=0 and abs(x.sum()-1)<1e-12
        drest=rest.astype(np.longdouble)@(x.astype(np.longdouble)-X.astype(np.longdouble))
        loss=np.longdouble(dt*lapse)*(sources['neutrino'][i]+sources['thermal_neutrino'][i])
        target=h0-drest-loss
        budget=max(2.,32*abs(np.spacing(float(target))),float(abs(drest+loss))*1e-8)
        # The existing inverse supplies Newton plus adjacent-representable lnT
        # search. Keep the original energy budget; no tolerance relaxation.
        best=old.thermal.g.s.enthalpy_inverse(query,lp,x,float(target),float(T+(-drest-loss)/cp),budget)
        _,lt,a=best
        error=np.longdouble(a[2])+np.longdouble(a[1])/a[0]-target
        assert abs(error)<=budget,(i,dt,float(error),budget)
        lr=float(np.log(np.longdouble(a[0])/rho))
        de=np.longdouble(rho)*(total0*np.expm1(np.longdouble(lr))+np.exp(np.longdouble(lr))*(drest+np.longdouble(a[2])-raw[2]))
        rows.append([lr/dt,float(de/dt),float((np.longdouble(lt)-T)/dt),float(error),budget,float(drest),float(loss)])
    return i,np.array(rows),calls


ns=dict(vars(old),OUT=OUT,__file__=__file__,initialize=initialize,cell=cell,old=old)
ns['block']=FunctionType(old.block.__code__,ns)
text=inspect.getsource(old.run).replace('files=[Path(__file__),','files=[Path(__file__),Path(old.__file__),Path(old.thermal.g.s.__file__),')
text=text.replace("claim='Compute", "preserved_failure='Original pilot cell 4809 had enthalpy residual 70.59375 erg/g against unchanged 64 erg/g budget. Use exact saved pressure, paired native baseline, and the existing adjacent-float inverse.',\n        claim='Compute")
exec(compile(text,__file__,'exec'),ns)

if __name__=='__main__':ns['run']()
