"""Resolve retarded light-cone cuts without adding physical cells or steps."""
from pathlib import Path
from types import FunctionType
import json
import signal
import sys
import time
import numpy as np
from numpy.polynomial import legendre as leg
import def_native_anisotropic_gr as base

OUT=base.OUT/'characteristic';C=base.C;LD=base.LD;write=base.write;sha=base.sha


class Response(base.Response):
    def setup(self,d,order):
        super().setup(d,order)
        self.xfaces=C*self.geo(d['edges']-self.model.m.RJ)[1]
        self.mid=(self.xfaces[:-1]+self.xfaces[1:])/2
        self.half=np.diff(self.xfaces)/2
        xn=(self.x.reshape(-1,order)-self.mid[:,None])/self.half[:,None]
        self.inverse=np.linalg.inv(leg.legvander(xn,order-1))
        # A degree(order-1) source times a quadratic time antiderivative is
        # integrated exactly by this rule between its characteristic cuts.
        self.gx,self.gw=leg.leggauss((order+3)//2)

    def segments(self,t):
        owners=[];cells=[];left=[];right=[]
        for i,x in enumerate(self.tx):
            cuts=np.unique(np.r_[x,x-C*(t-self.t[self.t<=t]),x+C*(t-self.t[self.t<=t])])
            cuts=cuts[(cuts>self.xfaces[0])&(cuts<self.xfaces[-1])]
            ids=np.searchsorted(self.xfaces,cuts,side='right')-1
            for j in np.unique(ids):
                edges=np.r_[self.xfaces[j],cuts[ids==j],self.xfaces[j+1]]
                owners.extend([i]*(len(edges)-1));cells.extend([j]*(len(edges)-1))
                left.extend(edges[:-1]);right.extend(edges[1:])
        return tuple(np.asarray(x,dtype=int if k<2 else float) for k,x in enumerate([owners,cells,left,right]))

    def propagate(self,source):
        poly=base.flow.green.polynomial(self.t,source);H=poly.antiderivative()
        hc=H.c.reshape(3,len(self.t)-1,-1,self.order)/self.dx.reshape(-1,self.order)
        coeff=np.einsum('cij,ptcj->ptci',self.inverse,hc)
        U=[];Ut=[];Ux=[];cols=np.arange(len(self.ids))[None,:]
        for t in self.t:
            ret=t-self.distance;cut=np.clip(ret,0,self.t[-1]);idx=np.clip(np.searchsorted(self.t,cut,side='right')-1,0,len(self.t)-2);dt=cut-self.t[idx]
            hv=np.zeros_like(ret)
            for co in H.c:hv=hv*dt+co[idx,cols]
            sv=poly.c[0][idx,cols]*dt+poly.c[1][idx,cols]
            hv[ret<=0]=0.;sv[ret<=0]=0.
            hv=hv.reshape(len(self.tx),-1,self.order).sum(2)
            dv=(sv*self.sign).reshape(len(self.tx),-1,self.order).sum(2)
            sv=sv.reshape(len(self.tx),-1,self.order).sum(2)
            owner,cell,lo,hi=self.segments(t)
            if len(cell):
                xx=(lo[:,None]+hi[:,None])/2+(hi-lo)[:,None]*self.gx/2
                rr=t-abs(self.tx[owner,None]-xx)/C
                jj=np.clip(np.searchsorted(self.t,np.clip(rr,0,self.t[-1]),side='right')-1,0,len(self.t)-2)
                tt=np.clip(rr,0,self.t[-1])-self.t[jj]
                vv=leg.legvander((xx-self.mid[cell,None])/self.half[cell,None],self.order-1)
                cc=np.array([np.sum(co[jj,cell[:,None]]*vv,axis=-1) for co in coeff])
                hh=(cc[0]*tt+cc[1])*tt+cc[2];ss=2*cc[0]*tt+cc[1]
                hh[rr<=0]=0.;ss[rr<=0]=0.
                ww=(hi-lo)[:,None]*self.gw/2
                affected=np.unique(np.column_stack([owner,cell]),axis=0)
                ii,jj=affected.T;hv[ii,jj]=0.;sv[ii,jj]=0.;dv[ii,jj]=0.
                np.add.at(hv,(owner,cell),np.sum(ww*hh,axis=1))
                np.add.at(sv,(owner,cell),np.sum(ww*ss,axis=1))
                np.add.at(dv,(owner,cell),np.sum(ww*ss*np.sign(self.tx[owner,None]-xx),axis=1))
            U.append(C/2*hv.sum(1,dtype=LD));Ut.append(C/2*sv.sum(1,dtype=LD));Ux.append(-.5*dv.sum(1,dtype=LD))
        return np.asarray(U,float),np.asarray(Ut,float),np.asarray(Ux,float)

    run=FunctionType(base.Response.run.__code__,dict(vars(base),OUT=OUT))


def check():
    # Exact box source: H(s)=s^2/2 for S(t,x)=t on [-1,1].
    # Include targets within cells, and wave fronts inside cells; the old
    # uncut Gauss rule is not exact for either of these cases.
    errors=[]
    for order in [4,8]:
        m=Response.__new__(Response);m.order=order;m.t=np.array([0.,.3,.7,1.])/C
        m.xfaces=np.array([-1.,0.,1.]);m.mid=np.array([-.5,.5]);m.half=np.array([.5,.5]);gx,gw=leg.leggauss(order)
        m.x=(m.mid[:,None]+m.half[:,None]*gx).ravel();m.ids=np.repeat(np.arange(2),order)
        m.dx=(m.half[:,None]*np.broadcast_to(gw,(2,order))).ravel();m.tx=np.array([-.75,0.,.3,2.])
        m.inverse=np.linalg.inv(leg.legvander(np.broadcast_to(gx,(2,order)),order-1));m.gx,m.gw=leg.leggauss((order+3)//2)
        m.distance=abs(m.tx[:,None]-m.x[None,:])/C;m.sign=np.sign(m.tx[:,None]-m.x[None,:])
        u,ut,ux=m.propagate((m.t*C)[:,None]*m.dx)
        exact=[]
        for t in m.t*C:
            primitive=lambda y:np.sign(y)*(t**3-max(t-abs(y),0.)**3)/3
            exact.append([(primitive(1-x)-primitive(-1-x))/4 for x in m.tx])
        errors.append(float(np.max(abs(u-exact))))
    assert max(errors)<3e-15,errors
    return dict(passed=True,exact_box_source_absolute_errors=errors,
        classification='Proven',scope='Polynomial quadrature self-check only; not a continuum GR certificate.')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    old=json.loads((base.OUT/'fields.json').read_text());assert not old['passed']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        failure_preserved=str(base.OUT/'fields.json'),failure=old['controls'],
        cause_candidate='Unsplit absolute-distance cusp and retarded stored-time knots cross existing source cells. Resolve these integration cuts at the same4/8 interpolation order and same physical cells/times.',
        method='Interpolate source density in optical distance within each existing cell. Split only cells containing the target or a retarded time knot. Integrate each smooth polynomial piece exactly; keep existing quadrature elsewhere.',
        claim='Does correcting characteristic quadrature close the original0.002 field gate without increasing physical resolution?',
        gates=json.loads((base.OUT/'plan.json').read_text())['controls'],
        budget=dict(seconds=90,CPU_threads=1,memory_GB=3,new_fluid_steps=0,new_native_calls=0),
        measured_basis='Original three field paths took12.43s. Extra split cells are at most35 per target/time. Halt90s; no automatic order/grid enlargement.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(base.__file__),base.OUT/'fields.json',base.OUT/'source-64.npz',base.OUT/'source-128.npz']}))
    write(OUT/'check.json',check())
    for steps in [64,128]:(OUT/f'source-{steps}.npz').write_bytes((base.OUT/f'source-{steps}.npz').read_bytes())
    (OUT/'sources.json').write_bytes((base.OUT/'sources.json').read_bytes())


def fields():
    # Reuse the original evaluator and gates; only the quadrature changes.
    for p,h in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert sha(p)==h,p
    fn=FunctionType(base.fields.__code__,dict(vars(base),OUT=OUT,Response=Response))
    fn()


if __name__=='__main__':globals()[sys.argv[1]]()
