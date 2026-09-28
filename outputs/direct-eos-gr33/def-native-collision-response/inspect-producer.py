"""Counterexample candidate: compensated radiation/heat/H response.

Reuse the actual moving collision kernel and corrected native EOS. Density,
velocity and advected inventories are prescribed by the accepted trajectory;
their *additional* mechanical response remains open, with its impulse saved.
"""
from pathlib import Path
import json
import signal
import sys
import time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu, LinearOperator, gmres
import def_native_compensated_transport as prior

flow=prior.flow; C=prior.C; write=prior.write; sha=prior.sha
OUT=prior.metric.OUT.parent/'def-native-collision-response'
AMPLITUDE=1e-26


def prepare():
    assert not OUT.exists(); OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='e71ed7eba',
        claim='Evolve actual moving absorption/emission/Thomson, native thermal energy and neutral H response together with the GR-driven photon transport, using separately scaled increments on all531 cells.',
        decision='Determine whether thermochemical feedback amplifies or screens the geometric photon response; save the paired material energy/species and momentum impulse for the remaining mechanical and GR closure.',
        reuse='Reuse the corrected64/128 backgrounds,17 lapse/source times,531 cells,8 angles,152 frequencies and original3.434ms horizon. No new native bank or nonlinear trajectory.',
        boundary='Same fixed deep incoming packet perturbation and external vacuum. Density, velocity and other advected inventories follow the saved background. The additional force is measured, not silently treated as a free fluid evolution.',
        method='Strang composition of shared implicit transport and a coupled local linear radiation/heat/H solve. Actual moving scattering is retained. A static Thomson plus absorption/heat/H Schur inverse preconditions the local solve. Time comparison decides this partition.',
        metric='Use delta_F=delta_N/measure-delta_log_volume*F0. Include the volume/time-rate change and canonical density/isentropic temperature response in collisions. Do not put frequency-box losses or gravitational work into gas heat.',
        derivatives='Centered relative temperature/neutral/density coefficient probes in the existing native interpolation model; compare1e-5 and5e-6 at representative times. This is sampled derivative evidence, not a uniform EOS derivative bound.',
        gates=dict(derivative=.0001,owner=1e-10,linear=1e-10,energy=1e-8,species=1e-8,time=.02,background_time=.02),
        budget=dict(inspect_seconds=25,bank_seconds=90,pilot_seconds=45,production_seconds=240,CPU_threads=1,memory_GB=3,new_native_bank_calls=0),
        production='Only64/128 on background128 and128 on background64, after a measured pilot with2x projected cost fits240s. Halt on budget or gate failure, preserve outputs; no automatic extra mesh, paths, horizon or relaxed tolerance.',
        limitations='Prescribed material motion and finite EOS/angle/frequency model; no full hydrodynamic/GR feedback or final physical charge.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(prior.__file__),Path(flow.__file__),prior.OUT/'result.json',flow.OUT/'coupled-128.npz',flow.OUT/'coupled-64.npz',prior.metric.OUT/'corrected/metric-128-g8.npz',prior.metric.OUT/'corrected/metric-64-g8.npz']}))


class Response(prior.Transport):
    def __init__(self,reference):
        super().__init__(reference); self.data=dict(np.load(flow.OUT/f'coupled-{reference}.npz'))
        self.keep=np.array([np.argmin(abs(self.data['snapshot_t']-t)) for t in self.t])
        b=self.model.bulk; self.scale=b.scale.copy(); self.a=np.r_[b.d['a'],self.model.m.a]
        self.volume=4*np.pi*self.W*self.a**3; self.nb=b.n
        self.Eweight=self.weights*self.E*self.scale
        self.Nweight=self.weights*self.scale
        self.cache={}

    def restore(self,k):
        m=self.model; b=m.bulk; f=m.flow; j=self.keep[k]; d=self.data
        m.h=d['snapshot_h'][j].copy(); m.Pi=d['snapshot_Pi'][j].copy(); m.mass=d['snapshot_mass'][j].copy(); m.set_material(self.t[k])
        self.theta=d['snapshot_theta'][j].copy(); self.eta=d['snapshot_eta'][j].copy()
        self.rho,self.beta,self.lt,self.y=f.primitive(d['snapshot_U'][j]); self.active=self.rho>=f.eos.floor
        self.bulk_beta=m.velocity(); self.bulk_x=b.eos.x.copy()

    def coefficients(self,dr=0.,dt=0.,dy=0.):
        """Existing bound-free owners, with relative density/T/neutral probes."""
        m=self.model; b=m.bulk; f=m.flow; nb=self.nb
        b.eos.x=self.bulk_x+(1+self.bulk_x)*np.expm1(dr)
        theta=self.theta+dt; eta=(1+self.eta)*np.exp(dy)-1
        p,u,ut,uy,pt,py,ne,*_=b.eos.gas(theta,eta)
        ab,em,*_=b.eos.radiation(theta,eta)
        df=b.eos.spectral['frequency'][0]*(1-b.eos.f)+b.eos.spectral['frequency'][1]*b.eos.f
        boost=-self.bulk_beta[:,None,None]*self.mu[None,:,None]
        absorb=b.factor[:,None,None]*np.exp(dr)*ab[:,None,:]*(1+boost*(1+df[:,0,None,:]))
        emit=b.factor[:,None,None]*np.exp(dr)*em[:,None,:]*(1+boost*(1+df[:,1,None,:]))
        sc=b.scfactor*ne
        energy=m.mass*b.d['a']*u
        neutral=m.mass*b.d['thermo'][:,4]*b.d['y0']*(1+eta)
        rho=m.mass/b.volume*np.exp(dr); kap=ne*6.6524587321e-25/rho
        n=self.n; E=np.zeros((n,self.q,self.nf)); loss=np.zeros_like(E)
        E[:nb]=emit; loss[:nb]=absorb-emit
        ss=np.zeros(n); ss[:nb]=sc
        gg=np.zeros(n); gg[:nb]=energy; hh=np.zeros(n); hh[:nb]=neutral
        rr=np.zeros(n); rr[:nb]=rho; vv=np.r_[self.bulk_beta,self.beta]; kk=np.zeros(n); kk[:nb]=kap
        idx=np.flatnonzero(self.active); rho=self.rho[idx]*np.exp(dr); lt=self.lt[idx]+dt; y=self.y[idx]*np.exp(dy)
        f.eos.y=y
        if len(idx):
            pp,uu,_,_,kap,_,_=f.eos.evaluate(rho,lt); rr0=rho*f.eos.rho0; a=self.a[nb+idx]
            gamma=1/np.sqrt(1-self.beta[idx]**2); D=gamma[:,None]*(1-self.beta[idx,None]*self.mu)
            ab,em=m.spectrum.coefficients(rr0,lt,y,self.E[None,None,:]*D[:,:,None]/a[:,None,None])
            factor=a[:,None,None]*C*rr0[:,None,None]*f.eos.nH*D[:,:,None]
            E[nb+idx]=factor*em; loss[nb+idx]=factor*(ab-em)
            ss[nb+idx]=a*C*rr0*kap
            # Exact lab energy derivative at prescribed density/velocity;
            # its opposite photon transfer is evolved, momentum is retained.
            pressure=pp*f.eos.rho0*C*C; eps=rr0*(f.eos.cx+uu)*C*C
            gg[nb+idx]=a*self.volume[nb+idx]*((eps+pressure)*gamma**2-pressure)
            hh[nb+idx]=self.volume[nb+idx]*gamma*rr0*f.eos.nH*y
            rr[nb+idx]=rr0; kk[nb+idx]=kap
        b.eos.x=self.bulk_x
        return dict(emit=E,loss=loss,sc=ss,energy=gg,neutral=hh,rho=rr,beta=vv,kap=kk)

    def scattering_matrix(self,c):
        """Exact sparse form of the existing positive two-node packet map."""
        n,q,nf=self.n,self.q,self.nf; v=c['beta']; gamma=1/np.sqrt(1-v*v)
        edges=(self.edges_mu[None,:]-v[:,None])/(1-v[:,None]*self.edges_mu[None,:]); dw=np.diff(edges)/2
        p2=(edges[:,:-1]**2+edges[:,:-1]*edges[:,1:]+edges[:,1:]**2-1)/2
        rate=c['sc'][:,None]*gamma[:,None]*(1-v[:,None]*self.mu)
        ids=np.arange(n*q*nf).reshape(n,q,nf); rows=[ids.ravel()]; cols=[ids.ravel()]; vals=[np.broadcast_to(-rate[:,:,None],ids.shape).ravel()]
        escaped=np.zeros((3,n,q,nf)); packet=self.w[None,:,None]*self.num*self.scale
        for j,muin in enumerate(self.mu):
            for k,muout in enumerate(self.mu):
                probability=dw[:,k]*(1+.5*p2[:,j]*p2[:,k])
                dest=self.E[None,:]*(1-v[:,None]*muin)/(1-v[:,None]*muout)
                high=np.searchsorted(self.extended,dest); low=high-1
                assert high.min()>0 and high.max()<len(self.extended)
                frac=(dest-self.extended[low])/(self.extended[high]-self.extended[low])
                for ind,w in [(low,1-frac),(high,frac)]:
                    inside=(ind>0)&(ind<=nf); destid=np.clip(ind-1,0,nf-1)
                    gain=rate[:,j,None]*probability[:,None]*w
                    value=gain*self.w[j]/self.w[k]*(self.num*self.scale)[None,:]/(self.num*self.scale)[destid]
                    rows.append(ids[:,k][np.arange(n)[:,None],destid][inside]); cols.append(ids[:,j][inside]); vals.append(value[inside])
                    lost=gain*(~inside)*packet[:,j]
                    escaped[0,:,j]+=lost; escaped[1,:,j]+=lost*self.extended[ind]; escaped[2,:,j]+=lost*self.extended[ind]*muout
        S=sparse.coo_matrix((np.concatenate(vals),(np.concatenate(rows),np.concatenate(cols))),shape=(ids.size,ids.size)).tocsr(); S.eliminate_zeros()
        return S,escaped* (4*np.pi*self.W)[None,:,None,None]

    def inspect(self):
        started=time.monotonic(); rows=[]
        for k in [0,8,16]:
            self.restore(k); c=self.coefficients(); then=time.monotonic(); S,esc=self.scattering_matrix(c)
            I=self.I[k]; X=I/self.scale; actual=(S@X.ravel()).reshape(X.shape)*self.scale
            expected,escape=self.model.scattering(I,c['rho'],c['beta'],c['kap'],self.a)
            weight=self.weights*self.E
            error=float(np.sum(abs(actual-expected)*weight)/max(np.sum(abs(expected)*weight),1.))
            actualesc=np.einsum('knqf,nqf->kn',esc,X)
            target=escape.T*self.volume
            # Owner escape is per proper volume, with Eref units.
            err=float(np.max(abs(actualesc-target))/max(np.max(abs(target)),1.))
            rows.append(dict(time=float(self.t[k]),active=int(self.active.sum()),matrix_nnz=int(S.nnz),matrix_seconds=time.monotonic()-then,owner_relative=error,escape_relative=err,max_loss=float(c['loss'].max()),max_sc=float(c['sc'].max())))
        row=dict(classification='Counterexample candidate',seconds=time.monotonic()-started,rows=rows,passed=max(max(r['owner_relative'],r['escape_relative']) for r in rows)<1e-10)
        write(OUT/'inspect.json',row); print(json.dumps(row),flush=True); assert row['passed']


def inspect():
    signal.signal(signal.SIGALRM,flow.old.optical.timeout); signal.alarm(25)
    Response(128).inspect(); signal.alarm(0)


if __name__=='__main__': globals()[sys.argv[1]]()
