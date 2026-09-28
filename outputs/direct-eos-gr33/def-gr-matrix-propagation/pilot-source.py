"""Same GR descriptor: exact source map and constrained matrix propagation.

Proven identity: u=g+x, g=K^-1 F, x-M x''=M g'', M=K^-1 D.
The correction lies in range(M), preserving the original algebraic rows.
Finite Krylov projection is a numerical candidate, never a changed EOS.
"""
from pathlib import Path
import argparse
import json
import resource
import signal
import time
import numpy as np
from scipy.sparse import coo_matrix, diags, vstack
from scipy.sparse.linalg import splu
from scipy.linalg import eig
import def_gr_laplace as prior

OUT=prior.OUT.parent/'def-gr-matrix-propagation'
write=prior.write


def source_map(problem):
    """Map the SAME full vector of integrated face energies to F+K H."""
    heat=problem.heat;bg=problem.bg;geom=problem.radiation.geometry;h=prior.task.h
    heat.d={k:heat.d[k] for k in heat.d.files}
    nf=len(heat.edges)
    def source(points):
        ids,f=heat.projection(points);n=len(ids);rows=np.repeat(np.arange(n),2)
        increment=coo_matrix((np.tile([1.,-1.],n),(rows,np.c_[ids,ids+1].ravel())),shape=(n,nf)).tocsr()
        enclosed=coo_matrix((np.c_[np.ones(n),f-1,-f].ravel(),(np.repeat(np.arange(n),3),np.c_[np.zeros(n,int),ids,ids+1].ravel())),shape=(n,nf)).tocsr()
        r=points['r'];N,a=geom.metric(r);A4=np.exp(-8*points['phi']**2)
        geo=h.gr.G*.1*geom.R**2/h.gr.C**4
        loss=diags(-geo/(heat.volumes[ids]*A4))@increment
        rho=np.interp(r,heat.rnative,heat.d['raw'][::-1,0])
        k=np.interp(r,heat.rnative,(heat.d['raw'][:,8]/heat.d['thermo'][:,5])[::-1])
        rr=diags(-k/(rho*geo))@loss
        J=diags(h.gr.G*1e-7/(h.gr.C**4*geom.R*N*a))@enclosed
        return rr,loss,J
    fn,_=prior.task.coupled.reactive.symbolic()
    ops,_,_,_=prior.task.coupled.operators(bg.mid,fn)
    mids=source(bg.mid);nodes=source(bg.nodes);n=len(bg.grid)-1;dx=np.diff(bg.grid)
    rowmaps=[]
    for j in range(4):
        rowmaps.append(sum(diags(dx*ops[:,j,k])@s for k,s in zip([4,5,8],mids)))
    jump=diags(bg.nodes['gamma'])@nodes[0]
    rowmaps[1]+=jump[1:]-jump[:-1]
    top=vstack(rowmaps).tocsr()[np.arange(4*n).reshape(4,n).T.ravel()]
    si=bg.surface_index;end={k:np.array([v[si]]) for k,v in bg.nodes.items()}
    with np.errstate(divide='ignore',invalid='ignore'):
        _,_,bracket,g=prior.task.coupled.operators(end,fn)
    zero=coo_matrix((1,nf)).tocsr()
    surface=sum(-bracket[k,0]/g[0]*a[si] for k,a in zip([4,5,8],nodes))
    forcing=vstack([top,zero,zero,surface,zero]).tocsc()
    r=bg.grid;ids=np.clip(np.searchsorted(heat.edges,r,side='right')-1,0,nf-2)
    f=(r-heat.edges[ids])/(heat.edges[ids+1]-heat.edges[ids]);n=len(r)
    interp=coo_matrix((np.c_[1-f,f].ravel(),(np.repeat(np.arange(n),2),np.c_[ids,ids+1].ravel())),shape=(n,nf)).tolil()
    factor=np.divide(1.,r**3,out=np.zeros_like(r),where=r>0)
    for i in np.flatnonzero(r<heat.edges[1]):interp.rows[i]=[1];interp.data[i]=[1.];factor[i]=heat.edges[1]**-3
    N,a=geom.metric(r);inside=(r<1)&(bg.nodes['p']>0)
    k=np.zeros(n);A4=np.exp(-8*bg.nodes['phi']**2)
    k[inside]=h.gr.G*1e-7*factor[inside]/(4*np.pi*h.gr.C**4*geom.R*N[inside]*a[inside]*A4[inside]*(bg.nodes['e'][inside]+bg.nodes['p'][inside]))
    hz=diags(k)@interp.tocsr();empty=coo_matrix(hz.shape).tocsr()
    lift=vstack([hz,empty,diags(r*bg.nodes['v'])@hz,empty]).tocsr()[np.arange(4*n).reshape(4,n).T.ravel()]
    total=(forcing+problem.K@lift).tocsc();total.eliminate_zeros()
    errors=[]
    for t in [1e-15,1e-8,1.]:
        energy=heat.faces(t)[1];actual=total@energy
        expected=problem.base(t)+problem.K@heat.lift(t,bg.nodes)
        errors.append(float(np.linalg.norm(actual-expected)/max(np.linalg.norm(expected),1e-100)))
    assert max(errors)<1e-12,errors
    return total,dict(classification='Counterexample candidate',map_relative_errors=errors,nnz=total.nnz)


class Propagator:
    def __init__(self,problem,source):
        self.problem=problem;self.L=source;self.K=problem.K;self.D=problem.D
        self.scale=np.asarray(abs(self.K).sum(1)).ravel()
        self.lu=splu(self.K.multiply((1/self.scale)[:,None]).tocsc())
        self.extended=self.K.astype(np.longdouble)
        heat=problem.heat;flux=np.zeros(len(heat.edges))
        flux[heat.face_ids]=np.sum(heat.amplitude,axis=1)*problem.radiation.geometry.tc
        self.slope=self.solve(self.L@flux)
        r=problem.bg.grid;self.advection=r*problem.bg.nodes['v']
        # Coordinate balancing only; the original K,D are used in every solve.
        y=self.slope.reshape(-1,4).copy();y[:,2]-=self.advection*y[:,0]
        self.units=np.maximum(np.max(abs(y),axis=0),1e-100)

    def solve(self,rhs,transpose=False):
        if transpose:
            return self.lu.solve(rhs,trans='T')/self.scale[:,None] if rhs.ndim==2 else self.lu.solve(rhs,trans='T')/self.scale
        z=self.lu.solve(rhs/self.scale[:,None] if rhs.ndim==2 else rhs/self.scale)
        for _ in range(2):
            defect=rhs.astype(np.longdouble)-self.extended@z.astype(np.longdouble)
            z+=self.lu.solve(np.asarray(defect/self.scale[:,None] if rhs.ndim==2 else defect/self.scale,float))
        return z

    def encode(self,x):
        y=x.reshape(-1,4).copy();y[:,2]-=self.advection*y[:,0]
        return (y/self.units).ravel()

    def decode(self,x):
        y=x.reshape(-1,4).copy()*self.units;y[:,2]+=self.advection*y[:,0]
        return y.ravel()

    def M(self,x):return self.encode(self.solve(self.D@self.decode(x)))

    def basis(self,size,shift_frequency=256.):
        sigma=shift_frequency**2
        A=(self.K-sigma*self.D).tocsc();scale=np.asarray(abs(A).sum(1)).ravel()
        lu=splu(A.multiply((1/scale)[:,None]).tocsc())
        def action(x):return self.encode(sigma*lu.solve((self.D@self.decode(x))/scale))
        q=self.M(self.encode(self.slope));q/=np.linalg.norm(q)
        Q=np.zeros((len(q),size));Q[:,0]=q
        for j in range(1,size):
            q=action(Q[:,j-1])
            for _ in range(2):q-=Q[:,:j]@(Q[:,:j].T@q)
            norm=np.linalg.norm(q);assert norm>1e-14,(j,norm)
            Q[:,j]=q/norm
        return Q


def prepare():
    assert not OUT.exists();OUT.mkdir()
    paths=[Path(__file__),Path(prior.__file__),Path(prior.task.__file__),Path(prior.task.coupled.__file__),
           prior.prior.OUT/'fine-bank.npz',prior.OUT/'result.json']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='f71d0558',
        claim='Keep the exact source and algebraic constraints by the static-plus-range(M) identity, then test source-driven rational Krylov propagation.',
        decision='First verify the full source-face map and measure one128-vector range basis. Require its projected finite spectrum to support a stable modal propagation before any production. Do not clip positive/complex modes or change K,D.',
        pilot=dict(basis_vectors=128,shift_frequency=256,hard_seconds=90,memory_GB=3,CPU_threads=1,new_EOS_calls=0),
        forecast='One sparse LU per shift plus128 static/range solves; cost unmeasured, stop at90s. No time, spatial or Fourier range increase.',
        bindings={str(p):prior.prior.digest(p) for p in paths}))
    signal.alarm(90);resource.setrlimit(resource.RLIMIT_AS,(int(3e9),int(3e9)));start=time.monotonic()
    p=prior.Problem(prior.prior.OUT/'fine-bank.npz');L,check=source_map(p)
    model=Propagator(p,L);setup=time.monotonic()-start;Q=model.basis(128)
    MQ=np.column_stack([model.M(q) for q in Q.T]);H=Q.T@MQ
    values=eig(H,right=False)
    np.savez_compressed(OUT/'pilot.npz',Q=Q,H=H,units=model.units,eigenvalues=values)
    result=dict(classification='Counterexample candidate',source_map=check,setup_seconds=setup,seconds=time.monotonic()-start,
        max_orthogonality_error=float(abs(Q.T@Q-np.eye(128)).max()),
        projected_positive_real_modes=int(np.sum(values.real>0)),projected_complex_modes=int(np.sum(abs(values.imag)>1e-10*abs(values))),
        eigenvalue_real_range=[float(values.real.min()),float(values.real.max())],units=model.units.tolist(),
        original_failure_resolved=False,full_dynamic_charge_solved=False)
    write(OUT/'pilot.json',result);signal.alarm(0);print('PILOT',json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare']);globals()[p.parse_args().action]()
