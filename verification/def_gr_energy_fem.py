"""Energy-consistent discretization of the proved continuous GR system.

Linear finite elements, two-point positive quadrature, physical pressure,
the same frozen heat poles and exact energy debit. This changes the old
midpoint spatial discretization, not the continuous equations or EOS.
"""
from pathlib import Path
import json
import signal
import time
import numpy as np
from scipy.sparse import coo_matrix,diags,vstack
from scipy.linalg import cholesky_banded
import def_gr_canonical as canonical

prior=canonical.prior
base=prior.prior
OUT=canonical.OUT;write=canonical.write


def source_points(heat,points):
    ids,f=heat.projection(points);n=len(ids);nf=len(heat.edges);geom=heat.geometry;h=base.task.h
    inc=coo_matrix((np.tile([1.,-1.],n),(np.repeat(np.arange(n),2),np.c_[ids,ids+1].ravel())),shape=(n,nf)).tocsr()
    enclosed=coo_matrix((np.c_[np.ones(n),f-1,-f].ravel(),(np.repeat(np.arange(n),3),np.c_[np.zeros(n,int),ids,ids+1].ravel())),shape=(n,nf)).tocsr()
    r=points['r'];N,a=geom.metric(r);A4=np.exp(-8*points['phi']**2);geo=h.gr.G*.1*geom.R**2/h.gr.C**4
    loss=diags(-geo/(heat.volumes[ids]*A4))@inc
    rho=np.interp(r,heat.rnative,heat.d['raw'][::-1,0]);ratio=np.interp(r,heat.rnative,(heat.d['raw'][:,8]/heat.d['thermo'][:,5])[::-1])
    rr=diags(-ratio/(rho*geo))@loss;J=diags(h.gr.G*1e-7/(h.gr.C**4*geom.R*N*a))@enclosed
    return rr,loss,J


def coefficients(bg,r):
    point=bg.sample(r);b=1-2*point['m']/r;C=point['N']*np.sqrt(b)
    n=len(r);result=np.zeros((n,36));inside=r<1
    if np.any(inside):
        rx=bg.saved['r']/bg.Rs;ids=np.clip(np.searchsorted(rx,r[inside],side='right')-1,0,len(rx)-2)
        ep=np.diff(bg.saved['e']*bg.Rs**2)[ids]/np.diff(rx)[ids]
        gp=np.diff(bg.saved['gamma'])[ids]/np.diff(rx)[ids]
        fn,_=canonical.symbolic()
        values=fn(*[point[k][inside] for k in ['r','m','p','e','phi','v','gamma']],np.exp(-8*point['phi'][inside]**2),C[inside],ep,gp)
        result[inside]=np.array([np.broadcast_to(v,ep.shape) for v in values]).T
    outside=~inside
    result[outside,7]=r[outside]**2*C[outside]          # B^-1 scalar
    result[outside,11]=-2*r[outside]**2*C[outside]*point['v'][outside]**2/b[outside]
    result[outside,15]=r[outside]**2/C[outside]       # W scalar
    result[outside,27]=-2*C[outside]*point['v'][outside]/b[outside]**2  # h_scalar,J
    assert np.all(np.isfinite(result))
    return result,point


def banded_spd(matrix):
    # Node-interleaved two-field linear elements have half-bandwidth three.
    diagonal=matrix.diagonal();assert np.all(diagonal>0)
    scale=1/np.sqrt(diagonal);A=(diags(scale)@matrix@diags(scale)).tocsc()
    band=np.zeros((4,len(diagonal)))
    for k in range(4):band[k,:len(diagonal)-k]=A.diagonal(-k)
    factor=cholesky_banded(band,lower=True,check_finite=True)
    return float(factor[0].min())


class Model:
    def __init__(self,bank,outer=2,coarse=False):
        self.original=base.Problem(bank,outer);p=self.original;heat=p.heat;bg=p.bg
        heat.d={k:heat.d[k] for k in heat.d.files}
        if coarse:
            si=bg.surface_index;indices=np.unique(np.r_[np.arange(0,si+1,2),si,np.arange(si+1,len(bg.grid),2),len(bg.grid)-1])
            self.grid=bg.grid[indices];self.surface_index=int(np.where(indices==si)[0][0])
        else:self.grid=bg.grid.copy();self.surface_index=bg.surface_index
        self.bg=bg;self.heat=heat;self.outer=outer;self.coarse=coarse
        n=len(self.grid);nf=self.surface_index+1;self.size=n+nf-1
        self.indices=np.full((n,2),-1,int)
        self.indices[:nf,0]=2*np.arange(nf);self.indices[:nf,1]=2*np.arange(nf)+1
        self.indices[nf:,1]=2*nf+np.arange(n-nf);self.indices[-1,1]=-1
        self.dofs=np.c_[self.indices[:-1],self.indices[1:]]
        dx=np.diff(self.grid);local=np.array([(1-1/np.sqrt(3))/2,(1+1/np.sqrt(3))/2])
        points=(self.grid[:-1,None]+dx[:,None]*local).ravel()
        data,self.points=coefficients(bg,points);nc=len(dx)
        self.canonical_data=data
        A=data[:,:4].reshape(-1,2,2);B=data[:,4:8].reshape(-1,2,2)
        C=data[:,8:12].reshape(-1,2,2);W=data[:,12:16].reshape(-1,2,2)
        gs=data[:,16:22].reshape(-1,2,3);hs=data[:,22:28].reshape(-1,2,3)
        shape=np.zeros((nc,2,2,4));derivative=shape.copy()
        shape[:,:,0,0]=1-local;shape[:,:,0,2]=local
        shape[:,:,1,1]=1-local;shape[:,:,1,3]=local
        derivative[:,:,0,0]=-1/dx[:,None];derivative[:,:,0,2]=1/dx[:,None]
        derivative[:,:,1,1]=-1/dx[:,None];derivative[:,:,1,3]=1/dx[:,None]
        shape=shape.reshape(-1,2,4);derivative=derivative.reshape(-1,2,4)
        covariant=derivative-np.einsum('nij,njk->nik',A,shape);weights=np.repeat(dx/2,2)
        stiffness=np.einsum('nia,nij,njb->nab',covariant,B,covariant)+np.einsum('nia,nij,njb->nab',shape,C,shape)
        inertia=np.einsum('nia,nij,njb->nab',shape,W,shape)
        rows=np.broadcast_to(self.dofs[:,:,None],(nc,4,4)).ravel();cols=np.broadcast_to(self.dofs[:,None,:],(nc,4,4)).ravel()
        valid=(rows>=0)&(cols>=0)
        def assemble(local_matrix):
            entries=np.sum((weights[:,None,None]*local_matrix).reshape(nc,2,4,4),axis=1).ravel()
            return coo_matrix((entries[valid],(rows[valid],cols[valid])),shape=(self.size,)*2).tocsc()
        self.K=assemble(stiffness);self.M=assemble(inertia)
        load=np.einsum('nia,nij,njk->nak',covariant,B,gs)-np.einsum('nia,nik->nak',shape,hs)
        rows=np.broadcast_to(np.repeat(self.dofs,2,axis=0)[:,:,None],load.shape).ravel()
        cols=np.broadcast_to((3*np.arange(len(points)))[:,None,None]+np.arange(3)[None,None,:],load.shape).ravel()
        valid=rows>=0
        test_map=coo_matrix(((weights[:,None,None]*load).ravel()[valid],(rows[valid],cols[valid])),shape=(self.size,3*len(points))).tocsc()
        src=source_points(heat,self.points)
        source=vstack(src).tocsr()[np.arange(3*len(points)).reshape(3,-1).T.ravel()]
        self.F=(test_map@source).tocsc()
        # The heat momentum lift has only the fluid displacement component
        # in Eulerian scalar coordinates. Its time history remains exact.
        fluid_r=self.grid[:nf];node=bg.sample(fluid_r);N,a=p.radiation.geometry.metric(fluid_r)
        ids=np.clip(np.searchsorted(heat.edges,fluid_r,side='right')-1,0,len(heat.edges)-2)
        f=(fluid_r-heat.edges[ids])/(heat.edges[ids+1]-heat.edges[ids])
        interp=coo_matrix((np.c_[1-f,f].ravel(),(np.repeat(np.arange(nf),2),np.c_[ids,ids+1].ravel())),shape=(nf,len(heat.edges))).tolil()
        inv=np.divide(1.,fluid_r**3,out=np.zeros_like(fluid_r),where=fluid_r>0)
        for i in np.flatnonzero(fluid_r<heat.edges[1]):interp.rows[i]=[1];interp.data[i]=[1.];inv[i]=heat.edges[1]**-3
        h=base.task.h;factor=np.zeros(nf);inside=(fluid_r<1)&(node['p']>0)
        factor[inside]=h.gr.G*1e-7*inv[inside]/(4*np.pi*h.gr.C**4*p.radiation.geometry.R*N[inside]*a[inside]*np.exp(-8*node['phi'][inside]**2)*(node['e'][inside]+node['p'][inside]))
        hz=(diags(factor)@interp.tocsr()).tocoo()
        self.H=coo_matrix((hz.data,(self.indices[hz.row,0],hz.col)),shape=(self.size,len(heat.edges))).tocsc()
        self.load=(self.F+self.K@self.H).tocsc()
        self.mass_pivot=banded_spd(self.M);self.shift_pivot=banded_spd(self.K+self.M)
        assert np.max(abs((self.K-self.K.T).data),initial=0)<1e-10*max(abs(self.K.data))


def main():
    assert not (OUT/'assembly.json').exists()
    write(OUT/'assembly-plan.json',dict(classification='Counterexample candidate',
        claim='Assemble the full canonical weak form with physical pressure, positive two-point element mass, exact heat source and momentum lift on the original grid.',
        changed='A different spatial discretization of the proved continuous equations, including the source and pressure reconstruction; not the rejected ad hoc change of D.',
        budget=dict(hard_seconds=60,new_EOS_calls=0,new_evolution_paths=0),source=prior.prior.prior.digest(Path(__file__)),canonical_source=prior.prior.prior.digest(Path(canonical.__file__))))
    signal.alarm(60);start=time.monotonic();m=Model(base.prior.OUT/'fine-bank.npz')
    result=dict(classification='Counterexample candidate',seconds=time.monotonic()-start,dofs=m.size,
        minimum_dr=float(np.diff(m.grid).min()),mass_cholesky_pivot=m.mass_pivot,
        stiffness_plus_mass_cholesky_pivot=m.shift_pivot,
        stiffness_nonsymmetry=float(np.max(abs((m.K-m.K.T).data),initial=0)/np.max(abs(m.K.data))),
        source_nnz=m.load.nnz,original_failure_resolved=False,full_dynamic_charge_solved=False)
    write(OUT/'assembly.json',result);signal.alarm(0);print('ASSEMBLY',json.dumps(result),flush=True)


if __name__=='__main__':main()
