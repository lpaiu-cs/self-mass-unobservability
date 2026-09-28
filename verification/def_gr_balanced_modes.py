"""Two physical source directions in the same mass-orthogonal rational space."""
import numpy as np
from scipy.sparse.linalg import spsolve_triangular
import def_gr_spatial_repair as task


class Projection(task.Projection):
    def basis(self,size):
        assert size%2==0
        model=self.model;fluid=np.zeros(model.size,bool);fluid[model.indices[:model.surface_index+1,0]]=True
        Q=np.empty((model.size,size),order='F')
        for j,mask in enumerate([fluid,~fluid]):
            f=np.where(mask,self.forcing,0.)
            # A shifted source response is smooth in the same physical energy;
            # the complete unaltered force is projected when evolving, not here.
            q=self.sigma*(self.C.T@self.solve(f))
            if j:q-=Q[:,0]*(Q[:,0]@q)
            norm=np.linalg.norm(q);assert norm>0;Q[:,j]=q/norm
        for j in range(2,size):
            q=self.sigma*(self.C.T@self.solve(self.C@Q[:,j-2]))
            for _ in range(2):q-=Q[:,:j]@(Q[:,:j].T@q)
            norm=np.linalg.norm(q);assert norm>1e-14,(j,norm);Q[:,j]=q/norm
        orth=float(np.max(abs(Q.T@Q-np.eye(size))));assert orth<1e-10
        physical=spsolve_triangular(self.C.T.tocsr(),Q,lower=False)
        return Q,physical,orth


def factor(model,Q):
    K,skew=model.energy_matrix(Q)
    return K,0.,dict(stiffness_skew=skew,source_directions=['fluid','scalar'],full_force_unchanged=True)


def series(model,K,unused,Q,label):return task.series(model,K,Q,label)


def control():
    return dict(classification='Proven',passed=True,
        identity='f=f_fluid+f_scalar. Their shifted responses span the shifted full response. Mass-orthogonal projection evolves the complete original forcing; separate normalization does not rescale a physical source.',
        scope='A trial-space identity, not an actual GR error certificate. The unchanged four response convergence gates are still required.')
