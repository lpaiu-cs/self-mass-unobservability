"""The same polynomial trial space in endpoint-and-bubble coordinates.

A constant scalar has endpoint coefficients one and bubble coefficients zero.
Its derivative uses only the exact opposite endpoint pair, even in microscopic
surface cells. No polynomial degree, physical equation or source is changed.
"""
import numpy as np
from scipy.sparse import coo_matrix,diags,eye
import def_gr_spatial_repair as task

NodalModel=task.Model


class Model(NodalModel):
    def evaluation(self,r):
        ids=np.clip(np.searchsorted(self.cells,r,side='right')-1,0,len(self.cells)-2)
        width=np.diff(self.cells)[ids];x=(r-self.cells[ids])/width
        shapes=np.column_stack([p(x) for p in self.polynomials]);derivatives=np.column_stack([p.deriv()(x) for p in self.polynomials])/width[:,None]
        shapes[:,0]=1-x;shapes[:,-1]=x;derivatives[:,0]=-1/width;derivatives[:,-1]=1/width
        nodes=ids[:,None]*self.degree+np.arange(self.degree+1)[None,:];rows=np.repeat(np.arange(len(r)),self.degree+1)
        values=[];slopes=[]
        for field in range(2):
            columns=self.indices[nodes,field].ravel();valid=columns>=0
            values.append(coo_matrix((shapes.ravel()[valid],(rows[valid],columns[valid])),shape=(len(r),self.size)).tocsr())
            slopes.append(coo_matrix((derivatives.ravel()[valid],(rows[valid],columns[valid])),shape=(len(r),self.size)).tocsr())
        return values,slopes

    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        # Convert the nodally supplied heat lift to endpoint/bubble coordinates.
        row=[];col=[];data=[]
        for i in range(1,self.degree):
            mid=np.arange(len(self.cells)-1)*self.degree+i;left=mid-i;right=left+self.degree
            x=self.local[i]
            for field in range(2):
                target=self.indices[mid,field]
                for nodes,value in [(left,-(1-x)),(right,-x)]:
                    source=self.indices[nodes,field];valid=(target>=0)&(source>=0)
                    row.extend(target[valid]);col.extend(source[valid]);data.extend(np.full(valid.sum(),value))
        self.to_hierarchical=eye(self.size,format='csc')+coo_matrix((data,(row,col)),shape=(self.size,self.size)).tocsc()
        self.H=(self.to_hierarchical@self.H).tocsc()
        d=self.data;B=d[:,4:8].reshape(-1,2,2);gs=d[:,16:22].reshape(-1,2,3);hs=d[:,22:28].reshape(-1,2,3)
        src=task.fem.source_points(self.heat,self.points)
        g=[sum(diags(gs[:,i,j])@src[j] for j in range(3)) for i in range(2)]
        h=[sum(diags(hs[:,i,j])@src[j] for j in range(3)) for i in range(2)]
        F=sum(self.cov[i].T@diags(self.weights*B[:,i,j])@g[j] for i in range(2) for j in range(2))
        F-=sum(self.V[i].T@diags(self.weights)@h[i] for i in range(2))
        self.load=(F+self.K@self.H).tocsc()
        # Compare the same constant on quadrature rows away from the fixed
        # outer endpoint: neither representation satisfies psi(outer)=0 here.
        _,oldD=NodalModel.evaluation(self,self.points['r'])
        flat=np.zeros(self.size);flat[self.indices[:-1,1]]=1.;hier=self.to_hierarchical@flat
        interior=self.points['r']<self.cells[-2]
        before=oldD[1]@flat;after=self.D[1]@hier
        self.constant_gradient_before=float(abs(before[interior]).max())
        self.constant_gradient_after=float(abs(after[interior]).max())
        assert self.constant_gradient_after==0.,self.constant_gradient_after


Projection=task.Projection


def factor(model,Q):
    K,skew=model.energy_matrix(Q)
    return K,0.,dict(stiffness_skew=skew,coordinates='endpoint and bubbles',
        constant_scalar_gradient_before=model.constant_gradient_before,constant_scalar_gradient_after=model.constant_gradient_after)


def series(model,K,unused,Q,label):return task.series(model,K,Q,label)


def control():
    import sympy as s
    x=s.symbols('x');errors=[]
    for degree in [1,2,4]:
        nodes,polys=task.basis(degree);samples=np.linspace(0,1,51)
        old=np.array([p(samples) for p in polys]);new=old.copy();new[0]=1-samples;new[-1]=samples
        P=np.eye(degree+1)
        for i in range(1,degree):P[i,0]=1-nodes[i];P[i,-1]=nodes[i]
        error=float(np.max(abs(P.T@old-new)));assert error<1e-12;errors.append(error)
    assert s.diff((1-x)+x,x)==0
    return dict(classification='Proven',passed=True,identity='Nodal coefficients equal endpoint-linear interpolation plus independent bubbles. This is an invertible coordinate change of the identical degree-p space. The endpoint pair differentiates a constant to exactly zero.',
        numerical_classification='Counterexample candidate',basis_equivalence_errors=errors,
        scope='Basis equivalence and constant derivative only; actual propagation/spatial gates are retained.')
