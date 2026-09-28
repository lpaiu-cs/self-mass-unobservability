"""Evaluate the identical projected energy through its square root.

If C+shift*W is pointwise positive, K_r+shift*M_r=X.T X. QR/SVD of
X avoids forming the high-condition-number normal matrix before eigenvalues.
This changes arithmetic only, not the thermal source or trial space.
"""
import numpy as np
from scipy.linalg.lapack import dgeqrf
from scipy.linalg import svd
import def_gr_spatial_repair as task


def factor(model,Q):
    n=len(model.weights);k=Q.shape[1]
    C=model.data[:,8:12].reshape(-1,2,2);W=model.data[:,12:16].reshape(-1,2,2)
    # Declare the shift from local coefficient positivity, independently of
    # the response. Vacuum has only the scalar channel; the fluid row is zero.
    diag=np.diagonal(W,axis1=1,axis2=2);inv=np.divide(1.,np.sqrt(diag),out=np.zeros_like(diag),where=diag>0)
    local=C*inv[:,:,None]*inv[:,None,:];minimum=float(np.linalg.eigvalsh(local).min())
    shift=max(1.,1-minimum);positive=C+shift*W
    X=np.empty((4*n,k),order='F')
    B=model.data[:,4:8].reshape(-1,2,2)
    assert np.max(abs(B[:,0,1]))==0 and np.max(abs(B[:,1,0]))==0
    for lo in range(0,n,1024):
        hi=min(lo+1024,n);ix=slice(lo,hi);v=np.array([a[ix]@Q for a in model.V]);d=np.array([a[ix]@Q for a in model.cov]);w=model.weights[ix]
        X[lo:hi]=np.sqrt(w*B[ix,0,0])[:,None]*d[0]
        X[n+lo:n+hi]=np.sqrt(w*B[ix,1,1])[:,None]*d[1]
        # A two-by-two scaled Cholesky retains tiny fluid/scalar cross terms.
        p=positive[ix];a=np.sqrt(np.maximum(p[:,0,0],0));b=np.divide(p[:,1,0],a,out=np.zeros_like(a),where=a>0)
        rem=p[:,1,1]-b*b;assert np.all(rem>0)
        X[2*n+lo:2*n+hi]=np.sqrt(w)[:,None]*(a[:,None]*v[0]+b[:,None]*v[1])
        X[3*n+lo:3*n+hi]=np.sqrt(w*rem)[:,None]*v[1]
    _,_,work,info=dgeqrf(X,lwork=-1,overwrite_a=True);assert info==0
    qr,_,_,info=dgeqrf(X,lwork=int(work[0]),overwrite_a=True);assert info==0
    R=np.triu(qr[:k,:]).copy();del X,qr
    return R,shift,dict(local_scaled_potential_minimum=minimum,energy_shift=shift)


def spectrum(R,shift):
    _,s,Vt=svd(R,full_matrices=False,check_finite=True,lapack_driver='gesvd')
    values=s*s-shift;order=np.argsort(values)
    return values[order],Vt[order].T


def series(model,R,shift,Q,label):
    saved=task.modes.eigh
    try:
        task.modes.eigh=lambda *args,**kwargs:spectrum(R,shift)
        return task.series(model,np.empty((0,0)),Q,label)
    finally:task.modes.eigh=saved


def control():
    import sympy as s
    a,b,c,d,v1,v2=s.symbols('a b c d v1 v2',real=True)
    L=s.Matrix([[a,0],[b,c]]);v=s.Matrix([v1,v2])
    assert s.expand((v.T*L*L.T*v)[0]-((L.T*v).T*(L.T*v))[0])==0
    rng=np.random.default_rng(76);U,_=np.linalg.qr(rng.normal(size=(32,8)));V,_=np.linalg.qr(rng.normal(size=(8,8)))
    singular=np.geomspace(1.,1e6,8);X=(U*singular)@V.T;R=np.linalg.qr(X,mode='r')
    vals,_=spectrum(R,0.);error=float(np.max(abs(vals/singular**2-1)))
    assert error<1e-9,error
    return dict(classification='Proven',gram_identity_checked=True,numerical_classification='Counterexample candidate',
        passed=True,known_spectrum_relative_error=error,scope='The Gram identity is exact. SVD has improved conditioning on the declared numerical control; actual GR projection and spatial errors remain independently gated.')
