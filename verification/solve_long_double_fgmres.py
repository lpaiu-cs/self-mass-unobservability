"""Flexible GMRES in long double with the existing double-precision preconditioner.

Counterexample candidate (numerical method). The coupled solves accept a linear
proposal only on the long-double true residual (vector 1e-14, physical moments
1e-13, per-component material gate). The existing correction is computed by
double-precision GMRES; near the end of the period the double operator and the
long-double operator differ enough that a correction can raise the true
residual (measured ratios 0.01..450 per refinement). Here the Arnoldi process,
the Hessenberg least squares and the update are carried in long double with the
same long-double operator used for acceptance. The preconditioner may stay in
double: flexible GMRES keeps the preconditioned vectors actually used, so their
rounding changes the search space, not the residual relation. Nothing here
changes an acceptance gate.
"""
import time
import numpy as np

LD=np.longdouble


def norm(v):
    return np.sqrt(np.dot(v,v))


def fgmres(op,precondition,b,restart=80,cycles=3,rtol=LD('1e-18'),log=None):
    """Return x with op(x)~b. op: long-double matvec; precondition: any matvec (double is fine)."""
    b=np.asarray(b,LD);x=np.zeros_like(b);bn=norm(b)
    if bn==0:return x
    for cycle in range(cycles):
        r=b-np.asarray(op(x),LD) if cycle else b.copy()
        beta=norm(r)
        if beta<=rtol*bn:break
        V=np.empty((restart+1,len(b)),LD);Z=np.empty((restart,len(b)),LD)
        V[0]=r/beta;H=np.zeros((restart+1,restart),LD);g=np.zeros(restart+1,LD);g[0]=beta
        cs=np.zeros(restart,LD);sn=np.zeros(restart,LD);k=0;start=time.monotonic()
        for j in range(restart):
            Z[j]=np.asarray(precondition(np.asarray(V[j],float)),LD)
            w=np.asarray(op(Z[j]),LD)
            for _ in range(2):  # modified Gram-Schmidt with one reorthogonalization pass
                for i in range(j+1):
                    c=np.dot(V[i],w);H[i,j]+=c;w-=c*V[i]
            hn=norm(w);H[j+1,j]=hn
            for i in range(j):
                a,d=H[i,j],H[i+1,j];H[i,j]=cs[i]*a+sn[i]*d;H[i+1,j]=-sn[i]*a+cs[i]*d
            d=np.sqrt(H[j,j]*H[j,j]+H[j+1,j]*H[j+1,j]);cs[j]=H[j,j]/d;sn[j]=H[j+1,j]/d
            H[j,j]=d;H[j+1,j]=0;g[j+1]=-sn[j]*g[j];g[j]=cs[j]*g[j];k=j+1
            if log is not None:log.append(dict(cycle=cycle,iteration=j+1,relative=float(abs(g[j+1])/bn),seconds=time.monotonic()-start))
            if abs(g[j+1])<=rtol*bn or hn==0:break
            V[j+1]=w/hn
        y=np.zeros(k,LD)
        for i in range(k-1,-1,-1):y[i]=(g[i]-np.dot(H[i,i+1:k],y[i+1:k]))/H[i,i]
        for i in range(k):x+=y[i]*Z[i]
        del V,Z
    return x


def self_check(n=300,seed=3):
    """Dense test (condition 1e10, consistent right-hand side): the true long-double residual falls below double precision."""
    rng=np.random.default_rng(seed);U,_=np.linalg.qr(rng.standard_normal((n,n)));W,_=np.linalg.qr(rng.standard_normal((n,n)))
    A=(U*np.logspace(0,-10,n))@W.T;ALD=A.astype(LD);b=ALD@rng.standard_normal(n).astype(LD)
    inverse=np.linalg.pinv(A)
    op=lambda v:ALD@np.asarray(v,LD);pre=lambda v:inverse@np.asarray(v,float)
    x=np.zeros(n,LD);history=[];inner=[]
    for _ in range(6):
        r=b-op(x);history.append(float(norm(r)/norm(b)))
        if history[-1]<1e-17:break
        x+=fgmres(op,pre,r,restart=40,cycles=2,log=inner)
    double=np.linalg.solve(A,np.asarray(b,float)).astype(LD);floor=float(norm(b-op(double))/norm(b))
    return dict(classification='Counterexample candidate',passed=history[-1]<1e-17,true_relative_history=history,
        double_direct_solve_relative=floor,inner_iterations=len(inner),
        scope='Synthetic dense check of the implementation only; no statement about the coupled operator.')
