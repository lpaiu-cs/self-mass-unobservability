"""The same unreduced resolvent, solved in its existing narrow band."""
import numpy as np
from scipy.linalg.lapack import zgbtrf,zgbtrs
import def_gr_spatial_direct as direct


class Problem(direct.Problem):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs);n=self.model.size;coo=self.K.tocoo()
        self.width=int(np.max(abs(coo.row-coo.col)));b=self.width
        self.kband=np.zeros((3*b+1,n),order='F');self.mband=self.kband.copy(order='F')
        for k in range(-b,b+1):
            j0=max(k,0);j1=min(n,n+k)
            self.kband[2*b-k,j0:j1]=self.K.diagonal(k)
            self.mband[2*b-k,j0:j1]=self.M.diagonal(k)
        self.costs=[]

    def transform(self,s):
        m=self.model;heat=m.heat;tc=m.original.radiation.geometry.tc;lam=heat.rates*tc
        E=np.zeros(len(heat.edges),complex);E[heat.face_ids]=np.sum(heat.amplitude*tc*lam/(s*s*(s+lam)),axis=1)
        rhs=m.load@E;A=(self.K+s*s*self.M).tocsc();absolute=abs(A);scale=np.asarray(absolute.sum(1)).ravel()
        band=np.array(self.kband+s*s*self.mband,order='F');n=m.size;b=self.width
        for k in range(-b,b+1):
            j0=max(k,0);j1=min(n,n+k)
            band[2*b-k,j0:j1]/=scale[j0-k:j1-k]
        lu,piv,info=zgbtrf(band,b,b,overwrite_ab=1);assert info==0
        def solve(value):
            answer,info=zgbtrs(lu,b,b,(value/scale).reshape(-1,1),piv,overwrite_b=1)
            assert info==0;return answer[:,0]
        u=solve(rhs);denom=absolute@abs(u)+abs(rhs)+1e-100
        error=float(np.max(abs(rhs-A@u)/denom))
        if error>1e-12:
            extended=A.astype(np.clongdouble)
            for _ in range(2):u+=solve(np.asarray(rhs.astype(np.clongdouble)-extended@u.astype(np.clongdouble),complex))
            error=float(np.max(abs(rhs.astype(np.clongdouble)-extended@u.astype(np.clongdouble))/(absolute@abs(u)+abs(rhs)+1e-100)))
        self.error=max(self.error,error);assert error<1e-9
        inertial=(self.mass_solve(rhs.real)+1j*self.mass_solve(rhs.imag))/(s*s)
        return u,inertial
