"""Cache invariant sparse data without changing the accepted GR equations."""
import numpy as np
from scipy.sparse import csc_matrix
import def_gr_full_weeks as go

OriginalProblem=go.Problem


class Problem(OriginalProblem):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        heat=self.model.heat;tc=np.longdouble(self.model.original.radiation.geometry.tc)
        self.lam=heat.rates.astype(np.longdouble)*tc
        self.numerator=heat.amplitude.astype(np.longdouble)*tc*self.lam
        self.absK=abs(self.Kx);self.absM=abs(self.Mx)
        pattern=(abs(self.K)+abs(self.M)).tocsc();pattern.sort_indices()
        self.indices=pattern.indices;self.indptr=pattern.indptr;self.columns=np.repeat(np.arange(self.model.size),np.diff(self.indptr))
        keys=self.columns*self.model.size+self.indices;assert np.all(np.diff(keys)>0)
        self.Kdata=np.zeros(len(keys));self.Mdata=self.Kdata.copy()
        for matrix,target in [(self.K,self.Kdata),(self.M,self.Mdata)]:
            a=matrix.tocoo();slots=np.searchsorted(keys,a.col.astype(np.int64)*self.model.size+a.row)
            assert np.array_equal(keys[slots],a.col.astype(np.int64)*self.model.size+a.row)
            target[slots]=a.data

    def transform(self,z):
        m=self.model;E=np.zeros(len(m.heat.edges),np.clongdouble)
        E[m.heat.face_ids]=np.sum(self.numerator/(z*z*(z+self.lam)),axis=1)
        rhs=self.load@E
        data=self.Kdata+complex(z*z)*self.Mdata
        matrix=csc_matrix((data,self.indices,self.indptr),shape=self.K.shape)
        scale=np.sqrt(abs(matrix.diagonal()));inverse=1/scale
        scaled=csc_matrix(((data*inverse[self.indices])*inverse[self.columns],self.indices,self.indptr),shape=self.K.shape)
        lu=go.splu(scaled);u=(lu.solve(np.asarray(rhs/scale,complex))/scale).astype(np.clongdouble)
        for _ in range(3):
            defect=rhs-self.Kx@u-z*z*(self.Mx@u)
            u+=(lu.solve(np.asarray(defect/scale,complex))/scale).astype(np.clongdouble)
        defect=rhs-self.Kx@u-z*z*(self.Mx@u);absolute=abs(u)
        error=float(np.max(abs(defect)/(self.absK@absolute+abs(z*z)*(self.absM@absolute)+abs(rhs)+1e-100)))
        self.error=max(self.error,error);assert self.error<1e-9
        return u
