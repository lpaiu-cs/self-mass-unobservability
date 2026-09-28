def jacobian(self,t,g):
    base=self.native(t,g);tangent,reset=self.selected[t]
    delta=np.maximum(np.max(abs(g),axis=0),LD('1e-100'))
    def action(q):
        reset(False);return self.native(t,q,tangent=tangent,metric=False)
    rows=[];cols=[];values=[];join=self.nb+2
    for component in range(4):
        for j in range(join):
            q=np.zeros_like(g);q[j,component]=delta[component];d=action(q)/delta[component]
            ii=np.flatnonzero(np.any(d!=0,axis=1))
            for r in range(4):
                rows.extend(4*ii+r);cols.extend([4*j+component]*len(ii));values.extend(d[ii,r])
        for color in range(5):
            jj=np.arange(join+color,self.n,5);q=np.zeros_like(g);q[jj,component]=delta[component]
            d=action(q)
            for offset in range(-2,3):
                ii=jj+offset;ok=(ii>=0)&(ii<self.n);ri,cj=ii[ok],jj[ok]
                for r in range(4):
                    rows.extend(4*ri+r);cols.extend(4*cj+component);values.extend(d[ri,r]/delta[component])
    J=sparse.coo_matrix((np.asarray(values,float),(rows,cols)),shape=(4*self.n,)*2).tocsr();J.eliminate_zeros()
    reset(False);offset=self.native(t,np.zeros_like(g),tangent=tangent)
    defect=(J@g.ravel()).reshape(self.n,4)+offset-base
    relative=np.sum(abs(defect)*self.units,axis=0)/np.maximum(np.sum(abs(base)*self.units,axis=0),LD('1e-290'))
    assert max(relative)<1e-12,('Selected native branch reconstruction',relative.tolist())
    return J,base
