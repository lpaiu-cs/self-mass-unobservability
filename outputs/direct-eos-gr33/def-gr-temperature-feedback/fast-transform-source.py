def transform_pair(self,z):
    m=self.model; tc=np.longdouble(m.original.radiation.geometry.tc)
    gain=np.sum(self.conductance*tc*self.lam/(z*(z+self.lam)),axis=1)
    E0=np.zeros(len(m.heat.edges),np.clongdouble)
    E0[m.heat.face_ids]=np.sum(self.numerator/(z*z*(z+self.lam)),axis=1)
    data=self.Kdata+complex(z*z)*self.Mdata
    matrix=csc_matrix((data,self.indices,self.indptr),shape=self.K.shape)
    scale=np.sqrt(abs(matrix.diagonal())); inv=1/scale
    scaled=csc_matrix(((data*inv[self.indices])*inv[self.columns],self.indices,self.indptr),shape=self.K.shape)
    lu=go.splu(scaled)
    solves=0
    def solve(E):
        nonlocal solves
        refinements=3 if solves==0 else 1
        solves+=1
        rhs=self.load@E
        u=(lu.solve(np.asarray(rhs/scale,complex))/scale).astype(np.clongdouble)
        for _ in range(refinements):
            defect=rhs-self.Kx@u-z*z*(self.Mx@u)
            u+=(lu.solve(np.asarray(defect/scale,complex))/scale).astype(np.clongdouble)
        defect=rhs-self.Kx@u-z*z*(self.Mx@u)
        error=float(np.max(abs(defect)/(self.absK@abs(u)+abs(z*z)*(self.absM@abs(u))+abs(rhs)+1e-100)))
        self.error=max(self.error,error); assert error<1e-9
        return u
    def feedback(u,E):
        out=np.zeros_like(E);out[m.heat.face_ids]=gain*(self.Gq@u+self.GE@E)
        return out
    u0=solve(E0); first=feedback(u0,E0); denom=max(float(abs(first).max()),1e-100)
    dE=np.zeros_like(E0); du=np.zeros_like(u0); term=first
    for iteration in range(6):
        contribution=solve(term); dE+=term; du+=contribution
        next_term=feedback(contribution,term)
        relative=float(abs(next_term).max()/denom)
        if iteration==0:self.second_to_first=max(self.second_to_first,relative)
        if relative<1e-11:break
        term=next_term
    else:raise AssertionError('Temperature feedback did not converge in six corrections')
    defect=dE-first-feedback(du,dE)
    residual=float(abs(defect).max()/denom)
    self.loop_residual=max(self.loop_residual,residual);assert residual<1e-10
    self.iterations=max(self.iterations,iteration+1)
    self.feedback_relative=max(self.feedback_relative,float(abs(dE).max()/max(abs(E0).max(),1e-100)))
    return u0-m.H@E0, du-m.H@dE, E0, dE
