def stage(self,h):
    # Same SDIRK elimination as Phase83/89, including an algebraic outgoing
    # face. Scaling and residuals use the actual new thermal/GR blocks.
    h=np.longdouble(h);K=self.K.astype(np.longdouble);M=self.M.astype(np.longdouble);H=self.H.astype(np.longdouble)
    F=self.F.astype(np.longdouble);Lq=self.Lq.astype(np.longdouble);LE=self.LE.astype(np.longdouble)
    gain=np.r_[h*h*self.lam/(1+h*self.lam),h]
    GG=K+M/(h*h);BB=F-M@H/(h*h)
    block=bmat([[GG,-BB@diags(self.energy_scale)],[-Lq,(diags(1/gain)-LE)@diags(self.energy_scale)]],format='csc')
    perm=self.permutation;AA=block[perm,:][:,perm].astype(float)
    row=np.asarray(abs(AA).max(axis=1).toarray()).ravel();aa=diags(1/row)@AA
    col=np.asarray(abs(aa).max(axis=0).toarray()).ravel();factor=condense_factor((aa@diags(1/col)).tocsc(),self)
    scale=np.sqrt(abs(GG.diagonal()));grfactor=splu((diags(1/scale)@GG@diags(1/scale)).astype(float).tocsc(),permc_spec='NATURAL')
    def invert(rhs):
        ans=np.empty(len(rhs));ans[perm]=factor.solve(np.asarray(rhs[perm]/row,float))/col
        return ans.astype(np.longdouble)
    def step(state,t):
        qr,wr,Er,fr=state
        rhs=M@(qr/(h*h)+wr/h+H@Er/(h*h))+self.photon_force(float(t))
        Eb=Er.copy();Eb[:-1]+=h*(fr[:-1]+h*self.lam*self.bg.tc*self.L0[:-1])/(1+h*self.lam)
        Eb[-1]+=h*self.bg.tc*self.L0[-1]
        vec=np.r_[rhs,Eb/gain];answer=invert(vec)
        for _ in range(3):answer+=invert(vec-block@answer)
        E=answer[self.size:]*self.energy_scale;target=rhs+BB@E
        q=(grfactor.solve(np.asarray(target/scale,float))/scale).astype(np.longdouble)
        for _ in range(3):q+=(grfactor.solve(np.asarray((target-K@q-M@q/(h*h))/scale,float))/scale).astype(np.longdouble)
        answer[:self.size]=q;defect=vec-block@answer
        error=float(max(abs(defect)/(abs(vec)+abs(block)@abs(answer)+1e-100)));self.max_error=max(self.max_error,error)
        assert error<1e-9,error
        targetflux=self.bg.tc*self.L0+Lq@q+LE@E
        f=np.r_[(fr[:-1]+h*self.lam*targetflux[:-1])/(1+h*self.lam),targetflux[-1]]
        heat_error=float(max(abs(E-Er-h*f))/max(max(abs(E)),max(abs(Er)),1e-100));self.max_heat_error=max(self.max_heat_error,heat_error)
        assert heat_error<1e-9,heat_error
        w=(q-qr+H@(E-Er))/h
        return q,w,E,f
    return step
