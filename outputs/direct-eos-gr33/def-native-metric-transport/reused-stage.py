def stage(self,h):
    h=np.longdouble(h);K=self.K.astype(np.longdouble);M=self.M.astype(np.longdouble)
    F=self.F.astype(np.longdouble);H=self.H.astype(np.longdouble)
    Lq=self.Lq.astype(np.longdouble);LE=self.LE.astype(np.longdouble)
    gain=np.r_[h*h*self.lam/(1+h*self.lam),h]
    GG=K+M/h**2;BB=F-M@H/h**2
    block=bmat([[GG,-BB@diags(self.energy_scale)],[-Lq,(diags(1/gain)-LE)@diags(self.energy_scale)]],format='csc')
    perm=self.permutation;AA=block[perm,:][:,perm].astype(float)
    row=np.asarray(abs(AA).max(axis=1).toarray()).ravel();aa=diags(1/row)@AA
    col=np.asarray(abs(aa).max(axis=0).toarray()).ravel()
    factor=condense_factor((aa@diags(1/col)).tocsc(),self)
    scale=np.sqrt(abs(GG.diagonal()))
    gr=splu((diags(1/scale)@GG@diags(1/scale)).astype(float).tocsc(),permc_spec='NATURAL')
    def invert(rhs):
        ans=np.empty(len(rhs));ans[perm]=factor.solve(np.asarray(rhs[perm]/row,float))/col
        return ans.astype(np.longdouble)
    def solve_once(state,t,correction,photons):
        qr,vr,er,dr=state;pred=er+h*dr
        rhs=M@(qr/h**2+vr/h)+t*self.F0+F@pred+photons
        vec=np.r_[rhs,t*self.LE0+LE@pred-dr+correction];answer=invert(vec)
        for _ in range(3):answer+=invert(vec-block@answer)
        delta=answer[self.size:]*self.energy_scale;target=rhs+BB@delta
        q=(gr.solve(np.asarray(target/scale,float))/scale).astype(np.longdouble)
        for _ in range(3):q+=(gr.solve(np.asarray((target-K@q-M@q/h**2)/scale,float))/scale).astype(np.longdouble)
        answer[:self.size]=q;defect=vec-block@answer
        error=float(max(abs(defect)/(abs(vec)+abs(block)@abs(answer)+1e-100)))
        self.max_error=max(self.max_error,error);assert error<1e-9,error
        e=pred+delta;targetflux=Lq@q+t*self.LE0+LE@e+correction
        d=np.r_[(dr[:-1]+h*self.lam*targetflux[:-1])/(1+h*self.lam),targetflux[-1]]
        # Local check in deviation coordinates: the large background
        # current cannot mask a failed small energy update.
        heat_error=float(max(abs(e-er-h*d)/(abs(e)+abs(er)+abs(h*d)+1e-100)))
        self.max_heat_error=max(self.max_heat_error,heat_error);assert heat_error<1e-9,heat_error
        return q,(q-qr)/h,e,d
    return self.close_stage(solve_once)
