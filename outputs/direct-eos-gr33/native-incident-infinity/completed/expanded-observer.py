def propagate(self,source):
    poly=base.flow.green.polynomial(self.t,source);H=poly.antiderivative()
    hc=H.c.reshape(3,len(self.t)-1,-1,self.order)/self.dx.reshape(-1,self.order)
    coeff=np.einsum('cij,ptcj->ptci',self.inverse,hc)
    U=[];Ut=[];Ux=[];cols=np.arange(len(self.ids))[None,:]
    for t in self.output_times:
        ret=t-self.distance;cut=np.clip(ret,0,self.t[-1]);idx=np.clip(np.searchsorted(self.t,cut,side='right')-1,0,len(self.t)-2);dt=cut-self.t[idx]
        hv=np.zeros_like(ret)
        for co in H.c:hv=hv*dt+co[idx,cols]
        sv=poly.c[0][idx,cols]*dt+poly.c[1][idx,cols]
        hv[ret<=0]=0.;sv[ret<=0]=0.
        hv=hv.reshape(len(self.tx),-1,self.order).sum(2)
        dv=(sv*self.sign).reshape(len(self.tx),-1,self.order).sum(2)
        sv=sv.reshape(len(self.tx),-1,self.order).sum(2)
        owner,cell,lo,hi=self.segments(t)
        if len(cell):
            xx=(lo[:,None]+hi[:,None])/2+(hi-lo)[:,None]*self.gx/2
            rr=t-abs(self.tx[owner,None]-xx)/C
            jj=np.clip(np.searchsorted(self.t,np.clip(rr,0,self.t[-1]),side='right')-1,0,len(self.t)-2)
            tt=np.clip(rr,0,self.t[-1])-self.t[jj]
            vv=leg.legvander((xx-self.mid[cell,None])/self.half[cell,None],self.order-1)
            cc=np.array([np.sum(co[jj,cell[:,None]]*vv,axis=-1) for co in coeff])
            hh=(cc[0]*tt+cc[1])*tt+cc[2];ss=2*cc[0]*tt+cc[1]
            hh[rr<=0]=0.;ss[rr<=0]=0.
            ww=(hi-lo)[:,None]*self.gw/2
            affected=np.unique(np.column_stack([owner,cell]),axis=0)
            ii,jj=affected.T;hv[ii,jj]=0.;sv[ii,jj]=0.;dv[ii,jj]=0.
            np.add.at(hv,(owner,cell),np.sum(ww*hh,axis=1))
            np.add.at(sv,(owner,cell),np.sum(ww*ss,axis=1))
            np.add.at(dv,(owner,cell),np.sum(ww*ss*np.sign(self.tx[owner,None]-xx),axis=1))
        U.append(C/2*hv.sum(1,dtype=LD));Ut.append(C/2*sv.sum(1,dtype=LD));Ux.append(-.5*dv.sum(1,dtype=LD))
    return np.asarray(U,float),np.asarray(Ut,float),np.asarray(Ux,float)

