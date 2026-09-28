def radiate(self,xb,I,u,theta,eta,h):
    if not self.feedback:return super().radiate(xb,I,u,theta,eta,h)
    b=self.bulk;beta=self.velocity();start_j=self.Pi.copy();initial_K=self.kinetic();trial_K=initial_K.copy();best=None;guess=xb*b.scale;gt=theta.copy();gy=eta.copy()
    for iteration in range(4):
        total,extra,escape,bound_extra=self.collision(guess,gt,gy,beta)
        en=b.d['num']*b.d['Einf'];a=b.d['a'];rho=self.mass/b.volume
        photon_force=np.einsum('iqf,q,f->i',total,b.w*b.mu,en)/a**4+escape[:,2]/a
        work=-(trial_K-initial_K)/(b.volume*h)
        b.extra=extra
        b.extra_u=-np.einsum('iqf,q,f->i',extra,b.w,en)/(a**4*rho)-escape[:,1]/(a*rho)+work/rho
        b.extra_y=np.einsum('iqf,q,f->i',bound_extra,b.w,b.d['num'])/(a**3*rho*b.d['thermo'][:,4]*b.d['y0'])
        result=super().radiate(xb,I,u,theta,eta,h);bb,aa,uu,tt,yy,port,it,err=result
        current,_,esc,_=self.collision(bb*b.scale,tt,yy,beta)
        force=-np.einsum('iqf,q,f->i',current,b.w*b.mu,en)/(a**4*C)-esc[:,2]/(a*C)
        next_j=start_j+h*b.volume*force
        self.Pi=next_j;new_K=self.kinetic();self.Pi=(start_j+next_j)/2;new_beta=self.velocity();change=float(max(abs(new_beta-beta)));self.Pi=start_j
        best=result,next_j,work,escape
        if iteration>0 and change<1e-13 and max(abs(gt-tt))<1e-11:break
        beta=new_beta;guess=bb*b.scale;gt=tt;gy=yy;trial_K=new_K
    else:raise AssertionError(('Interior moving-source iteration',change))
    self.maximum_inner_iterations=max(self.maximum_inner_iterations,iteration+1);self.Pi=best[1];self.max_frame=max(self.max_frame,float(max(abs(beta))));assert self.max_frame<1e-5
    self.radiation_work+=float(h*np.sum(b.volume*a*best[2]));self.deep_escape=getattr(self,'deep_escape',0.)+float(h*(b.volume@best[3][:,1]))
    return best[0]
