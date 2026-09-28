def primitive(m,k,z,field,bank):
    row=m.point(k);raw=m.raw(k,np.zeros_like(z),np.zeros_like(field),0.);model=m.model;b=model.bulk;f=model.flow;nb=m.nb
    q=row['Q'];active=row['active'];s=3*field[0]+field[2];db=np.divide(z[0],q[0],out=np.zeros(m.n),where=active)
    dy=np.divide(z[3],q[3],out=np.zeros(m.n),where=active)-db
    beta=bank['beta'];dr=db-s;dv=np.zeros(m.n);dt=np.zeros(m.n);dp=np.zeros(m.n)
    dh=np.r_[0.,-np.cumsum(z[0,:nb])];xi=np.zeros(nb)
    p,u,ut,uy,pt,py,*_=b.eos.gas(row['theta'],row['eta']);eta=(1+row['eta'])*dy[:nb]
    dv[:nb]=z[1,:nb]/(model.cx*q[0,:nb]*C*C)-beta[:nb]*db[:nb]
    du=np.asarray(z[2,:nb].astype(LD)/(m.a[:nb]*q[0,:nb]),float)
    du-=(u+.5*model.cx*C*C*beta[:nb]**2)*db[:nb]+model.cx*C*C*beta[:nb]*dv[:nb]
    dt[:nb]=(du-bank['dr_u'][:nb]*dr[:nb]-uy*eta+b.eos.inventory[1]*xi)/ut
    dp[:nb]=bank['dr_p'][:nb]*dr[:nb]+pt*dt[:nb]+py*eta-b.eos.inventory[0]*xi
    ids=np.flatnonzero(active[nb:])+nb;rho=bank['rho'][ids].astype(LD);v=beta[ids].astype(LD);W2=1/(1-v*v)
    pp=bank['p'][ids].astype(LD);uu=bank['u'][ids].astype(LD);H=rho*(model.cx*LD(C)**2+uu)+pp
    pr,pt,py=[bank[key+'_p'][ids].astype(LD) for key in ['dr','dt','dy']]
    Hr=rho*(model.cx*LD(C)**2+uu+bank['dr_u'][ids])+pr;Ht=rho*bank['dt_u'][ids]+pt;Hy=rho*bank['dy_u'][ids]+py
    DD=(db[ids]-s[ids]).astype(LD);Y=dy[ids].astype(LD)
    S=H*W2*v;E=H*W2-pp
    RS=z[1,ids].astype(LD)/m.V[ids]-S*s[ids]-W2*v*(Hr*DD+Hy*Y)
    W=np.sqrt(W2);wminus1=W2*v*v/(W+1);w2minus1=W2*v*v
    K0=rho*model.cx*LD(C)**2*W*wminus1+rho*uu*W2+pp*w2minus1
    Kr=rho*model.cx*LD(C)**2*W*wminus1+rho*(uu+bank['dr_u'][ids])*W2+pr*w2minus1
    Kt=rho*bank['dt_u'][ids]*W2+pt*w2minus1;Ky=rho*bank['dy_u'][ids]*W2+py*w2minus1
    Kv=rho*model.cx*LD(C)**2*v*W**3*(2*W-1)+2*(rho*uu+pp)*W2*W2*v
    RE=z[2,ids].astype(LD)/(m.a[ids]*m.V[ids])-K0*s[ids]-Kr*DD-Ky*Y
    ST=W2*v*Ht;SV=H*W2*W2*(1+v*v)-W2*W2*v*v*Hr
    ET=Kt;EV=Kv-Kr*W2*v
    det=ST*EV-SV*ET;assert np.all(det!=0)
    dt[ids]=np.asarray((RS*EV-SV*RE)/det,float);dv[ids]=np.asarray((ST*RE-RS*ET)/det,float);dr[ids]-=np.asarray(W2*v*dv[ids],float)
    dp[ids]=np.asarray(pr*dr[ids]+pt*dt[ids]+py*Y,float)
    return dict(dr=dr,dv=dv,dt=dt,dy=dy,dp=dp)
