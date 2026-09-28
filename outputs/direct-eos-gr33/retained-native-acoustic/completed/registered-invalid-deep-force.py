"""Stable finite HLL/acoustic differences, preserving the failed gross subtraction."""
from pathlib import Path
from types import MethodType
import inspect,sys,time,signal
import numpy as np
import sympy as sp
import def_retained_native_acoustic as native
prior=native.prior;OUT=native.OUT;DEST=OUT/'force-stable-final';LD=np.longdouble


def record_conserved(m):
    f=m.model.flow;old=f.conserved;f.acoustic_records=[]
    def conserved(self,rho,v,lt,y,a):
        U,F,thermo=old(rho,v,lt,y,a)
        if type(self.eos).__name__!='NativeEOS':return U,F,thermo
        p,u=thermo[:2];g=np.asarray(self.eos.parent.evaluate(rho,lt)[2],LD)
        c0=np.sqrt(g*p/np.maximum(np.asarray(rho,LD)*(self.eos.cx+u)+p,LD(1e-100)))
        assert np.all(np.isfinite(c0))
        self.acoustic_records.append(dict(U=U,F=F,v=np.asarray(v,LD),y=y,c0=c0,c1=thermo[-1]))
        return U,F,thermo
    f.conserved=MethodType(conserved,f)


def hll(left,right,scale,donor_left,donor_right):
    def speeds(i):
        l,r=left['v'],right['v'];cl,cr=left['c'+str(i)],right['c'+str(i)]
        return np.minimum(0,np.minimum((l-cl)/(1-l*cl),(r-cr)/(1-r*cr))),np.maximum(0,np.maximum((l+cl)/(1+l*cl),(r+cr)/(1+r*cr)))
    b,a=speeds(0);bn,an=speeds(1);da=an-a;db=bn-b;d=a-b;dn=an-bn;den=d*dn
    dw=np.divide(a*db-b*da,den,out=np.zeros_like(d),where=den>0)
    dj=np.divide(a*a*db-b*b*da+d*da*db,den,out=np.zeros_like(d),where=den>0)
    FL,FR,DU=left['F'],right['F'],right['U']-left['U']
    delta=(dw*(FL-FR)+dj*DU)*scale
    old=np.divide(a*FL-b*FR+a*b*DU,d,out=np.zeros_like(FL),where=d>0)*scale
    donor0=np.where(old[0]>=0,donor_left,donor_right)
    donor1=np.where(old[0]+delta[0]>=0,donor_left,donor_right)
    delta[3]=delta[0]*donor1+old[0]*(donor1-donor0)
    gross=(abs(dw)*(abs(FL-FR))+abs(dj)*abs(DU))*abs(scale)
    gross[3]=abs(delta[0]*donor1)+abs(old[0]*(donor1-donor0))
    # An independently subtracted full HLL flux agrees to its gross rounding.
    new=np.divide(an*FL-bn*FR+an*bn*DU,dn,out=np.zeros_like(FL),where=dn>0)*scale
    full=new-old;full[3]=new[0]*donor1-old[0]*donor0
    error=abs(full-delta);bound=64*np.finfo(LD).eps*(abs(new)+abs(old)+abs(delta))+1e-10*gross
    assert np.all(error<=bound+1e-100),float(np.max(error/np.maximum(bound,1e-300)))
    return delta,gross


def difference(m,u,theta,eta,K0):
    model=m.model;f=model.flow;geo=model.m;b=model.bulk
    records=f.acoustic_records;assert [len(r['v']) for r in records]==[512,513,513,1,1]
    y=records[0]['y'];ay=np.r_[f.incoming_y(y),y];by=np.r_[y,f.eos.y0]
    af,ae=hll(records[1],records[2],np.asarray(geo.af*geo.area,LD),ay,by)
    join,je=hll(records[3],records[4],np.asarray(geo.af[:1]*geo.area[:1],LD),records[3]['y'],np.array([y[0]]))
    af[:,0]=join[:,0];ae[:,0]=je[:,0]
    factors=LD(4*np.pi*geo.RJ**2*f.eos.rho0*prior.C)*np.array([1,prior.C**2,prior.C**2,f.eos.nH],LD)
    # Atmospheric returned momentum is converted from its original C-scaled
    # export in the SAME way as the existing whole-material owner.
    af=af*factors[:,None];ae=ae*factors[:,None]
    p=np.asarray(b.eos.gas(theta,eta)[0],LD);v=np.asarray(model.velocity()*prior.C,LD);u=np.asarray(u,LD)
    dp=p-model.f0['p0'];vl=np.r_[0.,v];vr=np.r_[v,0.];pl=np.r_[dp[0],dp];pr=np.r_[dp,dp[-1]]
    rho=np.asarray(model.face_rho,LD);s0=np.sqrt(np.asarray(K0,LD)/(model.cx*rho));s1=np.sqrt(np.asarray(model.face_K,LD)/(model.cx*rho))
    Z0=model.cx*rho*s0;Z1=model.cx*rho*s1;dZ=(np.asarray(model.face_K,LD)-K0)/(s1+s0)
    v0=(vl+vr)/2+(pl-pr)/(2*Z0);dv=-(pl-pr)*dZ/(2*Z0*Z1)
    ps0=(pl+pr)/2+Z0*(vl-vr)/2;dps=dZ*(vl-vr)/2
    v0[0]=0.;dv[0]=0.;ps0[0]=dp[0];dps[0]=0.
    area=np.asarray(model.area_gas,LD);mass0=area*rho*v0;dm=area*rho*dv;mass1=mass0+dm
    ul,ur=np.r_[u[0],u],np.r_[u,u[-1]];u0=np.where(mass0>=0,ul,ur);u1=np.where(mass1>=0,ul,ur)
    yy=np.asarray(b.d['thermo'][:,4]*b.d['y0']*(1+eta),LD);yl,yr=np.r_[yy[0],yy],np.r_[yy,yy[-1]]
    y0=np.where(mass0>=0,yl,yr);y1=np.where(mass1>=0,yl,yr);a=np.asarray(model.f0['af'],LD)
    e1=(a-model.m.a0)*model.cx*prior.C**2+a*(u1+(model.face_p+ps0+dps)/rho)
    de=dm*e1+mass0*a*(u1-u0+dps/rho)
    dpflux=area*(dps+model.cx*rho*dv*(2*v0+dv));dh=dm*y1+mass0*(y1-y0)
    df=np.array([dm,prior.C*dpflux,de,dh]);err=np.array([abs(dm),prior.C*area*(abs(dps)+abs(model.cx*rho*dv*(2*v0+dv))),abs(dm*e1)+abs(mass0*a*(u1-u0+dps/rho)),abs(dm*y1)+abs(mass0*(y1-y0))])
    df[:,-1]=af[:,0];err[:,-1]=ae[:,0]
    return np.c_[df[:,:-1],af],16*np.finfo(LD).eps*np.sum(np.c_[err[:,:-1],ae],axis=1)


def main():
    assert not DEST.exists();DEST.mkdir();start=time.monotonic()
    a,b,x,y=sp.symbols('a b da db');d=a-b;dn=d+x-y
    assert sp.factor((a+x)/dn-a/d-(a*y-b*x)/(d*dn))==0
    assert sp.factor((a+x)*(b+y)/dn-a*b/d-(a*a*y-b*b*x+d*x*y)/(d*dn))==0
    native.write(OUT/'stable-force-symbolic.json',dict(classification='Proven',passed=True,
        premise='Same conserved states and physical fluxes, changed HLL wave speeds; nonzero old/new speed spans.',
        weights='dw=(a*db-b*da)/((a-b)*(a+da-b-db)); dj=(a*a*db-b*b*da+(a-b)*da*db)/((a-b)*(a+da-b-db))',
        identity='delta_F=dw*(FL-FR)+dj*(UR-UL)',pressure_zero='Identical states or symmetric rest momentum give exact algebraic zero.'))
    source=inspect.getsource(native.force)
    source=source.replace("and not (OUT/'force.json').exists()", "and not (OUT/'force.json').exists()")
    source=source.replace("code=inspect.getsource(prior.hydro_sources)","code=inspect.getsource(prior.hydro_sources)\n    code=prior.replace(code,'    setup=time.monotonic()-start','    record_conserved(m)\\n    setup=time.monotonic()-start')\n    code=prior.replace(code,'            tick=time.monotonic();p=m.point(k);','            tick=time.monotonic();f.acoustic_records=[];p=m.point(k);')\n    code=prior.replace(code,\"                F=np.c_[df[:,:-1],af];G=np.r_[dg,ag]\",\"                F=np.c_[df[:,:-1],af];G=np.r_[dg,ag]\\n                stable_flux,stable_rounding=difference(m,m.d['snapshot_u'][m.ids[k]],nt,eta,baselineK)\")")
    source=source.replace("dF=F.astype(np.longdouble)-before[0].astype(np.longdouble)-previous['flux']",'dF=stable_flux')
    source=source.replace("dG=G.astype(np.longdouble)-before[1].astype(np.longdouble)-previous['gravity']",'dG=np.zeros_like(G,dtype=np.longdouble)')
    source=source.replace("    ns=dict(vars(prior),HYDRO=HYDRO,OLDHYDRO=prior.HYDRO,lookup=lookup,deep=deep,", "    code=prior.replace(code,'rounding=16*epsilon*np.sum(abs(F)+abs(before[0]),axis=1);rounding[1]+=16*epsilon*np.sum(abs(G)+abs(before[1]))','rounding=stable_rounding')\n    ns=dict(vars(prior),HYDRO=DEST,OLDHYDRO=prior.HYDRO,lookup=lookup,deep=deep,record_conserved=record_conserved,difference=difference,")
    source=source.replace("OUT/'force-plan.json'","OUT/'force-stable-plan.json'").replace("OUT/'expanded-native-force.py'","OUT/'expanded-stable-force.py'")
    ns=dict(vars(native),DEST=DEST,record_conserved=record_conserved,difference=difference)
    exec(compile(source,__file__,'exec'),ns);ns['force']()
    native.write(OUT/'force-repair.json',dict(classification='Counterexample candidate',passed=True,seconds=time.monotonic()-start,
        failed_force_preserved=True,new_native_calls=0,stable_delta_applied=True,source_sha256=native.sha(__file__)))


if __name__=='__main__':main()
