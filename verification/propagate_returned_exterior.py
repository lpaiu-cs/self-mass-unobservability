"""Continue the actually applied returned metric to exterior photon radii.

Counterexample candidate: same outgoing scalar continuation and signed Radau
emission as the accepted boundary. Exterior scalar scattering remains open.
"""
from pathlib import Path
import json, math, resource, sys, time
import numpy as np
import sympy as sp
from numpy.polynomial import Chebyshev, Polynomial
from scipy.integrate import solve_ivp
import propagate_exterior_vacuum as previous

p=previous.base; LD=p.LD; C=p.C; G=p.G
OUT=Path('native-returned-exterior262-work')
INPUT=p.ACTUAL; p.OUT=OUT
read,write,sha=p.read,p.write,p.sha


class Moments:
    """Exact moments of each declared linear emission/source interval."""
    def __init__(self,edges,a,b,T,degree):
        self.T=T; self.edges=np.asarray(edges)/T
        self.a=np.asarray(a,LD); self.b=np.asarray(b,LD)*LD(T)
        self.degree=degree
        lo=self.edges[:-1,None].astype(LD); hi=self.edges[1:,None].astype(LD)
        self.q=[]
        for k in range(degree+1):
            inc=(self.a-self.b*lo)*(hi**(k+1)-lo**(k+1))/(k+1)
            inc+=self.b*(hi**(k+2)-lo**(k+2))/(k+2)
            self.q.append(np.vstack([np.zeros((1,self.a.shape[1]),LD),np.cumsum(inc*LD(T),axis=0)]))

    def primitive(self,s,bins,k=0):
        v=np.clip(np.asarray(s)/self.T,0,self.edges[-1]).astype(LD)
        j=np.clip(np.searchsorted(self.edges,v,side='right')-1,0,len(self.edges)-2)
        lo=self.edges[j].astype(LD); a=self.a[j,bins]; b=self.b[j,bins]
        return self.q[k][j,bins]+LD(self.T)*((a-b*lo)*(v**(k+1)-lo**(k+1))/(k+1)+b*(v**(k+2)-lo**(k+2))/(k+2))

    def value(self,s,bins):
        v=np.asarray(s)/self.T; j=np.clip(np.searchsorted(self.edges,v,side='right')-1,0,len(self.edges)-2)
        result=self.a[j,bins]+self.b[j,bins]*(v-self.edges[j])
        return np.where((s>=0)&(s<=self.T),result,0)

    def convolution(self,t,lower,bins,coeff):
        # Kernel polynomial in age/T; integrate source moments analytically.
        end=np.maximum(t-lower,0); z=np.asarray(t)/self.T
        answer=np.zeros_like(z,dtype=LD)
        for j in range(len(coeff)):
            shifted=np.zeros_like(z,dtype=LD)
            for k in range(len(coeff)-1,j-1,-1):
                shifted=shifted*z+coeff[k]*math.comb(k,j)
            answer+=(-1)**j*shifted*self.primitive(end,bins,j)
        return np.where(t>=lower,answer,0)


class Metric(previous.Metric):
    def __init__(self,order,clock=128,degree=8):
        super().__init__(order); d=self.d; self.degree=degree
        self.field_path=INPUT/f'gr/fields-{clock}-g{order}.npz'
        self.metric_path=INPUT/f'metric/metric-{clock}-g{order}.npz'
        f=np.load(self.field_path); m=np.load(self.metric_path)
        self.clock=f['t']; self.U=f['U'][:,-1]; self.residual=m['asymptotic_mass_residual_cm']
        assert np.array_equal(self.clock,m['t']) and self.U[0]==0
        self.Udot=np.diff(self.U)/np.diff(self.clock)
        self.Rdot=np.diff(self.residual)/np.diff(self.clock)
        self.scalar=Moments(self.clock,self.U[:-1,None],self.Udot[:,None],d.T,degree)
        src=np.load(p.saved.full.source.saved(clock)); edges=src['actual_step_edges']
        e=p.saved.Emission(edges,src['accepted_angular_luminosity'][:2*(len(edges)-1)])
        self.emission=Moments(edges,e.a,e.b,d.T,degree)
        self.mu,self.mw,self.bins,self.rays,self.invariant=p.lapse.Lapse().rays(order,d.T)
        self.aw=self.mu*self.mw
        self.kernel_errors={}; self.inverse=[]
        self.coeff=self.fit(self.packet_kernel,'packet')
        self.qcoeff=self.fit(self.scalar_kernel,'scalar')
        self.dcoeff=np.arange(1,degree+1)[:,None]*self.coeff[1:]/d.T
        self.dqcoeff=np.arange(1,degree+1)[:,None]*self.qcoeff[1:]/d.T
        final=self.rays.sol(1).reshape(3,-1); self.reach=d.r0*final[0]
        for mu,end in zip(self.mu,self.reach):
            width=end-d.r0
            def rhs(s,y):
                r=d.r0+width*s; _,N,b,a,_=d.bg.metric(np.array([r/d.model.m.R]))
                direction=np.sqrt(1-(1-mu*mu)*(N[0]*d.r0/(d.z0['lapse'][0]*r))**2)
                return [width/(C*d.T*a[0]*direction)]
            sol=solve_ivp(rhs,[0,1],[0.],method='DOP853',rtol=2e-12,atol=2e-14,dense_output=True)
            assert sol.success; self.inverse.append(sol)
        self.scale=max(np.max(abs(m['delta_nu_faces'][:,-1])),np.max(abs(self.U*d.z0['Phi'][0])))
        assert self.scale>0

    def fit(self,function,name):
        nodes=(np.cos(np.pi*(np.arange(self.degree+1)+.5)/(self.degree+1))+1)/2
        values=function(nodes)[0]; coeff=[]
        for row in values.T:
            c=Chebyshev.fit(nodes,row,self.degree,domain=[0,1]).convert(kind=Polynomial).coef
            coeff.append(np.pad(c,(0,self.degree+1-len(c))))
        coeff=np.array(coeff).T
        query=np.linspace(0,1,257); exact,dot=function(query)
        value=np.polynomial.polynomial.polyval(query,coeff).T
        derivative=np.polynomial.polynomial.polyval(query,np.arange(1,len(coeff))[:,None]*coeff[1:]/self.d.T).T
        errors=[float(np.max(abs(v-w))/max(np.max(abs(w)),1e-290)) for v,w in [(value,exact),(derivative,dot)]]
        assert errors[0]<1e-10 and errors[1]<1e-7,(name,errors)
        self.kernel_errors[name]=errors
        return coeff

    def packet_kernel(self,s):
        d=self.d; y=self.rays.sol(np.asarray(s)).reshape(3,len(self.mu),-1)
        r=d.r0*y[0]; mu=y[1]; mass,N,b,a,Phi=d.bg.metric(r.ravel()/d.model.m.R)
        mass=(mass*d.model.m.R).reshape(r.shape); b=b.reshape(r.shape); a=a.reshape(r.shape); Phi=(Phi/d.model.m.R).reshape(r.shape)
        K=1/(r*a); nr=mass/(r*r*b)+r*Phi*Phi/2
        mdot=(1-mu*mu)*C*a*(1/r-nr)
        value=(1+mu*mu)*K
        dot=2*mu*mdot*K-(1+mu*mu)*K*C*a*mu/(r*b)
        return value.T,dot.T

    def scalar_kernel(self,s):
        d=self.d; r=d.inverse(C*d.T*np.asarray(s)); _,_,b,a,Phi=d.bg.metric(r/d.model.m.R); Phi/=d.model.m.R
        value=2*C*Phi*a/(r*b)
        dot=C*a*value*(-(2+1/b)/r+r*Phi*Phi)
        return value[:,None],dot[:,None]

    def wave(self,t,x,part='all'):
        t,x=np.broadcast_arrays(t,x); u=t-x/C
        j=np.clip(np.searchsorted(self.clock,u,side='right')-1,0,len(self.clock)-2)
        value=np.interp(u,self.clock,self.U,left=0,right=self.U[-1]); rate=np.where(u>=0,self.Udot[j],0)
        return value,rate,-rate/C

    def at(self,t,r):
        t,r=np.broadcast_arrays(t,r); shape=t.shape; t=t.ravel(); r=r.ravel(); d=self.d
        mass,N,b,a,Phi=d.bg.metric(r/d.model.m.R); mass*=d.model.m.R; Phi/=d.model.m.R
        x=d.optical(r); U,Ut,Ux=self.wave(t,x)
        zero=np.zeros(len(t),int); tau=x/C
        integral=self.scalar.convolution(t,tau,zero,self.qcoeff[:,0])
        integral_t=self.scalar.convolution(t,tau,zero,self.dqcoeff[:,0])+self.scalar.value(t-tau,zero)*np.polynomial.polynomial.polyval(tau/d.T,self.qcoeff[:,0])
        j=np.clip(np.searchsorted(self.clock,t,side='right')-1,0,len(self.clock)-2)
        residual=np.interp(t,self.clock,self.residual); residual_t=self.Rdot[j]
        outside=np.zeros(len(t),LD); flux=np.zeros(len(t),LD); kernel=np.zeros(len(t),LD); kernel_t=np.zeros(len(t),LD); pressure=np.zeros(len(t),LD)
        for k,mu0 in enumerate(self.mu):
            z=np.clip((r-d.r0)/(self.reach[k]-d.r0),0,1)
            age=d.T*self.inverse[k].sol(z)[0]; age=np.where(r>self.reach[k],2*d.T,age)
            ret=t-age; bins=np.full(len(t),self.bins[k],int)
            F=self.emission.primitive(np.maximum(ret,0),bins); L=self.emission.value(ret,bins)
            v=self.emission.convolution(t,age,bins,self.coeff[:,k])
            vt=self.emission.convolution(t,age,bins,self.dcoeff[:,k])+L*np.polynomial.polynomial.polyval(age/d.T,self.coeff[:,k])
            direction=np.sqrt(np.maximum(1-(1-mu0*mu0)*(N*d.r0/(d.z0['lapse'][0]*r))**2,0))
            weight=self.aw[k]; outside+=weight*F; flux+=weight*L; kernel+=weight*v; kernel_t+=weight*vt; pressure+=weight*L*direction
        fac=LD(G)/LD(C)**4; J=residual-fac*outside; Jt=residual_t-fac*flux; K=1/(r*a)
        lam=Phi*U+J*K; nu=Phi*U-integral-J*K-fac*kernel
        lt=Phi*Ut+Jt*K; nt=Phi*Ut-integral_t-Jt*K-fac*kernel_t
        nup=Phi*(Ux/a+U*(1/b-1)/r)+J/(r*r*N*b**1.5)+fac*pressure/(C*r*a*a)
        nr=mass/(r*r*b)+r*Phi**2/2; lr=r*Phi**2/2-mass/(r*r*b)
        return {k:np.asarray(v,float).reshape(shape) for k,v in dict(nu=nu,lam=lam,zeta=nu-lam,lt=lt,nt=nt,nup=nup,N=N,b=b,a=a,nr=nr,lr=lr).items()}


def self_check():
    # Independent Gauss integration of a discontinuous signed linear source.
    edges=np.array([0,.13,.6,1.]); a=np.array([[2.],[-3.],[1.]]); b=np.array([[.3],[.7],[-.4]])
    m=Moments(edges,a,b,1.,8); coeff=np.array([1.,-.3,.5,.07])
    gx,gw=np.polynomial.legendre.leggauss(8); errors=[]
    for t,lower in [(.08,0),(.4,.13),(.92,.11),(1.,.72)]:
        exact=0.
        for j in range(len(a)):
            lo=edges[j]; hi=min(edges[j+1],t-lower)
            if hi<=lo:continue
            te=lo+(hi-lo)*(gx+1)/2
            exact+=np.sum((hi-lo)*gw/2*(a[j,0]+b[j,0]*(te-lo))*np.polynomial.polynomial.polyval(t-te,coeff))
        actual=m.convolution(np.array([t]),np.array([lower]),np.array([0]),coeff)[0]
        errors.append(float(abs(actual-exact)))
    assert max(errors)<1e-13,errors
    r,a,b,c,J,F,L,mu,K=sp.symbols('r a b c J F L mu K',nonzero=True)
    # J_r=L/(c*a*mu), K_r=-K/(r*b), and the lower exterior
    # packet limit contributes +(1+mu^2)*K*L/(c*a*mu).
    derivative=-L/(c*a*mu)*K+J*K/(r*b)+(1+mu**2)*K*L/(c*a*mu)
    assert sp.simplify((derivative-J*K/(r*b)-L*mu/(c*r*a*a)).subs(K,1/(r*a)))==0
    return dict(classification='Proven',passed=True,polynomial_convolution_absolute=max(errors),
        constraint='C(r,t)=C_infinity(t)-G/c^4 sum angular_weight F(t-travel_time). lambda_mass=C/(r*a); nu_mass=-C/(r*a)-G/c^4 sum exterior E*(1+mu^2)/(r_packet*a_packet). Differentiating the moving lower limit supplies the radial photon pressure term.',
        scope='Declared piecewise-linear source times a polynomial kernel and exterior constraint identity; not continuum error certification.')


def prepare():
    assert not (OUT/'plan.json').exists()
    assert not OUT.exists() or all(v.name=='initialization' for v in OUT.iterdir())
    OUT.mkdir(exist_ok=True)
    assert read(INPUT/'result.json')['passed']
    files=[Path(__file__),Path(previous.__file__),Path(p.__file__),INPUT/'result.json']
    for n,q in [(128,8),(128,4),(64,8)]:
        files += [INPUT/f'gr/fields-{n}-g{q}.npz',INPUT/f'metric/metric-{n}-g{q}.npz',p.saved.full.source.saved(n)]
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Continue the actually applied returned metric from its boundary to every exterior radius, then propagate the same background photons through it.',
        decision='Supply the missing returned-metric photon stress on the accepted high/low history for the complete exterior boundary assembly.',
        representation='Same outgoing boundary U interpolation, actual signed Radau emission and homogeneous mass residual. Exact exterior linear constraints; polynomial approximations only to smooth static vacuum ray kernels, with independent value and derivative checks. No source smoothing, new physical grid or EOS/fluid replay.',
        gates=dict(boundary=.002,kernel_value=1e-10,kernel_derivative=1e-7,inverse=1e-9,energy_identity=.002,angular_invariant=1e-10,quadrature=.002,time=.02),
        budget=dict(check_seconds=600,pilot_seconds=1800,each_path_seconds=14400,threads_per_process=1,virtual_GiB=16),
        production='Start only after boundary and independent kernel checks and measured pilot. Reuse completed cohorts. Four original angular/emission/geometry controls; no tolerance or horizon enlargement.',
        remaining='Outgoing continuation matches the applied representation; exterior-generated/scattered scalar, reciprocal scalar/background-operator terms, matching inner work and new coupled boundary application remain required.',
        physical_final_charge_solved=False,bindings={str(f):sha(f) for f in dict.fromkeys(files)}))
    write(OUT/'symbolic.json',self_check())


def check():
    p.initialize(); m=Metric(8); d=m.d; raw=np.load(m.metric_path)
    t=m.clock; z=m.at(t,np.full(len(t),d.r0))
    old=raw['delta_nu_faces'][:,-1]; error=float(np.max(abs(z['nu']-old))/np.max(abs(old)))
    # Independent ray inversion, including near-grazing bins and early times.
    age=np.linspace(0,1,41); state=m.rays.sol(age).reshape(3,len(m.mu),-1)
    inverse=max(np.max(abs(sol.sol(np.clip((d.r0*state[0,k]-d.r0)/(m.reach[k]-d.r0),0,1))[0]-age)) for k,sol in enumerate(m.inverse))
    times=np.array([.213,.417,.713,.913])*d.T; radii=d.inverse(C*d.T*np.array([.013,.028,.071,.123]))
    dt=1e-7*d.T; dr=1e-7*C*d.T; mid=m.at(times,radii)
    before=m.at(times-dt,radii); after=m.at(times+dt,radii)
    left=m.at(times,radii-dr); right=m.at(times,radii+dr)
    derivatives={name:float(np.max(abs(value-mid[name]))/max(np.max(abs(mid[name])),1e-290)) for name,value in
        [('nt',(after['nu']-before['nu'])/(2*dt)),('lt',(after['lam']-before['lam'])/(2*dt)),('nup',(right['nu']-left['nu'])/(2*dr))]}
    result=dict(classification='Counterexample candidate',passed=error<.002 and inverse<1e-9,
        applied_boundary_relative=error,inverse_ray_age_fraction=float(inverse),kernel_controls=m.kernel_errors,
        independent_metric_derivative_relative=derivatives,
        returned_metric_scale=m.scale,source_knots=len(t),same_applied_outgoing_continuation=True,complete_physical_boundary=False)
    result['passed']=bool(result['passed'] and max(derivatives.values())<.002)
    write(OUT/'check.json',result); print(json.dumps(result),flush=True); assert result['passed'],result


def pilot():
    assert read(OUT/'check.json')['passed']; p.Metric=Metric; p.initialize(); m=p.Photons(8,8); rows=[]
    for cell in [0,7,15]:
        z,row=m.propagate(m.d.T,8,[cell]); np.savez_compressed(OUT/f'pilot-{cell}.npz',**z); rows.append(row)
        write(OUT/'pilot-progress.json',dict(rows=rows))
    upper=2*max(r['seconds']/r['packets'] for r in rows)*8*32*sum(range(1,17))+120
    result=dict(classification='Counterexample candidate',rows=rows,forecast_upper_seconds_per_full_path=upper,eligible=upper<14400)
    write(OUT/'pilot.json',result); print(json.dumps(result),flush=True)


if __name__=='__main__':
    action=sys.argv[1]; start=time.monotonic(); error=None
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,)*2); p.incident.native.deadline(1800 if action=='pilot' else 600)
    try:
        if action!='prepare':
            for f,h in read(OUT/'plan.json')['bindings'].items():assert sha(f)==h,f
        globals()[action]()
    except BaseException as exc:error=repr(exc); raise
    finally:
        if OUT.exists():write(OUT/f'{action}-receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
