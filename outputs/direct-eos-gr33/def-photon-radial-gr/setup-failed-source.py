"""Counterexample candidate: a finite spectral radiation patch coupled to radial GR.

Two adjacent material volumes and one photon-current face lie inside the real
outermost native cell. P1 closure and frozen transport geometry are explicit.
The existing full GR matrix is coupled in both directions through five ports;
no reduced mechanical basis or prescribed material-temperature history is used.
"""
from pathlib import Path
import argparse
import inspect
import json
import signal
import time
import resource
import numpy as np
from scipy import sparse
from scipy.linalg import lu_factor,lu_solve
from scipy.sparse.linalg import splu
import def_photon_population_coupled as ph
import def_gr_interface_patch as patch

go=patch.go;OUT=patch.OUT.parent/'def-photon-radial-gr';write=go.write


def geometry():
    b=np.load(ph.old.OUT/'bank.npz');d=np.load(go.task.fem.base.task.coupled.thermal.OUT/'coefficients.npz')
    R=float(go.task.fem.base.task.coupled.transport.Geometry().R)
    center=d['radius_cm'][0]/(100*R)
    half=.001*float(b['H'])/(float(d['A'][0]*d['metric'][0])*100*R)
    edges=center+half*np.array([-1.,0.,1.])
    assert d['faces_cm'][1]/(100*R)<edges[0]<edges[-1]<d['faces_cm'][0]/(100*R)
    return b,d,R,edges


class Background(patch.Background):
    def __init__(self,radiation,outer=2):
        super().__init__(radiation,outer)
        _,_,_,edges=geometry();extra=np.linspace(edges[0],edges[-1],65)
        self.grid=np.unique(np.r_[self.grid,extra]);self.surface_index=int(np.flatnonzero(self.grid==1)[0])
        self.nodes=self.sample(self.grid);self.mid=self.sample((self.grid[:-1]+self.grid[1:])/2)
        self.node_projection=radiation.projection(self.grid);self.mid_projection=radiation.projection(self.mid['r'])


class Mechanics:
    def __init__(self,degree=4):
        patch.base.Background=Background
        self.problem=p=go.Problem(degree);self.m=m=p.model
        b,d,R,self.edges=geometry();self.R=R;self.tc=m.original.radiation.geometry.tc
        self.geo=go.task.fem.base.task.h.gr.G*.1*R**2/go.task.fem.base.task.h.gr.C**4
        self.gr=go.task.fem.base.task.h.gr
        gx,gw=np.polynomial.legendre.leggauss(12)
        r=np.concatenate([(lo+hi)/2+(hi-lo)*gx/2 for lo,hi in zip(self.edges[:-1],self.edges[1:])])
        dr=np.concatenate([(hi-lo)*gw/2 for lo,hi in zip(self.edges[:-1],self.edges[1:])])
        point=m.bg.sample(r);N,a=m.original.radiation.geometry.metric(r);A4=np.exp(-8*point['phi']**2)
        weight=4*np.pi*(100*R)**3*dr*N*a*A4*r*r
        self.volumes=weight.reshape(2,-1).sum(1)
        avg=sparse.csr_matrix((weight/np.repeat(self.volumes,len(gx)),(np.repeat(np.arange(2),len(gx)),np.arange(len(r)))),shape=(2,len(r)))
        V,D=m.evaluation(r);mass,pressure,phi,v=[point[k] for k in ['m','p','phi','v']]
        bb=1-2*mass/r;alpha=-4*phi
        dlz=r*r*v*v/2-4*np.pi*r*r*A4*pressure/bb-mass/(r*bb)
        Dq=-sparse.diags(r)@D[0]-sparse.diags(3+dlz+3*alpha*r*v)@V[0]-sparse.diags(r*v+3*alpha)@V[1]
        self.Dq=(avg@Dq).tocsr()
        self.DS=(avg@(-sparse.diags(1/(r*bb))@self.sources(point)[2][:,:2])).toarray()
        self.gammaP=np.asarray(avg@(point['gamma']*point['p']/self.geo)).ravel()
        # Assemble all pressure, energy and enclosed-mass contributions together.
        data=m.data;B=data[:,4:8].reshape(-1,2,2);gs=data[:,16:22].reshape(-1,2,3);hs=data[:,22:28].reshape(-1,2,3)
        src=self.sources(m.points)
        g=[sum(sparse.diags(gs[:,i,j])@src[j] for j in range(3)) for i in range(2)]
        h=[sum(sparse.diags(hs[:,i,j])@src[j] for j in range(3)) for i in range(2)]
        F=sum(m.cov[i].T@sparse.diags(m.weights*B[:,i,j])@g[j] for i in range(2) for j in range(2))
        F-=sum(m.V[i].T@sparse.diags(m.weights)@h[i] for i in range(2))
        self.F=F.toarray().astype(np.longdouble)
        rf=np.array([self.edges[1]]);pf=m.bg.sample(rf);Nf,af=m.original.radiation.geometry.metric(rf)
        Af=np.exp(-2*pf['phi'][0]**2);self.luminosity_factor=float(4*np.pi*(100*R*rf[0])**2*Nf[0]**2*Af**4)
        self.face_volume=float(4*np.pi*(100*R*rf[0])**2*Nf[0]*af[0]*Af**4*(self.edges[2]-self.edges[0])*100*R/2)
        self.rate_clock=float(Af*Nf[0])*self.tc
        nf=m.surface_index+1;r=m.grid[:nf];point=m.bg.sample(r);N,a=m.original.radiation.geometry.metric(r)
        shape=np.interp(r,self.edges,[0.,1.,0.],left=0,right=0)
        coeff=np.zeros(nf);inside=(shape!=0)&(r>0)&(point['p']>0)
        coeff[inside]=self.gr.G*1e-7*self.tc*self.luminosity_factor*shape[inside]/(4*np.pi*self.gr.C**4*R*N[inside]*a[inside]*np.exp(-8*point['phi'][inside]**2)*(point['e'][inside]+point['p'][inside])*r[inside]**3)
        HQ=np.zeros(m.size);HQ[m.indices[:nf,0]]=coeff
        self.HQ=np.asarray(m.to_hierarchical@HQ,np.longdouble)
        value,_=m.evaluation(rf)
        self.Cv=(value[0]*(af[0]/Nf[0]*rf[0]*self.gr.C)).tocsr()
        self.error=0.

    def sources(self,point):
        r=point['r'];n=len(r);lo,mid,hi=self.edges
        mask=np.column_stack([(r>=lo)&(r<mid),(r>=mid)&(r<hi)]).astype(float)
        # Integrate exactly the same redshifted volume measure for J.
        gx,gw=np.polynomial.legendre.leggauss(12);enclosed=np.zeros((n,2))
        for j,(left,right) in enumerate(zip(self.edges[:-1],self.edges[1:])):
            upper=np.clip(r,left,right);x=(left+upper[:,None])/2+(upper-left)[:,None]*gx/2
            flat=x.ravel();bp=self.m.bg.sample(flat);N,a=self.m.original.radiation.geometry.metric(flat)
            val=(N*a*np.exp(-8*bp['phi']**2)*flat*flat).reshape(x.shape)
            enclosed[:,j]=4*np.pi*(100*self.R)**3*(upper-left)/2*(val@gw)
        N,a=self.m.original.radiation.geometry.metric(r)
        J=enclosed*self.gr.G*1e-7/(self.gr.C**4*self.R*(N*a)[:,None])
        rr=np.zeros((n,4));active=mask.sum(1)>0
        rr[active,2:]=-self.geo*mask[active]/(point['gamma'][active]*point['p'][active])[:,None]
        loss=np.c_[-self.geo*mask,np.zeros((n,2))]
        return [sparse.csr_matrix(rr),sparse.csr_matrix(loss),sparse.csr_matrix(np.c_[J,np.zeros((n,2))])]

    def factor(self,h):
        p=self.problem;z=1/h;matrix=(p.K+z*z*p.M).tocsc();scale=np.sqrt(abs(matrix.diagonal()));D=sparse.diags(1/scale)
        lu=splu((D@matrix@D).tocsc(),permc_spec='NATURAL')
        K=p.K.astype(np.longdouble);M=p.M.astype(np.longdouble)
        def solve(rhs):
            den=scale[:,None] if rhs.ndim==2 else scale
            q=(lu.solve(np.asarray(rhs/den,float))/den).astype(np.longdouble)
            for _ in range(3):q+=np.asarray(lu.solve(np.asarray((rhs-K@q-z*z*(M@q))/den,float))/den,np.longdouble)
            defect=rhs-K@q-z*z*(M@q)
            self.error=max(self.error,float(np.max(abs(defect)/(abs(K)@abs(q)+z*z*(abs(M)@abs(q))+abs(rhs)+1e-100))))
            assert self.error<1e-9
            return q
        ports=np.column_stack([self.F,-z*(M@self.HQ)])
        response=solve(ports)
        return solve,response,M


class Photons:
    def __init__(self,mechanics):
        self.gr=mech=mechanics;b,d,R,edges=geometry();self.b=dict(b)
        th=dict(np.load(ph.OUT/'thermo.npz'));mom=np.load(ph.OUT/'moments.npz')['fine']
        self.b['Cm']=th['Cf'];self.b.update(velocity_order=256,q_points=257)
        src=inspect.getsource(ph.photon_operator).replace('Operator(b,8,1,24)','Operator(b,2,0,24)')
        ns=dict(vars(ph));exec(compile(src,'radial-operator-source.py','exec'),ns);op=ns['photon_operator'](self.b)
        self.op=op;self.n=n=len(op.u);self.size=38+3*n
        A,B,weights=ph.coupled_blocks(op,th,mom);self.A=sparse.block_diag([A,A],format='csc')*mech.rate_clock
        BB=np.zeros((3*n,38));BB[:n,:19]=B;BB[n:2*n,19:]=B
        self.B=BB*mech.rate_clock
        diag=op.L.diagonal().reshape(n,2)
        self.diag=np.r_[diag[:,0],diag[:,0],diag[:,1]]*mech.rate_clock
        self.gains=[op.gain[0]*mech.rate_clock,op.gain[0]*mech.rate_clock,op.gain[1]*mech.rate_clock]
        self.off_bound=op.off_bound*mech.rate_clock
        c=float(b['c']);V=mech.volumes;Vf=mech.face_volume
        a=mech.tc*mech.luminosity_factor*c/np.sqrt(3*Vf*V)
        self.P=sparse.diags(self.diag,format='csc')
        ids=np.arange(n);rows=np.r_[ids,n+ids,2*n+ids,2*n+ids];cols=np.r_[2*n+ids,2*n+ids,ids,n+ids]
        vals=np.r_[np.full(n,a[0]),np.full(n,-a[1]),np.full(n,-a[0]),np.full(n,a[1])]
        self.P+=sparse.coo_matrix((vals,(rows,cols)),shape=(3*n,3*n)).tocsc()
        root=np.sqrt(op.Ci);T=float(b['T']);Cf=float(th['Cf']);U=th['U'];at=th['a'][-1]
        snap=np.load(ph.old.old.a.OUT/'eos-state.npz');st=np.load(ph.old.OUT/'station.npz')
        conversion=float(st['rho']*st['mass_scale'])*6.02214076e23
        def coords(i):return np.array([snap['number_fractions'][i,e,q+1:].sum() for e,q,_ in th['edges']])*conversion
        # Include derivative of the volume number density, not just fractions.
        arho=(coords(6)*np.exp(1e-4)-coords(5)*np.exp(-1e-4))/(2e-4)
        x0=th['coordinates'];brho=arho-x0;ub=U@brho;ua=U@arho;g=th['g']
        eos=snap['eos'][0];rho=float(b['material_rho']);pg=float(eos[1]);pT=pg*float(eos[6])/T
        latent=U.T@g;compT=(pg-rho*eos[9]+latent@brho)/np.sqrt(Cf)
        self.Hrho=np.zeros((self.size,2));self.Eq=np.zeros_like(self.Hrho)
        self.O=np.zeros((5,self.size));self.Drho=np.zeros((5,2))
        self.initial=np.zeros(self.size)
        # Frozen Taylor gradient from the two stored outer material states.
        theta=d['A']*d['N']*np.exp(d['lnT']);gradient=(theta[0]-theta[1])/(d['radius_cm'][0]-d['radius_cm'][1])
        centers=(edges[:-1]+edges[1:])/2;dt=gradient*(centers-edges[1])*100*R/(b['local_A']*b['local_N'])
        dt-=V@dt/V.sum();self.initial_temperature=dt
        for j in range(2):
            bath=slice(19*j,19*(j+1));radi=slice(38+n*j,38+n*(j+1));rv=np.sqrt(V[j])
            self.Hrho[bath,j]=rv*np.r_[compT,U@x0]
            self.Hrho[radi,j]=rv*root*T/3
            self.Eq[19*j+1:19*(j+1),j]=rv*ua
            self.O[j,bath]=weights/rv;self.O[j,radi]=root/rv
            self.O[j+2,bath]=np.r_[(pT+ub@g/T)/np.sqrt(Cf),-ub/T]/rv
            self.O[j+2,radi]=root/(3*rv)
            self.Drho[j,j]=-self.O[j]@self.Hrho[:,j]
            self.Drho[j+2,j]=pg*eos[5]+ub@ua/T-mech.gammaP[j]
            self.initial[bath]=rv*weights*dt[j];self.initial[radi]=rv*root*dt[j]
        self.O[4,38+2*n:]=c*root/np.sqrt(3*Vf)
        # v is in m/s; the photon library uses cm/s.
        self.Bv=np.zeros(self.size);self.Bv[38+2*n:]=np.sqrt(Vf)*root*T/(np.sqrt(3)*c)*100
        self.initial_outputs=self.O@self.initial
        self.energy=np.r_[np.sqrt(V[0])*weights,np.sqrt(V[1])*weights,np.sqrt(V[0])*root,np.sqrt(V[1])*root,np.zeros(n)]
        self.solver_error=0.;self.max_iterations=0
        radiation_energy=float(b['arad'])*T**4
        pressure_anchor=float(d['raw'][0,1])-radiation_energy/3-pg
        energy_anchor=rho*(float(d['raw'][0,2])-float(eos[2]))-radiation_energy
        self.input=dict(T_K=T,rho_g_cm3=rho,coordinate_radius_R=float(edges[1]),patch_edges_R=edges.tolist(),
            redshifted_volumes_cm3=V.tolist(),initial_temperature_offsets_K=dt.tolist(),photon_cells=n,
            pressure_native_difference=float(pg/(d['raw'][0,1]-float(b['arad'])*T**4/3)-1),
            heat_capacity_native_difference=float(b['Cm']/(rho*d['thermo'][0,3]/T-4*float(b['arad'])*T**3)-1),
            material_maxwell_relative=float(abs((pg-rho*eos[9])/(T*pT)-1)),
            free_energy_density_anchor=dict(constant_erg_cm3=-pressure_anchor,
                coefficient_erg_g=(energy_anchor+pressure_anchor)/rho,
                pressure_change_erg_cm3=pressure_anchor,energy_change_erg_cm3=energy_anchor,
                convention='Delta f_volume=a+b*rho independent of T. Delta p=-a and Delta(rho*u_lnrho)=-a cancel in the compression heat coefficient; pressure and heat-capacity tangents are unchanged.'))
        # An adjoint check tests every state, including the streaming face.
        ex,er=self.energy[:38],self.energy[38:]
        defect=np.r_[self.A.T@ex+self.B.T@er,self.B@ex+self.P.T@er+self.off(er)]
        scale=np.linalg.norm(self.P.T@er)+np.linalg.norm(self.B@ex)+np.linalg.norm(self.A@ex)
        self.checks=dict(total_energy_left_null=float(np.linalg.norm(defect)/scale),
            compression_cancel=float(np.max(abs(self.O[:2]@self.Hrho+self.Drho[:2]))),
            boost_enthalpy_relative=float(abs((self.O[4]@self.Bv)/(100*4*radiation_energy/3)-1)))
        assert self.checks['total_energy_left_null']<1e-10,self.checks
        assert self.checks['compression_cancel']<1e-10,self.checks
        assert self.checks['boost_enthalpy_relative']<1e-8,self.checks

    def off(self,rad):
        return np.concatenate([g@rad[i*self.n:(i+1)*self.n] for i,g in enumerate(self.gains)],axis=0)

    def collision(self,y):
        x,E=y[:38],y[38:]
        return np.concatenate([self.A@x+self.B.T@E,self.B@x+self.P@E+self.off(E)],axis=0)

    def inverse(self,h):
        factor=splu(sparse.eye(3*self.n,format='csc')+h*self.P)
        R=factor.solve(h*self.B);schur=lu_factor(np.eye(38)+h*self.A.toarray()-h*self.B.T@R)
        contraction=h*self.off_bound;assert contraction<.5,contraction
        def solve(rhs):
            x,E=rhs[:38],rhs[38:]
            def diagonal(right):
                z=factor.solve(right);xx=lu_solve(schur,x-h*self.B.T@z);return xx,z-R@xx
            xx,z=diagonal(E);norm=max(float(np.linalg.norm(rhs)),1e-100)
            for iteration in range(16):
                nx,nz=diagonal(E-h*self.off(z));error=np.linalg.norm(nz-z)*contraction/(1-contraction)/norm
                xx,z=nx,nz
                if error<1e-12:break
            else:raise AssertionError('photon inverse did not contract')
            y=np.concatenate([xx,z],axis=0);res=float(np.linalg.norm(y+h*self.collision(y)-rhs)/norm)
            self.solver_error=max(self.solver_error,res);self.max_iterations=max(self.max_iterations,iteration+1)
            assert res<1e-9,res
            return y
        return solve


class Coupled:
    def __init__(self,degree=4):
        self.m=Mechanics(degree);self.p=Photons(self.m)
        self.density_factor=lu_factor(np.eye(2)-self.m.DS@self.p.Drho[:2])
        # Eight local half-width light crossings resolve the initial exchange.
        self.duration=.008*float(self.p.b['H']/self.p.b['c'])/self.m.tc
        self.maximum_port_residual=0.

    def density(self,q,y):
        source=self.p.O[:2]@(y-self.p.initial)
        return lu_solve(self.density_factor,np.asarray(self.m.Dq@q+self.m.DS@source,float))

    def outputs(self,q,y,feedback=True):
        source=self.p.O@y-self.p.initial_outputs
        if not feedback:return source,np.asarray(self.m.Dq@q+self.m.DS@source[:2],float)
        rho=self.density(q,y)
        return source+self.p.Drho@rho,rho

    def stage_solver(self,h,feedback=True):
        m,p=self.m,self.p;grsolve,R,M=m.factor(h);phsolve=p.inverse(h)
        forcing=np.column_stack([p.Hrho+h*p.collision(p.Eq),-p.Bv])
        if not feedback:forcing*=0
        response=phsolve(forcing)
        T=p.O@response;T[:,:2]+=p.Drho if feedback else 0
        L=np.vstack([np.asarray(m.Dq@R,float)+np.c_[m.DS,np.zeros((2,3))],np.asarray(m.Cv@R,float)/h])
        small=lu_factor(np.eye(5)-T@L)
        def step(qrhs,wrhs,yrhs):
            oldrho=self.density(qrhs,yrhs);oldF=(p.O@yrhs)[4]
            oldv=float(m.Cv@(wrhs-m.HQ*oldF))
            rhs=yrhs-p.Hrho@oldrho+p.Bv*oldv if feedback else yrhs
            ybase=phsolve(rhs);base=p.O@ybase-p.initial_outputs
            qbase=grsolve(M@(qrhs/(h*h)+wrhs/h))
            offset=np.r_[np.asarray(m.Dq@qbase,float),float(m.Cv@(qbase-qrhs)/h)]
            ports=lu_solve(small,base+T@offset)
            motion=offset+L@ports;y=ybase+response@motion;q=qbase+R@ports
            w=(q-qrhs)/h+m.HQ*ports[4]
            actual,rho=self.outputs(q,y,feedback)
            residual=float(np.max(abs(actual-ports)/(abs(actual)+abs(ports)+np.max(abs(base))*1e-13+1e-100)))
            self.maximum_port_residual=max(self.maximum_port_residual,residual)
            # A norm criterion avoids promoting near-zero individual ports.
            assert np.linalg.norm(actual-ports)/max(np.linalg.norm(ports),1e-100)<1e-8
            return q,w,y
        return step

    def run(self,steps,label,feedback=True):
        start=time.monotonic();m,p=self.m,self.p;dt=self.duration/steps;gamma=1-1/np.sqrt(2)
        stage=self.stage_solver(gamma*dt,feedback);q=np.zeros(m.m.size,np.longdouble);w=q.copy();y=p.initial.copy()
        history=[];balance=0.;ph_scale=max(np.linalg.norm(p.initial),1.)
        def read(t):
            ports,rho=self.outputs(q,y,feedback);velocity=w-m.HQ*ports[4]
            native_v=m.problem.speed*(m.m.nativeV[0]@velocity);scalar=m.m.nativeV[1]@q
            weights=m.m.original.weights
            energy=m.volumes@ports[:2];scale=max(np.sum(abs(m.volumes*ports[:2])),ph_scale,1.)
            return dict(tau=float(t),velocity_mass_RMS_m_s=float(np.sqrt(weights@native_v**2)),
                scalar_mass_RMS=float(np.sqrt(weights@scalar**2)),face_velocity_m_s=float(m.Cv@velocity),
                luminosity_erg_s=float(m.luminosity_factor*ports[4]),cell_energy_erg=(m.volumes*ports[:2]).tolist(),
                density_perturbation=rho.tolist(),pressure_increment_erg_cm3=ports[2:4].tolist(),
                exchange_balance=float(abs(energy)/scale)),ports
        history.append(read(0)[0])
        for j in range(steps):
            q1,w1,y1=stage(q,w,y)
            factor=(1-gamma)/gamma
            q,w,y=stage(q+factor*(q1-q),w+factor*(w1-w),y+factor*(y1-y))
            row,ports=read((j+1)*dt);history.append(row);balance=max(balance,row['exchange_balance'])
        np.savez_compressed(OUT/(label+'.npz'),q=q,w=w,y=y,initial=p.initial,grid=m.m.grid,ports=ports,
            volume=m.volumes,density=self.density(q,y),source_energy=p.O[:2]@(y-p.initial),source_pressure=p.O[2:4]@(y-p.initial))
        row=dict(classification='Counterexample candidate',steps=steps,degree=m.m.degree,feedback=feedback,history=history,
            seconds=time.monotonic()-start,GR_residual=m.error,photon_residual=p.solver_error,
            exchange_balance=balance,maximum_port_component_residual=self.maximum_port_residual,
            photon_iterations=p.max_iterations,structural_checks=p.checks)
        write(OUT/(label+'.json'),row);print('RADIAL PHOTON GR',label,json.dumps({k:v for k,v in row.items() if k!='history'}),flush=True)
        return row


def symbolic():
    import sympy as s
    rho,a,b,p,ur,c=s.symbols('rho a b p ur c',nonzero=True)
    df=a+b*rho;dp=s.expand(rho*s.diff(df,rho)-df)
    dur=s.simplify(rho*s.diff(df/rho,rho))
    assert dp==-a and s.simplify(dp-rho*dur)==0
    E,v=s.symbols('E v')
    assert s.expand(E*v+E*v/3-4*E*v/3)==0
    VL,VR,L=s.symbols('VL VR L',positive=True)
    assert s.simplify(VL*(-L/VL)+VR*(L/VR))==0
    return dict(classification='Proven',passed=True,
        identity='Affine free-energy density anchors static p/e without changing p-rho*u_lnrho. Pair luminosity telescopes with the same redshifted volumes. The isotropic radiation boost carries (E+P)v=4Ev/3; replacing LTE photons must subtract this momentum as well.',
        scope='Algebraic identities of the declared local P1 matched-tangent model, not global opacity or atmosphere closure.')


def prepare():
    assert not OUT.exists();OUT.mkdir();signal.alarm(120);start=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(int(6e9),int(6e9)))
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='f470378f3',
        claim='Actually solve a conservative two-volume spectral photon/material patch and unreduced radial GR together, including photon pressure, stored ion energy, heat momentum and GR compression/velocity feedback.',
        scope='A local P1 finite-volume candidate inside the actual outer native cell. All24134 photon frequencies and18 retained ionization edges are retained. Two reservoirs exchange through one physical interior face; omitted surrounding photon channels are not certified zero. Frozen transport geometry and coefficients; no whole-atmosphere or full nonlinear claim.',
        thermodynamics='Use the Phase70 MHD chemical Hessian and radiation-free material heat capacity. Radiation energy/pressure are counted once. The small local MHD/background value mismatch is recorded: local free-energy constants anchor static pressure/energy to the saved background, while the new native tangent and its difference from the reference GR tangent enter the actual source ports. This is an explicit local matched-tangent candidate, not a new global native-EOS background.',
        domain='Patch half-width0.001 of the measured proper temperature scale;64 fixed mechanical subintervals; coordinate duration0.008*H/c, eight local half-width light crossings. The local frozen AN factor is retained in collision rates. Native outer-cell coordinates and no new EOS states.',
        method='Full GR5-port Schur coupling at every SDIRK2 stage. P1 source currents and spectral collisions share the same energy ledger. Compare32/64/128 time steps, then fixed p2 spatial and one-way controls at128 only if time passes. The one-way counterfactual prescribes thermal energy/pressure/current without GR compression, velocity or reference-adiabat subtraction; it is a code control, not another conservative physical star.',
        gates=dict(time_relative=.02,time_order=1.5,spatial_relative=.02,linear_residual=1e-9,energy_balance=1e-8,photon_residual=1e-9),
        budget=dict(pilot_seconds=120,production_seconds=600,total_seconds=900,CPU_threads=1,memory_GB=6,new_EOS_calls=0,automatic_expansion=False),
        references=['https://academic.oup.com/mnras/article/367/4/1739/1747567','verification/def_photon_matter_split.py'],
        bindings={str(path):go.task.digest(path) for path in [Path(__file__),Path(ph.__file__),Path(patch.__file__),Path(go.__file__),
            ph.OUT/'thermo.npz',ph.OUT/'moments.npz',ph.old.OUT/'bank.npz',patch.OUT/'result.json']}))
    write(OUT/'symbolic.json',symbolic())
    c=Coupled();write(OUT/'input.json',dict(classification='Counterexample candidate',**c.p.input,
        structural_checks=c.p.checks,coordinate_duration_seconds=c.duration*c.m.tc,GR_degrees_of_freedom=c.m.m.size))
    row=c.run(8,'pilot');elapsed=time.monotonic()-start
    forecast=1.4*(row['seconds']*480/8+3*(elapsed-row['seconds']))
    write(OUT/'pilot.json',dict(classification='Counterexample candidate',seconds=elapsed,path_seconds=row['seconds'],forecast_seconds=forecast,
        memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
        assumption='Eight-stage-count pilot scaled linearly to480 steps plus three setup costs and40percent margin; factorization reuse makes this conservative. No automatic time/radial/angular expansion.'))
    signal.alarm(0);print('RADIAL PILOT',elapsed,forecast,flush=True)


def run():
    assert not (OUT/'result.json').exists();plan=json.loads((OUT/'plan.json').read_text())
    for path,h in plan['bindings'].items():assert go.task.digest(Path(path))==h,path
    assert json.loads((OUT/'pilot.json').read_text())['forecast_seconds']<600
    signal.alarm(600);resource.setrlimit(resource.RLIMIT_AS,(int(6e9),int(6e9)));start=time.monotonic();c=Coupled()
    rows=[c.run(n,'fine-'+str(n)) for n in [32,64,128]]
    fields=['velocity_mass_RMS_m_s','scalar_mass_RMS','face_velocity_m_s','luminosity_erg_s'];comparisons={}
    for key in fields:
        a,b,d=[np.array([h[key] for h in row['history']]) for row in rows];scale=max(abs(d).max(),1e-100)
        first=np.max(abs(a-b[::2]))/scale;last=np.max(abs(b-d[::2]))/scale
        comparisons[key]=dict(previous=float(first),last=float(last),order=float(np.log2(first/last)))
    time_pass=all(r['last']<.02 and r['order']>1.5 for r in comparisons.values())
    spatial_pass=False;accepted_rows=list(rows)
    if time_pass:
        control=c.run(128,'one-way',False);fine=rows[-1];other=Coupled(2).run(128,'spatial-p2')
        for key in fields:
            x=np.array([h[key] for h in fine['history']]);norm=max(abs(x).max(),1e-100)
            comparisons[key]['spatial']=float(np.max(abs(x-np.array([h[key] for h in other['history']])))/norm)
            comparisons[key]['feedback_effect']=float(np.max(abs(x-np.array([h[key] for h in control['history']])))/norm)
        spatial_pass=all(r['spatial']<.02 for r in comparisons.values());accepted_rows+=[other]
    result=dict(classification='Counterexample candidate',actual_radial_photon_material_GR_evolved=True,
        bidirectional_compression_velocity_feedback=True,time_passed=time_pass,spatial_passed=spatial_pass,
        passed=bool(time_pass and spatial_pass and max(r['exchange_balance'] for r in accepted_rows)<1e-8),
        comparisons=comparisons,seconds=time.monotonic()-start,full_atmosphere_or_global_opacity=False,
        full_nonlinear_evolution=False,full_dynamic_charge_solved=False)
    write(OUT/'result.json',result);signal.alarm(0);print('RADIAL RESULT',json.dumps(result),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','run']);globals()[parser.parse_args().action]()
