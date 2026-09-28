"""Counterexample candidate: resolved anisotropic metric/scalar increments.

Evolve the retarded scalar response and solve its mass/pressure-volume feedback
using actual saved material/photon sources. Spatial material/packet transport
has not yet been re-evolved with these increments; no full-GR flag is earned.
"""
from pathlib import Path
from types import FunctionType
import json
import signal
import sys
import time
import numpy as np
import sympy as s
import def_native_projected_evolution as flow
import def_native_feedback_gr as previous

OUT=flow.OUT.parent/'def-native-anisotropic-gr';C=flow.C;G=flow.G
write=flow.write;sha=flow.sha;LD=np.longdouble


def symbolic():
    r,m,N,Phi,A,alpha,beta,Eg,Pg,Er,Pr,R4,Kg,f,fr,J,eF,pF,tF=s.symbols(
        'r m N Phi A alpha beta Eg Pg Er Pr R4 Kg f fr J eF pF tF',nonzero=True)
    b=1-2*m/r;E=Eg+Er;P=Pg+Pr;T=Eg-3*Pg;H=Eg+Pg
    mr=4*s.pi*r*r*A**4*E+r*r*b*Phi**2/2
    nu=m/(r*r*b)+4*s.pi*r*A**4*P/b+r*Phi**2/2
    lam=(mr/r-m/r**2)/b
    Phir=-(2/r+nu-lam)*Phi+4*s.pi*A**4*alpha*T/b
    dm=r*r*b*Phi*f+J;dl=dm/(r*b)
    de=eF-3*alpha*H*f-H*dl-4*alpha*Er*f-(Er+Pr)*dl
    dp=pF-3*alpha*Kg*f-Kg*dl-4*alpha*Pr*f-(3*Pr-R4)*dl
    dt=tF-(H-3*Kg)*(3*alpha*f+dl)
    Jr=4*s.pi*r*r*A**4*eF-(nu+lam)*J
    derived=sum(s.diff(dm,x)*dx for x,dx in [(r,1),(m,mr),(Phi,Phir),(f,fr),(J,Jr)])
    target=4*s.pi*r*r*A**4*(de+4*alpha*E*f)-r*Phi**2*dm+r*r*b*Phi*fr
    assert s.simplify(derived-target)==0
    dnu=(1+8*s.pi*r*r*A**4*P)*dm/(r*r*b*b)+4*s.pi*r*A**4*(dp+4*alpha*P*f)/b+r*Phi*fr
    dlprime=sum(s.diff(dl,x)*dx for x,dx in [(r,1),(m,mr),(Phi,Phir),(f,fr),(J,Jr)])
    K=2*N*N*Phi/(r*b)*(1+4*s.pi*r*r*A**4*(P-E))-8*s.pi*N*N*A**4*alpha*T/b
    full=r*N*N*b*Phi*(dnu-dlprime)-8*s.pi*r*N*N*A**4*alpha*T*dl
    full-=4*s.pi*r*N*N*A**4*(alpha*dt+(beta+4*alpha**2)*T*f)
    Dphi=4*s.pi*r*N*N*A**4*(r*Phi*alpha*(3*(H-Kg)+4*(Er-Pr))+3*alpha**2*(H-3*Kg))
    Dl=4*s.pi*r*N*N*A**4*(r*Phi*(E+P-Kg-3*Pr+R4)+alpha*(H-3*Kg))
    fc=K*r*r*b*Phi+16*s.pi*alpha*r*r*N*N*A**4*Phi*(P-E)-4*s.pi*r*N*N*A**4*(beta+4*alpha**2)*T+Dphi+Dl*r*Phi
    forcing=(K+Dl/(r*b))*J-4*s.pi*r*N*N*A**4*(alpha*tF+r*Phi*(eF-pF))+fc*f
    assert s.simplify(full-forcing)==0
    # Derive photon variations from a fixed canonical momentum packet.
    mu,l,h=s.symbols('mu l h');dn=-3*alpha*f-l;depacket=-alpha*f-mu**2*l
    dmu2=-2*mu**2*(1-mu**2)*l
    assert s.expand((dn+depacket)*mu**2+dmu2+4*alpha*f*mu**2+(3*mu**2-mu**4)*l)==0
    assert s.simplify((nu+lam)-r*Phi**2-4*s.pi*r*A**4*(E+P)/b)==0
    # Recover the old gas-only isotropic closure as a controlled limit.
    Dgas=4*s.pi*r*N*N*A**4*(r*Phi*(H-Kg)+alpha*(H-3*Kg))
    assert s.simplify(Dl.subs({Er:0,Pr:0,R4:0})-Dgas)==0
    assert s.simplify(Dphi.subs({Er:0,Pr:0,R4:0})-3*alpha*Dgas)==0
    return dict(classification='Proven',passed=True,
        premise='First variation on the momentarily scalar-balanced polar-areal initial slice, Pi_phi=0. Gas coordinate baryons, entropy and ionic inventory and photon canonical momenta are held fixed for the metric substep. Additional matter/packet transport is the prescribed source, not certified feedback.',
        gas='delta_Eg=eFg-(Eg+Pg)*(3*alpha*f+dlambda); delta_Pg=pFg-Kg*(3*alpha*f+dlambda).',
        photons='delta_Er=eFr-4*alpha*Er*f-(Er+Pr)*dlambda; delta_Pr=pFr-4*alpha*Pr*f-(3*Pr-R4)*dlambda. R4=int(E*mu^4); trace remains exactly zero.',
        mass='delta_m=r^2*b*Phi*f+J; J_prime+(nu_prime+lambda_prime)*J=4*pi*r^2*A^4*eF. This identity survives anisotropic radiation with its canonical-momentum response.',
        wave='U_tt/c^2-U_xx+Veff*U=Keff*J-4*pi*r*N^2*A^4*[alpha*tF+r*Phi*(eF-pF)], U=r*f, dx=dr/(N*sqrt(b)). The implemented Dphi,Dlambda and Veff are verified symbolically.',
        precision='Metric increments are independent state variables. Their being below the ulp of the full background does not make them physically zero.',
        limits='Not a nonlinear Einstein evolution, hydrodynamic stability theorem, continuum angular/spectral certificate or full response to exterior/deep-core sources.')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='c57efa4b6',
        claim='Derive the missing anisotropic conserved-phase-space metric closure and evolve separately stored retarded scalar/mass increments using the completed corrected-background trajectories. Establish the precision and actual forcing for the next spatial transport feedback.',
        decision='Do the coupled scalar potential/volume response and actual photon radial stresses preserve the represented direct sign? Determine whether ordinary background addition would erase the computed metric increments before attempting any new fluid run.',
        reuse='Completed Phase12164/128 matter/photon paths, initial20-point constraint panels and native bulk moduli. No fluid trajectory, EOS bank or wider physical horizon.',
        model='Momentarily scalar-balanced initial operator; canonical photon radial deformation retains its actual fourth angular moment. Gas responds isentropically at fixed coordinate baryons/inventory. Actual additional energy, radial stress and trace are saved forcing. Solve compact retarded Volterra potential feedback separately from the background, with zero initial increments.',
        limits='Spatial material/photon transport is still prescribed. Fields omit unrepresented exterior/deep-core forcing; finite cell reconstruction and saved time interpolation remain. No old GR lower bound is inherited.',
        budget=dict(source_seconds=25,field_seconds=90,max_potential_iterations=4,readout_seconds=30,CPU_threads=1,memory_GB=3,new_fluid_steps=0,new_native_calls=0),
        controls=dict(sources=[64,128],radial_orders=[4,8],field_times=17,field_locations='actual531 centers plus outer face',
            identity=1e-12,direct_reproduction=1e-9,time=.02,quadrature=.002,iteration=1e-8,contraction=.01),
        forecast='Saved source export/readout previously took9.3s. Scalar field evaluation uses about38million retarded source evaluations per order8 pass and blocks one time at a time; max four potential iterates. Stop at90s; no automatic mesh/time or iteration enlargement.',
        references=['https://arxiv.org/abs/gr-qc/9707041,2.23-2.26','https://arxiv.org/abs/gr-qc/0201064,233-237'],
        bindings={str(p):sha(p) for p in [Path(__file__),Path(flow.__file__),Path(previous.__file__),flow.OUT/'source-64.npz',flow.OUT/'source-128.npz',flow.OUT/'coupled-128.npz',flow.INPUT/'balanced-20.npz']}))
    write(OUT/'symbolic.json',symbolic())


def sources():
    FunctionType(previous.sources.__code__,dict(vars(previous),OUT=OUT,flow=flow))()


class Response:
    def __init__(self):
        self.model=model=flow.Coupled();self.bg=bg=model.m.bg;self.geo=flow.Geometry(model.m);b=model.bulk;f=model.flow
        inv=np.linalg.inv(np.polynomial.legendre.legvander(bg.q.x,19)).astype(LD)
        for k in ['gas_energy','gas_pressure','photon_energy','photon_radial_pressure']:bg.co[k]=bg.z['panel_'+k].astype(LD)@inv.T
        rho,v,lt,y=f.primitive(f.initial);f.eos.y=y;_,_,gam,*_=f.eos(rho,lt)
        self.gamma=np.r_[model.mech.K[0]/model.f0['p0'],gam]
        photons=np.concatenate([b.initial,model.initial_I]);mu4=(b.edges_mu[1:]**5-b.edges_mu[:-1]**5)/(5*np.diff(b.edges_mu))
        self.ratio4=np.einsum('iqf,q,f->i',photons,b.w*mu4,b.d['num']*b.d['Einf'])/np.einsum('iqf,q,f->i',photons,b.w,b.d['num']*b.d['Einf'])
        self.faces=bg.edges[-(b.n+model.n+1):]

    def coeff(self,r):
        z=self.bg.fields(r);ids=np.clip(np.searchsorted(self.faces,r,side='right')-1,0,len(self.gamma)-1)
        factor=G/C**4;Eg=z['gas_energy']*factor;Pg=z['gas_pressure']*factor;Er=z['photon_energy']*factor;Pr=z['photon_radial_pressure']*factor
        outside=r>self.bg.rend
        for value in [Eg,Pg,Er,Pr]:value[outside]=0
        Kg=self.gamma[ids]*Pg;R4=self.ratio4[ids]*Er;N=z['lapse'];b=z['b'];Phi=z['Phi'];alpha=-4*z['phi'];A4=np.exp(-8*z['phi']**2)
        E=Eg+Er;P=Pg+Pr;T=Eg-3*Pg;H=Eg+Pg
        K=2*N*N*Phi/(r*b)*(1+4*np.pi*r*r*A4*(P-E))-8*np.pi*N*N*A4*alpha*T/b
        Dphi=4*np.pi*r*N*N*A4*(r*Phi*alpha*(3*(H-Kg)+4*(Er-Pr))+3*alpha**2*(H-3*Kg))
        Dl=4*np.pi*r*N*N*A4*(r*Phi*(E+P-Kg-3*Pr+R4)+alpha*(H-3*Kg))
        fc=K*r*r*b*Phi+16*np.pi*alpha*r*r*N*N*A4*Phi*(P-E)-4*np.pi*r*N*N*A4*(-4+4*alpha**2)*T+Dphi+Dl*r*Phi
        V=N*N*(2*z['mass']/r**3+4*np.pi*A4*(P-E))-fc/r
        return dict(z,Eg=Eg,Pg=Pg,Er=Er,Pr=Pr,Kg=Kg,R4=R4,A4=A4,alpha=alpha,V=V,K=K+Dl/(r*b))

    def setup(self,d,order):
        edges=d['edges'];q=flow.initial.Quadrature(edges,order);rJ=q.r.ravel();w,delay,a,B,re=self.geo(rJ-self.model.m.RJ);r=re*self.model.m.R
        self.order=order;self.ids=np.repeat(np.arange(len(edges)-1),order);self.r=r;self.x=C*delay;self.z=z=self.coeff(r)
        measure=(rJ*rJ*B).reshape(q.r.shape);norm=measure@q.w
        self.weights=(measure*q.w/norm[:,None]).ravel()
        self.dx=(q.h[:,None]*q.w*B.reshape(q.r.shape)/a.reshape(q.r.shape)).ravel()
        apart=measure*a.reshape(q.r.shape)
        mean_a=(apart@q.w)/norm
        partial_a=(apart@q.Q.T)/norm[:,None]
        rest=np.asarray(d['baryon_g'],LD)*LD(d['cx'])*LD(C)**2
        energy=rest+d['gas_nonrest_energy_erg']+d['photon_energy_erg']
        killing=energy*mean_a;prefix=np.cumsum(killing,axis=1,dtype=LD)-killing
        bracket=-d['inner_cumulative_energy_erg'][:,None]+prefix[:,self.ids]+energy[:,self.ids]*partial_a.ravel()
        self.J=np.asarray(G/C**4*np.sqrt(z['b'])/z['lapse']*bracket,float)
        trace=np.asarray(rest+d['nonrest_trace_erg'],float);stress=d['metric_stress_erg']
        self.source=-G/C**4*self.weights*(w*trace[:,self.ids]+a*z['Phi']*stress[:,self.ids])+self.dx*z['K']*self.J
        self.direct=-G/C**4*self.weights*w*trace[:,self.ids]
        self.t=d['t'];self.data=d
        self.targets=np.r_[d['radius'],edges[-1]];self.tx=C*self.geo(self.targets-self.model.m.RJ)[1]
        self.distance=abs(self.tx[:,None]-self.x[None,:])/C
        self.sign=np.sign(self.tx[:,None]-self.x[None,:])
        self.tr=self.geo.metric(self.targets-self.model.m.RJ)[3]*self.model.m.R
        self.tz=self.coeff(self.tr)

    def propagate(self,source):
        poly=flow.green.polynomial(self.t,source);H=poly.antiderivative();U=[];Ut=[];Ux=[];cols=np.arange(len(self.ids))[None,:]
        for t in self.t:
            ret=t-self.distance;cut=np.clip(ret,0,self.t[-1]);idx=np.clip(np.searchsorted(self.t,cut,side='right')-1,0,len(self.t)-2);dt=cut-self.t[idx]
            value=np.zeros_like(ret)
            for co in H.c:value=value*dt+co[idx,cols]
            value[ret<=0]=0.;U.append(C/2*np.sum(value,axis=1,dtype=LD))
            value=poly.c[0][idx,cols]*dt+poly.c[1][idx,cols];value[ret<=0]=0.
            Ut.append(C/2*np.sum(value,axis=1,dtype=LD));Ux.append(-.5*np.sum(value*self.sign,axis=1,dtype=LD))
        return np.asarray(U,float),np.asarray(Ut,float),np.asarray(Ux,float)

    def run(self,steps,order):
        d=dict(np.load(OUT/f'source-{steps}.npz'));start=time.monotonic();self.setup(d,order)
        free=self.propagate(self.source);U=free[0];history=[]
        eta=float(C*self.t[-1]/2*np.sum(self.dx*abs(self.z['V']),dtype=LD));assert eta<.01
        # Piecewise-constant cell field, with a piecewise-linear stored time
        # representation. Its interpolation norm is <=1; no unstable slope
        # amplification is hidden in the Volterra contraction.
        for iteration in range(4):
            potential=-self.dx*self.z['V']*U[:,:-1][:,self.ids]
            add=self.propagate(potential);new=free[0]+add[0]
            err=float(np.max(abs(new-U))/max(np.max(abs(new)),1e-300));history.append(err);U=new
            if err<1e-8:break
        else:raise AssertionError(('Retarded potential iteration',history))
        Ut=free[1]+add[1];Ux=free[2]+add[2]
        z=self.tz;r=self.tr;f=U/r;fr=Ux/(z['lapse']*np.sqrt(z['b'])*r)-U/r**2
        # Interpolate the independently integrated J within the same cell;
        # the outer value is the final quadrature extrapolation here and is
        # not promoted to an exact outer mass constraint.
        J=np.array([np.interp(self.tx,self.x,row) for row in self.J])
        dm=r*r*z['b']*z['Phi']*f+J;dl=dm/(r*z['b']);volume=3*z['alpha']*f+dl
        de_g=-(z['Eg']+z['Pg'])*volume;dp_g=-z['Kg']*volume
        de_r=-4*z['alpha']*z['Er']*f-(z['Er']+z['Pr'])*dl
        dp_r=-4*z['alpha']*z['Pr']*f-(3*z['Pr']-z['R4'])*dl
        # All arrays below are increments, not background+increment sums.
        np.savez_compressed(OUT/f'fields-{steps}-g{order}.npz',t=self.t,radius_E=r,U=U,U_t=Ut,U_x=Ux,
            delta_phi=f,delta_Phi=fr,delta_mass_cm=dm,delta_lambda=dl,delta_log_proper_volume=volume,
            delta_gas_energy_geom=de_g,delta_gas_radial_pressure_geom=dp_g,
            delta_photon_energy_geom=de_r,delta_photon_radial_pressure_geom=dp_r,
            direct_and_mass_stress_U=free[0],potential_U=add[0],J_center_interpolated=J,potential_iterations=history)
        # At the outer face all represented sources lie on one side; the
        # original outgoing Green readout must be recovered exactly.
        direct=self.propagate(self.direct)[0][:,-1];M=float(d['M_cm'])
        row=dict(classification='Counterexample candidate',steps=steps,order=order,seconds=time.monotonic()-start,
            maximum_delta_phi=float(np.max(abs(f))),maximum_delta_mass_cm=float(np.max(abs(dm))),
            maximum_delta_lambda=float(np.max(abs(dl))),maximum_delta_log_proper_volume=float(np.max(abs(volume))),
            maximum_mass_increment_over_background_ulp=float(np.max(abs(dm)/np.spacing(z['mass']))),
            maximum_scalar_increment_over_background_ulp=float(np.max(abs(f)/np.spacing(z['phi']))),
            compact_potential_contraction=eta,potential_iteration_relative=history,
            endpoint_direct=-float(direct[-1])/M,endpoint_compact_with_metric=-float(U[-1,-1])/M,
            compact_potential_change=-float(add[0][-1,-1])/M,
            spatial_material_photon_feedback=False,full_GR_feedback=False,final_charge_solved=False)
        write(OUT/f'fields-{steps}-g{order}.json',row);print(json.dumps(row),flush=True);return row


def fields():
    assert not (OUT/'fields.json').exists();assert json.loads((OUT/'sources.json').read_text())['passed']
    for p,h in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert sha(p)==h,p
    signal.signal(signal.SIGALRM,flow.old.optical.timeout);signal.alarm(90);start=time.monotonic();model=Response();rows=[]
    for steps,order in [(128,8),(64,8),(128,4)]:rows.append(model.run(steps,order))
    a=np.load(OUT/'fields-128-g8.npz');b=np.load(OUT/'fields-64-g8.npz');c=np.load(OUT/'fields-128-g4.npz')
    errors={key:float(np.max(abs(a['U']-other['U']))/max(np.max(abs(a['U'])),1e-300)) for key,other in [('time',b),('quadrature',c)]}
    original=np.load(flow.OUT/'wave-128.npz');errors['direct_reproduction']=float(abs(rows[0]['endpoint_direct']/original['direct'][-1]-1))
    passed=errors['time']<.02 and errors['quadrature']<.002 and errors['direct_reproduction']<1e-9
    row=dict(classification='Counterexample candidate',passed=bool(passed),controls=errors,paths=rows,seconds=time.monotonic()-start,
        represented_anisotropic_GR_scalar_increment_evolved=True,conserved_volume_feedback_in_wave_operator=True,
        spatial_material_photon_feedback=False,unrepresented_exterior_source_enclosed=False,full_GR_feedback=False,final_charge_solved=False)
    write(OUT/'fields.json',row);print(json.dumps(row),flush=True);signal.alarm(0)


if __name__=='__main__':globals()[sys.argv[1]]()
