"""Conserved-volume metric feedback and a bounded scalar-wave resolvent.

Counterexample candidate: linear adiabatic metric response at fixed coordinate
baryon inventory, with the saved nonlinear gas as a prescribed source. Material
displacement feedback, scattered photons and a physical tail remain separate.
"""
from pathlib import Path
import argparse
import inspect
import json
import signal
import time
import numpy as np
import mpmath as mp
import sympy as s
import def_native_metric_charge as old

prior=old.prior
C,G=old.C,old.G
OUT=old.evolution.OUT.parent/'def-native-conserved-wave'
write=prior.write


def symbolic():
    r,m,N,Phi,A,alpha,E,P,gamma,f,J,eF,pF=s.symbols('r m N Phi A alpha E P gamma f J eF pF',nonzero=True)
    b=1-2*m/r;H=E+P;volume=(3*alpha+r*Phi)*f+J/(r*b)
    de=eF-H*volume;dp=pF-gamma*P*volume
    original=4*s.pi*r*r*A**4*(de+H*(3*alpha+r*Phi)*f)-r*Phi**2*J
    conserved=4*s.pi*r*r*A**4*eF-(r*Phi**2+4*s.pi*r*A**4*H/b)*J
    assert s.simplify(original-conserved)==0
    mr=4*s.pi*r*r*A**4*E+r*r*b*Phi**2/2
    nu=m/(r*r*b)+4*s.pi*r*A**4*P/b+r*Phi**2/2
    lam=(mr/r-m/r**2)/b
    assert s.simplify(nu+lam-r*Phi**2-4*s.pi*r*A**4*H/b)==0
    Be=-4*s.pi*N*N*A**4*(r*alpha+r*r*Phi)
    Bp=4*s.pi*N*N*A**4*(3*r*alpha+r*r*Phi)
    Dv=4*s.pi*N*N*A**4*r*(alpha*(H-3*gamma*P)+r*Phi*(H-gamma*P))
    assert s.simplify(Be*de+Bp*dp-(Be*eF+Bp*pF+Dv*volume))==0
    return dict(classification='Proven',passed=True,
        premise='First-order metric perturbations around the declared static isotropic background. Coordinate baryon inventory, entropy and composition are fixed for this algebraic response; additional material displacement and heating are source terms.',
        volume='dlnV=(3*alpha+r*Phi)*f+J/(r*b); dE=eF-(E+P)*dlnV; dP=pF-Gamma1*P*dlnV.',
        mass='J_prime+(nu_prime+lambda_prime)*J=4*pi*r^2*A^4*eF. Thus J=sqrt(b)/N times the cumulative prescribed Killing-energy perturbation in geometric units. The f-dependent source cancels only after the conserved-volume matter response is included.',
        wave='U_tt/c^2-U_xx+Veff*U=S_forced; Veff=Vfixed-Dv*(3*alpha+r*Phi)/r; KJ_eff=KJ+Dv/(r*b).',
        Dv='4*pi*N^2*A^4*r*[alpha*(E+P-3*Gamma1*P)+r*Phi*(E+P-Gamma1*P)].',
        resolvent='For eta=(c*Delta_u/2)*integral|Veff|dx<1, the retarded Green operator is a contraction. ||U-U0||<=eta*M0/(1-eta); ||U-U0-KU0||<=eta^2*M0/(1-eta), if ||U0||<=M0. This is a bound for the stated potential operator, not for missing physical source equations.')


def coeff(bg,r):
    z=bg.sample(np.asarray(r)/bg.R);b=1-2*z['m']/(r/bg.R);N=z['N'];A=np.exp(-2*z['phi']**2)
    Phi=z['v']/bg.R;alpha=-4*z['phi'];E=z['e']/bg.R**2;P=z['p']/bg.R**2;gam=z['gamma'];mass=z['m']*bg.R
    Kj=2*N*N*Phi/(r*b)*(1+4*np.pi*r*r*A**4*(P-E))-8*np.pi*N*N/b*alpha*A**4*(E-3*P)
    fc=Kj*r*r*b*Phi+16*np.pi*N*N*alpha*r*r*A**4*Phi*(P-E)-4*np.pi*N*N*r*A**4*(-4+4*alpha**2)*(E-3*P)
    Vf=N*N*(2*mass/r**3+4*np.pi*A**4*(P-E))-fc/r
    Dv=4*np.pi*N*N*A**4*r*(alpha*(E+P-3*gam*P)+r*Phi*(E+P-gam*P))
    return Vf-Dv*(3*alpha+r*Phi)/r,Kj+Dv/(r*b),b,N


def readout(n,channel,direction,utimes,raw,packets):
    base=prior.Geometry if channel=='trace' else old.MetricGeometry
    class Geometry(base):
        def __call__(self,x):
            w,d,a,B,re=super().__call__(x)
            return w,direction*d,a,B,re
    env=dict(vars(prior),OUT=OUT/f'{channel}-{direction}',prior=type('Inputs',(),dict(Flow=staticmethod(old.adapter),OUT=old.evolution.OUT)),
        Geometry=Geometry,packets=packets,channel=channel,direction=direction)
    acoustic=inspect.getsource(prior.acoustic)
    if direction==-1:
        assert acoustic.count('(1+cs)')==4
        acoustic=acoustic.replace('(1+cs)','(1-cs)')
    if channel=='stress':acoustic=acoustic.replace('u+p/rr-3*gamma*p/rr','u+p/rr-gamma*p/rr')
    exec(compile(acoustic,__file__,'exec'),env)
    source=inspect.getsource(prior.readout)
    if channel=='stress':source=source.replace('rho*u-3*p','rho*u-p')
    anchor='    nr=np.array(nr);nr-=nr[0]'
    assert source.count(anchor)==1
    source=source.replace(anchor,anchor+'\n    if direction==1: packets[(n,channel)]=dict(times=t,baryon=np.asarray(baryon,float),nonrest=nr,weights=w,cx=m.eos.cx)')
    exec(compile(source,__file__,'exec'),env)
    return env['readout'](n,'linear',utimes,raw)


def mass_wave(n,utimes):
    m=old.adapter(n);geom=prior.Geometry(m);d=np.load(old.evolution.OUT/f'cells-{n}.npz');h=d['history'];t=np.r_[0,h[:,0]]
    states=np.concatenate([d['initial'][None],d['snapshots']]).astype(np.longdouble);scale=np.longdouble(4*np.pi*m.RJ**2*m.eos.rho0)
    delta=(states-states[0])*m.vol*scale;pb=np.r_[0,h[:,3]].astype(np.longdouble)*scale;pk=np.r_[0,h[:,4]-h[:,5]].astype(np.longdouble)*scale
    energy=delta[:,2]+np.longdouble(m.a0*m.eos.cx)*delta[:,0]
    prefix=-pk[:,None]-np.longdouble(m.a0*m.eos.cx)*pb[:,None]+np.cumsum(energy,axis=1)-energy/2
    _,delay,a,B,re=geom(m.x);r=re*m.R;V,kj,b,N=coeff(m.bg,r)
    J=np.asarray(G/C**2*np.sqrt(b)/N*prefix,float);hj=prior.polynomial(t,J).antiderivative();waves=[]
    for sign in [1,-1]:
        waves.append(np.array([C/2*np.sum(m.dx*B/a*kj*prior.paired(hj,tt+sign*delay),dtype=np.longdouble) for tt in utimes],float))
    M0J=C/2*float(t[-1])*np.sum(m.dx*B/a*abs(kj)*np.max(abs(J),axis=0),dtype=np.longdouble)
    previous=np.load(old.OUT/f'mass-source-{n}.npz')['J_forced_cm']
    change=float(np.max(abs(J-previous))/max(np.max(abs(J)),1e-100))
    np.savez_compressed(OUT/f'mass-{n}.npz',times=t,radius_cm=r,J_cm=J,coefficient_cm_minus2=kj,potential_cm_minus2=V,
        u_seconds=utimes,outgoing_cm=waves[0],incoming_cm=waves[1])
    return m,np.array(waves),float(M0J),pb,change


def sample_H(poly,t):
    return np.where(t<=poly.x[0],0,poly(np.clip(t,poly.x[0],poly.x[-1])))


def q_of_r(bg,r):
    r=np.asarray(r);gx,gw=np.polynomial.legendre.leggauss(16)
    rr=r[:,None]+(bg.R-r[:,None])*(gx+1)/2
    z=bg.sample(rr.ravel()/bg.R);cc=(z['N']*np.sqrt(1-2*z['m']/(rr.ravel()/bg.R))).reshape(rr.shape)
    return -(bg.R-r)/2*(1/cc@gw)


def potential_correction(m,utimes,waves,cells,order):
    bg=m.bg;geom=prior.Geometry(m);lo,hi=geom(np.array([-28000.,120000.]))[1]*C
    duration=float(utimes[-1]-utimes[0]);qmin=lo-C*duration/2
    edges=np.linspace(qmin,lo,cells+1);gx,gw=np.polynomial.legendre.leggauss(4)
    q=(edges[:-1,None]+np.diff(edges)[:,None]*(gx+1)/2).ravel();weights=(np.diff(edges)[:,None]*gw/2).ravel()
    r=bg.R+q*bg.Nb*np.sqrt(1-2*bg.mu)
    for _ in range(3):
        z=bg.sample(r/bg.R);cc=z['N']*np.sqrt(1-2*z['m']/(r/bg.R));r-=(q_of_r(bg,r)-q)*cc
    qerror=float(max(abs(q_of_r(bg,r)-q)));assert qerror<1e-4
    V=coeff(bg,r)[0];Hin=prior.polynomial(utimes,waves[1]).antiderivative();Hout=prior.polynomial(utimes,waves[0]).antiderivative()
    inner=np.array([-C/2*np.sum(weights*V*sample_H(Hin,u+2*q/C),dtype=np.longdouble) for u in utimes],float)
    # Compactify the exact vacuum integral. No outer-radius cutoff is used.
    rhi=geom.metric(np.array([120000.]))[3][0]*m.R;gx,gw=np.polynomial.legendre.leggauss(order);zmax=bg.R/rhi
    zz=(gx+1)*zmax/2;rr=bg.R/zz;mm,nn,bb,cc,phi=bg.metric(rr/bg.R);phi/=bg.R;mm*=bg.R
    vv=2*nn*nn*(mm/rr**3-phi*phi);Iout=float(np.sum(gw*zmax/2*vv*bg.R/(zz*zz*cc)))
    outer=-C/2*Iout*Hout(utimes)
    return inner+outer,dict(cells=cells,exterior_order=order,q_inverse_error_cm=qerror,exterior_potential_integral_per_cm=Iout,
        endpoint_inner_cm=float(inner[-1]),endpoint_outer_cm=float(outer[-1]),compact_q_cm=[float(lo),float(hi)],qmin_cm=float(qmin))


def bound(m,utimes,packets,n,rows,M0J,pb,quadrature):
    # All rounding in the bound arithmetic is outward interval arithmetic.
    # These constants describe the frozen saved coefficient/source model.
    mp.iv.dps=40;iv=mp.iv;I=iv.mpf
    def up(z):return float(np.nextafter(float(z.b),np.inf))
    bg=m.bg;R=I(bg.R);M=I(bg.M*bg.R);T=I(float(utimes[-1]-utimes[0]));qlo,qhi=quadrature['compact_q_cm'];qmin=I(qlo)-I(C)*T/2
    rmin=R+qmin;idx=max(0,int(np.searchsorted(bg.d['radius_cm'],float(rmin.a)))-1)
    E=I(float(max(bg.d['energy_cgs'][idx:])))*I(G)/I(C)**4
    P=I(float(max(bg.d['pressure_cgs'][idx:])))*I(G)/I(C)**4
    phi=I(float(max(abs(bg.d['phi_prime_cm'][idx:]))));alpha=4*I(float(max(abs(bg.d['phi'][idx:]))));gamma=I(float(max(bg.d['gamma1'][idx:])))
    H=E+P;gp=gamma*P
    Vbound=2*M/rmin**3+4*iv.pi*H+2*phi**2*(1+4*iv.pi*R**2*H)
    Vbound+=8*iv.pi*alpha*R*phi*(E+3*P)+16*iv.pi*alpha*R*phi*H+4*iv.pi*(4+4*alpha**2)*(E+3*P)
    Vbound+=4*iv.pi*(alpha*(H+3*gp)+R*phi*(H+gp))*(3*alpha+R*phi)
    bmin=1-2*M/R;Cs=I(bg.Nb)*iv.sqrt(bmin);K=I(bg.K*bg.R)
    Ivac=(M/R**2+2*K*K/(3*Cs**2*R**3))/iv.sqrt(bmin)
    integral=-qmin*Vbound+Ivac;eta=I(C)*T/2*integral
    assert up(eta)<1
    M0=I(M0J);source_time=I(float(np.load(old.evolution.OUT/f'cells-{n}.npz')['history'][-1,0]));portTV=I(float(np.sum(abs(np.diff(pb)),dtype=np.longdouble)))
    norms={}
    for channel in ['trace','stress']:
        packet=packets[(n,channel)];wb=abs(packet['weights']);mb=np.max(abs(packet['baryon']),axis=0);en=np.max(abs(packet['nonrest']),axis=0)
        ab=I(0);an=I(0)
        for w,bv,ev in zip(wb,mb,en):ab+=I(float(w))*I(float(bv));an+=I(float(w))*I(float(ev))
        body_w=alpha/rmin if channel=='trace' else phi
        ab+=body_w*portTV;an+=body_w*portTV*I(abs(rows[channel]['bulk']['trace_nonrest_erg_per_g']))
        z=I(G)/(2*I(C))*I(packet['cx'])*source_time*ab+I(G)/(2*I(C)**3)*source_time*an
        norms[channel]=up(z);M0+=z
    # Small input-reduction guard covers the binary64/long-double maxima and
    # port total-variation reductions before the interval arithmetic begins.
    M0*=I('1.000000001')
    full=eta/(1-eta)*M0;remainder=eta**2/(1-eta)*M0
    central=I(C)*T/2*(I(qhi)-I(qlo))*Vbound*M0
    row=dict(classification='Counterexample candidate',conditional=True,
        theorem='Retarded-Green contraction and Neumann remainder under the stated coefficient and source bounds.',
        potential_absolute_integral_bound_per_cm=up(integral),contraction_bound=up(eta),free_field_norm_bound_cm=up(M0),free_field_component_bounds_cm=norms,
        all_orders_change_bound_cm=up(full),higher_Born_terms_bound_cm=up(remainder),omitted_compact_first_Born_bound_cm=up(central),
        all_orders_change_normalized_bound=up(full/M),higher_Born_normalized_bound=up(remainder/M),compact_normalized_bound=up(central/M),
        coefficient_enclosure='On the causal inner shell use extrema of the frozen piecewise-linear background, N<=1, A<=1, m<=M and r>=R+qmin. Outside use the exact positive-mass vacuum, C>=C_surface and Phi=K/(C*r^2).',
        source_enclosure='Piecewise-linear saved source histories, the declared local acoustic port profile and prescribed forcing J. The norm uses absolute source envelopes, not cancellation in the observed charge.',
        exclusions='Additional material displacement feedback, nonlinear metric terms, physical EOS/chemistry, unknown dilute-tail stress and causal scattered-photon forcing are not enclosed by this potential bound.',
        first_Born_quadrature_certified=False,all_physical_charge_enclosed=False)
    write(OUT/f'bound-{n}.json',row)
    return row


def prepare():
    assert not OUT.exists();OUT.mkdir()
    for channel in ['trace','stress']:
        for direction in [1,-1]:(OUT/f'{channel}-{direction}').mkdir()
    paths=[Path(__file__),Path(old.__file__),Path(prior.__file__),old.evolution.old.OUT/'runtime-columns.npz',prior.old.prior.OUT/'background.npz']
    paths += [old.evolution.OUT/f'cells-{n}.npz' for n in [896,1792]]
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='2282eb421',
        claim='Apply conserved-volume mass/scalar feedback, compute outgoing and incoming prescribed waves from the saved conservative release, and apply the curved-space scalar potential with an all-orders contraction bound.',
        closure='Linear adiabatic metric-volume response at fixed coordinate baryon inventory. No claim that independently prescribed fluid displacement or scattered photons are solved.',
        seconds=120,CPU_threads=1,memory_GB=2,native_calls=0,new_fluid_steps=0,
        measured_basis='Phase102 four direct component readouts plus two source integrations took24.07seconds. Eight directional readouts, two conserved source integrations and potential quadratures are budgeted at120seconds; the first complete grid supplies a stop-before-second-grid forecast.',
        settings=dict(saved_gas_grids=[1792,896],potential_inner_cells=[64,128],vacuum_quadrature=[32,64],Born_terms_computed=1),
        gates=dict(grid_wave=.02,first_Born_quadrature=.001,all_orders_potential_bound_over_free_wave=.02,manufactured_bound=True),
        stop='No new EOS, gas evolution, duration, spatial grids or extra Born terms on failure. Retain failed results and reassess before changing the plan.',
        bindings={str(p.relative_to(prior.old.ROOT)):prior.old.photons.digest(p) for p in paths},symbolic=symbolic()))


def run():
    assert not (OUT/'result.json').exists();start=time.monotonic();signal.alarm(120)
    plan=json.loads((OUT/'plan.json').read_text())
    for p,h in plan['bindings'].items():assert prior.old.photons.digest(prior.old.ROOT/p)==h,p
    original=np.load(old.OUT/'components.npz')['u_seconds'];m=old.adapter(1792);delay=prior.Geometry(m)(np.array([-28000.,120000.]))[1]
    utimes=np.unique(np.r_[np.linspace(-1.1*max(abs(delay)),0,9),original]);raw=np.load(prior.OUT/'native-bulk.npz')['raw'];packets={};results={};waveforms={}
    for n in [1792,896]:
        began=time.monotonic();directional=np.zeros((2,len(utimes)));rows={}
        for channel in ['trace','stress']:
            for i,direction in enumerate([1,-1]):
                alpha,row=readout(n,channel,direction,utimes,raw,packets);directional[i]+=-alpha*m.bg.M*m.R
                if direction==1:rows[channel]=row
        model,jwaves,M0J,pb,jchange=mass_wave(n,utimes);directional+=jwaves
        coarse,cq=potential_correction(model,utimes,directional,64,32);fine,fq=potential_correction(model,utimes,directional,128,64)
        bnd=bound(model,utimes,packets,n,rows,M0J,pb,fq);M=model.bg.M*model.R
        quaderr=float(max(abs(coarse-fine))/max(max(abs(fine)),1e-100));peak=float(max(abs(directional[0])))
        potential_bound=bnd['all_orders_change_bound_cm']/peak
        estimate=directional[0]+fine;waveforms[n]=-estimate/M
        row=dict(classification='Counterexample candidate',cells=n,seconds=time.monotonic()-began,
            first_Born_quadrature_relative=quaderr,potential_all_orders_bound_over_free_wave=potential_bound,
            endpoint_free_normalized=float(-directional[0,-1]/M),endpoint_first_Born_normalized=float(-fine[-1]/M),
            endpoint_corrected_estimate_normalized=float(-estimate[-1]/M),conserved_J_change_relative=jchange,
            initial_outgoing_history_extension_peak=float(max(abs(directional[0,utimes<0]))/peak),
            coarse_quadrature=cq,fine_quadrature=fq,passed=quaderr<.001 and potential_bound<.02)
        np.savez_compressed(OUT/f'wave-{n}.npz',u_seconds=utimes,outgoing_free_cm=directional[0],incoming_free_cm=directional[1],
            first_Born_cm=fine,corrected_estimate_cm=estimate,normalized_estimate=-estimate/M)
        for channel in ['trace','stress']:np.savez_compressed(OUT/f'{channel}-source-{n}.npz',**packets[(n,channel)])
        write(OUT/f'wave-{n}.json',row);results[n]=row;print('WAVE',n,json.dumps(row),flush=True)
        if n==1792:
            elapsed=time.monotonic()-start;forecast=row['seconds']*1.3+5
            write(OUT/'measured-budget.json',dict(first_grid_seconds=row['seconds'],elapsed_seconds=elapsed,second_grid_forecast_seconds=forecast,remaining_seconds=120-elapsed,assumption='The half-size saved source with identical quadrature is allowed130percent of the measured first-grid cost plus5seconds.'))
            assert elapsed+forecast<120,'Remaining wave work exceeds budget'
    grid=float(max(abs(waveforms[896]-waveforms[1792]))/max(abs(waveforms[1792])))
    # Independently known delta-potential solution tests the bound beyond a
    # vanishing-potential limit: U=t^2-k*integral(U), k*T=0.2.
    mp.mp.dps=60;k=mp.mpf('.2');exact=2/k-2*(-mp.expm1(-k))/k**2;born=1-k/3
    manufactured=abs(exact-born)<=k*k/(1-k)
    result=dict(classification='Counterexample candidate',passed=all(r['passed'] for r in results.values()) and grid<.02 and bool(manufactured),
        source_grid_relative=grid,endpoint_corrected_estimate_normalized=results[1792]['endpoint_corrected_estimate_normalized'],
        endpoint_first_Born_normalized=results[1792]['endpoint_first_Born_normalized'],
        potential_all_orders_bound_over_free_wave=results[1792]['potential_all_orders_bound_over_free_wave'],
        seconds=time.monotonic()-start,manufactured_resolvent_bound_passed=bool(manufactured),
        conserved_volume_metric_feedback_applied=True,prescribed_source_potential_bounded_to_all_orders=True,
        full_material_displacement_feedback=False,causal_scattered_photon_metric_source_applied=False,
        final_charge_solved=False,full_goal_complete=False)
    write(OUT/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
