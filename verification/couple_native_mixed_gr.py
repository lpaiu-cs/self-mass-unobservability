"""Counterexample candidate: apply the missing reciprocal compact GR source.

Use the existing physical material increments once. This computes a specified
bilinear source on the initial canonical operator, not nonlinear transport.
"""
from pathlib import Path
import json, resource, sys, time
import numpy as np
import sympy as sp
from numpy.polynomial import legendre as leg
import couple_native_interior_incident as previous

OUT=Path('native-mixed-gr161-work')
wave=previous.wave; inf=previous.infinity; inc=inf.incident
C=inf.C; G=inf.G; LD=np.longdouble
read,write,sha=previous.read,previous.write,previous.sha
CAPS=dict(prepare=45,apply=240,audit=60)
SETTINGS=[(128,8), (128,4), ('128-17',8)]
FACTOR=2.**63


def symbolic():
    e,d=sp.symbols('e d');r,b,a,N,A,alpha,beta,E,P,T=sp.symbols('r b a N A alpha beta E P T',nonzero=True)
    fb,fi,lb,li,nb,ni,eb,ei,pb,pi,tb,ti,frb,fri,ftb,fti,c=sp.symbols('fb fi lb li nb ni eb ei pb pi tb ti frb fri ftb fti c')
    f=e*fb+d*fi; l=e*lb+d*li; n=e*nb+d*ni
    def cross(v):return sp.diff(v,e,d).subs({e:0,d:0})
    common=sp.exp(2*l+4*alpha*f+2*beta*f*f)
    kB=2*lb+4*alpha*fb;kI=2*li+4*alpha*fi
    k=kB*kI+4*beta*fb*fi
    for sign,value,vb,vi in [(-1,E,eb,ei),(1,P,pb,pi)]:
        exact=sign*sp.exp(2*l)/(2*r*b)+4*sp.pi*r*A/b*common*(value+e*vb+d*vi)
        exact+=r/2*(e*frb+d*fri)**2+r/(2*c*c*a*a)*(e*ftb+d*fti)**2
        known=sign*2*lb*li/(r*b)+4*sp.pi*r*A/b*(value*k+kB*vi+kI*vb)+r*frb*fri+r/(c*c*a*a)*ftb*fti
        assert sp.simplify(cross(exact)-known)==0
    hb=2*nb+4*alpha*fb;hi=2*ni+4*alpha*fi
    h=(alpha*(hb*hi+4*beta*fb*fi)+beta*(hb*fi+hi*fb))*T
    h+=(alpha*hb+beta*fb)*ti+(alpha*hi+beta*fi)*tb
    exact=sp.exp(2*n+4*alpha*f+2*beta*f*f)*(alpha+beta*f)*(T+e*tb+d*ti)
    assert sp.simplify(cross(exact)-h)==0
    # Static-gradient cross term after the two first metric actions are removed.
    zb,zi,zrb,zri,L0,Phi=sp.symbols('zb zi zrb zri L0 Phi')
    exact=sp.exp(2*(e*zb+d*zi))*(L0+(e*zrb+d*zri)*Phi)
    assert sp.simplify(cross(exact)-4*zb*zi*L0-2*(zb*zri+zi*zrb)*Phi)==0
    ql,qn=sp.symbols('ql qn')
    ef=b*ql/(4*sp.pi*r*A);pf=b*qn/(4*sp.pi*r*A)
    assert sp.simplify(4*sp.pi*r*r*A*ef-r*b*ql)==0
    assert sp.simplify(-4*sp.pi*r*N*N*A*r*Phi*(ef-pf)-r*N*N*b*Phi*(qn-ql))==0
    return dict(classification='Proven',passed=True,
        scope='Mixed derivatives of the declared polar Einstein constraints and scalar trace factor, alpha_prime=beta=-4; the two known first variations are held independent. Not a closure of second-order material evolution.',
        canonical='Recover physical dE,dP,dT by subtracting the initial canonical response ONCE from the saved forcing. Its actual noncanonical material response is already present.',
        mass='Jx_prime+(nu0_prime+lambda0_prime)*Jx=r*b*q_lambda. Jx=sqrt(b)/N integral N*r*sqrt(b)*q_lambda dr, with zero additional inner port.',
        scalar='Keff*Jx+r*a^2*Phi*(q_nu-q_lambda)+S_trace+S_static_gradient+S_incident_metric[U_background]. The previous background-metric action on the primary is not added again.',
        physical_mass='delta_m_mixed=r*b*delta_lambda_mixed-2*r*b*lambda_background*lambda_incident. The auxiliary constraint Jx is not a new matter energy inventory.')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    files=[Path(__file__),Path(previous.__file__),Path(wave.__file__),Path(wave.base.__file__),Path(inf.__file__),Path(inc.__file__),
           wave.base.flow.INPUT/'balanced-20.npz',previous.OUT/'applied-charge.npz',inf.BACKGROUND,
           inf.BEFORE/'expanded-incident-source.py',inf.BEFORE/'gr/source-128-reference-128.npz',
           inf.SELF/'sweep-2/gr/source-128-reference-128.npz',inf.ROOT/'def-native-global-scalar-closure/bound.json']
    for n,q in SETTINGS:
        files += [previous.FIELDS/f'source-{n}.npz',previous.FIELDS/f'fields-{n}-g{q}.npz',previous.OUT/f'operator-{n}-g{q}.npz']
    for q in [4,8]:
        files += [inc.FIELDS/f'born-g{q}.npz']
        for folder,stem in [(inf.SELF/'fields','fields'),(inf.SELF/'metric','metric'),
                            (inf.SELF/'sweep-2/returned-fields','fields'),(inf.SELF/'sweep-2/returned-metric','metric')]:
            files.append(folder/f'{stem}-128-g{q}.npz')
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='5e6501e94',
        claim='Actually propagate the omitted reciprocal compact bilinear Einstein-scalar source and retain the induced field for material/photon return.',
        decision='Does this missing gravitational coupling materially change the stored body charge, and which new field must be returned to transport? A failed control remains failed.',
        background='The exact same Phase160 current-source composition and saved4/8 fields. No background or EOS replay.',
        incident='Existing primary evaluated analytically at source knots, saved first Born, Phase156 actual matter response, and Phase157 physical body and normalized self-return. Divide stored self-return by2^63 exactly. Restore physical stress from initial canonical terms once.',
        representation='Bilinear nodal source on existing531 cells and33 knots, linear source-time interpolation and inverse-Legendre spatial density.17 knots are a source interpolation control, not independent fluid evolution. No new driving waveform is supplied to matter. Higher product interpolation error remains explicit.',
        derivatives='Compute background Uxx by differentiating the exact retarded piecewise-linear source including its saved potential iterate; use Uxx=Utt/c^2-S. No finite differences of tiny sampled waveforms.',
        terms=['Mixed lambda/nu constraints including scalar radial and time stress, propagated through Jx and the scalar stress source.',
               'Mixed conformal trace factor and static-gradient source.',
               'Incident metric acting on the generated background scalar.'],
        exclusions=['Exterior mixed source, mixed physical ADM/arrived-mass readout and interface completion.',
                    'Transport return of the newly generated field and complete mixed canonical material Hessian.',
                    'Previous failed upstream controls, body/input interpolation error, complete EOS/boundary/nonlinear/observational closure.'],
        gates=dict(quadrature=.002,cadence=.02,algebra=1e-10,observer=1e-10,potential=.01),
        budgets=CAPS,total_seconds=sum(CAPS.values()),CPU_threads=1,virtual_GiB=3,new_fluid_steps=0,new_EOS_roots=0,
        pilot='Measure the first fine route including analytic second derivative and propagate. Twice that cost for each remaining route plus15s must fit the remaining240s. Reuse it.',
        stop='Do not automatically enlarge time interval, mesh, source clock, iteration budget or gates. Preserve raw outputs and stop after registered controls.',
        reference='https://arxiv.org/abs/gr-qc/9707041; mixed identities are derived here from the existing repository convention, not imported as a result from that paper.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)}))
    write(OUT/'symbolic.json',symbolic());(OUT/'registered-producer.py').write_bytes(Path(__file__).read_bytes())


def interp(t,clock,v):
    j=np.clip(np.searchsorted(clock,t,side='right')-1,0,len(clock)-2)
    w=((np.asarray(t)-clock[j])/np.diff(clock)[j]).reshape((-1,)+(1,)*(v.ndim-1))
    return (1-w)*v[j]+w*v[j+1]


def node_field(m,field,t,key):
    value=interp(t,field['t'],field[key])
    return np.array([np.interp(m.r,field['radius_E'],v) for v in value])


def physical(m,d,t,f,lam):
    z=m.z
    rest=np.asarray(d['baryon_g'],LD)*LD(d['cx'])*LD(C)**2
    energy=rest+d['gas_nonrest_energy_erg']+d['photon_energy_erg']
    values=[energy,energy-d['metric_stress_erg'],rest+d['nonrest_trace_erg']]
    e,p,tr=[np.asarray(interp(t,d['t'],v)/d['volume']*LD(G)/LD(C)**4,float)[:,m.ids] for v in values]
    volume=3*z['alpha']*f+lam;H=z['Eg']+z['Pg']
    return (e-H*volume-4*z['alpha']*z['Er']*f-(z['Er']+z['Pr'])*lam,
            p-z['Kg']*volume-4*z['alpha']*z['Pr']*f-(3*z['Pr']-z['R4'])*lam,
            tr-(H-3*z['Kg'])*volume)


def second_derivative(m,source):
    """Exact Uxx of the declared retarded linear-time polynomial density."""
    density=source/m.dx;slopes=np.diff(density,axis=0)/np.diff(m.t)[:,None]
    co=np.einsum('cij,tcj->tci',m.inverse,slopes.reshape(len(m.t)-1,-1,m.order))
    initial=np.einsum('cij,cj->ci',m.inverse,density[0].reshape(-1,m.order))
    def evaluate(coeff,x):
        ids=np.clip(np.searchsorted(m.xfaces,x,side='right')-1,0,len(m.mid)-1)
        vv=leg.legvander((x-m.mid[ids])/m.half[ids],m.order-1)
        value=np.sum(coeff[ids]*vv,axis=-1)
        return np.where((x>=m.xfaces[0])&(x<=m.xfaces[-1]),value,0.)
    result=[];cols=np.arange(len(m.ids))[None,:]
    for t in m.t:
        ret=t-m.distance;idx=np.clip(np.searchsorted(m.t,ret,side='right')-1,0,len(m.t)-2)
        value=slopes[idx,cols]*m.dx;value[ret<=0]=0
        cell_values=value.reshape(len(m.tx),-1,m.order).sum(2)
        owner,cell,lo,hi=m.segments(t)
        if len(cell):
            xx=(lo[:,None]+hi[:,None])/2+(hi-lo)[:,None]*m.gx/2
            rr=t-abs(m.tx[owner,None]-xx)/C
            jj=np.clip(np.searchsorted(m.t,rr,side='right')-1,0,len(m.t)-2)
            vv=leg.legvander((xx-m.mid[cell,None])/m.half[cell,None],m.order-1)
            value=np.sum(co[jj,cell[:,None]]*vv,axis=-1);value[rr<=0]=0
            affected=np.unique(np.column_stack([owner,cell]),axis=0);ii,jj=affected.T;cell_values[ii,jj]=0
            np.add.at(cell_values,(owner,cell),np.sum((hi-lo)[:,None]*m.gw/2*value,axis=1))
        utt_over_c2=cell_values.sum(1,dtype=LD)/(2*C)
        ti=np.clip(np.searchsorted(m.t,t,side='right')-1,0,len(m.t)-2)
        nodal=density[ti]+(t-m.t[ti])*slopes[ti]
        coeff=np.einsum('cij,cj->ci',m.inverse,nodal.reshape(-1,m.order))
        fronts=(evaluate(initial,m.tx-C*t)+evaluate(initial,m.tx+C*t))/2
        result.append(utt_over_c2+fronts-evaluate(coeff,m.tx))
    return np.asarray(result,float)


def fields_and_stress(m,d,field,operator,qorder):
    t=d['t'];r=m.r;z=m.z;a=z['lapse']*np.sqrt(z['b']);q=wave.base.flow.initial.Quadrature(d['edges'],qorder)
    jac=1/(np.exp(-2*z['phi']**2)*(1+z['alpha']*r*z['Phi']))
    B=dict(f=node_field(m,field,t,'delta_phi'),fr=node_field(m,field,t,'delta_Phi'),
           ft=node_field(m,field,t,'U_t')/r,U=node_field(m,field,t,'U'),Ut=node_field(m,field,t,'U_t'),Ux=node_field(m,field,t,'U_x'))
    B['lam']=r*z['Phi']*B['f']+m.J/(r*z['b']);B['zet']=operator['zeta'];B['zr']=operator['zeta_prime']
    B['nu']=B['lam']+B['zet'];B['de'],B['dp'],B['dt']=physical(m,d,t,B['f'],B['lam'])
    assert len(field['potential_iterations'])==1
    bg_source=m.source-m.dx*z['V']*field['direct_and_mass_stress_U'][:,:-1][:,m.ids]
    bx=second_derivative(m,bg_source)
    B['Uxx']=np.array([np.interp(r,field['radius_E'],v) for v in bx])
    data={k:v.copy() for k,v in np.load(inf.BEFORE/'gr/source-128-reference-128.npz').items()}
    ds=np.load(inf.SELF/'sweep-2/gr/source-128-reference-128.npz')
    for k in previous.current.KEYS:data[k]=data[k]+ds[k]/LD(FACTOR)
    fs={k:v.copy() for k,v in np.load(inf.SELF/f'fields/fields-128-g{qorder}.npz').items()}
    extra=np.load(inf.SELF/f'sweep-2/returned-fields/fields-128-g{qorder}.npz')
    for k in ['U','U_t','U_x','delta_phi','delta_Phi']:fs[k]+=extra[k]/FACTOR
    cs=previous.retained.constraints.centers(m,data,fs,qorder)
    I=dict(f=node_field(m,fs,t,'delta_phi'),fr=node_field(m,fs,t,'delta_Phi'),ft=node_field(m,fs,t,'U_t')/r)
    J=interp(t,data['t'],m.J);I['lam']=r*z['Phi']*I['f']+J/(r*z['b'])
    met=np.load(inf.SELF/f'metric/metric-128-g{qorder}.npz')
    metx=np.load(inf.SELF/f'sweep-2/returned-metric/metric-128-g{qorder}.npz')
    outer=interp(t,met['t'],met['delta_nu_faces'][:,-1]+metx['delta_nu_faces'][:,-1]/FACTOR-cs['delta_lambda'][:,-1])
    stress=np.asarray(interp(t,data['t'],data['metric_stress_erg'])/data['volume']*LD(G)/LD(C)**4,float)[:,m.ids]
    kernel_l=2/(r*z['b'])+4*np.pi*r*z['A4']/z['b']*(3*z['Pg']-z['Eg']-z['Kg']-z['Er']+z['R4'])
    kernel_f=4*np.pi*r*z['A4']/z['b']*z['alpha']*(7*z['Pg']-z['Eg']-3*z['Kg'])
    I['zr']=kernel_l*I['lam']+kernel_f*I['f']-4*np.pi*r*z['A4']/z['b']*stress
    I['zet']=np.array([(edge-(q.integrate(v*jac)[1][-1]-q.integrate(v*jac)[0])).ravel() for edge,v in zip(outer,I['zr'])])
    slopes=np.diff(I['zet'],axis=0)/np.diff(t)[:,None]
    I['zt']=np.r_[slopes,slopes[-1:]]
    driver=inc.Driver(qorder);x=m.x-m.xfaces[-1];kernel=kernel_l*r*z['Phi']+kernel_f
    for it,now in enumerate(t):
        U,Ut,Ux=driver.wave(now,x);f=U/r;ft=Ut/r;fr=Ux/(a*r)-U/r**2
        ue,uet,_=driver.wave(now,driver.xout)
        def speed(value,exterior_value):
            ext=np.sum(driver.eq.h*((2*driver.zout['Phi']*exterior_value/(driver.rout*driver.zout['b'])).reshape(driver.eq.r.shape)@driver.eq.w),dtype=LD)
            integ,faces=q.integrate(kernel*value*jac)
            return np.asarray(-ext-(faces[-1]-integ),float).ravel()
        I['f'][it]+=f;I['ft'][it]+=ft;I['fr'][it]+=fr;I['lam'][it]+=r*z['Phi']*f
        I['zet'][it]+=speed(f,ue);I['zt'][it]+=speed(ft,uet);I['zr'][it]+=kernel*f
    I['nu']=I['lam']+I['zet'];I['de'],I['dp'],I['dt']=physical(m,data,t,I['f'],I['lam'])
    m.setup(d,qorder)
    return B,I,q,jac,bg_source,bx


def sources(m,B,I,q,jac):
    z=m.z;r=m.r;a=z['lapse']*np.sqrt(z['b']);alpha=z['alpha'];beta=-4
    E=z['Eg']+z['Er'];P=z['Pg']+z['Pr'];T=z['Eg']-3*z['Pg']
    kb=2*B['lam']+4*alpha*B['f'];ki=2*I['lam']+4*alpha*I['f'];cross=kb*ki+4*beta*B['f']*I['f']
    scalar=r*B['fr']*I['fr']+r/(C*C*a*a)*B['ft']*I['ft']
    radial=2/(r*z['b'])*B['lam']*I['lam'];fac=4*np.pi*r*z['A4']/z['b']
    ql=-radial+fac*(E*cross+kb*I['de']+ki*B['de'])+scalar
    qn=radial+fac*(P*cross+kb*I['dp']+ki*B['dp'])+scalar
    J=np.array([(np.sqrt(z['b'])/z['lapse']*q.integrate(a*r*row*jac)[0].ravel()) for row in ql])
    hb=2*B['nu']+4*alpha*B['f'];hi=2*I['nu']+4*alpha*I['f']
    trace=(alpha*(hb*hi+4*beta*B['f']*I['f'])+beta*(hb*I['f']+hi*B['f']))*T
    trace+=(alpha*hb+beta*B['f'])*I['dt']+(alpha*hi+beta*I['f'])*B['dt']
    trace=-4*np.pi*r*z['lapse']**2*z['A4']*trace
    gradient=4*B['zet']*I['zet']*r*4*np.pi*z['lapse']**2*z['A4']*alpha*T
    gradient+=2*a*a*r*z['Phi']*(B['zet']*I['zr']+I['zet']*B['zr'])
    geom=z['lapse']**2*(2*z['mass']/r**3+4*np.pi*z['A4']*(P-E))
    reciprocal=2*I['zet']*B['Uxx']+a*I['zr']*B['Ux']+(-2*I['zet']*geom-a*a*I['zr']/r)*B['U']+I['zt']*B['Ut']/C**2
    parts=dict(trace=trace,static_gradient=gradient,reciprocal=reciprocal,
               constraint=z['K']*J+r*a*a*z['Phi']*(qn-ql))
    return parts,J,ql,qn


def route(n,order):
    mark=time.monotonic();d=dict(np.load(previous.FIELDS/f'source-{n}.npz'));f=np.load(previous.FIELDS/f'fields-{n}-g{order}.npz')
    op=np.load(previous.OUT/f'operator-{n}-g{order}.npz');m=wave.Response();m.setup(d,order)
    B,I,q,jac,bg_source,bx=fields_and_stress(m,d,f,op,order)
    parts,J,ql,qn=sources(m,B,I,q,jac);source=m.dx*sum(parts.values())
    free=m.propagate(source);U=free[0];history=[]
    contraction=float(C*m.t[-1]/2*np.sum(m.dx*abs(m.z['V']),dtype=LD));assert contraction<.01
    for _ in range(4):
        potential=-m.dx*m.z['V']*U[:,:-1][:,m.ids];add=m.propagate(potential);new=free[0]+add[0]
        history.append(float(np.max(abs(new-U))/max(np.max(abs(new)),1e-290)));U=new
        if history[-1]<1e-8:break
    else:raise AssertionError(('potential iteration',history))
    Ut=free[1]+add[1];Ux=free[2]+add[2];r=m.tr;z=m.tz
    np.savez_compressed(OUT/f'field-{n}-g{order}.npz',t=m.t,radius_E=r,U=U,U_t=Ut,U_x=Ux,
        delta_phi=U/r,delta_Phi=Ux/(r*z['lapse']*np.sqrt(z['b']))-U/r**2,
        background_Uxx=bx,source=source,potential_source=potential,source_dx=m.dx,source_x=m.x-m.xfaces[-1],
        xfaces=m.xfaces-m.xfaces[-1],inverse=m.inverse,J_mixed=J,q_lambda=ql,q_nu=qn,
        **{'density_'+k:v for k,v in parts.items()})
    clock=np.load(previous.OUT/'applied-charge.npz')['t'];origin=m.xfaces[-1]
    m.x-=origin;m.xfaces-=origin;m.mid-=origin
    charges={k:-inf.observe(m,m.dx*v,clock,3*C*clock[-1])/float(d['M_cm']) for k,v in parts.items()}
    charges['potential']=-inf.observe(m,potential,clock,3*C*clock[-1])/float(d['M_cm'])
    charge=sum(charges.values());other=-inf.observe(m,source+potential,clock,4*C*clock[-1])/float(d['M_cm'])
    observer_error=float(np.max(abs(charge-other))/max(np.max(abs(charge)),1e-290));assert observer_error<1e-10
    np.savez_compressed(OUT/f'charge-{n}-g{order}.npz',t=clock,total=charge,other_observer=other,**charges)
    return dict(source=n,order=order,seconds=time.monotonic()-mark,endpoint=float(charge[-1]),
        components={k:float(v[-1]) for k,v in charges.items()},maximum_phi=float(np.max(abs(U/r))),
        observer_relative=observer_error,compact_potential_contraction=contraction,potential_iteration=history)


def apply():
    previous.initialize();start=time.monotonic();rows=[route(*SETTINGS[0])]
    remaining=4*rows[0]['seconds']+15
    write(OUT/'pilot.json',dict(first_seconds=rows[0]['seconds'],remaining_upper_seconds=remaining,eligible=remaining<CAPS['apply']-(time.monotonic()-start)))
    assert remaining<CAPS['apply']-(time.monotonic()-start)
    rows += [route(*setting) for setting in SETTINGS[1:]]
    fine=np.load(OUT/'charge-128-g8.npz')['total'];scale=max(np.max(abs(fine)),1e-290)
    errors={k:float(np.max(abs(np.load(OUT/f'charge-{n}-g{q}.npz')['total']-fine))/scale)
            for k,(n,q) in zip(['quadrature','cadence'],SETTINGS[1:])}
    old=np.load(previous.OUT/'applied-charge.npz');epsilon=np.load(inf.BACKGROUND)['epsilon'].astype(LD)
    contribution=fine.astype(LD)/(1-epsilon);combined=old['body_plus_selected_return']+contribution
    np.savez_compressed(OUT/'applied-charge.npz',t=old['t'],old_body=old['old_body'],previous_selected=old['body_plus_selected_return'],
        reciprocal_compact_return=contribution,body_plus_selected_returns=combined,direct=old['direct'])
    result=dict(classification='Counterexample candidate',passed=errors['quadrature']<.002 and errors['cadence']<.02,
        controls=errors,endpoint_reciprocal_compact_return=float(contribution[-1]),endpoint_previous_selected=float(old['body_plus_selected_return'][-1]),
        endpoint_combined_selected=float(combined[-1]),return_over_original_body=float(contribution[-1]/old['old_body'][-1]),rows=rows,
        reciprocal_compact_source_applied=True,new_wave_applied_to_material_photons=False,exterior_mixed_source_closed=False,
        mixed_ADM_normalization_closed=False,full_goal_complete=False,seconds=time.monotonic()-start)
    write(OUT/'result.json',result);print(json.dumps(result),flush=True)
    assert result['passed'],errors


def repair():
    failed=read(OUT/'apply-receipt.json');assert failed['error']=="AssertionError('Nonzero initial source needs wavefront terms.')"
    assert not (OUT/'repair-plan.json').exists()
    write(OUT/'repair-plan.json',dict(classification='Counterexample candidate',
        failure='The actual initial source is not identically zero. The guarded derivative path stopped before propagation; do not erase or zero that initial source.',
        correction='Add the exact initial fronts [S0(x-c*t)+S0(x+c*t)]/2 to Uxx. Keep the same background, pulse, clocks, grids and acceptance gates.',
        original_sha256=sha(OUT/'registered-producer.py'),producer_sha256=sha(__file__),
        audit_sha256=sha(Path(__file__).with_name('audit_native_mixed_gr.py')),
        remaining_apply_seconds=CAPS['apply']-failed['seconds']))
    from audit_native_mixed_gr import derivative_check
    write(OUT/'derivative-check.json',derivative_check())
    inc.native.deadline(CAPS['apply']-failed['seconds']);apply()


if __name__=='__main__':
    action=sys.argv[1];assert action in ['prepare','apply','repair'];receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));inc.native.deadline(CAPS.get(action,CAPS['apply']))
    started=time.monotonic();cpu=time.process_time();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():
                assert sha(OUT/'registered-producer.py' if action=='repair' and p==str(Path(__file__)) else p)==h,p
        globals()[action]()
    except Exception as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(action=action,seconds=time.monotonic()-started,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
