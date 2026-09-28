"""Resolve interpolation cuts in the same mixed source before transport.

Counterexample candidate. No extra physical clock, cells or source history.
"""
from pathlib import Path
import copy, json, time
import numpy as np
from numpy.polynomial import legendre as leg
from scipy.special import betainc
from scipy.integrate import quad
import return_native_mixed_gr as run

mixed=run.mixed;C=run.C;LD=run.LD;OUT=run.OUT/'resolved';DEST=run.OUT/'corrected-fields'


def pulse_propagate(m,coefficient,amplitude,duration):
    """Green integral of K(t,x)*g((t+x/c)/D), exact declared polynomials.

    K is the ORIGINAL inverse-Legendre spatial / linear-time coefficient.
    Additional nodes integrate its product; they are not new source samples.
    """
    t=m.t;order=m.order;origin=m.xfaces[-1];oldfaces=m.xfaces-origin
    oldmid=m.mid-origin;oldhalf=m.half;targets=m.tx-origin
    faces=np.unique(np.r_[oldfaces,-C*t,C*(duration-t)])
    faces=faces[(faces>=oldfaces[0])&(faces<=0)]
    mid=(faces[:-1]+faces[1:])/2;half=np.diff(faces)/2
    degree=order+10;gx,gw=leg.leggauss(degree);inverse=np.linalg.inv(leg.legvander(gx,degree-1))
    x=(mid[:,None]+half[:,None]*gx).ravel();weights=(half[:,None]*gw).ravel();nc=len(mid)
    owners=np.clip(np.searchsorted(oldfaces,x,side='right')-1,0,len(oldmid)-1)
    vc=leg.legvander((x-oldmid[owners])/oldhalf[owners],order-1)
    original=np.einsum('cij,tcj->tci',m.inverse,coefficient.reshape(len(t),-1,order))
    values=np.einsum('tij,ij->ti',original[:,owners],vc)
    kc=(values.reshape(len(t),nc,degree)@inverse.T)[...,:order]
    slopes=np.diff(values,axis=0)/np.diff(t)[:,None]
    def primitive(y):
        y=np.clip(y,0,1)
        return 128/315*betainc(5,5,y),64/315*betainc(6,5,y)
    def interval(A,B,start,end,xx):
        p0,p1=primitive((start+xx/C)/duration);q0,q1=primitive((end+xx/C)/duration)
        return amplitude*duration*((A-B*(start+xx/C))*(q0-p0)+B*duration*(q1-p1))
    increments=interval(values[:-1],slopes,t[:-1,None],t[1:,None],x)
    H=np.r_[np.zeros((1,len(x))),np.cumsum(increments,axis=0)]
    hc=H.reshape(len(t),nc,degree)@inverse.T
    def evaluate(ret,xx,cell=None):
        s=np.clip(ret,0,t[-1]);j=np.clip(np.searchsorted(t,s,side='right')-1,0,len(t)-2)
        if cell is None:
            cols=np.arange(len(x))[None,:];A=values[j,cols];B=slopes[j,cols];initial=H[j,cols]
        else:
            cell=np.broadcast_to(cell,xx.shape);v=leg.legvander((xx-mid[cell])/half[cell],degree-1)
            A=np.sum(kc[j,cell]*v[...,:order],axis=-1)
            B=np.sum((kc[j+1,cell]-kc[j,cell])*v[...,:order],axis=-1)/(t[j+1]-t[j])
            initial=np.sum(hc[j,cell]*v,axis=-1)
        value=initial+interval(A,B,t[j],s,xx)
        y=(s+xx/C)/duration;yy=np.clip(y,0,1)
        source=(A+B*(s-t[j]))*amplitude*256*yy**4*(1-yy)**4
        source[(ret<=0)|(y<=0)|(y>=1)]=0.;value[ret<=0]=0.
        return value,source
    outputs=[np.zeros((len(t),len(targets))) for _ in range(3)]
    for it,now in enumerate(t):
        for begin in range(0,len(targets),64):
            tx=targets[begin:begin+64];ret=now-abs(tx[:,None]-x)/C
            hv,sv=evaluate(ret,x[None,:]);dv=sv*np.sign(tx[:,None]-x)
            hsum=(hv*weights).reshape(len(tx),nc,degree).sum(2)
            ssum=(sv*weights).reshape(len(tx),nc,degree).sum(2);dsum=(dv*weights).reshape(len(tx),nc,degree).sum(2)
            owners=[];cells=[];left=[];right=[]
            for i,target in enumerate(tx):
                cuts=np.unique(np.r_[target,target-C*(now-t[t<=now]),target+C*(now-t[t<=now]),
                                      (target-C*now)/2,(target+C*(duration-now))/2])
                cuts=cuts[(cuts>faces[0])&(cuts<faces[-1])];ids=np.searchsorted(faces,cuts,side='right')-1
                for j in np.unique(ids):
                    edges=np.unique(np.r_[faces[j],cuts[ids==j],faces[j+1]])
                    owners.extend([i]*(len(edges)-1));cells.extend([j]*(len(edges)-1));left.extend(edges[:-1]);right.extend(edges[1:])
            if cells:
                owner=np.asarray(owners);cell=np.asarray(cells);lo=np.asarray(left);hi=np.asarray(right)
                xx=(lo[:,None]+hi[:,None])/2+(hi-lo)[:,None]*gx/2;rr=now-abs(tx[owner,None]-xx)/C
                hh,ss=evaluate(rr,xx,cell[:,None]);ww=(hi-lo)[:,None]*gw/2
                affected=np.unique(np.column_stack([owner,cell]),axis=0);ii,jj=affected.T
                hsum[ii,jj]=0.;ssum[ii,jj]=0.;dsum[ii,jj]=0.
                np.add.at(hsum,(owner,cell),np.sum(ww*hh,axis=1));np.add.at(ssum,(owner,cell),np.sum(ww*ss,axis=1))
                np.add.at(dsum,(owner,cell),np.sum(ww*ss*np.sign(tx[owner,None]-xx),axis=1))
            for out,value,factor in zip(outputs,[hsum,ssum,dsum],[C/2,C/2,-.5]):out[it,begin:begin+len(tx)]=factor*value.sum(1,dtype=LD)
    return tuple(outputs)


def check():
    from types import SimpleNamespace
    errors=[]
    for order in [4,8]:
        x,w=leg.leggauss(order);m=SimpleNamespace(order=order,t=np.array([0.,.2,.5,.8,1.])/C,
            xfaces=np.array([-1.,0.]),mid=np.array([-.5]),half=np.array([.5]),tx=np.array([-.75,-.3,.1]),
            inverse=np.linalg.inv(leg.legvander(x,order-1))[None])
        actual=pulse_propagate(m,np.ones((len(m.t),order)),1.,.5/C)
        for it,T in enumerate(C*m.t):
            for j,X in enumerate(m.tx):
                def f(y,which):
                    upper=T-abs(X-y);lo=max(0.,-y);hi=min(upper,.5-y)
                    if hi<=lo:return 0.
                    if which==0:return .25*(128/315)*(betainc(5,5,np.clip(2*(hi+y),0,1))-betainc(5,5,np.clip(2*(lo+y),0,1)))
                    z=2*(upper+y)
                    val=128*z**4*(1-z)**4 if 0<z<1 else 0.
                    return val if which==1 else -np.sign(X-y)*val
                points=[v for v in [X,(X-T)/2,(X+.5-T)/2,-T,.5-T] if -1<v<0]
                for which in range(3):
                    expected=quad(lambda y:f(y,which),-1,0,points=points,epsabs=1e-12,epsrel=1e-12)[0]
                    value=actual[which][it,j]/(C if which==1 else 1.)
                    errors.append(abs(value-expected))
    assert max(errors)<1e-10, max(errors)
    return dict(classification='Proven',passed=True,exact_pulse_box_max_absolute=max(errors),
        scope='Polynomial coefficient times the exact original pulse, with internal and external targets. Not physical source error certification.')


def constraint_context(m,d,field,order):
    di={k:v.copy() for k,v in np.load(mixed.inf.BEFORE/'gr/source-128-reference-128.npz').items()}
    extra=np.load(mixed.inf.SELF/'sweep-2/gr/source-128-reference-128.npz')
    for k in mixed.previous.current.KEYS:di[k]=di[k]+extra[k]/LD(mixed.FACTOR)
    fi={k:v.copy() for k,v in np.load(mixed.inf.SELF/f'fields/fields-128-g{order}.npz').items()}
    fx=np.load(mixed.inf.SELF/f'sweep-2/returned-fields/fields-128-g{order}.npz')
    for k in ['delta_phi','delta_Phi','U_t']:fi[k]+=fx[k]/mixed.FACTOR
    values=[]
    for source,fields in [(d,field),(di,fi)]:
        m.setup(source,order);J=mixed.interp(d['t'],source['t'],m.J)
        co=np.einsum('cij,tcj->tci',m.inverse,J.reshape(len(d['t']),-1,order))
        rest=np.asarray(source['baryon_g'],LD)*LD(source['cx'])*LD(C)**2
        energy=rest+source['gas_nonrest_energy_erg']+source['photon_energy_erg']
        raw=[np.asarray(mixed.interp(d['t'],source['t'],v)/source['volume']*LD(run.G)/LD(C)**4,float) for v in [energy,energy-source['metric_stress_erg']]]
        fv={k:mixed.interp(d['t'],fields['t'],fields[k]) for k in ['delta_phi','delta_Phi','U_t']}
        values.append(dict(J=co,raw=raw,field=fv,radius=fields['radius_E']))
    m.setup(d,order)
    return values


def constraints(m,context,driver,it,now,r,x,z,ids):
    radius=np.asarray(r);a=z['lapse']*np.sqrt(z['b']);A=z['alpha'];Phi=z['Phi'];H=z['Eg']+z['Pg']
    vander=leg.legvander((x-m.mid[ids])/m.half[ids],m.order-1);states=[]
    for k,ctx in enumerate(context):
        f=np.interp(radius,ctx['radius'],ctx['field']['delta_phi'][it])
        fr=np.interp(radius,ctx['radius'],ctx['field']['delta_Phi'][it]);ft=np.interp(radius,ctx['radius'],ctx['field']['U_t'][it])/radius
        if k:
            U,Ut,Ux=driver.wave(now,x-m.xfaces[-1]);f+=U/radius;fr+=Ux/(a*radius)-U/radius**2;ft+=Ut/radius
        J=np.sum(ctx['J'][it,ids]*vander,axis=-1);lam=radius*Phi*f+J/(radius*z['b']);volume=3*A*f+lam
        E=ctx['raw'][0][it,ids]-H*volume-4*A*z['Er']*f-(z['Er']+z['Pr'])*lam
        P=ctx['raw'][1][it,ids]-z['Kg']*volume-4*A*z['Pr']*f-(3*z['Pr']-z['R4'])*lam
        states.append(dict(f=f,fr=fr,ft=ft,lam=lam,E=E,P=P))
    B,I=states;kb=2*B['lam']+4*A*B['f'];ki=2*I['lam']+4*A*I['f'];cross=kb*ki-16*B['f']*I['f']
    scalar=radius*B['fr']*I['fr']+radius/(C*C*a*a)*B['ft']*I['ft'];radial=2/(radius*z['b'])*B['lam']*I['lam'];fac=4*np.pi*radius*z['A4']/z['b']
    ql=-radial+fac*((z['Eg']+z['Er'])*cross+kb*I['E']+ki*B['E'])+scalar
    qn=radial+fac*((z['Pg']+z['Pr'])*cross+kb*I['P']+ki*B['P'])+scalar
    return ql,qn


def resolved_constraints(m,d,field,order):
    context=constraint_context(m,d,field,order);driver=mixed.inc.Driver(order)
    q=run.wave.base.flow.initial.Quadrature(d['edges'],order);targets=np.r_[q.r.ravel(),d['radius'],d['edges'][-1]]
    rtargets=np.r_[m.r,m.tr];zt=m.coeff(rtargets);gx,gw=leg.leggauss(order);rows=[];qn=[];checks=[]
    for it,now in enumerate(d['t']):
        cuts=[targets,d['edges']]
        for xx in [-C*now,C*(driver.D-now)]:
            if m.xfaces[0]-m.xfaces[-1]<xx<0:cuts.append(m.geo.physical(driver.inverse(np.array([xx])))[0])
        edges=np.unique(np.concatenate(cuts));rJ=(edges[:-1,None]+np.diff(edges)[:,None]*(gx+1)/2).ravel()
        _,delay,_,_,re=m.geo(rJ-m.model.m.RJ);r=re*m.model.m.R;x=C*delay;z=m.coeff(r)
        ids=np.clip(np.searchsorted(d['edges'],rJ,side='right')-1,0,len(d['radius'])-1)
        ql,_=constraints(m,context,driver,it,now,r,x,z,ids)
        jac=1/(np.exp(-2*z['phi']**2)*(1+z['alpha']*r*z['Phi']))
        integrand=z['lapse']*np.sqrt(z['b'])*r*ql*jac
        whole=np.diff(edges)/2*(integrand.reshape(-1,order)@gw);prefix=np.r_[LD(0),np.cumsum(whole,dtype=LD)]
        rows.append(np.asarray(prefix[np.searchsorted(edges,targets)]*np.sqrt(zt['b'])/zt['lapse'],float))
        rc=m.tr[:-1];xc=m.tx[:-1];zc=m.coeff(rc);ic=np.arange(len(rc))
        qn.append(constraints(m,context,driver,it,now,rc,xc,zc,ic)[1])
        if it in [0,len(d['t'])//2,len(d['t'])-1]:
            old=constraints(m,context,driver,it,now,m.r,m.x,m.z,m.ids)[0]
            checks.append(old)
    return np.asarray(rows),np.asarray(qn),np.asarray(checks)


def route(n,order):
    start=time.monotonic();d=dict(np.load(mixed.previous.FIELDS/f'source-{n}.npz'));bg=np.load(mixed.previous.FIELDS/f'fields-{n}-g{order}.npz')
    op=np.load(mixed.previous.OUT/f'operator-{n}-g{order}.npz');old=dict(np.load(mixed.OUT/f'field-{n}-g{order}.npz'))
    m=run.wave.Response();m.setup(d,order);B,I,q,jac,_,_=mixed.fields_and_stress(m,d,bg,op,order)
    z=m.z;r=m.r;A=z['alpha'];beta=-4.;T=z['Eg']-3*z['Pg'];H3=z['Eg']+z['Pg']-3*z['Kg']
    hb=2*B['nu']+4*A*B['f'];hi=2*r*z['Phi']+4*A
    coef=(A*(hb*hi+4*beta*B['f'])+beta*(hb+hi*B['f']))*T
    coef+=(A*hb+beta*B['f'])*(-H3*(3*A+r*z['Phi']))+(A*hi+beta)*B['dt']
    coef=-4*np.pi*z['lapse']**2*z['A4']*coef
    driver=mixed.inc.Driver(order);sample=coef*mixed.inc.ETA*driver.r0*mixed.inc.pulse((d['t'][:,None]+(m.x-m.xfaces[-1])/C)/driver.D)
    old_primary=m.propagate(m.dx*sample)
    mark=time.monotonic();resolved=pulse_propagate(m,coef,mixed.inc.ETA*driver.r0,driver.D);pulse_seconds=time.monotonic()-mark
    jj,qn,identity=resolved_constraints(m,d,bg,order);nn=len(m.r);newJ=jj[:,:nn]
    ks=[0,len(d['t'])//2,len(d['t'])-1]
    error=float(np.max(abs(identity-old['q_lambda'][ks]))/max(np.max(abs(old['q_lambda'][ks])),1e-290));assert error<1e-10,error
    mass_source=m.dx*m.z['K']*(newJ-old['J_mixed']);mass=m.propagate(mass_source)
    delta=[a-b+c for a,b,c in zip(resolved,old_primary,mass)]
    potential=-m.dx*m.z['V']*delta[0][:,:-1][:,m.ids];extra=m.propagate(potential)
    fields=[old[k]+dv+v for k,dv,v in zip(['U','U_t','U_x'],delta,extra)]
    U,Ut,Ux=fields;rr=m.tr;zz=m.tz
    old.update(U=U,U_t=Ut,U_x=Ux,delta_phi=U/rr,delta_Phi=Ux/(rr*zz['lapse']*np.sqrt(zz['b']))-U/rr**2,
        J_mixed=newJ,J_centers=jj[:,nn:],q_nu_centers=qn,
        correction_U=delta[0]+extra[0],exact_primary_coefficient=coef,corrected_mass_source=mass_source,
        corrected_potential_source=potential)
    np.savez_compressed(DEST/f'field-{n}-g{order}.npz',**old)
    return dict(source=n,order=order,seconds=time.monotonic()-start,pulse_seconds=pulse_seconds,
        original_constraint_node_reproduction=error,maximum_corrected_phi=float(np.max(abs(U/rr))),
        old_endpoint_phi=float(np.load(mixed.OUT/f'field-{n}-g{order}.npz')['delta_phi'][-1,-1]),
        new_endpoint_phi=float(U[-1,-1]/rr[-1]))


def main():
    assert not OUT.exists();OUT.mkdir();DEST.mkdir();mixed.previous.initialize();start=time.monotonic()
    run.write(OUT/'check.json',check());rows=[route(*mixed.SETTINGS[0])]
    forecast=4*rows[0]['seconds']+15;run.write(OUT/'pilot.json',dict(first_seconds=rows[0]['seconds'],remaining_upper_seconds=forecast,eligible=forecast<480-(time.monotonic()-start)))
    assert forecast<480-(time.monotonic()-start)
    rows += [route(*setting) for setting in mixed.SETTINGS[1:]]
    run.write(OUT/'result.json',dict(classification='Counterexample candidate',rows=rows,
        source_clock_and_physical_cells_unchanged=True,exact_primary_time_dependence_applied=True,
        constraint_integrals_split_at_existing_centers_and_pulse_fronts=True,full_goal_complete=False))
