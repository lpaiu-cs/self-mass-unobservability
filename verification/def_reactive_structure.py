"""Initial reactive quasi-static free-surface mass and scalar readout.

This is a tangent of equilibria, not a time-dependent transport evolution.
Delta m remains independent: the adiabatic momentum constraint is insufficient
when entropy and composition change.
"""
from pathlib import Path
import json
import time
import numpy as np
import sympy as sp
import mpmath as mp
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
import def_reactive_energy as chemical
import def_free_surface_response_normalized as surface

h=surface.h
OUT=chemical.OUT/'structure'


def symbolic():
    r,m,p,e,phi,v,beta=sp.symbols('r m p e phi v beta',real=True)
    b=1-2*m/r;A4=sp.exp(2*beta*phi**2);alpha=beta*phi
    M=4*sp.pi*r*r*A4*e+r*r*b*v*v/2
    g=m/(r*r*b)+4*sp.pi*r*A4*p/b+r*v*v/2+alpha*v
    F=4*sp.pi*A4/b*(alpha*(e-3*p)+r*v*(e-p))-2*(r-m)/(r*r*b)*v
    z,eta,f,V,d,rr,ee,ga=sp.symbols('z eta f V d rr ee ga')
    xi=r*z;dl=(d-m*z/r)/b;xip=-eta/ga-rr-2*z-dl-3*alpha*f
    inc=[xi,r*d,p*eta,(e+p)*eta/ga+ee,f,V]
    DG=lambda G:sum(sp.diff(G,t)*u for t,u in zip([r,m,p,e,phi,v],inc))
    dp_prime=-g*((e+p)*eta/ga+ee+p*eta+(e+p)*xip)-(e+p)*DG(g)
    bracket=g*(2*z+dl+3*alpha*f+rr-ee/(e+p)+eta)-DG(g)
    assert sp.simplify(dp_prime/p+(e+p)*g*eta/p-((e+p)*bracket/p-g*eta))==0
    fn=sp.lambdify([r,m,p,e,phi,v,beta],[M,g,F,*[sp.diff(G,t) for G in [M,g,F] for t in [r,m,p,e,phi,v]]],'numpy',cse=True)
    return fn,dict(classification='Proven',passed=True,
        equations=['Delta ln rho=eta/Gamma1+rho_ref; Delta e=w*eta/Gamma1+e_ref',
            'xi_prime=-eta/Gamma1-rho_ref-2*zeta-Delta lambda-3*alpha*f',
            'Delta p_prime=-g*(Delta e+Delta p+w*xi_prime)-w*Delta g',
            'f_prime=V+xi_prime*Phi; V_prime=Delta F+xi_prime*F',
            'Delta m_prime=Delta M+xi_prime*M'],
        source_identity='At fixed pressure, d(e/rho+p/rho)=-loss: e_ref=w*rho_ref-rho*loss. Nuclear rest-to-thermal conversion must not be counted as added total energy.',
        scope='Lagrangian variation of static field equations at fixed baryon inventory. No assumption of a stationary thermal background or a dynamic transport solution.')


def coefficients(bg,rho_ref,loss,fn):
    r,m,p,e,phi,v,N,ga=[bg[k] for k in ['r','m','p','e','phi','v','N','gamma']]
    b=1-2*m/r;alpha=-4*phi;w=e+p
    val=fn(r,m,p,e,phi,v,-4.);M,g,F=val[:3];dM,dg,dF=val[3:9],val[9:15],val[15:21]
    # Last column is the known source, not another unknown.
    basis=np.zeros((5,6,len(r)));basis[:,:5,:]=np.eye(5)[:,:,None]
    z,eta,f,V,d=basis;xi=r*z
    rr=np.zeros_like(z);rr[-1]=rho_ref
    le=np.zeros_like(z);le[-1]=loss
    dl=(d-m*z/r)/b;xip=-eta/ga-rr-2*z-dl-3*alpha*f
    de=w*(eta/ga+rr)-le
    inc=[xi,r*d,p*eta,de,f,V]
    DG=lambda der:sum(q*u for q,u in zip(der,inc))
    # Cancel w*rho_ref analytically against e_ref: retain the actual neutrino
    # energy loss even when subtracting two large thermochemical terms loses it.
    with np.errstate(divide='ignore',invalid='ignore'):
        bracket=g*(2*z+dl+3*alpha*f+eta)-DG(dg)
        bracket[-1]+=np.divide(g*loss,w,out=np.zeros_like(w),where=w!=0)
        result=np.array([(xip-z)/r,w/p*bracket-g*eta,V+xip*v,DG(dF)+xip*F,(DG(dM)+xip*M-d)/r])
    return np.moveaxis(result,-1,0),bracket,M,F


def exterior_rows(saved):
    mp.mp.dps=60
    R=float(saved['grid'][-1]);mu=mp.mpf(float(saved['m'][-1]/R));q=mp.mpf(float(R*saved['v'][-1]))
    exact=surface.old.exterior.q.exterior.exact
    funcs=[lambda a,c:a+c*c*exact(a,c)[0],lambda a,c:-c*exact(a,c)[1],lambda a,c:c*exact(a,c)[2]]
    values=np.array([float(f(mu,q)) for f in funcs])
    jac=np.array([[float(mp.diff(lambda a:f(a,q),mu)),float(mp.diff(lambda c:f(mu,c),q))] for f in funcs])
    return values,jac


def solve(refinement,probe,fn,ext,drive=0.,outer='zero'):
    saved=np.load(surface.OUT/'background.npz');old=saved['grid']
    grid=np.sort(np.concatenate([old]+[old[:-1]+j*np.diff(old)/refinement for j in range(1,refinement)]))
    bg=surface.namespace['sample'](grid);r=bg['r'];R0=float(saved['R'])
    d=np.load(chemical.paired.old.thermal.OUT/'coefficients.npz');s=np.load(chemical.paired.old.thermal.OUT/'sources.npz')
    a=np.load(chemical.OUT/'forcing.npz')['rows'];rx=d['radius_cm'][::-1]/(100*R0)
    rates=a[::-1,probe,0];geo=h.gr.G*.1*R0*R0/h.gr.C**4
    losses=(d['raw'][:,0]*d['A']*d['N']*(s['neutrino']+s['thermal_neutrino']))[::-1]*geo
    if drive:rates=rates*0;losses=losses*0
    rho_ref=np.interp(r,rx,rates,right=0 if outer=='zero' else rates[-1])
    loss=np.interp(r,rx,losses,right=0 if outer=='zero' else losses[-1])
    # Any continuation vanishes at the true surface; this avoids a spurious
    # finite volumetric source in vacuum. The alternative is a sensitivity test.
    if outer=='continued':
        mask=r>rx[-1];shape=bg['e'][mask]/np.interp(rx[-1],saved['r'],saved['e'])
        loss[mask]=losses[-1]*shape
    A,_,_,_=coefficients(bg,rho_ref,loss,fn);n=len(r);dx=np.diff(grid);I=np.eye(5)
    left=-I[None]-dx[:,None,None]*A[:,:,:5]/2;right=I[None]-dx[:,None,None]*A[:,:,:5]/2
    rows=[];cols=[];values=[]
    for shift,block in [(0,left),(5,right)]:
        rows.extend((5*np.arange(n)[:,None,None]+np.arange(5)[None,:,None]+np.zeros((1,1,5),int)).ravel())
        cols.extend((5*np.arange(n)[:,None,None]+np.arange(5)[None,None,:]+shift+np.zeros((1,5,1),int)).ravel())
        values.extend(block.ravel())
    rhs=np.zeros(5*(n+1));rhs[:5*n]=(dx[:,None]*A[:,:,5]).ravel()
    end={k:np.array([saved[k][-1]]) for k in ['r','m','p','e','phi','v','N','gamma']}
    _,bracket,M,F=coefficients(end,np.zeros(1),np.zeros(1),fn)
    R=grid[-1];v=saved['v'][-1];alpha=-4*saved['phi'][0]
    ev,jac=ext
    # Eulerian endpoint increments: delta mu=d-zeta*M, delta q=R*V-R^2*zeta*F.
    phiinfty=np.array([-R*v-jac[2,0]*M[0]-jac[2,1]*R*R*F[0],0,1,jac[2,1]*R,jac[2,0]])
    bcs=[(0,[3,1/saved['gamma'][0],3*alpha,0,0],-rates[0]),
         (0,[0,0,0,1,0],0),(0,[0,0,0,0,1],0),
         (5*n,bracket[:5,0]/float(fn(R,saved['m'][-1],0,0,saved['phi'][-1],v,-4.)[1]),0),
         (5*n,phiinfty,drive)]
    for j,(offset,row,target) in enumerate(bcs):
        rows.extend([5*n+j]*5);cols.extend(offset+np.arange(5));values.extend(row);rhs[5*n+j]=target
    matrix=coo_matrix((values,(rows,cols)),shape=(len(rhs),len(rhs))).tocsc()
    scale=np.asarray(abs(matrix).sum(1)).ravel()
    y=spsolve(matrix.multiply((1/scale)[:,None]).tocsc(),rhs/scale).reshape(n+1,5)
    residual=float(np.max(abs(matrix@y.ravel()-rhs)/(np.asarray(abs(matrix)@abs(y.ravel())).ravel()+abs(rhs)+1e-100)))
    z,eta,f,V,dd=y[-1];dmu=dd-z*M[0];dq=R*V-R*R*z*F[0]
    dm,dcharge,dshift=jac@np.array([dmu,dq]);mass=R*ev[0];charge=R*ev[1]
    mass_rate=R*dm;charge_rate=R*dcharge;normalized=charge_rate/mass-charge*mass_rate/mass**2
    row=dict(classification='Counterexample candidate',refinement=refinement,probe=probe,outer=outer,drive=drive,
        linear_residual=residual,surface_radial_fraction_rate=float(z),surface_m_per_second=float(z*R*R0),
        ADM_geom_m_per_second=float(mass_rate*R0),scalar_charge_geom_m_per_second=float(charge_rate*R0),
        normalized_charge_rate=float(normalized),asymptotic_scalar_residual=float(f-R*z*v+dshift-drive))
    name=f"{'adiabatic' if drive else 'reactive'}-{probe}-{refinement}-{outer}"
    np.savez_compressed(OUT/(name+'.npz'),grid=grid,response=y)
    h.write(OUT/(name+'.json'),row)
    return row,y


def main():
    assert not OUT.exists();OUT.mkdir();began=time.monotonic();fn,proof=symbolic()
    paths=[Path(__file__),chemical.OUT/'forcing.npz',chemical.OUT/'result.json',surface.OUT/'background.npz',
           chemical.paired.old.thermal.OUT/'sources.npz',chemical.paired.old.thermal.OUT/'coefficients.npz']
    h.write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='b3552307',
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths},symbolic=proof,
        claim='Connect native reaction plus chemical energy to baryon-preserving quasi-static free-surface displacement, independent mass constraint and exact Just scalar/mass readout at fixed asymptotic scalar.',
        missing='This tangent omits inertia and evolving heat transport. It is neither a physical trajectory nor a measured orbital charge.',
        outer_source='Reference source is zero beyond the original native inventory. Compare constant fractional density source and density-scaled volumetric loss outside it; both are declared continuations, not physical certification.',
        gates=dict(linear_residual=1e-9,adiabatic_control_relative=.001,spatial_displacement=.02,spatial_charge=.02,probe_readout=.01,mass_first_law=.02),
        budget=dict(hard_timeout_seconds=60,native_calls=0,new_time_steps=0,maximum_solves=6,refinements=[1,2],automatic_expansion=False)))
    saved=np.load(surface.OUT/'background.npz');ext=exterior_rows(saved)
    control,y=solve(1,1,fn,ext,drive=1.)
    original=np.load(surface.OUT/'harmonic-0-grid-1.npz')['response'].real
    control_error=float(np.max(abs(y[:,0]-original[:,0]))/np.max(abs(original[:,0])))
    rows=[]
    for probe in [0,1]:
        for refinement in [1,2]:rows.append(solve(refinement,probe,fn,ext)[0])
    outer=solve(2,1,fn,ext,outer='continued')[0]
    fine=rows[-1];coarse=rows[-2];wide=rows[1]
    fields=['surface_radial_fraction_rate','normalized_charge_rate']
    spatial={k:abs(fine[k]-coarse[k])/max(abs(fine[k]),1e-100) for k in fields}
    probes={k:abs(fine[k]-wide[k])/max(abs(fine[k]),1e-100) for k in fields}
    exterior_control={k:abs(fine[k]-outer[k])/max(abs(fine[k]),1e-100) for k in fields}
    d=np.load(chemical.paired.old.thermal.OUT/'coefficients.npz');s=np.load(chemical.paired.old.thermal.OUT/'sources.npz')
    power=float(d['dm']@(d['A']**2*d['N']**2*(s['neutrino']+s['thermal_neutrino'])))
    expected=-power*h.gr.G*.1/h.gr.C**4
    mass_error=abs(fine['ADM_geom_m_per_second']/expected-1)
    result=dict(classification='Counterexample candidate',symbolic=proof,rows=rows,adiabatic_control=control,
        adiabatic_displacement_relative_difference=control_error,spatial_relative_difference=spatial,
        probe_relative_difference=probes,outer_source_control_relative_difference=exterior_control,
        expected_neutrino_ADM_geom_m_per_second=expected,mass_first_law_relative_difference=mass_error,
        all_gates_passed=bool(control_error<.001 and max(spatial.values())<.02 and max(probes.values())<.01 and mass_error<.02 and max(r['linear_residual'] for r in rows+[control,outer])<1e-9),
        seconds=time.monotonic()-began,physical_time_evolution=False,thermal_transport_evolved=False,full_dynamic_charge_solved=False)
    h.write(OUT/'result.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
