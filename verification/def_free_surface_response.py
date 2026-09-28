"""Finite-temperature adiabatic fluid + scalar/metric + outgoing exterior.

This is the mechanical block of the full thermal problem, not an orbit-long
stationary-response assertion. The cold gas completion is explicit and tested
against the saved native overlap; no failed native value is extrapolated.
"""
from pathlib import Path
import argparse
import json
import time
import numpy as np
import sympy as sp
from scipy.optimize import brentq
from scipy.interpolate import PchipInterpolator
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
import def_material_surface_enclosure as surface
import def_radiative_exterior as exterior

h=surface.h
OUT=h.molecular.g.OUT/'def-free-surface-response'


def symbolic():
    r,m,p,en,phi,v,beta=sp.symbols('r m p en phi v beta',real=True)
    b=1-2*m/r;A4=sp.exp(2*beta*phi**2);alpha=beta*phi
    mr=4*sp.pi*r*r*A4*en+r*r*b*v*v/2
    nr=m/(r*r*b)+4*sp.pi*r*A4*p/b+r*v*v/2
    g=nr+alpha*v
    F=4*sp.pi*A4/b*(alpha*(en-3*p)+r*v*(en-p))-2*(r-m)/(r*r*b)*v
    xi,f=sp.symbols('xi f');dm=-4*sp.pi*r*r*A4*(en+p)*xi+r*r*b*v*(f-xi*v)
    Dm=r*r*b*v*f-(4*sp.pi*r*r*A4*p+r*r*b*v*v/2)*xi
    assert sp.simplify(dm+xi*mr-Dm)==0
    # The Lagrangian equations cancel advected equilibrium gradients exactly.
    ep,pp,gp,xip,de,dp,dg,w,omega,N=sp.symbols('ep pp gp xip de dp dg w omega N')
    euler=w*omega**2*xi/(N*N*b)-g*(de+dp-xi*(ep+pp))-w*(dg-xi*gp)
    lag=euler+xip*(-w*g)+xi*(-(ep+pp)*g-w*gp)
    expected=w*omega**2*xi/(N*N*b)-g*(de+dp+w*xip)-w*dg
    assert sp.simplify(lag-expected)==0
    variables=[r,m,p,en,phi,v]
    derivatives=[g,F,*[sp.diff(g,z) for z in variables],*[sp.diff(F,z) for z in variables]]
    fn=sp.lambdify([*variables,beta],derivatives,'numpy',cse=True)
    return fn,dict(classification='Proven',passed=True,
        equations=['Delta m=r^2*b*Phi*Delta phi-(4*pi*r^2*A^4*p+r^2*b*Phi^2/2)*xi',
            'xi_prime=-Delta p/(Gamma1*p)-2*xi/r-Delta lambda-3*alpha*Delta phi',
            'Delta p_prime=w*(omega^2*xi/(N^2*b)+g*(2*xi/r+Delta lambda+3*alpha*Delta phi)-Delta g)-g*Delta p',
            'Delta phi_prime=Delta Phi+xi_prime*Phi',
            'Delta Phi_prime=sum(F_z*Delta z)+xi*F_r+xi_prime*F-omega^2*(Delta phi-xi*Phi)/(N^2*b)'],
        scope='Adiabatic Delta s=Delta X=0 on a static mechanical background, finite-temperature Gamma1. Not the omitted thermal/reactive perturbation block.')


def build():
    assert not OUT.exists();OUT.mkdir();start=time.monotonic()
    files=[Path(__file__),surface.OUT/'result.json',h.OUT/'absolute-shoot/background-0.001.npz',
        h.OUT/'absolute-shoot/native-audit.npz',h.OUT/'absolute-shoot/lapse.npz',
        surface.atmosphere.parent.OUT/'extended-cold/progress.json',h.molecular.model.OUT/'compiled-level-inputs.npz',
        exterior.s.OUT/'companion-benchmark.json']
    _,proof=symbolic()
    h.write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='Phase49 checkpoint',
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in files},symbolic=proof,
        claim='Actually solve the finite-temperature adiabatic material/scalar/metric response with regular centre, free material surface and dynamic outgoing exterior. This supplies a necessary mechanical block, not a replacement for thermal/nonstationary response.',
        cold_completion='Neutral ideal mixture + the same 302-level H2 sum + photons, native composition and molecular ground energy. Entropy additive constant set by the last native state; check four additional native overlap states. Excluded nonideal/excited terms and the gas-branch premise remain physical uncertainties.',
        gates=dict(cold_pressure_relative=1e-8,cold_thermal_energy_relative=1e-8,cold_Gamma1_relative=1e-7,linear_residual=1e-9,spatial_response_relative=.02),
        frequencies_harmonics=[0,1,2,3],refinements=[1,2],
        budget=dict(build_hard_timeout_seconds=30,pilot_hard_timeout_seconds=30,remaining_hard_timeout_seconds=90,native_calls=0,new_time_steps=0,automatic_expansion=False),
        decision='Use the coupled result to identify the response scale and whether mechanical amplification is present. Do not turn adiabatic frozen-background output into measured charge, full radiation hydrodynamics or a stationary orbital result.'))
    data,_=h.inputs();body=h.Structure(.001);bg=np.load(h.OUT/'absolute-shoot/background-0.001.npz');lap=np.load(h.OUT/'absolute-shoot/lapse.npz');native=np.load(h.OUT/'absolute-shoot/native-audit.npz')['raw']
    states=bg['states'];R=bg['faces'][0,0]*body.R;rf=bg['faces'][::-1,0]*body.R/R
    # All stored cells and native finite-temperature Gamma1 are retained.
    rc=states[::-1,0]*body.R/R;pc=np.exp(states[::-1,2]);raw=native[::-1];rho=raw[:,0]
    cx=data['CX'][::-1];e=rho*((h.gr.C*100)**2*cx+raw[:,2])
    phic=.001*(1+body.mu*states[::-1,3]);vc=.001*body.mu*states[::-1,4]*R/body.R
    geo_p=h.gr.G*.1*R*R/h.gr.C**4
    centre=body.state(bg['parameters'][0],len(data['dm'])-1)
    fields=dict(r=np.r_[0,rc],m=np.r_[0,states[::-1,1]*body.B/R],p=np.r_[centre[0]*R*R,pc*geo_p],
        e=np.r_[centre[1]*R*R,e*geo_p],phi=np.r_[.001*(1+body.mu*bg['parameters'][3]),phic],
        v=np.r_[0,vc],N=np.exp(np.r_[lap['nu_faces'][-1],lap['nu_mid'][::-1]]),gamma=np.r_[raw[0,4],raw[:,4]])
    rows=json.loads((surface.atmosphere.parent.OUT/'extended-cold/progress.json').read_text())['rows'];a=np.array([z['raw'] for z in rows]);T=np.exp([z['lnT'] for z in rows]);cx0=float(data['CX'][0]);X=data['X'][0]
    # Exact constants are those in the archived native mod_free_eos_constants.
    k=1.380649e-16;NA=6.02214076e23;planck=6.62607015e-27;c=2.99792458e10
    arad=8*np.pi**5*k**4/(15*planck**3*c**3);gas=k*NA
    nuclei=X/h.molecular.g.c.A;hydrogen=nuclei[h.molecular.g.c.Z==1].sum();Rg=gas*(nuclei.sum()-hydrogen/2);Rh=gas*hydrogen/2
    levels=np.load(h.molecular.model.OUT/'compiled-level-inputs.npz');energy=levels['kelvin'];weights=levels['weights'];u0=-cx0*a[-1,11]
    def neutral(rho,T):
        xx=energy/T;ww=weights*np.exp(-xx);prob=ww/ww.sum();mean=prob@xx;variance=prob@((xx-mean)**2)
        pr=arad*T**4/3;p=rho*Rg*T+pr;u=1.5*Rg*T+Rh*T*mean+3*pr/rho
        entropy=Rg*(1.5*np.log(T)-np.log(rho))+Rh*(np.log(ww.sum())+mean)+4*pr/(rho*T)
        cvT=1.5*Rg*T+Rh*T*variance+12*pr/rho;slope=(p/rho+3*pr/rho)/cvT
        gamma=(rho*Rg*T+(rho*Rg*T+4*pr)*slope)/p
        return p,u,entropy,gamma
    sref=neutral(a[-1,0],T[-1])[2];offset=a[-1,3]-sref;controls=[]
    for j in range(-5,0):
        p,u,s,gamma=neutral(a[j,0],T[j]);actual=a[j,2]-u0
        controls.append(dict(index=j,temperature_K=float(T[j]),pressure_relative=abs(p/a[j,1]-1),thermal_energy_relative=abs(u/actual-1),Gamma1_relative=abs(gamma/a[j,4]-1),entropy_difference=float(s+offset-a[j,3])))
    accepted=all(v['pressure_relative']<1e-8 and v['thermal_energy_relative']<1e-8 and v['Gamma1_relative']<1e-7 for v in controls)
    h.write(OUT/'cold-overlap.json',dict(classification='Counterexample candidate',passed=accepted,rows=controls,entropy_offset=float(offset),scope='Finite native overlap, not a uniform physical error certificate.'))
    assert accepted,controls
    cold=[];lr0=np.log(a[-1,0])
    for t in np.geomspace(T[-1],1e-5,65)[1:]:
        f=lambda lr:neutral(np.exp(lr),t)[2]-sref
        lr=brentq(f,lr0-100,lr0,xtol=2e-13);rh=np.exp(lr);p,u,_,ga=neutral(rh,t);lr0=lr
        cold.append((rh,p,u+u0,ga,t))
    material=[(z[0],z[1],z[2],z[4],tt) for z,tt in zip(a,T)]+cold+[(0,0,u0,5/3,0)]
    mp=surface.mp;mp.mp.dps=65;M=mp.mpf(float(bg['faces'][0,1]*body.B));RR=mp.mpf(float(R));mu=M/RR
    pb=mp.mpf(float(.001*(1+body.mu*bg['faces'][0,3])));qr=mp.mpf(float(.001*body.mu*bg['faces'][0,4]*R/body.R));Nb=mp.exp(mp.mpf(float(lap['nu_faces'][0])));Q=Nb*mp.sqrt(1-2*mu)*qr
    C2=mp.mpf(float(h.gr.C*100))**2;hb=mp.mpf(cx0)+(mp.mpf(float(a[0,2]))+mp.mpf(float(a[0,1]))/mp.mpf(float(a[0,0])))/C2
    outer=[]
    for rh,p,u,ga,t in material:
        H=mp.mpf(cx0)+(mp.mpf(float(u))+(mp.mpf(float(p))/mp.mpf(float(rh)) if rh else 0))/C2
        r,z,ph,m=surface.just_surface(mp,mu,qr,pb,mp.log(hb/H),RR);x=float(r/RR);N=float(Nb*mp.exp(z));v=float(Q)/(N*np.sqrt(1-2*float(m/RR)/x)*x*x)
        outer.append([x,float(m/RR),p*geo_p,rh*((h.gr.C*100)**2*cx0+u)*geo_p,float(ph),v,N,ga])
    outer=np.array(outer)
    # A radius below machine separation does not define a new usable cell.
    keep=np.r_[True,np.diff(outer[:,0])>2e-13];keep[-1]=False;outer=np.r_[outer[keep],outer[-1:]]
    for j,key in enumerate(fields):fields[key]=np.r_[fields[key],outer[:,j]]
    assert np.all(np.diff(fields['r'])>0)
    grid=np.r_[rf,outer[1:,0]]
    np.savez_compressed(OUT/'background.npz',**fields,grid=grid,R=R,source_dm=bg['dm'],source_X=bg['X'],outer_start=len(rc)+1)
    h.write(OUT/'background.json',dict(classification='Counterexample candidate',original_cells=len(rc),nodes=len(grid),cold_completion_overlap_passed=True,
        radius_m=float(R*grid[-1]),seconds=time.monotonic()-start,native_calls=0,static_atmosphere_reference='Exact Just test-atmosphere from Phase48; material-source geometric error bounded there, not a dynamic response bound.',thermal_stationarity=False))
    print('BACKGROUND',len(grid),'seconds',time.monotonic()-start,flush=True)


def sample(grid):
    saved=np.load(OUT/'background.npz');r=saved['r'];x=(grid[:-1]+grid[1:])/2
    out={k:np.interp(x,r,saved[k]) for k in ['m','p','e','phi','v','N','gamma']};out['r']=x
    # Regular centre powers, not a spurious central point mass or scalar cusp.
    centre=x<r[1]
    out['m'][centre]=saved['m'][1]*(x[centre]/r[1])**3
    out['v'][centre]=saved['v'][1]*x[centre]/r[1]
    return out


def operators(bg,omega):
    fn,_=symbolic();x,m,p,en,phi,v,N,ga=[bg[k] for k in ['r','m','p','e','phi','v','N','gamma']]
    beta=-4.;b=1-2*m/x;A4=np.exp(2*beta*phi*phi);alpha=beta*phi;w=en+p
    deriv=fn(x,m,p,en,phi,v,beta);g,F=deriv[:2];dg,dF=deriv[2:8],deriv[8:]
    n=np.size(x);basis=np.eye(4)[:,:,None]+np.zeros((4,4,n));z,eta,f,V=basis;xi=x*z
    Dm=x*x*b*v*f-(4*np.pi*x*x*A4*p+x*x*b*v*v/2)*xi
    Dl=(Dm/x-m*xi/(x*x))/b;drho=eta/ga;xip=-drho-2*z-Dl-3*alpha*f
    increments=[xi,Dm,p*eta,w*drho,f,V]
    Dg=sum(j*q for j,q in zip(dg,increments));DF=sum(j*q for j,q in zip(dF,increments))
    bracket=omega*omega/(N*N*b)*xi+g*(2*z+Dl+3*alpha*f)-Dg+g*eta
    scalar=DF+xip*F-omega*omega/(N*N*b)*(f-xi*v)
    with np.errstate(divide='ignore',invalid='ignore'):
        result=np.array([(xip-z)/x,w/p*bracket-g*eta,V+xip*v,scalar])
    return np.moveaxis(result,-1,0),bracket,g,F


def solve(harmonic,refinement):
    saved=np.load(OUT/'background.npz');old=saved['grid'];grid=np.sort(np.concatenate([old]+[old[:-1]+j*np.diff(old)/refinement for j in range(1,refinement)]))
    R=float(saved['R']);benchmark=json.loads((exterior.s.OUT/'companion-benchmark.json').read_text());omega=harmonic*2*np.pi*R/(h.gr.C*benchmark['leading_drive']['period_seconds']);began=time.monotonic()
    A,_,_,_=operators(sample(grid),omega);n=len(grid)-1;step=np.diff(grid);I=np.eye(4)
    left=-I[None,:,:]-step[:,None,None]*A/2;right=I[None,:,:]-step[:,None,None]*A/2
    rr=[];cc=[];vv=[]
    for shift,block in [(0,left),(4,right)]:
        rows=np.arange(4*n).reshape(n,4,1)+np.zeros((1,1,4),int);cols=4*np.arange(n)[:,None,None]+np.arange(4)[None,None,:]+shift+np.zeros((1,4,1),int)
        rr.extend(rows.ravel());cc.extend(cols.ravel());vv.extend(block.ravel())
    rhs=np.zeros(4*(n+1),complex);endbg={k:np.array([saved[k][-1]]) for k in ['r','m','p','e','phi','v','N','gamma']}
    _,bracket,g,F=operators(endbg,omega);Rs=grid[-1];m=endbg['m'][0];Phi=endbg['v'][0];phi=endbg['phi'][0];flux=Rs*Phi
    wave=exterior.outgoing(m/Rs,flux,omega*Rs);Z=wave['impedance'];hw=wave['h'];Fmetric=endbg['N'][0]*np.sqrt(1-2*m/Rs)
    drive=np.exp(-1j*omega*Rs)/(Fmetric*hw)
    centre=np.array([3*saved['gamma'][0],1,3*saved['gamma'][0]*(-4*saved['phi'][0]),0])
    boundary=[(0,centre,0),(0,[0,0,0,1],0),(4*n,bracket[:,0]/g[0],0),
        (4*n,[-Rs*Rs*F[0]+Z*Rs*Phi,0,-Z,Rs],drive)]
    for j,(offset,row,value) in enumerate(boundary):
        rr.extend([4*n+j]*4);cc.extend(offset+np.arange(4));vv.extend(row);rhs[4*n+j]=value
    matrix=coo_matrix((np.asarray(vv,complex),(rr,cc)),shape=(len(rhs),len(rhs))).tocsc()
    scale=np.asarray(abs(matrix).sum(1)).ravel();scaled=matrix.multiply((1/scale)[:,None]).tocsc();answer=spsolve(scaled,rhs/scale)
    residual=float(np.max(abs(matrix@answer-rhs)/(np.asarray(abs(matrix)@abs(answer)).ravel()+abs(rhs)+1e-100)))
    y=answer.reshape(n+1,4);euler_surface=y[-1,2]-Rs*y[-1,0]*Phi
    name=f'harmonic-{harmonic}-grid-{refinement}';np.savez_compressed(OUT/(name+'.npz'),grid=grid,response=y)
    value=dict(classification='Counterexample candidate',harmonic=harmonic,refinement=refinement,cells=n,omega_R0_over_c=omega,
        seconds=time.monotonic()-began,linear_residual=residual,surface_displacement_over_radius=[float(y[-1,0].real),float(y[-1,0].imag)],
        Eulerian_surface_scalar=[float(euler_surface.real),float(euler_surface.imag)],surface_Lagrangian_scalar=[float(y[-1,2].real),float(y[-1,2].imag)],
        h=[float(hw.real),float(hw.imag)],impedance=[float(Z.real),float(Z.imag)],
        radiative_current_relative_error=wave['current_relative_error'],normalization='Unit regular incident monopole; exp(-i omega t), incident phase at the surface is exp(-i omega R/c). At zero frequency use its regular Wronskian limit.',
        finite_temperature_adiabatic_only=True,thermal_stationarity=False,full_dynamic_charge_solved=False)
    h.write(OUT/(name+'.json'),value);print(name,'seconds',value['seconds'],'residual',residual,'xi',value['surface_displacement_over_radius'],flush=True)
    assert residual<1e-9
    return value


def pilot():
    assert not (OUT/'pilot.json').exists();row=solve(1,1)
    h.write(OUT/'pilot.json',dict(classification='Counterexample candidate',result=row,forecast_remaining_seconds=row['seconds']*14,scope='Finest grid and other harmonics unmeasured; no automatic expansion.'))


def run():
    assert not (OUT/'result.json').exists();pilot=json.loads((OUT/'pilot.json').read_text());assert pilot['forecast_remaining_seconds']<85,pilot
    results=[]
    for harmonic in [0,1,2,3]:
        pair=[]
        for refinement in [1,2]:
            record=OUT/f'harmonic-{harmonic}-grid-{refinement}.json'
            pair.append(json.loads(record.read_text()) if record.exists() else solve(harmonic,refinement))
        a,b=[np.load(OUT/f'harmonic-{harmonic}-grid-{refinement}.npz') for refinement in [1,2]]
        ratio=float(np.max(abs(b['response'][::2,0]-a['response'][:,0]))/np.max(abs(b['response'][:,0])))
        results.append(dict(harmonic=harmonic,solves=pair,relative_displacement_grid_difference=ratio,spatial_gate_passed=ratio<.02))
    h.write(OUT/'result.json',dict(classification='Counterexample candidate',coupled_adiabatic_solutions_completed=True,rows=results,spatial_gate_passed=all(q['spatial_gate_passed'] for q in results),native_calls=0,new_time_steps=0,physical_cold_EOS_certified=False,thermal_response_included=False,thermal_stationarity=False,full_dynamic_charge_solved=False))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['build','pilot','run']);globals()[p.parse_args().action]()
