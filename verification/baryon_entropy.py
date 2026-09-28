"""Request23: Lagrangian, composition/entropy preserving thermal TOV projection.

Counterexample candidate. Preserve each original cell's baryon fraction and
FreeEOS reference entropy. Hydrostatic projection, not a stellar evolution run.
"""
from pathlib import Path
import json, sys, time
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import least_squares
from scipy.integrate import solve_ivp
from thermal_wd import mesa
from thermal_restart import sha
import gr_mass as gr

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/baryon-entropy23'
OLD=ROOT/'outputs/thermal-closure22'


def save(name,obj):
    OUT.mkdir(exist_ok=True,parents=True)
    (OUT/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')


class Material:
    def __init__(self):
        self.h,self.d=mesa(OLD/'selected.data.gz');d=self.d
        self.path=gr.ThermalPath();self.lp=self.path.lp;self.lt=d['logT']*np.log(10.)
        self.baryon_g=float(sum(d['dm']));self.B=self.baryon_g*.001*gr.G/gr.C**2
        self.R=float(self.h['photosphere_r'])*gr.RSUN
        self.dm=d['dm']/self.baryon_g
        # Sum from the nearer endpoint: avoid subtracting a nearly total mass
        # when resolving the exceedingly light photospheric cells.
        self.outer=np.r_[0.,np.cumsum(self.dm)]
        self.inner=np.r_[np.cumsum(self.dm[::-1])[::-1],0.]
        self.outer[-1]=1.;self.inner[0]=1.
        self.split=int(np.argmin(abs(self.outer-.5)))
        self.cx=np.array([self.path(p)[1] for p in self.lp])
        self.eps=np.array([self.path(p)[2] for p in self.lp])


def prepare():
    assert not OUT.exists();m=Material()
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='729aca0',
        previous_manifest_sha256=sha(OLD/'manifest.json'),profile_sha256=sha(OLD/'selected.data.gz'),
        closure='Original cell baryon masses and renormalized non-F isotope fractions are fixed. Each cell has constant composition and entropy per baryon gram, defined by FreeEOS at its original P,T. Pressure, density, temperature and radius adjust by TOV. Piecewise-constant material data are the declared model, not a smoothed stellar evolution solution.',
        entropy_reference='FreeEOS entropy at the source P,T,X, multiplied by C_X. This is not certification of the original MESA blended-EOS absolute entropy.',
        primary='Fixed original baryon mass; central pressure, optical radius and gravitational mass determined by centre/surface matching.',
        secondary='Only if required, separately register and solve a family with all cell baryon masses scaled together to match gravitational mass. Preserve the primary result; the scaled family is not the same original material.',
        boundary='Original nonzero optical pressure; append the previously declared finite gamma=5/3 mathematical atmosphere and report its additional baryon budget separately. No atmosphere transport is solved.',
        entropy_table=dict(logP_offset=[-1.,1.],initial_points=17,refinement_points=33),
        checks=['EOS entropy inversion and independent finite-difference heat capacity',
            'Cell-preserving inventories and entropy','RK4 step subdivision 1,2,4',
            'Direct EOS evaluation independent of the lookup tables','Symbolic baryon-coordinate TOV identities'],
        tolerances=dict(interface=1e-8,GR_mass_relative=1e-8,entropy_relative=1e-10),
        not_solved=['reactions','thermal_transport','rotation','stability','fluid_scalar_response','observational_inference']))
    print('Material cells',len(m.lp),'baryon_g',m.baryon_g,'split',m.split,flush=True)


def invert(eos,lp,entropy,eps,guess):
    lt=float(guess)
    for iteration in range(15):
        a=eos(1,lp,lt,eps)
        # dh/dlnT|P / T = cp = ds/dlnT|P, at fixed composition.
        cp=(a[10]-a[1]/a[0]*a[8])/np.exp(lt)
        assert cp>0 and np.isfinite(cp)
        error=a[3]-entropy
        if abs(error)<2e-12*max(abs(entropy),1.): return a,lt,iteration+1
        lt-=np.clip(error/cp,-.15,.15)
    raise ValueError(('entropy inversion',lp,lt,error,entropy))


def tables(points):
    m=Material();eos=gr.EOS();offset=np.linspace(-1.,1.,points)
    reference=np.array([eos(1,p,t,eps) for p,t,eps in zip(m.lp,m.lt,m.eps)])
    entropy=reference[:,3];data=np.empty((len(m.lp),points,4));iterations=[]
    start=time.monotonic()
    for i in range(len(m.lp)):
        a=reference[i];cp=(a[10]-a[1]/a[0]*a[8])/np.exp(m.lt[i])
        nabla=-a[1]/a[0]/np.exp(m.lt[i])*a[8]/cp
        for j,dx in enumerate(offset):
            a,lt,n=invert(eos,m.lp[i]+dx,entropy[i],m.eps[i],m.lt[i]+nabla*dx)
            data[i,j]=[np.log(a[0]),lt,a[2]/gr.C**2*1e-4,a[3]];iterations.append(n)
        if i%500==0: print('EOS cell',i,'/',len(m.lp),'seconds',time.monotonic()-start,flush=True)
    np.savez_compressed(OUT/f'adiabats-{points}.npz',offset=offset,values=data,reference=reference)
    save(f'table-{points}.json',dict(classification='Counterexample candidate',cells=len(m.lp),points=points,
        max_Newton_iterations=int(max(iterations)),elapsed_s=time.monotonic()-start,
        max_relative_entropy_residual=float(np.max(abs(data[:,:,3]/entropy[:,None]-1)))))
    print('EOS table complete',points,time.monotonic()-start,flush=True)


def table17(): tables(17)
def table33(): tables(33)


class Solver:
    def __init__(self,points=17,subdivision=1,direct=False,method='RK4'):
        self.mat=Material();m=self.mat;self.sub=subdivision;self.direct=direct;self.method=method
        data=np.load(OUT/f'adiabats-{points}.npz');self.ref=data['reference'];self.offset=data['offset']
        self.coef=PchipInterpolator(self.offset,data['values'][:,:,:3],axis=1).c
        self.eos=gr.EOS();self.calls=0;self.max_entropy_error=0.

    def state(self,lp,i):
        self.calls+=1;m=self.mat
        dx=lp-m.lp[i]
        if not -1<=dx<=1: raise ValueError(('pressure outside declared EOS grid',i,dx))
        j=min(len(self.offset)-2,max(0,int(np.searchsorted(self.offset,dx)-1)));z=dx-self.offset[j]
        a=self.coef[:,j,i,:];v=((a[0]*z+a[1])*z+a[2])*z+a[3]
        if self.direct:
            real,lt,_=invert(self.eos,lp,self.ref[i,3],m.eps[i],v[1])
            self.max_entropy_error=max(self.max_entropy_error,abs(real[3]/self.ref[i,3]-1))
            rho=real[0];u=real[2]*1e-4/gr.C**2
        else: rho=np.exp(v[0]);u=v[2]
        p=gr.G*np.exp(lp)*.1/gr.C**4
        b=gr.G*(rho/m.cx[i])*1000/gr.C**2;e=gr.G*rho*1000/gr.C**2*(1+u)
        return p,e,b

    def rhs(self,x,y,i,B,outer):
        m=self.mat;r=y[0]*m.R;mass=y[1]*B;lp=y[2]
        p,e,b=self.state(lp,i);f=1-2*mass/r
        if r<=0 or mass<=0 or f<=0: raise ValueError(('invalid structure',r,mass,f))
        db=B*np.exp(x)*(-1 if outer else 1)
        rB=np.sqrt(f)/(4*np.pi*r*r*b)
        massB=(e/b)*np.sqrt(f)
        lpB=-(e+p)*(mass+4*np.pi*r**3*p)/(4*np.pi*r**4*b*np.sqrt(f)*p)
        return np.array([rB/m.R,massB/B,lpB])*db

    def step(self,a,b,y,i,B,outer):
        if self.method=='DOP853':
            result=solve_ivp(lambda x,y:self.rhs(x,y,i,B,outer),(a,b),y,
                method='DOP853',rtol=2e-12,atol=[1e-13,1e-14,1e-12],max_step=.05)
            assert result.success,result.message
            return result.y[:,-1]
        n=max(self.sub,int(np.ceil(abs(b-a)/(.12/self.sub))))
        h=(b-a)/n
        for k in range(n):
            x=a+k*h
            k1=self.rhs(x,y,i,B,outer);k2=self.rhs(x+h/2,y+h*k1/2,i,B,outer)
            k3=self.rhs(x+h/2,y+h*k2/2,i,B,outer);k4=self.rhs(x+h,y+h*k3,i,B,outer)
            y=y+h*(k1+2*k2+2*k3+k4)/6
        return y

    def branches(self,parameters,Bscale=1.,record=False):
        m=self.mat;B=m.B*Bscale;pc,rs,ms=parameters
        rs=m.R*np.exp(rs);ms=gr.TARGET*np.exp(ms)
        p,e,b=self.state(pc,len(m.lp)-1);r0=1.
        q0=4*np.pi*b*r0**3/(3*B)
        central_drop=2*np.pi/3*(e+p)*(e+3*p)*r0*r0
        yc=np.array([r0/m.R,(e/b)*q0,pc+np.log1p(-central_drop/p)])
        p,e,b=self.state(m.lp[0],0);f=1-2*ms/rs
        rB=np.sqrt(f)/(4*np.pi*rs*rs*b);massB=e/b*np.sqrt(f)
        lpB=-(e+p)*(ms+4*np.pi*rs**3*p)/(4*np.pi*rs**4*b*np.sqrt(f)*p)
        w0=min(m.dm[0]*1e-6,1e-8/abs(lpB*B))
        yo=np.array([(rs-rB*B*w0)/m.R,(ms-massB*B*w0)/B,m.lp[0]-lpB*B*w0])
        inner=[(q0,*yc)];outer=[(w0,*yo)]
        for i in range(len(m.lp)-1,m.split-1,-1):
            low=q0 if i==len(m.lp)-1 else m.inner[i+1]
            yc=self.step(np.log(low),np.log(m.inner[i]),yc,i,B,False)
            if record: inner.append((m.inner[i],*yc))
        for i in range(m.split):
            low=w0 if i==0 else m.outer[i]
            yo=self.step(np.log(low),np.log(m.outer[i+1]),yo,i,B,True)
            if record: outer.append((m.outer[i+1],*yo))
        return yc-yo,np.array(inner),np.array(outer)


def match(points=17,subdivision=1,family='fixed',direct=False,initial=None):
    solver=Solver(points,subdivision,direct);m=solver.mat;trials=[]
    prefix=family+('-direct' if direct else '')
    previous=OUT/f'{family}-{points}-{subdivision//2}.json'
    supplied=initial
    if supplied is not None: initial=supplied
    elif previous.exists(): initial=json.loads(previous.read_text())['parameters']
    elif (OUT/'fixed-17-1.json').exists(): initial=json.loads((OUT/'fixed-17-1.json').read_text())['parameters']
    else: initial=[m.lp[-1]+.005,0.,.0007]
    def objective(x):
        args=x if family=='fixed' else [x[0],x[1],0.]
        error,_,_=solver.branches(args,1. if family=='fixed' else np.exp(x[2]))
        row=dict(parameters=list(map(float,x)),residual=error.tolist());trials.append(row)
        print('match',prefix,points,subdivision,row,flush=True)
        return error
    fit=least_squares(objective,initial,diff_step=2e-5,xtol=1e-11,ftol=1e-11,gtol=1e-11,
        bounds=([m.lp[-1]-.4,-.3,-.02],[m.lp[-1]+.4,.3,.02]),max_nfev=40)
    args=fit.x if family=='fixed' else [fit.x[0],fit.x[1],0.]
    Bscale=1. if family=='fixed' else np.exp(fit.x[2])
    error,inner,outer=solver.branches(args,Bscale,record=True)
    save(f'{prefix}-{points}-{subdivision}.json',dict(classification='Counterexample candidate',
        parameters=fit.x.tolist(),interface_residual=error.tolist(),success=bool(fit.success),
        baryon_mass_g=float(m.baryon_g*Bscale),baryon_scale=float(Bscale),photospheric_radius_m=float(m.R*np.exp(fit.x[1])),
        photospheric_mass_GM_solar=float(.197536385307*np.exp(args[2])),
        target_mass_relative_residual=float(np.expm1(args[2])),EOS_calls=solver.calls,
        direct_EOS=direct,max_direct_entropy_error=solver.max_entropy_error))
    save(f'trials-{prefix}-{points}-{subdivision}.json',trials)
    np.savez_compressed(OUT/f'{prefix}-structure-{points}-{subdivision}.npz',inner=inner,outer=outer)
    assert fit.success and max(abs(error))<1e-8,error


def pilot(): match(17,1)


def refine():
    for points,sub in [(17,2),(17,4),(33,2),(33,4)]: match(points,sub)


def mass_family_plan():
    fixed=json.loads((OUT/'fixed-33-4.json').read_text())
    assert abs(fixed['target_mass_relative_residual'])>1e-8
    save('mass-family-plan.json',dict(classification='Counterexample candidate',
        fixed_material_result_sha256=sha(OUT/'fixed-33-4.json'),
        reason='The fixed original baryon inventory determines a gravitational mass different from the specified numerical target. This is not an observational exclusion.',
        family='Scale all cell baryon masses by a single positive factor, retaining their isotope fractions and entropy per baryon gram. Fit central pressure, optical radius and this scale, while imposing the target photospheric gravitational mass.',
        bounds='Baryon scale within exp(-0.02) to exp(0.02); no extra temperature or entropy parameter.',
        primary_unchanged=True,normalized_material_profile_preserved=True,absolute_original_inventory_preserved=False,
        atmosphere='Append and report the tiny atmosphere after the photospheric mass fit; require its mass contribution below the declared 1e-8 tolerance.',
        optical='Evaluate retained-luminosity Teff after the fit; Teff/radius is not a fit target in this family, but retained luminosity and material template remain external assumptions.'))


def mass_family():
    assert (OUT/'mass-family-plan.json').exists()
    fixed=json.loads((OUT/'fixed-33-4.json').read_text())['parameters']
    initial=[fixed[0],fixed[1],-fixed[2]]
    match(33,2,'scaled',initial=initial);match(33,4,'scaled')


def eos_checks():
    m=Material();eos=gr.EOS();table=np.load(OUT/'adiabats-33.npz');rows=[]
    for i in np.linspace(0,len(m.lp)-1,41).astype(int):
        for dx in [-.08125,.0375]:
            lp=m.lp[i]+dx;s0=table['reference'][i,3]
            a,lt,n=invert(eos,lp,s0,m.eps[i],m.lt[i]+.2*dx)
            cp=(a[10]-a[1]/a[0]*a[8])/np.exp(lt)
            ap=eos(1,lp,lt+1e-4,m.eps[i]);am=eos(1,lp,lt-1e-4,m.eps[i])
            fd=(ap[3]-am[3])/2e-4
            rows.append(dict(cell=int(i),offset=dx,relative_entropy_error=float(a[3]/s0-1),cp_relative_FD_error=float(fd/cp-1)))
    assert max(abs(x['relative_entropy_error']) for x in rows)<1e-10
    assert max(abs(x['cp_relative_FD_error']) for x in rows)<1e-5
    save('eos-controls.json',dict(classification='Proven',states=len(rows),rows=rows,
        interpretation='Sampled entropy inversion and independent entropy finite-difference heat capacity checks; not a global EOS error certificate.'))
    print('PASS EOS controls',max(abs(x['cp_relative_FD_error']) for x in rows),flush=True)


def direct_audit(write=True):
    result={}
    for family in ['fixed','scaled']:
        row=json.loads((OUT/f'{family}-33-4.json').read_text());x=row['parameters']
        args=x if family=='fixed' else [x[0],x[1],0.]
        solver=Solver(33,1,True,'DOP853');start=time.monotonic()
        error,_,_=solver.branches(args,row['baryon_scale'])
        result[family]=dict(interface_residual=error.tolist(),
            entropy_residual=solver.max_entropy_error,EOS_calls=solver.calls,elapsed_s=time.monotonic()-start,
            passed=bool(max(abs(error))<1e-8))
        print('direct audit',family,result[family],flush=True)
        if write:
            save('direct-audit.json',dict(classification='Proven',method='DOP853 with direct FreeEOS entropy inversion; independent of tabulated EOS and fixed-step RK4',results=result))
    assert all(row['passed'] for row in result.values())


def recheck_direct(): direct_audit(False)


def symbolic():
    import sympy as s
    r,m,p,e,b=s.symbols('r m p e b',positive=True);f=1-2*m/r
    Bprime=4*s.pi*r*r*b/s.sqrt(f)
    Pprime=-(e+p)*(m+4*s.pi*r**3*p)/(r*r*f)
    expected=[s.sqrt(f)/(4*s.pi*r*r*b),e/b*s.sqrt(f),
        -(e+p)*(m+4*s.pi*r**3*p)/(4*s.pi*r**4*b*s.sqrt(f)*p)]
    original=[1,4*s.pi*r*r*e,Pprime/p]
    assert all(s.simplify(a/Bprime-z)==0 for a,z in zip(original,expected))
    cx,T,cp=s.symbols('cx T cp',positive=True);satom=s.Function('s')(T)
    assert s.diff(cx*satom,T)==cx*s.diff(satom,T)
    deltaM,X,S,scale=s.symbols('deltaM X S scale',positive=True)
    assert s.simplify((scale*deltaM*X)/(scale*deltaM)-X)==0
    assert s.simplify((scale*deltaM*S)/(scale*deltaM)-S)==0
    save('symbolic-checks.json',dict(classification='Proven',baryon_coordinate_TOV=True,
        fixed_composition_entropy_normalization=True,uniform_mass_scaling_preserves_specific_not_total_inventories=True,
        limitation='Identities and boundary formulation only; no evolutionary, stability or dynamical scalar claim.'))
    print('PASS baryon-coordinate TOV and material identities')


def uniform_control():
    from types import SimpleNamespace
    R=1.;e=.03;b=.02;a=8*np.pi*e/3;surface=np.sqrt(1-a*R*R)
    def baryon(r): return 2*np.pi*b/a**1.5*(np.arcsin(np.sqrt(a)*r)-np.sqrt(a)*r*np.sqrt(1-a*r*r))
    B=baryon(R);errors=[]
    fake=SimpleNamespace(mat=SimpleNamespace(R=R),
        state=lambda lp,i:(gr.G*np.exp(lp)*.1/gr.C**4,e,b))
    for radius in np.linspace(.05,.95,19):
        f=1-a*radius*radius;z=np.sqrt(f);p=e*(z-surface)/(3*surface-z)
        mass=4*np.pi*e*radius**3/3;Br=baryon(radius)
        dp=-2*e*surface*a*radius/(z*(3*surface-z)**2)
        dr=Br/(4*np.pi*radius*radius*b/np.sqrt(f))
        expected=np.array([dr/R,4*np.pi*radius*radius*e*dr/B,dp/p*dr])
        value=Solver.rhs(fake,np.log(Br/B),np.array([radius/R,mass/B,np.log(p*gr.C**4/(gr.G*.1))]),0,B,False)
        errors.append(float(max(abs(value/expected-1))))
    assert max(errors)<1e-12
    save('uniform-control.json',dict(classification='Proven',
        control='Analytic constant-energy and constant-baryon-density Schwarzschild interior, including the independent analytic proper-volume integral.',
        samples=len(errors),max_relative_RHS_error=max(errors)))
    print('PASS analytic uniform-density baryon coordinate control',max(errors))


def atmosphere(row):
    mat=Material();a=gr.EOS()(1,mat.lp[0],mat.lt[0],mat.eps[0])
    ps=gr.G*np.exp(mat.lp[0])*.1/gr.C**4;rho=gr.G*a[0]*1000/gr.C**2
    bs=rho/mat.cx[0];u0=a[2]*1e-4/gr.C**2-1.5*ps/rho
    r=row['photospheric_radius_m'];mass=row['photospheric_mass_GM_solar']*gr.GM_SUN/gr.C**2
    B=row['baryon_mass_g']*.001*gr.G/gr.C**2
    def rhs(z,y):
        r,m,b=y;p=ps*z**2.5;den=rho*z**1.5
        e=den*(1+u0)+1.5*p;f=1-2*m/r
        dr=-2.5*(ps/rho)*r*r*f/((1+u0+2.5*ps/rho*z)*(m+4*np.pi*r**3*p))
        return [dr,4*np.pi*r*r*e*dr,4*np.pi*r*r*bs*z**1.5/np.sqrt(f)*dr]
    sol=solve_ivp(rhs,(1.,0.),[r,mass,B],method='DOP853',rtol=1e-12,atol=[1e-5,1e-12,1e-12],max_step=.02)
    assert sol.success
    R,M,NB=sol.y[:,-1];teff=float(mat.h['Teff'])*np.sqrt(mat.R/r)
    logg=float(np.log10(M*gr.C**2/r**2/np.sqrt(1-2*M/r)*100))
    result=dict(vacuum_radius_m=float(R),ADM_mass_GM_solar=float(M*gr.C**2/gr.GM_SUN),
        total_baryon_mass_g=float(NB*gr.C**2/gr.G*1000),
        extra_atmosphere_baryon_fraction=float((NB-B)/B),
        extra_atmosphere_mass_fraction=float((M-mass)/mass),
        atmospheric_thickness_m=float(R-r),conditional_Teff_K=float(teff),surface_logg_cgs=logg,
        optical_cut_passed=bool(15500<=teff<=16100 and 5.67<=logg<=5.97))
    assert abs(result['extra_atmosphere_mass_fraction'])<1e-8
    assert abs(result['extra_atmosphere_baryon_fraction'])<1e-8
    return result


def evidence():
    m=Material();isos=json.loads((ROOT/'outputs/thermal-restart19/network-diagnosis.json').read_text())['profile_isotopes']
    keep=[k for k in isos if not k.startswith('f')];norm=sum(m.d[k] for k in keep)
    X=np.array([m.d[k]/norm for k in keep]);dm=m.d['dm']
    assert np.max(abs(X.sum(axis=0)-1))<1e-14
    ref=np.load(OUT/'adiabats-33.npz')['reference'];sB=m.cx*ref[:,3]
    initial=X@dm;S=float(np.dot(dm,sB));atomic=float(np.dot(dm,m.cx));results={}
    for family in ['fixed','scaled']:
        row=json.loads((OUT/f'{family}-33-4.json').read_text());scale=row['baryon_scale']
        final=atmosphere(row);new=X@(dm*scale)
        error=np.max(abs(new/initial/scale-1))
        assert error<1e-14
        atomic_GM=atomic*.001*gr.G/gr.GM_SUN*scale
        results[family]=dict(**row,**final,
            original_cell_inventory_relative_error_after_scaling=float(error),
            isotope_inventory_g={k:dict(original=float(a),interior=float(b)) for k,a,b in zip(keep,initial,new)},
            original_interior_entropy_erg_K=S,interior_entropy_erg_K=S*scale,
            neutral_atomic_rest_mass_GM_solar=float(atomic_GM),
            ADM_minus_neutral_atomic_rest_mass_GM_solar=final['ADM_mass_GM_solar']-atomic_GM,
            source_baryon_mass_GM_solar=float(m.baryon_g*.001*gr.G/gr.GM_SUN),
            exact_original_interior_inventories_preserved=family=='fixed',
            temperature_or_entropy_parameter_fitted=False,
            limitation='Cell inventories/entropy are exact for the declared piecewise material template. Mathematical atmosphere adds separately reported tiny material. EOS/reference-template errors and actual reaction/heat evolution are not certified.')
    save('family-evidence.json',dict(classification='Counterexample candidate',results=results,
        source_material_cell_closure_error=float(np.max(abs(m.inner+m.outer-1)))))
    for k,v in results.items(): print(k,{a:v[a] for a in ['ADM_mass_GM_solar','baryon_scale','conditional_Teff_K','surface_logg_cgs','extra_atmosphere_baryon_fraction','ADM_minus_neutral_atomic_rest_mass_GM_solar']},flush=True)


def finalize():
    data=json.loads((OUT/'family-evidence.json').read_text())['results']
    audit=json.loads((OUT/'direct-audit.json').read_text())['results']
    assert all(v['passed'] and v['entropy_residual']<1e-10 for v in audit.values())
    fixed=data['fixed'];scaled=data['scaled']
    assert abs(fixed['ADM_mass_GM_solar']/.197536385307-1)>1e-8
    assert abs(scaled['ADM_mass_GM_solar']/.197536385307-1)<1e-8
    assert fixed['optical_cut_passed'] and scaled['optical_cut_passed']
    comparisons={}
    for a,b in [('fixed-17-1','fixed-17-2'),('fixed-17-2','fixed-17-4'),
                ('fixed-17-4','fixed-33-4'),('fixed-33-2','fixed-33-4'),('scaled-33-2','scaled-33-4')]:
        x=json.loads((OUT/(a+'.json')).read_text());y=json.loads((OUT/(b+'.json')).read_text())
        comparisons[a+' -> '+b]={k:float(y[k]/x[k]-1) for k in
            ['photospheric_radius_m','photospheric_mass_GM_solar','baryon_mass_g']}
    save('refinement-summary.json',dict(classification='Counterexample candidate',comparisons=comparisons,
        interpretation='Finite EOS-table and ODE-step refinements; not physical mesh or continuum certification. Original material cells are the declared piecewise model.'))
    save('gates.json',dict(classification='Counterexample candidate',
        original_interior_baryon_composition_and_reference_entropy_preserved_in_primary=True,
        primary_matches_numerical_gravitational_mass_target=False,
        separate_uniform_baryon_scale_matches_GR_target=True,
        scaled_family_preserves_specific_material_profile=True,
        scaled_family_preserves_original_total_inventory=False,
        additional_temperature_or_entropy_parameter_fitted=False,
        conditional_optical_cut_passes_without_optical_fit=True,
        direct_EOS_independent_integrator_passed=True,
        atmosphere_material_budget_reported=True,
        actual_MESA_absolute_entropy_reference_certified=False,
        full_22_isotope_reactive_microphysics=False,
        GR_thermal_composition_evolution_solved=False,
        stability_or_fluid_scalar_dynamics_certified=False,
        complete_nonlinear_observational_inference=False,final_PDF_or_ZIP_updated=False,
        classification_detail='Theorem progress: material-coordinate TOV and entropy/inventory identities. Loophole progress: conservative declared-material GR projection and separately scaled mass-matched candidate; no actual thermal evolution or dynamic observable established.'))


def maintain():
    import shutil
    prior=json.loads((OLD/'manifest.json').read_text())['sha256']
    additions={
        'model-definition':'분류: Counterexample candidate. Request23은 원래 5735개 구역의 바리온 질량, 불소 제외 재규격화 조성과 FreeEOS 기준 엔트로피를 물질 좌표에 고정하고 TOV를 다시 풀었다. 원래 내부 물질을 보존한 해는 지정 중력질량 목표보다 0.0698907% 높다. 별도 모형족에서 모든 구역 질량을 0.999301574445배로 조정하면 목표 GR 질량을 맞춘다. 이는 원래 총 물질을 보존한 변환과 다른 모형족이며 온도·엔트로피 조정 인자는 없다.',
        'observable-targets':'분류: Counterexample candidate. 별도 바리온 질량 조정 후보의 ADM 질량은 0.197536385307002 GM_sun 단위, 광학 반지름은 68819648.589 m, 유지한 광도를 조건으로 한 Teff는 15834.371 K, logg는 5.743136이다. 이 단계에서 광학량을 적합하지 않았으나 물질 템플릿과 광도가 외부 조건이므로 독립적인 실제 항성 예측이나 새로운 동적 관측량으로 세지 않는다.',
        'adiabatic-limit':'분류: Proven. 고유 바리온 질량을 독립 변수로 한 TOV 식과 고정 조성에서의 엔트로피 정규화를 기호 검산했다. 조성과 단위 바리온 질량당 엔트로피를 유지한 균일 질량 배율은 비율을 보존하지만 총 원소 재고와 총 엔트로피는 그 배율만큼 바꾼다. 정확한 정적 계량의 순 열유속 no-go와 기존 평탄 구각 scalar 오차 상계의 적용 경계는 그대로다.',
        'nonadiabatic-regime':'분류: Counterexample candidate. 직접 FreeEOS 엔트로피 역산과 독립 DOP853 적분은 두 모형의 연결 잔차를 1e-8 이내로 재확인했다. 이 검사는 반응·열수송을 끈 지정 물질의 정수압 재구성이다.\n\n분류: Conjectural. 다음 경계는 조성이 변할 때의 정지질량·핵반응·내부에너지 기준을 중복 없이 연결하고, 새 상태의 반응·손실·전도·복사·대류와 GR 열·조성 진화를 계산하는 것이다. 안정성 및 유체·metric·scalar 응답은 그 다음이며 정적 재구성을 궤도 완화나 관측 추론으로 확대하지 않는다.',
        'failure-ledger-dynamic-chi':'분류: Counterexample candidate. Request22에서 실패한 압력 고정 조성의 보존 해석은 원래 후보의 실패로 보존한다. 새 바리온 좌표 해는 원래 내부 재고와 선택한 EOS 엔트로피를 보존하지만 목표 중력질량과 일치하지 않는다. 별도로 등록한 균일 바리온 질량 조정 모형은 원래 총 재고를 약 0.0698426% 줄이는 다른 후보이며, 이를 원래 물질의 보존 변환이라고 부르지 않는다.\n\n분류: Conjectural. 보존한 엔트로피는 원래 P,T에서 정의한 FreeEOS 값으로 원래 MESA 전체 EOS의 절대 엔트로피 인증이 아니다. 수학적 외곽 대기의 추가 바리온, 구역별 상수 템플릿, 불소 생략과 열·조성 진화 미완료를 계속 명시한다.'}
    dest=OUT/'request22-notes';dest.mkdir(exist_ok=False);bindings={}
    for rel in ['docs/'+k+'.md' for k in additions]+['paper/revision-manifest.json']:
        assert sha(ROOT/rel)==prior[rel],rel
        snap=dest/Path(rel).name;shutil.copy2(ROOT/rel,snap)
        bindings[rel]=dict(snapshot=snap.relative_to(ROOT).as_posix(),sha256=prior[rel],historical_manifest=rel.startswith('paper/'))
    save('historical-note-bindings.json',bindings)
    revision=json.loads((ROOT/'paper/revision-manifest.json').read_text())
    for stem,body in additions.items():
        rel='docs/'+stem+'.md'
        with (ROOT/rel).open('ab') as f:
            f.write(('\n\n## Request 23 바리온 좌표와 엔트로피 보존 GR 재구성\n\n'+body+
                '\n\n세부 근거: [한글 보고서](../notes/REQUEST23_BARYON_ENTROPY_KO.md).\n').encode())
        revision['sha256'][rel]=sha(ROOT/rel)
    revision['request23_supporting_note_update']=dict(evidence_manifest='outputs/baryon-entropy23/manifest.json',
        historical_notes='outputs/baryon-entropy23/historical-note-bindings.json',
        status='물질·기준 엔트로피 보존 TOV 재구성 및 별도 바리온 배율 GR 질량 후보; 실제 열 진화 미완료',
        artifact_status='원고 PDF와 ZIP은 Request 12 동결본이다.')
    (ROOT/'paper/revision-manifest.json').write_text(json.dumps(revision,ensure_ascii=False,indent=2)+'\n')


def seal():
    import scipy,sympy
    save('provenance.json',dict(classification='Proven',before_task_checkpoint='729aca0',
        previous_manifest_sha256=sha(OLD/'manifest.json'),interpreter=sys.executable,
        versions={m.__name__:m.__version__ for m in [np,scipy,sympy]},
        runtime_provenance='outputs/gr-mass21/provenance.json',
        source_profile='outputs/thermal-closure22/selected.data.gz',producer='verification/baryon_entropy.py',
        entropy_source='outputs/gr-mass21/sources/src/mod_free_eos.f90: fixed-composition entropy and enthalpy derivatives',
        primary_reference=dict(title='Althaus et al. (2022), Structure and evolution of ultra-massive white dwarfs in general relativity',
            url='https://www.aanda.org/articles/aa/pdf/2022/12/aa44604-22.pdf',
            use='Primary reference for distinguishing conserved rest-mass sequences and time-dependent gravitational mass; current equations also independently derived.'),
        record_note='EOS_calls in solver records counts state-evaluation requests, including table lookups; a direct entropy inversion may require multiple low-level FreeEOS calls. Initial 17-point pilot predates optional record metadata. Timing fields are run records, not deterministic outputs.'))
    files=[p for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json']
    files += [ROOT/'verification/baryon_entropy.py',ROOT/'notes/REQUEST23_BARYON_ENTROPY_KO.md',ROOT/'.gitattributes']
    files += [ROOT/k for k in json.loads((OUT/'historical-note-bindings.json').read_text())]
    save('manifest.json',dict(classification='Proven',scope='Declared material and entropy preserving TOV projection plus separately scaled GR mass family',
        sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(files)}))


def verify():
    stages=['validated-variational','remaining-levers15','nbody-readout16','nonzero-drive17',
            'thermal-wd18','thermal-restart19','thermal-robustness20','gr-mass21','thermal-closure22','baryon-entropy23']
    histories={f'outputs/{a}/manifest.json':f'outputs/{b}/historical-note-bindings.json' for a,b in zip(stages,stages[1:])}
    count=0
    for label in ['outputs/research-remediation/manifest.json',*histories,'paper/revision-manifest.json','outputs/baryon-entropy23/manifest.json']:
        old=json.loads((ROOT/histories[label]).read_text()) if label in histories else {}
        for name,expected in json.loads((ROOT/label).read_text())['sha256'].items():
            path=ROOT/name
            if name in old:
                bind=old[name];path=ROOT/bind['snapshot'];assert expected==bind['sha256']
                if bind.get('historical_manifest'):
                    before=json.loads(path.read_text());after=json.loads((ROOT/name).read_text())
                    for k,value in before.items():
                        if k!='sha256': assert after[k]==value,k
                    for k,value in before['sha256'].items():
                        if k not in old: assert after['sha256'][k]==value,k
                else: assert (ROOT/name).read_bytes().startswith(path.read_bytes())
            assert sha(path)==expected,(label,name);count+=1
    plan=json.loads((OUT/'mass-family-plan.json').read_text())
    assert plan['fixed_material_result_sha256']==sha(OUT/'fixed-33-4.json')
    assert all(x['passed'] for x in json.loads((OUT/'direct-audit.json').read_text())['results'].values())
    gates=json.loads((OUT/'gates.json').read_text())
    assert gates['original_interior_baryon_composition_and_reference_entropy_preserved_in_primary']
    assert gates['separate_uniform_baryon_scale_matches_GR_target']
    for k in ['primary_matches_numerical_gravitational_mass_target','scaled_family_preserves_original_total_inventory',
              'additional_temperature_or_entropy_parameter_fitted','GR_thermal_composition_evolution_solved',
              'complete_nonlinear_observational_inference','final_PDF_or_ZIP_updated']: assert not gates[k],k
    provenance=json.loads((ROOT/'outputs/gr-mass21/provenance.json').read_text())
    for key in ['runtime_sha256','module_sha256']:
        for path,digest in provenance[key].items(): assert sha(path)==digest,path
    print('PASS:',count,'현재·역사 SHA, 원래 물질 보존 해와 별도 질량 조정의 경계')


if __name__=='__main__': globals()[sys.argv[1]]()
