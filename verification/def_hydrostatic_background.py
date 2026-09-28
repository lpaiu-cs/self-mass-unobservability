"""Fixed-inventory finite-temperature DEF hydrostatic reconstruction.

The optical boundary retains its prescribed finite pressure. It is not a
vacuum material surface or a stationary radiating atmosphere. The solved
profile is a mechanical initial state; thermal evolution is still required.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import argparse
import hashlib
import json
import time
import numpy as np
import mpmath as mp
import sympy as sp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import root
import gr_molecular_conservative_initial as initial
import gr_scalar_nonlinear_exterior as exterior

molecular=initial.molecular
gr=molecular.g.c.gr
ROOT=Path(__file__).resolve().parents[1]
OLD=molecular.OUT
OUT=molecular.g.OUT/'def-hydrostatic-background'


def digest(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d): Path(p).write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n')


def inputs():
    return dict(np.load(OLD/'molecular-state-17-8.npz')),dict(np.load(OLD/'molecular-adiabats-17.npz'))


def coverage(data,table):
    rho=np.exp(data['lnd']);p=np.exp(data['logP'])
    h=data['CX']*(gr.C*100)**2+data['u_W']+p/rho
    drop=abs(np.diff(data['nu_faces']))*h*rho/p
    # This is only a proposed lookup domain, never an EOS extrapolation.
    span=np.maximum(.04,1.25*drop+.015)
    return [(int(i),np.arange(-int(np.ceil(span[i]/.005)),int(np.ceil(span[i]/.005))+1)*.005)
            for i in np.where(span>.040000001)[0]]


def symbolic():
    r,m,p,en,beta,phi,v=sp.symbols('r m p en beta phi v',nonzero=True)
    b=1-2*m/r
    mr=4*sp.pi*r*r*en+r*r*b*v*v/2
    nr=m/(r*r*b)+4*sp.pi*r*p/b+r*v*v/2
    lr=(mr/r-m/r**2)/b
    wave=-(2/r+nr-lr)*v+4*sp.pi*beta*phi*(en-3*p)/b
    source=4*sp.pi/b*(beta*phi*(en-3*p)+r*v*(en-p))-2*(r-m)/(r*r*b)*v
    assert sp.simplify(wave-source)==0
    rho,A,B=sp.symbols('rho A B',positive=True)
    drdB=sp.sqrt(b)/(4*sp.pi*r*r*A**3*rho)
    assert sp.simplify(drdB*(4*sp.pi*r*r*A**3*rho/sp.sqrt(b)))==1
    return dict(classification='Proven',passed=True,scope='DEF equations (3.6a-f), log lapse is half the source nu; fixed Jordan baryon mass coordinate. No thermal stationarity assertion.')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    data,table=inputs();targets=coverage(data,table)
    count=sum(np.sum(abs(x)>.040000001) for _,x in targets)
    paths=[Path(__file__),OLD/'molecular-state-17-8.npz',OLD/'molecular-adiabats-17.npz',
           OLD/'GR-8.json',OLD/'candidate.py',ROOT/'verification/gr_molecular_conservative_initial.py',
           ROOT/'verification/baryon_entropy.py',ROOT/'verification/gr_scalar_nonlinear_exterior.py']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='972195f',
        bindings={p.relative_to(ROOT).as_posix():digest(p) for p in paths},symbolic=symbolic(),
        source='https://arxiv.org/pdf/gr-qc/9602056 (3.6a-f); finite-temperature fixed shell entropy/composition EOS here, not the source cold barotrope.',
        beta=-4,phi_infinity=[0,.001],subdivision=4,cells=len(data['dm']),
        material='Every saved shell baryon mass, 26-isotope fraction and EOS entropy is retained. Radius, mass, pressure and scalar profile are solved together; no subtracted hydrostatic reference force.',
        boundary='Prescribed original finite Jordan photospheric pressure. Just map normalizes the truncated metric/scalar only; no solved physical atmosphere or vacuum material junction.',
        table_extension_cells=len(targets),table_extension_roots=int(count),
        gates=dict(interface=1e-8,GR_reproduction_relative=1e-8,native_logrho=1e-7,native_entropy_over_cv=1e-7),
        budget=dict(pilot_seconds=90,table_seconds=180,shoot_seconds=600,native_audit_seconds=180,
                    workers=4,maximum_objective_calls_per_background=45,automatic_expansion=False),
        decision='Close the nonzero mechanical-background gap if matching and native controls pass; retain finite-pressure/heat limitations. Do not infer stationarity over an orbit or execute an orbit.'))


def eos_worker(task):
    i,offsets,tag=task;data,table=inputs();eos=molecular.model.EOS()
    stats=dict(calls=0,evaluations=0,maximum_score=0.,label=tag)
    inv=molecular.inverse(stats,OUT/tag)
    ref=table['reference'][i];lp=np.log(ref[1]);entropy=ref[3]
    poly=PchipInterpolator(table['offset'],table['values'][i],axis=0)
    slope=poly.derivative()(0)[1];rows=[];start=time.monotonic()
    for dx in offsets:
        a,t,_=inv(eos,lp+float(dx),entropy,data['X'][i],float(table['lnT'][i]+slope*dx))
        rows.append([np.log(a[0]),t,a[2]*1e-4/gr.C**2])
    return dict(cell=i,offsets=list(map(float,offsets)),values=rows,statistics=stats,seconds=time.monotonic()-start)


def pilot():
    assert not (OUT/'pilot.json').exists()
    data,table=inputs();targets=coverage(data,table)
    jobs=[(i,[float(x[0]),float(x[-1])],f'pilot-{i:04}') for i,x in targets[:3]]
    start=time.monotonic()
    with ProcessPoolExecutor(max_workers=4) as pool: rows=list(pool.map(eos_worker,jobs))
    write(OUT/'pilot.json',dict(classification='Counterexample candidate',rows=rows,seconds=time.monotonic()-start))
    print(json.dumps(dict(seconds=time.monotonic()-start,rows=rows)),flush=True)


def table_build():
    assert not (OUT/'extended-table.npz').exists()
    plan=json.loads((OUT/'plan.json').read_text());pilot=json.loads((OUT/'pilot.json').read_text())
    data,table=inputs();targets=coverage(data,table)
    rate=sum(r['seconds'] for r in pilot['rows'])/6
    projection=rate*plan['table_extension_roots']/4
    assert projection<120,('Table budget exceeded',projection)
    write(OUT/'table-budget.json',dict(classification='Counterexample candidate',measured_inverse_seconds=rate,
          projected_seconds=projection,hard_timeout_seconds=180,scope='Outer-cell pilot; other states unmeasured. No automatic expansion.'))
    jobs=[(i,[float(v) for v in x if abs(v)>.040000001],f'table-{i:04}') for i,x in targets]
    start=time.monotonic()
    with ProcessPoolExecutor(max_workers=4) as pool: rows=list(pool.map(eos_worker,jobs))
    arrays={}
    for row in rows:
        i=row['cell'];xx=np.r_[table['offset'],row['offsets']];yy=np.r_[table['values'][i],row['values']]
        order=np.argsort(xx);arrays[f'x{i}']=xx[order];arrays[f'y{i}']=yy[order]
    np.savez_compressed(OUT/'extended-table.npz',**arrays)
    write(OUT/'table-result.json',dict(classification='Counterexample candidate',rows=rows,seconds=time.monotonic()-start,
          physical_EOS_certified=False,continuous_table_error_certified=False))
    print('TABLE',len(rows),'seconds',time.monotonic()-start,flush=True)


class Structure(molecular.g.c.be.Solver):
    """Same two baryon branches/RK4; five actual DEF unknowns, guarded EOS table."""
    def __init__(self,phi0,sub=4):
        self.data,tab=inputs();d=self.data
        self.sub=sub;self.method='RK4';self.phi0=phi0;self.beta=-4.;self.calls=0
        self.R=float(d['radius_faces_m'][0]);self.M=float(d['mass_faces_geom'][0]);self.mu=self.M/self.R
        self.B=float(d['dm'].sum())*.001*gr.G/gr.C**2;self.f=d['dm']/d['dm'].sum()
        self.outer=np.r_[0,np.cumsum(self.f)];self.inner=np.r_[np.cumsum(self.f[::-1])[::-1],0]
        self.outer[-1]=self.inner[0]=1.;self.split=int(np.argmin(abs(self.outer-.5)))
        self.ref=tab['reference'];self.lp=np.log(self.ref[:,1]);self.offset=tab['offset']
        self.coef=PchipInterpolator(self.offset,tab['values'],axis=1).c
        ext=np.load(OUT/'extended-table.npz');self.extra={}
        for key in ext.files:
            if key.startswith('x'):
                i=int(key[1:]);xx=ext[key];self.extra[i]=(xx,PchipInterpolator(xx,ext[f'y{i}'],axis=0).c)

    def state(self,lp,i):
        self.calls+=1;dx=lp-self.lp[i]
        if i in self.extra: xx,co=self.extra[i]
        else: xx,co=self.offset,self.coef[:,:,i]
        if not xx[0]<=dx<=xx[-1]:raise ValueError(('No EOS extrapolation',i,float(dx),float(xx[0]),float(xx[-1])))
        j=min(len(xx)-2,max(0,int(np.searchsorted(xx,dx)-1)));h=dx-xx[j];a=co[:,j]
        v=((a[0]*h+a[1])*h+a[2])*h+a[3]
        b=gr.G*np.exp(v[0])*1000/gr.C**2
        return gr.G*np.exp(lp)*.1/gr.C**4,b*(self.data['CX'][i]+v[2]),b,v

    def radial(self,y,i):
        r=y[0]*self.R;m=y[1]*self.B;phi=self.phi0*(1+self.mu*y[3]);v=self.phi0*self.mu*y[4]/self.R
        p,en,rho,_=self.state(y[2],i);A=np.exp(self.beta*phi*phi/2);pE,enE=A**4*p,A**4*en
        b=1-2*m/r
        if r<=0 or m<=0 or b<=0:raise ValueError(('Invalid radial state',i,r,m,b))
        drdB=np.sqrt(b)/(4*np.pi*r*r*A**3*rho)
        nr=m/(r*r*b)+4*np.pi*r*pE/b+r*v*v/2
        zr=4*np.pi/b*(self.beta*(1+self.mu*y[3])*(enE-3*pE)*self.R/self.mu+r*y[4]*(enE-pE))-2*(r-m)/(r*r*b)*y[4]
        return np.array([1/self.R,(4*np.pi*r*r*enE+r*r*b*v*v/2)/self.B,
            -(en+p)/p*(nr+self.beta*phi*v),y[4]/self.R,zr])*drdB

    def rhs(self,x,y,i,B,outer):return self.radial(y,i)*B*np.exp(x)*(-1 if outer else 1)

    def branches(self,parameters,record=False):
        pc,lr,lm,uc,j=parameters;R=self.R*np.exp(lr);M=self.M*np.exp(lm);r0=1.
        p,en,rho,_=self.state(pc,len(self.f)-1);phic=self.phi0*(1+self.mu*uc);A=np.exp(self.beta*phic**2/2)
        pE,enE=A**4*p,A**4*en;Phi=4*np.pi/3*self.beta*(1+self.mu*uc)*(enE-3*pE)*r0
        q0=4*np.pi*A**3*rho*r0**3/(3*self.B)
        drop=2*np.pi/3*(enE+3*pE)*r0*r0+self.beta*phic*self.phi0*Phi*r0/2
        yc=np.array([r0/self.R,4*np.pi/3*enE*r0**3/self.B,pc-(en+p)/p*drop,uc+Phi*r0/(2*self.mu),Phi*self.R/self.mu])
        G=float(exterior.exact(mp.mpf(M/R),mp.mpf(self.phi0*self.mu*j))[2])
        ys=np.array([R/self.R,M/self.B,float(self.data['boundary_logP']),-j*G,j*self.R/R])
        derivative=self.radial(ys,0);w0=min(self.f[0]*1e-6,1e-8/abs(derivative[2]*self.B))
        yo=ys-derivative*self.B*w0
        inner=[(q0,*yc)];outer=[(w0,*yo)]
        for i in range(len(self.f)-1,self.split-1,-1):
            low=q0 if i==len(self.f)-1 else self.inner[i+1]
            yc=self.step(np.log(low),np.log(self.inner[i]),yc,i,self.B,False)
            if record:inner.append((self.inner[i],*yc))
        for i in range(self.split):
            low=w0 if i==0 else self.outer[i]
            yo=self.step(np.log(low),np.log(self.outer[i+1]),yo,i,self.B,True)
            if record:outer.append((self.outer[i+1],*yo))
        return yc-yo,np.asarray(inner),np.asarray(outer),ys


def shoot_pilot():
    assert not (OUT/'shoot-pilot.json').exists();s=Structure(0);d=s.data
    x=np.array([json.loads((OLD/'GR-8.json').read_text())['parameters'][0],0,0,s.beta*d['nu_faces'][-1]/s.mu,-4])
    start=time.monotonic();err,*_=s.branches(x);seconds=time.monotonic()-start
    write(OUT/'shoot-pilot.json',dict(classification='Counterexample candidate',seconds=seconds,parameters=x.tolist(),residual=err.tolist(),lookups=s.calls))
    print('SHOOT PILOT',seconds,err,flush=True)


def run():
    assert not (OUT/'result.json').exists();plan=json.loads((OUT/'plan.json').read_text())
    for rel,h in plan['bindings'].items():assert digest(ROOT/rel)==h,rel
    pilot=json.loads((OUT/'shoot-pilot.json').read_text());estimate=pilot['seconds']*40
    assert estimate<500,('Shooting budget exceeded',estimate)
    write(OUT/'shoot-budget.json',dict(classification='Counterexample candidate',single_branch_seconds=pilot['seconds'],
          projected_40_branches_seconds=estimate,hard_timeout_seconds=600,maximum_calls_per_background=45,automatic_expansion=False))
    x=np.array(pilot['parameters']);results=[];start=time.monotonic()
    for phi0 in plan['phi_infinity']:
        s=Structure(phi0);history=[]
        def objective(v):
            if len(history)>=45:raise RuntimeError('Shooting call cap')
            error,*_=s.branches(v);history.append(dict(parameters=v.tolist(),residual=error.tolist()))
            write(OUT/f'progress-{phi0}.json',dict(classification='Counterexample candidate',history=history,seconds=time.monotonic()-start))
            return error
        fit=root(objective,x,method='hybr',options=dict(xtol=2e-10,eps=1e-10,maxfev=45));x=fit.x
        error,inner,outer,ys=s.branches(x,record=True)
        assert abs(error).max()<plan['gates']['interface'],(fit.message,error)
        n=len(s.f);faces=np.zeros((n+1,5));faces[0]=ys;faces[1:s.split+1]=outer[1:,1:];faces[s.split:n]=inner[1:,1:][::-1]
        faces[n]=[0,0,x[0],x[3],0]
        # Midpoints use the same material branch, not a face interpolation.
        states=[];thermo=[]
        for i in range(n):
            outside=i<s.split
            if outside:begin=outer[i,0];mid=(s.outer[i]+s.outer[i+1])/2;value=outer[i,1:]
            else:j=n-1-i;begin=inner[j,0];mid=(s.inner[i]+s.inner[i+1])/2;value=inner[j,1:]
            value=s.step(np.log(begin),np.log(mid),value,i,s.B,outside)
            states.append(value);thermo.append(s.state(value[2],i)[3])
        states=np.asarray(states);thermo=np.asarray(thermo)
        np.savez_compressed(OUT/f'background-{phi0}.npz',faces=faces,states=states,thermo=thermo,dm=s.data['dm'],X=s.data['X'],entropy=s.ref[:,3],parameters=x)
        R=ys[0]*s.R;M=ys[1]*s.B;q=phi0*s.mu*x[4];ext=exterior.exact(mp.mpf(M/R),mp.mpf(q));ADM=M+R*q*q*float(ext[0])
        row=dict(classification='Counterexample candidate',phi_infinity=phi0,parameters=x.tolist(),interface=error.tolist(),
             objective_calls=len(history),EOS_lookups=s.calls,radius_m=R,mass_geom_m=M,ADM_geom_m=ADM,
             alpha_over_phi_infinity=s.mu*x[4]*R*float(ext[1])/ADM,table_sha256=digest(OUT/'extended-table.npz'),
             fixed_inventory=True,thermal_stationarity=False,physical_atmosphere=False)
        if phi0==0:
            row['GR_reproduction_relative']=max(abs(R/s.R-1),abs(M/s.M-1),float(abs(states[:,2]-s.data['logP']).max()))
            assert row['GR_reproduction_relative']<plan['gates']['GR_reproduction_relative'],row
        results.append(row);write(OUT/f'result-{phi0}.json',row);print('BACKGROUND',json.dumps(row),flush=True)
    write(OUT/'result.json',dict(classification='Counterexample candidate',mechanical_matching_passed=True,rows=results,seconds=time.monotonic()-start,native_audit_passed=False,full_dynamic_charge_solved=False))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','pilot','table_build','shoot_pilot','run'])
    globals()[p.parse_args().action]()
