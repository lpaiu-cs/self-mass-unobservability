"""Counterexample candidate: actual native cells resolve the material junction.

Split only the outermost deep cell; retain fifteen original native cells.
All shared material/photon evolution is inherited from the corrected Phase118
owner. This is a local radial comparison, not global spatial certification.
"""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import inspect
import json
import signal
import sys
import textwrap
import time
import numpy as np
from scipy.interpolate import PchipInterpolator
import def_native_material_centered as prior

feedback=prior.previous;old=prior.old;chem=feedback.chem;angular=old.deep.prior
C=prior.C;write=prior.write;sha=prior.sha;replace=prior.replace
OUT=prior.OUT.parent.parent/'def-native-boundary-layer'


def prepare():
    assert not OUT.exists();OUT.mkdir();(OUT/'thermal-support').mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='bf060cff9',
        claim='Replace the last34km center-to-material-face reconstruction by actual native radiative/reactive moving cells, then decide its influence on the conditional final charge.',
        grid='Keep the original first15 native cells exactly. Replace only the last cell by four physical cells with additional faces at-4000,-2000,-1000m and the same shared face at-400m. Their native centers are-36375,-3000,-1500,-700m. The nearest center-face gap becomes300m. Atmosphere remains512 cells; same152 frequencies,8 angles and3.434ms.',
        reuse='Existing15-cell geometry, frozen local inventories, nonlinear thermal/rate tables and density/frequency derivatives. Existing completed16-cell64/128 trajectories are the comparator; never repeat them. All common photon, central material and shared HLL source methods are reused.',
        native='Only four new local native inventory anchors, nonlinear H/temperature support identical to the old bank, same density/frequency derivatives and independent withheld native checks. Endpoint derivative reference states use interpolated original frozen references, not new freely adjusted EOS parameters.',
        decision='Apply the new bank to actual19-cell64/128 coupled paths. Compare direct charge, direct-plus-photon mass, shared material flux and surface spectrum against saved16-cell paths. A failed2percent local radial comparison stops automatic refinement; it does not become a new accepted physical sign.',
        gates=dict(native=.002,energy=1e-8,baryon=1e-10,source=1e-9,time_direct=.02,time_total=.02,space_direct=.02,space_total=.02,space_surface_spectrum=.02,positive_photons=True),
        budget=dict(native_calls=1800,native_seconds=75,pilot_steps_each=2,pilot_seconds=35,production_seconds=500,saved_readout_seconds=60,CPU_threads=1,memory_GB=3),
        forecast='The corrected16-cell two-path dispatch took323.75s. Increasing only bulk16 to19 cells is expected to add cost to the gas Schur solve; actual later stiffness and new short cells are unmeasured. Measure two steps of both clocks, require1.8x remaining forecast plus10s within500s, then continue saved prefixes only.',
        stop='Stop on native support, positivity, conservation, forecast or hard cap. No extra radial/time/frequency/angle path or relaxed gate. This local refinement does not settle the remaining fifteen-cell radial or initial Einstein-constraint error.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(prior.__file__),Path(prior.base.__file__),Path(feedback.__file__),Path(old.__file__),Path(angular.__file__),prior.OUT/'coupled-64.npz',prior.OUT/'coupled-128.npz',prior.OUT/'result.json']}))


# Rebind only bank paths and constructor owners; retain the actual equations.
thermal_init=textwrap.dedent(inspect.getsource(angular.ThermalTable.__init__))
thermal_init=replace(thermal_init,"prior.OUT/'bank-16-8.npz'","OUT/'geometry.npz'")
tns=dict(vars(angular),OUT=OUT);exec(compile(thermal_init,__file__,'exec'),tns)


class Thermal(angular.ThermalTable):
    __init__=tns['__init__']


angular_init=FunctionType(angular.Model.__init__.__code__,dict(vars(angular),OUT=OUT,ThermalTable=Thermal),argdefs=angular.Model.__init__.__defaults__)
unsplit_init=textwrap.dedent(inspect.getsource(old.deep.Model.__init__))
unsplit_init=replace(unsplit_init,'super().__init__(angles)','angular_init(self,angles)')
uns=dict(vars(old.deep),angular_init=angular_init);exec(compile(unsplit_init,__file__,'exec'),uns)


class Table(feedback.Table):
    __init__=FunctionType(feedback.Table.__init__.__code__,dict(vars(feedback),OUT=OUT,motion=SimpleNamespace(OUT=OUT)))


class Bulk(feedback.Bulk):
    def __init__(self,angles=8):
        uns['__init__'](self,angles);self.eos=Table(self.eos)
        self.extra=np.zeros_like(self.initial);self.extra_u=np.zeros(self.n);self.extra_y=np.zeros(self.n)


def mechanics(b):
    n=b.n;d=b.d;average=np.zeros((n,n+1));average[np.arange(n),np.arange(n)]=.5;average[np.arange(n),np.arange(1,n+1)]=.5
    return SimpleNamespace(xi=average/(4*np.pi*d['r']**2*d['B']*d['rho'])[:,None],K=np.load(OUT/'native.npz')['K'])


base_coupled_init=FunctionType(old.Coupled.__init__.__code__,dict(vars(old),deep=SimpleNamespace(**dict(vars(old.deep),Model=Bulk))),argdefs=old.Coupled.__init__.__defaults__)
feedback_init=textwrap.dedent(inspect.getsource(feedback.Coupled.__init__))
feedback_init=replace(feedback_init,'super().__init__(448,8)','base_coupled_init(self,448,8)')
feedback_init=replace(feedback_init,'self.mech=motion.Mechanics(448)','self.mech=mechanics(b)')
fns=dict(vars(feedback),OUT=OUT,Bulk=Bulk,base_coupled_init=base_coupled_init,mechanics=mechanics,motion=SimpleNamespace(OUT=OUT))
exec(compile(feedback_init,__file__,'exec'),fns)
join_init=textwrap.dedent(inspect.getsource(prior.base.Coupled.__init__))
join_init=replace(join_init,'super().__init__()','feedback_init(self)')
jns=dict(vars(prior.base),feedback_init=fns['__init__']);exec(compile(join_init,__file__,'exec'),jns)


class Coupled(prior.Coupled):
    __init__=jns['__init__']
    run=FunctionType(prior.Coupled.run.__code__,dict(prior.Coupled.run.__globals__,OUT=OUT),argdefs=prior.Coupled.run.__defaults__)


def bank():
    assert not (OUT/'bank.json').exists();start=time.monotonic();signal.signal(signal.SIGALRM,old.optical.timeout);signal.alarm(75)
    bg=chem.prior.Background();baseline=prior.Coupled();bd=baseline.bulk.d;n=19;keep=15
    original=np.load(chem.prior.OUT/'bank-16-8.npz');oldthermal=np.load(angular.OUT/'thermal-support/bank.npz')
    oldnative=np.load(feedback.motion.OUT/'native.npz');oldspec=np.load(feedback.OUT/'bank.npz');ref=np.load(feedback.motion.OUT/'forcing-896.npz')
    # Last unshortened edge is the surface, as required by the common photon
    # constructor. Its exact same final400m truncation is retained below.
    edges=np.r_[original['edges'][:-1],bg.R+np.array([-400000.,-200000.,-100000.,0.])]
    actual_edges=edges.copy();actual_edges[-1]=baseline.m.rf[0]
    r=(actual_edges[:-1]+actual_edges[1:])/2;r[:keep]=original['r'][:keep]
    d=dict(bg.sample(r));zf=bg.sample(edges)
    d.update(edges=edges,face_a=zf['a'],face_B=zf['B'],face_T=zf['T'],face_L=bg.luminosity(zf),Einf=original['Einf'],num=original['num'],cx=original['cx'])
    for key in ['r','a','B','rho','T','phi']:d[key][:keep]=original[key][:keep]
    for key in ['coeff','thermo','raw','y0','target']:
        d[key]=np.zeros((n,)+original[key].shape[1:]);d[key][:keep]=original[key][:keep]
    raw=np.zeros((n,)+oldthermal['raw'].shape[1:]);rates=np.zeros((n,)+oldthermal['rates'].shape[1:]);raw[:keep]=oldthermal['raw'][:keep];rates[:keep]=oldthermal['rates'][:keep]
    mechanics_raw=np.zeros((2,n,5,21));mechanics_raw[:,:keep]=oldnative['raw'][:,:keep]
    spectral={key:np.zeros((2,n)+oldspec[key].shape[2:]) for key in ['density','frequency']}
    for key in spectral:spectral[key][:,:keep]=oldspec[key][:,:keep]
    offsets=oldthermal['offsets'];ratios=oldthermal['ratios'];h=1e-4;native=chem.old.Native(cap=1800);derivative=[]
    ref_theta=np.vstack([np.interp(r,original['r'],ref['theta'][it]) for it in [0,-1]])
    ref_eta=np.vstack([np.interp(r,original['r'],ref['eta'][it]) for it in [0,-1]])
    for j in range(keep,n):
        rho,T,a=d['rho'][j],d['T'][j],d['a'][j];lt=np.log(T)
        eq=native.ion.snapshot(np.log(rho),lt,np.zeros(318));target=eq['number_fractions'];d['target'][j]=target;d['y0'][j]=target[0,0]/target[0,:2].sum()
        chem.setup(native,d,j);y=native.y0;states=[native.state(0.,lt+dt,y*(1+dy)) for dt,dy in [(0,0),(h,0),(-h,0),(0,h),(0,-h)]]
        rs=np.array([s['raw'] for s in states]);spec=np.array([chem.prior.coefficients(native,s,d['Einf']/a) for s in states]);ct=(spec[1]-spec[2])/(2*h);cy=(spec[3]-spec[4])/(2*h)
        bb=1/np.expm1(d['Einf']/(a*chem.prior.K*T));d['coeff'][j]=[spec[0,0],ct[1]-ct[0]*bb,cy[1]-cy[0]*bb,ct[0],cy[0]]
        dt=(rs[1]-rs[2])/(2*h);dy=(rs[3]-rs[4])/(2*h)
        d['thermo'][j]=[dt[2],dy[2],dt[1],dy[1],native.nH,rs[0,13]/1.66053906660e-24,dt[13]/1.66053906660e-24,dy[13]/1.66053906660e-24]
        d['raw'][j]=rs[0];derivative.append(abs(dt[2]/(T*dt[3])-1))
        for iy,ratio in enumerate(ratios):
            for it,off in enumerate(offsets):
                yy=y*ratio;s=native.state(0.,lt+off,yy);chi,em=chem.prior.coefficients(native,s,d['Einf']/a)
                raw[j,iy,it]=s['raw'];rates[j,iy,it]=[(chi+em)/yy,em/(1-yy)]
        for it in [0,1]:
            yy=y*(1+ref_eta[it,j]);tt=lt+ref_theta[it,j]
            mm=[native.state(x,tt+dt,yy) for x,dt in [(0,0),(h,0),(-h,0),(0,h),(0,-h)]]
            mechanics_raw[it,j]=[s['raw'] for s in mm]
            v=[]
            for s in [mm[2],mm[0],mm[1]]:
                chi,em=chem.prior.coefficients(native,s,d['Einf']/a);v.append(np.array([chi+em,em]))
            freq=[]
            for off in [-h,h]:
                chi,em=chem.prior.coefficients(native,mm[0],d['Einf']/a*np.exp(off));freq.append(np.array([chi+em,em]))
            mask=v[1]>1e-280
            spectral['density'][it,j]=np.divide(v[2]-v[0],2*h*v[1],out=np.zeros_like(v[1]),where=mask)
            spectral['frequency'][it,j]=np.divide(freq[1]-freq[0],2*h*v[1],out=np.zeros_like(v[1]),where=mask)
    d['native_calls']=native.ion.calls;d['derivative_errors']=np.array(derivative);np.savez_compressed(OUT/'geometry.npz',**d)
    np.savez_compressed(OUT/'thermal-support/bank.npz',raw=raw,rates=rates,offsets=offsets,ratios=ratios,native_calls=native.ion.calls)
    pr=(mechanics_raw[:,:,1,1]-mechanics_raw[:,:,2,1])/(2*h);ur=(mechanics_raw[:,:,1,2]-mechanics_raw[:,:,2,2])/(2*h)
    pt=(mechanics_raw[:,:,3,1]-mechanics_raw[:,:,4,1])/(2*h);ut=(mechanics_raw[:,:,3,2]-mechanics_raw[:,:,4,2])/(2*h)
    K=pr+pt*(mechanics_raw[:,:,0,1]/d['rho']-ur)/ut;assert np.min(ut)>0 and np.min(K)>0
    np.savez_compressed(OUT/'native.npz',raw=mechanics_raw,K=K,pr=pr,ur=ur,pt=pt,ut=ut,native_calls=native.ion.calls)
    np.savez_compressed(OUT/'bank.npz',**spectral,native_calls=native.ion.calls)
    eos=Thermal();p0,u0,*_=eos.gas(np.zeros(n),np.zeros(n));geo=feedback.prior.green.Geometry(baseline.m)
    _,af,Bf,_,_=geo.metric(actual_edges-baseline.m.RJ);_,a,B,re,phi=geo.metric(r-baseline.m.RJ);z=baseline.m.bg.sample(re)
    A=np.exp(-2*phi**2);bb=1-2*z['m']/re;alpha=-4*phi;den=1+alpha*re*z['v']
    ap=a*(z['m']/(re*re*bb)+4*np.pi*re*A**4*z['p']/bb+re*z['v']**2/2+alpha*z['v'])/(baseline.m.R*A*den)
    lt=PchipInterpolator(bg.r,np.log(bg.d['temperature_K']));Prad=float(bg.env['Prad'][0]*3/bg.env['T'][0]**4)*d['T']**4/3
    rho_prime=d['rho']*PchipInterpolator(bg.r,np.log(bg.d['density_cgs']))(r,1)
    pgprime=-(d['rho']*(float(d['cx'])*C*C+u0)+p0+4*Prad)*ap/a-4*Prad*lt(r,1)
    np.savez_compressed(OUT/'forcing-896.npz',r=r,edges=actual_edges,volume=4*np.pi*r*r*B*np.diff(actual_edges),a=a,B=B,ap=ap,af=af,Bf=Bf,
        rho=d['rho'],rho_face=bg.sample(actual_edges)['rho'],p0=p0,u0=u0,cx=d['cx'],pgprime=pgprime,rho_prime=rho_prime,uprime=PchipInterpolator(r,u0)(r,1),
        initial_support=a/B*(4*Prad*lt(r,1)+4*Prad*ap/a),theta=ref_theta,eta=ref_eta)
    checks=[]
    for j in range(keep,n):
        chem.setup(native,d,j);theta=.09;eta=-.95;s=native.state(0.,np.log(d['T'][j])+theta,d['y0'][j]*(1+eta));chi,em=chem.prior.coefficients(native,s,d['Einf']/d['a'][j])
        pp,uu,*_=eos.gas(np.full(n,theta),np.full(n,eta));aa,ee,*_=eos.radiation(np.full(n,theta),np.full(n,eta));mask=(chi+em>1e-280)&(em>1e-280)
        checks.append(dict(cell=j,constitutive=float(max(abs(pp[j]/s['raw'][1]-1),abs(uu[j]/s['raw'][2]-1))),spectral=float(max(np.max(abs(aa[j,mask]/(chi+em)[mask]-1)),np.max(abs(ee[j,mask]/em[mask]-1))))))
    for key in ['r','a','B','rho','T','phi','raw','y0','thermo','target']:assert np.array_equal(d[key][:keep],original[key][:keep]),key
    assert np.array_equal(raw[:keep],oldthermal['raw'][:keep]) and np.array_equal(rates[:keep],oldthermal['rates'][:keep])
    row=dict(classification='Counterexample candidate',passed=bool(max(derivative)<.002 and max(max(x['constitutive'],x['spectral']) for x in checks)<.002),
        cells=n,new_native_cells=4,retained15_initial_banks_bitwise=True,nearest_center_face_gap_m=float((actual_edges[-1]-r[-1])/100),
        native_derivative_relative=max(derivative),controls=checks,native_calls=native.ion.calls,seconds=time.monotonic()-start)
    write(OUT/'bank.json',row);signal.alarm(0);print(json.dumps(row),flush=True);assert row['passed']


def pilot():
    assert json.loads((OUT/'bank.json').read_text())['passed'];assert not (OUT/'pilot.json').exists();start=time.monotonic();signal.signal(signal.SIGALRM,old.optical.timeout);signal.alarm(35);rows=[]
    model=Coupled();assert model.bulk.n==19 and model.bulk.d['edges'][-1]==model.m.rf[0]
    assert abs(model.bulk.area[-1]/model.area[0]-1)<1e-14
    zero=np.zeros(19);p=model.bulk.eos.gas(zero,zero)[0];assert np.max(abs(p/model.f0['p0']-1))<1e-12
    for steps in [64,128]:
        row=(model if steps==64 else Coupled()).run(steps,f'pilot-{steps}',2);rows.append(row)
        if not row['passed']:break
    forecast=None if len(rows)!=2 or not all(r['passed'] for r in rows) else sum(r['seconds']/2*(r['steps']-2) for r in rows)
    upper=None if forecast is None else 1.8*forecast+10
    write(OUT/'pilot.json',dict(classification='Counterexample candidate',paths=rows,forecast_seconds=forecast,upper_seconds=upper,
        eligible=bool(upper is not None and upper<500),seconds=time.monotonic()-start));signal.alarm(0)


def production():
    assert json.loads((OUT/'pilot.json').read_text())['eligible'];assert not (OUT/'production.json').exists()
    for path,value in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert sha(path)==value,path
    start=time.monotonic();signal.signal(signal.SIGALRM,old.optical.timeout);signal.alarm(500);rows=[]
    for steps in [64,128]:
        row=Coupled().run(steps,f'coupled-{steps}',restart=f'pilot-{steps}');rows.append(row)
        if not row['passed']:break
    write(OUT/'production.json',dict(classification='Counterexample candidate',passed=bool(len(rows)==2 and all(r['passed'] for r in rows)),paths=rows,seconds=time.monotonic()-start,
        actual_local_native_refinement=True,global_radial_error_certified=False,full_GR_feedback=False,final_charge_solved=False));signal.alarm(0)


if __name__=='__main__':globals()[sys.argv[1]]()
