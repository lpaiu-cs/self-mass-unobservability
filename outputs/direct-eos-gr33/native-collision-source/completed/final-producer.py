"""Counterexample candidate: integrate the known paired collision forcing.

Keep the actual Radau stages and physical angular ports. Only the known
collision of the exact photon source primitive gets moment-fitted rates.
"""
from pathlib import Path
import inspect,json,math,resource,shutil,sys,time
import numpy as np
import sympy as sp
import solve_native_radau_transfer as radau

OUT=Path('native-collision-source169-work');OLD=radau.OUT;original=radau.prior
read,write,sha=radau.read,radau.write,radau.sha
AMP=radau.AMP;LD=np.longdouble;C=original.C
CAPS=dict(prepare=30,inspect=60,check=90,pilot=240,finish=30);TOTAL=450
FITS=[];raw_stages=radau.stages


def wave_moment(d,lo,hi,x,m):
    """Integral (time-lo)^m wave(time,x), primary plus declared linear Born."""
    x=np.asarray(x,LD);D=LD(d.D);dt=LD(hi-lo);a=(LD(lo)+x/C)/D;b=a+dt/D
    coeff=np.r_[np.zeros(4,dtype=LD),LD(256)*np.array([1,-4,6,-4,1],LD)]
    def integral(power):
        value=np.zeros_like(x)
        for k in range(power+1):
            poly=np.polynomial.polynomial.polyint(np.r_[np.zeros(k,dtype=LD),coeff])
            value+=math.comb(power,k)*(-a)**(power-k)*(np.polynomial.polynomial.polyval(np.clip(b,0,1),poly)-np.polynomial.polynomial.polyval(np.clip(a,0,1),poly))
        return D**(power+1)*value
    pulse=lambda v:np.where((v>0)&(v<1),np.polynomial.polynomial.polyval(np.clip(v,0,1),coeff),LD(0))
    U=integral(m);Ut=dt**m*pulse(b)-(pulse(a) if m==0 else m*integral(m-1))
    result=[LD(original.inf.incident.ETA)*LD(d.r0)*v for v in [U,Ut,Ut/C]]
    for i,key in enumerate(['U','Ut','Ux']):
        value=np.zeros_like(d.born[key][0],dtype=LD)
        for j in range(len(d.times)-1):
            left=max(lo,d.times[j]);right=min(hi,d.times[j+1])
            if right<=left:continue
            slope=(d.born[key][j+1].astype(LD)-d.born[key][j])/(d.times[j+1]-d.times[j])
            y=d.born[key][j].astype(LD)+(left-d.times[j])*slope;h=LD(right-left);w=LD(left-lo)
            for k in range(m+1):value+=math.comb(m,k)*w**(m-k)*(y*h**(k+1)/(k+1)+slope*h**(k+2)/(k+2))
        result[i]+=np.interp(np.asarray(x,float),d.tx,np.asarray(value,float))
    return result


def lift_moments(model,lo,hi):
    """Integrals of H and (time-lo)H, not a newly sampled forcing history."""
    j=np.clip(np.searchsorted(model.t,lo,side='right')-1,0,len(model.t)-2)
    assert hi<=model.t[j+1]+1e-18
    dt=LD(hi-lo);width=model.t[j+1]-model.t[j];f=(lo-model.t[j])/width
    I=(1-f)*model.I[j]+f*model.I[j+1];It=(model.I[j+1]-model.I[j])/width
    H=model.lift(lo)[0];d=model.redshift_driver;wave=d.wave;fields=[]
    try:
        for power in range(4):
            d.wave=lambda now,x:wave_moment(d,lo,hi,x,power)
            d.cache={};fields.append(d.at(hi))
    finally:d.wave=wave;d.cache={}
    M=[]
    for k in range(3):
        a=original.frequency_source(model,I,fields[k],model.drive_scale)
        b=original.frequency_source(model,It,fields[k+1],model.drive_scale)
        assert max(a[2],b[2])<1e-12
        M.append((a[0].astype(LD)+b[0])/(model.scale*AMP))
    return dt*H+dt*M[0]-M[1],dt*dt/2*H+(dt*dt*M[0]-M[2])/2


def collision_fit(model,t,h,cs):
    assert np.count_nonzero(model.motion)==0 and np.count_nonzero(model.energy_offset)==0
    integrals=[];zero=np.zeros((model.n,2));errors=[]
    for fraction in radau.RK_C:
        K0,K1=lift_moments(model,t,t+h*fraction)
        # C(time) is affine on this existing canonical background interval.
        # The two actual Radau coefficient points reconstruct it exactly.
        a=model.collision(cs[0],LD(1.5)*(K0-K1/h),zero)
        b=model.collision(cs[1],LD(1.5)*K1/h-LD(.5)*K0,zero)
        integrals.append([a[k]+b[k] for k in [0,3,2]])
    inverse=np.linalg.inv(radau.RK_A)/h
    for j in range(2):
        H=model.lift(t+h*radau.RK_C[j])[0];p,_,e,b=model.collision(cs[j],H,zero)
        errors.append(max(float(np.max(abs(cs[j][key]-target))/max(np.max(abs(target)),1e-290))
                          for key,target in [('q',p),('qb',b),('qe',e)]))
    assert max(errors)<1e-12,('Known collision source is not exactly C(t)H(t)',errors)
    for j in range(2):
        for k,key in enumerate(['q','qb','qe']):cs[j][key]=sum(inverse[j,l]*integrals[l][k] for l in range(2))
    moment=0.
    for i in range(2):
        for k,key in enumerate(['q','qb','qe']):
            got=h*sum(radau.RK_A[i,j]*cs[j][key] for j in range(2));target=integrals[i][k]
            moment=max(moment,float(np.max(abs(got-target))/max(np.max(abs(target)),1e-290)))
    assert moment<1e-12,moment
    FITS.append(dict(t=float(t),h=float(h),known_source_relative=max(errors),fitted_moment_relative=moment))
    return cs


def initialize(fitted=True):
    radau.OUT=OUT
    if fitted:
        source=inspect.getsource(raw_stages);anchor="    mechanical=cs[0]['mechanical'].copy()"
        assert source.count(anchor)==1;source=source.replace(anchor,"    cs=collision_fit(m,t,h,cs)\n"+anchor)
        ns=dict(raw_stages.__globals__,collision_fit=collision_fit);exec(compile(source,__file__,'exec'),ns)
        radau.stages=ns['stages'];(OUT/'fitted-stage-producer.py').write_text(source)
    base_initialize()


def prepare():
    assert not OUT.exists() and not read(OLD/'result.json')['passed'];OUT.mkdir()
    files=[Path(__file__),Path(radau.__file__),Path(original.__file__),OLD/'result.json',OLD/'error-profile.json']
    for s in [0,1]:
        for folder in ['photons','material']:(OUT/f'sweep-{s}/{folder}').mkdir(parents=True)
    for p in (OLD/'sweep-0').rglob('*.npz'):
        shutil.copyfile(p,OUT/p.relative_to(OLD));files.append(p)
    for name in ['normalization.json','photon-conservation-plan.json','lift-plan.json','primitive-plan.json']:
        shutil.copyfile(OLD/name,OUT/name);files.append(OLD/name)
    a=sp.Matrix([[sp.Rational(5,12),-sp.Rational(1,12)],[sp.Rational(3,4),sp.Rational(1,4)]])
    h,K0,K1=sp.symbols('h K0 K1',nonzero=True);c0,c1=sp.symbols('c0 c1')
    assert sp.expand(c0*sp.Rational(3,2)*(K0-K1/h)+c1*(sp.Rational(3,2)*K1/h-K0/2)-((3*c0-c1)*K0/2+3*(c1-c0)*K1/(2*h)))==0
    J=sp.Matrix(sp.symbols('J1 J2'));assert h*a*(a.inv()*J/h)==J
    z=sp.symbols('z')
    for power in range(11):
        M=[h**(power+j+1)/sp.Integer(power+j+1) for j in range(3)]
        H=z**(power+1)/sp.Integer(power+1)
        assert sp.simplify(h*M[0]-M[1]-sp.integrate(H,(z,0,h)))==0
        assert sp.simplify((h*h*M[0]-M[2])/2-sp.integrate(z*H,(z,0,h)))==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        identities=['integral H = dt*H(lo)+dt*M0-M1','integral (time-lo)*H = dt^2*H(lo)/2+(dt^2*M0-M2)/2',
                    'C(t) affine: integral C*H=C1*1.5*(K0-K1/h)+C2*(1.5*K1/h-K0/2)',
                    'q_effective=A_inverse*J/h gives h*A*q_effective=J; paired collision moments and spectral ghosts use the same rates'],
        scope='Known-source quadrature identities in the declared finite linear model; no global PDE or EOS error bound.'))
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='39341e99f',
        claim='Resolve the cell15-dominated material time mismatch by integrating the known H collision source over each actual Radau stage interval, then applying the same paired source to both photon and gas equations.',
        decision='User criterion: scientific value is whether the final-charge conclusion survives on the SAME coupled solution after its dominant error is resolved. Source/prefix gates are intermediate prerequisites only. A passed unchanged4/8-step pair supports continuing the corrected coupled solution and final-charge readout with measured budgets; it does not establish that conclusion. If it fails, preserve it and do not expand clocks, horizon or schemes automatically.',
        method='Exact moments0..3 of compact polynomial and separately stored linear Born U/Ut/Ux; combine with original linear photon history and affine collision coefficients. Fit q/qb/qe to the two Radau cumulative integrals. Unknown stages, operator, physical angular samples, stream forcing, gates and physical inputs stay unchanged.',
        restrictions='This increment has zero prior motion and extra metric. Assert its known collision source equals C(t)H(t). Do not silently reuse this specialization for a later nonzero mechanical sweep.',
        gates=dict(time=.02,energy_H=1e-8,stage=1e-12,physical_stage_moment=1e-13,port=1e-12,source=1e-12,moment_derivative=1e-5),
        budget=dict(actions=CAPS,total_action_seconds=TOTAL,CPU_threads=1,virtual_GiB=3),
        forecast='The actual168pair used84.84action seconds and2.412GB peakRSS. Added analytic source moments can increase work and live arrays;240s is a bounded pilot cap, not measured production speed. Reuse both stage coefficient objects and release each clock model. No full-horizon computation admitted.',
        inputs='Same531cells8angles152frequencies64/128 clocks,17background knots,primary plus finite Born and original0.214652ms prefix. No native roots or new background.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)}))


def inspect_input():
    initialize(False);m=original.c.Response(128);d=m.redshift_driver;cell=15;end=m.t[-1]/16
    rows=[]
    for k in [14,15,16,17,18]:rows.append(dict(cell=k,radius_cm=float(m.r[k]),scalar_arrival_seconds=float(-d.xc[k]/C),
        arrival_over_prefix=float(-d.xc[k]/C/end)))
    c=m.local(end/2);rate=c['loss'][cell];x=m.lift(end)[0]
    p,g,e,b=m.collision(c,x,np.zeros((m.n,2)))
    result=dict(classification='Counterexample candidate',rows=rows,prefix_seconds=float(end),deep_cells=m.nb,
        collision_absorption_range=[float(rate.min()),float(rate.max())],
        last_deep_cell=bool(cell==m.nb-1),source_is_paired_collision=True,
        limits='Arrival and saved source localization are inputs to the repair, not an isolated proof of the unique error cause.')
    write(OUT/'inspection.json',result);print(json.dumps(result),flush=True)


def check():
    initialize(False);m=original.c.Response(128);d=m.redshift_driver;h=m.t[-1]/64;rows=[]
    for lo,hi in [(0.,h),(2*h,3*h),(3*h,4*h)]:
        for power in [0,1]:
            new=wave_moment(d,lo,hi,d.xc,power);old=original.integrated_wave(d,lo,hi,d.xc,power)
            error=max(float(np.max(abs(a-b))/max(np.max(abs(b)),1e-290)) for a,b in zip(new,old))
            assert error<1e-12,error
        end=hi-h*.1;eps=h*1e-5
        a=lift_moments(m,lo,end+eps);b=lift_moments(m,lo,end-eps);H=m.lift(end)[0]
        errors=[float(np.max(abs((a[j]-b[j])/(2*eps)-H*(end-lo)**j))/max(np.max(abs(H*(end-lo)**j)),1e-290)) for j in range(2)]
        cs=[m.local(lo+h*f) for f in radau.RK_C];collision_fit(m,lo,h,cs)
        rows.append(dict(lo=float(lo),hi=float(hi),derivative_time=float(end),moment_derivative_relative=errors))
        assert max(errors)<1e-5,rows[-1]
    write(OUT/'source-check.json',dict(classification='Counterexample candidate',passed=True,rows=rows,fits=FITS))
    print(json.dumps(rows),flush=True)


def pilot():
    assert read(OUT/'source-check.json')['passed']
    # Avoid recursion: pilot's initializer must invoke the captured Radau owner.
    radau.initialize=lambda:initialize()
    radau.pilot()
    write(OUT/'actual-source-fits.json',dict(classification='Counterexample candidate',rows=FITS,passed=True))


base_initialize=radau.initialize


if __name__=='__main__':
    action=sys.argv[1];receipt=OUT/f'{action}-receipt.json';assert action in CAPS and not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
            assert sum(read(p)['seconds'] for p in OUT.glob('*-receipt.json'))+CAPS[action]<=TOTAL
        if action=='inspect':inspect_input()
        elif action=='finish':radau.OUT=OUT;radau.finish()
        else:globals()[action]()
    except Exception as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(action=action,seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
