"""Counterexample candidate: integrate known collision AND stream forcing.

Keep physical stage samples. Store their separate known-source quadrature
correction; a corrected quadrature rate is not a pointwise physical packet.
"""
from pathlib import Path
import gc,inspect,json,resource,shutil,sys,time
import numpy as np
import sympy as sp
import integrate_native_collision_source as prior

OUT=Path('native-known-forcing170-work');OLD=prior.OUT
radau=prior.radau;original=prior.original
read,write,sha=prior.read,prior.write,prior.sha;LD=prior.LD;AMP=prior.AMP
CAPS=dict(prepare=30,check=90,pilot=240,finish=30);TOTAL=390
FITS=[]


def fit(m,t,h,cs,sources,ledgers):
    assert np.count_nonzero(m.motion)==0 and np.count_nonzero(m.energy_offset)==0
    inverse=np.linalg.inv(radau.RK_A)/h;integrals=[];HK=[]
    stream=lambda x:(m.A@x.reshape(m.n*m.q,m.nf)).reshape(x.shape)
    for fraction in radau.RK_C:
        K0,K1=prior.lift_moments(m,t,t+h*fraction);HK.append(K0)
        a=m.collision(cs[0],LD(1.5)*(K0-K1/h),np.zeros((m.n,2)))
        b=m.collision(cs[1],LD(1.5)*K1/h-LD(.5)*K0,np.zeros((m.n,2)))
        integrals.append([a[k]+b[k] for k in [0,3,2]])
    identity=0.;m.stage_lift_delta={};known_stream=[]
    for j,fraction in enumerate(radau.RK_C):
        now=t+h*fraction;H=m.lift(now)[0];p,_,e,b=m.collision(cs[j],H,np.zeros((m.n,2)))
        for key,target in [('q',p),('qb',b),('qe',e)]:
            identity=max(identity,float(np.max(abs(cs[j][key]-target))/max(np.max(abs(target)),1e-290)))
        base=stream(H);identity=max(identity,float(np.max(abs(sources[j]-base))/max(np.max(abs(base)),1e-290)))
        effective=sum(inverse[j,i]*HK[i] for i in range(2));delta=effective-H
        m.stage_lift_delta[now]=delta
        known_stream.append(stream(effective));sources[j]+=stream(delta)
        ledgers[j][3:]+=m.port(delta*m.scale)
        for k,key in enumerate(['q','qb','qe']):cs[j][key]=sum(inverse[j,i]*integrals[i][k] for i in range(2))
    assert identity<1e-10,('Known forcing specialization',identity)
    moment=0.;balance=0.
    for i in range(2):
        for k,key in enumerate(['q','qb','qe']):
            got=h*sum(radau.RK_A[i,j]*cs[j][key] for j in range(2));target=integrals[i][k]
            moment=max(moment,float(np.max(abs(got-target))/max(np.max(abs(target)),1e-290)))
        got=h*sum(radau.RK_A[i,j]*known_stream[j] for j in range(2));target=stream(HK[i])
        moment=max(moment,float(np.max(abs(got-target))/max(np.max(abs(target)),1e-290)))
        packet=target*m.scale;norm=np.array([np.sum(abs(packet)*m.weights,dtype=LD),np.sum(abs(packet)*m.weights*m.E,dtype=LD)])
        balance=max(balance,float(np.max(abs(m.moments(packet)-m.port(HK[i]*m.scale))/np.maximum(norm,LD('1e-290')))))
    assert max(moment,balance)<1e-12,(moment,balance)
    FITS.append(dict(t=float(t),h=float(h),known_source_relative=identity,fitted_moment_relative=moment,stream_port_balance_relative=balance))
    return cs,sources,ledgers


def initialize():
    prior.OUT=OUT;radau.OUT=OUT
    source=inspect.getsource(prior.raw_stages);anchor="    mechanical=cs[0]['mechanical'].copy()"
    assert source.count(anchor)==1;source=source.replace(anchor,"    cs,sources,ls=fit(m,t,h,cs,sources,ls)\n"+anchor)
    namespace=dict(prior.raw_stages.__globals__,fit=fit);exec(compile(source,__file__,'exec'),namespace)
    radau.stages=namespace['stages'];prior.base_initialize()
    (OUT/'fitted-stage-producer.py').write_text(source)
    Parent=original.c.Response
    class Response(Parent):
        def __init__(self,n):super().__init__(n);self.angular_correction=[]
        def boundary_ports(self,t,x):
            actual=super().boundary_ports(t,x);sample=np.asarray(self.angular[-1]).copy()
            effective=super().boundary_ports(t,x+self.stage_lift_delta[t])
            correction=np.asarray(self.angular.pop())-sample;assert self.angular_times.pop()==t
            self.angular_correction.append(correction)
            return effective
    runner=(OUT/'sweep-1/expanded-radau-run.py').read_text()
    runner=runner.replace('accepted_angular_luminosity=self.angular,',
        'accepted_angular_luminosity=self.angular,accepted_angular_quadrature_correction=self.angular_correction,')
    # Checkpoint and final output both retain actual samples and the correction.
    runner=runner.replace('accepted_angular_luminosity=self.angular)',
        'accepted_angular_luminosity=self.angular,accepted_angular_quadrature_correction=self.angular_correction)')
    anchor="self.angular=list(z['accepted_angular_luminosity']);"
    assert runner.count(anchor)==1;runner=runner.replace(anchor,anchor+"self.angular_correction=list(z['accepted_angular_quadrature_correction']);")
    namespace=dict(Parent.run.__globals__);exec(compile(runner,__file__,'exec'),namespace)
    Response.run=namespace['run'];original.c.Response=Response
    (OUT/'sweep-1/expanded-known-forcing-run.py').write_text(runner)


def prepare():
    assert not OUT.exists() and not read(OLD/'result.json')['passed'];OUT.mkdir()
    files=[Path(__file__),Path(prior.__file__),Path(radau.__file__),OLD/'result.json',OLD/'pressure-profile.json',OLD/'source-check.json']
    for s in [0,1]:
        for folder in ['photons','material']:(OUT/f'sweep-{s}/{folder}').mkdir(parents=True)
    for p in (OLD/'sweep-0').rglob('*.npz'):
        shutil.copyfile(p,OUT/p.relative_to(OLD));files.append(p)
    for name in ['normalization.json','photon-conservation-plan.json','lift-plan.json','primitive-plan.json']:
        shutil.copyfile(OLD/name,OUT/name);files.append(OLD/name)
    a=sp.Matrix([[sp.Rational(5,12),-sp.Rational(1,12)],[sp.Rational(3,4),sp.Rational(1,4)]])
    h,L,P=sp.symbols('h L P',nonzero=True);J=sp.Matrix(sp.symbols('J1 J2'));fitted=a.inv()*J/h
    assert h*a*(L*fitted)==L*J and h*a*(P*fitted)==P*J
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        identity='For known H and fixed linear stream/port maps, Heff=A_inverse*integral(H)/h integrates both maps consistently at Radau stage endpoints. Collision uses the separate exact integral C(t)H. Actual point values H(tj) and quadrature rates Heff_j remain distinct.',
        scope='Known-forcing quadrature identity only, not stiff/global error certification.'))
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='873a09ad13d0dd9df748d8ab6ea871e119f21068',
        claim='Test whether consistent integration of all known H forcing resolves the pressure mismatch in the same actual coupled photon/gas prefix.',
        cause_evidence='169pressure is the difference of large E/H contributions. Collision-only fitting leaves streamH and its boundary contribution point-sampled; the combined forcing quadrature is inconsistent in accuracy. This is a targeted candidate, not a demonstrated unique cause.',
        method='Reuse169analytic moments. Fit streamH and both boundary port moments to integralH together with C(t)H. Keep unknown coupled Radau stages and actual angular samples; save the separate known-source angular quadrature correction. No correction is added to a mass or final charge.',
        decision='Only a pair satisfying all unchanged gates permits considering full coupled continuation and final-charge readout. User criterion remains final-charge stability on the SAME corrected coupled solution; intermediate gates do not satisfy it.',
        gates=dict(time=.02,energy_H=1e-8,stage=1e-12,physical_stage_moment=1e-13,port=1e-12,source=1e-10,integral=1e-12),
        restrictions='Same zero previous motion/extra metric increment,531cells8angles152frequencies,64/128 clocks,4/8steps,0.214652ms. Future nonzero mechanical sweeps must generalize this asserted specialization.',
        consumers='Physical angular samples alone no longer integrate radial_ports. Consume saved weights times (actual sample plus known-source quadrature correction), and use the same representation in retarded kernels with its original time control.',
        budget=dict(actions=CAPS,total_action_seconds=TOTAL,CPU_threads=1,virtual_GiB=3),
        forecast='169actual pair128.03s; reuse the same analytic moments with a few added sparse stream/port operations. Estimated100-170s subject to contention;240s hard cap. No production or extra pair admitted.',
        stop='Any gate or cap ends this trial. Preserve169failure; no automatic finer clock, longer interval, method ladder, relaxed pressure gate or full production.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)}))


def check():
    initialize();m=original.c.Response(128);h=m.t[-1]/64
    for t in [0.,2*h,3*h]:
        cs=[];ss=[];ls=[]
        for c in radau.RK_C:
            now=t+h*c;cs.append(m.local(now));q,l,_=m.source(now);ss.append(q/(m.scale*AMP));ls.append(l/AMP)
        fit(m,t,h,cs,ss,ls)
        m.angular=[];m.angular_times=[];m.angular_correction=[]
        for c in radau.RK_C:
            now=t+h*c;port=m.boundary_ports(now,np.zeros_like(m.I[0]))
            actual=np.asarray(m.angular[-1]);correction=np.asarray(m.angular_correction[-1]);flux=(actual+correction)@(np.arange(1,8,2)/32)
            assert abs(flux-port[1,1]*AMP)/max(abs(flux),1e-290)<1e-12
    write(OUT/'source-check.json',dict(classification='Counterexample candidate',passed=True,rows=FITS))
    print(json.dumps(FITS),flush=True)


def pilot():
    assert read(OUT/'source-check.json')['passed'];initialize()
    source=inspect.getsource(radau.pilot).replace('    initialize();','    pass;')
    source=source.replace("p['accepted_angular_luminosity']@", "(p['accepted_angular_luminosity']+p['accepted_angular_quadrature_correction'])@")
    source=source.replace("row['actual_angular_quadrature_relative']=error;", "row['actual_angular_quadrature_relative']=error;row['angular_quadrature_includes_separate_known_correction']=True;")
    namespace=dict(radau.pilot.__globals__);exec(compile(source,__file__,'exec'),namespace)
    (OUT/'pilot-producer.py').write_text(source)
    try:namespace['pilot']()
    finally:write(OUT/'actual-source-fits.json',dict(classification='Counterexample candidate',rows=FITS))


def finish():
    radau.OUT=OUT;radau.finish()
    result=read(OUT/'result.json');result.pop('export_repair',None)
    result['final_charge_conclusion']='unadjudicated on the corrected coupled solution'
    result['angular_representation']='Actual physical samples plus a separately stored known-source quadrature correction; raw point samples are preserved.'
    write(OUT/'result.json',result)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
            assert sum(read(p)['seconds'] for p in OUT.glob('*-receipt.json'))+CAPS[action]<=TOTAL
        globals()[action]()
    except Exception as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(action=action,seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
