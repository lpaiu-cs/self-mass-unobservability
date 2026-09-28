"""Counterexample candidate: close native neutral transport inside radiation.

Reuse the frozen175 histories and native tangent. Before any new evolution,
identify the actual missing neutral rate and test a sparse stage representation.
"""
from pathlib import Path
import inspect,json,resource,shutil,sys,time
import numpy as np
import sympy as sp
from scipy import sparse
import return_native_pressure_reciprocity as previous

OUT=Path('native-neutral-coupled177-work');BEFORE=previous.OUT
read,write,sha=previous.read,previous.write,previous.sha
LD=np.longdouble;AMP=previous.AMP
CAPS=dict(prepare=10,localize=30,assemble=60,pilot=240)
BASE_PATHS=previous.paths
radau=previous.front.radau


def paths(sweep):return OUT/f'sweep-{sweep}/photons',OUT/f'sweep-{sweep}/material'


def flux(m,k,z):return previous.engine.tangent(m,k,z,1.)[0][3]


def face_matrix(m,k,component):
    """Five colors separate every column of the four-cell face stencil."""
    q=m.point(k)['Q'];units=np.maximum(abs(q[component]),1.).astype(LD)
    rows=[];cols=[];values=[];odd=0.
    for color in range(5):
        z=np.zeros((4,m.n),LD);z[component,color::5]=units[color::5]
        a,b=flux(m,k,z),flux(m,k,-z)
        odd=max(odd,float(np.sum(abs(a+b),dtype=LD)/max(np.sum(abs(a),dtype=LD),LD('1e-290'))))
        for offset in [-2,-1,0,1]:
            ii=np.arange(m.n+1);jj=ii+offset;mask=(jj>=0)&(jj<m.n)&(jj%5==color)
            ii,jj=ii[mask],jj[mask];rows.extend(ii);cols.extend(jj);values.extend(a[ii]/units[jj])
    matrix=sparse.coo_matrix((np.asarray(values,float),(rows,cols)),shape=(m.n+1,m.n)).tocsr()
    matrix.eliminate_zeros();return matrix,odd


def assemble():
    assert read(OUT/'localize-receipt.json')['error'] is None
    write(OUT/'stage-plan.json',dict(classification='Conjectural',
        evidence='At saved knots1,8,16 the native neutral-rate change is dominated by the H state, with B/S/E changes far smaller; native component additivity is within1.7e-15 on those directions.',
        repair='Move T_H(B,S,Eref,H) into the actual Radau photon/thermal/H stage. Assemble native face-flux E/H columns with five nonoverlapping stencil colors; retain lagged B/S, and keep the original Etilde mechanical source. No post-hoc replacement of evolved H.',
        controls='Both signs of every colored column, all17 saved old/new states and independent signed oscillatory E/H directions must reproduce native face flux within0.2percent. Form cell divergence from one shared face matrix. Check the native flux again on EVERY accepted new Radau stage.',
        integration='Same RadauIIA2,64/128 base clocks and fixed local bisections. H affine source uses its own actual stage time; only the retained Etilde mechanical rate stays first-stage-held. Store actual H transport stages and use their Radau-weighted value in collision/number ledgers.',
        budgets=CAPS,CPU_threads=1,virtual_GiB=3,full_horizon_authorized=False,
        admission='Assemble within60s; then one equal-horizon4/8macro-step pilot within240s. Original2percent time,1e-8 balances,1e-12 stage and0.2percent native-H gates. Preserve175 failure. No extra full sweep or clock refinement on failure.',
        forecast='175 pilots took135.1s including startup; native material constructor about4.5s each. Sparse transport adds work not yet measured;240s is a hard pilot cap. Full-period cost must be measured separately before admission.',
        bindings={str(p):sha(p) for p in [Path(__file__),OUT/'plan.json',OUT/'localization.json',
            Path(radau.__file__),Path(previous.engine.face.feedback.__file__)]}))
    previous.initialize();m=previous.run.c.Material(128,128);states=[];rows=[]
    for s in [0,1]:
        d=np.load(BASE_PATHS(s)[1]/'steps-128-reference-128.npz')
        states.append(d['history_scaled'][[np.argmin(abs(d['t']-t)) for t in m.t]].copy())
    for k in range(len(m.t)):
        matrices=[];odd=[]
        for component in [2,3]:
            matrix,error=face_matrix(m,k,component);matrices.append(matrix);odd.append(error)
            sparse.save_npz(OUT/f'face-{k}-{component}.npz',matrix)
        errors=[]
        for state in [states[0][k],states[1][k]]:
            fixed=state.copy();fixed[2:]=0
            exact=flux(m,k,state);estimate=flux(m,k,fixed)+matrices[0]@state[2]+matrices[1]@state[3]
            errors.append(float(np.sum(abs(estimate-exact),dtype=LD)/max(np.sum(abs(exact),dtype=LD),LD('1e-290'))))
        for sign in [-1,1]:
            state=np.zeros((4,m.n),LD);q=m.point(k)['Q'];phase=np.arange(m.n)+1.
            state[2]=sign*np.maximum(abs(q[2]),1.)*np.sin(phase)
            state[3]=sign*np.maximum(abs(q[3]),1.)*np.cos(phase*.71)
            exact=flux(m,k,state);estimate=matrices[0]@state[2]+matrices[1]@state[3]
            errors.append(float(np.sum(abs(estimate-exact),dtype=LD)/max(np.sum(abs(exact),dtype=LD),LD('1e-290'))))
        rows.append(dict(k=k,odd_relative=odd,native_face_relative=errors,nnz=[a.nnz for a in matrices]))
        write(OUT/'assembly-result.json',dict(classification='Counterexample candidate',passed=False,completed_knots=k+1,rows=rows))
        assert max(odd+errors)<.002,rows[-1]
    files=list(OUT.glob('face-*.npz'))
    write(OUT/'assembly-result.json',dict(classification='Counterexample candidate',passed=True,rows=rows,
        matrices={str(p):sha(p) for p in files},scope='Native face representation on declared controls; actual coupled-stage verification remains mandatory.'))
    for s in [0,1]:
        for p in paths(s):p.mkdir(parents=True)
    for old,new in zip(BASE_PATHS(1),paths(0)):
        for n in [64,128]:
            p=old/f'steps-{n}-reference-128.npz';shutil.copyfile(p,new/p.name);assert sha(p)==sha(new/p.name)
    for name in ['normalization.json','photon-conservation-plan.json']:
        shutil.copyfile(BEFORE/name,OUT/name)


def stage_source():
    source=inspect.getsource(radau.stages)
    def change(old,new):
        nonlocal source
        assert source.count(old)==1,old;source=source.replace(old,new)
    change("mechanical=cs[0]['mechanical'].copy()\n    for c in cs:c['mechanical']=mechanical",
        "for c in cs:c['mechanical'][:,0]=cs[0]['mechanical'][:,0]")
    change("m.gas(c['q'],c['qb'],c['qe'])+mechanical", "m.gas(c['q'],c['qb'],c['qe'])+c['mechanical']")
    change('    result=[]\n', '    result=[];transport=[]\n')
    change("xx,gg=m.unpack(vv);p,q,e,_=m.collision(cs[j],xx,gg,True)",
        "xx,gg=m.unpack(vv);p,q,e,_=m.collision(cs[j],xx,gg,True)\n        mech=np.asarray(cs[j]['mechanical'],np.longdouble).copy();mech[:,1]+=m.neutral_linear(cs[j],gg)\n        m.check_neutral(cs[j],gg,mech[:,1],h*RK_B[j]);transport.append(mech)")
    change('return result,mechanical', 'return result,sum(b*v for b,v in zip(RK_B,transport))')
    return source


def initialize_coupled():
    global Response
    previous.OUT=OUT;previous.paths=paths;previous.initialize();Parent=previous.run.c.Response
    class Coupled(Parent):
        def __init__(self,n):
            super().__init__(n);self.neutral_material=self.material
            self.material.face_owner_error=0.
            self.material.deep_tangent=previous.original.reuse.deep.deep_tangent.__get__(self.material,type(self.material))
            self.neutral_faces={k:[sparse.load_npz(OUT/f'face-{k}-{j}.npz') for j in [2,3]] for k in range(17)}
            self.neutral_times=[];self.neutral_weights=[];self.neutral_rates=[];self.neutral_error=0.
        def neutral_fixed(self,t):
            k=np.clip(np.searchsorted(self.t,t,side='right')-1,0,15);w=(t-self.t[k])/(self.t[k+1]-self.t[k])
            z=((1-w)*self.motion[k]+w*self.motion[k+1]).astype(LD);z[2:]=0
            z[2]=(1-w)*self.energy_offset[k]+w*self.energy_offset[k+1]
            return k,w,z
        def local(self,t):
            c=super().local(t);k,w,z=self.neutral_fixed(t)
            c['neutral_faces']=[(1-w)*self.neutral_faces[k][j]+w*self.neutral_faces[k+1][j] for j in [0,1]]
            f=sum(weight*flux(self.neutral_material,i,z) for i,weight in [(k,1-w),(k+1,w)] if weight)
            c['mechanical'][:,1]=-np.diff(f)/self.nu;c['neutral_time']=t
            return c
        def neutral_linear(self,c,g):
            f=c['neutral_faces'][0]@(g[:,0]*self.eu)+c['neutral_faces'][1]@(g[:,1]*self.nu)
            return -np.diff(f)/self.nu
        def collision(self,c,x,g,source=False):
            p,q,e,b=super().collision(c,x,g,source);q[:,1]+=self.neutral_linear(c,g)
            return p,q,e,b
        def check_neutral(self,c,g,rate,weight):
            t=c['neutral_time'];k,w,z=self.neutral_fixed(t);z[2]+=g[:,0]*self.eu;z[3]=g[:,1]*self.nu
            native=sum(v*flux(self.neutral_material,i,z) for i,v in [(k,1-w),(k+1,w)] if v)
            exact=-np.diff(native);actual=rate*self.nu
            error=float(np.sum(abs(actual-exact),dtype=LD)/max(np.sum(abs(exact),dtype=LD),LD('1e-290')))
            self.neutral_error=max(self.neutral_error,error);assert error<.002,('Actual native H stage',t,error)
            self.neutral_times.append(t);self.neutral_weights.append(weight);self.neutral_rates.append(actual.copy()*AMP)
    Response=Coupled
    ns=dict(radau.stages.__globals__);exec(compile(stage_source(),__file__,'exec'),ns)
    source=(OUT/'sweep-1/expanded-corrected-run.py').read_text()
    source=source.replace('transfer=np.zeros((self.n,2));','transfer=np.zeros((self.n,2),dtype=LD);')
    namespace=dict(Parent.run.__globals__,stages=ns['stages'],LD=LD)
    exec(compile(source,__file__,'exec'),namespace);Coupled.run=namespace['run']
    (OUT/'expanded-neutral-stages.py').write_text(stage_source());(OUT/'expanded-neutral-run.py').write_text(source)


def pilot():
    assembled=read(OUT/'assembly-result.json');assert assembled['passed']
    for p,h in assembled['matrices'].items():assert sha(p)==h,p
    initialize_coupled();rows=[];histories=[];start=time.monotonic()
    for n in [64,128]:
        mark=time.monotonic();m=Response(n);label=f'pilot-{n}'
        row=m.run(n,label,n//16);path=paths(1)[0]/f'{label}.npz';data=dict(np.load(path))
        times=np.array(m.neutral_times);weights=np.array(m.neutral_weights);rates=np.array(m.neutral_rates)
        assert np.max(abs(times-data['accepted_angular_times']))<1e-18
        assert np.max(abs(weights-data['accepted_angular_quadrature_weights']))<1e-18
        data.update(native_neutral_stage_times=times,native_neutral_stage_weights=weights,native_neutral_stage_rates=rates)
        integral=np.sum(weights[:,None].astype(LD)*rates,axis=0,dtype=LD)
        direct=data['moments'][-1,2].astype(LD)-data['collision_transfer'][-1,:,1].astype(LD)
        error=float(np.sum(abs(direct-integral),dtype=LD)/max(np.sum(abs(integral),dtype=LD),LD('1e-290')))
        row.update(native_neutral_stage_relative=m.neutral_error,neutral_partition_relative=error,
            native_H_transport_in_Radau_operator=True,lagged_H_mechanical_source_used=False,worker_seconds=time.monotonic()-mark)
        row['passed']=row['passed'] and m.neutral_error<.002 and error<.002
        np.savez_compressed(path,**data);write(path.with_suffix('.json'),row);rows.append(row)
        histories.append(data['moments'][:,[0,1,2,3,5,6]]);assert row['passed'],row
    norm=np.maximum(np.max(np.sum(abs(histories[1]),axis=2),axis=0),LD('1e-290'))
    errors=(np.max(np.sum(abs(histories[0]-histories[1]),axis=2),axis=0)/norm).astype(float).tolist()
    result=dict(classification='Counterexample candidate',passed=max(errors)<.02,rows=rows,time_comparison=errors,
        seconds=time.monotonic()-start,full_horizon_completed=False,free_material_reciprocal_closure_verified=False,
        final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'pilot-result.json',result);print(json.dumps(result),flush=True);assert result['passed'],result


def prepare():
    assert not OUT.exists();OUT.mkdir()
    failed=read(BEFORE/'block-result.json');assert not failed['passed']
    files=[Path(__file__),Path(previous.__file__),Path(previous.engine.__file__),Path(previous.engine.face.__file__),
        BEFORE/'block-result.json',BEFORE/'result.json']
    files.extend(p/f'steps-128-reference-128.npz' for s in [0,1] for p in previous.paths(s))
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='04bf81e3d',
        claim='Resolve the failed native noncollisional H feedback in the same photon/material equation before any GR charge readout.',
        failure='175 full histories pass their individual time/conservation tests but M_H and dM_H differ20.7percent and42.1percent; preserve the0.2percent rejection.',
        method='Reuse the native analytic conserved-to-primitive/HLL/donor tangent. On three saved states decompose the actual neutral flux difference and test additivity before selecting a sparse coupled-stage representation. A separately evolved gas H must not be overwritten to fake agreement.',
        decision='Only a representation that reproduces the native neutral flux, including branch effects, may enter an actual joint stage solve. Failure requires a source-level repair or explicitly nonlinear stage; no automatic further waveform sweep.',
        budgets=CAPS,CPU_threads=1,virtual_GiB=3,new_photon_or_material_steps=0,
        scope='Source-level preparation for a real coupled solve; neither final charge progress by itself nor a uniform derivative certificate.',
        gates=dict(native_representation=.002,conservation=1e-12),
        bindings={str(p):sha(p) for p in files}))
    H,C,M,T=sp.symbols('H C M T',cls=sp.Function);t=sp.symbols('t')
    assert sp.simplify(sp.diff(H(t)-C(t),t).subs(sp.diff(H(t),t),sp.diff(C(t),t)+T(t))-T(t))==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        scope='For Hdot=Cdot+T_H, the noncollisional inventory M=H-C satisfies Mdot=T_H. This algebra does not establish a sparse native derivative or a coupled numerical solution.'))


def localize():
    previous.initialize();m=previous.run.c.Material(128,128);rows=[];states=[]
    for s in [0,1]:
        d=np.load(previous.paths(s)[1]/'steps-128-reference-128.npz')
        ids=[np.argmin(abs(d['t']-t)) for t in m.t];states.append(d['history_scaled'][ids].copy())
    norm=lambda a:np.sum(abs(a),dtype=LD)
    for k in [1,8,16]:
        old,new=states[0][k],states[1][k]
        def transport(z):return -np.diff(previous.engine.tangent(m,k,z,1.)[0][3])
        oldrate,newrate=transport(old),transport(new);difference=newrate-oldrate
        changes=[]
        for j in range(4):
            z=old.copy();z[j]=new[j];changes.append(transport(z)-oldrate)
        zero=np.zeros_like(new);parts=[]
        for j in range(4):
            z=zero.copy();z[j]=new[j];parts.append(transport(z))
        scale=max(norm(newrate),LD('1e-290'))
        rows.append(dict(k=k,t=float(m.t[k]),deep_cells=m.nb,largest_changed_rate_cell=int(np.argmax(abs(difference))),
            neutral_rate_L1_physical=float(AMP*norm(newrate)),difference_L1_physical=float(AMP*norm(difference)),
            changed_component_L1_physical=[float(AMP*norm(v)) for v in changes],component_order=['B','S','Eref','H'],
            difference_decomposition_relative=float(norm(sum(changes)-difference)/max(norm(difference),LD('1e-290'))),
            zero_based_additivity_relative=float(norm(sum(parts)-newrate)/scale),
            neutral_only_rate_L1_physical=float(AMP*norm(parts[3])),
            conservation_relative=float(abs(np.sum(newrate,dtype=LD))/scale)))
    write(OUT/'localization.json',dict(classification='Counterexample candidate',rows=rows,no_evolution_steps=True))
    print(json.dumps(rows),flush=True)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS
    receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():
                target=OUT/'initial-producer.py' if Path(p)==Path(__file__) else p
                assert sha(target)==h,p
            if action=='pilot':
                for p,h in read(OUT/'stage-plan.json')['bindings'].items():assert sha(p)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
