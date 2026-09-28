"""Evaluate actual finite outer Gauss nodes at three preregistered states."""
from concurrent.futures import ProcessPoolExecutor,as_completed
from decimal import Decimal,localcontext
from fractions import Fraction as F
from pathlib import Path
from types import ModuleType
import gzip,json,math,shlex,shutil,subprocess,sys,time
import numpy as np
from mpmath import iv
import gr_self_regular_domain as domain
import gr_polarization_infrared as infrared

highq=domain.highq;cusp=domain.cusp;ROOT=domain.ROOT;OUT=highq.OUT.parent/'gr-outer-self-pilot';CACHE=highq.moments.native.g.CACHE/'outer-self-pilot'
I=domain.I;low=domain.low;high=domain.high;sha=domain.sha;SMOOTH=highq.OUT.parent/'gr-smooth-response-defined'

def save(name,value):(OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    domain.verify();infrared.verify();assert not OUT.exists() and not CACHE.exists();OUT.mkdir();CACHE.mkdir()
    cpp=(SMOOTH/'effective-source.cpp').read_text();changes={'MF::setDefaultPrecision(128)':'MF::setDefaultPrecision(256)',
        'MI beta=binary(bs),Q=binary(qs)':'MI beta=binary(bs),Q=encoded(qs,qs)',
        'encoded("1e-13","1e-13")':'encoded("1e-18","1e-18")'}
    for a,b in changes.items():assert cpp.count(a)==1;cpp=cpp.replace(a,b)
    (OUT/'finite-source.cpp').write_text(cpp);shutil.copy2(SMOOTH/'rule.tsv',OUT/'rule.tsv')
    text=(ROOT/'verification/gr_logarithmic_response_certificate.py').read_text();cchanges={
        "Q=I(DATA['Q'][i,j])":"Q=I(DATA['Q_exact'][i,j])",
        "Q_hex=float(DATA['Q'][i,j]).hex()":"Q_exact=str(DATA['Q_exact'][i,j])"}
    for a,b in cchanges.items():assert text.count(a)==1;text=text.replace(a,b)
    (OUT/'cusp-source.py').write_text(text)
    files=[ROOT/'verification/gr_outer_self_pilot.py',OUT/'finite-source.cpp',OUT/'cusp-source.py',OUT/'rule.tsv',
        domain.OUT/'manifest.json',domain.OUT/'result.json',SMOOTH/'manifest.json',SMOOTH/'rule.json',cusp.OUT/'plan.json',cusp.rule.OUT/'result.json',
        highq.OUT/'manifest.json',infrared.OUT/'manifest.json',cusp.ROOTS/'result.json',cusp.THERMO/'states.npz',highq.moments.OUT/'inputs.npz',
        highq.neutral.density.ionic.OUT/'inputs.npz',highq.neutral.density.ionic.OUT/'constants.json']
    save('plan.json',dict(classification='Proven',checkpoint='62fdec7',bits=256,node_binary_bits=80,positions=[0,1603,3205],preflight_nodes=16,
        cusp_radius_ratio=128,cusp_component_budget='1e-18',finite_remainder_budget='1e-18',physical_field_budget='2e-7',processes=3,
        bindings={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        interval_headers=json.loads((highq.moments.native.OUT/'plan.json').read_text())['interval_headers'],
        finite_source_substitutions=changes,cusp_source_substitutions=cchanges,
        node_binding='Use the certified 16-point interval Gauss rule on each exact dyadic Q panel. Round its interval midpoint to an exact 80-bit dyadic Q. The native momentum evaluator reads that exact decimal dyadic at 256 bits. Enclose the displacement to the true algebraic Gauss node using the joint Q Cauchy bound M/(R-half_width-displacement).',
        momentum='Reuse the certified complete cusp, finite complement and positive tail formulas. Narrow the cusp to R/128 and reduce the finite analytic remainder to 1e-18. The true eta interval is retained. Both finite and cusp sources are bound after explicit substitutions; exact dyadic coverage checks remain active.',
        assembly='Propagate all six full response intervals through G at the exact Q. Add the Gauss-node displacement before the positive weighted sum. Convert the raw integral with the certified neutral matrix and exact inventory, add the pre-certified outer remainder, then add infrared and high-Q field intervals without overlap.',
        boundary='Three selected actual states only. A successful pilot does not certify all 3206 states, the full correlated EOS, GR evolution or observations.'))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert sha(ROOT/rel)==digest,rel
    for rel,digest in p['interval_headers'].items():assert sha(Path(rel))==digest,rel
    return p


def build():
    bindings();native=highq.moments.native;flags=shlex.split(subprocess.check_output([str(native.CAPD/'build-request15-mp/bin/capd-config'),'--cflags','--libs'],text=True))
    cmd=['g++',str(OUT/'finite-source.cpp'),f'-I{native.DEPS}/include',f'-I{native.DEPS}/include/x86_64-linux-gnu',f'-L{native.DEPS}/lib/x86_64-linux-gnu',*flags,'-O2','-o',str(CACHE/'finite')]
    done=subprocess.run(cmd,capture_output=True,text=True);(OUT/'build.log').write_text(done.stdout+done.stderr);save('build.json',dict(command=cmd,returncode=done.returncode));assert done.returncode==0,done.stderr
    linked=subprocess.check_output(['ldd',str(CACHE/'finite')],text=True);(OUT/'linked-libraries.txt').write_text(linked)
    libs={word:sha(Path(word)) for line in linked.splitlines() for word in line.split() if word.startswith('/') and Path(word).is_file()}
    save('runtime.json',dict(binary=str(CACHE/'finite'),binary_sha256=sha(CACHE/'finite'),libraries=libs))


def runtime():
    bindings();r=json.loads((OUT/'runtime.json').read_text());assert sha(Path(r['binary']))==r['binary_sha256']
    for path,digest in r['libraries'].items():assert sha(Path(path))==digest,path
    return r['binary']


def dyadic_midpoint(value,bits):
    lo,hi=cusp.endpoints(cusp.interval_text(value));mid=(lo+hi)/2;assert mid>0
    step=domain.dyadic_below(mid)/2**(bits-1);q=((mid/step+F(1,2)).numerator//(mid/step+F(1,2)).denominator)*step
    with localcontext() as ctx:
        ctx.prec=512;text=str(Decimal(q.numerator)/Decimal(q.denominator))
    assert F(text)==q;return q,text,max(abs(q-lo),abs(hi-q))


def nodes(schedule,row,plan):
    rule=json.loads((SMOOTH/'rule.json').read_text())['nodes'];out=[]
    for a,b in schedule['panels']:
        lo,hi=F(a),F(b);center=(lo+hi)/2;h=(hi-lo)/2;qD=F(row['screening_disk_radius']);s=F(row['strip_radius_lower'])
        R=max(qD if center<=2*qD else min(s,center/2),center/256)
        for point in rule:
            exact=I(lo)+I(hi-lo)*cusp.interval(point['node']);q,text,error=dyadic_midpoint(exact,plan['node_binary_bits']);assert error<R-h
            out.append(dict(Q=str(q),Q_decimal=text,weight=cusp.interval_text(I(hi-lo)*cusp.interval(point['weight'])),
                shift_factor=str(high(I(error)/I(R-h-error)))))
    assert len(out)==schedule['nodes'];return out


def evaluate(position,limit=None):
    plan=bindings();binary=runtime();iv.prec=plan['bits'];began=time.perf_counter();meta=json.loads((domain.OUT/'result.json').read_text())
    schedule=next(x for x in meta['schedules'] if x['position']==position);row=meta['rows'][position];all_nodes=nodes(schedule,row,plan)
    selected=all_nodes if limit is None else all_nodes[:limit];prefix=f'cell-{position:04d}' if limit is None else f'preflight-{position:04d}'
    state,ions,roots,data,B,_=highq.context();iv.prec=plan['bits'];root=roots[position]
    # Recreate alpha/pi at the requested precision after the shared context loader.
    const=json.loads((highq.neutral.density.ionic.OUT/'constants.json').read_text())['native_binary64'];B=I(const['alpha'])/iv.pi
    obj=ModuleType('outer_cusp');exec(compile((OUT/'cusp-source.py').read_text(),str(OUT/'cusp-source.py'),'exec'),obj.__dict__)
    obj.PLAN=json.loads((cusp.OUT/'plan.json').read_text());obj.PLAN.update(bits=plan['bits'],window_radius_ratio=plan['cusp_radius_ratio'],per_component_absolute_error_budget=plan['cusp_component_budget'])
    obj.DATA=dict(cells=np.array([root['cell']]),beta=np.array([data['beta'][position]]),Sref=np.array([data['Sref'][position]]),root_intervals=np.array([root['root']]),Q_exact=np.array([[x['Q'] for x in selected]]))
    obj.RULES=[dict(norm=I(r['norm_exact']),nodes=[(cusp.interval(n['node']),cusp.interval(n['weight'])) for n in r['nodes']]) for r in json.loads((cusp.rule.OUT/'result.json').read_text())['rules']]
    input_path=OUT/f'{prefix}-inputs.tsv';cusp_path=OUT/f'{prefix}-cusp.jsonl.gz';assert not input_path.exists() and not cusp_path.exists()
    with input_path.open('x') as inputs,gzip.open(cusp_path,'xt',encoding='utf-8',compresslevel=1) as stream:
        for j,node in enumerate(selected):
            window=obj.certified_window(0,j);window['position']=position;stream.write(json.dumps(window)+'\n');elo,ehi=root['root'][1:-1].split(',')
            inputs.write(' '.join(map(str,[j,root['cell'],j,float(data['beta'][position]).hex(),node['Q_decimal'],float(data['Sref'][position]).hex(),
                elo.strip(),ehi.strip(),float(F(window['window_half_width'])).hex(),float(data['pcut'][position]).hex()]))+'\n')
            if j%128==0:print('OUTER CUSP',position,j,len(selected),flush=True)
    finite_path=OUT/f'{prefix}-finite.jsonl';assert not finite_path.exists()
    done=subprocess.run([binary,str(input_path),str(OUT/'rule.tsv'),str(finite_path),'0',str(len(selected))],capture_output=True,text=True)
    (OUT/f'{prefix}-native.log').write_text(done.stdout+done.stderr);assert done.returncode==0,done.stderr
    beta=I(data['beta'][position]);Sref=I(data['Sref'][position]);scale=I(state['scale'][position]);eta_hi=cusp.endpoints(root['root'])[1];tc=64+max(eta_hi,F(0));T=I(tc)
    thermal=[iv.exp(I(eta_hi-tc))*sum(I(F(math.factorial(n),math.factorial(k)))*T**k for k in range(n+1)) for n in range(6)]
    raw=[I(0)]*6;point_max=[F(0)]*6;point_path=OUT/f'{prefix}-points.jsonl.gz';finite_panels=0
    with gzip.open(cusp_path,'rt',encoding='utf-8') as cs,finite_path.open() as fs,gzip.open(point_path,'xt',encoding='utf-8',compresslevel=1) as points:
        for j,(cl,fl,node) in enumerate(zip(cs,fs,selected,strict=True)):
            window=json.loads(cl);finite=json.loads(fl);assert finite['case']==j and window['z_index']==j;finite_panels+=finite['panels']
            Q=I(node['Q']);h=I(window['window_half_width']);L=iv.log((2+h)/h);tail=[]
            for k in range(3):tail.append((beta*(thermal[k]+beta*thermal[k+1])+L/(2*Q)*iv.sqrt(beta/(2*T))*(thermal[k]+3*beta*thermal[k+1]+3*beta**2*thermal[k+2]+beta**3*thermal[k+3]+Q*Q*(thermal[k]+beta*thermal[k+1])))/Sref)
            bounds=[tail[0],tail[0],tail[0],tail[1],tail[1],tail[2]+tail[1]]
            response=[cusp.interval(w)+iv.mpf([I(a).a,I(b).b])+cusp.symmetric(high(t)) for w,(a,b),t in zip(window['enclosures'],finite['enclosures'],bounds)]
            X=Q*Q/(B*Sref);den=X+response[0];assert low(response[0])>0
            R=response;G=[R[0]/den,X*R[1]/den**2,X*(R[2]/den**2-2*R[1]**2/den**3),X*R[3]/den**2,
                X*(R[4]/den**2-2*R[1]*R[3]/den**3),X*(R[5]/den**2-2*R[3]**2/den**3)]
            shift=I(node['shift_factor']);G=[x+cusp.symmetric(high(shift*I(m))) for x,m in zip(G,meta['derivative_majorants'])]
            weight=cusp.interval(node['weight']);raw=[old+weight*x/scale for old,x in zip(raw,G)]
            point_max=[max(old,high(x)-low(x)) for old,x in zip(point_max,G)]
            points.write(json.dumps(dict(node=j,Q_exact=node['Q'],response=list(map(cusp.interval_text,response)),G=list(map(cusp.interval_text,G))))+'\n')
    with gzip.open(finite_path.with_suffix('.jsonl.gz'),'xb',compresslevel=1) as stream:stream.write(finite_path.read_bytes())
    finite_path.unlink()
    result=dict(classification='Proven',position=position,cell=root['cell'],nodes=len(selected),finite_momentum_panels=finite_panels,
        partial_raw_integrals=list(map(cusp.interval_text,raw)),G_maximum_width=list(map(str,point_max)),elapsed_seconds=time.perf_counter()-began,
        complete_finite_interval=limit is None,full_physical_EOS_certified=False)
    if limit is None:
        coeff=[cusp.interval(root['derivatives'][name]) for name in highq.neutral.density.FIELDS];counts=[I(x)/I(a) for x,a in zip(ions['X'][position],ions['A'])]
        Z2=sum(x*I(z)**2 for x,z in zip(counts,ions['Z']))/sum(counts);pref=2*B*Z2*scale/beta
        remainder=F(schedule['field_remainder_upper']);middle=[-pref*sum(c*x for c,x in zip(line,raw))+cusp.symmetric(remainder) for line in highq.neutral.matrix(coeff)]
        ir=json.loads((infrared.OUT/'result.json').read_text())['rows'][position];hq=next(x for x in highq.records() if x['position']==position)
        total=[x+cusp.interval(a)+cusp.interval(b) for x,a,b in zip(middle,ir['field_enclosures'],hq['field_enclosures'])];texts=list(map(cusp.interval_text,total));approx=[];errors=[]
        for text in texts:
            lo,hi=cusp.endpoints(text);point=float((lo+hi)/2);approx.append(point);errors.append(max(abs(F.from_float(point)-lo),abs(hi-F.from_float(point))))
        reference=np.asarray(state['fine_fields'])[:,position];difference=float(np.max(np.abs(np.array(approx)-reference)))
        result.update(passed=max(errors)<F(plan['physical_field_budget']),middle_field_enclosures=list(map(cusp.interval_text,middle)),
            complete_field_enclosures=texts,field_approximation=approx,field_error_upper=list(map(str,errors)),
            finite_reference_control=dict(classification='Counterexample candidate',maximum_absolute_difference=difference,passed=difference<float(plan['physical_field_budget']),
                scope='Compare with the frozen independent numerical thermodynamic table at its stored neutral midpoint/coefficients. This is a finite numerical control, not a rigorous reference enclosure.'))
    save(prefix+'-result.json',result);return result


def preflight():
    plan=bindings()
    for i in plan['positions']:
        r=evaluate(i,plan['preflight_nodes']);print('OUTER PREFLIGHT',i,r['elapsed_seconds'],r['finite_momentum_panels'],flush=True)


def run():
    p=bindings()
    with ProcessPoolExecutor(max_workers=p['processes']) as pool:
        jobs=[pool.submit(evaluate,i) for i in p['positions']]
        for future in as_completed(jobs):
            r=future.result();print('COMPLETE OUTER PILOT',r['position'],r['passed'],r['finite_reference_control'],flush=True)


def finalize():
    runtime();plan=bindings();results=[json.loads((OUT/f'cell-{i:04d}-result.json').read_text()) for i in plan['positions']]
    assert all(r['passed'] and r['finite_reference_control']['passed'] for r in results)
    save('manifest.json',dict(sha256={path.relative_to(ROOT).as_posix():sha(path) for path in OUT.iterdir() if path.is_file()}));verify()


def verify():
    p=bindings();runtime()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    for i in p['positions']:
        r=json.loads((OUT/f'cell-{i:04d}-result.json').read_text());assert r['passed'] and r['complete_finite_interval'] and r['finite_reference_control']['passed']
        for text,point,bound in zip(r['complete_field_enclosures'],r['field_approximation'],r['field_error_upper'],strict=True):
            lo,hi=cusp.endpoints(text);assert max(abs(lo-F.from_float(point)),abs(hi-F.from_float(point)))==F(bound)<F(p['physical_field_budget'])
    print('PASS three actual complete self-integral pilots and seven fields; remaining states/full EOS/GR/observations open',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
