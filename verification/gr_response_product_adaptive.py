"""Bound discarded Legendre degrees before timing a separate response evaluator."""
from decimal import Decimal, localcontext
from fractions import Fraction as F
from pathlib import Path
import json, re, statistics, subprocess, sys
import sympy as sp
from mpmath import iv
import gr_outer_product_pilot as previous

ROOT=previous.ROOT;OUT=previous.OUT.parent/'gr-response-product-adaptive'
CACHE=previous.CACHE.parent/'response-product-adaptive';sha=previous.sha
g=previous.g;cusp=previous.cusp


def save(name,value):(OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def replace_once(text,old,new):
    assert text.count(old)==1,old
    return text.replace(old,new)


def source():
    text=(previous.OUT/'effective-source.cpp').read_text()
    start=text.index('std::vector<MI> log_moments(');end=text.index('\nstruct Panel',start)
    function=text[start:end]
    function=re.sub(r'\bN\b','degree',function)
    function=replace_once(function,'log_moments(MI z){','log_moments(MI z,int degree=N){')
    function=replace_once(function,' std::vector<MI> J(degree+1);MI az=abs(z);',
        ' if(degree<0||degree>N)throw std::runtime_error("moment degree");\n std::vector<MI> J(degree+1);MI az=abs(z);\n if(degree==0){J[0]=xlogabs(z+1)-xlogabs(z-1)-2;return J;}')
    text=text[:start]+function+text[end:]
    text=replace_once(text,'std::array<Six,N+1> product{};','std::array<Six,N+1> product{};std::array<MI,N+2> tail_norm{};')
    anchor='   }panels.push_back(panel);'
    text=replace_once(text,anchor,'''    Six suffix{};
    for(int n=N;n>=0;--n){MI maximum_tail(0);for(int k=0;k<6;++k){suffix[k]+=abs(b[n][k]);if(suffix[k].rightBound()>maximum_tail.rightBound())maximum_tail=MI(suffix[k].rightBound());}panel.tail_norm[n]=maximum_tail;}
   }panels.push_back(panel);''')
    text=replace_once(text,'  int j;std::string qtext;while(qs>>j>>qtext){MI Q=encoded(qtext,qtext);exact(Q);Six sum{};',
        '''  std::array<unsigned long long,N+1> histogram{};MI largest_omission(0);
  int j;std::string qtext;while(qs>>j>>qtext){MI Q=encoded(qtext,qtext);exact(Q);if(Q.leftBound()<MF(0))throw std::runtime_error("negative Q");Six sum{};MI omission(0);''')
    text=replace_once(text,'''    auto plus=log_moments(-(c+Q)/h),minus=log_moments((Q-c)/h);const auto&b=panel.product;
    MI scale=h/(2*Q);for(int n=0;n<=N;++n){MI moment=scale*(plus[n]-minus[n]);for(int k=0;k<6;++k)sum[k]+=b[n][k]*moment;}''',
        '''    MI zp=-(c+Q)/h,zm=(Q-c)/h,scale=h/(2*Q),budget=encoded("5e-19","5e-19")*(panel.b-panel.a)/Pcut;
    std::array<MI,2> ratio{},bound{};int kside=0;
    for(MI z:{zp,zm}){MI az=abs(z);if(az.leftBound()>=MF(1.25)){ratio[kside]=1/MI(az.leftBound());bound[kside]=2*ratio[kside]/(1-ratio[kside]);}++kside;}
    int degree=N;MI remainder(0);
    for(int n=1;n<=N+1;++n){MI weight(0);for(int side=0;side<2;++side){MI term(4);if(ratio[side].rightBound()>MF(0)){MI distant=bound[side]/n;if(distant.rightBound()<term.rightBound())term=MI(distant.rightBound());bound[side]*=ratio[side];}weight+=term;}
     MI candidate=scale*panel.tail_norm[n]*weight;
     if(candidate.rightBound()<=budget.leftBound()){degree=n-1;remainder=MI(candidate.rightBound());break;}}
    ++histogram[degree];omission+=remainder;
    auto plus=log_moments(zp,degree),minus=log_moments(zm,degree);const auto&b=panel.product;
    for(int n=0;n<=degree;++n){MI moment=scale*(plus[n]-minus[n]);for(int k=0;k<6;++k)sum[k]+=b[n][k]*moment;}''')
    text=replace_once(text,'   for(int k=0;k<6;++k)sum[k]+=symmetric(Pcut*maximum[k]);',
        '''   if(omission.rightBound()>encoded("1e-18","1e-18").leftBound())throw std::runtime_error("kernel remainder budget");
   if(omission.rightBound()>largest_omission.rightBound())largest_omission=MI(omission.rightBound());
   for(int k=0;k<6;++k)sum[k]+=symmetric(Pcut*maximum[k])+symmetric(omission);''')
    old='  if(!out)throw std::runtime_error("response output");std::cout<<"{\\"seconds\\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<"}\\n";'
    new='''  if(!out)throw std::runtime_error("response output");std::cout<<"{\\"seconds\\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<",\\"kernel_remainder_upper\\":\\""<<endpoint(largest_omission.rightBound(),true)<<"\\",\\"degree_histogram\\":[";for(int n=0;n<=N;++n){if(n)std::cout<<',';std::cout<<histogram[n];}std::cout<<"]}\\n";'''
    text=replace_once(text,old,new)
    text=replace_once(text,'auto J=log_moments(fraction(z));out<<',
        '''auto J=log_moments(fraction(z));for(int d=0;d<=N;++d){auto prefix=log_moments(fraction(z),d);for(int n=0;n<=d;++n)if(prefix[n].leftBound()>J[n].rightBound()||J[n].leftBound()>prefix[n].rightBound())throw std::runtime_error("degree prefix overlap");}out<<''')
    assert 'MI Q=encoded(qtext,qtext);exact(Q);' in text
    return text


def decimal_exact(q):
    q=F(q)
    assert q.denominator & (q.denominator-1)==0
    with localcontext() as ctx:
        ctx.prec=2000
        text=format(Decimal(q.numerator)/Decimal(q.denominator),'f')
    assert F(text)==q
    return text


def prepare():
    previous.verify();assert not OUT.exists() and not CACHE.exists();OUT.mkdir();CACHE.mkdir()
    (OUT/'effective-source.cpp').write_text(source())
    files=[ROOT/'verification/gr_response_product_adaptive.py',OUT/'effective-source.cpp',previous.OUT/'manifest.json',
        previous.OUT/'runtime.json',previous.product.OUT/'inputs.tsv',previous.product.OUT/'rule.tsv',
        g.rule.OUT/'controls.json',g.rule.OUT/'tail.json',previous.product.OUT/'result.json']
    selections=[]
    for i in [0,1603,3205]:
        qold=previous.product.OUT/f'q-{i:04d}.tsv';lines=[]
        for line in qold.read_text().splitlines():
            j,q=line.split();lines.append(f'{j} {decimal_exact(F.from_float(float.fromhex(q)))}\n')
        qpath=OUT/f'control-q-{i:04d}.tsv';qpath.write_text(''.join(lines));files += [qold,qpath]
        qouter=previous.OUT/f'cell-{i:04d}-q.tsv';allrows=qouter.read_text().splitlines();indices=[k*(len(allrows)-1)//47 for k in range(48)]
        assert len(set(indices))==48
        target=OUT/f'benchmark-q-{i:04d}.tsv';target.write_text(''.join(allrows[j]+'\n' for j in indices))
        selections.append(dict(position=i,indices=indices,total_nodes=len(allrows)));files += [qouter,target,previous.product.OUT/f'h-{i:04d}.tsv']
    p=previous.bindings()
    save('plan.json',dict(classification='Proven',checkpoint='6c29d985',positions=[0,1603,3205],native_bits=128,
        kernel_budget='1e-18',allocated_kernel_budget='5e-19',response_midpoint_budget='1e-15',physical_field_budget='2e-7',
        selections=selections,repeats=3,order='For each position and repeat, baseline/adaptive if (position_index+repeat) is even, adaptive/baseline otherwise. No warmup exclusions.',
        bindings={f.relative_to(ROOT).as_posix():sha(f) for f in files},interval_headers=p['interval_headers'],
        target='Retain the exact H polynomial and all original interpolation/tail errors. Bound discarded pH Legendre degrees panel by panel before evaluation. A zero tail at d=32 is the unconditional fallback. Preserve the exact 80-bit Q parser.',
        theorem='For n>=1, |Jn(z)|<=4 on the whole real line. If |z|>=5/4, |Jn(z)|<=2/[n*|z|^n*(1-1/|z|)]. Orthogonality removes the constant and powers below n. Bound the omitted panel response by h/(2Q)*max_k sum_(n>d)|b_n,k| times the sum of the two (d+1) kernel bounds. Coefficient magnitudes and all arithmetic are outward intervals.',
        boundary='These are finite response controls and timings at three frozen states, not a new all-state integral certificate or physical EOS/GR/observation closure. The existing full-table producer remains untouched.'))
    x=sp.symbols('x');checks=0
    for n in range(1,33):
        poly=sp.Poly(sp.legendre(n,x),x)
        for k in range(n):
            integral=sum(c*sp.Rational(1-(-1)**(degree[0]+k+1),degree[0]+k+1) for degree,c in poly.terms())
            assert integral==0;checks+=1
    assert sp.Rational(9,4)<sum(sp.Rational(1,sp.factorial(k)) for k in range(3))
    assert 5<sum(sp.Rational(2)**k/sp.factorial(k) for k in range(4))
    save('proof.json',dict(classification='Proven',passed=True,exact_orthogonality_checks=checks,
        near='For |z|<=5/4, split |log|z-x|| into negative and positive parts. Their integrals are <=2 and <=2 log(9/4). Since e>5/2>9/4, the sum is <4. Use |Pn(x)|<=1 from the previously certified Legendre representation.',
        far_uniform='For |z|>=5/4 subtract log|z| using orthogonality. Then |log(1-x/z)|<=log 5 and 2 log 5<4 since e^2>19/3>5. This proves a universal bound even for intervals crossing 5/4.',
        far_degree='The uniformly absolutely convergent log series has no contributions for powers k<n. Bound each integral by2, 1/k<=1/n and the remaining geometric series. Its bound decreases with n, so the first omitted degree bounds every discarded moment.',
        allocation='The panel bound uses a suffix sum of outward absolute coefficient bounds, maximized over six jets. Its upper endpoint must be <= the lower endpoint of 5e-19*(b-a)/P. Exact panel coverage implies the mathematical sum <=5e-19. The computed outward sum must additionally pass1e-18 and is explicitly added to every response component.',
        recurrence='Every requested prefix retains the whole continued-fraction terminal interval [0,rho]. Changing the finite depth changes tightness only. Endpoint and near-endpoint formulas are unchanged; d=0 uses the continuous analytic J0 directly.'))
    print('PREPARED adaptive degree proof and frozen paired benchmark',flush=True)


def bindings():
    p=json.loads((OUT/'plan.json').read_text());assert (OUT/'effective-source.cpp').read_text()==source()
    for rel,digest in p['bindings'].items():assert sha(ROOT/rel)==digest,rel
    for path,digest in p['interval_headers'].items():assert sha(Path(path))==digest,path
    return p


def engine():
    obj=previous.engine();obj.OUT=OUT;obj.CACHE=CACHE;obj.bindings=bindings;return obj


def build():engine().build()


def native_call(binary,position,qpath,prefix,mode='response',hpath=None):
    target=OUT/(prefix+('.tsv' if mode=='build' else '.jsonl'));assert not target.exists()
    cmd=[binary,mode,str(previous.product.OUT/'inputs.tsv'),str(previous.product.OUT/'rule.tsv'),str(target),str(position)]
    if mode=='response':cmd += [str(hpath or previous.product.OUT/f'h-{position:04d}.tsv'),str(qpath)]
    done=subprocess.run(cmd,capture_output=True,text=True);(OUT/(prefix+'.log')).write_text(done.stdout+done.stderr)
    assert done.returncode==0,(prefix,done.stderr)
    return target,json.loads(done.stdout)


def overlap(a,b):return max(F(a[0]),F(b[0]))<=min(F(a[1]),F(b[1]))


def compare_files(left,right):
    count=0
    for a,b in zip(left.read_text().splitlines(),right.read_text().splitlines(),strict=True):
        a=json.loads(a);b=json.loads(b)
        assert (a['position'],a['cell'],a['z_index'])==(b['position'],b['cell'],b['z_index'])
        for x,y in zip(a['finite_response'],b['finite_response'],strict=True):assert overlap(x,y),(a['position'],a['z_index'],x,y);count+=1
    return count


def point_errors(path,qpath):
    qvalues={int(j):F(q) for j,q in (line.split() for line in qpath.read_text().splitlines())}
    tails=json.loads((g.rule.OUT/'tail.json').read_text())['rows'];maximum=F(0)
    for row in map(json.loads,path.read_text().splitlines()):
        q=qvalues[row['z_index']]
        for (a,b),t in zip(row['finite_response'],tails[row['position']]['upper_coefficients'],strict=True):
            tail=F(t)/max(F(1),q*q);lo=F(a)-tail;hi=F(b)+tail;mid=F.from_float(float((lo+hi)/2))
            maximum=max(maximum,abs(lo-mid),abs(hi-mid))
    assert maximum<F(bindings()['response_midpoint_budget']),str(maximum)
    return str(maximum)


def controls():
    p=bindings();obj=engine();binary=obj.runtime();obj.moment_controls();moments=list(map(json.loads,(OUT/'native-moments.jsonl').read_text().splitlines()))
    assert len(moments)==11 and all(x['passed'] for x in moments if x['z'].startswith('interval'))
    old=json.loads((previous.product.OUT/'result.json').read_text())['results'];rows=[]
    for i in p['positions']:
        h,_=native_call(binary,i,None,f'h-replay-{i:04d}',mode='build');assert h.read_bytes()==(previous.product.OUT/f'h-{i:04d}.tsv').read_bytes()
        target,timing=native_call(binary,i,OUT/f'control-q-{i:04d}.tsv',f'control-{i:04d}')
        assert F(timing['kernel_remainder_upper'])<=F(p['kernel_budget'])
        count=compare_files(target,previous.product.OUT/f'response-{i:04d}.jsonl')
        error=point_errors(target,OUT/f'control-q-{i:04d}.tsv')
        # The old total intervals already include the same entire H tail.
        reference=next(x for x in old if x['position']==i)
        assert count==len(reference['rows'])*6==72
        rows.append(dict(position=i,bitwise_H_replay=True,components=count,maximum_midpoint_error=error,**timing))
        print('ADAPTIVE CONTROL',i,count,float(F(error)),timing,flush=True)
    save('controls.json',dict(classification='Counterexample candidate',passed=True,rows=rows,independent_moment_components=54,endpoint_components=66,prefix_overlap_components=5049,
        boundary='Prefix overlaps and response overlaps are necessary finite compatibility controls. The rigorous truncation argument is separate in proof.json. H byte identity retains the earlier independent H controls.'))


def benchmark():
    p=bindings();assert json.loads((OUT/'controls.json').read_text())['passed'];new=engine().runtime();old=previous.engine().runtime();rows=[]
    for index,i in enumerate(p['positions']):
        trials=[]
        for repeat in range(p['repeats']):
            order=['baseline','adaptive'] if (index+repeat)%2==0 else ['adaptive','baseline'];trial=dict(repeat=repeat,order=order)
            paths={}
            for name in order:
                path,timing=native_call(old if name=='baseline' else new,i,OUT/f'benchmark-q-{i:04d}.tsv',f'benchmark-{i:04d}-{repeat}-{name}')
                paths[name]=path;trial[name]=timing
            assert compare_files(paths['baseline'],paths['adaptive'])==48*6
            trial['maximum_midpoint_error']=point_errors(paths['adaptive'],OUT/f'benchmark-q-{i:04d}.tsv')
            assert F(trial['adaptive']['kernel_remainder_upper'])<=F(p['kernel_budget'])
            trials.append(trial);print('ADAPTIVE PAIRED TIMING',i,repeat,trial['baseline']['seconds'],trial['adaptive']['seconds'],flush=True)
        a=statistics.median(t['baseline']['seconds'] for t in trials);b=statistics.median(t['adaptive']['seconds'] for t in trials)
        rows.append(dict(position=i,baseline_median_seconds=a,adaptive_median_seconds=b,median_ratio=a/b,trials=trials))
    save('benchmark.json',dict(classification='Counterexample candidate',passed=True,rows=rows,components=3*3*48*6,
        scope='Three paired repeats at48 frozen exact outer nodes per state, under shared host load on one inherited CPU core. No whole3206-state or physical inference performance claim.'))


def finalize():
    bindings();engine().runtime()
    for name in ['proof.json','controls.json','moment-controls.json','benchmark.json']:assert json.loads((OUT/name).read_text())['passed'],name
    save('manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    p=bindings();engine().runtime()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    for name in ['proof.json','controls.json','moment-controls.json','benchmark.json']:assert json.loads((OUT/name).read_text())['passed'],name
    r=json.loads((OUT/'controls.json').read_text());assert sum(x['components'] for x in r['rows'])==216 and r['prefix_overlap_components']==5049
    for x in r['rows']:assert x['bitwise_H_replay'] and F(x['maximum_midpoint_error'])<F(p['response_midpoint_budget'])
    r=json.loads((OUT/'benchmark.json').read_text());assert r['components']==2592
    for x in r['rows']:
        assert len(x['trials'])==p['repeats']
        for t in x['trials']:assert F(t['maximum_midpoint_error'])<F(p['response_midpoint_budget']) and F(t['adaptive']['kernel_remainder_upper'])<=F(p['kernel_budget'])
    print('PASS rigorous discarded-degree bound,54 independent moments,5049 prefix overlaps,216 response controls and2592 paired exact-Q components',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
