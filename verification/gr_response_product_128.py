"""Shorter rigorous continued fractions and Q-independent coefficient reuse."""
from types import ModuleType
import json,sys
import gr_response_product_defined as previous

g=previous.original;OUT=previous.OUT.parent/'gr-response-product-128'


def sources():
    driver,_=previous.source();cpp,_=previous.cpp_source()
    dchanges={"'gr-response-product-defined'":"'gr-response-product-128'","'response-product-defined'":"'response-product-128'",
        'response_bits=256':'response_bits=128',
        'Evaluate the analytic logarithmic product moments at 256 bits.':'Evaluate the analytic logarithmic product moments at 128 bits.',
        'uses 160 extra levels.':'uses 16 extra levels.'}
    coefficient='''    for(int n=0;n<N;++n)for(int k=0;k<6;++k){b[n][k]+=c*panel.c[n][k];b[n+1][k]+=h*MI(n+1)/(2*n+1)*panel.c[n][k];if(n>0)b[n-1][k]+=h*MI(n)/(2*n+1)*panel.c[n][k];}'''
    cchanges={
        'MF::setDefaultPrecision(mode=="build"?128:256);':'MF::setDefaultPrecision(128);',
        'constexpr int padding=160;':'constexpr int padding=16;',
        'struct Stored{MI a,b;bool skip;Six error;Coeff c;};':'struct Stored{MI a,b,h,center;bool skip;Six error;Coeff c;std::array<Six,N+1> product{};};',
        'panels.push_back(panel);':'''panel.h=(panel.b-panel.a)/2;panel.center=(panel.a+panel.b)/2;exact(panel.h);exact(panel.center);
   if(!panel.skip){MI h=panel.h,c=panel.center;auto&b=panel.product;
'''+coefficient+'''
   }panels.push_back(panel);''',
        'MI h=(panel.b-panel.a)/2,c=(panel.a+panel.b)/2;exact(h);exact(c);':'MI h=panel.h,c=panel.center;',
        'std::array<Six,N+1>b{};\n'+coefficient:'const auto&b=panel.product;',
        '    for(int n=0;n<=N;++n)for(int k=0;k<6;++k)sum[k]+=h/(2*Q)*b[n][k]*(plus[n]-minus[n]);':
        '    MI scale=h/(2*Q);for(int n=0;n<=N;++n){MI moment=scale*(plus[n]-minus[n]);for(int k=0;k<6;++k)sum[k]+=b[n][k]*moment;}'}
    for old,new in dchanges.items():assert driver.count(old)==1;driver=driver.replace(old,new)
    # Apply the coefficient removal before its relocation adds the same statement.
    order=[key for key in cchanges if key!='panels.push_back(panel);']+['panels.push_back(panel);']
    for old in order:assert cpp.count(old)==1,old;cpp=cpp.replace(old,cchanges[old])
    return driver,cpp,dict(driver=dchanges,native=cchanges)


def module(text):
    obj=ModuleType('gr_response_product_128_candidate');exec(compile(text,str(OUT/'candidate.py'),'exec'),obj.__dict__);return obj


def prepare():
    previous.verify();assert not OUT.exists();text,cpp,changes=sources();obj=module(text);obj.prepare()
    (OUT/'candidate.py').write_text(text);(OUT/'effective-source.cpp').write_text(cpp)
    p=json.loads((OUT/'plan.json').read_text());p.update(checkpoint='a7b87ee',substitutions=changes,
        change='Use 128-bit outward arithmetic and 16 extra continued-fraction levels. The distant ratio is still the whole [0,rho], so validity holds independently of depth; the original 1e-15 response interval gate decides whether the width suffices. Precompute pH Legendre coefficients once per panel and share the scalar product moment across all six components. The 128-bit H build and its error budget are unchanged.',
        retained_reference='The original 2e-9 numerical vector comparison is retained only as a finite diagnostic; its convergence flags were not enforced. Independent cold-convolution audit is a separate stage. No physical or old reference verdict is promoted.',
        extra_controls='All H tables must replay bitwise. All 216 response intervals must overlap the original certified intervals. Same 54 moment and 66 endpoint-crossing controls. Report actual whole 12-Q timings, not only moment-kernel timing.')
    files=[g.ROOT/'verification/gr_response_product_128.py',OUT/'candidate.py',OUT/'effective-source.cpp',previous.OUT/'manifest.json',previous.OUT/'result.json']
    files+=list(previous.OUT.glob('h-*.tsv'))
    p['bindings'].update({x.relative_to(g.ROOT).as_posix():g.sha(x) for x in files});obj.save('plan.json',p)


def run(name):
    text,cpp,_=sources();assert text==(OUT/'candidate.py').read_text() and cpp==(OUT/'effective-source.cpp').read_text();obj=module(text);obj.bindings();getattr(obj,name)()


def controls():
    text,_,_=sources();obj=module(text);obj.bindings();old=json.loads((previous.OUT/'result.json').read_text())['results'];new=json.loads((OUT/'result.json').read_text())['results'];checks=[]
    for a,b in zip(old,new,strict=True):
        i=a['position'];assert i==b['position'];assert (OUT/f'h-{i:04d}.tsv').read_bytes()==(previous.OUT/f'h-{i:04d}.tsv').read_bytes()
        for x,y in zip(a['rows'],b['rows'],strict=True):
            for u,v in zip(x['enclosures'],y['enclosures'],strict=True):
                lo,hi=g.cusp.endpoints(u);left,right=g.cusp.endpoints(v);assert max(lo,left)<=min(hi,right),(i,x['z_index'])
        checks.append(dict(position=i,bitwise_H_replay=True,components=72,old_seconds=a['response_seconds'],new_seconds=b['response_seconds'],speedup=a['response_seconds']/b['response_seconds']))
    rows=[json.loads(x) for x in (OUT/'native-moments.jsonl').read_text().splitlines()];intervals=[x for x in rows if x['z'].startswith('interval')]
    assert len(intervals)==2 and all(x['passed'] for x in intervals)
    obj.save('reuse-controls.json',dict(classification='Counterexample candidate',passed=True,checks=checks,endpoint_components=66,
        scope='Three completed 12-Q jobs under the current shared host load. A finite timing comparison, not a whole-table EOS speed claim.'))


def finalize():controls();run('finalize');verify()


def verify():
    run('verify');assert json.loads((OUT/'reuse-controls.json').read_text())['passed']


if __name__=='__main__':
    if sys.argv[1] in ['prepare','controls','finalize','verify']:globals()[sys.argv[1]]()
    else:run(sys.argv[1])
