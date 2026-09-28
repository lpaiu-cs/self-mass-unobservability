"""Preserve the first failure; handle serialized boundaries and rational log coordinates."""
from types import ModuleType
import json,sys
import sympy as sp
import gr_response_product_native as original

OUT=original.OUT.parent/'gr-response-product-defined'


def cpp_source():
    before=(original.ROOT/'verification/gr_response_product_native.cpp').read_text()
    oldnear='''  MI q0=log(abs((z+1)/(z-1)))/2;std::vector<MI>P(N+2),W(N+2);P[0]=1;P[1]=z;W[0]=0;W[1]=1;
  for(int n=1;n<=N;++n){P[n+1]=((2*n+1)*z*P[n]-n*P[n-1])/(n+1);W[n+1]=((2*n+1)*z*W[n]-n*W[n-1])/(n+1);}
  J[0]=(z+1)*log(abs(z+1))-(z-1)*log(abs(z-1))-2;
  for(int n=1;n<=N;++n)J[n]=2*((P[n+1]-P[n-1])*q0-W[n+1]+W[n-1])/(2*n+1);'''
    newnear='''  std::vector<MI>P(N+2),D(N+2),W(N+2);P[0]=1;P[1]=z;D[0]=0;D[1]=1;W[0]=0;W[1]=1;
  for(int n=1;n<=N;++n){P[n+1]=((2*n+1)*z*P[n]-n*P[n-1])/(n+1);D[n+1]=((2*n+1)*(P[n]+z*D[n])-n*D[n-1])/(n+1);W[n+1]=((2*n+1)*z*W[n]-n*W[n-1])/(n+1);}
  MI left=xlogabs(z+1),right=xlogabs(z-1),B=((z-1)*left-(z+1)*right)/2;J[0]=left-right-2;
  for(int n=1;n<=N;++n)J[n]=2*D[n]*B/(n*(n+1))-2*(W[n+1]-W[n-1])/(2*n+1);'''
    helpers='''MI hull(const MI&a,const MI&b){return MI(a.leftBound()<b.leftBound()?a.leftBound():b.leftBound(),a.rightBound()>b.rightBound()?a.rightBound():b.rightBound());}
MI xlogpoint(const MF&x){if(x==MF(0))return MI(0);MI p(x);return p*log(abs(p));}
MI xlogabs(const MI&x){
 MI value=hull(xlogpoint(x.leftBound()),xlogpoint(x.rightBound())),e=exp(MI(-1));
 if(x.leftBound()<=MF(0)&&x.rightBound()>=MF(0))value=hull(value,MI(0));
 if(x.leftBound()<=e.rightBound()&&x.rightBound()>=e.leftBound())value=hull(value,-e);
 MI negative=-e;if(x.leftBound()<=negative.rightBound()&&x.rightBound()>=negative.leftBound())value=hull(value,e);return value;
}
'''
    controls='''  for(int sign:{-1,1}){MI center(sign),eps=MI(1)/power(MI(2),80),lo=center-eps,hi=center+eps,box(lo.leftBound(),hi.rightBound());auto enclosure=log_moments(box);
   for(MI sample:{lo,center,hi}){auto point_values=log_moments(sample);for(int n=0;n<=N;++n)if(enclosure[n].leftBound()>point_values[n].leftBound()||enclosure[n].rightBound()<point_values[n].rightBound())throw std::runtime_error("endpoint-crossing interval control");}
   out<<"{\\\"z\\\":\\\"interval"<<sign<<"\\\",\\\"passed\\\":true,\\\"moments\\\":[";
   for(int n=0;n<=N;++n){if(n)out<<',';out<<"[\\\""<<endpoint(enclosure[n].leftBound(),false)<<"\\\",\\\""<<endpoint(enclosure[n].rightBound(),true)<<"\\\"]";}out<<"]}\\n";}
'''
    # Keep ordinary C++ JSON escaping; this Python string is written once to the frozen source.
    controls=controls.replace('\\\\','\\')
    changes={
        'auto s=endpoint(x.leftBound(),false,200);MI y=encoded(s,s);':'auto s=endpoint(x.leftBound(),false,200);std::istringstream tokens(s);tokens>>s;MI y=encoded(s,s);',
        'std::vector<MI> log_moments(MI z){':helpers+'std::vector<MI> log_moments(MI z){',
        ' exact(z);std::vector<MI> J(N+1);MI az=abs(z);':' std::vector<MI> J(N+1);MI az=abs(z);',
        'if(az.leftBound()==MF(1)){':'if(az.leftBound()==MF(1)&&az.rightBound()==MF(1)){',
        oldnear:newnear,
        '  if(!out)throw std::runtime_error("moment output");return 0;':controls+'  if(!out)throw std::runtime_error("moment output");return 0;'}
    after=before
    for a,b in changes.items():assert after.count(a)==1;after=after.replace(a,b)
    reverse=after
    for a,b in reversed(list(changes.items())):assert reverse.count(b)==1;reverse=reverse.replace(b,a)
    assert reverse==before;return after,changes


def source():
    before=(original.ROOT/'verification/gr_response_product_native.py').read_text()
    changes={"OUT=rule.OUT.parent/'gr-response-product-native';native=rule.highq.moments.native;CACHE=native.g.CACHE/'response-product-native'":
        "OUT=rule.OUT.parent/'gr-response-product-defined';native=rule.highq.moments.native;CACHE=native.g.CACHE/'response-product-defined'",
        "str(ROOT/'verification/gr_response_product_native.cpp')":"str(OUT/'effective-source.cpp')",
        'All Q coordinates and panel-derived z coordinates must be exact dyadics.':'Q and panel boundaries remain exact dyadics; derived rational z coordinates are enclosed outward.'}
    after=before
    for a,b in changes.items():assert after.count(a)==1;after=after.replace(a,b)
    reverse=after
    for a,b in reversed(list(changes.items())):reverse=reverse.replace(b,a)
    assert reverse==before;return after,changes


def module(text):
    obj=ModuleType('gr_response_product_defined_candidate');exec(compile(text,str(OUT/'candidate.py'),'exec'),obj.__dict__);return obj


def prepare():
    assert not OUT.exists();text,changes=source();cpp,cchanges=cpp_source();obj=module(text);obj.prepare()
    (OUT/'candidate.py').write_text(text);(OUT/'effective-source.cpp').write_text(cpp)
    p=json.loads((OUT/'plan.json').read_text());p.update(checkpoint='1fa59fc',driver_substitutions=changes,native_substitutions=cchanges,
        correction='MPFR formats exact zero with trailing whitespace, which its direct string constructor rejected during the first boundary self-check. Trim the serialized scalar token before the exact round-trip check. Also remove the unjustified exact-dyadic assertion for z=(Q-c)/h: a quotient of dyadics can be rational. Keep exact Q/panel coordinates and enclose z outward. Rewrite near-endpoint logarithmic moments with the continuous x*log(abs(x)) function so intervals crossing +/-1 remain valid. Scientific tolerances, states, H formula, quadrature schedules and reference centers are unchanged.')
    files=[original.ROOT/'verification/gr_response_product_defined.py',OUT/'candidate.py',OUT/'effective-source.cpp',original.OUT/'plan.json',original.OUT/'h-0000-native.log']
    p['bindings'].update({x.relative_to(original.ROOT).as_posix():original.sha(x) for x in files});obj.save('plan.json',p)
    z=sp.symbols('z')
    for n in range(1,33):assert sp.expand(n*(n+1)*(sp.legendre(n+1,z)-sp.legendre(n-1,z))-(2*n+1)*(z*z-1)*sp.diff(sp.legendre(n,z),z))==0
    obj.save('stability-proof.json',dict(classification='Proven',passed=True,
        continuous_form='Let B=((z-1)*[(z+1)log|z+1|]-(z+1)*[(z-1)log|z-1|])/2=(z^2-1)Q0. The exact identity P_(n+1)-P_(n-1)=(2n+1)(z^2-1)Pn_prime/[n(n+1)] rewrites Jn without any division by z+/-1. Verified as exact polynomials through n=32.',
        interval_extension='x*log|x| is continuous at zero. Its only non-endpoint extrema are at +/-exp(-1), with values -/+exp(-1). Taking the hull of endpoints, zero when included, and any overlapping critical-point interval gives an outward range, including intervals crossing zero.',
        coordinate='Q and panel endpoints are exact binary inputs, but z=(Q-c)/h need not be dyadic. The interval formulas enclose that rational coordinate directly. An exact +/-1 shortcut is used only for a singleton endpoint interval.'))


def run(name):
    text,_=source();cpp,_=cpp_source();assert text==(OUT/'candidate.py').read_text() and cpp==(OUT/'effective-source.cpp').read_text();obj=module(text);obj.bindings();getattr(obj,name)()


def moment_controls():
    run('moment_controls');rows=[json.loads(x) for x in (OUT/'native-moments.jsonl').read_text().splitlines()];intervals=[x for x in rows if x['z'].startswith('interval')]
    assert len(intervals)==2 and all(x['passed'] for x in intervals)
    (OUT/'interval-controls.json').write_text(json.dumps(dict(classification='Counterexample candidate',passed=True,endpoint_intervals=2,moment_components=66),indent=2)+'\n')


def verify():
    run('verify');assert json.loads((OUT/'stability-proof.json').read_text())['passed'] and json.loads((OUT/'interval-controls.json').read_text())['passed']


if __name__=='__main__':
    if sys.argv[1]=='prepare':prepare()
    elif sys.argv[1]=='moment_controls':moment_controls()
    elif sys.argv[1]=='verify':verify()
    else:run(sys.argv[1])
