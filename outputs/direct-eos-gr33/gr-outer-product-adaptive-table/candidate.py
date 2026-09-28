"""Bind reusable H responses to exact outer nodes and true-neutral field bounds."""
from fractions import Fraction as F
from pathlib import Path
import gzip,json,math,subprocess,sys,time
import numpy as np
from mpmath import iv
import gr_response_product_128 as product
import gr_outer_self_pilot as outer
import gr_polarization_neutral_refinement as neutral_shift

g=product.g;ROOT=g.ROOT;OUT=product.OUT.parent/'gr-outer-product-pilot';CACHE=g.native.g.CACHE/'outer-product-pilot'
I=outer.I;low=outer.low;high=outer.high;cusp=outer.cusp;sha=outer.sha


def save(name,value):(OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    product.verify();outer.domain.verify();neutral_shift.verify();assert not OUT.exists() and not CACHE.exists();OUT.mkdir();CACHE.mkdir()
    cpp=(product.OUT/'effective-source.cpp').read_text();old='MI Q=binary(qtext);Six sum{};';new='MI Q=encoded(qtext,qtext);exact(Q);Six sum{};'
    assert cpp.count(old)==1;(OUT/'effective-source.cpp').write_text(cpp.replace(old,new))
    files=[ROOT/'verification/gr_outer_product_pilot.py',ROOT/'verification/gr_outer_self_pilot.py',OUT/'effective-source.cpp',
        product.OUT/'manifest.json',product.OUT/'result.json',product.OUT/'inputs.tsv',product.OUT/'rule.tsv',
        outer.domain.OUT/'manifest.json',outer.domain.OUT/'result.json',outer.SMOOTH/'rule.json',neutral_shift.OUT/'manifest.json',neutral_shift.OUT/'result.json',
        outer.highq.OUT/'manifest.json',outer.infrared.OUT/'manifest.json',cusp.ROOTS/'result.json',cusp.THERMO/'states.npz',
        outer.highq.moments.OUT/'inputs.npz',outer.highq.neutral.density.ionic.OUT/'inputs.npz',outer.highq.neutral.density.ionic.OUT/'constants.json']
    files+=list(product.OUT.glob('h-*.tsv'))
    save('plan.json',dict(classification='Proven',checkpoint='d56386d',positions=[0,1603,3205],bits=256,native_bits=128,node_binary_bits=80,
        preflight_nodes=16,response_midpoint_budget='1e-15',physical_field_budget='2e-7',
        bindings={x.relative_to(ROOT).as_posix():sha(x) for x in files},interval_headers=json.loads((product.OUT/'plan.json').read_text())['interval_headers'],
        source_change={old:new},
        scope='Three whole outer-integral pilots with stored eta/neutral coefficients, separately enclosed displacement to the true neutral model, and certified IR/high-Q endpoints. Retain the original direct-momentum outer run as an independent path. No whole-table or physical EOS/GR claim.',
        node_binding='Reuse the original exact dyadic panel schedules and 80-bit dyadic Q rounding. The native evaluator parses Q as an exact decimal dyadic, rejects any non-singleton interval and never rounds it to binary64.',
        center_domain='For any true root in the certified eta interval let delta bound its distance from the stored eta. The original eta-polydisk still contains the disk about the stored center of radius h_eta-delta>0. Apply Cauchy there with unchanged tau/Q domains. Recompute the six derivative majorants and old-neutral-matrix physical envelope; rescale the original outer remainder by the exact ratio of new to original maximum physical envelopes.',
        neutral_correction='The earlier center/coefficient error proof integrates nonnegative all-Q envelopes. The same whole-integral upper bound therefore bounds their integral over this fixed middle subinterval. Add that bound to the stored-center/stored-coefficient middle fields only; IR and high-Q fields already use true roots. Endpoints, beta, scale, B and inventories are held fixed in differentiation.',
        remainder='No new outer quadrature tolerance. Actual point intervals, exact-node displacement, rescaled analytic outer remainder, neutral displacement and both endpoints enter the final field intervals. Keep the original 2e-7 field gate.'))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert sha(ROOT/rel)==digest,rel
    for path,digest in p['interval_headers'].items():assert sha(Path(path))==digest,path
    return p


def engine():
    text,_,_=product.sources();obj=product.module(text);obj.OUT=OUT;obj.CACHE=CACHE;obj.bindings=bindings;return obj


def build():engine().build()


def evaluate(position,limit=None):
    p=bindings();binary=engine().runtime();began=time.perf_counter();state,ions,roots,data,_,_=outer.highq.context();iv.prec=p['bits']
    meta=json.loads((outer.domain.OUT/'result.json').read_text());row=meta['rows'][position];schedule=full_schedule(position,row,state,data)
    selected=outer.nodes(schedule,row,p);selected=selected if limit is None else selected[:limit];prefix=('cell' if limit is None else 'preflight')+f'-{position:04d}'
    qpath=OUT/f'{prefix}-q.tsv';assert not qpath.exists();qpath.write_text(''.join(f"{j} {x['Q_decimal']}\n" for j,x in enumerate(selected)))
    save(prefix+'-nodes.json',dict(nodes=selected))
    target=OUT/f'{prefix}-native.jsonl';assert not target.exists()
    cmd=[binary,'response',str(FULL_INPUTS),str(product.OUT/'rule.tsv'),str(target),str(position),str(CURRENT_H),str(qpath)]
    done=subprocess.run(cmd,capture_output=True,text=True);(OUT/f'{prefix}-native.log').write_text(done.stdout+done.stderr);assert done.returncode==0,done.stderr
    root=roots[position];center=F.from_float(float(state['eta'][position]));a,b=cusp.endpoints(root['root']);delta=max(abs(a-center),abs(b-center))
    dp=json.loads((outer.domain.OUT/'plan.json').read_text());he=F(dp['eta_radius'])-delta;ht=F(dp['tau_radius']);assert he>0
    majorants=[I(2*math.factorial(n)*math.factorial(m))/(I(he)**n*I(ht)**m) for n,m in outer.domain.ORDERS]
    beta=I(data['beta'][position]);Sref=I(data['Sref'][position]);scale=I(state['scale'][position]);const=json.loads((outer.highq.neutral.density.ionic.OUT/'constants.json').read_text())['native_binary64'];B=I(const['alpha'])/iv.pi
    coeff=[I(x) for x in state['neutral_coefficients'][:,position]];matrix=outer.highq.neutral.matrix(coeff);counts=[I(x)/I(a) for x,a in zip(ions['X'][position],ions['A'])]
    Z2=sum(x*I(z)**2 for x,z in zip(counts,ions['Z']))/sum(counts);pref=2*B*Z2*scale/beta
    envelopes=[pref*sum(abs(c)*m for c,m in zip(line,majorants))/scale for line in matrix]
    ratio=I(max(map(high,envelopes)))/I(max(map(F,row['physical_integrand_majorants'])));remainder=high(I(schedule['field_remainder_upper'])*ratio)
    tail=json.loads((g.rule.OUT/'tail.json').read_text())['rows'][position];raw=[I(0)]*6;maximum=F(0);point_path=OUT/f'{prefix}-points.jsonl.gz'
    assert not point_path.exists()
    with target.open() as native,gzip.open(point_path,'xt',compresslevel=1) as stream:
        for j,(line,node) in enumerate(zip(native,selected,strict=True)):
            rec=json.loads(line);assert rec['position']==position and rec['cell']==root['cell'] and rec['z_index']==j;Q=I(node['Q']);R=[]
            for (a,b),t in zip(rec['finite_response'],tail['upper_coefficients'],strict=True):
                value=iv.mpf([I(a).a,I(b).b])+cusp.symmetric(high(I(t)/I(max(F(1),low(Q*Q)))))
                left,right=low(value),high(value);mid=F.from_float(float((left+right)/2));maximum=max(maximum,abs(left-mid),abs(right-mid));R.append(value)
            assert low(R[0])>0;X=Q*Q/(B*Sref);den=X+R[0]
            G=[R[0]/den,X*R[1]/den**2,X*(R[2]/den**2-2*R[1]**2/den**3),X*R[3]/den**2,
                X*(R[4]/den**2-2*R[1]*R[3]/den**3),X*(R[5]/den**2-2*R[3]**2/den**3)]
            G=[x+cusp.symmetric(high(I(node['shift_factor'])*m)) for x,m in zip(G,majorants)]
            raw=[x+cusp.interval(node['weight'])*y/scale for x,y in zip(raw,G)]
            stream.write(json.dumps(dict(node=j,Q_exact=node['Q'],response=list(map(cusp.interval_text,R)),G=list(map(cusp.interval_text,G))))+'\n')
    result=dict(classification='Proven',position=position,cell=root['cell'],nodes=len(selected),complete_finite_interval=limit is None,
        response_midpoint_error_upper=str(maximum),point_budget_passed=maximum<F(p['response_midpoint_budget']),
        partial_raw_integrals=list(map(cusp.interval_text,raw)),eta_center_distance_upper=str(delta),eta_Cauchy_radius=str(he),derivative_majorants=list(map(cusp.interval_text,majorants)),
        outer_field_remainder_upper=str(remainder),native_seconds=json.loads(done.stdout)['seconds'],elapsed_seconds=time.perf_counter()-began,
        physical_EOS_certified=False)
    if limit is None:
        correction=json.loads((neutral_shift.OUT/'result.json').read_text())['records'][position];assert correction['cell']==root['cell']
        middle=[-pref*sum(c*x for c,x in zip(line,raw))+cusp.symmetric(remainder+high(cusp.interval(correction['field_error_upper'][field]))) for line,field in zip(matrix,outer.highq.FIELDS)]
        ir=json.loads((outer.infrared.OUT/'result.json').read_text())['rows'][position];hq=next(x for x in outer.highq.records() if x['position']==position)
        total=[x+cusp.interval(a)+cusp.interval(b) for x,a,b in zip(middle,ir['field_enclosures'],hq['field_enclosures'])];texts=list(map(cusp.interval_text,total));points=[];errors=[]
        for text in texts:
            a,b=cusp.endpoints(text);point=float((a+b)/2);points.append(point);errors.append(max(abs(a-F.from_float(point)),abs(b-F.from_float(point))))
        difference=float(np.max(abs(np.array(points)-state['fine_fields'][:,position])))
        result.update(passed=result['point_budget_passed'] and max(errors)<F(p['physical_field_budget']),middle_field_enclosures=list(map(cusp.interval_text,middle)),
            complete_field_enclosures=texts,field_approximation=points,field_error_upper=list(map(str,errors)),
            finite_reference_control=dict(classification='Counterexample candidate',maximum_absolute_difference=difference,passed=difference<float(p['physical_field_budget']),
                scope='Frozen numerical table at its stored eta/coefficients and binary64 constants. Its adaptive convergence flags were not checked, so this is only a finite comparison.'))
    save(prefix+'-result.json',result)
    # The native plain output is retained while the pilot is active; do not rewrite the original direct run.
    print('PRODUCT OUTER',position,len(selected),float(maximum),result.get('passed'),result['native_seconds'],flush=True);return result


def preflight():
    p=bindings();results=[evaluate(i,p['preflight_nodes']) for i in p['positions']];assert all(x['point_budget_passed'] for x in results)
    save('preflight-manifest.json',dict(sha256={x.relative_to(ROOT).as_posix():sha(x) for x in OUT.iterdir() if x.is_file()}));verify_preflight()


def verify_preflight():
    p=bindings();engine().runtime()
    for rel,digest in json.loads((OUT/'preflight-manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    for i in p['positions']:
        r=json.loads((OUT/f'preflight-{i:04d}-result.json').read_text());assert r['point_budget_passed'] and r['nodes']==16 and not r['complete_finite_interval']
    print('PASS 48 actual 80-bit Q product preflight nodes, shared-domain center shift and neutral assembly definitions',flush=True)


def run():
    verify_preflight()
    for i in bindings()['positions']:evaluate(i)


def finalize():
    p=bindings();engine().runtime()
    for i in p['positions']:
        r=json.loads((OUT/f'cell-{i:04d}-result.json').read_text());assert r['passed'] and r['finite_reference_control']['passed']
    save('manifest.json',dict(sha256={x.relative_to(ROOT).as_posix():sha(x) for x in OUT.iterdir() if x.is_file()}));verify()


def verify():
    verify_preflight();p=bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    for i in p['positions']:
        r=json.loads((OUT/f'cell-{i:04d}-result.json').read_text());assert r['passed'] and r['complete_finite_interval'] and r['finite_reference_control']['passed']
        for text,point,error in zip(r['complete_field_enclosures'],r['field_approximation'],r['field_error_upper'],strict=True):
            a,b=cusp.endpoints(text);assert max(abs(a-F.from_float(point)),abs(b-F.from_float(point)))==F(error)<F(p['physical_field_budget'])
    print('PASS three product self-integral pilots including center/neutral, quadrature, IR and infinite high-Q errors; full EOS/GR remain open',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
