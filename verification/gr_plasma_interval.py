"""Bind native interval quadrature to frozen actual plasma states and root residuals."""
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction as F
import json, shlex, subprocess, sys
import numpy as np
from mpmath import iv
import verify_plasma_moment_certificate as audit
from interval_records import interval_text, exact_endpoint

g=audit.g;OUT=g.OUT/'gr-plasma-interval';CACHE=g.CACHE/'plasma-interval'
CAPD=g.CACHE.parent/'nutimo_pilot/request14_capd'
DEPS=g.CACHE.parent/'nutimo_pilot/request15_deps/usr'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')
def upper(x): return exact_endpoint(x._mpi_[1])
def lower(x): return exact_endpoint(x._mpi_[0])
def I(x):
    if isinstance(x,F):return iv.mpf(x.numerator)/x.denominator
    return iv.mpf(str(x))


def prepare():
    assert not OUT.exists() and not CACHE.exists();OUT.mkdir();CACHE.mkdir();audit.verify()
    data=np.load(audit.certificate.dispersion.plasma.OUT/'stellar-fermi-plasma.npz')
    (OUT/'states.tsv').write_text(''.join(f'{i} {float(e).hex()} {float(b).hex()}\n' for i,(e,b) in enumerate(zip(data['eta'],data['beta']))))
    paths=[g.ROOT/'verification/gr_plasma_interval.py',g.ROOT/'verification/gr_plasma_interval.cpp',audit.OUT/'manifest.json',OUT/'states.tsv']
    # Reuse the installed interval arithmetic, never its old application runtime.
    paths += [CAPD/'build-request15-mp/libcapd.a',CAPD/'build-request15-mp/bin/capd-config',DEPS/'include/mpfr.h']
    headers={str(p):g.c.sha(p) for directory in [CAPD/'capdAlg/include/capd/intervals',CAPD/'capdAlg/include/capd/multiPrec'] for p in directory.rglob('*') if p.is_file()}
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='c4ff34b',
        bindings={str(p):g.c.sha(p) for p in paths},interval_headers=headers,
        panels=4096,bits=128,controls=[0,1175,2972,3043,4352,5734],cells=5735,block_size=128,processes=4,
        actual_inputs='Exact saved binary64 eta and beta, transferred in hexadecimal. Not an enclosure of physical or native EOS parameter error.',
        method='Existing CAPD MPFR interval operations; generated Gauss nodes/weights verified through degree seven. Evaluate the new declared electron rule and add its continuous moment remainder plus omitted positrons. Endpoint decimal output uses directed MPFR rounding at 60 digits.',
        root_error_tolerance='1e-10',adaptive_moment_normalized_error_tolerance='1e-10',
        control='At six actual cells validate native interval moment enclosures against prior adaptive moments at relative tolerance 1e-10 with floor one, and certify all registered frozen transverse/longitudinal candidate roots at absolute h error 1e-10. If controls pass, execute all 5735 cells with four independent block workers.',
        scope='Pointwise certified leading-order model roots at exact saved eta,beta,q inputs. No physical parameter uncertainty, full state derivative certificate, higher-order plasma effects or native EOS/GR replacement.'))


def bindings():
    plan=json.loads((OUT/'plan.json').read_text())
    for key in ['bindings','interval_headers']:
        for name,digest in plan[key].items():
            from pathlib import Path
            assert g.c.sha(Path(name))==digest,name
    return plan


def build():
    bindings();flags=shlex.split(subprocess.check_output([str(CAPD/'build-request15-mp/bin/capd-config'),'--cflags','--libs'],text=True))
    cmd=['g++',str(g.ROOT/'verification/gr_plasma_interval.cpp'),f'-I{DEPS}/include',f'-I{DEPS}/include/x86_64-linux-gnu',
        f'-L{DEPS}/lib/x86_64-linux-gnu',*flags,'-o',str(CACHE/'moments')]
    p=subprocess.run(cmd,capture_output=True,text=True);(OUT/'build.log').write_text(p.stdout+p.stderr)
    save('build.json',dict(classification='Counterexample candidate',command=cmd,returncode=p.returncode))
    assert p.returncode==0,p.stderr
    libraries=subprocess.check_output(['ldd',str(CACHE/'moments')],text=True);(OUT/'linked-libraries.txt').write_text(libraries)
    from pathlib import Path
    libs={}
    for line in libraries.splitlines():
        for word in line.split():
            if word.startswith('/') and Path(word).is_file():libs[word]=g.c.sha(Path(word))
    save('runtime.json',dict(binary=str(CACHE/'moments'),binary_sha256=g.c.sha(CACHE/'moments'),libraries=libs))


def evaluate(first,count,label):
    plan=json.loads((OUT/'plan.json').read_text());runtime=json.loads((OUT/'runtime.json').read_text())
    assert g.c.sha(CACHE/'moments')==runtime['binary_sha256']
    target=OUT/(label+'.jsonl');assert not target.exists()
    cmd=[str(CACHE/'moments'),str(OUT/'states.tsv'),str(target),str(plan['panels']),str(first),str(count)]
    p=subprocess.run(cmd,capture_output=True,text=True);(OUT/(label+'.log')).write_text(p.stdout+p.stderr)
    assert p.returncode==0,p.stderr
    return label


def analyze(labels):
    plan=json.loads((OUT/'plan.json').read_text());iv.dps=60
    cert=json.loads((audit.certificate.OUT/'uniform-certificate.json').read_text())
    tails=json.loads((audit.certificate.OUT/'result.json').read_text())
    errors=[]
    for j in range(26):
        row=next(r for r in cert['rows'] if r['component']==j and r['eta_order']==r['beta_order']==0)
        omitted=I(row['electron_tail'])+I(row['omitted_positron'])
        errors.append(((I(row['panels'])/plan['panels'])**8*(I(row['error'])-omitted)+omitted).b)
    frozen=dict(np.load(audit.certificate.dispersion.OUT/'stellar-dispersion.npz'))
    states=dict(np.load(audit.certificate.dispersion.plasma.OUT/'stellar-fermi-plasma.npz'))
    rows=[];max_error=F(0);max_width=F(0);total_seconds=0.
    for label in labels:
        for line in (OUT/(label+'.jsonl')).read_text().splitlines():
            raw=json.loads(line);i=raw['cell'];assert raw['panels']==plan['panels'] and raw['bits']==plan['bits']
            total_seconds+=raw['seconds'];J=[]
            for value,error in zip(raw['moments'],errors):
                # Parse exported directed decimal endpoints before adding the remainder.
                a,b=map(F,value);assert a<=b;max_width=max(max_width,b-a)
                J.append(iv.mpf([I(a).a,I(b).b])+iv.mpf([-error,error]))
            denominator=J[0]-J[1]/3;assert lower(denominator)>0
            beta=I(F.from_float(float(states['beta'][i])));eta=I(F.from_float(float(states['eta'][i])))
            scale=I(F.from_float(float(states['dimensionless_density_scale'][i])))
            factor=scale/(iv.sqrt(I(2))*beta*iv.sqrt(beta))*(iv.exp(-eta) if float(states['eta'][i])<0 else I(1))
            moment_error=F(0)
            for j in range(26):
                previous=I(F.from_float(float(frozen['raw_moments'][j,i])))*factor
                error=upper(abs(previous-J[j]))/max(F(1),upper(abs(previous)))
                moment_error=max(moment_error,error)
            co=[I(1)]+[(J[j]/(2*j+1)-J[j+1]/(2*j+3))/denominator for j in range(1,25)]
            tests=[]
            for kind,count in [('T',6),('L',5)]:
                coefficients=co if kind=='T' else [(2*j+1)*c for j,c in enumerate(co)]
                endpoint=sum(coefficients);tail=I(tails['transverse_uniform_tail' if kind=='T' else 'longitudinal_uniform_tail']).b
                for k in range(count):
                    h=I(F.from_float(float(frozen[f'{kind}-{k}-h'][i])))
                    q=I(F.from_float(float(frozen[f'{kind}-{k}-q'][i])))
                    z=q*q/(q*q+h) if kind=='T' else q*q/h
                    assert lower(z)>=0 and upper(z)<=1
                    if kind=='L':assert upper(q*q)<lower(endpoint), 'timelike branch not certified'
                    value=I(0)
                    for c in coefficients[::-1]:value=value*z+c
                    residual=h-value-iv.mpf([I(0),(tail*z**25).b])
                    error=upper(abs(residual));max_error=max(max_error,error)
                    tests.append(dict(kind=kind,index=k,residual=interval_text(residual),h_error_upper_rational=str(error),
                        passed=error<F(plan['root_error_tolerance'])))
            rows.append(dict(cell=i,tests=tests,adaptive_moment_error_upper_rational=str(moment_error),
                adaptive_moments_passed=moment_error<F(plan['adaptive_moment_normalized_error_tolerance'])))
    assert len({r['cell'] for r in rows})==len(rows)
    return dict(classification='Proven',rows=rows,cells=len(rows),passed=all(r['adaptive_moments_passed'] and all(t['passed'] for t in r['tests']) for r in rows),
        maximum_h_error_upper_rational=str(max_error),maximum_raw_interval_width_rational=str(max_width),
        display_only=dict(maximum_h_error=float(max_error),maximum_raw_interval_width=float(max_width),total_compute_seconds=total_seconds),
        exact_saved_parameter_inputs=True,physical_EOS_certified=False,full_GR_evolution=False)


def controls():
    plan=bindings();labels=[]
    for cell in plan['controls']:labels.append(evaluate(cell,1,f'control-{cell}'))
    result=analyze(labels);save('controls.json',result);assert result['passed']
    save('controls-manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    print('PASS actual plasma interval controls',result['cells'],result['display_only'],flush=True)


def run():
    plan=bindings();assert json.loads((OUT/'controls.json').read_text())['passed']
    labels=[]
    with ProcessPoolExecutor(max_workers=plan['processes']) as pool:
        futures=[pool.submit(evaluate,i,min(plan['block_size'],plan['cells']-i),f'block-{i:04}') for i in range(0,plan['cells'],plan['block_size'])]
        for future in as_completed(futures):
            labels.append(future.result());save('progress.json',dict(completed_blocks=labels,total_blocks=len(futures)))
            print('PLASMA INTERVAL BLOCK',len(labels),'/',len(futures),flush=True)
    result=analyze(sorted(labels));save('result.json',result);assert result['cells']==plan['cells'] and result['passed']
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'result.json').read_text())['passed']
    print('PASS all saved plasma-state interval root certificates; physical/GR closure remains open',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
