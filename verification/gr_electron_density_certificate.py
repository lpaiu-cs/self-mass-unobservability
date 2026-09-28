"""Neutral ideal-electron roots and implicit derivatives at actual ion inventories."""
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction as F
from pathlib import Path
import json, shlex, subprocess, sys
import numpy as np
import mpmath as mp
from mpmath import iv
import sympy as sp
from scipy.integrate import quad_vec
from scipy.special import expit
from interval_records import exact_endpoint, interval_text
import gr_polarization_self_scaled as self_energy
import gr_ionic_hamiltonian_bound as ionic
import gr_plasma_interval as native
import verify_fermi_uniform as uniform

g=ionic.g;OUT=g.OUT/'gr-electron-density-certificate';CACHE=g.CACHE/'electron-density-certificate'
ORDER=[(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)]
FIELDS=['eta_lnn','eta_lnT','eta_lnn_lnn','eta_lnn_lnT','eta_lnT_lnT']


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2,default=self_energy.original.finite.previous.original.previous.scalar)+'\n')
def low(x):return exact_endpoint(x._mpi_[0])
def high(x):return exact_endpoint(x._mpi_[1])
def I(x):
    if isinstance(x,F):return iv.mpf(x.numerator)/x.denominator
    if isinstance(x,(float,np.floating)):return iv.mpf(float(x))
    return iv.mpf(x)
def read_interval(pair):
    a,b=map(F,pair);assert a<=b
    return iv.mpf([I(a).a,I(b).b])


def prepare():
    assert not OUT.exists() and not CACHE.exists();OUT.mkdir();CACHE.mkdir()
    self_energy.verify();uniform.verify();p=json.loads((native.OUT/'plan.json').read_text())
    paths=[g.ROOT/'verification/gr_electron_density_certificate.py',g.ROOT/'verification/gr_electron_density_certificate.cpp',
        g.ROOT/'verification/interval_records.py',self_energy.OUT/'manifest.json',
        self_energy.original.finite.OUT/'states.npz',ionic.OUT/'inputs.npz',ionic.OUT/'constants.json',
        uniform.OUT/'manifest.json',uniform.OUT/'uniform-certificate.json',native.OUT/'plan.json']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='2565748',
        bindings={str(x):g.c.sha(x) for x in paths},interval_headers=p['interval_headers'],
        arithmetic_sources={str(x):g.c.sha(x) for x in Path(mp.__file__).parent.rglob('*.py')},
        panels=8192,bits=128,root_half_width_power=-40,processes=3,block_size=128,cells=3206,
        control_positions=[0,801,1603,2404,3205],candidate_Newton_steps=3,candidate_tolerance=2e-13,
        control_momentum_digits=60,finite_comparison_gate='1e-7',fields=FIELDS,
        roots='Solve N(eta,beta)=N_A*rho*sum_j X_j*Z_j/A_j for each actual saved positive-ion inventory, treating the binary rho,X,A,Z and constants as exact model parameters. N is the ideal relativistic electron-only number density. The assumption that each isotope is fully ionized is explicit; this does not certify actual ionization or interacting chemical potentials.',
        exact_rule='Use the earlier uniformly certified composite four-point Gauss rule after t=s^2, s in [0,10], at 8192 panels. Evaluate D and its first/second eta,beta derivatives with CAPD/MPFR outward arithmetic. Add the earlier continuous quadrature plus infinite-tail error, also uniformly over a center +/- 2^-40 eta box.',
        root_certificate='Enclose D(center)-target, certify a strictly positive D_eta throughout the box, and prove residual_upper < slope_lower*2^-40. Mean value theorem gives existence, uniqueness and an interval Newton enclosure of the neutral root. Enclose both first and all second derivatives with respect to ln(number) and ln(T) at this root.',
        strict_scope='Complete root and derivative enclosures for these declared ideal-electron/fully-ionized exact input states. Report widths; no target physical EOS tolerance is inferred. Floating predictors and old eta values are not treated as certified roots.',
        physical_EOS_certified=False,native_EOS_replaced=False))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for group in ['bindings','interval_headers','arithmetic_sources']:
        for name,digest in p[group].items():assert g.c.sha(Path(name))==digest,name
    native.bindings()
    return p


def momentum(eta,beta,tolerance):
    """Six derivatives of N, divided by a fixed sqrt(2)*beta^(3/2)."""
    def f(w):
        p2=2*beta*w*w;gamma=np.sqrt(1+p2);t=2*w*w/(gamma+1);q=expit(eta-t);v=expit(t-eta)
        base=2*w*w*q
        return base*np.array([np.ones_like(q),v,v*(1-2*q),t*v,t*v*(1-2*q),(t*t*(1-2*q)-t)*v])
    return quad_vec(f,0,np.inf,epsabs=tolerance,epsrel=tolerance,norm='max',limit=1200)[0]


def candidates():
    p=bindings();iv.prec=128
    star=dict(np.load(ionic.OUT/'inputs.npz'));a=dict(np.load(self_energy.original.finite.OUT/'states.npz'))
    assert np.array_equal(star['cells'],a['cells']);c=json.loads((ionic.OUT/'constants.json').read_text())['native_binary64']
    targets=[];densities=[]
    for i in range(len(a['cells'])):
        n=I(c['N_A'])*I(star['rho'][i])*sum(I(x)*I(z)/I(m) for x,z,m in zip(star['X'][i],star['Z'],star['A']))
        beta=I(a['beta'][i]);target=n/(I(c['number_prefactor_cm3'])*iv.sqrt(2)*beta*iv.sqrt(beta))
        targets.append(dict(position=i,cell=int(a['cells'][i]),target=interval_text(target),number_density_cm3=interval_text(n)))
        densities.append(float((low(target)+high(target))/2))
    target=np.array(densities);eta=a['eta'].copy();steps=[]
    for _ in range(p['candidate_Newton_steps']):
        values=momentum(eta,a['beta'],p['candidate_tolerance']);delta=(values[0]-target)/values[1];eta-=delta
        steps.append(float(np.max(abs(delta))))
    np.savez_compressed(OUT/'candidates.npz',cells=a['cells'],eta_original=a['eta'],eta_center=eta,beta=a['beta'],target=target,
        momentum=momentum(eta,a['beta'],p['candidate_tolerance']))
    save('targets.json',dict(classification='Proven',input_scope='Exact binary inventory and constants with outward interval evaluation; full ionization is an explicit model assumption.',rows=targets))
    (OUT/'states.tsv').write_text(''.join(f'{i} {int(cell)} {float(e).hex()} {float(b).hex()}\n' for i,(cell,e,b) in enumerate(zip(a['cells'],eta,a['beta']))))
    save('candidate-result.json',dict(classification='Counterexample candidate',maximum_Newton_steps=steps,
        maximum_eta_shift=float(np.max(abs(eta-a['eta']))),certified_roots=False))
    save('candidate-manifest.json',dict(sha256={x.relative_to(g.ROOT).as_posix():g.c.sha(x) for x in OUT.iterdir() if x.name in ['candidates.npz','targets.json','states.tsv','candidate-result.json']}))


def build():
    bindings();flags=shlex.split(subprocess.check_output([str(native.CAPD/'build-request15-mp/bin/capd-config'),'--cflags','--libs'],text=True))
    cmd=['g++',str(g.ROOT/'verification/gr_electron_density_certificate.cpp'),f'-I{native.DEPS}/include',f'-I{native.DEPS}/include/x86_64-linux-gnu',f'-L{native.DEPS}/lib/x86_64-linux-gnu',*flags,'-o',str(CACHE/'density')]
    done=subprocess.run(cmd,capture_output=True,text=True);(OUT/'build.log').write_text(done.stdout+done.stderr)
    save('build.json',dict(command=cmd,returncode=done.returncode));assert done.returncode==0,done.stderr
    linked=subprocess.check_output(['ldd',str(CACHE/'density')],text=True);(OUT/'linked-libraries.txt').write_text(linked)
    libs={word:g.c.sha(Path(word)) for line in linked.splitlines() for word in line.split() if word.startswith('/') and Path(word).is_file()}
    save('runtime.json',dict(binary=str(CACHE/'density'),binary_sha256=g.c.sha(CACHE/'density'),libraries=libs))


def partial_errors(beta):
    cert=json.loads((uniform.OUT/'uniform-certificate.json').read_text());fields={}
    for row in cert['rows']:
        if row['nu'] not in ['0.5','1.5']:continue
        old=read_interval(row['uniform_error_upper'][1:-1].split(','));tail=read_interval(row['tail_upper'][1:-1].split(','))
        # Positive overestimate: scaling the entire old bound and adding its tail again.
        fields[(row['nu'],row['eta_order'],row['beta_order'])]=((I(row['panels'])/8192)**8*old+tail).b
    answer=[]
    for k,l in ORDER:
        error=fields[('0.5',k,l)]+beta*fields[('1.5',k,l)]
        if l:error+=l*fields[('1.5',k,l-1)]
        answer.append(error.b)
    return answer


def inverse_derivatives(target,values,beta):
    D,De,Db,Dee,Deb,Dbb=values
    Nt=I(3)/2*D+beta*Db;Net=I(3)/2*De+beta*Deb;Ntt=I(9)/4*D+4*beta*Db+beta**2*Dbb
    en=target/De;et=-Nt/De
    return [en,et,(target-Dee*en**2)/De,-en*(Dee*et+Net)/De,-(Ntt+2*Net*et+Dee*et**2)/De]


def symbolic():
    t,b,e=sp.symbols('t b e',positive=True);w=1+b*t/2;H=sp.sqrt(w)*(1+b*t)
    assert sp.simplify(sp.diff(H,b)-(t*sp.sqrt(w)+t*(1+b*t)/(4*sp.sqrt(w))))==0
    assert sp.simplify(sp.diff(H,b,2)-(t*t/(2*sp.sqrt(w))-t*t*(1+b*t)/(16*w**sp.Rational(3,2))))==0
    D=sp.Function('D')(e,b);physical=b**sp.Rational(3,2)*D
    first=b*sp.diff(physical,b);second=b*sp.diff(first,b)
    assert sp.simplify(first/b**sp.Rational(3,2)-(sp.Rational(3,2)*D+b*sp.diff(D,b)))==0
    assert sp.simplify(second/b**sp.Rational(3,2)-(sp.Rational(9,4)*D+4*b*sp.diff(D,b)+b*b*sp.diff(D,b,2)))==0
    N,E,EE,T,ET,TT=sp.symbols('N E EE T ET TT',nonzero=True)
    en=N/E;et=-T/E;enn=(N-EE*en**2)/E;ent=-en*(EE*et+ET)/E;ett=-(TT+2*ET*et+EE*et**2)/E
    for value in [E*en-N,E*et+T,EE*en**2+E*enn-N,EE*en*et+ET*en+E*ent,TT+2*ET*et+EE*et**2+E*ett]:assert sp.simplify(value)==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        density='N/n_pref=sqrt(2)*beta^(3/2)*D. D=F_(1/2)+beta*F_(3/2). Exact beta derivatives include the product-rule terms.',
        thermal='After dividing each derivative by the same fixed sqrt(2)*beta^(3/2), N_tau=(3/2)D+beta*D_beta, N_etatau=(3/2)D_eta+beta*D_etabeta, N_tautau=(9/4)D+4*beta*D_beta+beta^2*D_betabeta, tau=ln(T).',
        implicit='eta_n=N/N_eta, eta_tau=-N_tau/N_eta; eta_nn=(N-N_etaeta*eta_n^2)/N_eta; eta_ntau=-eta_n*(N_etaeta*eta_tau+N_etatau)/N_eta; eta_tautau=-(N_tautau+2*N_etatau*eta_tau+N_etaeta*eta_tau^2)/N_eta; n denotes ln(number density).',
        root_proof='A positive lower bound m on D_eta throughout [c-delta,c+delta], and |D(c)-target|<m*delta, imply a unique root inside this bracket by the mean value theorem. c+(target-D(c))/D_eta(box) encloses that root. All derivatives are enclosed over the original bracket, so they enclose their values at the neutral root.',
        limitations='The exact input inventory is fully ionized by assumption. No nonideal chemical potential, full EOS or GR closure is inferred.'))


def independent_momentum(eta,beta):
    eta,beta=mp.mpf(float(eta)),mp.mpf(float(beta));normal=mp.sqrt(2)*beta**mp.mpf('1.5')
    top=max(mp.mpf(2),eta);points=[mp.mpf(0)]+[mp.sqrt(beta*t*(2+beta*t)) for t in [1,top,top+8,top+40]]+[mp.inf]
    def f(p,j):
        gamma=mp.sqrt(1+p*p);t=p*p/(beta*(gamma+1));q=1/(1+mp.exp(t-eta));v=1-q
        factors=[1,v,v*(1-2*q),t*v,t*v*(1-2*q),(t*t*(1-2*q)-t)*v]
        return p*p*q*factors[j]/normal
    return np.array([float(mp.quad(lambda p:f(p,j),points)) for j in range(6)])


def evaluate(first,count,label):
    runtime=json.loads((OUT/'runtime.json').read_text());assert g.c.sha(CACHE/'density')==runtime['binary_sha256']
    target=OUT/(label+'.jsonl');assert not target.exists()
    done=subprocess.run([str(CACHE/'density'),str(OUT/'states.tsv'),str(target),'8192',str(first),str(count)],capture_output=True,text=True)
    (OUT/(label+'.log')).write_text(done.stdout+done.stderr);assert done.returncode==0,done.stderr
    return label


def analyze(labels):
    p=json.loads((OUT/'plan.json').read_text());iv.prec=128
    a=dict(np.load(OUT/'candidates.npz'));targets=json.loads((OUT/'targets.json').read_text())['rows'];rows=[];widths=[];scores=[]
    radius=I(1)/2**40
    for label in labels:
        for line in (OUT/(label+'.jsonl')).read_text().splitlines():
            raw=json.loads(line);i=raw['position'];assert raw['cell']==int(a['cells'][i]) and raw['panels']==8192 and raw['bits']==128
            beta=I(a['beta'][i]);errors=partial_errors(beta)
            values=[read_interval(pair)+iv.mpf([-e,e]) for pair,e in zip(raw['box_D_partials'],errors)]
            center=read_interval(raw['center_D'])+iv.mpf([-errors[0],errors[0]])
            target=read_interval(targets[i]['target'][1:-1].split(','));slope=low(values[1]);assert slope>0
            residual=high(abs(center-target));bracket_pass=residual<slope*low(radius)
            root=I(a['eta_center'][i])+(target-center)/values[1]
            assert low(root)>=F.from_float(float(a['eta_center'][i]))-low(radius) and high(root)<=F.from_float(float(a['eta_center'][i]))+low(radius)
            derivatives=inverse_derivatives(target,values,beta)
            m=a['momentum'][:,i];ref=[m[0]/m[1],-m[3]/m[1]]
            ref += [(m[0]-m[2]*ref[0]**2)/m[1],-ref[0]*(m[2]*ref[1]+m[4])/m[1],-(m[5]+2*m[4]*ref[1]+m[2]*ref[1]**2)/m[1]]
            score=max(high(abs(v-I(float(r))))/max(F(1),F.from_float(abs(float(r)))) for v,r in zip(derivatives,ref))
            scores.append(score);width=high(root)-low(root);widths.append(width)
            rows.append(dict(position=i,cell=raw['cell'],root=interval_text(root),density_slope_lower=str(slope),
                root_error_from_center_upper=str(residual/slope),root_width=str(width),bracket_passed=bracket_pass,
                derivatives={name:interval_text(v) for name,v in zip(FIELDS,derivatives)},
                finite_momentum_comparison_upper=str(score),finite_passed=score<F(p['finite_comparison_gate'])))
    assert len({r['position'] for r in rows})==len(rows)
    return dict(classification='Proven',cells=len(rows),passed=all(r['bracket_passed'] and r['finite_passed'] for r in rows),rows=rows,
        maximum_root_width=str(max(widths)),maximum_finite_momentum_comparison=str(max(scores)),
        display_only=dict(maximum_root_width=float(max(widths)),maximum_finite_momentum_comparison=float(max(scores))),
        scope='Neutral roots and first/second implicit derivatives of the declared ideal electron EOS at exact fully ionized ion inventories. Full physical EOS and native equilibrium remain open.')


def controls():
    p=bindings();symbolic();mp.mp.dps=p['control_momentum_digits']
    for rel,digest in json.loads((OUT/'candidate-manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    labels=[evaluate(i,1,f'control-{i}') for i in p['control_positions']]
    r=analyze(labels);a=dict(np.load(OUT/'candidates.npz'));checks=[]
    for i in p['control_positions']:
        ref=independent_momentum(a['eta_center'][i],a['beta'][i]);score=float(np.max(abs(a['momentum'][:,i]-ref)/np.maximum(1,abs(ref))))
        checks.append(dict(classification='Counterexample candidate',cell=int(a['cells'][i]),reference=ref.tolist(),score=score,passed=score<float(p['finite_comparison_gate'])))
    r['independent_momentum_controls']=checks;r['passed']=r['passed'] and all(x['passed'] for x in checks)
    save('controls.json',r);assert r['passed'];print('PASS neutral root controls',r['display_only'],flush=True)


def run():
    p=bindings();assert json.loads((OUT/'controls.json').read_text())['passed'];labels=[]
    with ProcessPoolExecutor(max_workers=p['processes']) as pool:
        work=[pool.submit(evaluate,i,min(p['block_size'],p['cells']-i),f'block-{i:04d}') for i in range(0,p['cells'],p['block_size'])]
        for done in as_completed(work):
            labels.append(done.result());save('progress.json',dict(completed_blocks=sorted(labels)));print('DENSITY BLOCKS',len(labels),flush=True)
    r=analyze(sorted(labels));save('result.json',r);assert r['passed'] and r['cells']==p['cells']
    save('manifest.json',dict(sha256={x.relative_to(g.ROOT).as_posix():g.c.sha(x) for x in OUT.iterdir() if x.is_file()}))
    verify()


def verify():
    p=bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    runtime=json.loads((OUT/'runtime.json').read_text());assert g.c.sha(Path(runtime['binary']))==runtime['binary_sha256']
    for path,digest in runtime['libraries'].items():assert g.c.sha(Path(path))==digest,path
    r=json.loads((OUT/'result.json').read_text());assert r['passed'] and r['cells']==p['cells']==3206
    assert json.loads((OUT/'symbolic.json').read_text())['passed'] and json.loads((OUT/'controls.json').read_text())['passed']
    print('PASS3206 neutral ideal-electron roots and implicit density/temperature Hessians; full physical EOS remains open',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
