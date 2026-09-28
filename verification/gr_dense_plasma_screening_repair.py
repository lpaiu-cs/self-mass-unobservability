"""Correct the isolated quotient derivative without changing the fitted free energy."""
import ctypes,json,subprocess,sys
import mpmath as mp
import numpy as np
import sympy as sp
import gr_dense_plasma as d
import gr_dense_plasma_derivative_probe as probe

g=d.g;OUT=g.OUT/'gr-dense-plasma-screening-repair';CACHE=g.CACHE/'dense-plasma-screening-repair'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists() and not CACHE.exists();OUT.mkdir();CACHE.mkdir();d.verify()
    original=(d.OUT/'potekhin-chabrier-eos22.f').read_text()
    old='COR0DXX=(U0DXX-(2.*U0DX*D0DX+U0*D0DXX)/D0+2.*(D0DX/D0)**2)/D0'
    new='COR0DXX=(U0DXX-(2.*U0DX*D0DX+U0*D0DXX)/D0+\n     +  2.*U0*(D0DX/D0)**2)/D0'
    assert original.count(old)==1;changed=original.replace(old,new);assert changed.replace(new,old)==original
    (OUT/'eos22-screening.f').write_text(changed)
    paths=[g.ROOT/'verification/gr_dense_plasma_screening_repair.py',probe.OUT/'manifest.json',d.OUT/'manifest.json',OUT/'eos22-screening.f']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='e9500e5',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        substitution=dict(old=old,new=new),
        problem='FSCRliq8 second x derivative of 1+U0/D0 omits U0 in the final quotient-rule term. This feeds FDXX and PDRSCR, while free energy, U, P, CV and PDT are unchanged.',
        controls=[2972,3043,4352,5734],charges=[1,2,6,20],
        mp_controls='60-decimal free energy and its independent log n/log T derivatives; preserve binary32 rounding of unsuffixed Fortran decimal literals before their promotion to double precision.',
        independent_score_tolerance=1e-10,all_unaffected_outputs_bitwise=True,
        finite_log_steps=[2e-4,1e-4],finite_derivative_tolerance=1e-5,
        policy='Retain the original failed 1e-5 mixture derivative gate. Repeat that same gate after one quotient-rule correction in a separate library. No refitting, physical input change, native EOS change or GR restart.'))


def bindings():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    return plan


def build():
    bindings();command=['gfortran','-O2','-fPIC','-shared','-std=legacy',str(OUT/'eos22-screening.f'),'-o',str(CACHE/'pc.so')]
    result=subprocess.run(command,capture_output=True,text=True)
    save('build.json',dict(command=command,returncode=result.returncode,stdout=result.stdout,stderr=result.stderr));assert result.returncode==0
    save('runtime.json',dict(sha256={str(CACHE/'pc.so'):g.c.sha(CACHE/'pc.so')},
        inherited_runtime_sha256=json.loads((d.OUT/'runtime.json').read_text())['sha256']))


class Provider(d.Provider):
    def __init__(self):self.lib=ctypes.CDLL(str(CACHE/'pc.so'))


def free_mp(rs,ge,z):
    # A decimal real literal without d0 is rounded as Fortran default real.
    q=lambda s:mp.mpf(float(np.float32(s)))
    X=q('.0140047')/rs;logz=mp.log(z);z13=mp.exp(logz/3)
    cdh=z/q('1.73205')*(mp.sqrt(z+1)**3-mp.sqrt(z)**3-1)
    ctf=z*z*q('.2513')*(z13-1+q('.2')/mp.sqrt(z13))
    p01=q('1.11')*mp.exp(q('.475')*logz);p03=q('.2')+q('.078')*logz**2;power=q('1.16')+q('.08')*logz
    tx=ge**power
    cor1=1+(z-1)/9*rs**3/(1+6*rs**2)*(1+1/(q('.001')*z*z+2*ge))
    cor0=1+q('.78')*mp.sqrt(ge/z)*rs**3/(ge*z**3+21*rs**3)
    h1=(1+X*X/5)/(1+q('.18')/mp.sqrt(mp.sqrt(z))*X+(q('.2')+q('.37')/mp.sqrt(z))*X*X)
    up=cdh*mp.sqrt(ge)+p01*ctf*tx*cor0*h1
    den=1+(p03*mp.sqrt(ge)+p01/rs*tx*cor1)/mp.sqrt(1+X*X)
    return -up/den*ge


def independent(rs,ge,z):
    mp.mp.dps=60;rs,ge,z=[mp.mpf(float(v)) for v in [rs,ge,z]]
    f=lambda s,t:free_mp(rs*mp.exp(-s/3),ge*mp.exp(s/3-t),z)
    F=f(0,0);U=-mp.diff(f,(0,0),(0,1));P=mp.diff(f,(0,0),(1,0))
    return np.array(list(map(float,[F,U,P,U-F,U-mp.diff(f,(0,0),(0,2)),
        P+mp.diff(f,(0,0),(1,1)),P+mp.diff(f,(0,0),(2,0))])))


def controls():
    plan=bindings();p=Provider();old=d.Provider();a=np.load(d.OUT/'stellar-comparison.npz');rows=[]
    x=sp.symbols('x');u=sp.Function('u')(x);v=sp.Function('v')(x)
    correct=(sp.diff(u,x,2)-(2*sp.diff(u,x)*sp.diff(v,x)+u*sp.diff(v,x,2))/v+2*u*(sp.diff(v,x)/v)**2)/v
    assert sp.simplify(sp.diff(u/v,x,2)-correct)==0
    source=(sp.diff(u,x,2)-(2*sp.diff(u,x)*sp.diff(v,x)+u*sp.diff(v,x,2))/v+2*(sp.diff(v,x)/v)**2)/v
    assert sp.simplify(source-correct-2*(1-u)*sp.diff(v,x)**2/v**3)==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        correction='(U/D)_xx=U_xx/D-(2U_x D_x+U D_xx)/D^2+2U D_x^2/D^3.',
        original_error='The distributed COR0DXX differs by 2*(1-U0)*D0_x^2/D0^3. Only the density second derivative chain is changed.'))
    for i in plan['controls']:
        k=int(np.flatnonzero(a['cells']==i)[0]);rs,ge=a['parameters'][k,:2]
        for z in plan['charges']:
            before=d.seven(old.call('fscrliq8',[rs,ge,float(z)],6));after=d.seven(p.call('fscrliq8',[rs,ge,float(z)],6))
            ref=independent(rs,ge,z);same=np.array_equal(before[:6],after[:6]);score=float(np.max(abs(after-ref)/np.maximum(1,abs(ref))))
            rows.append(dict(cell=i,Z=z,before=before.tolist(),after=after.tolist(),independent=ref.tolist(),
                unaffected_bitwise=same,corrected_score=score,passed=bool(same and score<plan['independent_score_tolerance'])))
    save('controls.json',dict(classification='Counterexample candidate',rows=rows,passed=all(r['passed'] for r in rows),
        maximum_score=max(r['corrected_score'] for r in rows)))
    assert all(r['passed'] for r in rows)
    print('SCREENING MP controls',len(rows),max(r['corrected_score'] for r in rows),flush=True)


def run():
    plan=bindings();assert json.loads((OUT/'controls.json').read_text())['passed'];p=Provider()
    state,_=d.state_data();saved=np.load(d.OUT/'stellar-comparison.npz');values=[];rows=[]
    for k,i in enumerate(saved['cells']):
        r,t,X=state['lnd'][i],state['lnT'][i],state['X'][i];b=d.mixture(p,r,t,X)['corrected'];values.append(b)
        assert np.array_equal(b[:6],saved['corrected'][k,:6]),int(i)
        if i in plan['controls']:
            for h in plan['finite_log_steps']:
                rm=d.mixture(p,r-h,t,X)['corrected'];rp=d.mixture(p,r+h,t,X)['corrected']
                tm=d.mixture(p,r,t-h,X)['corrected'];tp=d.mixture(p,r,t+h,X)['corrected']
                finite=np.array([-(tp[0]-tm[0])/(2*h),(rp[0]-rm[0])/(2*h),b[1]+(tp[1]-tm[1])/(2*h),
                    b[2]+(tp[2]-tm[2])/(2*h),b[2]+(rp[2]-rm[2])/(2*h)])
                expected=b[[1,2,4,5,6]];score=float(np.max(abs(finite-expected)/np.maximum(1,abs(expected))))
                rows.append(dict(cell=int(i),step=h,finite=finite.tolist(),expected=expected.tolist(),score=score,passed=score<plan['finite_derivative_tolerance']))
        if len(values)%512==0:print('SCREENING REPLAY',len(values),'/',len(saved['cells']),flush=True)
    values=np.array(values);np.savez_compressed(OUT/'corrected-stellar.npz',cells=saved['cells'],values=values,
        PDR_change=values[:,6]-saved['corrected'][:,6])
    save('result.json',dict(classification='Counterexample candidate',completed=True,cells=len(values),
        all_unaffected_outputs_bitwise=True,finite_controls=rows,finite_derivatives_passed=all(r['passed'] for r in rows),
        maximum_finite_score=max(r['score'] for r in rows),max_PDR_change=float(np.max(abs(values[:,6]-saved['corrected'][:,6]))),
        original_derivative_failure_preserved=True,physical_EOS_certified=False,native_EOS_replaced=False,full_GR_evolution=False))
    assert all(r['passed'] for r in rows)
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    bindings();d.verify()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    runtime=json.loads((OUT/'runtime.json').read_text())
    for path,digest in {**runtime['sha256'],**runtime['inherited_runtime_sha256']}.items():assert g.c.sha(path)==digest,path
    assert json.loads((OUT/'result.json').read_text())['finite_derivatives_passed']
    print('PASS source quotient correction and unchanged-gate stellar derivative replay',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
