"""Request 13: isolated live engine accuracy audit. Run inside the existing WSL.

These comparisons are diagnostics, not rigorous enclosures of the continuum ODE.
"""
import argparse
import difflib
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'outputs/research-remediation'
RUNTIME = Path.home()/'work/nutimo_pilot'
SRC = Path(os.environ.get('TIMING13_SOURCE',str(RUNTIME/'nutimo_request13/src')))
PAR = 'parfile-planetGR-max-bestfit'
TIM = '0337_20211005-sorted-sun5deg-res25microsec-58631_58780_clipped.tim'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def build():
    assert not SRC.exists(), 'Do not overwrite an existing build'
    shutil.copytree(RUNTIME/'nutimo_request12/src', SRC)
    patch = []
    for name in ['AllTheories3Bodies.cpp', 'Fittriple-compute.cpp']:
        path = SRC/name
        old = path.read_text()
        new = old
        if name.startswith('All'):
            # Check the index BEFORE dereferencing it (short circuit order matters).
            a = 'interv.second > tpos[tn] and tn < nt'
            assert new.count(a) == 1
            new = new.replace(a, 'tn < nt and interv.second > tpos[tn]')
            a = 'interv.second > tneg[tn] and tn >= 0'
            assert new.count(a) == 1
            new = new.replace(a, 'tn >= 0 and interv.second > tneg[tn]')
        else:
            a = 'value_type errtoe = pow(dix, -14);'
            assert new.count(a) == 1
            new = new.replace(a, 'const char* e13 = getenv("TIMING13_ERRTOE");\n    value_type errtoe = e13 ? strtold(e13, NULL) : pow(dix, -14);')
        path.write_text(new)
        patch += difflib.unified_diff(old.splitlines(True), new.splitlines(True),
                                    fromfile='request12/'+name, tofile='request13/'+name)
    (OUT/'runtime.patch').write_text(''.join(patch))
    command = json.loads((ROOT/'outputs/research-completion/runtime12/build.json').read_text())['command']
    with (OUT/'build.log').open('w') as log:
        result = subprocess.run(command, cwd=SRC, stdout=log, stderr=subprocess.STDOUT)
    assert result.returncode == 0, 'See build.log'
    (OUT/'build.json').write_text(json.dumps(dict(command=command, library_sha256=sha(SRC/'libFittriplecpp.so'),
        source_sha256={n:sha(SRC/n) for n in ['AllTheories3Bodies.cpp','Fittriple-compute.cpp']}),indent=2)+'\n')


def live(label, tolerance='1e-16', mesh=250, errtoe='1e-14'):
    run = RUNTIME/('run_request13_'+label)
    if not run.exists():
        shutil.copytree(RUNTIME/'run_request12', run)
        p = run/PAR
        text = p.read_text()
        text,n = re.subn(r'(?m)^integration_tolerance\s+\S+', 'integration_tolerance '+tolerance,text)
        assert n == 1
        text,n = re.subn(r'(?m)^interpsteps\s+\S+', 'interpsteps '+str(mesh),text)
        assert n == 1
        p.write_text(text)
    else:
        text = (run/PAR).read_text()
        assert float(re.search(r'(?m)^integration_tolerance\s+(\S+)',text)[1]) == float(tolerance)
        assert int(re.search(r'(?m)^interpsteps\s+(\S+)',text)[1]) == mesh
    os.chdir(run)
    for key in ['SEPDYN_A','SEPDYN_W','SEPDYN_PH','SEPDYN_TAU']:
        os.environ.pop(key,None)
    os.environ['TIMING13_ERRTOE'] = errtoe
    sys.path.insert(0,str(SRC))
    import python_Fittriple_interface as pfi
    start = time.monotonic()
    fit = pfi.PyFittriple(PAR,TIM)
    fit.Compute_lnposterior()
    r = fit.Get_time_residuals().copy()
    base = np.load(ROOT/'request10_external/baseline_planetGR.npz',allow_pickle=True)
    assert len(r) == len(base['res']) and np.isfinite(r).all()
    report = dict(label=label,tolerance=tolerance,mesh=mesh,errtoe_days=errtoe,
        source_sha256=sha(SRC/'AllTheories3Bodies.cpp'),compute_sha256=sha(SRC/'Fittriple-compute.cpp'),
        library_sha256=sha(SRC/'libFittriplecpp.so'),interface_sha256=sha(pfi.__file__),
        par_sha256=sha(PAR),tim_sha256=sha(TIM),baseline_seconds=time.monotonic()-start,
        archived_max_difference_us=float(np.max(abs(r-base['res']))))
    (OUT/(label+'-preflight.json')).write_text(json.dumps(report,indent=2)+'\n')
    np.savez(OUT/(label+'-baseline.npz'),res=r)
    print('PREFLIGHT',report,flush=True)
    return fit,base,r


def stable_build():
    source=RUNTIME/'nutimo_request13/src'
    target=RUNTIME/'nutimo_request13_stable/src'
    assert not target.exists()
    shutil.copytree(source,target)
    p=target/'Fittriple-compute.cpp'; old=p.read_text(); new=old
    a='value_type * delay ;'
    assert new.count(a)==1
    new=new.replace(a,a+'\n    vector<value_type> emission_delay(ntoas, zero);')
    a='toes[i] = Te1 - delay_einstein_i( Te1, toa_in_interp[i] - marginmin, toa_in_interp[i] + marginsup ) ;'
    assert new.count(a)==1
    new=new.replace(a,a+'''\n            emission_delay[i] = delays_interpolation(ta-delay0, toa_in_interp[i]-marginmin, toa_in_interp[i]+marginsup)
                + delay_geom_i(ta-delay0, ntisaround/2-marginmin, ntisaround/2+marginmin)
                + delay_einstein_i(Te1, toa_in_interp[i]-marginmin, toa_in_interp[i]+marginsup);''')
    a='if (fractional == 0)'
    assert new.count(a)==1
    new=new.replace(a,'''if (fractional == 0 && remove_mean
            && strcmp(parameters.specialcase,"RN_PL") != 0
            && strcmp(parameters.specialcase,"Circum-ternary_Kepler") != 0)
        {
            const value_type D=parameters.DopplerF;
            const value_type T0=toas[0]*D-parameters.treference;
            const value_type d0=emission_delay[0];
            for (i=0; i<ntoas; ++i) {
                const value_type Ti=toas[i]*D-parameters.treference;
                const value_type di=emission_delay[i];
                value_type r=fmal(spinfreq,toas[i]-toas[0],-static_cast<value_type>(turns[i]));
                r += undemi*spinfreq1*(toas[i]-toas[0])*(toas[i]+toas[0]-deux*parameters.treference/D);
                r -= (spinfreq/D+spinfreq1*Ti/(D*D))*di
                    -(spinfreq/D+spinfreq1*T0/(D*D))*d0;
                r += undemi*spinfreq1/(D*D)*(di*di-d0*d0);
                residuals[i]=r*spinPeriodmicrosec;
                mean+=residuals[i]*weights[i];
            }
        }
        else if (fractional == 0)''')
    p.write_text(new)
    (OUT/'stable-residual.patch').write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),
        fromfile='request13/Fittriple-compute.cpp',tofile='request13_stable/Fittriple-compute.cpp')))
    command=json.loads((OUT/'build.json').read_text())['command']
    with (OUT/'stable-build.log').open('w') as log:
        result=subprocess.run(command,cwd=target,stdout=log,stderr=subprocess.STDOUT)
    assert result.returncode==0
    (OUT/'stable-build.json').write_text(json.dumps(dict(command=command,source_sha256=sha(p),
        library_sha256=sha(target/'libFittriplecpp.so')),indent=2)+'\n')


def audit(args):
    fit,base,r = live(args.label,args.tolerance,args.mesh,args.errtoe)
    meta = json.loads((ROOT/'request10_external/finite_jacobian_v2_meta.json').read_text())
    scales = base['scales'][base['fmap']].astype(float)
    names = [str(base['names'][j]) for j in base['fmap']]
    assert names == meta['columns']
    for j in map(int,args.columns.split(',')):
        for fraction in [1.,.5]:
            path = OUT/f'{args.label}-jac-{j:02d}-{fraction}.npz'
            if path.exists():
                continue
            h = meta['abs_steps'][j]*fraction
            results = []
            for sign in [1,-1]:
                delta = np.zeros(28); delta[j] = sign*h/scales[j]
                fit.Set_fitted_parameter_relativeshifts(delta)
                fit.Compute_lnposterior()
                results.append(fit.Get_time_residuals().copy())
            plus,minus = results
            assert np.isfinite(plus).all() and np.isfinite(minus).all()
            assert max(abs(plus-r).max(),abs(minus-r).max()) < 1000
            np.savez(path,plus=plus,minus=minus,h=h,dcol=(plus-minus)/(2*h),name=names[j])
            print('JAC',args.label,j,fraction,flush=True)
    fit.Set_fitted_parameter_relativeshifts(np.zeros(28))
    fit.Compute_lnposterior()
    recovery = float(np.max(abs(fit.Get_time_residuals()-r)))
    assert recovery < 1e-7
    (OUT/(args.label+'-done.json')).write_text(json.dumps(dict(recovery_us=recovery))+'\n')
    print('DONE',args.label,flush=True)


def gradient(args):
    fit,base,r=live('gradient_'+args.label+'_'+str(args.worker))
    checkpoint=OUT/('nonlinear-'+args.label+'.npz')
    saved=np.load(checkpoint)
    z=saved['scaled_parameters']
    meta=json.loads((ROOT/'request10_external/finite_jacobian_v2_meta.json').read_text())
    h=np.asarray(meta['abs_steps'])
    scales=base['scales'][base['fmap']].astype(float)
    for j in range(args.worker,28,args.workers):
        pair=[]
        for sign in [1.,-1.]:
            trial=z.copy(); trial[j]+=sign*.5
            fit.Set_fitted_parameter_relativeshifts(h*trial/scales)
            fit.Compute_lnposterior(0)
            pair.append(fit.Get_time_residuals().copy())
        plus,minus=pair
        assert np.isfinite(plus).all() and np.isfinite(minus).all()
        np.savez(OUT/f'gradient-{args.label}-{j:02d}.npz',plus=plus,minus=minus,dcol=plus-minus,
            checkpoint_sha256=sha(checkpoint),scaled_parameters=z)
        print('GRADIENT',args.label,j,flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('mode',choices=['build','audit','gradient','stable_build'])
    p.add_argument('--label',default='control')
    p.add_argument('--tolerance',default='1e-16')
    p.add_argument('--mesh',type=int,default=250)
    p.add_argument('--errtoe',default='1e-14')
    p.add_argument('--columns',default='21,22,27')
    p.add_argument('--worker',type=int,default=0)
    p.add_argument('--workers',type=int,default=1)
    a = p.parse_args()
    OUT.mkdir(exist_ok=True,parents=True)
    {'build':lambda:build(),'audit':lambda:audit(a),'gradient':lambda:gradient(a),'stable_build':lambda:stable_build()}[a.mode]()
