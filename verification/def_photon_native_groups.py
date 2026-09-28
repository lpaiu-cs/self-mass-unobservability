"""Native frequency-group absorption on the current GR surface state.

Groups retain provider line integrals; sampled spectra remain separately
identified and are not rescaled to reproduce those integrals.
"""
from pathlib import Path
from types import FunctionType
import argparse
import json
import shutil
import time
import numpy as np
from bs4 import BeautifulSoup
from scipy.integrate import quad
import lanl_tops_stellar_probe as retrieval
import lanl_tops_stellar_audit as audit
import lanl_tops_cutoff_reader as reader
import def_ionic_structure_transport as model

ex=model.ex;h=model.h;OUT=model.OUT.parent/'def-photon-native-groups'


def prepare():
    assert not OUT.exists();OUT.mkdir()
    form=h.ROOT/'outputs/lanl-tops-form.html';shutil.copy2(form,OUT/form.name)
    defaults=retrieval.fields(BeautifulSoup(form.read_bytes(),'html.parser').find('form'))
    data,physical=model.base.inputs();c=model.base.thermal.g.c
    Y=data['X'][0]/c.A;mass=float(Y@c.W);rho=float(np.exp(data['lnd'][0])*mass)
    parts=[]
    for z in sorted(set(c.Z.astype(int))):
        mask=c.Z==z;amount=float(Y[mask].sum())
        if amount:parts.append(f'{amount:.12g} {z} {float(Y[mask]@c.W[mask]/amount):.10g}')
    rows=[]
    for count in [257,1000]:
        for T in ['.0015','.002']:
            name=f'surf{count}T'+T.replace('.','')
            assert name.isalnum() and len(name)<=15
            fields=dict(defaults,lib='new',fractype='atomic',isotope='isotope',mixture=' '.join(parts),mixname=name,
                tgrid='specific',temps=T,rgrid='specific',dens=f'{rho:.12g}',
                datype='contgrup',plasnu='off',egrid='range',egplow='1e-8',egphigh='1',ngpengs=str(count),espace='log')
            rows.append(dict(name=name,cell=0,boundaries=count,plasma_cutoff='off',fields=fields))
    ex.write(OUT/'requests.json',rows)
    paths=[Path(__file__),Path(retrieval.__file__),Path(audit.__file__),reader.OUT/'corrected-header.py',
        model.base.thermal.OUT/'coefficients.npz',OUT/form.name,OUT/'requests.json']
    ex.write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='0f6b09ed',
        claim='Supply source-resolved photon absorption groups and scattering spectra at the exact present surface mixture/density, resolving the previous sampled-spectrum Planck integral mismatch without mean calibration.',
        source='https://aphysics2.lanl.gov/static/opacdocs/opac-help.html',
        actual_state=dict(cell=0,T_K=float(physical['T'][0]),T_keV=float(model.base.k*physical['T'][0]/1.602176634e-16),rho_neutral_g_cm3=rho,neutral_mass_per_baryon=mass),
        requests='Four single-temperature queries, two bracketing native temperatures and 257/1000 logarithmic boundaries. Every submit/results pair has its own cookie session. No paired-temperature request and no hyphenated display name.',
        acceptance='Exact printed input intervals and material identities; no warnings or flat density substitutes. Positive finite groups, reject the exact 1e10 sentinel. Recombine native Planck/Rosseland groups at relative 1e-3 and compare group counts. Preserve the independent sampled-spectrum discrepancy.',
        gates=dict(group_mean_relative=.001,input_printing_score=1.,spectrum_sum_printing_interval=True),
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths},
        budget=dict(network_queries=4,requests_per_query=2,per_request_timeout_seconds=90,query_process_hard_seconds=190,CPU_workers=1,new_native_EOS_calls=0,new_stellar_steps=0,automatic_retry=False),
        limits='Provider LTE populations, bracketing temperatures and scalar opacity data are not a complete angular/energy redistribution kernel or whole-star thermal closure. No opacity/EOS model-error certificate.'))


def fetch(index):
    plan=json.loads((OUT/'plan.json').read_text())
    for p,sha in plan['bindings'].items():assert h.digest(h.ROOT/p)==sha,p
    row=json.loads((OUT/'requests.json').read_text())[index];assert not list(OUT.glob(row['name']+'-*'))
    env=dict(vars(retrieval),OUT=OUT);env['save']=lambda name,value:ex.write(OUT/name,value)
    function=FunctionType(retrieval.fetch.__code__,env,argdefs=retrieval.fetch.__defaults__)
    start=time.monotonic()
    try:result=function(row)
    except Exception as exc:
        ex.write(OUT/f"{row['name']}-failure.json",dict(classification='Counterexample candidate',error=repr(exc),seconds=time.monotonic()-start));raise
    result['seconds']=time.monotonic()-start;ex.write(OUT/f"{row['name']}-retrieval.json",result);print('RETRIEVAL',result,flush=True)


def analyze():
    assert not (OUT/'result.json').exists();results=[]
    header=reader.namespace()['header']
    for request in json.loads((OUT/'requests.json').read_text()):
        name=request['name'];path=OUT/(name+'-table.txt')
        common=header(path,dict(request,fields=dict(request['fields'],datype='gray')))
        data=audit.parse(path);assert len(data['temperatures'])==1;token=data['temperatures'][0];T=float(token)
        lines=[s.strip() for s in path.read_text().splitlines() if s.strip()]
        j=next(j for j,s in enumerate(lines) if s.startswith('Photon grid'));count=int(lines[j].split()[-2]);assert count==request['boundaries']-1
        printed=[];j+=1
        while not lines[j].startswith('Rosseland'):printed.extend(lines[j].split());j+=1
        bounds=np.geomspace(float(request['fields']['egplow']),float(request['fields']['egphigh']),count+1)
        echo=[audit.score(audit.F(float(x)),y) for x,y in zip(bounds[:-1],printed)]
        second=json.loads((OUT/(name+'-results-request.json')).read_text())
        assert float(second['egplow'])==bounds[0] and float(second['egphigh'])==bounds[-1]
        j=next(j for j,s in enumerate(lines) if s.startswith('Energy') and 'density =' in s)
        group=np.array([s.split() for s in lines[j+1:j+1+count]],float)
        assert group.shape==(count,3) and np.array_equal(group[:,0],np.array(printed,float))
        spectrum=data['spectra'][token];failures=[]
        for n,t in enumerate(data['tokens'][token][1]):
            lt,ut=audit.interval(t[1]);la,ua=audit.interval(t[2]);ls,us=audit.interval(t[3])
            if lt>ua+us or ut<la+ls:failures.append(n)
        weights=[]
        for lo,hi in zip(bounds[:-1]/T,bounds[1:]/T):
            def p(u):return 15/np.pi**4*u**3*np.exp(-u)/(-np.expm1(-u))
            def r(u):return 15/(4*np.pi**4)*u**4*np.exp(-u)/(-np.expm1(-u))**2
            weights.append([quad(p,lo,hi,epsabs=1e-14,epsrel=1e-11)[0],quad(r,lo,hi,epsabs=1e-14,epsrel=1e-11)[0]])
        weights=np.array(weights);native=np.array(common['means'][token][1:3])[::-1]
        reconstructed=np.array([weights[:,0]@group[:,2]/weights[:,0].sum(),weights[:,1].sum()/np.sum(weights[:,1]/group[:,1])])
        error=abs(reconstructed/native-1)
        finite=np.isfinite(group).all() and (group>0).all() and (group[:,1:]!=1e10).all()
        finite=finite and np.isfinite(spectrum).all() and (np.diff(spectrum[:,0])>0).all() and (spectrum[:,1:]>0).all() and np.ptp(spectrum[:,2])>0
        passed=bool(common['input_passed'] and max(echo)<=1 and finite and not failures and error.max()<.001)
        np.savez_compressed(OUT/f'{name}.npz',bounds_keV=bounds,groups=group,weights=weights,spectrum=spectrum,gray=np.array(common['means'][token]),T_keV=T)
        results.append(dict(name=name,passed=passed,groups=count,T_keV=T,input_max_score=max(common['input_max_score'],max(echo)),
            mean_relative=error.tolist(),native_Planck_Rosseland=native.tolist(),reconstructed=reconstructed.tolist(),spectrum_sum_failures=failures))
    refinements=[max(abs(np.array(a['reconstructed'])/np.array(b['reconstructed'])-1)) for a,b in zip(results[:2],results[2:])]
    result=dict(classification='Counterexample candidate',passed=all(r['passed'] for r in results) and max(refinements)<.001,
        checks=results,group_count_recombination_difference=refinements,photon_collision_kernel_identified=False,whole_star_heat_closed=False,full_dynamic_charge_solved=False)
    ex.write(OUT/'result.json',result);print('RESULT',result,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','fetch','analyze']);parser.add_argument('--index',type=int)
    args=parser.parse_args();fetch(args.index) if args.action=='fetch' else globals()[args.action]()
