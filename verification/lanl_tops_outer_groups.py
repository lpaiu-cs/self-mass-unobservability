"""Retrieve native outer group integrals after monochromatic quadrature failed."""
from types import FunctionType
from decimal import Decimal
from fractions import Fraction as F
import json, shutil, sys
import numpy as np
from scipy.integrate import quad
import lanl_tops_outer_audit as spectral
import lanl_tops_cutoff_reader as reader

g=spectral.g;original=spectral.retrieval.original;OUT=g.OUT/'lanl-tops-outer-groups'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();spectral.verify();reader.verify()
    shutil.copy2(original.OUT/'lanl-tops-form.html',OUT/'lanl-tops-form.html')
    base=next(x for x in json.loads((spectral.retrieval.control.audit.retrieval.OUT/'requests.json').read_text()) if x['cell']==0)
    rows=[]
    for count in [16,32]:
        name=f'outer{count}groups';energies=' '.join(format(x,'.6g') for x in np.geomspace(1e-8,1,count+1))
        assert len(energies)<450
        fields=dict(base['fields'],mixname=name,datype='groups',plasnu='off',egrid='specific',energies=energies)
        rows.append(dict(name=name,cell=0,fields=fields))
    save('requests.json',rows)
    paths=[g.ROOT/'verification/lanl_tops_outer_groups.py',spectral.OUT/'manifest.json',
           reader.OUT/'manifest.json',OUT/'requests.json',OUT/'lanl-tops-form.html']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='f305247',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        reason='Both separately retrieved outer spectra pass all input and arithmetic checks but fail the frozen Planck quadrature gate. Preserve both failures. Retrieve native Planck/Rosseland group averages at the same exact outer mixture/density and both temperatures with cutoff OFF.',
        grid='16 and 32 logarithmic groups with explicitly submitted six-digit boundaries from 1e-8 to 1 keV. These are provider group integrals, not reconstructions from the failed output spectra.',
        acceptance='Check exact input printing intervals, every boundary and final upper edge, absence of warnings and strictly positive finite group means. Weight native Planck groups by their exact submitted Planck intervals and Rosseland groups by the corresponding derivative weights. Recombine all groups and compare each gray mean and the 16/32 recombinations at fixed relative 1e-3. Do not replace prior failed spectral gates.',
        relative_tolerance=1e-3,physical_LTE_atmosphere_certified=False))


def namespace():
    env=dict(vars(original),OUT=OUT)
    for name in ['save','bindings','fetch','run','verify']:
        fn=getattr(original,name);env[name]=FunctionType(fn.__code__,env,argdefs=fn.__defaults__)
    return env


def run():namespace()['run']()


def analyze():
    namespace()['verify']();plan=json.loads((OUT/'plan.json').read_text());rows=[];allvalues={}
    header=reader.namespace()['header'];audit=spectral.audit
    for request in json.loads((OUT/'requests.json').read_text()):
        name=request['name'];path=OUT/(name+'-table.txt')
        # Reuse only the common composition/temperature/density part of the reader.
        common=header(path,dict(request,fields=dict(request['fields'],datype='gray')))
        lines=[x.strip() for x in path.read_text().splitlines() if x.strip()];fields=request['fields']
        bounds=np.array(fields['energies'].split(),float);count=len(bounds)-1
        second=json.loads((OUT/(name+'-results-request.json')).read_text())
        assert F(Decimal(second['egplow']))==F(Decimal(fields['energies'].split()[0]))
        assert F(Decimal(second['egphigh']))==F(Decimal(fields['energies'].split()[-1]))
        i=next(i for i,x in enumerate(lines) if x.startswith('Photon grid'));assert int(lines[i].split()[-2])==count
        i+=1;printed=[]
        while not lines[i].startswith('Rosseland'):printed.extend(lines[i].split());i+=1
        assert len(printed)==count
        echo=[audit.score(F(Decimal(x)),y) for x,y in zip(fields['energies'].split()[:-1],printed)]
        groups={}
        for i,line in enumerate(lines):
            if line.startswith('Energy') and 'density =' in line:
                T,rho=line.split('=')[1].split();echo.append(audit.score(F(Decimal(fields['dens'])),rho))
                group=np.array([x.split() for x in lines[i+1:i+1+count]],float)
                assert group.shape==(count,3) and np.array_equal(group[:,0],np.array(printed,float))
                groups[T]=group
        assert set(groups)==set(common['means'])
        for T,group in groups.items():
            temperature=float(T);weights=[]
            for lo,hi in zip(bounds[:-1]/temperature,bounds[1:]/temperature):
                def p(x):return 15/np.pi**4*x**3*np.exp(-x)/(-np.expm1(-x))
                def r(x):return 15/(4*np.pi**4)*x**4*np.exp(-x)/(-np.expm1(-x))**2
                weights.append([quad(p,lo,hi,epsabs=1e-13,epsrel=1e-11)[0],quad(r,lo,hi,epsabs=1e-13,epsrel=1e-11)[0]])
            weights=np.array(weights);values=np.array([weights[:,0]@group[:,2]/weights[:,0].sum(),weights[:,1].sum()/np.sum(weights[:,1]/group[:,1])])
            native=np.array(common['means'][T][1:3])[::-1];error=abs(values/native-1)
            finite=bool(np.all(np.isfinite(group)) and np.all(group>0) and np.all(group[:,1:]<1e10))
            passed=bool(common['input_passed'] and max(echo)<=1 and finite and error.max()<plan['relative_tolerance'])
            rows.append(dict(name=name,groups=count,T_keV=temperature,input_max_score=max(max(echo),common['input_max_score']),
                native_gray_Planck_Rosseland=native.tolist(),recombined_Planck_Rosseland=values.tolist(),
                relative_errors=error.tolist(),total_weight=weights.sum(axis=0).tolist(),passed=passed))
            allvalues[count,T]=values;np.savez_compressed(OUT/f'groups-{count}-{T}.npz',bounds_keV=bounds,groups=group,weights=weights)
    comparisons=[dict(T_keV=float(T),relative_difference=abs(allvalues[16,T]/allvalues[32,T]-1).tolist()) for count,T in allvalues if count==16]
    passed=all(x['passed'] for x in rows) and max(y for x in comparisons for y in x['relative_difference'])<plan['relative_tolerance']
    save('analysis.json',dict(classification='Counterexample candidate',completed=True,rows=rows,refinement=comparisons,all_passed=bool(passed),
        original_spectral_quadrature_failure_erased=False,physical_LTE_atmosphere_certified=False,
        monochromatic_shape_certified=False,full_GR_evolution=False))
    save('analysis-manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify()


def verify():
    namespace()['verify']()
    for rel,digest in json.loads((OUT/'analysis-manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'analysis.json').read_text())['completed']
    print('PASS native outer group audit bindings; original spectral failures preserved',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
