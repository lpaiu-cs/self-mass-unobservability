"""Audit gray/group outputs and test an explicitly exploratory cutoff model."""
from decimal import Decimal
from fractions import Fraction as F
import json, re, shutil, sys
import numpy as np
from scipy.constants import physical_constants, hbar, epsilon_0, elementary_charge, m_e, Avogadro
from scipy.integrate import simpson, trapezoid
import lanl_tops_boundary_control as retrieval

g=retrieval.g;audit=retrieval.audit;OUT=g.OUT/'lanl-tops-cutoff-audit'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();retrieval.verify();audit.verify()
    for name in ['lanl-opacity-tables2025.pdf','lanl-opacity-tables2025.txt']:
        shutil.copy2(g.ROOT/'outputs'/name,OUT/name)
    paths=[g.ROOT/'verification/lanl_tops_cutoff_audit.py',retrieval.OUT/'manifest.json',
           audit.OUT/'manifest.json']+list(OUT.iterdir())
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='361c8e0',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        source_url='https://www.osti.gov/servlets/purl/2570872',
        input_gate='Exact decimal printing intervals for all composition fractions, temperatures and densities; same material IDs as prior accepted mixtures. Reject warnings. The table reports group lower edges; validate all submitted boundaries including the final upper edge against the actual returned request form.',
        group_gate='Compare cutoff-OFF group means with positive log-interpolation of the already retrieved spectra, both Simpson and trapezoid quadrature, each normalized over its own requested energy interval. Relative tolerance 1e-3, fixed before reconstruction. Do not widen any failed gate.',
        group_relative_tolerance=1e-3,
        cutoff_diagnostic='Exploratory, after observing ON/OFF gray values: compute the nonrelativistic electron plasma energy using returned free electrons per atom, the submitted neutral mass density and mean atomic mass. Compare normalized integrals above this independently computed threshold to both central ON tables. Do not fit a threshold to the gray data. Report errors with no physical acceptance claim; the service source and dense-plasma dispersion remain distinct from this diagnostic.',
        outer='A successful gray output is a separate output-path result. The earlier frequency request remains failed, and low-density LTE physical validity is not asserted.',
        physical_acceptance=False))


def header(path,request):
    lines=[x.strip() for x in path.read_text().splitlines() if x.strip()]
    fields=request['fields'];tokens=fields['mixture'].split()
    amounts={int(tokens[i+1]):(F(Decimal(tokens[i])),F(Decimal(tokens[i+2]))) for i in range(0,len(tokens),3)}
    nt=sum(x[0] for x in amounts.values());mt=sum(x[0]*x[1] for x in amounts.values())
    start=next(i for i,x in enumerate(lines) if x.startswith('No. Fraction'))+1;echo=[];ids={}
    while not lines[start].startswith('Temperature grid'):
        number,mass,Z,symbol,mat=lines[start].split();Z=int(Z);n,w=amounts[Z]
        echo.extend([audit.score(n/nt,number),audit.score(n*w/mt,mass)]);ids[Z]=int(mat);start+=1
    assert set(ids)==set(amounts)
    expected=next(x for x in json.loads((audit.OUT/'result.json').read_text())['rows'] if x.get('accepted'))['material_ids']
    assert ids=={x['Z']:x['returned'] for x in expected}
    T=lines[start+1].split();rho=lines[start+3].split();assert len(rho)==1
    Ts=sorted(F(Decimal(x)) for x in fields['temps'].split());assert len(Ts)==len(T)
    echo.extend(audit.score(a,b) for a,b in zip(Ts,sorted(T,key=Decimal)))
    density=F(Decimal(fields['dens']));echo.append(audit.score(density,rho[0]));means={}
    for i,line in enumerate(lines):
        if line.startswith('Density') and 'T=' in line:
            t=line.split('T=')[1].strip();v=lines[i+1].split();assert len(v)==5
            echo.append(audit.score(density,v[0]));means[t]=list(map(float,v))
    assert set(means)==set(T)
    groups=[];energies=[]
    if fields['datype']=='groups':
        i=next(i for i,x in enumerate(lines) if x.startswith('Photon grid'))
        count=int(re.search(r'(\d+) points',lines[i]).group(1));i+=1
        while not lines[i].startswith('Rosseland'):energies.extend(lines[i].split());i+=1
        assert len(energies)==count
        submitted=fields['energies'].split();assert len(submitted)==count+1
        second=json.loads((retrieval.OUT/(request['name']+'-results-request.json')).read_text())
        assert second['energies'].split()==submitted
        echo.extend(audit.score(F(Decimal(a)),b) for a,b in zip(submitted[:-1],energies))
        i=next(i for i,x in enumerate(lines) if x.startswith('Energy') and 'density =' in x)
        tg,rg=lines[i].split('=')[1].split();echo.append(audit.score(density,rg))
        echo.append(audit.score(Ts[0],tg));groups=np.array([x.split() for x in lines[i+1:]],float)
        assert groups.shape==(count,3) and np.all(groups>0) and np.all(np.isfinite(groups))
        assert np.array_equal(groups[:,0],np.array(energies,float))
    passed=max(echo)<=1 and not any('warning' in x.lower() for x in lines)
    return dict(input_passed=passed,input_max_score=max(echo),means=means,groups=np.asarray(groups),
                mass_sum=float(mt),number_sum=float(nt))


def integrate(a,T,lo,hi):
    assert a[0,0]<=lo<hi<=a[-1,0]
    interior=a[(a[:,0]>lo)&(a[:,0]<hi)]
    def endpoint(E):return np.r_[E,np.exp([np.interp(np.log(E),np.log(a[:,0]),np.log(a[:,j])) for j in [1,2,3]])]
    b=np.vstack([endpoint(lo),interior,endpoint(hi)]);x=b[:,0]/T
    wP=np.zeros_like(x);wR=np.zeros_like(x);mask=x<700;y=x[mask];em=np.exp(-y);den=-np.expm1(-y)
    wP[mask]=y**3*em/den;wR[mask]=y**4*em/(den*den)
    values=[];weights=[]
    for method in [simpson,trapezoid]:
        p=method(wP,x=x);r=method(wR,x=x);weights.append([p*15/np.pi**4,r*15/(4*np.pi**4)])
        values.append([method(wP*b[:,2],x=x)/p,r/method(wR/b[:,1],x=x)])
    return np.array(values),np.array(weights)


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    requests=json.loads((retrieval.OUT/'requests.json').read_text());parsed={};rows=[]
    for request in requests:
        name=request['name'];data=header(retrieval.OUT/(name+'-table.txt'),request);parsed[name]=data
        rows.append(dict(name=name,input_passed=data['input_passed'],input_max_score=data['input_max_score'],means=data['means']))
    off=audit.parse(audit.retrieval.OUT/'cell-5734-off-table.txt');a=off['spectra']['1.7500E+00'];T=1.75
    request=next(x for x in requests if x['name']=='coregroupsoff');bounds=np.array(request['fields']['energies'].split(),float)
    group_controls=[]
    for row,lo,hi in zip(parsed['coregroupsoff']['groups'],bounds[:-1],bounds[1:]):
        values,weights=integrate(a,T,lo,hi);native=row[1:][::-1];error=abs(values/native-1)
        mutual=abs(values[0]/values[1]-1)
        group_controls.append(dict(bounds_keV=[lo,hi],native_Planck_Rosseland=native.tolist(),
            integrated_Planck_Rosseland=values.tolist(),relative_error=error.tolist(),mutual_difference=mutual.tolist(),
            passed=bool(error.max()<plan['group_relative_tolerance'] and mutual.max()<plan['group_relative_tolerance'])))
    on=parsed['coregroupson']['groups'];goff=parsed['coregroupsoff']['groups']
    sentinel=np.all(on[:,1:]==1e10,axis=1);assert np.any(sentinel) and np.any(~sentinel)
    cutoff_rows=[];central=parsed['coregroupsoff'];W=central['mass_sum']/central['number_sum']
    density=float(request['fields']['dens'])
    constants=dict(hbar_J_s=hbar,epsilon0_SI=epsilon_0,e_C=elementary_charge,m_e_kg=m_e,Avogadro=Avogadro,
                   atomic_mass_constant_kg=physical_constants['atomic mass constant'][0])
    save('physical-constants.json',dict(classification='Imported from prior work',source='installed scipy.constants CODATA',values=constants))
    for t,a in off['spectra'].items():
        rho,R,P,zbar,z2=off['means'][t];zbar=float(zbar)
        # Number density uses a neutral atom of mass W times the SI atomic-mass constant.
        ne_SI=density*1000/W/physical_constants['atomic mass constant'][0]*zbar
        Ep=hbar*np.sqrt(ne_SI*elementary_charge**2/(epsilon_0*m_e))/(1000*elementary_charge)
        values,weights=integrate(a,float(t),Ep,a[-1,0])
        target=np.array(audit.parse(audit.retrieval.OUT/'cell-5734-on-table.txt')['means'][t][1:3],float)[::-1]
        cutoff_rows.append(dict(T_keV=float(t),nonrelativistic_plasma_energy_keV=float(Ep),electron_density_m3=float(ne_SI),
            native_ON_Planck_Rosseland=target.tolist(),conditional_tail_Planck_Rosseland=values.tolist(),
            relative_error_to_ON=abs(values/target-1).tolist(),vacuum_tail_Planck_Rosseland_weight=weights.tolist(),
            fitted_to_gray_means=False,physical_or_provider_formula_certified=False))
    save('result.json',dict(classification='Counterexample candidate',completed=True,rows=rows,
        all_input_gates_passed=all(x['input_passed'] for x in rows),group_controls=group_controls,
        all_group_controls_passed=all(x['passed'] for x in group_controls),
        cutoff_sentinel_lower_edges_keV=on[sentinel,0].tolist(),
        groups_above_cutoff_bitwise_equal=bool(np.array_equal(on[~sentinel],goff[~sentinel])),
        cutoff_exploratory_controls=cutoff_rows,
        outer_ON_OFF_gray_equal=parsed['outergrayon']['means']==parsed['outergrayoff']['means'],
        outer_frequency_failure_erased=False,physical_LTE_atmosphere_certified=False,full_GR_evolution=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify()


def verify():
    for name,key in [('plan.json','bindings'),('manifest.json','sha256')]:
        for rel,digest in json.loads((OUT/name).read_text())[key].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'result.json').read_text())['completed']
    print('PASS TOPS cutoff diagnostic bindings; inspect finite gates, no physical cutoff certificate',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
