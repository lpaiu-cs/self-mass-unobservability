"""Audit returned stellar TOPS spectra against the actual submitted inputs."""
from decimal import Decimal
from fractions import Fraction as F
import json, re, shutil, sys
import numpy as np
from scipy.integrate import simpson
import lanl_tops_name_runner as retrieval

g=retrieval.g;OUT=g.OUT/'lanl-tops-stellar-audit'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();retrieval.verify()
    for name in ['lanl-opac-faq.html','lanl-opac-physics.html']:
        shutil.copy2(g.ROOT/'outputs'/name,OUT/name)
    paths=[g.ROOT/'verification/lanl_tops_stellar_audit.py',retrieval.OUT/'manifest.json',
           retrieval.original.OUT/'manifest.json',OUT/'lanl-opac-faq.html',OUT/'lanl-opac-physics.html']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='faefcb8',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        input_echo='Check exact decimal printing intervals for every submitted number fraction, mass fraction, density and table temperature. No composition refit or tolerance based on a small missing species. Report all returned material identities.',
        material_ids='The returned integer material ID equals the reversal of the last four digits of the documented n1xxxxx filename convention in the observed Al-new/Al-old/He-new controls. Check this inferred interface mapping for every element; it is not a native source or physical opacity certificate.',
        spectra='Require positive strictly increasing photon energies, finite nonnegative total/absorption/scattering, nonconstant absorption, both requested temperatures and density. Check total=absorption+scattering with exact decimal output intervals. Reject all actual WARNING messages and any flat boundary spectrum.',
        finite_mean_control='For the explicitly cutoff-OFF central control only, integrate the returned absorption/total spectra against vacuum Planck/Rosseland weights with Simpson and trapezoid quadrature. Compare both to the printed gray means and to each other at relative tolerance 1e-3. These are finite spectral consistency tests; unsampled spectra and tails have no physical error bound.',
        mean_relative_tolerance=1e-3,
        cutoff_control='Compare the central cutoff-ON/OFF tables at identical inputs. Keep gray-mean changes separate from changes in the frequency-resolved arrays; do not silently insert cutoff-adjusted means into an unchanged vacuum-photon EOS.',
        acceptance='Input echo, material-ID mapping, warning absence and spectrum arithmetic are separate from mean consistency, physical plasma/EOS error, temperature interpolation and actual evolution. The outer-cell server failure stays failed.'))


def interval(token):
    value=Decimal(token);half=F(10)**value.as_tuple().exponent/2
    return F(value)-half,F(value)+half


def score(value,token):
    lo,hi=interval(token);return float(abs(F(value)-(lo+hi)/2)/((hi-lo)/2))


def parse(path):
    lines=path.read_text().splitlines();lines=[x.strip() for x in lines if x.strip()]
    counts=list(map(int,re.findall(r'=\s*(\d+)',lines[0])));assert len(counts)==3
    start=next(i for i,line in enumerate(lines) if line.startswith('No. Fraction'))+1;composition=[]
    while not lines[start].startswith('Temperature grid'):
        tokens=lines[start].split();assert len(tokens)==5
        composition.append(dict(number=tokens[0],mass=tokens[1],Z=int(tokens[2]),symbol=tokens[3],matid=int(tokens[4])))
        start+=1
    temperatures=lines[start+1].split();assert lines[start+2].startswith('Density grid')
    densities=lines[start+3].split();means={};spectra={};spectrum_tokens={};i=start+4
    while i<len(lines):
        if lines[i].startswith('Density') and 'T=' in lines[i]:
            T=lines[i].split('T=')[1].strip();tokens=lines[i+1].split();assert len(tokens)==5
            means[T]=tokens;i+=2
        elif lines[i].startswith('T (keV), Density'):
            T,rho=lines[i].split('=')[1].split();assert lines[i+1]=='Photon energy(keV), Total, Absorp, Scatt (cm**2/gm)'
            i+=2;tokens=[]
            while i<len(lines) and len(lines[i].split())==4 and re.fullmatch(r'[+\-.\dEe\s]+',lines[i]):
                tokens.append(lines[i].split());i+=1
            assert T not in spectra;spectra[T]=np.array(tokens,float);spectrum_tokens[T]=(rho,tokens)
        else:i+=1
    assert counts==[len(temperatures),len(densities),len(composition)]
    assert set(means)==set(spectra)==set(temperatures)
    return dict(composition=composition,temperatures=temperatures,densities=densities,
        means=means,spectra=spectra,tokens=spectrum_tokens,warnings=[x for x in lines if 'warning' in x.lower()])


def means_from_spectrum(a,T):
    x=a[:,0]/T;weightP=np.zeros_like(x);weightR=np.zeros_like(x);mask=x<700
    y=x[mask];em=np.exp(-y);den=-np.expm1(-y)
    weightP[mask]=15/np.pi**4*y**3*em/den
    weightR[mask]=15/(4*np.pi**4)*y**4*em/(den*den)
    P=a[:,2]*weightP;R=weightR/a[:,1]
    return np.array([[simpson(P,x=x),1/simpson(R,x=x)],
                     [np.trapz(P,x=x),1/np.trapz(R,x=x)]])


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    retrieval.verify();requests=json.loads((retrieval.OUT/'requests.json').read_text())
    results=json.loads((retrieval.OUT/'retrieval-result.json').read_text())['rows']
    selected={x['cell']:x for x in json.loads((retrieval.OUT/'selected-inputs.json').read_text())['rows']}
    assert len(requests)==len(results)
    from bs4 import BeautifulSoup
    material_text=BeautifulSoup((retrieval.original.OUT/'lanl-avmat.html').read_bytes(),'html.parser').get_text(' ',strip=True)
    expected_ids={int(z):int(new[-4:][::-1]) for z,symbol,new,old in re.findall(r'(\d+)\s+([A-Za-z]+)\s+(n\d+)\s+(n\d+)',material_text)}
    assert expected_ids[13]==9173 and expected_ids[2]==4675
    allrows=[];parsed={};finite=[]
    for request,result in zip(requests,results):
        assert request['name']==result['name'];name=request['name']
        if not result['retrieved']:
            allrows.append(dict(name=name,cell=request['cell'],retrieved=False,accepted=False));continue
        data=parse(retrieval.OUT/(name+'-table.txt'));parsed[name]=data
        fields=request['fields'];tokens=fields['mixture'].split();assert len(tokens)%3==0
        amounts={int(tokens[i+1]):(F(Decimal(tokens[i])),F(Decimal(tokens[i+2]))) for i in range(0,len(tokens),3)}
        number_sum=sum(x[0] for x in amounts.values());mass_sum=sum(x[0]*x[1] for x in amounts.values())
        echo=[];ids=[]
        assert {x['Z'] for x in data['composition']}==set(amounts)
        for row in data['composition']:
            number,mass=amounts[row['Z']]
            echo.extend([score(number/number_sum,row['number']),score(number*mass/mass_sum,row['mass'])])
            ids.append(dict(Z=row['Z'],returned=row['matid'],expected=expected_ids[row['Z']],passed=row['matid']==expected_ids[row['Z']]))
        submitted_T=sorted(F(Decimal(x)) for x in fields['temps'].split())
        returned_T=sorted(data['temperatures'],key=lambda x:Decimal(x));assert len(submitted_T)==len(returned_T)
        echo.extend(score(a,b) for a,b in zip(submitted_T,returned_T))
        assert len(data['densities'])==1;density=F(Decimal(fields['dens']))
        echo.append(score(density,data['densities'][0]));spectrum_rows=[]
        for T,a in data['spectra'].items():
            rho,tokens=data['tokens'][T];echo.append(score(density,rho));echo.append(score(density,data['means'][T][0]))
            sum_failures=[]
            for i,t in enumerate(tokens):
                lt,ut=interval(t[1]);la,ua=interval(t[2]);ls,us=interval(t[3])
                if lt>ua+us or ut<la+ls:sum_failures.append(i)
            physical_array=bool(np.all(np.isfinite(a)) and np.all(a[:,0]>0) and np.all(np.diff(a[:,0])>0)
                and np.all(a[:,1]>0) and np.all(a[:,2:]>=0) and np.ptp(a[:,2])>0)
            spectrum_rows.append(dict(T_keV=float(T),points=len(a),energy_range_keV=[float(a[0,0]),float(a[-1,0])],
                finite_positive_nonflat=physical_array,sum_printing_interval_failures=sum_failures))
            np.savez_compressed(OUT/f'{name}-{T}.npz',spectrum=a,gray=np.array(data['means'][T],float))
            if request['plasma_cutoff']=='off':
                values=means_from_spectrum(a,float(T));gray=np.array(data['means'][T][1:3],float)[::-1]
                errors=abs(values/gray-1);mutual=abs(values[0]/values[1]-1)
                finite.append(dict(name=name,T_keV=float(T),native_Planck_Rosseland=gray.tolist(),
                    Simpson_trapezoid_Planck_Rosseland=values.tolist(),relative_errors=errors.tolist(),
                    quadrature_mutual_difference=mutual.tolist(),
                    passed=bool(errors.max()<plan['mean_relative_tolerance'] and mutual.max()<plan['mean_relative_tolerance'])))
        accepted=bool(max(echo)<=1 and all(x['passed'] for x in ids) and not data['warnings']
            and all(x['finite_positive_nonflat'] and not x['sum_printing_interval_failures'] for x in spectrum_rows))
        allrows.append(dict(name=name,cell=request['cell'],retrieved=True,accepted=accepted,
            maximum_exact_input_echo_score=max(echo),material_ids=ids,warnings=data['warnings'],spectra=spectrum_rows,
            submitted_neutral_mass_per_baryon=float(mass_sum),
            opacity_conversion_to_per_baryon_gram=selected[request['cell']]['neutral_mass_per_baryon'],
            submitted_mean_atomic_mass=float(mass_sum/number_sum)))
    on=parsed['cell-5734-on'];off=parsed['cell-5734-off'];cutoff=[]
    for T,a in on['spectra'].items():
        b=off['spectra'][T];assert a.shape==b.shape
        cutoff.append(dict(T_keV=float(T),frequency_arrays_bitwise_equal=bool(np.array_equal(a,b)),
            gray_on_Rosseland_Planck=list(map(float,on['means'][T][1:3])),
            gray_off_Rosseland_Planck=list(map(float,off['means'][T][1:3]))))
    save('result.json',dict(classification='Counterexample candidate',completed=True,rows=allrows,
        retrieved_queries=sum(x['retrieved'] for x in allrows),accepted_queries=sum(x['accepted'] for x in allrows),
        cutoff_off_finite_mean_controls=finite,cutoff_off_finite_mean_controls_passed=all(x['passed'] for x in finite),
        central_cutoff_comparison=cutoff,material_id_convention_inferred=True,
        physical_EOS_or_opacity_certified=False,temperature_interpolation_certified=False,
        scattering_kernel_complete=False,outer_cell_resolved=False,full_GR_evolution=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify()


def verify():
    for name,key in [('plan.json','bindings'),('manifest.json','sha256')]:
        for rel,digest in json.loads((OUT/name).read_text())[key].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'result.json').read_text())['completed']
    print('PASS TOPS returned-state audit bindings; inspect acceptance and finite mean gates',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
