"""Freeze and fetch small TOPS spectral queries at actual stellar states.

Counterexample candidate. Retrieval and physical acceptance are separate. The
service may substitute density-boundary values; never infer success from HTTP 200.
"""
from http.cookiejar import CookieJar
import json, shutil, sys, traceback, urllib.parse, urllib.request
import numpy as np
from bs4 import BeautifulSoup
import op_planck_reader_runner as op

g=op.g;reference=op.original.reference;OUT=g.OUT/'lanl-tops-stellar-probe'
BASE='https://aphysics2.lanl.gov/apps/'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def fields(form):
    data={}
    for element in form.find_all('input'):
        kind=element.get('type','').lower()
        if not element.get('name') or kind in ['submit','reset']:continue
        if kind in ['radio','checkbox'] and not element.has_attr('checked'):continue
        data[element['name']]=element.get('value','')
    for element in form.find_all('select'):
        option=element.find('option',selected=True) or element.find('option')
        data[element['name']]=option.get('value',option.get_text(strip=True))
    return data


def prepare():
    assert not OUT.exists();OUT.mkdir();op.verify()
    copied=['lanl-tops-form.html','lanl-avmat.html','lanl-opac-help.html','lanl-opactemps.html',
        'lanl-tops-al1keV-request.json','lanl-tops-al1keV-response.html','lanl-tops-al1keV-table-request.json',
        'lanl-tops-al1keV-table.html','lanl-tops-HeFlow.request.json','lanl-tops-HeFlow.summary.html',
        'lanl-tops-HeFlow.table-request.json','lanl-tops-HeFlow.table.html',
        'lanl-tops-AlOldFlow.request.json','lanl-tops-AlOldFlow.summary.html',
        'lanl-tops-AlOldFlow.table-request.json','lanl-tops-AlOldFlow.table.html']
    for name in copied:shutil.copy2(g.ROOT/'outputs'/name,OUT/name)
    form=BeautifulSoup((OUT/'lanl-tops-form.html').read_bytes(),'html.parser').find('form')
    defaults=fields(form);temperatures=np.array([float(x.get_text(strip=True)) for x in form.find('select',attrs={'name':'tlow'}).find_all('option')])
    state=dict(np.load(reference.OUT/'reference-state.npz'));requests=[];inputs=[]
    for i in [0,1175,2972,3043,4352,5734]:
        X=state['X'][i];Y=X/g.c.A;rhoB=np.exp(state['lnd'][i]);T=np.exp(state['lnT'][i])
        neutral_mass_per_baryon=float(Y@g.c.W);rho=rhoB*neutral_mass_per_baryon
        temperature_keV=T*1.380649e-16/(1.602176634e-9)
        j=int(np.searchsorted(temperatures,temperature_keV));assert 0<j<len(temperatures)
        pair=temperatures[j-1:j+1];parts=[];elements=[]
        for z in sorted(set(g.c.Z.astype(int))):
            selected=g.c.Z==z;amount=float(Y[selected].sum());assert amount>0
            mean_mass=float(Y[selected]@g.c.W[selected]/amount)
            parts.append(f'{amount:.12g} {z} {mean_mass:.10g}')
            elements.append(dict(Z=int(z),number_per_baryon=amount,mean_neutral_atomic_mass=mean_mass))
        mixture=' '.join(parts);assert len(mixture)<=450,(i,len(mixture))
        input_row=dict(cell=i,rho_B_g_cm3=float(rhoB),T_K=float(T),T_keV=float(temperature_keV),
            queried_temperatures_keV=pair.tolist(),queried_neutral_mass_density_g_cm3=float(rho),
            neutral_mass_per_baryon=neutral_mass_per_baryon,elements=elements,mixture=mixture)
        inputs.append(input_row)
        for cutoff in (['on','off'] if i==5734 else ['on']):
            name=f'cell-{i:04}-{cutoff}';data=dict(defaults,lib='new',fractype='atomic',isotope='isotope',
                mixture=mixture,mixname=name,tgrid='specific',temps=' '.join(format(x,'.12g') for x in pair),
                rgrid='specific',dens=format(rho,'.12g'),datype='cont',plasnu=cutoff)
            requests.append(dict(name=name,cell=i,plasma_cutoff=cutoff,fields=data))
    np.savez_compressed(OUT/'stellar-inputs.npz',**state)
    save('selected-inputs.json',dict(classification='Counterexample candidate',rows=inputs))
    save('requests.json',requests)
    paths=[g.ROOT/'verification/lanl_tops_stellar_probe.py',reference.OUT/'manifest.json',op.OUT/'manifest.json']+list(OUT.iterdir())
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='faefcb8',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths if p.is_file()},
        source_urls=[BASE,'https://aphysics2.lanl.gov/static/opacdocs/avmat.html',
            'https://aphysics2.lanl.gov/static/opacdocs/opac-help.html','https://aphysics2.lanl.gov/static/opacdocs/opactemps.html'],
        purpose='Query the new/ATOMIC option with Li,Be,B,F and actual high-density stellar mixtures absent from OP. Six selected cells, two tabulated temperatures bracketing each actual T, one exact requested neutral-mass density. The centre has an additional explicit plasma-cutoff-off diagnostic. Each request stays below the service six-state frequency-output limit.',
        normalization='Keep isotope number abundances Y_i=X_i/A_i. Group number Y_z=sum Y_i and mean neutral mass W_z=sum(Y_i W_i)/Y_z; send atomic fractions Y_z and isotopic weight W_z, rho_neutral=rho_B*sum(Y_i W_i). This preserves each element number density in the declared mass convention. A returned opacity per neutral gram converts to per baryon gram by multiplying sum(Y_i W_i). Grouping masses does not supply isotope-resolved line physics.',
        query_rounding='12 significant digits in atom-number fractions and density, 10 in mean neutral isotope masses; preserve these actual submitted strings and compare returned compositions to them before acceptance.',
        workflow='A fresh /submit calculation and its returned numerical-results form for EVERY query, in one cookie session. Changing only /results parameters is not a fresh physical calculation. Capture both requests and responses. Never count HTTP 200 alone as a passed scientific query.',
        acceptance='Follow-up analysis must validate all returned temperatures, densities, composition, material identities and WARNINGS. Reject density-boundary substitutions and the documented flat surrogate frequency spectra. Distinguish the provider plasma cutoff convention from a physical in-medium photon EOS or exact scattering kernel.',
        physical_EOS_or_opacity_certified=False,stellar_evolution_completed=False))


def bindings():
    for rel,digest in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel


def fetch(row):
    name=row['name'];opener=urllib.request.build_opener(urllib.request.HTTPCookieProcessor(CookieJar()))
    form=BeautifulSoup((OUT/'lanl-tops-form.html').read_bytes(),'html.parser').find('form')
    url=urllib.parse.urljoin(BASE,form['action']);assert url=='https://aphysics2.lanl.gov/submit'
    with opener.open(url,urllib.parse.urlencode(row['fields']).encode(),timeout=90) as response:
        raw=response.read();response_url=response.url
    (OUT/(name+'-summary.html')).write_bytes(raw)
    soup=BeautifulSoup(raw,'html.parser')
    results=next(f for f in soup.find_all('form') if f.find('input',attrs={'name':'output','value':'tabcol'}))
    second=fields(results);save(name+'-results-request.json',second)
    url=urllib.parse.urljoin(response_url,results['action']);assert url=='https://aphysics2.lanl.gov/results'
    with opener.open(url,urllib.parse.urlencode(second).encode(),timeout=90) as response:raw=response.read()
    (OUT/(name+'-table.html')).write_bytes(raw)
    soup=BeautifulSoup(raw,'html.parser');code=soup.find('code');assert code is not None
    text=code.get_text('\n',strip=True).replace('\xa0',' ')
    (OUT/(name+'-table.txt')).write_text(text)
    return dict(name=name,cell=row['cell'],retrieved=True,
        warnings_present='warning' in text.lower(),physical_acceptance_not_yet_evaluated=True)


def run():
    bindings();rows=[]
    for row in json.loads((OUT/'requests.json').read_text()):
        assert not list(OUT.glob(row['name']+'-*'))
        try:result=fetch(row)
        except Exception:
            result=dict(name=row['name'],cell=row['cell'],retrieved=False,traceback=traceback.format_exc())
        rows.append(result);save('retrieval-progress.json',dict(classification='Counterexample candidate',rows=rows))
        print('TOPS RETRIEVAL',result,flush=True)
    save('retrieval-result.json',dict(classification='Counterexample candidate',completed=True,rows=rows,
        physical_acceptance_not_yet_evaluated=True,full_GR_evolution=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify()


def verify():
    bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'retrieval-result.json').read_text())['completed']
    print('PASS TOPS stellar request/response bindings; acceptance needs separate analysis',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
