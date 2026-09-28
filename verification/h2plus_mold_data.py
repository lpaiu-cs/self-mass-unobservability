"""Audit all 423 retrieved H2+ level keys; keep model and physical errors distinct."""
from decimal import Decimal, localcontext
import json, shutil, sys, xml.etree.ElementTree as ET
import numpy as np
import mpmath as mp
import h2plus_spectral_data as previous

g=previous.g; base=previous.previous
DATA=g.OUT/'gr-h2plus-mold-source'; OUT=g.OUT/'gr-h2plus-mold-audit'
NMAX=[35,34,33,31,30,28,27,25,24,22,20,19,17,15,13,11,9,6,3,1]


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists() and not DATA.exists();OUT.mkdir();DATA.mkdir();previous.verify()
    names=['mold-index-http.html','mold-combo-ajax.js','mold-xsams.xml','mold-xsams-headers.txt',
           'mold-paper.html','zammit2017.pdf','curtin31657507-metadata.json','nist-factors2022.pdf']
    for name in names:shutil.copy2(g.ROOT/'outputs'/('h2plus-'+name),DATA/name)
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='99c551b18735ad5a73f350cb864f3e414a0c96e3',
        discovery='The HTTP XSAMS response and its 423 keyed energy states were inspected before this audit plan. This is a confirming source audit, not a blind preregistered discovery.',
        source=dict(classification='Imported from prior work',retrieved_local_date='2026-09-11',
            database='http://servo.aob.rs/mold/',paper='https://arxiv.org/abs/1603.08200',
            query='http://servo.aob.rs/mold/tap/sync?REQUEST=doQuery&LANG=VSS2&FORMAT=XSAMS&QUERY=select%20%2A%20where%20InchiKey%3D%27ZZIJOQHRUPVPQC-UHFFFAOYSA-N%27',
            XSAMS_physics_DOI='10.1051/0004-6361:20077206',
            comparison='Zammit et al. 2017, ApJ 851,64, DOI 10.3847/1538-4357/aa9712, PDF page 7 Tables 3 and 4 visually inspected. Nmax support and selected Q values manually transcribed. The full 423 energies here are a different MOL-D calculation, not a Zammit energy table.',
            comparison_pdf='https://ndownloader.figshare.com/files/62672668',
            comparison_pdf_metadata='https://api.figshare.com/v2/articles/31657507',
            units='XSAMS energies are au relative to one common auxiliary reference labelled ionized. Q uses E(v,N)-E(0,0), invariant to that common origin. The auxiliary label alone is not an absolute dissociation/ionization-energy certification.',
            transport='MOL-D HTTPS connection failed; the article and live UI advertise HTTP. HTTP source bytes and response headers are SHA-pinned, without authenticated-transport claims. No TLS verification was bypassed.',
            uncertainty='No level-specific Accuracy entries were supplied in this response. Digit count, agreement with another calculation and matching level support do not certify physical energy errors.'),
        expected_Nmax_by_v=NMAX,temperature_K=[1000,2000,3150,5040,8400,12600,16800,25200,100000,1000000,32000000],
        published_Q_Zammit2017=['26.1','66.2','142.1','361.2','1024.0','2043.1','3009.2','4546.3'],
        published_Q_comparison='Raw difference only; no fitted normalization, rejection threshold or claim of physical uncertainty.',
        hartree_J='4.3597447222060e-18',hartree_standard_uncertainty_J='4.8e-30',
        constants_source='https://physics.nist.gov/cuu/pdf/factors_2022.pdf',
        constants_scope='CODATA 2022 central Hartree value treated as an exact declared conversion in this model. Its quoted standard uncertainty is not a hard error bound and is not included in the interval certificate.',
        uniform_T_K=[999999,1000001],interval_digits=70,maximum_logT_derivative=10,
        printed_level_halfwidth_cm_inverse='0.00000003',
        interval_scope='All 423 retrieved positive fixed-spectrum terms per electron-spin state, normalized nuclear weights (2N+1)/4 even N and 3(2N+1)/4 odd N. The excitation printing box +/-3e-8 cm^-1 encloses two half-last-digit au bins after the declared exact conversion, not physical uncertainty. Electron spin 2 adds ln2 only to logQ. No excited-electronic/continuum/plasma occupation or implicit chemical-equilibrium EOS error is certified.',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in [
            g.ROOT/'verification/h2plus_mold_data.py',previous.OUT/'manifest.json',*[DATA/name for name in names]]}))


def parse():
    root=ET.parse(DATA/'mold-xsams.xml').getroot();molecules=root.findall('.//{*}Molecule')
    assert len(molecules)==1 and molecules[0].find('.//{*}InChI').text=='1S/H2/h1H/q+1'
    states=root.findall('.//{*}MolecularState');ids={s.attrib['stateID'] for s in states}
    assert len(ids)==len(states);rows={};origins=set();auxiliary=[]
    for state in states:
        node=state.find('.//{*}StateEnergy')
        if node is None:
            assert state.attrib['auxillary']=='true';auxiliary.append(state.attrib['stateID']);continue
        assert state.attrib['auxillary']=='false' and node.find('{*}Value').attrib['units']=='au'
        key=tuple(int(state.find('.//{*}'+n).text) for n in ['v','J']);assert key not in rows
        value=node.find('{*}Value').text;energy=Decimal(value);assert energy.is_finite() and energy<0
        rows[key]=dict(v=key[0],N=key[1],energy_au=value,state_id=state.attrib['stateID'])
        origins.add(node.attrib['energyOrigin'])
    assert len(auxiliary)==1 and origins==set(auxiliary)
    expected={(v,N) for v,n in enumerate(NMAX) for N in range(n+1)}
    assert set(rows)==expected and len(rows)==423
    refs=set();cross_sections=root.findall('.//{*}AbsorptionCrossSection');samples=0
    for section in cross_sections:
        ref=section.find('.//{*}StateRef').text;assert ref in ids and ref not in refs;refs.add(ref)
        assert section.find('.//{*}SpeciesRef').text==molecules[0].attrib['speciesID']
        assert [s.text for s in section.findall('{*}SourceRef')]==['BmolD-1','BmolD-6']
        lists=[section.find('{*}'+n+'/{*}DataList') for n in ['X','Y']]
        arrays=[np.fromstring(s.text,sep=' ') for s in lists]
        assert all(len(a)==int(s.attrib['count']) and np.all(np.isfinite(a)) for a,s in zip(arrays,lists))
        assert len(arrays[0])==len(arrays[1]) and np.all(np.diff(arrays[0])>0) and np.all(arrays[1]>=0)
        samples+=len(arrays[0])
    assert refs==ids-set(auxiliary) and len(cross_sections)==423
    return rows,dict(actual_energy_states=len(rows),auxiliary_reference_states=len(auxiliary),
        state_support_matches_all_Nmax=True,cross_sections=len(cross_sections),cross_section_samples=samples,
        all_state_references_resolve=True,physical_energy_accuracy_entries=len(root.findall('.//{*}Accuracy')))


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    levels,report=parse();rows=[]
    with localcontext() as ctx:
        ctx.prec=90
        factor=Decimal(plan['hartree_J'])/(Decimal('6.62607015e-34')*299792458*100)
        E0=Decimal(levels[0,0]['energy_au']);width0=Decimal(10)**E0.as_tuple().exponent/2
        max_printing=Decimal(0)
        for key,row in sorted(levels.items()):
            energy=Decimal(row['energy_au']);excitation=(energy-E0)*factor
            assert excitation>=0
            width=(width0+Decimal(10)**energy.as_tuple().exponent/2)*factor
            max_printing=max(max_printing,width);assert width<Decimal(plan['printed_level_halfwidth_cm_inverse'])
            row['excitation_cm_inverse']=str(excitation);rows.append((*key,str(excitation)))
    for v,n in enumerate(NMAX):
        assert all(Decimal(levels[v,N+1]['energy_au'])>Decimal(levels[v,N]['energy_au']) for N in range(n))
    old=json.loads((previous.OUT/'levels.json').read_text())['rows'];oldmap={(r['v'],r['N']):r for r in old}
    assert set(oldmap)<set(levels) and len(set(levels)-set(oldmap))==86
    comparison=[dict(v=v,N=N,mold_minus_Babb_excitation_cm_inverse=float(Decimal(levels[v,N]['excitation_cm_inverse'])-Decimal(r['excitation_cm_inverse']))) for (v,N),r in sorted(oldmap.items())]
    E=np.array([float(e) for v,N,e in rows]);weights=np.array([(2*N+1)*(1 if N%2==0 else 3)/2 for v,N,e in rows])
    missing=np.array([(v,N) not in oldmap for v,N,e in rows]);c2=6.62607015e-34*299792458*100/1.380649e-23
    comparisons=[]
    for j,T in enumerate(plan['temperature_K']):
        x=c2*E/T;w=weights*np.exp(-x);p=w/w.sum();mean=p@x;published=float(plan['published_Q_Zammit2017'][j]) if j<8 else None
        comparisons.append(dict(T_K=T,Q_with_electron_spin=float(w.sum()),DlogQ=float(mean),
            internal_Cv_over_kB=float(p@((x-mean)**2)),new_86_fraction_of_Q=float(w[missing].sum()/w.sum()),
            published_Q=published,Q_over_published_minus_one=float(w.sum()/published-1) if published else None))
    boxes=base.interval_audit(rows,base.coefficients(),plan)
    subset=boxes.pop('fixed_348_level_model');boxes['fixed_423_level_model']=subset;save('interval.json',boxes)
    mp.mp.dps=90;c2mp=mp.mpf('6.62607015e-34')*299792458*100/mp.mpf('1.380649e-23')
    args=[(mp.mpf(e),mp.mpf((2*N+1)*(1 if N%2==0 else 3))/4) for v,N,e in rows]
    def logQ(t):return mp.log(mp.fsum(w*mp.exp(-c2mp*e/mp.exp(t)) for e,w in args))
    controls=[]
    for T in [999999,1000000,1000001]:
        for n,box in enumerate(subset['logQ_derivatives']):
            value=mp.diff(logQ,mp.log(T),n);lo,hi=[mp.mpf(tuple(x)) for x in box['binary_endpoints']]
            assert lo<=value<=hi,(T,n);controls.append(dict(T_K=T,order=n,value=mp.nstr(value,90),contained=True))
    save('independent-controls.json',dict(classification='Counterexample candidate',rows=controls,all_contained=True))
    save('levels.json',dict(classification='Imported from prior work',rows=list(levels.values())))
    save('overlap-comparison.json',dict(classification='Counterexample candidate',rows=comparison,
        boundary='Different calculations and energy origins. Differences of excitation energies are model discrepancies, not certified physical error bars.'))
    save('result.json',dict(classification='Counterexample candidate',completed=True,**report,
        recovered_missing_level_keys=[list(k) for k in sorted(set(levels)-set(oldmap))],
        max_excitation_printing_scenario_halfwidth_cm_inverse=str(max_printing),
        maximum_absolute_overlap_excitation_difference_cm_inverse=max(abs(r['mold_minus_Babb_excitation_cm_inverse']) for r in comparison),
        comparisons=comparisons,independent_controls=33,all_independent_controls_contained=True,
        physical_EOS_certified=False,native_EOS_replaced=False,full_GR_evolution=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify();print('MOLD FULL SUPPORT',report,comparisons[-2],flush=True)


def verify():
    for name,key in [('plan.json','bindings'),('manifest.json','sha256')]:
        for rel,digest in json.loads((OUT/name).read_text())[key].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    _,audit=parse();result=json.loads((OUT/'result.json').read_text())
    assert result['completed'] and result['all_independent_controls_contained']
    assert audit['actual_energy_states']==423 and len(result['recovered_missing_level_keys'])==86
    print('PASS actual MOL-D 423-level support, provenance, positive model and independent jet controls',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
