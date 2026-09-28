"""Extract the actual published H2+ support; do not call a subset a full spectrum."""
from decimal import Decimal
import gzip, io, json, sys
import numpy as np
import mpmath as mp
import molecular_partition_data as previous

g=previous.g;DATA=g.OUT/'gr-h2plus-spectral-source';OUT=g.OUT/'gr-h2plus-spectral-audit'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();previous.verify()
    raw=gzip.decompress((DATA/'table1.dat.gz').read_bytes());assert len(raw.splitlines())==284168
    files={'ReadMe.txt':'https://cdsarc.cds.unistra.fr/ftp/J/ApJS/216/21/ReadMe',
        'index.html':'https://cdsarc.cds.unistra.fr/ftp/J/ApJS/216/21/',
        'table1.dat.gz':'https://cdsarc.cds.unistra.fr/ftp/J/ApJS/216/21/table1.dat.gz',
        'paper.html':'https://arxiv.org/html/1412.2606v1',
        'errata.html':'https://lweb.cfa.harvard.edu/~babb/errata.html'}
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='870e0c8',
        source=dict(classification='Imported from prior work',author='James F. Babb',year=2015,
            doi='10.1088/0067-0049/216/1/21',cds='J/ApJS/216/21',retrieved_date='2026-09-11',urls=files,
            erratum='ApJS 237,20 (2018), author errata: the prefactor in Eq.9 is 7.80e-26. The original 1.475e-20 is not used to compute a radiative-association rate.',
            parent_claimed_bound_levels=423,
            table_selection='Only dipole-matrix-element squared values greater than 1e-6 are listed. A transition-selected table need not contain every bound level. The level count must be measured.',
            energies='The E(v,N) column is the positive magnitude of binding to the H(1s)+H+ dissociation limit in cm^-1, not the excitation energy above (0,0). Excitation term value is D00-DvN.',
            approximation='Born-Oppenheimer potential; adiabatic, nonadiabatic, relativistic and radiative corrections omitted. The reported count is the parent claim, not an independently certified complete spectrum.'),
        expected_table_rows=284168,temperature_K=[1000,5000,20000,100000,999999,1000000,1000001,32000000],
        uniform_T_K=[999999,1000001],interval_digits=70,maximum_logT_derivative=10,
        printed_level_halfwidth_cm_inverse='0.01',
        interval_scope='The positive fixed subset sum per electron-spin state, with normalized even/odd nuclear weights 1/4 and 3/4. Each excitation energy is a difference of two printed binding energies, so +/-0.01 cm^-1 covers a scenario of +/-0.005 for each printed entry. This is not a physical energy-uncertainty bound. A constant electron-spin degeneracy 2 adds ln(2) to log Q; derivatives of orders >=1 are unchanged. Missing levels, plasma occupation and the complete EOS are excluded.',
        source_gzip_crc_verified=True,
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in [
            g.ROOT/'verification/h2plus_spectral_data.py',previous.OUT/'manifest.json',*[DATA/name for name in files]]}))


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    text=gzip.decompress((DATA/'table1.dat.gz').read_bytes()).decode('ascii')
    table=np.loadtxt(io.StringIO(text));assert table.shape==(plan['expected_table_rows'],6)
    assert np.all(np.isfinite(table)) and np.all(table[:,2:]>=0)
    levels={};counts={}
    for line in text.splitlines():
        v,N=int(line[:2]),int(line[3:5]);binding=Decimal(line[17:25].strip());key=(v,N)
        if key in levels:assert levels[key]==binding,key
        else:levels[key]=binding
        counts[key]=counts.get(key,0)+1
    D00=levels[0,0];assert D00==max(levels.values())
    rows=[(v,N,str(D00-binding)) for (v,N),binding in sorted(levels.items())]
    assert all(Decimal(E)>=0 for v,N,E in rows)
    support=[dict(v=v,N=N,binding_cm_inverse=str(levels[v,N]),excitation_cm_inverse=E,
        tabulated_transition_rows=counts[v,N]) for v,N,E in rows]
    save('levels.json',dict(classification='Imported from prior work',rows=support))
    energy=np.array([float(E) for v,N,E in rows]);nuclear=np.array([(2*N+1)*(1 if N%2==0 else 3)/4 for v,N,E in rows])
    c2=6.62607015e-34*299792458*100/1.380649e-23;native=previous.native_function();comparisons=[]
    for T in plan['temperature_K']:
        x=c2*energy/T;weights=nuclear*np.exp(-x);p=weights/weights.sum();q1=p@x
        comparisons.append(dict(T_K=T,subset_Q_with_electron_spin=float(2*weights.sum()),
            subset_DlogQ=float(q1),subset_internal_Cv_over_kB=float(p@((x-q1)**2)),
            native_H2plus_logQ_jets=native(T)[1].tolist()))
    boxes=previous.interval_audit(rows,previous.coefficients(),plan)
    subset=boxes.pop('fixed_348_level_model');boxes['fixed_actual_subset_model']=subset
    boxes['number_of_included_levels']=len(rows);save('interval.json',boxes)
    mp.mp.dps=90
    energies=[mp.mpf(E) for v,N,E in rows];weights=[mp.mpf((2*N+1)*(1 if N%2==0 else 3))/4 for v,N,E in rows]
    factor=mp.mpf('6.62607015e-34')*299792458*100/mp.mpf('1.380649e-23')
    def logQ(t):return mp.log(mp.fsum(w*mp.exp(-factor*E/mp.exp(t)) for E,w in zip(energies,weights)))
    controls=[]
    for T in [999999,1000000,1000001]:
        for n,box in enumerate(subset['logQ_derivatives']):
            value=mp.diff(logQ,mp.log(T),n);low,high=[mp.mpf(tuple(x)) for x in box['binary_endpoints']]
            assert low<=value<=high,(T,n)
            controls.append(dict(T_K=T,order=n,contained=True,value=mp.nstr(value,90)))
    save('independent-controls.json',dict(classification='Counterexample candidate',rows=controls,all_contained=True))
    actual_count=len(rows);claimed=plan['source']['parent_claimed_bound_levels']
    save('result.json',dict(classification='Counterexample candidate',completed=True,
        actual_transition_rows=len(table),actual_unique_levels=actual_count,
        all_repeated_level_energies_identical=True,v_range=[min(v for v,N in levels),max(v for v,N in levels)],
        N_range=[min(N for v,N in levels),max(N for v,N in levels)],
        levels_per_N={str(N):sum(n==N for v,n in levels) for N in sorted({N for v,N in levels})},
        subset_shortfall_relative_to_parent_423_claim=claimed-actual_count,
        all_parent_claimed_levels_supplied=actual_count==claimed,
        binding_ground_cm_inverse=str(D00),comparisons=comparisons,
        no_full_spectrum_or_physical_EOS_certificate=True,no_native_EOS_replacement=True,no_GR_evolution=True))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify();print('H2PLUS SUPPORT',actual_count,'of parent claim',claimed,'D00',D00,comparisons[-3],flush=True)


def verify():
    for name,key in [('plan.json','bindings'),('manifest.json','sha256')]:
        for rel,digest in json.loads((OUT/name).read_text())[key].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'result.json').read_text())['completed']
    assert json.loads((OUT/'independent-controls.json').read_text())['all_contained']
    print('PASS actual H2+ spectral support and conditional jet controls; missing levels remain explicit',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
