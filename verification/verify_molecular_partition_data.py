"""Independently differentiate the saved spectral sum inside its interval box."""
import json, sys
import mpmath as mp
import molecular_partition_data as original

g=original.g;OUT=g.OUT/'gr-molecular-partition-independent'


def run():
    assert not OUT.exists();OUT.mkdir();original.verify();mp.mp.dps=90
    record=json.loads((original.OUT/'interval.json').read_text());levels=[]
    for line in (original.DATA/'eleroyh2.dat').read_text().splitlines():
        if line.strip():
            v,J,E=line.split();J=int(J)
            levels.append((mp.mpf(E),mp.mpf((2*J+1)*(1 if J%2==0 else 3))/4))
    c2=mp.mpf('6.62607015e-34')*299792458*100/mp.mpf('1.380649e-23')
    def logQ(lnT):return mp.log(mp.fsum(w*mp.exp(-c2*E/mp.exp(lnT)) for E,w in levels))
    rows=[]
    for T in [999999,1000000,1000001]:
        for order,box in enumerate(record['fixed_348_level_model']['logQ_derivatives']):
            value=mp.diff(logQ,mp.log(T),order)
            low,high=[mp.mpf(tuple(endpoint)) for endpoint in box['binary_endpoints']]
            assert low<=value<=high,(T,order,mp.nstr(value,30),low,high)
            rows.append(dict(T_K=T,order=order,value_90_digits=mp.nstr(value,90),contained=True))
    for row in record['rounded_coefficient_Taylor_model']:
        low,high=[mp.mpf(tuple(endpoint)) for endpoint in row['internal_Cv_over_kB']['binary_endpoints']]
        assert low<=high<0
    source=json.loads((original.OUT/'source-column-audit.json').read_text())
    assert source['min_Cp_minus_Cv']>20.78 and source['max_Cp_minus_Cv']<20.79
    result=dict(classification='Counterexample candidate',completed=True,rows=rows,all_contained=True,
        scope='90-digit numerical differentiation of the directly summed spectrum at three temperatures independently checks all 0--10 log-temperature jets. This is an implementation control, not a new continuum or physical error proof; the declared interval computation owns its mathematical enclosure. Source Cv is not silently treated as total ideal-molecule Cv.')
    (OUT/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    paths=[g.ROOT/'verification/verify_molecular_partition_data.py',original.OUT/'manifest.json',OUT/'result.json']
    (OUT/'manifest.json').write_text(json.dumps(dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths}),indent=2)+'\n')
    verify()


def verify():
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'result.json').read_text())['all_contained']
    print('PASS 33 independently differentiated spectral jets inside their declared interval enclosures',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
