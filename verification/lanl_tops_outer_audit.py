"""Accept separately retrieved outer spectra only after input and mean checks."""
import json, sys
import numpy as np
import lanl_tops_outer_spectra as retrieval
import lanl_tops_cutoff_reader as reader

g=retrieval.g;audit=retrieval.control.audit;OUT=g.OUT/'lanl-tops-outer-audit'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();retrieval.verify();reader.verify()
    paths=[g.ROOT/'verification/lanl_tops_outer_audit.py',retrieval.OUT/'manifest.json',reader.OUT/'manifest.json']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='f305247',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        acceptance='Same exact decimal input, finite positive/nonflat spectrum and total=absorption+scattering printing interval gates as earlier. Gray means must be bitwise equal to the separately retrieved outer gray table at each temperature. Vacuum mean quadrature against these ON spectra is a finite diagnostic because gray ON/OFF agree at this outer state. Require Simpson/trapezoid means and mutual differences below the same 1e-3 relative gate; preserve any failure.',
        mean_relative_tolerance=1e-3,
        boundary='Successful single-temperature queries do not erase failed paired queries or certify an LTE atmosphere. No parameter changes and no high-density cutoff reinterpretation.'))


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    requests=json.loads((retrieval.OUT/'requests.json').read_text());rows=[]
    result=json.loads((retrieval.OUT/'retrieval-result.json').read_text())['rows'];header=reader.namespace()['header']
    gray=next(x for x in json.loads((reader.OUT/'result.json').read_text())['rows'] if x['name']=='outergrayon')['means']
    for request,outcome in zip(requests,result):
        assert request['name']==outcome['name'];name=request['name']
        if not outcome['retrieved']:rows.append(dict(name=name,retrieved=False,accepted=False));continue
        path=retrieval.OUT/(name+'-table.txt');common=header(path,request);data=audit.parse(path)
        assert len(data['spectra'])==1;t,a=next(iter(data['spectra'].items()));rho,tokens=data['tokens'][t]
        nonnegative=bool(np.all(np.isfinite(a)) and np.all(a[:,0]>0) and np.all(np.diff(a[:,0])>0)
            and np.all(a[:,1]>0) and np.all(a[:,2:]>=0) and np.ptp(a[:,2])>0)
        failures=[]
        for i,tok in enumerate(tokens):
            lt,ut=audit.interval(tok[1]);la,ua=audit.interval(tok[2]);ls,us=audit.interval(tok[3])
            if lt>ua+us or ut<la+ls:failures.append(i)
        density_score=audit.score(audit.F(audit.Decimal(request['fields']['dens'])),rho)
        gray_equal=common['means'][t]==gray[t];values=audit.means_from_spectrum(a,float(t))
        native=np.array(gray[t][1:3])[::-1];error=abs(values/native-1);mutual=abs(values[0]/values[1]-1)
        finite=bool(error.max()<plan['mean_relative_tolerance'] and mutual.max()<plan['mean_relative_tolerance'])
        accepted=bool(common['input_passed'] and density_score<=1 and nonnegative and not failures and gray_equal and finite)
        np.savez_compressed(OUT/(name+'.npz'),spectrum=a,gray=np.array(gray[t]))
        rows.append(dict(name=name,retrieved=True,accepted=accepted,points=len(a),input_passed=common['input_passed'],
            input_max_score=max(density_score,common['input_max_score']),finite_positive_nonflat=nonnegative,
            sum_printing_interval_failures=failures,gray_matches_independent_request=gray_equal,
            integrated_Planck_Rosseland=values.tolist(),relative_errors=error.tolist(),mutual_difference=mutual.tolist(),
            finite_mean_gate_passed=finite))
    save('result.json',dict(classification='Counterexample candidate',completed=True,rows=rows,
        accepted_queries=sum(x['accepted'] for x in rows),original_pair_failures_preserved=True,
        physical_LTE_atmosphere_certified=False,full_stellar_coverage=False,full_GR_evolution=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify()


def verify():
    for name,key in [('plan.json','bindings'),('manifest.json','sha256')]:
        for rel,digest in json.loads((OUT/name).read_text())[key].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'result.json').read_text())['completed']
    print('PASS outer spectral audit bindings; inspect finite acceptance, no LTE atmosphere certificate',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
