def run():
    plan=bindings();d.verify();a=dict(np.load(d.OUT/'stellar-comparison.npz'));state,native=d.state_data();p=d.Provider();rows=[]
    for k,i in enumerate(a['cells']):
        X=state['X'][i];count=X/g.c.A;x=count/count.sum();z=g.c.Z
        rs,ge=a['parameters'][k,:2];ideal=np.zeros(7);rawql=np.zeros(7);domain=[]
        for j in np.flatnonzero(x>0):
            gamma=ge*z[j]**(5/3);theta=gamma/np.sqrt(rs)*np.sqrt(3/(1822.88848*g.c.A[j]))/z[j]**(7/6)
            R=3*(gamma/theta)**2;fid=1.5*np.log(theta*theta/gamma)-float(np.float32(1.323515))
            ideal+=x[j]*np.array([fid,1.5,1.,1.5-fid,1.5,1.,1.])
            rawql+=x[j]*d.seven(p.call('liqubc',[float(R),float(theta)],6))
            domain.append([float(R),float(theta),float(gamma)])
        groupdifference=a['raw_groups'][k,1]-a['raw_groups'][k,0]-ideal
        expected=rawql.copy();expected[6]=rawql[5]
        error=float(np.max(abs(groupdifference-expected)/np.maximum(1,abs(expected))))
        rawtotal=a['raw_groups'][k,1]+a['mixing'][k];correcttotal=a['corrected'][k]+ideal
        dropped=(x>0)&(x<1e-7)
        rows.append(dict(cell=int(i),return_group_identity_score=error,
            quantum_F_omitted_from_first_group=float(a['quantum'][k,0]),
            raw_quantum_PDR_error=float(rawql[6]-a['quantum'][k,6]),
            PDR2_uses_PDTQL_error=float(rawql[5]-rawql[6]),
            total_ion_nonideal_PDR_difference=float(rawtotal[6]-correcttotal[6]),
            total_ion_nonideal_PDR_relative_difference=float((rawtotal[6]-correcttotal[6])/max(1e-10,abs(correcttotal[6]))),
            trace_isotopes_below_MELANGE_cutoff=int(dropped.sum()),
            trace_number_fraction=float(x[dropped].sum()),trace_baryon_fraction=float(X[dropped].sum()),
            R_range=[float(np.min(domain,axis=0)[0]),float(np.max(domain,axis=0)[0])],
            theta_range=[float(np.min(domain,axis=0)[1]),float(np.max(domain,axis=0)[1])],
            gamma_range=[float(np.min(domain,axis=0)[2]),float(np.max(domain,axis=0)[2])]))
    passed=all(r['return_group_identity_score']<plan['identity_score_tolerance'] for r in rows)
    save('return-audit.json',dict(classification='Counterexample candidate',cells=len(rows),source_identity_passed=passed,rows=rows,
        maximum_identity_score=max(r['return_group_identity_score'] for r in rows),
        max_quantum_F_omitted=max(r['quantum_F_omitted_from_first_group'] for r in rows),
        max_total_ion_nonideal_PDR_relative_difference=max(abs(r['total_ion_nonideal_PDR_relative_difference']) for r in rows),
        diagnosis='The original excess group excludes liquid quantum corrections. Second minus first minus ideal equals LIQUBC except PDR uses PDTQL. Exact free-energy derivatives fix this separate mixed evaluator; source originals stay unchanged.'))
    assert passed
    save('result.json',dict(classification='Counterexample candidate',completed=True,source_contract_audit_passed=passed,
        source_result_sha256=g.c.sha(d.OUT/'manifest.json'),
        calibration_replay_passed=json.loads((OUT/'calibration.json').read_text())['calibration_replay_passed'],
        physical_EOS_certified=False,full_GR_evolution=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}));verify()
