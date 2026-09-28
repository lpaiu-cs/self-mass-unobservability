def run():
    plan=bindings();assert json.loads((OUT/'controls.json').read_text())['passed'];p=Provider()
    state,_=d.state_data();saved=dict(np.load(d.OUT/'stellar-comparison.npz'));values=[];rows=[]
    for k,i in enumerate(saved['cells']):
        r,t,X=state['lnd'][i],state['lnT'][i],state['X'][i];b=d.mixture(p,r,t,X)['corrected'];values.append(b)
        assert np.array_equal(b[[0,1,3,4]],saved['corrected'][k,[0,1,3,4]]),int(i)
        if i in plan['controls']:
            for h in plan['finite_log_steps']:
                rm=d.mixture(p,r-h,t,X)['corrected'];rp=d.mixture(p,r+h,t,X)['corrected']
                tm=d.mixture(p,r,t-h,X)['corrected'];tp=d.mixture(p,r,t+h,X)['corrected']
                finite=np.array([-(tp[0]-tm[0])/(2*h),(rp[0]-rm[0])/(2*h),b[1]+(tp[1]-tm[1])/(2*h),
                    b[2]+(tp[2]-tm[2])/(2*h),b[2]+(rp[2]-rm[2])/(2*h)])
                expected=b[[1,2,4,5,6]];score=float(np.max(abs(finite-expected)/np.maximum(1,abs(expected))))
                rows.append(dict(cell=int(i),step=h,finite=finite.tolist(),expected=expected.tolist(),score=score,passed=score<plan['finite_derivative_tolerance']))
        if len(values)%512==0:print('SCREENING REPLAY',len(values),'/',len(saved['cells']),flush=True)
    values=np.array(values);np.savez_compressed(OUT/'corrected-stellar.npz',cells=saved['cells'],values=values,
        PDR_change=values[:,6]-saved['corrected'][:,6])
    save('result.json',dict(classification='Counterexample candidate',completed=True,cells=len(values),
        all_unaffected_outputs_bitwise=True,finite_controls=rows,finite_derivatives_passed=all(r['passed'] for r in rows),
        maximum_finite_score=max(r['score'] for r in rows),max_PDR_change=float(np.max(abs(values[:,6]-saved['corrected'][:,6]))),
        original_derivative_failure_preserved=True,physical_EOS_certified=False,native_EOS_replaced=False,full_GR_evolution=False))
    assert all(r['passed'] for r in rows)
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}));verify()
