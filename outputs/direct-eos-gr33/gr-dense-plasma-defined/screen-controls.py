def controls():
    plan=bindings();p=Provider();old=d.Provider();a=np.load(d.OUT/'stellar-comparison.npz');rows=[]
    x=sp.symbols('x');u=sp.Function('u')(x);v=sp.Function('v')(x)
    correct=(sp.diff(u,x,2)-(2*sp.diff(u,x)*sp.diff(v,x)+u*sp.diff(v,x,2))/v+2*u*(sp.diff(v,x)/v)**2)/v
    assert sp.simplify(sp.diff(u/v,x,2)-correct)==0
    source=(sp.diff(u,x,2)-(2*sp.diff(u,x)*sp.diff(v,x)+u*sp.diff(v,x,2))/v+2*(sp.diff(v,x)/v)**2)/v
    assert sp.simplify(source-correct-2*(1-u)*sp.diff(v,x)**2/v**3)==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        correction='(U/D)_xx=U_xx/D-(2U_x D_x+U D_xx)/D^2+2U D_x^2/D^3.',
        original_error='The distributed COR0DXX differs by 2*(1-U0)*D0_x^2/D0^3. Only the density second derivative chain is changed.'))
    for i in plan['controls']:
        k=int(np.flatnonzero(a['cells']==i)[0]);rs,ge=a['parameters'][k,:2]
        for z in plan['charges']:
            before=d.seven(old.call('fscrliq8',[rs,ge,float(z)],6));after=d.seven(p.call('fscrliq8',[rs,ge,float(z)],6))
            ref=independent(rs,ge,z);same=np.array_equal(before[[0,1,3,4]],after[[0,1,3,4]]);score=float(np.max(abs(after-ref)/np.maximum(1,abs(ref))))
            rows.append(dict(cell=i,Z=z,before=before.tolist(),after=after.tolist(),independent=ref.tolist(),
                unaffected_bitwise=same,corrected_score=score,passed=bool(same and score<plan['independent_score_tolerance'])))
    save('controls.json',dict(classification='Counterexample candidate',rows=rows,passed=all(r['passed'] for r in rows),
        maximum_score=max(r['corrected_score'] for r in rows)))
    assert all(r['passed'] for r in rows)
    print('SCREENING MP controls',len(rows),max(r['corrected_score'] for r in rows),flush=True)
