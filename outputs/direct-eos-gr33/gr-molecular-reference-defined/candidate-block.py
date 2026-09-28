def block(start):
    plan=bindings();state=dict(np.load(g.OUT/'initial-state-17-4.npz'));eos=model.EOS();old=g.EOS()
    stop=min(start+plan['block_size'],plan['cells']);rows=[];values=[];previous=[];molecules=[]
    assert not (OUT/f'block-{start:04}.json').exists()
    for i in range(start,stop):
        r,t,X=state['lnd'][i],state['lnT'][i],state['X'][i]
        snap=eos.snapshot(r,t,X);report=model.previous.switch.molecular.p.c.s.check(snap,X,r,eos)
        a=snap['eos'];b=old(2,r,t,X);assert a.shape==b.shape==(21,)
        passed=not report['missing_nonzero_elements'] and report['inventory_error']<plan['inventory_tolerance'] and report['charge_error']<plan['charge_tolerance'] and a[0]>0 and a[1]>0 and a[10]>0
        rows.append(dict(cell=i,**report,passed=bool(passed)));values.append(a);previous.append(b)
        molecules.append(snap['molecular_H_fractions'])
    target=OUT/f'block-{start:04}.npz'
    np.savez_compressed(target,cells=np.arange(start,stop),new=np.array(values),old=np.array(previous),molecules=np.array(molecules))
    record=dict(classification='Counterexample candidate',start=start,stop=stop,rows=rows,all_passed=all(r['passed'] for r in rows),
        plan_sha256=g.c.sha(OUT/'plan.json'),output_sha256=g.c.sha(target))
    save(f'block-{start:04}.json',record);print('NEW MOLECULAR REFERENCE',start,stop,record['all_passed'],flush=True)
    return record
