def analyze():
    assert not (OUT/'result.json').exists();plan=json.loads((OUT/'plan.json').read_text());header=native.reader.namespace()['header'];rows=[[],[]]
    for request in json.loads((OUT/'requests.json').read_text()):
        path=OUT/(request['name']+'-table.txt');common=header(path,dict(request,fields=dict(request['fields'],datype='gray')))
        assert common['input_passed'];lines=[s.strip() for s in path.read_text().splitlines() if s.strip()]
        j=next(j for j,s in enumerate(lines) if s.startswith('Energy') and 'density =' in s);tokens=[s.split() for s in lines[j+1:j+1000]]
        group=np.array(tokens,float);assert group.shape==(999,3) and (group>0).all() and np.isfinite(group).all() and (group[:,1:]!=1e10).all()
        edge=np.geomspace(float(request['fields']['egplow']),float(request['fields']['egphigh']),1000)
        assert max(native.audit.score(native.audit.F(float(x)),t[0]) for x,t in zip(edge[:-1],tokens))<=1
        second=json.loads((OUT/(request['name']+'-generated-results-request.json')).read_text())
        assert float(second['egplow'])==edge[0] and float(second['egphigh'])==edge[-1]
        rows[request['temperature_index']].append((edge,group[:,1:]))
    results=[]
    for j,blocks in enumerate(rows):
        old=np.load(previous.OUT/f'loss-{j}.npz');original=np.load(native.OUT/(['surf1000T0015.npz','surf1000T002.npz'][j]));T=float(old['T_keV']);lo,hi=plan['replace_indices']
        parts=[previous.bounds(edge,group,T,original['spectrum']) for edge,group in blocks]
        edge=np.concatenate([old['bounds_keV'][:lo+1]]+[a[0][1:] for a in blocks]+[old['bounds_keV'][hi+1:]])
        groups=np.concatenate([old['groups'][:lo]]+[b for _,b in blocks]+[old['groups'][hi:]])
        w=np.concatenate([old['weights'][:lo]]+[p['w'] for p in parts]+[old['weights'][hi:]])
        inverse=np.concatenate([old['resolvent_uncertainty'][:lo]]+[p['inverse'] for p in parts]+[old['resolvent_uncertainty'][hi:]])
        energy=np.concatenate([old['energy_response_uncertainty'][:lo]]+[p['energy'] for p in parts]+[old['energy_response_uncertainty'][hi:]])
        assert len(edge)==len(groups)+1 and np.all(np.diff(edge)>0)
        dc=w[:,1]@(1/groups[:,0]);means=np.array([w[:,0]@groups[:,1]/w[:,0].sum(),w[:,1].sum()/dc]);comparison=abs(means/original['gray'][[2,1]]-1)
        inv=float(w[:,1]@inverse/dc);en=float(w[:,1]@energy);passed=bool(comparison.max()<.001 and inv<.01 and en<.01)
        np.savez_compressed(OUT/f'loss-{j}.npz',bounds_keV=edge,groups=groups,weights=w,loss_rate_per_opacity=groups[:,0],
            resolvent_uncertainty=inverse,energy_response_uncertainty=energy,T_keV=T,rho_neutral_g_cm3=original['gray'][0],passed=passed)
        results.append(dict(T_keV=T,groups=len(groups),passed=passed,gray_relative=comparison.tolist(),
            uniform_complex_resolvent_relative_to_DC=inv,uniform_energy_response_absolute=en))
    result=dict(classification='Counterexample candidate',passed=all(r['passed'] for r in results),checks=results,
        local_loss_propagator_constructed=True,angular_gain_kernel_complete=False,material_energy_exchange_closed=False,whole_star_heat_closed=False,full_dynamic_charge_solved=False)
    ex.write(OUT/'result.json',result);print('RESULT',result,flush=True)
