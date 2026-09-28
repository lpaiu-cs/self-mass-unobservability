"""Reuse broad opacity groups; resolve only the remaining line-dominated band."""
from pathlib import Path
import argparse
import json
import time
import numpy as np
import def_photon_loss_bounds as previous

native=previous.native;ex=previous.ex;h=previous.h
OUT=previous.OUT.parent/'def-photon-line-resolved'


def prepare():
    assert not OUT.exists();OUT.mkdir()
    old=np.load(previous.OUT/'loss-0.npz');edge=old['bounds_keV']
    lo=int(np.searchsorted(edge,.0077)-1);hi=int(np.searchsorted(edge,.0135))
    windows=np.geomspace(edge[lo],edge[hi],9);requests=[]
    old_requests=json.loads((previous.OUT/'requests.json').read_text())[:2]
    for i,(a,b) in enumerate(zip(windows[:-1],windows[1:])):
        for j,row in enumerate(old_requests):
            name=f'phLine{i}T{j}';fields=dict(row['fields'],mixname=name,egplow=f'{a:.17g}',egphigh=f'{b:.17g}')
            requests.append(dict(name=name,cell=0,boundaries=1000,window=i,temperature_index=j,fields=fields))
    ex.write(OUT/'requests.json',requests);(OUT/'lanl-tops-form.html').write_bytes((native.OUT/'lanl-tops-form.html').read_bytes())
    paths=[Path(__file__),Path(previous.__file__),previous.OUT/'plan.json',previous.OUT/'result.json',
        previous.OUT/'loss-0.npz',previous.OUT/'loss-1.npz',OUT/'requests.json',OUT/'lanl-tops-form.html']
    forecast=16*max(r['seconds'] for r in json.loads((previous.OUT/'progress.json').read_text()))
    ex.write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='0f6b09ed',
        claim='Resolve the remaining line-dominated uncertainty using focused native groups while reusing every other spectral interval.',
        reassessment='The broad moment enclosure failed: maximum uniform loss/energy bounds 1.769/2.556 percent and gray Rosseland recombination 0.1043 percent. Over 87 percent of each uncertainty sum is in its largest 100 bins, concentrated at 7.8 to 13 eV. Replace only 7.7 to 13.5 eV with eight finer windows; no gray normalization or acceptance change.',
        replace_indices=[lo,hi],window_edges_keV=windows.tolist(),groups_per_window=999,
        gates=dict(native_gray_recombination=.001,uniform_complex_resolvent_relative_to_DC=.01,uniform_energy_response_absolute=.01),
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths},
        budget=dict(queries=16,forecast_network_seconds=forecast,hard_network_seconds=60,per_query_hard_seconds=20,CPU_workers=1,new_native_EOS_calls=0,new_stellar_steps=0,automatic_expansion=False),
        limits='Same declared positive-opacity moment assumptions; local loss propagator only, not a full scattering/emission operator or whole-star photon evolution.'))
    ex.write(OUT/'symbolic.json',previous.symbolic());print('FORECAST',forecast,flush=True)


def fetch():
    oldout=previous.OUT;previous.OUT=OUT;start=time.monotonic()
    try:previous.fetch()
    finally:previous.OUT=oldout
    assert time.monotonic()-start<60


def analyze():
    assert not (OUT/'result.json').exists();plan=json.loads((OUT/'plan.json').read_text());header=native.reader.namespace()['header'];rows=[[],[]]
    for request in json.loads((OUT/'requests.json').read_text()):
        path=OUT/(request['name']+'-table.txt');common=header(path,dict(request,fields=dict(request['fields'],datype='gray')))
        assert common['input_passed'];lines=[s.strip() for s in path.read_text().splitlines() if s.strip()]
        j=next(j for j,s in enumerate(lines) if s.startswith('Energy') and 'density =' in s);tokens=[s.split() for s in lines[j+1:j+1000]]
        group=np.array(tokens,float);assert group.shape==(999,3) and (group>0).all() and np.isfinite(group).all() and (group[:,1:]!=1e10).all()
        edge=np.geomspace(float(request['fields']['egplow']),float(request['fields']['egphigh']),1000)
        assert max(native.audit.score(native.audit.F(float(x)),t[0]) for x,t in zip(edge[:-1],tokens))<=1
        second=json.loads((OUT/(request['name']+'-results-request.json')).read_text())
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


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','fetch','analyze']);globals()[parser.parse_args().action]()
