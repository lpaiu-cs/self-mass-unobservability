"""Positive loss-resolvent surrogate enclosed by native opacity moments.

The frequency distribution is not inferred uniquely. Arithmetic/harmonic
moments and a declared scattering floor bound its possible response.
This is a loss propagator, not an angular redistribution or energy closure.
"""
from pathlib import Path
from types import FunctionType
import argparse
import json
import time
import numpy as np
from numpy.polynomial.legendre import leggauss
import def_photon_native_groups as native

model=native.model;ex=model.ex;h=model.h
OUT=native.OUT.parent/'def-photon-loss-bounds'


def weights(edges,T,order=24):
    gx,gw=leggauss(order);lo=edges[:-1]/T;hi=edges[1:]/T
    u=(hi[:,None]+lo[:,None])/2+(hi-lo)[:,None]*gx/2
    em=np.exp(-u);den=-np.expm1(-u);P=15/np.pi**4*u**3*em/den
    R=P*u/(4*den)
    return np.c_[P@gw,R@gw]*(hi-lo)[:,None]/2


def maxima(f,a,c,power):
    # The logarithmic derivative is strictly decreasing from power to power-3.
    lo=np.full_like(a,-80.);hi=np.full_like(a,80.)
    for _ in range(90):
        middle=(lo+hi)/2;r=np.exp(middle)
        positive=power-r/(f+r)-r/(a+r)-r/(c+r)>0
        lo=np.where(positive,middle,lo);hi=np.where(positive,hi,middle)
    rl=np.nextafter(np.exp(lo),0);rh=np.nextafter(np.exp(hi),np.inf)
    return rh**power/((f+rl)*(a+rl)*(c+rl))


def bounds(edges,group,T,spectrum,tail=False):
    w=weights(edges,T);tight=weights(edges,T,48)
    assert np.max(abs(w-tight))<1e-12
    # Five-digit printed opacities: relative 5e-5 encloses half a last-place unit.
    Hlo=1/(group[:,0]*(1+5e-5));Hhi=1/(group[:,0]*(1-5e-5));H=1/group[:,0]
    if tail:return dict(w=w,H=H,inverse=2*Hhi,energy=np.full(len(H),2.))
    selected=(spectrum[:,0]>=edges[0])&(spectrum[:,0]<=edges[-1])
    near=np.r_[np.flatnonzero(selected),np.searchsorted(spectrum[:,0],edges[[0,-1]])+np.array([-1,0])]
    scattering=spectrum[np.clip(near,0,len(spectrum)-1),3]
    f=float(scattering.min()*(1-5e-5));smax=float(scattering.max()*(1+5e-5))
    u=edges[1:]/T;q=u/(-np.expm1(-u))
    M=group[:,1]*(1+5e-5)*w[:,0]/w[:,1]*q/4+smax
    assert np.all(M*Hhi>=1),'Inconsistent native moment bounds'
    Hlo=np.maximum(Hlo,1/M) # Remove impossible harmonic moments, not data.
    a=1/Hhi;c=(M-f)/(1-f*Hlo)
    assert np.all(a>f) and np.all(c>f)
    coefficient=M-a
    # For any admissible H and M, the real-axis Jensen/Radau gap is bounded
    # by coefficient*r/((f+r)*(a+r)*(c+r)). Cauchy gives <=2*sqrt(2) times
    # this gap on the imaginary axis. Multiply by r for energy-response error.
    inverse=2*np.sqrt(2)*coefficient*maxima(f,a,c,1)+np.maximum(Hhi-H,H-Hlo)
    energy=2*np.sqrt(2)*coefficient*maxima(f,a,c,2)+np.maximum(abs(1/Hlo-group[:,0])/(1/Hlo+group[:,0]),abs(a-group[:,0])/(a+group[:,0]))
    return dict(w=w,H=H,inverse=np.minimum(inverse,2*Hhi),energy=np.minimum(energy,2.),floor=f,arithmetic_upper=M)


def symbolic():
    import sympy as sp
    f,a,c,r=sp.symbols('f a c r',positive=True)
    v=(a-f)*(c-a);p=(a-f)/(c-f)
    radau=(1-p)/(f+r)+p/(c+r)
    assert sp.factor(radau-1/(a+r)-v/((f+r)*(a+r)*(c+r)))==0
    x=sp.symbols('x',positive=True)
    assert sp.factor(1/(x+r)-1/(a+r)+(x-a)/(a+r)**2-(x-a)**2/((a+r)**2*(x+r)))==0
    return dict(classification='Proven',passed=True,
        assumptions='Positive opacity k>=f>0 under a normalized spectral measure, harmonic moment H=E(1/k), arithmetic moment M=E(k).',
        reweight='dnu=dmu/(k H); a=1/H, Var_nu=M/H-a^2, c=(M-f)/(1-f H). The Radau nodes f,c match the first two nu moments.',
        real_gap='For r>0, 0 <= 1/(a+r)-E[1/(k+r)] <= r*(M-a)/((f+r)*(a+r)*(c+r)).',
        complex_gap='For s=i*w, the absolute error is bounded by 2*sqrt(2) times the real-axis gap at r=abs(w). This follows from abs(k+i*w)>=(k+abs(w))/sqrt(2) and (a+r)^2/(a^2+r^2)<=2.',
        uniform='Maximize r^p/((f+r)*(a+r)*(c+r)), p=1 or 2, at its unique positive derivative root. p=1 bounds the loss resolvent; p=2 bounds 1-s*resolvent. The probability-weighted sum of per-group bounds is uniform in real temporal frequency.',
        scope='Exact conditional moment inequalities; tabulated inputs, rounded arithmetic, scattering interpolation and continuum physical-model error remain separate.')


def prepare():
    assert not OUT.exists();OUT.mkdir();assert json.loads((native.OUT/'result.json').read_text())['passed']
    data=np.load(native.OUT/'surf1000T0015.npz');coarse=data['bounds_keV'];endpoints=coarse[[435,847]]
    windows=np.geomspace(*endpoints,17);original=json.loads((native.OUT/'requests.json').read_text());requests=[]
    for index,(lo,hi) in enumerate(zip(windows[:-1],windows[1:])):
        for j,row in enumerate(original[2:]):
            name=f'phW{index:02}T{j}';fields=dict(row['fields'],mixname=name,datype='groups',egplow=f'{lo:.17g}',egphigh=f'{hi:.17g}')
            requests.append(dict(name=name,cell=0,boundaries=1000,window=index,temperature_index=j,fields=fields))
    ex.write(OUT/'requests.json',requests);(OUT/'lanl-tops-form.html').write_bytes((native.OUT/'lanl-tops-form.html').read_bytes())
    paths=[Path(__file__),Path(native.__file__),native.OUT/'plan.json',native.OUT/'result.json',OUT/'requests.json',OUT/'lanl-tops-form.html']
    paths += [native.OUT/(row['name']+'.npz') for row in original[2:]]
    times=[json.loads((native.OUT/(r['name']+'-retrieval.json')).read_text())['seconds'] for r in original]
    ex.write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='0f6b09ed',
        claim='Supply a positive photon total-loss response surrogate with explicit opacity-distribution uncertainty, without treating mean opacities as a uniquely identified kernel.',
        reassessment='Even 999 logarithmic group Planck coefficients used as constant absorptions shift the Rosseland flux coefficient by 14 to 19 percent. Moment bounds on the original broad log grid also remain wide. Retrieve native group moments only in the thermally relevant energy window and retain coarse tail data with distribution-independent tail bounds.',
        window_edges_keV=windows.tolist(),coarse_tail_split_indices=[435,847],groups_per_window=999,
        method='One harmonic-loss pole per group. Bound the unknown distribution via native arithmetic/harmonic moments, monotone Bprime/B across the bin, printed opacity intervals and the minimum of the tabulated log-PCHIP scattering curve over the whole queried window. No mean-fitted spectral amplitudes.',
        gates=dict(native_gray_recombination=.001,uniform_complex_resolvent_relative_to_DC=.01,uniform_energy_response_absolute=.01),
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths},
        budget=dict(queries=32,forecast_network_seconds=32*max(times),hard_network_seconds=180,per_query_hard_seconds=20,CPU_workers=1,new_native_EOS_calls=0,new_stellar_steps=0,automatic_expansion=False),
        limits='Uniform frequency bounds concern a local total-loss resolvent under declared positive-opacity/scattering assumptions. No angular/Compton gain kernel, general opacity model certificate, actual-temperature interpolation or whole-star thermal closure.'))
    ex.write(OUT/'symbolic.json',symbolic());print('NETWORK FORECAST',32*max(times),flush=True)


def fetch():
    plan=json.loads((OUT/'plan.json').read_text())
    for p,sha in plan['bindings'].items():assert h.digest(h.ROOT/p)==sha,p
    start=time.monotonic();results=[]
    env=dict(vars(native.retrieval),OUT=OUT);env['save']=lambda name,value:ex.write(OUT/name,value)
    function=FunctionType(native.retrieval.fetch.__code__,env,argdefs=native.retrieval.fetch.__defaults__)
    for row in json.loads((OUT/'requests.json').read_text()):
        assert not list(OUT.glob(row['name']+'-*'));began=time.monotonic()
        try:result=function(row)
        except Exception as exc:
            ex.write(OUT/'failure.json',dict(error=repr(exc),name=row['name'],seconds=time.monotonic()-start));raise
        result['seconds']=time.monotonic()-began;results.append(result);ex.write(OUT/'progress.json',results)
        print('QUERY',len(results),row['name'],result['seconds'],flush=True)
        assert result['seconds']<20 and time.monotonic()-start<180,'frozen network budget'
    ex.write(OUT/'retrieval.json',dict(seconds=time.monotonic()-start,queries=len(results),passed=all(not r['warnings_present'] for r in results)))


def analyze():
    assert not (OUT/'result.json').exists();plan=json.loads((OUT/'plan.json').read_text());header=native.reader.namespace()['header']
    rows=[[],[]]
    for request in json.loads((OUT/'requests.json').read_text()):
        path=OUT/(request['name']+'-table.txt');common=header(path,dict(request,fields=dict(request['fields'],datype='gray')))
        assert common['input_passed'];lines=[s.strip() for s in path.read_text().splitlines() if s.strip()]
        j=next(j for j,s in enumerate(lines) if s.startswith('Energy') and 'density =' in s);count=request['boundaries']-1
        tokens=[s.split() for s in lines[j+1:j+1+count]];group=np.array(tokens,float);assert group.shape==(count,3)
        assert np.isfinite(group).all() and (group>0).all() and (group[:,1:]!=1e10).all()
        edges=np.geomspace(float(request['fields']['egplow']),float(request['fields']['egphigh']),count+1)
        assert max(native.audit.score(native.audit.F(float(x)),t[0]) for x,t in zip(edges[:-1],tokens))<=1
        second=json.loads((OUT/(request['name']+'-results-request.json')).read_text())
        assert float(second['egplow'])==edges[0] and float(second['egphigh'])==edges[-1]
        rows[request['temperature_index']].append((edges,group[:,1:]))
    results=[]
    for j,blocks in enumerate(rows):
        original=np.load(native.OUT/(['surf1000T0015.npz','surf1000T002.npz'][j]));T=float(original['T_keV']);coarse=original['bounds_keV'];split=plan['coarse_tail_split_indices']
        blocks=[(coarse[:split[0]+1],original['groups'][:split[0],1:])]+blocks+[(coarse[split[1]:],original['groups'][split[1]:,1:])]
        for a,b in zip(blocks[:-1],blocks[1:]):assert abs(a[0][-1]/b[0][0]-1)<1e-14
        edges=np.concatenate([blocks[0][0]]+[a[0][1:] for a in blocks[1:]]);groups=np.concatenate([b for _,b in blocks]);parts=[]
        for n,(edge,group) in enumerate(blocks):parts.append(bounds(edge,group,T,original['spectrum'],tail=n in [0,len(blocks)-1]))
        w=np.concatenate([p['w'] for p in parts]);H=np.concatenate([p['H'] for p in parts]);inverse=np.concatenate([p['inverse'] for p in parts]);energy=np.concatenate([p['energy'] for p in parts])
        dc=w[:,1]@H;means=np.array([w[:,0]@groups[:,1]/w[:,0].sum(),w[:,1].sum()/dc]);gray=original['gray'][[2,1]]
        comparison=abs(means/gray-1);inverse_bound=float(w[:,1]@inverse/dc);energy_bound=float(w[:,1]@energy)
        passed=bool(comparison.max()<.001 and inverse_bound<.01 and energy_bound<.01)
        np.savez_compressed(OUT/f'loss-{j}.npz',bounds_keV=edges,groups=groups,weights=w,
            loss_rate_per_opacity=groups[:,0],resolvent_uncertainty=inverse,energy_response_uncertainty=energy,T_keV=T,
            rho_neutral_g_cm3=original['gray'][0],passed=passed)
        results.append(dict(T_keV=T,groups=len(groups),passed=passed,gray_relative=comparison.tolist(),
            uniform_complex_resolvent_relative_to_DC=inverse_bound,uniform_energy_response_absolute=energy_bound))
    result=dict(classification='Counterexample candidate',passed=all(r['passed'] for r in results),checks=results,
        local_loss_propagator_constructed=True,angular_gain_kernel_complete=False,material_energy_exchange_closed=False,
        whole_star_heat_closed=False,full_dynamic_charge_solved=False)
    ex.write(OUT/'result.json',result);print('RESULT',result,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','fetch','analyze'])
    globals()[parser.parse_args().action]()
