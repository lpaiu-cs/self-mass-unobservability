"""Independent identities, native differences and retained-rate envelope."""
import json
import re
import signal
import time
import numpy as np
import sympy as s
import def_native_ion_closure as task
from def_native_ion_history import provider


def main():
    out=task.OUT;assert not (out/'audit.json').exists();begin=time.monotonic();signal.alarm(45)
    task.write(out/'audit-plan.json',dict(classification='Counterexample candidate',seconds=45,EOS_calls=180,
        controls='Independent two-step fixed-population thermodynamic differences at three accepted adiabat states. Analytic native first-law and retained RR comparison principle. Independent published hydrogen fit and interval event ceiling.',
        gates=dict(first_law=1e-4,step_contrast=1e-4,rate_implementation=1e-6),
        coding_repair='Preserve the first audit TypeError from abs(list), before the first control was written. Use numpy.abs. Reserve its worst-case129 state evaluations and9.1s from the unused interrupted adiabat allocation; this single audit restart does not enlarge the combined allocated call/time budget or change gates.',
        boundary='The registered49node adiabat did not complete. Auditing its accepted prefix does not change that failure or certify omitted reaction channels.',source_sha256=task.sha(__file__)))
    rho,T=s.symbols('rho T',positive=True);F=s.Function('F')(rho,T)
    P=rho**2*s.diff(F,rho);S=-s.diff(F,T);U=F+T*S
    assert s.simplify(s.diff(U,T)-T*s.diff(S,T))==0
    assert s.simplify(s.diff(U,rho)-T*s.diff(S,rho)-P/rho**2)==0
    t=s.symbols('t',real=True);rate=s.Function('rate')(t);q=s.Function('q')(t);A=s.Function('A')(t)
    assert s.simplify(s.diff(s.exp(A)*q,t).subs({s.diff(A,t):rate,s.diff(q,t):-rate*q}))==0
    symbolic=dict(classification='Proven',passed=True,
        fixed_composition='At fixed chemical coordinates u=F-T F_T, s=-F_T and P=rho^2 F_rho imply u_T=T s_T and u_rho=T s_rho+P/rho^2.',
        comparison='If dx/dtau >= -R(tau)x, x>=0, and R<=Rmax, then x(tau)>=x0 exp(-integral Rmax dtau). Nonnegative ionization cannot weaken this lower bound. Missing neutralization channels cannot be omitted from a claim about full physical R.',
        interval='For log-linear rho,T on each saved interval, decreasing alpha_RR(T), ne<=rho*ne_per_rho and dproper_time/dt<=1, integral R dtau <= sum dt*max(rho endpoints)*alpha_RR(min(T endpoints))*ne_per_rho.')
    task.write(out/'symbolic.json',symbolic)
    cache=np.load(out/'second-adiabat-native-states.npz');initial=cache['number_fractions'][0];lr0=cache['lrho'][0]
    ion=task.Ions(cap=180);controls=[]
    for x in [0.,-2.,-4.]:
        j=np.flatnonzero(abs(cache['lrho']-(lr0+x))<1e-12)[-1]
        lr=float(cache['lrho'][j]);lt=float(cache['logT'][j]);fields=cache['fields'][j]
        base,_,_=ion.constrain(lr,lt,initial,fields,target_molecules=cache['molecular_H_fractions'][0],tolerance=1e-12)
        a=base['eos'];derivatives=[];laws=[]
        for h in [2e-4,1e-4]:
            neighbors=[]
            for dr,dt in [(h,0),(-h,0),(0,h),(0,-h)]:
                b,_,_=ion.constrain(lr+dr,lt+dt,initial,fields,target_molecules=cache['molecular_H_fractions'][0],tolerance=1e-12)
                neighbors.append(b['eos'])
            ar,br,at,bt=neighbors;radial=(ar-br)/(2*h);thermal=(at-bt)/(2*h)
            laws.append([float((np.exp(lt)*thermal[3]-thermal[2])/thermal[2]),
                         float((np.exp(lt)*radial[3]-radial[2]+a[1]/a[0])/(a[1]/a[0]))])
            gamma=radial[1]/a[1]+thermal[1]/a[1]*(a[1]/a[0]-radial[2])/thermal[2]
            derivatives.append([float(thermal[2]),float(gamma)])
        step=float(max(abs(np.array(derivatives[0])/derivatives[1]-1)))
        controls.append(dict(log_density_ratio=x,T=float(np.exp(lt)),first_law_residuals=laws,
            cvT=derivatives[-1][0],gamma1=derivatives[-1][1],step_relative=step,
            passed=bool(np.max(np.abs(laws))<1e-4 and step<1e-4 and derivatives[-1][0]>0 and derivatives[-1][1]>1)))
    # Independent evaluation of the four published H fit parameters.
    source=(out/'rrfit.f').read_text();m=re.search(r'data\(rnew\(i,\s*1,\s*1\),i=1,4\)/([^/]+)',source,re.I)
    assert m;pars=np.array([float(v) for v in m.group(1).replace('&','').split(',')]);a,b,t0,t1=pars
    temps=np.array([800.,1e4,1e5]);root=np.sqrt(temps/t0)
    expected=a/(root*(1+root)**(1-b)*(1+np.sqrt(temps/t1))**(1+b));rr=provider()
    rate_error=float(max(abs(rr(temps)/expected-1)));assert rate_error<1e-6
    rows=[];historical=json.loads((out/'history.json').read_text())
    for n in [1792,896]:
        d=np.load(out/f'material-history-{n}.npz');temp=d['T'];ne=d['electron_ceiling'];dt=np.diff(d['t'])[:,None]
        envelope=np.sum(dt*np.maximum(ne[:-1],ne[1:])*rr(np.minimum(temp[:-1],temp[1:])),axis=0)
        path=next(p for p in historical['paths'] if p['cells']==n)
        detail=[]
        for j,old in enumerate(path['rows']):
            lower=old['initial_H_ion']*np.exp(-envelope[j]);ratio=old['required_log_neutralization']/envelope[j]
            assert lower>old['LTE_final_H_ion']
            detail.append(dict(quantile=old['outside_baryon_quantile'],retained_RR_event_ceiling=float(envelope[j]),
                retained_RR_lower_H_ion=float(lower),required_over_ceiling=float(ratio)))
        rows.append(dict(cells=n,rows=detail))
    task.write(out/'rate-envelope.json',dict(classification='Counterexample candidate',passed=True,paths=rows,
        source_model='Piecewise log-linear saved material histories; original monotone total spontaneous RR fit and maximum fully ionized electron inventory. This excludes no omitted neutralization process.',
        fit_parameters=pars.tolist(),independent_fit_relative=rate_error))
    # Recover the accepted prefix; never mark the original49-node target passed.
    count=json.loads((out/'second-adiabat-failure.json').read_text())['nodes'];density=np.linspace(0,-6,49)[:count]
    indices=[np.flatnonzero(abs(cache['lrho']-(lr0+x))<1e-12)[-1] for x in density]
    arr=cache['eos'][indices];temperatures=np.exp(cache['logT'][indices]);eq=np.load(task.native.OUT/'fine.npz');comparisons=[]
    for x in [-1.,-2.,-4.]:
        i=int(np.argmin(abs(density-x)));j=int(np.argmin(abs(eq['log_density_ratio']-x)))
        comparisons.append(dict(log_density_ratio=x,fixed_T=float(temperatures[i]),LTE_T=float(eq['T'][j]),
            fixed_to_LTE_pressure=float(arr[i,1]/eq['raw'][j,1]),fixed_H_ion=float(arr[i,14]),LTE_H_ion=float(eq['raw'][j,14])))
    np.savez_compressed(out/'accepted-adiabat-prefix.npz',log_density_ratio=density,T=temperatures,raw=arr,
        fields=cache['fields'][indices],number_fractions=cache['number_fractions'][indices])
    task.write(out/'adiabat-prefix.json',dict(classification='Counterexample candidate',accepted_nodes=count,
        planned_nodes=49,planned_domain_completed=False,minimum_T=float(min(temperatures)),last_log_density_ratio=float(density[-1]),
        comparisons=comparisons,physical_flow_or_tail_completed=False))
    result=dict(classification='Counterexample candidate',passed=all(c['passed'] for c in controls),controls=controls,
        EOS_calls=ion.calls,seconds=time.monotonic()-begin,source_sha256=task.sha(__file__),
        fixed_population_first_law_checked=True,retained_RR_incompatible_with_LTE_path=True,
        full_adiabat_passed=False,physical_chemistry_closed=False,final_charge_solved=False,full_goal_complete=False)
    task.write(out/'audit.json',result);ion.save('audit-native-states.npz');signal.alarm(0)
    print(json.dumps(result),flush=True);assert result['passed']


if __name__=='__main__':main()
