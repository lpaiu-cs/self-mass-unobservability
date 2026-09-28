"""Conditional variational bound for isotope translational mass terms.

Proven only for the declared model change. No derivative or plasma-error bound
is inferred from a free-energy bound, and no old EOS verdict is overwritten.
"""
import json, sys
import numpy as np
import sympy as sp
from mpmath import iv
import direct_eos_gr as g


def run():
    out=g.OUT/'isotope-mass-bound';out.mkdir(parents=True,exist_ok=True)
    assert not (out/'result.json').exists();iv.dps=60
    def save(name,value): (out/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')
    save('plan.json',dict(classification='Proven',checkpoint='932888a',
        comparison='Two variational free-energy models with finite infima on exactly the same feasible populations. Change only the monatomic charge-dependent isotope translational masses after aligning neutral-atom mass entropy gauges. Keep all electronic, molecular and nonideal terms fixed.',
        masses='Treat the saved nuclear/reference masses and the source-constant census electron mass as declared exact binary64 inputs. No experimental mass-uncertainty or electronic isotope-shift claim.',
        bound='With delta_iq=(3/2)log[(1-q me/M_i)/(1-q me/R_i)], a_i=min(0,delta_iZ), b_i=max(0,delta_iZ), the minimized free-energy change per baryon gram lies in [-NA kB T sum_i(X_i/A_i)b_i, -NA kB T sum_i(X_i/A_i)a_i].',
        scope='All admissible ion fractions for each specified T,X, including a non-monatomic pool unchanged between these two models. Not all isotope/molecular physics, not pressure/entropy/derivative errors and not full EOS certification.'))
    M,R,me,q=sp.symbols('M R me q',positive=True)
    delta=sp.Rational(3,2)*sp.log((1-q*me/M)/(1-q*me/R))
    expected=sp.Rational(3,2)*me*(M-R)/((M-q*me)*(R-q*me))
    assert sp.simplify(sp.diff(delta,q)-expected)==0
    assert sp.simplify(delta.subs(q,0))==0 and sp.simplify(delta.subs(R,M))==0
    # Exact one-nucleus/two-charge partition controls. The minimizer changes;
    # a mere composition gauge would give zero and fails the positive control.
    controls=[]
    for shift in ['-0.5','-0.0001','0','0.0001','0.5']:
        v=iv.mpf(shift);answer=-iv.ln((1+iv.exp(v))/2)
        lower=-max(v,iv.mpf(0));upper=-min(v,iv.mpf(0))
        assert answer.a>=lower.a and answer.b<=upper.b
        if shift=='0': assert answer.a==0 and answer.b==0
        else: assert answer.b<0 or answer.a>0
        controls.append(dict(shift=shift,free_energy_over_kT=str(answer)))
    data=dict(np.load(g.OLD/'initial-state.npz'));c=g.c;d=g.d
    census=json.loads((d.OUT/'isotope-no-go.json').read_text())
    electron=iv.mpf(census['finite_source_constant_electron_mass_amu'])
    weights=json.loads((d.OUT/'model-data.json').read_text())['atomic_weights']
    lower=[];upper=[];species=[]
    for k,name in enumerate(c.NAMES):
        i=int(np.flatnonzero(d.CHARGES==c.Z[k])[0]);mass=iv.mpf(float(c.W[k]));ref=iv.mpf(weights[i])
        charge=int(c.Z[k]);assert min(mass.a,ref.a)>(charge*electron).b
        end=iv.mpf('1.5')*iv.ln((1-charge*electron/mass)/(1-charge*electron/ref))
        a=iv.mpf(min(iv.mpf(0).a,end.a));b=iv.mpf(max(iv.mpf(0).b,end.b))
        lower.append(a);upper.append(b)
        species.append(dict(species=name,delta_min=str(a),delta_max=str(b)))
    rows=[];total_width=iv.mpf(0);maximum=iv.mpf(0);worst=None
    for i,(t,x,dm) in enumerate(zip(data['lnT'],data['X'],data['dm'])):
        prefactor=iv.mpf(c.NA)*iv.mpf(c.KB)*iv.exp(iv.mpf(float(t)))
        lo=-prefactor*sum(iv.mpf(float(v))/int(a)*b for v,a,b in zip(x,c.A,upper))
        hi=-prefactor*sum(iv.mpf(float(v))/int(a)*b for v,a,b in zip(x,c.A,lower))
        assert lo.b<=0 and hi.a>=0
        width=(hi-lo).b;total_width+=iv.mpf(float(dm))*width
        if width>maximum: maximum=width;worst=i
        rows.append(dict(cell=i,lower_erg_per_baryon_g=str(lo),upper_erg_per_baryon_g=str(hi)))
    save('result.json',dict(classification='Proven',passed=True,
        symbolic_monotonicity=True,zero_mass_difference_control=True,partition_positive_controls=controls,
        proof='For M,R>Z me, delta_iq has constant derivative sign and is between its q=0 and q=Z values. Every feasible normalized charge population therefore gives the stated linear free-energy interval. If L<=F_new(z)-F_old(z)<=U for all z in the same feasible set, taking infima yields L<=inf F_new-inf F_old<=U. This also holds with an unchanged reservoir assigned delta=0. No uniqueness or numerical minimization is needed.',
        species=species,rows=rows,worst_width_cell=worst,maximum_interval_width_erg_per_baryon_g=str(maximum),
        sum_cell_baryon_mass_times_width_erg=str(total_width),
        input_sha256={str(p.relative_to(g.ROOT)):c.sha(p) for p in [g.OLD/'initial-state.npz',
            d.OUT/'isotope-no-go.json',d.OUT/'model-data.json',out/'plan.json']},
        pressure_or_derivative_bound=False,full_physical_EOS_certified=False))
    print('ISOTOPE MASS BOUND',len(rows),'states; maximum width',str(maximum),'at',worst,flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
