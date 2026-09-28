"""Proven conditional boundary: hard removal of an ideal-mixing species."""
import json, sys
import numpy as np
import sympy as sp
import direct_eos_gr as g

OUT=g.OUT/'molecular-switch-identities'


def run():
    assert not OUT.exists();OUT.mkdir()
    n,x,kT,C=sp.symbols('n x kT C',positive=True)
    f=lambda y:y*(sp.log(y)-1)
    change=kT*(f(n-2*x)+f(x)-f(n))+C*x
    assert sp.limit(change/x,x,0,dir='+')==-sp.oo
    assert sp.limit(change,x,0,dir='+')==0
    assert sp.simplify(sp.diff(change,x,2)-kT*(1/x+4/(n-2*x)))==0
    result=dict(classification='Proven',passed=True,
        theorem='Let a removed molecular coordinate x enter the same free energy through kT*x*(ln x-1), kT>0. Suppose at a restricted equilibrium there is a feasible creation direction, every depleted donor has positive population, and every other term has a finite directional derivative there. Then [F(x)-F(0)]/x tends to minus infinity as x decreases to zero. Thus x=0 is not a local minimum of the full free energy, and the attained full minimum is strictly below the restricted minimum.',
        continuity_boundary='If the full and restricted minima are continuous in T near T0 and the hard switch uses the full minimum below T0 and restricted minimum at/above T0, the minimized free energy has a strictly positive jump at T0. Uniform classical C1/C2 derivative certification across that switch is impossible without changing the model or restricting the domain.',
        assumptions='Positive donor populations, a feasible nuclear/charge-conserving direction, finite nonideal/partition directional derivatives, attained minima and continuity of both minima. These are conditional assumptions, not silently certified for every native branch.',
        symbolic_control='H2 creation from two positive H atoms: F(x)-F(0)=kT[f(n-2x)+f(x)-f(n)]+C*x, f(y)=y(ln y-1); limit of the quotient is -infinity and curvature is kT[1/x+4/(n-2x)].',
        physical_EOS_certified=False)
    (OUT/'identity.json').write_text(json.dumps(result,indent=2)+'\n')
    data=dict(np.load(g.OUT/'molecular-switch/boundary.npz'));v=data['values'];t=np.log(1e6)+data['steps']
    # The same stored binary64 outputs evaluated with extended precision.
    # This reduces cancellation but is not an interval for native EOS error.
    ld=np.longdouble
    delta=(v[:,1,2].astype(ld)-v[:,3,2].astype(ld))-np.exp(t.astype(ld))*(v[:,1,3].astype(ld)-v[:,3,3].astype(ld))
    rows=[dict(cell=int(i),log_step=float(h),restricted_minus_retained_free_energy_erg_g=float(d))
        for i,h,d in zip(data['cells'],data['steps'],delta)]
    (OUT/'finite-free-energy.json').write_text(json.dumps(dict(classification='Counterexample candidate',rows=rows,
        computed_from_saved_native_outputs=True,continuous_or_physical_error_bound=False),indent=2)+'\n')
    paths=[g.ROOT/'verification/molecular_switch_identity.py',g.OUT/'molecular-switch/boundary.npz',*OUT.glob('*.json')]
    (OUT/'manifest.json').write_text(json.dumps(dict(classification='Proven',
        sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths}),indent=2)+'\n')
    verify()


def verify():
    entries=json.loads((OUT/'manifest.json').read_text())['sha256']
    for rel,digest in entries.items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'identity.json').read_text())['passed']
    print('PASS MOLECULAR REMOVAL BOUNDARY',len(entries),'artifact SHA',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
