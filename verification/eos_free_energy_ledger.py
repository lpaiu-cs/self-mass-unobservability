"""Source-bound bookkeeping for the declared fixed-temperature EOS Hessian."""
import json, sys
from pathlib import Path
import sympy as sp
import excitation_scaled as scaled

g=scaled.g; OUT=g.OUT/'free-energy-ledger'
SOURCE=g.d.CACHE/'full-integral-source/src'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def run():
    assert not OUT.exists(); OUT.mkdir()
    n=sp.symbols('n0:3',positive=True); T,k=sp.symbols('T k',positive=True)
    coefficients=[sp.Function('a'+str(i))(T) for i in range(3)]
    lnw=[sp.Function('lnw'+str(i))(T) for i in range(3)]
    ideal=k*T*sum(v*(sp.log(v)-1) for v in n)
    linear=sum(a*v for a,v in zip(coefficients,n))
    PL=-k*T*sum(v*w for v,w in zip(n,lnw))
    assert sp.hessian(ideal,n)==sp.diag(*[k*T/v for v in n])
    assert sp.hessian(linear+PL,n)==sp.zeros(3,3)
    assert sp.diff(PL,n[0],T)==-k*(T*sp.diff(lnw[0],T)+lnw[0])
    cr,s0,u0,sum1,sum2=sp.symbols('cr s0 u0 sum1 sum2')
    spi=s0+cr*sum1; upi=u0+cr*T*sum2; free_pl=cr*T*(sum2-sum1)
    free_pi=upi-T*spi-free_pl
    assert sp.expand(free_pi-(u0-T*s0))==0
    assert sp.expand(free_pi+free_pl-(u0-T*s0+free_pl))==0
    weights=[sp.symbols('w'+str(i)) for i in range(3)]
    slopes=[sp.symbols('d'+str(i)) for i in range(3)]
    s1=sum(v*(w+d) for v,w,d in zip(n,weights,slopes));s2=sum(v*d for v,d in zip(n,slopes))
    assert sp.expand(s2-s1+sum(v*w for v,w in zip(n,weights)))==0
    # Species-independent radiation terms have no internal-composition Hessian.
    radiation=sp.Function('F_rad')(T); assert sp.hessian(radiation,n)==sp.zeros(3,3)
    excerpts={}
    for name,ranges in {
        'free_eos_detailed.f90':[(3167,3229),(3496,3513)],
        'eos_calc.f90':[(390,403),(2434,2444)]}.items():
        path=SOURCE/name;lines=path.read_text().splitlines()
        excerpts[name]=dict(source=str(path),sha256=g.c.sha(path),
            excerpts=[dict(first=a,last=b,text='\n'.join(lines[a-1:b])) for a,b in ranges])
    detailed=(SOURCE/'free_eos_detailed.f90').read_text()
    for line in ['     spi = spi + cr*sumpl1','     upi = upi + cr*t*sumpl2',
        '  free_pi = upi -t*spi - free_pl','  free_excited = uexcited -t*sexcited',
        '       free_excited + free_pl + free_ion']:
        assert line in detailed,line
    save('source-ledger.json',dict(classification='Imported from prior work',sources=excerpts,
        terms=[
            dict(term='radiation',fixed_T_species_Hessian='zero',condition='Radiation depends on T; the species redistribution holds T and volume fixed.'),
            dict(term='ideal atoms, ions and molecules',fixed_T_species_Hessian='kT*diag(1/n)',condition='Positive species support; fixed-T translational masses, statistical weights and internal binding/rotation/vibration terms are linear in species populations. Count molecular nuclear inventory with multiplicity two.'),
            dict(term='Planck-Larkin ground-state correction',fixed_T_species_Hessian='zero',condition='The temperature-only ln occupation weights are multiplied linearly by each population. The source adds PL to spi/upi and then subtracts it from the pressure-ionization bookkeeping column, leaving one PL term in the total.'),
            dict(term='electron plus nonlinear exchange',fixed_T_species_Hessian='charge pullback of the Legendre-compressibility coefficient',condition='Consistent interior branch and positive grand-canonical compressibility. Previously traced returned-input coefficients are frozen; their continuous evaluation is not certified.'),
            dict(term='Coulomb',fixed_T_species_Hessian='three-moment Hessian pullback',condition='Declared DH/OCP model, fixed T, frozen coefficients.'),
            dict(term='MDH pressure ionization',fixed_T_species_Hessian='nine-moment polynomial Hessian pullback',condition='Declared radii and quad=10. Ground-state MDH trace precedes PL additions.'),
            dict(term='density-dependent excitation',fixed_T_species_Hessian='-(D.T*C+C.T*D+C.T*W*C)',condition='Differentiable fixed branch and linear moments; actual scaled component derivatives are separately captured.'),
        ],
        scope='A source-bound algebraic ledger. It is not a measured total-free-energy or equilibrium-residual replay and is not a continuous/physical derivative certificate.'))
    save('symbolic.json',dict(classification='Proven',passed=True,
        identities=['Hessian_n of ideal translational term is kT*diag(1/n).',
            'Fixed-temperature linear binding/statistical/PL terms have zero species Hessian.',
            'Source PL bookkeeping leaves exactly one PL term; the traced ground-state MDH term excludes it.',
            'Zero species Hessian does not imply zero temperature gradient or mixed derivative.'],
        conditional_coverage='For the declared fixed-temperature, fixed-volume, positive-support, differentiable finite model with exactly the listed free-energy terms, all nonzero species-Hessian contributions are the ideal, electron/exchange, Coulomb, MDH and excitation terms.',
        unresolved=['Native derivative accuracy for partition approximations and cutoffs','Omitted-species and branch boundaries','Total gradient and coupled stationarity residual','Uniform curvature and derivative bounds over the relevant continuous domain','Physical model discrepancy and calibration','Self-consistent GR evolution and nonlinear observation closure']))
    paths=[OUT/'source-ledger.json',OUT/'symbolic.json',g.ROOT/'verification/eos_free_energy_ledger.py']
    save('manifest.json',dict(classification='Proven',sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        source_sha256={v['source']:v['sha256'] for v in excerpts.values()}))
    print('PASS EOS FREE ENERGY LEDGER: symbolic bookkeeping and source anchors',flush=True)


def verify():
    m=json.loads((OUT/'manifest.json').read_text())
    for rel,digest in m['sha256'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    for path,digest in m['source_sha256'].items(): assert g.c.sha(Path(path))==digest,path
    print('PASS EOS FREE ENERGY LEDGER',len(m['sha256']),'artifact SHA',flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
