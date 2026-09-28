"""Connect EOS-resolved H/He II weights to reversible optical line factors.

Occupation ratios require an explicit optical survival model. PL thermodynamic
weights alone do not certify that interpretation; unsupported metal levels fail
the state-set contract instead of being silently populated in opacity only.
"""
import json
import subprocess
from pathlib import Path
import numpy as np
import sympy as sp
import def_photon_eos_populations as p
import def_photon_eos_levels as levels


def run():
    target=p.OUT/'level-balance.json';assert not target.exists()
    s=np.load(p.OUT/'eos-state.npz');q=np.load(levels.OUT/'level-sums.npz')
    inventory=np.load(p.a.previous.OUT/'inventory.npz');T=float(q['T']);n=np.arange(1,11.)
    source=levels.OUT/'constants.f90'
    source.write_text('''program native_constants
  use mod_free_eos_constants, only: c2,rydberg,electron_mass,h_mass
  implicit none
  write(*,'(4es26.17)') c2,rydberg,electron_mass,h_mass
end program
''')
    parent=p.SOURCE.parent.parent/'build/src';binary=levels.CACHE/'constants'
    subprocess.run(['gfortran','-I'+str(parent),str(source),'-o',str(binary)],check=True,capture_output=True,timeout=30)
    constants=np.fromstring(subprocess.check_output([str(binary)]).decode(),sep=' ')
    assert len(constants)==4
    np.savetxt(levels.OUT/'constants.txt',constants,fmt='%.17e')
    c2,rydberg,electron,h_mass=constants
    # Use the actual thermodynamic excitation correction as the common
    # ground-relative normalization, never an empirical line-strength factor.
    output={};rows=[]
    for z,tag,element,charge in [(1,(1,0),0,0),(2,(2,1),1,1)]:
        j=list(map(tuple,s['ids'][0])).index(tag)
        L=float(s['value'][0,j,3]/s['scale'][0,j])
        tails=q['cumulative'][z-1,0];summands=tails-np.r_[tails[1:],0.]
        fractions=np.r_[1.,np.expm1(L)*summands[1:]/tails[1]]/np.exp(L)
        assert abs(fractions.sum()-1)<1e-14 and np.all(fractions>0)
        populations=fractions*inventory['ni'][element,charge]
        # Infer relative weights algebraically from the native partition.
        # Match qstar_calc's exact reduced-mass convention, including 4*m_H
        # for He II. This is not the separate nuclear-inventory mass law.
        mass=h_mass*(1 if z==1 else 4)
        binding=c2*rydberg*z*z/(1+electron/mass)/n**2
        excitation=binding[0]-binding;g=2*n*n
        relative_w=fractions/fractions[0]*g[0]/g*np.exp(excitation/T)
        residuals=[];naive=[];probabilities=[]
        for lower in range(9):
            for upper in range(lower+1,10):
                u=(excitation[upper]-excitation[lower])/T
                availability=relative_w[upper]/relative_w[lower]
                assert 0<availability<=1+1e-12
                # B_ul=1; g_l B_lu=g_u B_ul; corrected upward availability.
                absorption=populations[lower]*(g[upper]/g[lower])*availability-populations[upper]
                emission=populations[upper]
                residuals.append(abs(emission*np.expm1(u)/absorption-1))
                naive_abs=populations[lower]*(g[upper]/g[lower])-populations[upper]
                naive.append(abs(emission*np.expm1(u)/naive_abs-1))
                probabilities.append(availability)
        rows.append(dict(species=tag,excited_fraction=float(1-fractions[0]),
            population_sum_relative=float(populations.sum()/inventory['ni'][element,charge]-1),
            conditional_Kirchhoff_relative=max(residuals),unmodified_Einstein_relative=max(naive),
            transitions=len(residuals),availability_range=[min(probabilities),max(probabilities)]))
        output['fractions_Z'+str(z)]=fractions;output['populations_Z'+str(z)]=populations
        output['relative_weights_Z'+str(z)]=relative_w
    np.savez_compressed(p.OUT/'level-populations.npz',**output)
    x,A,wl,wu,gl,gu=sp.symbols('x A wl wu gl gu',positive=True)
    lower=A*gl*wl;upper=A*gu*wu*sp.exp(-x)
    chi=lower*gu/gl*wu/wl-upper;j=upper
    assert sp.simplify(j/chi-1/(sp.exp(x)-1))==0
    # Missing upper state at finite T: positive upward absorption and zero
    # spontaneous emission cannot satisfy Kirchhoff. No numerical tolerance
    # can supply the missing state or its free-energy contribution.
    assert sp.simplify(chi.subs(wu,0))==0 and j.subs(wu,0)==0
    metal_species=[(int(Z),int(q)) for Z,ion in zip(p.a.previous.inventory_reader.g.d.CHARGES,inventory['ni'])
        if Z>2 for q,value in enumerate(ion[:int(Z)]) if value>0]
    assert all(0<=charge<Z for Z,charge in metal_species)
    p.a.write(target,dict(classification='Counterexample candidate',passed=all(r['conditional_Kirchhoff_relative']<1e-10 for r in rows),
        rows=rows,unsupported_excited_metal_species=metal_species,
        same_EOS_full_optical_model_ready=False,physical_optical_survival_certified=False,
        exact_statement='For n_i proportional to g_i w_i exp(-E_i/kT), an upward availability factor w_j/w_i and unmodified reverse coefficient give Kirchhoff exactly. At w_j=0 the same model removes BOTH bound-bound directions; keeping positive absorption with zero upper population violates LTE.',
        scientific_status='The algebra is Proven under the declared survival assumption. Applying thermodynamic PL weights as optical survival is Conjectural and is not inferred from matching sums. The current metal EOS has no excited-state partition correction, so an atomic level free-energy extension is required before claiming common-state full opacity.',
        sources=['https://academic.oup.com/mnras/article/335/2/499/1046898','https://adsabs.harvard.edu/pdf/1988ApJ...331..794H'],
        boundary='Provider work and a conditional optical contract, not a calibrated opacity or a new coupled/stellar evolution. Exact native hydrogenic energy constants are used; cross sections and optical-survival interpretation still require physical closure.'))
    print('LEVEL BALANCE',rows,'unsupported metal species',len(metal_species))


if __name__=='__main__':run()
