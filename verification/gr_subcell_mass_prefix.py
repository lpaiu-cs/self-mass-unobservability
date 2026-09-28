"""Audit common radial masses from unchanged cell energies using exact dyadic sums."""
from fractions import Fraction as F
from pathlib import Path
import json,sys
import numpy as np
import sympy as sp
import gr_metric_coupled_subcell as metric
import gr_balanced_conservation as balance

ROOT=metric.ROOT;OUT=metric.OUT.parent/'gr-subcell-mass-prefix';sha=metric.sha;ld=np.longdouble


def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2)+'\n')


def exact(value):return F(*value.as_integer_ratio())


def prepare():
    assert not OUT.exists();prior=metric.bindings();balance.verify();OUT.mkdir()
    files=[Path(__file__),metric.OUT/'plan.json',balance.OUT/'manifest.json']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='900aecc6',
        bindings=dict(prior['bindings'],**{p.relative_to(ROOT).as_posix():sha(p) for p in files}),
        cells=5735,nodes=[8,16],
        method='Evaluate each unchanged reference shell energy from its stored mixed-precision native rho,u,C_X and positive coordinate weights in 64-significand-bit long double. Freeze each resulting geometric shell mass as its exact dyadic value. Sum these dyadics from centre to surface exactly; compare every common face with the original saved mass as an exact dyadic. Do not replace any saved profile or rescale energy.',
        boundary='Exact accumulation is exact only for the rounded shell masses. Native EOS, subcell quadrature, constants, continuum and physical errors remain. The new metric inverse recalls the original native library, whereas this audit uses the frozen mixed-precision reference outputs. These two evaluator paths are not identified.'))


def bindings():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert sha(ROOT/rel)==digest,rel
    return plan


def symbolic():
    k,m0=sp.symbols('k m0',positive=True);energy=sp.symbols('e0:3');target=sp.symbols('t0:3')
    actual=[m0+k*sum(energy[:i]) for i in range(4)]
    supplied=[m0+k*sum(target[:i]) for i in range(4)]
    for i in range(3):
        assert sp.expand(actual[i+1]-actual[i]-k*energy[i])==0
        assert sp.expand((actual[i]-supplied[i]).subs(dict(zip(energy,target))))==0
    return dict(classification='Proven',passed=True,
        argument='For any finite ordered cells with prescribed target coordinate energies T_i and centre mass m_0, the radial constraint m_(i+1)-m_i=k*T_i uniquely supplies m_i=m_0+k*sum_(j<i)T_j by induction. If every local primitive solve reproduces its T_i and the other target moments using that supplied inner mass, all common mass faces agree. Conversely any global solution with those energies must use these same inner masses. Thus known conserved energies eliminate the apparent cross-cell mass-prefix dependence from local primitive inversion; local nodal metric feedback remains.',
        conditions='Fixed radius chart, shared positive conversion k=G/c^4, prescribed energies and centre mass. Does not establish local root existence/uniqueness, lapse reconstruction, Q evolution, a valid chart or continuous GR error.',
        imported_shared_flux='The previously proved common-face conservative energy flux and telescoping identity remain in gr-balanced-conservation; they are not counted as a new proof here.')


def calculate():
    plan=bindings();reference=json.loads((metric.OUT/'plan.json').read_text())
    grid=np.load(metric.g.OUT/'gr-increment-structure/path-4.npz')
    saved=[exact(v*ld(100)) for v in grid['mass_geom_m'].astype(ld)]
    assert len(saved)==5736 and saved[-1]==0 and saved[0]>0
    c=ld(metric.g.c.gr.C)*100;factor=ld(metric.g.c.gr.G)*1000/c**4;rows=[];prefixes=[]
    for number in plan['nodes']:
        masses={}
        for rel in sorted(set(reference['cell_sources'].values())):
            source=ROOT/rel;path=source.with_name(source.stem+f'-nodes-{number}.npz')
            data=metric.node_file(str(path))
            for j,index in enumerate(data['cells']):
                index=int(index);eos=data['eos'][j].astype(ld);rho=eos[:,0]
                rest=ld(data['C_X'][j])*c*c;eps=rho*(rest+eos[:,2])
                mass=factor*np.sum(eps*data['coordinate_weights_cm3'][j].astype(ld))
                assert np.isfinite(mass) and mass>0 and index not in masses;masses[index]=exact(mass)
        assert sorted(masses)==list(range(5735))
        prefix=[F(0)]*5736
        for i in range(5734,-1,-1):prefix[i]=prefix[i+1]+masses[i]
        jumps=[saved[i+1]+masses[i]-saved[i] for i in range(5735)]
        defects=[a-b for a,b in zip(prefix,saved)]
        assert all(prefix[i]-prefix[i+1]==masses[i] for i in range(5735))
        assert sum(jumps)==prefix[0]-saved[0]
        peak=max(range(5736),key=lambda i:abs(defects[i]))
        # A changed single cell must change every enclosing face and the total.
        altered=masses[0]/1024;assert sum(jumps)+altered!=prefix[0]-saved[0]
        rows.append(dict(nodes=number,shell_masses_geom_cm=list(map(str,(masses[i] for i in range(5735)))),
            prefix_faces_geom_cm=list(map(str,prefix)),maximum_face_defect_index=peak,
            maximum_face_defect_geom_cm=str(abs(defects[peak])),
            maximum_face_defect_over_total=str(abs(defects[peak])/saved[0]),
            total_mass_difference_geom_cm=str(defects[0]),total_mass_relative_difference=str(defects[0]/saved[0]),
            maximum_local_face_jump_over_total=str(max(map(abs,jumps))/saved[0]),
            sum_local_face_jumps_equals_total_defect=True,exact_prefix_identities_passed=True,
            altered_cell_negative_control_passed=True))
        prefixes.append(prefix)
    difference=max(abs(a-b) for a,b in zip(*prefixes))/saved[0]
    return dict(classification='Counterexample candidate',completed=True,cells=5735,rows=rows,
        reference_surface_mass_geom_cm=str(saved[0]),maximum_8_16_face_difference_over_total=str(difference),
        original_profile_changed=False,native_EOS_error_certified=False,continuous_errors_certified=False,full_GR_evolution=False)


def run():
    assert not (OUT/'result.json').exists();save('symbolic.json',symbolic());result=calculate();save('result.json',result)
    save('manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in OUT.iterdir() if p.is_file()}));verify()
    for r in result['rows']:print('PREFIX',r['nodes'],float(F(r['total_mass_relative_difference'])),float(F(r['maximum_face_defect_over_total'])),flush=True)


def verify():
    bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    result=json.loads((OUT/'result.json').read_text());assert result==calculate()
    assert json.loads((OUT/'symbolic.json').read_text())==symbolic()
    print('PASS exact dyadic mass-prefix accumulation and finite reference comparison at both node counts',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
