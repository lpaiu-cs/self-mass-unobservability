"""Independent exact arithmetic of saved complete GR shell/face inventories."""
import json,sys
from fractions import Fraction as F
import numpy as np
import gr_increment_structure as original

g=original.g;OUT=g.OUT/'gr-increment-structure-audit'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def fraction(x):return F(*x.as_integer_ratio())


def run():
    assert not OUT.exists();OUT.mkdir();original.verify()
    plan=json.loads((original.OUT/'plan.json').read_text());rows=[]
    reference=dict(np.load(g.OUT/'reference-state.npz'));previous=None
    for sub in plan['subdivisions']:
        data=dict(np.load(original.OUT/f'path-{sub}.npz'))
        record=json.loads((original.OUT/f'path-{sub}.json').read_text())
        for old,new in [('dm','original_baryon_g'),('X','original_X'),('s_B','original_reference_entropy')]:
            assert np.array_equal(reference[old],data[new]),(sub,old)
        assert np.all(np.diff(data['radius_m'])<0) and np.all(np.diff(data['mass_geom_m'])<0)
        shell=list(map(fraction,data['shell_mass_geom_m']));mass=list(map(fraction,data['mass_geom_m']))
        local=[mass[i]-mass[i+1]-shell[i] for i in range(len(shell))]
        total=sum(shell,F());gap=total-(mass[0]-mass[-1])
        assert sum(local,F())==-gap
        error=data['inner'][-1,1:]-data['outer'][-1,1:]
        assert float(abs(error).max())==record['interface_max']
        assert (float(abs(error).max())<plan['interface_tolerance'])==record['interface_passed']
        worst=max(range(len(shell)),key=lambda i:abs(local[i]/shell[i]))
        row=dict(classification='Proven',subdivision=sub,cells=len(shell),same_baryon_species_entropy=True,
            positive_ordered_faces=True,interface_passed=record['interface_passed'],
            total_shell_mass_geom_m=float(total),total_shell_minus_face_mass_exact=str(gap),
            total_relative_face_difference_exact=str(gap/(mass[0]-mass[-1])),
            total_relative_face_difference=float(gap/(mass[0]-mass[-1])),
            worst_local_face_subtraction_cell=worst,worst_local_face_subtraction_relative=float(local[worst]/shell[worst]),
            first_three_face_subtraction_relative=[float(local[i]/shell[i]) for i in range(3)],
            interpretation='Exact identities of the saved binary numbers, including branch matching/seed/rounding discrepancies. Independently stored shell increments are the cell masses; subtracting large face masses reintroduces cancellation.')
        if previous is not None:
            delta=max(abs(fraction(a)-fraction(b)) for a,b in zip(data['faces'].flat,previous.flat))
            row['finite_endpoint_maximum_exact']=str(delta)
            row['finite_endpoint_maximum']=float(delta)
            row['finite_endpoint_passed']=delta<F(plan['finite_endpoint_refinement_tolerance'])
            assert row['finite_endpoint_passed']==record['finite_endpoint_refinement_passed']
        previous=data['faces'];rows.append(row)
    result=dict(classification='Proven',completed=True,rows=rows,
        scope='Independent saved-value inventory, ordering, matching and exact summation audit. No direct physical EOS, continuum or actual GR time evolution certification.')
    save('result.json',result)
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in [
        g.ROOT/'verification/gr_increment_structure_audit.py',original.OUT/'manifest.json',
        g.OUT/'reference-state.npz',OUT/'result.json']}))
    verify();print('EXACT FULL STRUCTURE',rows,flush=True)


def verify():
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'result.json').read_text())['completed']
    print('PASS exact complete GR cell inventory audit',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
