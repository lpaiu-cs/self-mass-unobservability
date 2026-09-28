"""Correct only the stored pressure-mode EOS jet; preserve the measured heat drive."""
from types import ModuleType
import json,sys
import gr_subcell_heat_drive as original
import gr_eos_derivative_coordinates as coordinates

ROOT=original.ROOT;OUT=original.OUT.parent/'gr-subcell-heat-drive-coordinates';sha=original.sha


def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2)+'\n')


def source():
    before=(ROOT/'verification/gr_subcell_heat_drive.py').read_text()
    changes={
        'import opacity_tables as opacity':'import opacity_tables as opacity\nimport gr_eos_derivative_coordinates as coordinates',
        'slope=(P/rho-eos[:,9])/eos[:,10];gamma=eos[:,5]+eos[:,6]*slope':
        'cr,ct,ur,cv=coordinates.pressure_to_density(eos[:,7],eos[:,8],eos[:,9],eos[:,10])\n        slope=(P/rho-ur)/cv;gamma=cr+ct*slope'}
    after=before
    for a,b in changes.items():assert after.count(a)==1;after=after.replace(a,b)
    reverse=after
    for a,b in reversed(list(changes.items())):assert reverse.count(b)==1;reverse=reverse.replace(b,a)
    assert reverse==before;return after,changes


def prepare():
    assert not OUT.exists();original.verify_preflight();coordinates.verify();plan=original.bindings();OUT.mkdir()
    text,changes=source();(OUT/'candidate.py').write_text(text)
    files=[ROOT/'verification/gr_subcell_heat_drive_coordinates.py',OUT/'candidate.py',
        coordinates.OUT/'manifest.json',ROOT/'verification/gr_eos_derivative_coordinates.py',
        original.OUT/'preflight-manifest.json',original.OUT/'preflight.json',
        ROOT/'verification/opacity_cubic.py',ROOT/'verification/native_opacity.py']
    plan['bindings'].update({p.relative_to(ROOT).as_posix():sha(p) for p in files})
    plan.update(classification='Counterexample candidate',checkpoint='668bc739',substitutions=changes,
        correction='The stored subcell EOS calls use pressure mode1. Convert the four density/energy derivatives to density coordinates before using the isentropic-jet identity. The directly differentiated temperature/lapse profile, opacity, prescribed Q, both node sets, all inputs and every other formula are unchanged.',
        gate='All twelve preflight cell/node-count records must retain every non-jet field bit for bit. Verify the correction with the separately bound native two-coordinate controls; report the remaining spectral/jet and constitutive discrepancies without fitting Q or a tolerance.')
    save('plan.json',plan)
    for name in ['symbolic.json','controls.json']:(OUT/name).write_bytes((original.OUT/name).read_bytes())


def engine():
    text,_=source();assert text==(OUT/'candidate.py').read_text()
    name='gr_subcell_heat_drive_coordinates_candidate';obj=ModuleType(name);obj.__file__=str(OUT/'candidate.py')
    sys.modules[name]=obj;exec(compile(text,obj.__file__,'exec'),obj.__dict__);obj.OUT=OUT;return obj


def comparison():
    before=json.loads((original.OUT/'preflight.json').read_text());after=json.loads((OUT/'preflight.json').read_text())
    assert before['evaluations_completed'] and after['evaluations_completed'];rows=[]
    changed={'isentropic_jet_Fourier_moment_over_enthalpy','maximum_abs_jet_QF_over_enthalpy','maximum_spectral_jet_difference_over_enthalpy'}
    for a,b in zip(before['rows'],after['rows'],strict=True):
        assert a['cell']==b['cell'] and a['finite_8_16_moment_difference_over_enthalpy']==b['finite_8_16_moment_difference_over_enthalpy']
        for x,y in zip(a['rows'],b['rows'],strict=True):
            assert {k:v for k,v in x.items() if k not in changed}=={k:v for k,v in y.items() if k not in changed}
            target=y['spectral_Fourier_moment_over_enthalpy'];assert target!=0
            rows.append(dict(cell=b['cell'],nodes=y['nodes'],spectral_heat_drive_unchanged=True,
                old_jet_relative_moment_difference=abs(x['isentropic_jet_Fourier_moment_over_enthalpy']/target-1),
                corrected_jet_relative_moment_difference=abs(y['isentropic_jet_Fourier_moment_over_enthalpy']/target-1)))
    assert len(rows)==12
    save('comparison.json',dict(classification='Counterexample candidate',passed=True,rows=rows,
        all_nonjet_fields_identical=True,original_pressure_mode_misinterpretation_preserved=True,full_GR_evolution=False))
    files=[OUT/'preflight-manifest.json',OUT/'comparison.json',coordinates.OUT/'manifest.json']
    save('repair-manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in files}));verify_preflight()


def preflight():engine().preflight();comparison()


def verify_preflight():
    engine().verify_preflight();coordinates.verify()
    for rel,digest in json.loads((OUT/'repair-manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    r=json.loads((OUT/'comparison.json').read_text());assert r['passed'] and r['all_nonjet_fields_identical']
    print('PASS pressure-coordinate jet correction with every directly measured heat-drive field unchanged',flush=True)


def run():verify_preflight();engine().run()


def verify():verify_preflight();engine().verify()


if __name__=='__main__':globals()[sys.argv[1]]()
