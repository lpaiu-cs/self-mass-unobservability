"""Preserve and localize the Request33 full-output history-control failure."""
import json, shutil, subprocess, sys
import numpy as np
import direct_eos_gr as g

NAMES=['rho','P','u','s','gamma1','chi_rho','chi_T','rho_P','rho_T','u_P','u_T',
    'H2_reference','eta','rmue','fh2','fhe2','fhe3','xmu1','xmu3','lambda','gamma_e','sound2']


def probe():
    data=dict(np.load(g.OUT/'reference-state.npz'));records=[];raw=[]
    for evaluator in ['ct','tight','full']:
        eos=g.d.EOS(evaluator)
        for i,j in [(0,89),(2972,5734),(5734,0)]:
            r,t,x=data['lnd'][i],data['lnT'][i],data['X'][i]
            a=eos(2,r,t,x);b=eos(2,r,t,x)
            eos(1,np.log(a[1]),t,x);z=eos(2,r,t,x)
            eos(2,data['lnd'][j],data['lnT'][j],data['X'][j]);w=eos(2,r,t,x)
            values=np.array([a,b,z,w]);raw.append(values)
            differences={}
            for k,name in enumerate(NAMES):
                if np.any(values[:,k]!=values[0,k]):
                    differences[name]=dict(values=values[:,k].tolist(),
                        maximum_absolute_change=float(abs(values[:,k]-values[0,k]).max()),
                        maximum_relative_change=float(abs(values[:,k]-values[0,k]).max()/max(abs(values[0,k]),1.)))
            record=dict(evaluator=evaluator,cell=i,remote_cell=j,differences=differences,
                ordered_history=['initial','repeat','after pressure inverse','after remote cell'])
            records.append(record);print('EOS HISTORY',record,flush=True)
    np.savez_compressed(g.OUT/'history-probe.npz',values=np.array(raw))
    g.save('history-probe.json',dict(classification='Counterexample candidate',rows=records,
        original_inverse_failure_preserved=True,physical_or_continuous_certificate=False))


def repair_contract():
    assert not (g.OUT/'output-contract-repair.json').exists()
    original=json.loads((g.OUT/'plan.json').read_text())
    before=subprocess.check_output(['git','show','4b861ce:verification/direct_eos_gr.py'],cwd=g.ROOT)
    assert g.e.hashlib.sha256(before).hexdigest()==original['inputs_sha256']['verification/direct_eos_gr.py']
    (g.OUT/'before-output-contract-direct_eos_gr.py').write_bytes(before)
    for name in ['plan.json','inverse-control.json']:
        shutil.copy2(g.OUT/name,g.OUT/('before-output-contract-'+name))
    source=g.d.CACHE/'full-integral-source/src/coulomb.f90';text=source.read_text()
    assert text.count('gamma_e = dc_lambda/lambda/gamma_e_const')==1
    assert text.index("error stop 'coulomb: the diffraction correction is disabled.'")<text.index('gamma_e = dc_lambda/lambda/gamma_e_const')
    assert 'if(if_dc.eq.1) then' in text
    target=g.OUT/'sources/coulomb.f90';target.parent.mkdir(exist_ok=True);shutil.copy2(source,target)
    rows=json.loads((g.OUT/'history-probe.json').read_text())['rows']
    assert all(set(r['differences'])=={'gamma_e'} for r in rows)
    original['inputs_sha256']['verification/direct_eos_gr.py']=g.c.sha(g.ROOT/'verification/direct_eos_gr.py')
    original['output_contract']=dict(legacy_bridge_slots=22,defined_native_slots=list(range(20))+[21],
        returned_slots=21,removed_output='gamma_e: unassigned by the disabled diffraction branch',
        effect='Drop the undefined diagnostic only. Keep all thermodynamic, electron and sound-speed results and the original inverse/history criteria.',
        sound_squared_return_index=20)
    g.save('plan.json',original)
    g.save('output-contract-repair.json',dict(classification='Counterexample candidate',
        preserved_original_failure=True,source_sha256=g.c.sha(target),
        before_plan_sha256=g.c.sha(g.OUT/'before-output-contract-plan.json'),after_plan_sha256=g.c.sha(g.OUT/'plan.json'),
        undefined_output='gamma_e',physical_EOS_or_tolerance_changed=False,
        previous_pressure_energy_entropy_electron_derivative_tests_used_this_output=False))
    print('REPAIRED OUTPUT CONTRACT: 21 defined results; undefined gamma_e omitted',flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
