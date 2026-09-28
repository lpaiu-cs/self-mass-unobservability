"""Provider return-contract audit and source-table PIMC calibration replay."""
import json, sys
import numpy as np
import mpmath as mp
import gr_dense_plasma as d

g=d.g;OUT=g.OUT/'gr-dense-plasma-audit'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();d.bindings()
    # Table 2 is already used in the published fit; it is not a held-out test.
    table=[dict(R=500,gamma=v,value=u) for v,u in zip(
        [90,101,112,123,134,145,156,167,178],[6.344,6.951,7.565,8.185,8.811,9.438,10.072,10.707,11.343])]
    save('source-table.json',dict(classification='Imported from prior work',
        source='https://arxiv.org/abs/2112.04822',table=2,pdf_page=9,
        quantity='(U-U_M)/T per ion, includes classical kinetic energy 3/2; theta=Gamma*sqrt(3/R).',
        Madelung='U_M/T=-0.895929256*Gamma; bcc constant from frozen EOS22 FHARM12.',
        entries=table,decimal_last_place=.001,statistical_error_bars_provided_in_table=False,
        independent_holdout=False))
    paths=[g.ROOT/'verification/gr_dense_plasma_audit.py',d.OUT/'plan.json',d.OUT/'runtime.json',
        d.OUT/'quantum-controls.json',d.OUT/'symbolic.json',OUT/'source-table.json']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='e9500e5',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        identity_score_tolerance=1e-10,mp_classical_score_tolerance=1e-11,
        calibration_displayed_thermal_relative_tolerance=.01,
        policy='Freeze before evaluating the PIMC table. The 1% gate is a declared numerical fit-replay diagnostic relative to displayed thermal energy, not an author uncertainty interval or physical certification. Preserve any failure. No coefficients are fitted here.',
        purpose='Audit original EOSFI22 excess/ideal groups and the missing quantum density derivative on all saved near-fully-ionized states; independently evaluate classical OCP energies at 60 decimal digits and compare the sourced quantum fit with its original nine PIMC calibration values.',
        scope='A finite provider consistency and calibration replay. No universal error envelope, many-body proof, independent holdout, or EOS/GR replacement.'))


def mp_classical(gamma):
    mp.mp.dps=60;G=mp.mpf(str(gamma));a1=mp.mpf('-.907347');a2=mp.mpf('.62849')
    a3=-mp.mpf('.8660254038')-a1/mp.sqrt(a2)
    u=G**mp.mpf('1.5')*(a1/mp.sqrt(a2+G)+a3/(1+G))
    return u+mp.mpf('.004500')*G*G/(170+G)-mp.mpf('.000084')*G*G/(mp.mpf('.0037')+G*G)


def calibration():
    plan=bindings();p=d.Provider();table=json.loads((OUT/'source-table.json').read_text());rows=[]
    for r in table['entries']:
        R=r['R'];gamma=r['gamma'];theta=gamma*np.sqrt(3/R)
        cl=p.call('fition9',[float(gamma)],6)[1];ind=float(mp_classical(gamma))
        cq=d.quantum(R,theta)[1];raw=p.call('liqubc',[float(R),float(theta)],6)[1]
        predicted=1.5+cl+cq+.895929256*gamma;delta=predicted-r['value'];score=abs(delta)/r['value']
        row=dict(**r,theta=float(theta),classical_provider=float(cl),classical_independent=ind,
            classical_score=float(abs(cl-ind)/max(1,abs(ind))),quantum_derived=float(cq),quantum_provider=float(raw),
            predicted=float(predicted),difference=float(delta),relative_displayed_thermal_difference=float(score),
            displayed_rounding_half_width=.0005,passed=bool(score<plan['calibration_displayed_thermal_relative_tolerance']))
        rows.append(row)
    assert all(r['classical_score']<plan['mp_classical_score_tolerance'] for r in rows)
    save('calibration.json',dict(classification='Counterexample candidate',rows=rows,
        calibration_replay_passed=all(r['passed'] for r in rows),
        maximum_relative_displayed_thermal_difference=max(r['relative_displayed_thermal_difference'] for r in rows),
        maximum_absolute_difference=max(abs(r['difference']) for r in rows),
        no_refitting=True,independent_holdout=False,physical_error_certificate=False))
    print('PIMC calibration replay',max(r['relative_displayed_thermal_difference'] for r in rows),flush=True)


def bindings():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    return plan


def run():
    plan=bindings();d.verify();a=np.load(d.OUT/'stellar-comparison.npz');state,native=d.state_data();p=d.Provider();rows=[]
    for k,i in enumerate(a['cells']):
        X=state['X'][i];count=X/g.c.A;x=count/count.sum();z=g.c.Z
        rs,ge=a['parameters'][k,:2];ideal=np.zeros(7);rawql=np.zeros(7);domain=[]
        for j in np.flatnonzero(x>0):
            gamma=ge*z[j]**(5/3);theta=gamma/np.sqrt(rs)*np.sqrt(3/(1822.88848*g.c.A[j]))/z[j]**(7/6)
            R=3*(gamma/theta)**2;fid=1.5*np.log(theta*theta/gamma)-1.323515
            ideal+=x[j]*np.array([fid,1.5,1.,1.5-fid,1.5,1.,1.])
            rawql+=x[j]*d.seven(p.call('liqubc',[float(R),float(theta)],6))
            domain.append([float(R),float(theta),float(gamma)])
        groupdifference=a['raw_groups'][k,1]-a['raw_groups'][k,0]-ideal
        expected=rawql.copy();expected[6]=rawql[5]
        error=float(np.max(abs(groupdifference-expected)/np.maximum(1,abs(expected))))
        rawtotal=a['raw_groups'][k,1]+a['mixing'][k];correcttotal=a['corrected'][k]+ideal
        dropped=(x>0)&(x<1e-7)
        rows.append(dict(cell=int(i),return_group_identity_score=error,
            quantum_F_omitted_from_first_group=float(a['quantum'][k,0]),
            raw_quantum_PDR_error=float(rawql[6]-a['quantum'][k,6]),
            PDR2_uses_PDTQL_error=float(rawql[5]-rawql[6]),
            total_ion_nonideal_PDR_difference=float(rawtotal[6]-correcttotal[6]),
            total_ion_nonideal_PDR_relative_difference=float((rawtotal[6]-correcttotal[6])/max(1e-10,abs(correcttotal[6]))),
            trace_isotopes_below_MELANGE_cutoff=int(dropped.sum()),
            trace_number_fraction=float(x[dropped].sum()),trace_baryon_fraction=float(X[dropped].sum()),
            R_range=[float(np.min(domain,axis=0)[0]),float(np.max(domain,axis=0)[0])],
            theta_range=[float(np.min(domain,axis=0)[1]),float(np.max(domain,axis=0)[1])],
            gamma_range=[float(np.min(domain,axis=0)[2]),float(np.max(domain,axis=0)[2])]))
    passed=all(r['return_group_identity_score']<plan['identity_score_tolerance'] for r in rows)
    save('return-audit.json',dict(classification='Counterexample candidate',cells=len(rows),source_identity_passed=passed,rows=rows,
        maximum_identity_score=max(r['return_group_identity_score'] for r in rows),
        max_quantum_F_omitted=max(r['quantum_F_omitted_from_first_group'] for r in rows),
        max_total_ion_nonideal_PDR_relative_difference=max(abs(r['total_ion_nonideal_PDR_relative_difference']) for r in rows),
        diagnosis='The original excess group excludes liquid quantum corrections. Second minus first minus ideal equals LIQUBC except PDR uses PDTQL. Exact free-energy derivatives fix this separate mixed evaluator; source originals stay unchanged.'))
    assert passed
    save('result.json',dict(classification='Counterexample candidate',completed=True,source_contract_audit_passed=passed,
        source_result_sha256=g.c.sha(d.OUT/'manifest.json'),
        calibration_replay_passed=json.loads((OUT/'calibration.json').read_text())['calibration_replay_passed'],
        physical_EOS_certified=False,full_GR_evolution=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    bindings();d.verify()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    result=json.loads((OUT/'result.json').read_text());assert result['completed'] and result['source_contract_audit_passed']
    assert g.c.sha(d.OUT/'manifest.json')==result['source_result_sha256']
    print('PASS source return audit and immutable PIMC calibration diagnostics',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
