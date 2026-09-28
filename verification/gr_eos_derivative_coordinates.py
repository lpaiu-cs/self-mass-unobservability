"""Keep pressure-mode and density-mode FreeEOS derivatives distinct."""
from fractions import Fraction as F
from pathlib import Path
import json,sys
import numpy as np
import sympy as sp
import gr_metric_coupled_subcell as metric
import gr_subcell_audit as old

ROOT=metric.ROOT;OUT=metric.OUT.parent/'gr-eos-derivative-coordinates';sha=metric.sha;ld=np.longdouble


def pressure_to_density(rho_P,rho_T,u_P,u_T):
    """Inputs use (lnP,lnT); return chi_rho,chi_T,u_lnrho,cv*T."""
    assert np.all(rho_P>0)
    return 1/rho_P,-rho_T/rho_P,u_P/rho_P,u_T-u_P*rho_T/rho_P


def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2)+'\n')


def symbolic():
    rp=sp.symbols('rp',positive=True)
    rt,up,ut,p,rho,dr,dt=sp.symbols('rt up ut P rho dr dt',nonzero=True)
    cr,ct,ur,cv=pressure_to_density(rp,rt,up,ut)
    assert sp.simplify(up*(dr-rt*dt)/rp+ut*dt-ur*dr-cv*dt)==0
    density_slope=(p/rho-ur)/cv;density_ad=density_slope/(cr+ct*density_slope)
    pressure_ad=(p/rho*rp-up)/(ut-p/rho*rt)
    assert sp.simplify(density_ad-pressure_ad)==0
    # Exercise the same numerical/scalar implementation in exact arithmetic.
    assert pressure_to_density(F(2),F(-3),F(5),F(7))==(F(1,2),F(3,2),F(5,2),F(29,2))
    return dict(classification='Proven',passed=True,
        coordinates='For rho_P=(partial ln rho/partial ln P)_T>0 and rho_T=(partial ln rho/partial ln T)_P, chi_rho=1/rho_P, chi_T=-rho_T/rho_P, u_lnrho=u_P/rho_P and cv*T=u_T-u_P*rho_T/rho_P. The transformation is a chain rule, not a physical EOS calibration.',
        entropy_slope='Under the differentiable Gibbs premise, nabla_ad=[(P/rho)*rho_P-u_P]/[u_T-(P/rho)*rho_T], equal to the density-coordinate expression. A pressure-mode u_T is generally neither cv*T nor cp*T.',
        applicability='The stored full/subcell reference nodes were obtained in pressure-inversion mode1. The metric-coupled inverse recalls mode2 and already uses density-mode derivatives. Archived fixed-metric inverse equations/Jacobians also recall mode2; only its fixed initial residual scale and the subsequent frozen-factor interpretation used the pressure-mode slot as cv*T.')


def prepare():
    assert not OUT.exists();prior=metric.bindings();old.verify();OUT.mkdir()
    files=[Path(__file__),metric.OUT/'plan.json',old.OUT/'manifest.json',old.inverse.OUT/'manifest.json',
        ROOT/'verification/gr_subcell_heat_drive.py',metric.OUT.parent/'gr-subcell-heat-drive/preflight-manifest.json',
        ROOT/'verification/audit_structured_enthalpy.py',ROOT/'outputs/gr-subcell-heat-drive33-mode-inspect.log']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='668bc739',
        bindings=dict(prior['bindings'],**{p.relative_to(ROOT).as_posix():sha(p) for p in files}),
        runtime_sha256=metric.runtime(),native_cells=[0,2688,5734],native_node_indices=[0,8,15],native_coordinate_gate=1e-8,
        audit='Check the pressure-mode tags and converted positive heat capacity on every stored 8/16-node reference. Compare the chain-rule-converted mode1 derivatives with fresh mode2 calls at nine actual states. Reassess the four archived exact factors and twelve saved inverse residuals under the corrected physical heat-capacity scale, retaining every original file and original verdict.',
        scope='Coordinate interpretation and finite native controls. No physical EOS, continuous derivative, exact native-Gibbs, global inverse or finite GR certificate.'))
    save('symbolic.json',symbolic())


def bindings():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert sha(ROOT/rel)==digest,rel
    for path,digest in plan['runtime_sha256'].items():assert sha(Path(path))==digest,path
    return plan


def run():
    plan=bindings();assert not (OUT/'result.json').exists()
    reference=json.loads((metric.OUT/'plan.json').read_text());state=dict(np.load(metric.g.OUT/'initial-state-17-4.npz'))
    counts={8:0,16:0};minimum=ld(np.inf);ratios=[]
    for n in [8,16]:
        for rel in sorted(set(reference['cell_sources'].values())):
            path=ROOT/rel;data=metric.node_file(str(path.with_name(path.stem+f'-nodes-{n}.npz')));a=data['eos'].astype(ld)
            assert np.all(a[...,5]==1) and np.all(a[...,6]==0)
            *_,cv=pressure_to_density(a[...,7],a[...,8],a[...,9],a[...,10])
            assert np.all(np.isfinite(cv)) and np.all(cv>0)
            counts[n]+=cv.size;minimum=min(minimum,np.min(cv));ratios.extend([float(np.min(a[...,10]/cv)),float(np.max(a[...,10]/cv))])
    assert counts=={8:5735*8,16:5735*16}
    eos=metric.g.EOS();native=[]
    for i in plan['native_cells']:
        path=ROOT/reference['cell_sources'][str(i)];data=metric.node_file(str(path.with_name(path.stem+'-nodes-16.npz')));j=list(data['cells']).index(i)
        for k in plan['native_node_indices']:
            p=eos(1,float(data['logP'][j,k]),float(data['lnT'][j,k]),state['X'][i])
            d=eos(2,float(np.log(p[0])),float(data['lnT'][j,k]),state['X'][i])
            assert p[5]==1 and p[6]==0 and d[7]==1 and d[8]==0
            converted=np.array(pressure_to_density(*map(ld,p[7:11])))
            target=d[[5,6,9,10]].astype(ld);scale=np.array([max(abs(target[0]),1),max(abs(target[1]),1),max(abs(target[3]),1),max(abs(target[3]),1)],dtype=ld)
            error=float(np.max(abs(converted-target)/scale));assert error<=plan['native_coordinate_gate'],(i,k,error)
            native.append(dict(cell=i,node=k,relative_coordinate_error=error,passed=True,
                pressure_mode_derivatives=p[5:11].tolist(),density_mode_derivatives=d[5:11].tolist()))
    factors=[];recoveries=[];p=json.loads((old.inverse.OUT/'plan.json').read_text())
    original_factors={r['cell']:r for r in json.loads((old.OUT/'positive-factors.json').read_text())['rows']}
    for i in p['cells']:
        data=dict(np.load(old.inverse.reference.OUT/f'cell-{i}-tolerance-1-nodes-16.npz'));a=data['eos'].astype(ld)
        assert np.all(a[:,5]==1) and np.all(a[:,6]==0)
        *_,cv=pressure_to_density(a[:,7],a[:,8],a[:,9],a[:,10]);weights=data['coordinate_weights_cm3'].astype(ld)
        capacity=np.sum(a[:,0]*cv*weights);assert capacity>0
        rational=F(0)
        for row,cw in zip(data['eos'],data['coordinate_weights_cm3']):
            *_,value=pressure_to_density(*(F(float(x)) for x in row[7:11]));assert value>0
            rational+=F(float(row[0]))*value*F(float(cw))
        f=original_factors[i];assert F(f['B'])>0 and F(f['enthalpy_proper'])>0 and rational>0
        factors.append(dict(cell=i,original_named_Cv_coordinate=f['Cv_coordinate'],corrected_Cv_coordinate=str(rational),
            B_unchanged=f['B'],enthalpy_unchanged=f['enthalpy_proper'],corrected_exact_factor_positive=True,
            original_over_corrected_capacity=float(F(f['Cv_coordinate'])/rational)))
        denom=F(*capacity.as_integer_ratio())
        for case in range(3):
            d=dict(np.load(old.inverse.OUT/f'cell-{i}-case-{case}.npz'))
            energy=abs(F(*d['recovered_moments'][1].as_integer_ratio())-F(*d['target'][1].as_integer_ratio()))/denom
            score=max(energy,abs(F(float(d['residual'][1]))))
            recoveries.append(dict(cell=i,case=case,corrected_scaled_residual=str(score),passed=score<=F(str(p['scaled_residual_tolerance']))))
    result=dict(classification='Counterexample candidate',completed=True,pressure_mode_nodes=counts,
        minimum_converted_cvT_display=float(minimum),original_slot_over_cvT_range=[min(ratios),max(ratios)],
        native_controls=native,corrected_factors=factors,corrected_recovery_residuals=recoveries,
        all_corrected_residual_gates_passed=all(r['passed'] for r in recoveries),
        original_files_changed=False,physical_EOS_certified=False,continuous_derivative_error_certified=False,full_GR_evolution=False)
    save('result.json',result);save('manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());assert r['completed'] and len(r['native_controls'])==9
    assert r['pressure_mode_nodes']=={'8':45880,'16':91760} and len(r['corrected_factors'])==4 and len(r['corrected_recovery_residuals'])==12
    assert json.loads((OUT/'symbolic.json').read_text())==symbolic()
    print('PASS pressure/density derivative-coordinate audit; corrected archived residual gates:',r['all_corrected_residual_gates_passed'],flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
