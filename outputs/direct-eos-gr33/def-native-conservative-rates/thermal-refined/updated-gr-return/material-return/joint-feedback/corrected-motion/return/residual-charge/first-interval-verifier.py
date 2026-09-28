"""Independently audit the fine input-box certificate with 40-digit intervals."""
from pathlib import Path
import hashlib,json,signal,time
import mpmath as mp
import numpy as np

OUT=Path('outputs/direct-eos-gr33/def-native-conservative-rates/thermal-refined/updated-gr-return/material-return/joint-feedback/corrected-motion/return/residual-charge')


def main():
    result=json.loads((OUT/'enclosure.json').read_text());assert result['passed']
    remaining=int(90-result['seconds']);assert remaining>0 and not (OUT/'interval-audit.json').exists()
    paths=[Path(__file__),Path('verification/def_native_residual_charge_envelope.py'),OUT/'enclosure.json',OUT/'enclosure-inputs-128-128.npz',OUT/'source-128-128.npz',OUT/'residual-128-128.npz']
    plan=dict(classification='Counterexample candidate',claim='Recompute every fine-path pressure,source-polynomial and potential bound independently using40-digit outward intervals from exact binary inputs. The original producer used binary64 nextafter arithmetic.',
        scope='Same declared frozen coefficient map and source box. No EOS or coupled-gain certificate. The coarser/background paths retain their original binary64 enclosure only.',
        remaining_original_enclosure_seconds=remaining,bindings={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    (OUT/'interval-audit-plan.json').write_text(json.dumps(plan,indent=2)+'\n')
    def timeout(*_):raise TimeoutError('Original90s enclosure budget remainder')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(remaining);start=time.monotonic()
    mp.iv.dps=40;I=mp.iv.mpf
    def B(x):
        n,d=float(x).as_integer_ratio();return I(n)/I(d)
    def checked(stored,exact):assert B(stored).a>=exact.b,(stored,str(exact))
    data=np.load(OUT/'enclosure-inputs-128-128.npz');source=np.load(OUT/'source-128-128.npz');r=result['paths'][1]
    box=data['box'];coeff=data['pressure_coeff'];n=box.shape[1];energies=[];traces=[];stresses=[]
    for j in range(n):
        e=B(box[0,j])/B(data['reference_lapse'][j]);h=B(box[1,j]);energies.append(e)
        ts=[];ss=[]
        for k in range(17):
            pp=[abs(B(coeff[k,p,j,0]))*B(box[0,j])+abs(B(coeff[k,p,j,1]))*h for p in [0,1]]
            ts.append((e+pp[1]+2*pp[0]).b);ss.append((e+pp[1]).b)
        traces.append(max(ts));stresses.append(max(ss))
    prefix=I(0);free_integral=I(0);potential_integral=I(0);checks=0
    for j in range(n):
        node=[]
        for q in range(8):
            i=8*j+q
            J=abs(B(data['frozen_mass_coefficient'][i]))*(prefix+energies[j]*abs(B(data['frozen_partial_lapse'][j,q])))
            v=(abs(B(data['frozen_direct_coefficient'][i]))*traces[j]+abs(B(data['frozen_stress_coefficient'][i]))*stresses[j]+abs(B(data['frozen_mass_source_coefficient'][i]))*J)
            node.append(v/B(data['frozen_dx'][i]))
        spatial=[];potential=[]
        for p in range(8):
            sc=sum((abs(B(data['frozen_inverse'][j,p,q]))*node[q] for q in range(8)),I(0))
            vc=sum((abs(B(data['frozen_inverse'][j,p,q]))*abs(B(data['frozen_potential'][8*j+q])) for q in range(8)),I(0))
            checked(data['source_coefficient_bounds'][j,p],sc);checked(data['potential_coefficient_bounds'][j,p],vc)
            spatial.append(sc);potential.append(vc);checks+=2
        measure=B(data['cell_positive_quadrature_measure'][j])
        free_integral+=measure*sum(spatial,I(0));potential_integral+=measure*sum(potential,I(0))
        prefix+=energies[j]*abs(B(data['frozen_mean_lapse'][j]))
    factor=B(29979245800.)*B(source['t'][-1])/2
    free=factor*free_integral;eta=factor*potential_integral;bound=free/(1-eta)/B(source['M_cm'])
    checked(r['free_norm_cm'],free);checked(r['compact_potential_contraction'],eta);checked(r['normalized_box_bound'],bound)
    # Report unmixed E/H fractions even where the stored pressure map is zero.
    residual=np.load(OUT/'residual-128-128.npz')['residual'];zero=np.all(coeff==0,axis=(1,3))
    fraction=np.max(np.sum(abs(residual)*zero[:,None,:],axis=2),axis=0)/np.maximum(np.max(np.sum(abs(residual),axis=2),axis=0),1e-300)
    assert np.max(fraction)<1e-12
    row=dict(classification='Proven',passed=True,scope='Independent40-digit interval verification of all848? fine spatial/potential coefficient bounds; the exact frozen binary inputs are premises, not a certified physical EOS.',
        coefficient_inequalities_checked=checks,pressure_knots=17,physical_cells=n,
        normalized_box_interval=[str(bound.a),str(bound.b)],potential_norm_interval=[str(eta.a),str(eta.b)],
        zero_pressure_map_residual_E_H_fractions=fraction.tolist(),seconds=time.monotonic()-start,
        original_enclosure_plus_independent_audit_seconds=result['seconds']+time.monotonic()-start,
        physical_EOS_certified=False,coupled_fixed_point_verified=False,final_charge_solved=False)
    row['scope']='Independent40-digit interval verification of every fine spatial/potential coefficient bound. Exact frozen binary inputs and the declared linear reconstruction are premises; physical EOS and dynamical error remain outside the certificate.'
    (OUT/'interval-audit.json').write_text(json.dumps(row,indent=2)+'\n');signal.alarm(0);print(json.dumps(row),flush=True)


if __name__=='__main__':main()
