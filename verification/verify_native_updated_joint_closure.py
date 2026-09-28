"""Audit common saved times and all material endpoints, without evolution."""
from pathlib import Path
import json
import numpy as np
import sympy as sp
import verify_native_updated_joint_feedback as run


def main():
    out=run.OUT;assert json.loads((out/'audit.json').read_text())['passed']
    # Proven only for the fixed-coefficient, invertible energy coordinate.
    et,B,k,se,sb,c1,c2,h,g=sp.symbols('et B k se sb c1 c2 h g')
    E=et+k*B
    assert sp.expand(E.subs({et:et+h*(se-k*sb),B:B+h*sb})-E-h*se)==0
    assert sp.expand(h*((1-g)*c1+g*c2)+h*((1-g)*(-c1)+g*(-c2)))==0
    checks=[];source_error=0.;histories={}
    for n,r in run.matter.PATHS:
        p=np.load(run.photon_path(n,r));m=np.load(out/f'steps-{n}-reference-{r}.npz')
        ids=[int(np.argmin(abs(m['t']-t))) for t in p['t']]
        assert np.max(abs(p['t']-m['t'][ids]))<1e-18 and len(m['t'])==n+1
        actual=m['history_scaled'][ids][:,[2,3]]*run.matter.AMP
        declared=p['moments'][:,[1,2]]
        norm=np.maximum(np.max(np.sum(abs(actual),axis=2),axis=0),1.)
        residual=np.sum(abs(declared-actual),axis=2)/norm
        checks.append(dict(steps=n,reference=r,common_saved_times=len(ids),material_endpoints=len(m['t']),
            energy_H_relative=residual.max(0).tolist(),
            maximum_at_time=p['t'][np.argmax(residual,axis=0)].tolist()))
        histories[n,r]=m['history_scaled']
        d=np.load(run.GR/f'source-{n}-reference-{r}.npz')
        stress=np.load(out/f'stress-{n}-reference-{r}.npz')['material']
        rest=d['baryon_g'].astype(np.longdouble)*np.longdouble(d['cx'])*np.longdouble(run.prior.C)**2
        total=rest+d['gas_nonrest_energy_erg']
        errors=[total-stress[:,0],d['nonrest_trace_erg']+rest-(stress[:,0]-stress[:,1]-2*stress[:,3]),
            d['nonrest_stress_erg']+rest-(stress[:,0]-stress[:,1]),d['pressure_volume_erg']-stress[:,3],
            d['metric_stress_erg']-(total+d['photon_energy_erg']-stress[:,1]-d['photon_radial_pressure_erg'])]
        source_error=max(source_error,float(max(np.max(abs(x)) for x in errors)/max(np.max(abs(stress)),1.)))
    assert source_error<1e-12
    fine=histories[128,128]
    def compare(a,b):
        return (np.max(np.sum(abs(a-b),axis=2),axis=0)/np.maximum(np.max(np.sum(abs(b),axis=2),axis=0),1.)).tolist()
    comparisons=dict(time_at_all_65_common_endpoints=compare(histories[64,128],fine[::2]),
        background_at_all_129_endpoints=compare(histories[128,64],fine))
    assert max(v for values in comparisons.values() for v in values)<.02
    # The old and new normalized residuals use different material norms.
    # Compare their absolute numerators at the same17 times as well.
    absolute=[]
    for photon_folder,material_folder in [(run.prior.run.RESPONSE,run.matter.OUT),(run.joint.OUT,out)]:
        p=np.load(photon_folder/'steps-128-reference-128.npz')
        m=np.load(material_folder/'steps-128-reference-128.npz')
        times=np.linspace(0,m['t'][-1],17)
        pi=[int(np.argmin(abs(p['t']-t))) for t in times]
        mi=[int(np.argmin(abs(m['t']-t))) for t in times]
        assert max(np.max(abs(p['t'][pi]-times)),np.max(abs(m['t'][mi]-times)))<1e-18
        actual=m['history_scaled'][mi][:,[2,3]]*run.matter.AMP
        absolute.append(np.max(np.sum(abs(p['moments'][pi][:,[1,2]]-actual),axis=2),axis=0))
    result=dict(classification='Counterexample candidate',passed=True,source_identity_relative=source_error,
        common_saved_time_waveform_residual=checks,material_comparisons=comparisons,
        fine_17time_absolute_energy_erg_H_count_residual_before=absolute[0].tolist(),
        fine_17time_absolute_energy_erg_H_count_residual_after=absolute[1].tolist(),
        absolute_residual_after_over_before=(absolute[1]/absolute[0]).tolist(),
        symbolic=dict(classification='Proven',passed=True,
            scope='With fixed k, E=Etilde+k*B transfers the same baryon source exactly. Opposite collision transfers cancel under the same SDIRK weights. Neither identity bounds EOS derivatives or continuum integration error.'),
        residual_scope='Only17 canonical photon/material times, plus saved pilot endpoint where present. All65/129 material endpoints are compared separately. No missing photon-time reconstruction, continuous-time supremum, contraction constant or full coupled-error bound.',
        coupled_fixed_point_verified=False,final_charge_solved=False,
        bindings={str(p):run.sha(p) for p in [Path(__file__),out/'audit.json',out/'sources.json',run.joint.OUT/'audit.json']})
    run.write(out/'closure-audit.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
