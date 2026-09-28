"""Trace the actual conserved energy into the frozen GR source, before repair."""
from pathlib import Path
import json,signal,time
import numpy as np
import def_native_characteristic_gr as gr

flow=gr.base.flow;OUT=gr.OUT.parent.parent/'def-native-energy-source-reconciliation'
C=flow.C;LD=np.longdouble;write=flow.write;sha=flow.sha


def main():
    assert not OUT.exists();OUT.mkdir();start=time.monotonic()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='414e242ad',previous_goal_turn='progress',
        claim='Locate the0.154percent mismatch in actual conserved-state energy, rest-mass restoration, photon subtraction, subcell lapse, or sparse-time port quadrature before changing any source.',
        decision='Repair the actual failing conversion and propagate the corrected source to the charge. A fitted energy offset or relaxed tolerance is not a repair.',
        reuse='Saved64/128 original states and their exact integrated boundary/escape/discard ledgers. No fluid steps, EOS roots or physical resolution changes in this audit.',
        budget_seconds=35,CPU_threads=1,memory_GB=2,
        stop='If missing stage ledgers require replay, design and time that bounded replay separately; do not infer an unrecorded history from the endpoint conservation residual.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(flow.__file__),Path(gr.base.__file__),flow.OUT/'coupled-64.npz',flow.OUT/'coupled-128.npz',gr.base.OUT/'source-64.npz',gr.base.OUT/'source-128.npz']}))
    signal.signal(signal.SIGALRM,flow.old.optical.timeout);signal.alarm(35);m=flow.Coupled();b=m.bulk;f=m.flow;geo=flow.Geometry(m.m);rows=[]
    for steps in [64,128]:
        z=np.load(flow.OUT/f'coupled-{steps}.npz');d=np.load(gr.base.OUT/f'source-{steps}.npz');k=-1
        m.Pi=z['Pi'];m.h=z['h'];m.mass=m.mass0-m.h[1:]+m.h[:-1];m.set_material(flow.old.END)
        dm=-np.diff(m.h).astype(LD);a=b.d['a'].astype(LD);a0=LD(m.m.a0);cx=LD(m.cx);scale=LD(m.gas_scale)
        deep=a*(m.mass.astype(LD)*(z['u'].astype(LD)-b.u0)+dm*b.u0+m.kinetic())+(a-a0)*cx*LD(C)**2*dm
        deltaU=z['U'].astype(LD)-f.initial.astype(LD)
        atmo=deltaU[2]*m.m.vol*scale;db=deltaU[0]*m.m.vol*scale/LD(C)**2
        photons=np.r_[np.sum((z['bulk_I'].astype(LD)-b.initial)*b.photon_energy_weight,axis=(1,2)),
            np.sum((z['I'].astype(LD).sum(0)-m.initial_I)*m.energy_weight,axis=(1,2))]
        reduced=np.r_[deep,atmo]+photons;full=reduced+a0*cx*LD(C)**2*np.r_[dm,db]
        rest=d['baryon_g'].astype(LD)*LD(d['cx'])*LD(C)**2
        exported=(rest+d['gas_nonrest_energy_erg']+d['photon_energy_erg'])*d['a']
        photon_exported=d['photon_energy_erg'][-1]*d['a']
        discard=float(z['discard'][2]*scale);rest_discard=float(z['discard'][0]*scale*a0*cx)
        escape=float(z['ledger'][5]+z['scalar_deep_escape']);exact=float(z['boundary'][0]);sampled=float(d['inner_cumulative_energy_erg'][-1]-d['outer_cumulative_energy_erg'][-1]);norm=float(d['inner_cumulative_energy_erg'][-1])
        q=flow.initial.Quadrature(d['edges'],8);_,_,aa,BB,_=geo(q.r.ravel()-m.m.RJ)
        measure=q.r*q.r*BB.reshape(q.r.shape);mean=((measure*aa.reshape(q.r.shape))@q.w)/(measure@q.w)
        full_difference=exported[-1]-full
        row=dict(steps=steps,scale_inner_erg=norm,exact_integrated_net_photon_port=exact,sampled_net_photon_port=sampled,
            sparse_port_error=(exact-sampled)/norm,conserved_reduced_energy=float(np.sum(reduced,dtype=LD)),
            conserved_reduced_balance=float((np.sum(reduced,dtype=LD)+discard+escape-exact)/norm),
            restored_baryon_energy=float(a0*cx*LD(C)**2*np.sum(np.r_[dm,db],dtype=LD)),
            restored_baryon_energy_relative=float(a0*cx*LD(C)**2*np.sum(np.r_[dm,db],dtype=LD)/norm),
            baryon_with_discard_g=float(np.sum(np.r_[dm,db],dtype=LD)+LD(z['discard'][0])*scale/LD(C)**2),
            total_source_conversion_error=float(np.sum(full_difference,dtype=LD)/norm),
            photon_subtraction_error=float(np.sum(photon_exported-photons,dtype=LD)/norm),
            subcell_mean_lapse_change=float(np.sum((rest[-1]+d['gas_nonrest_energy_erg'][-1]+d['photon_energy_erg'][-1])*(mean-d['a']),dtype=LD)/norm),
            actual_bulk_cx=m.cx,exported_cx=float(d['cx']),atmosphere_cx=float(f.eos.cx),discard_energy_erg=discard,rest_discard_energy_erg=rest_discard,
            spectral_escape_erg=escape,exported_port_error=float((np.sum(exported[-1],dtype=LD)-sampled)/norm))
        rows.append(row);np.savez_compressed(OUT/f'components-{steps}.npz',full_conserved_Killing=full,exported_Killing=exported[-1],conversion_difference=full_difference,
            photon_difference=photon_exported-photons,baryon_g=np.r_[dm,db]);print(json.dumps(row),flush=True)
    result=dict(classification='Counterexample candidate',paths=rows,seconds=time.monotonic()-start,source_repaired=False,final_charge_solved=False,full_goal_complete=False)
    write(OUT/'diagnosis.json',result);signal.alarm(0)


if __name__=='__main__':main()
