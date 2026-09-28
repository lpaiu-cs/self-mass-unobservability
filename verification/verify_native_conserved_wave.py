"""Audit premises and interval meaning of the saved potential response."""
from pathlib import Path
import hashlib
import json
import time
import numpy as np
import mpmath as mp
import def_native_conserved_wave as task


def main():
    out=task.OUT;assert not (out/'audit.json').exists();start=time.monotonic()
    # Recover the exact unexecuted preflight source, whose sole subsequent
    # change included the reference tail point in the absolute source norm.
    initial=json.loads((out/'preflight-plan.json').read_text());raw=Path(task.__file__).read_bytes();newline=b'\r\n' if b'\r\n' in raw else b'\n'
    block="""        # The centered rest readout is equivalent to placing the unresolved
        # reference mass at x=0. Include that mathematical point source in
        # the norm; its physical displacement and stress remain unclosed.
        ab+=body_w*(portTV+I(rows[channel]['maximum_unresolved_baryon_g']))
        an+=body_w*portTV*I(abs(rows[channel]['bulk']['trace_nonrest_erg_per_g']))"""
    before="        ab+=body_w*portTV;an+=body_w*portTV*I(abs(rows[channel]['bulk']['trace_nonrest_erg_per_g']))"
    reconstruct=raw.replace(block.encode().replace(b'\n',newline),before.encode())
    expected=initial['bindings']['verification/def_native_conserved_wave.py']
    assert hashlib.sha256(reconstruct).hexdigest()==expected
    (out/'preflight-source.py').write_bytes(reconstruct)
    plan=json.loads((out/'plan.json').read_text())
    for p,h in plan['bindings'].items():assert task.prior.old.photons.digest(task.prior.old.ROOT/p)==h,p
    model=task.old.adapter(1792);bg=model.bg;wave=np.load(out/'wave-1792.npz');bound=json.loads((out/'bound-1792.json').read_text());row=json.loads((out/'wave-1792.json').read_text())
    qmin=row['fine_quadrature']['qmin_cm'];rmin=bg.R+qmin;idx=max(0,int(np.searchsorted(bg.d['radius_cm'],rmin))-1)
    d=bg.d;region=slice(idx,None)
    assert np.all(d['energy_cgs'][region]>=0) and np.all(d['pressure_cgs'][region]>=0)
    assert np.all((d['lapse'][region]>0)&(d['lapse'][region]<=1))
    assert np.all(d['mass_geom_cm'][region]>=0) and max(d['mass_geom_cm'][region])<=bg.M*bg.R
    assert np.all(d['gamma1'][region]>0)
    # Independent dense coefficient quadrature only checks the implementation.
    # The enclosure itself uses node extrema and the analytic vacuum bound.
    r=np.linspace(rmin,bg.R,4097);V=task.coeff(bg,r)[0];z=bg.sample(r/bg.R);cc=z['N']*np.sqrt(1-2*z['m']/(r/bg.R))
    inner=float(np.trapz(abs(V)/cc,r))
    x,w=np.polynomial.legendre.leggauss(128);zz=(x+1)/2;rr=bg.R/zz;mm,nn,bb,cv,phi=bg.metric(rr/bg.R)
    vv=2*nn*nn*(mm*bg.R/rr**3-(phi/bg.R)**2);outer=float(np.sum(w/2*abs(vv)*bg.R/(zz*zz*cv)))
    actual_integral=inner+outer
    assert actual_integral<bound['potential_absolute_integral_bound_per_cm']
    assert max(abs(wave['outgoing_free_cm']))<bound['free_field_norm_bound_cm']
    assert max(abs(wave['incoming_free_cm']))<bound['free_field_norm_bound_cm']
    mp.mp.dps=60;controls=[]
    for k in [mp.mpf('-.2'),mp.mpf('.2')]:
        exact=2/k-2*(-mp.expm1(-k))/k**2;born=1-k/3
        assert abs(exact-born)<=k*k/(1-abs(k))
        controls.append(dict(point_potential_strength=float(k),exact=float(exact),first_Born=float(born),remainder=float(abs(exact-born)),bound=float(k*k/(1-abs(k)))))
    M=bg.M*bg.R;free=-float(wave['outgoing_free_cm'][-1])/M;radius=bound['all_orders_change_normalized_bound']
    interval=[float(np.nextafter(free-radius,-np.inf)),float(np.nextafter(free+radius,np.inf))]
    assert interval[1]<0
    result=dict(classification='Counterexample candidate',passed=True,seconds=time.monotonic()-start,
        source_bindings_verified=True,preflight_source_recovered_exactly=True,
        coefficient_premises_checked=True,numerical_absolute_potential_integral_per_cm=actual_integral,
        analytic_integral_enclosure_margin=float(bound['potential_absolute_integral_bound_per_cm']/actual_integral-1),
        manufactured_positive_and_negative_potential=controls,
        declared_1792_source_endpoint_potential_interval=interval,
        interval_scope='All orders of the stated scalar potential, including the algebraic fixed-coordinate-baryon adiabatic volume response. The reference tail point implicit in the centered rest readout is included in the norm. This interval is not a continuum/EOS/whole-star/observational uncertainty interval.',
        archived_metadata_unit_correction='Phase102 mass-source NPZ key scalar_source_coefficient_cm3 contains KJ with units cm^-2. Numerical use was correct; Phase103 uses coefficient_cm_minus2. Preserve the archived key and bytes.',
        compact_first_Born_not_evaluated=True,first_Born_quadrature_difference_is_not_a_rigorous_quadrature_error_bound=True,
        full_goal_complete=False)
    task.write(out/'audit.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
