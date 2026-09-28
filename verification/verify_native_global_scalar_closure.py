"""Conditional global first-variation bound on the current saved background.

Prove the algebra separately; do not rename a frozen-source enclosure as
an EOS, source-continuum, waveform fixed-point or nonlinear certificate.
"""
from pathlib import Path
import json,signal,time
import numpy as np
import mpmath as mp
import def_native_global_scalar_closure as task

OUT=task.OUT;write=task.write;sha=task.sha;C=task.C;G=task.G;LD=np.longdouble


def main():
    assert not (OUT/'bound.json').exists();start=time.monotonic()
    inputs=[Path(__file__),Path(task.__file__),OUT/'result.json',OUT/'check.json',
        task.compact.OUT/'fields-128-g8.npz',task.compact.base.OUT/'source-128.npz',
        task.flow.INPUT/'balanced-20.npz',task.flow.OUT/'coupled-128.npz',task.flow.OUT/'coupled-64.npz']
    write(OUT/'bound-plan.json',dict(classification='Counterexample candidate',
        claim='Bound global potential scattering, arbitrary nonnegative outward exterior forcing and the compact mass-integration-constant effect of the actual port mismatch on this initial operator. Exclude direct deep forcing only for this observer cone.',
        decision='Does the conditional scalar sign survive these remaining GR terms on the current background? Do not inherit the old background bound or certify the source error from an energy residual.',
        method='Exact binary Legendre coefficient envelopes; interval arithmetic for the declared operator. Positive-mass static vacuum inequalities outside all initial photons. Absolute source norm, no cancellation-based norm. Zero initial increments make the initial charge exactly zero.',
        gates=dict(contraction=1.,conditional_lower=0.,seconds=40),CPU_threads=1,new_EOS_calls=0,
        measured_basis='Exterior four-path production9.21s; one background load about3s, algebra/envelopes only. No evolution or new physical resolution.',
        bindings={str(p):sha(p) for p in inputs}))
    signal.signal(signal.SIGALRM,task.flow.old.optical.timeout);signal.alarm(40)
    assert json.loads((OUT/'result.json').read_text())['passed']
    m=task.Exterior();rsp=m.response;bg=m.bg;d=dict(np.load(task.compact.base.OUT/'source-128.npz'))
    rsp.setup(d,8);mp.iv.dps=40;iv=mp.iv;I=iv.mpf
    def binary(v):
        a,b=v.as_integer_ratio();return I(a)/I(b)
    def up(v):return float(np.nextafter(float(v.b),np.inf))
    def lo(v):return float(np.nextafter(float(v.a),-np.inf))
    cc=binary(float(C));gg=binary(float(G));T=binary(m.T);M=binary(m.M);K=binary(abs(m.K));r0=binary(m.r0)
    # With m>=0,N<=1 the past cone cannot extend farther inward than c*T.
    rmin=r0-cc*T;selected=bg.edges[1:]>=lo(rmin);bounds={}
    for key,co in bg.co.items():
        lower=[];upper=[];absolute=[]
        for row in co[selected]:
            tail=sum((abs(binary(v)) for v in row[1:]),I(0));first=binary(row[0])
            lower.append(lo(first-tail));upper.append(up(first+tail));absolute.append(up(abs(first)+tail))
        bounds[key]=dict(lower=min(lower),upper=max(upper),absolute=max(absolute))
    assert bounds['mass']['lower']>0 and bounds['lapse']['upper']<1
    assert lo(rmin)>rsp.faces[0]
    Nmin=I(bounds['lapse']['lower']);massmax=I(bounds['mass']['upper']);bmin=1-2*massmax/rmin
    assert lo(Nmin)>0 and lo(bmin)>0 and up(bmin)<=1
    Phi=I(bounds['Phi']['absolute']);alpha=4*I(bounds['phi']['absolute']);factor=gg/cc**4
    Eg,Pg,Er,Pr=[I(bounds[k]['absolute'])*factor for k in ['gas_energy','gas_pressure','photon_energy','photon_radial_pressure']]
    Kg=I(float(np.max(abs(rsp.gamma))))*Pg;R4=I(float(np.max(abs(rsp.ratio4))))*Er
    E=Eg+Er;P=Pg+Pr;H=Eg+Pg
    kj=2*Phi/(rmin*bmin)*(1+4*iv.pi*r0*r0*(E+P))+8*iv.pi*alpha*(Eg+3*Pg)/bmin
    dp=4*iv.pi*r0*(r0*Phi*alpha*(3*(H+Kg)+4*(Er+Pr))+3*alpha**2*(H+3*Kg))
    dl=4*iv.pi*r0*(r0*Phi*(E+P+Kg+3*Pr+R4)+alpha*(H+3*Kg))
    fc=kj*r0*r0*Phi+16*iv.pi*alpha*r0*r0*Phi*(E+P)+4*iv.pi*r0*(4+4*alpha**2)*(Eg+3*Pg)+dp+dl*r0*Phi
    V=2*massmax/rmin**3+4*iv.pi*(E+P)+fc/rmin;ke=kj+dl/(rmin*bmin)
    b0=1-2*M/r0;N0=binary(m.N0);c0=N0*iv.sqrt(b0)
    # Exact vacuum: m<=M, N<=1, c=N*sqrt(b)>=c0 and Phi=K/(c*r^2).
    Ivac=(M/r0**2+2*K*K/(3*c0*c0*r0**3))/iv.sqrt(b0)
    length=(r0-rmin)/(Nmin*iv.sqrt(bmin));eta=cc*T/2*(length*V+Ivac)
    assert up(eta)<1
    # Bound every spatial polynomial and time knot, not only the final trace.
    density=rsp.source.reshape(len(rsp.t),-1,8)/rsp.dx.reshape(-1,8)
    source_co=np.einsum('cij,tcj->tci',rsp.inverse,density)
    np.savez_compressed(OUT/'enclosure-inputs.npz',source_legendre=source_co,
        **{key:value[selected] for key,value in bg.co.items()})
    integral=I(0)
    for j in range(len(rsp.half)):
        if rsp.faces[j+1]<lo(rmin):continue
        maximum=max(up(sum((abs(binary(v)) for v in row),I(0))) for row in source_co[:,j])
        integral+=I(maximum)*2*binary(rsp.half[j])
    compact_norm=cc*T/2*integral
    ext=np.load(OUT/'exterior-128-a8-r8.npz');emitted=binary(float(ext['emitted_energy_erg']));e=gg/cc**4*emitted
    kappa=M/(r0*b0)+K*K/(2*c0*c0*r0*r0);assert up(kappa)<1
    ext_stress=e*K*iv.pi/(4*c0*c0*r0*iv.sqrt(1-kappa))
    ext_mass=cc*T/2*e*K/(c0*b0*r0*r0)
    # Use the same finite-volume lapse average as Response.setup, not the
    # center-lapse approximation of the older port diagnostic.
    q=task.flow.initial.Quadrature(d['edges'],8);_,_,a,B,_=rsp.geo(q.r.ravel()-rsp.model.m.RJ)
    measure=q.r*q.r*B.reshape(q.r.shape);mean_a=((measure*a.reshape(q.r.shape))@q.w)/(measure@q.w)
    rest=np.asarray(d['baryon_g'],LD)*LD(d['cx'])*LD(C)**2
    energy=rest+d['gas_nonrest_energy_erg']+d['photon_energy_erg'];ek=np.asarray((energy*mean_a).sum(1,dtype=LD),float)
    dt=np.diff(d['t']);lum=d['inner_luminosity']-d['outer_luminosity']
    mismatch=ek-d['inner_cumulative_energy_erg']+d['outer_cumulative_energy_erg']
    a2=-np.diff(lum)/(2*dt);a1=np.diff(ek)/dt-lum[:-1];a0=mismatch[:-1]
    root=np.divide(-a1,2*a2,out=np.zeros_like(a1),where=a2!=0);valid=(root>0)&(root<dt)
    port=max(np.max(abs(mismatch)),np.max(abs((a2*root*root+a1*root+a0)[valid]),initial=0))
    # A constant-in-r adjustment of J may enforce pairing, but is not a
    # replacement for the actual trace/stress/energy discretization error.
    port_effect=cc*T/2*length*ke*gg/cc**4/Nmin*binary(float(port))
    potential=eta/(1-eta)*(compact_norm+ext_stress+ext_mass+port_effect)
    fields=np.load(task.compact.OUT/'fields-128-g8.npz');free=-float(fields['direct_and_mass_stress_U'][-1,-1])/m.M
    assert not np.any(fields['direct_and_mass_stress_U'][0]) and ext['normalized_exterior'][0]==0
    arbitrary=(ext_stress+ext_mass+port_effect+potential)/M
    scalar_lower=I(free)-arbitrary;scalar_upper=I(free)+arbitrary
    epsilon=gg/cc**4*binary(float(ext['arrived_energy_erg'][-1]))/M;alpha0=K/M
    nominal=(I(free)+binary(float(ext['normalized_exterior'][-1]))+alpha0*epsilon)/(1-epsilon)
    conditional_error=(port_effect+potential)/M/(1-epsilon)
    samples=rsp.coeff(np.linspace(lo(rmin),m.r0,2049));assert np.max(abs(samples['V']))<=up(V)
    assert np.max(abs(samples['K']))<=up(ke)
    assert abs(ext['normalized_exterior'][-1])<=up((ext_stress+ext_mass)/M)
    row=dict(classification='Counterexample candidate',passed=lo(scalar_lower)>0,
        causal_radius_min_cm=lo(rmin),represented_inner_radius_cm=float(rsp.faces[0]),deep_direct_source_excluded=True,
        frozen_polynomial_bounds=bounds,potential_absolute_per_cm2=up(V),mass_source_coefficient_absolute_per_cm2=up(ke),
        global_potential_contraction=up(eta),compact_absolute_free_source_norm_cm=up(compact_norm),
        all_orders_potential_normalized_bound=up(potential/M),
        arbitrary_outward_exterior_stress_bound=up(ext_stress/M),arbitrary_outward_exterior_mass_bound=up(ext_mass/M),
        exact_mean_lapse_port_mismatch_erg=float(port),port_over_inner_energy=float(port/np.max(abs(d['inner_cumulative_energy_erg']))),
        compact_mass_constant_port_effect_bound=up(port_effect/M),
        free_compact_endpoint=free,arbitrary_outward_conditional_scalar_interval=[lo(scalar_lower),up(scalar_upper)],
        nominal_scalar_plus_arrived_photon_mass=lo(nominal),
        conditional_GR_only_interval=[lo(nominal-conditional_error),up(nominal+conditional_error)],
        nominal_exterior_to_compact=float(ext['normalized_exterior'][-1]/free),
        potential_coefficient_dense_check=True,arbitrary_angle_bound_contains_computed_exterior=True,
        initial_increment_exactly_zero=True,continuous_mass_port_extrema_checked=True,
        scope='The declared saved linear source, finite spatial/time representation and exact vacuum idealization. Legendre envelopes use exact binary coefficients and outward40-digit interval arithmetic; this does not enclose their reconstruction from continuum physics or numerical vacuum/background errors.',
        port_limit='Only the compact mass-integration-constant consequence of the measured source/port mismatch is bounded. No unpaired instantaneous exterior mass field is silently propagated; exterior packets are paired with equal body debit. Trace/stress and constitutive discretization errors remain open.',
        deep_limit='Zero initial increments/no incoming wave and this3.434ms null-infinity retarded endpoint only. Does not exclude changes of initial data, inner mass constraints, or later-time deep response.',
        full_source_error_enclosed=False,coupled_fixed_point_verified=False,nonlinear_GR=False,final_charge_solved=False,full_goal_complete=False,
        seconds=time.monotonic()-start)
    write(OUT/'bound.json',row);signal.alarm(0);print(json.dumps(row),flush=True)


if __name__=='__main__':main()
