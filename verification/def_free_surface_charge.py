"""Paired outgoing mechanical contribution and thermal-background eligibility."""
from pathlib import Path
import json
import time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
import def_free_surface_response_normalized as repaired

old=repaired.old
h=old.h
OUT=repaired.OUT


def clamped(harmonic,refinement):
    full=np.load(OUT/f'harmonic-{harmonic}-grid-{refinement}.npz');grid=full['grid'];bg=repaired.namespace['sample'](grid)
    row=json.loads((OUT/f'harmonic-{harmonic}-grid-{refinement}.json').read_text());omega=row['omega_R0_over_c']
    fn,_=old.symbolic();r,m,p,en,phi,v,N,ga=[bg[k] for k in ['r','m','p','e','phi','v','N','gamma']];b=1-2*m/r;alpha=-4*phi
    D=fn(r,m,p,en,phi,v,-4.);Fder=D[8:];compression=-(r*v+3*alpha)
    coefficient=Fder[1]*r*r*b*v+Fder[2]*ga*p*compression+Fder[3]*(en+p)*compression+Fder[4]-omega*omega/(N*N*b)
    n=len(r);A=np.zeros((n,2,2));A[:,0,1]=1;A[:,1,0]=coefficient;A[:,1,1]=Fder[5]
    dx=np.diff(grid);I=np.eye(2);blocks=[-I[None,:,:]-dx[:,None,None]*A/2,I[None,:,:]-dx[:,None,None]*A/2]
    rr=[];cc=[];vv=[]
    for offset,block in zip([0,2],blocks):
        rr.extend((2*np.arange(n)[:,None,None]+np.arange(2)[None,:,None]+np.zeros((1,1,2),int)).ravel())
        cc.extend((2*np.arange(n)[:,None,None]+np.arange(2)[None,None,:]+offset+np.zeros((1,2,1),int)).ravel());vv.extend(block.ravel())
    saved=np.load(OUT/'background.npz');Rs=grid[-1];Nf=saved['N'][-1];mf=saved['m'][-1];Fmetric=Nf*np.sqrt(1-2*mf/Rs)
    Z=complex(*row['impedance']);hw=complex(*row['h']);drive=np.exp(-1j*omega*Rs)/(Fmetric*hw)
    rr.extend([2*n,2*n+1,2*n+1]);cc.extend([1,2*n,2*n+1]);vv.extend([1,-Z,Rs])
    rhs=np.zeros(2*(n+1),complex);rhs[-1]=drive;matrix=coo_matrix((vv,(rr,cc)),shape=(len(rhs),len(rhs))).tocsc()
    scale=np.asarray(abs(matrix).sum(1)).ravel();answer=spsolve(matrix.multiply((1/scale)[:,None]).tocsc(),rhs/scale).reshape(n+1,2)
    residual=float(np.max(abs(matrix@answer.ravel()-rhs)/(np.asarray(abs(matrix)@abs(answer.ravel())).ravel()+abs(rhs)+1e-100)))
    assert residual<1e-9
    full_scalar=complex(*row['Eulerian_surface_scalar']);contrast=(full_scalar-answer[-1,0])/hw
    background=json.loads((h.OUT/'absolute-shoot/result-0.001.json').read_text());M=background['ADM_geom_m'];R=float(saved['R'])*Rs
    # Same incoming field: the outgoing difference avoids subtracting O(1/omega)
    # incident/outgoing amplitudes. The clamped body needs an external support.
    charge=-R/M*contrast
    np.savez_compressed(OUT/f'clamped-{harmonic}-{refinement}.npz',grid=grid,response=answer)
    return dict(harmonic=harmonic,refinement=refinement,linear_residual=residual,
        clamped_scalar_surface=[float(answer[-1,0].real),float(answer[-1,0].imag)],
        outgoing_full_minus_clamped=[float(contrast.real),float(contrast.imag)],
        delta_alpha_per_unit_incident_scalar=[float(charge.real),float(charge.imag)],
        scope='Clamped control sets xi=0 but retains baryon-volume, conformal, adiabatic EOS and scalar/metric response. The external holding force is a defined comparison, not a second freely evolving star.')


def run():
    assert not (OUT/'charge-result.json').exists();started=time.monotonic()
    paths=[Path(__file__),Path(repaired.__file__),Path(old.__file__),OUT/'background.npz',OUT/'result.json',h.OUT/'absolute-shoot/lapse.npz',h.OUT/'absolute-shoot/background-0.001.npz']
    h.write(OUT/'charge-plan.json',dict(classification='Counterexample candidate',bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths},
        claim='Extract same-incoming outgoing full-minus-clamped material motion contribution from the coupled solutions; separately test necessary thermal equilibrium, without claiming a complete orbital charge.',
        paired_subtraction='Delta outgoing=(Delta phi_surface-xi*Phi_surface contrast)/h. Static mass-denominator terms cancel in this same-background, same-asymptotic-field mechanical contrast.',
        clamped_control='xi=0, delta ln rho=-(r*Phi+3*alpha)*delta phi; same adiabatic native Gamma1 and scalar/metric constraints. Momentum is held by a prescribed support; this is an explicit decomposition.',
        thermal_test='At zero velocity and zero initial heat flux, positive steady conductive/radiative transport requires d ln(T_J*A*N)/dr=0. Test this necessary condition on the saved fine background; no cooling-time bound inferred.',
        gates=dict(linear_residual=1e-9,mechanical_charge_grid_relative=.02),budget=dict(hard_timeout_seconds=30,native_calls=0,time_steps=0)))
    rows=[clamped(k,j) for k in [0,1,2,3] for j in [1,2]];comparisons=[]
    benchmark=json.loads((old.exterior.s.OUT/'companion-benchmark.json').read_text());drive=benchmark['leading_drive']['maximum_scalar_excursion_from_inverse_semimajor_reference'];phi0=.001
    for k in [0,1,2,3]:
        pair=[r for r in rows if r['harmonic']==k];a,b=[complex(*r['delta_alpha_per_unit_incident_scalar']) for r in pair]
        static=complex(*next(r for r in rows if r['harmonic']==0 and r['refinement']==2)['delta_alpha_per_unit_incident_scalar'])
        relative=abs(a-b)/max(abs(b),1e-30)
        comparisons.append(dict(harmonic=k,mechanical_charge_gain=[b.real,b.imag],relative_grid_difference=relative,
            passed=relative<.02,benchmark_maximum_drive_scaled_delta_alpha_over_phi0=float(abs(b)*drive/phi0),
            static_subtracted_maximum_drive_scaled=float(abs(b-static)*drive/phi0),
            scaling_scope='Maximum declared scalar excursion used as a common scale, not the individual Fourier harmonic amplitude and not a certified physical charge bound.'))
    bg=np.load(h.OUT/'absolute-shoot/background-0.001.npz');lap=np.load(h.OUT/'absolute-shoot/lapse.npz');body=h.Structure(.001)
    phi=.001*(1+body.mu*bg['states'][:,3]);lnTheta=bg['thermo'][:,1]-2*phi*phi+lap['nu_mid']
    span=float(np.ptp(lnTheta));thermal=dict(classification='Counterexample candidate',zero_heat_steady_Tolman_condition_passed=span<1e-10,
        redshifted_temperature_log_span=span,maximum_to_minimum_ratio=float(np.exp(span)),
        consequence='The zero-heat mechanical background is not a stationary solution of the positive-conductivity heat model. These adiabatic frequency solves cannot be promoted to an orbit-long thermal response without evolving the background or bounding its effect.',
        cooling_time_or_dynamic_error_bound=False)
    h.write(OUT/'thermal-eligibility.json',thermal)
    result=dict(classification='Counterexample candidate',coupled_free_surface_adiabatic_response=True,rows=rows,comparisons=comparisons,
        mechanical_charge_spatial_gate_passed=all(r['passed'] for r in comparisons),thermal_background=thermal,
        seconds=time.monotonic()-started,actual_orbital_harmonics_applied=False,restricted_derivative_comparator_removed=False,
        complete_thermal_response=False,physical_radiative_atmosphere=False,joint_errors_certified=False,full_dynamic_charge_solved=False)
    h.write(OUT/'charge-result.json',result);print(json.dumps(dict(comparisons=comparisons,thermal=thermal,seconds=result['seconds'])),flush=True)


if __name__=='__main__':run()
