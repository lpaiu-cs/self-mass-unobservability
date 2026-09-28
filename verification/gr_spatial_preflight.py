"""Locate spatial/background and surface limits before a GR time integrator.

These are diagnostics of saved point samples. They do not promote a midpoint
profile to a conservative cell-average state or rescue failed time refinement.
"""
import json,sys
from fractions import Fraction as F
import numpy as np
import gr_heat_entropy_closure as heat

g=heat.g;OUT=g.OUT/'gr-spatial-preflight'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def five_point(x,y):
    result=[];conditions=[];n=len(x)
    for i in range(n):
        start=max(0,min(i-2,n-5));idx=np.arange(start,start+5)
        offset=x[idx]-x[i];scale=max(abs(offset));z=offset/scale
        matrix=np.vander(z,5,increasing=True)
        weights=np.linalg.solve(matrix.T,np.array([0.,1.,0.,0.,0.]))
        result.append((weights@(y[idx]/y[i]-1))*y[i]/scale)
        conditions.append(np.linalg.cond(matrix))
    return np.array(result),np.array(conditions)


def run():
    assert not OUT.exists();OUT.mkdir();heat.verify()
    paths=[g.ROOT/'verification/gr_spatial_preflight.py',heat.OUT/'manifest.json',heat.OUT/'root-brackets.json',
        heat.OUT/'initial-rates.npz',g.OUT/'initial-state-17-4.npz',g.OUT/'initial-structure-17-4.json',
        g.OUT/'gr-microphysics/auxiliaries.npz',g.OUT/'gr-heat-initial-constraints/initial-data.npz',
        g.OUT/'gr-heat-relaxation-boundary/thresholds.npz',g.OUT/'gr-radiative-boundary/diagnostics.npz',
        g.OUT/'gr-caloric-increment/path-1.npz',g.OUT/'gr-caloric-increment/path-2.npz',
        g.OUT/'gr-caloric-increment/path-2.json',g.OUT/'gr-nonlinear-thermal/duration.json']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='c8b8dad',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        diagnostic_defined_after_time_failure=True,
        scope='Compare independently differenced point pressure to the TOV force, point values to shell integrals, and the existing failed maximum temperature gate to spatial/optical diagnostics. No new pass tolerance, claimed grid convergence, interval derivative bound or physical surface completion.',
        full_GR_evolution=False,physical_EOS_certified=False))
    x=np.array([0.,.03,.2,.7,.9,1.1,1.7,2.]);y=2+.3*x+.1*x*x-.04*x**3+.001*x**4
    derivative,condition=five_point(x,y);expected=.3+.2*x-.12*x*x+.004*x**3
    control=float(abs(derivative-expected).max());assert control<1e-10
    state,aux=heat.old.char.prior.micro.inputs();eos=aux['eos'];dm=state['dm'];rho=np.exp(state['lnd'])
    initial=dict(np.load(g.OUT/'gr-heat-initial-constraints/initial-data.npz'))
    transport=dict(np.load(g.OUT/'gr-heat-relaxation-boundary/thresholds.npz'))
    optical=dict(np.load(g.OUT/'gr-radiative-boundary/diagnostics.npz'))
    r=state['r_mid_m']*100;rf=state['radius_faces_m']*100;mf=state['mass_faces_geom']*100
    c=g.c.gr.C*100;G=g.c.gr.G*1000;N=initial['lapse'];a=initial['radial_metric_factor']
    w=transport['enthalpy_density'];P=eos[:,1];m=initial['m_geom_cm'];energy=w-P
    gravity=c*c*(m+4*np.pi*r**3*G*P/c**4)*a/r**2
    exact_proper_dp=-w*gravity/c**2
    dp3=np.gradient(P,r,edge_order=2)/a;dp5,conditioning=five_point(r,P);dp5/=a
    imbalance=np.array([-c*N*(dp-exact_proper_dp)/w for dp in [dp3,dp5]])
    relative=np.array([abs(dp/exact_proper_dp-1) for dp in [dp3,dp5]])
    # A midpoint value is not a shell average. Factor the cubic difference
    # to avoid cancellation between two large almost equal radii.
    dr=rf[:-1]-rf[1:];V=4*np.pi/3*dr*(rf[:-1]**2+rf[:-1]*rf[1:]+rf[1:]**2)
    baryon_guess=rho*a*V;mass_guess=G*energy*V/c**4;mass_difference=mf[:-1]-mf[1:]
    resolved=mass_difference>64*np.maximum(np.spacing(abs(mf[:-1])),np.spacing(abs(mf[1:])))
    mass_ratio=np.full(len(r),np.nan);mass_ratio[resolved]=mass_guess[resolved]/mass_difference[resolved]-1
    speed=np.zeros((2,len(r)))
    for row in json.loads((heat.OUT/'root-brackets.json').read_text())['rows']:
        speed[row['model'],row['cell']]=max(float(abs(F(x))) for pair in row['bounds'] for x in pair)
    crossing=dr[None,:]/(c*N[None,:]/a[None,:]*speed)
    path1=dict(np.load(g.OUT/'gr-caloric-increment/path-1.npz'))
    path2=dict(np.load(g.OUT/'gr-caloric-increment/path-2.npz'))
    error=abs(path2['lnT_shift']-path1['lnT_shift']);bad=error>1e-4;worst=int(np.argmax(error))
    previous=json.loads((g.OUT/'gr-caloric-increment/path-2.json').read_text())
    assert abs(float(error.max())-previous['time_refinement_logT_difference'])<1e-18
    duration=json.loads((g.OUT/'gr-nonlinear-thermal/duration.json').read_text())['coordinate_seconds']
    full_rates=dict(np.load(heat.OUT/'initial-rates.npz'))['values'];heat_acc=np.zeros((2,len(r)))
    for row in full_rates:heat_acc[int(row[1]),int(row[0])]=row[2]
    top=[]
    for i in np.argsort(error)[-8:][::-1]:
        top.append(dict(cell=int(i),time_logT_difference=float(error[i]),rho_cgs=float(rho[i]),
            temperature_K=float(np.exp(state['lnT'][i])),baryon_cell_fraction=float(dm[i]/dm.sum()),
            radius_cm=float(r[i]),radial_cell_width_cm=float(dr[i]),
            optical_depth_from_stored_boundary=float(optical['optical_depth_from_saved_boundary'][i]),
            Rosseland_Knudsen=float(optical['Rosseland_transport_Knudsen'][i]),
            pressure_balance_relative_3point=float(relative[0,i]),pressure_balance_relative_5point=float(relative[1,i])))
    np.savez_compressed(OUT/'diagnostics.npz',pressure_imbalance_vdot=imbalance,
        pressure_derivative_relative_error=relative,pressure_stencil_condition=conditioning,
        point_to_cell_baryon_relative_difference=baryon_guess/dm-1,
        point_to_cell_geometric_mass_relative_difference=mass_ratio,geometric_mass_difference_resolved=resolved,
        characteristic_cell_crossing_coordinate_seconds=crossing,
        time_refinement_logT_difference=error,heat_only_initial_vdot=heat_acc)
    result=dict(classification='Counterexample candidate',completed=True,cells=len(r),
        manufactured_quartic_derivative_absolute_error=control,
        maximum_pressure_stencil_condition=float(conditioning.max()),
        maximum_pressure_derivative_relative_errors=[float(v.max()) for v in relative],
        maximum_abs_pressure_imbalance_vdot=[float(abs(v).max()) for v in imbalance],
        maximum_abs_heat_only_initial_vdot=[float(abs(v).max()) for v in heat_acc],
        mass_fraction_5point_pressure_acceleration_larger_than_heat_only=[float(dm@(abs(imbalance[1])>abs(v))/dm.sum()) for v in heat_acc],
        maximum_abs_point_to_cell_baryon_relative_difference=float(abs(baryon_guess/dm-1).max()),
        total_point_to_cell_baryon_relative_difference=float((baryon_guess.sum()-dm.sum())/dm.sum()),
        geometric_mass_cells_above_64ulp_difference=int(resolved.sum()),
        maximum_resolved_point_to_cell_mass_relative_difference=float(np.nanmax(abs(mass_ratio))),
        minimum_pointwise_characteristic_cell_crossing_seconds=[float(v.min()) for v in crossing],
        duration_over_minimum_cell_crossing=[float(duration/v.min()) for v in crossing],
        time_gate_failed_cells=int(bad.sum()),time_gate_failed_baryon_fraction=float(dm@bad/dm.sum()),
        largest_time_errors=top,outermost_midpoint_pressure_dyn_cm2=float(P[0]),
        conclusions='Pointwise radial differences and midpoint-to-average substitutions must not be read as the physical initial acceleration or as a conservative GR state. A well-balanced spatial discretization and consistent cell-average reconstruction must be validated before finite evolution. Failure confined to any small mass region still fails the original maximum-norm time gate.',
        physical_boundary_or_spatial_error_certified=False,full_GR_evolution=False)
    save('result.json',result)
    save('manifest.json',dict(classification='Counterexample candidate',
        sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    print('GR SPATIAL PREFLIGHT',result,flush=True);verify()


def verify():
    for rel,digest in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'result.json').read_text())['completed']
    print('PASS spatial preflight bindings; diagnostics are not a GR evolution pass',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
