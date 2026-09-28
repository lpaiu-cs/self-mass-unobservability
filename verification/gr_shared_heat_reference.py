"""Shared material heat faces and a regular centre for the native GR reference."""
import json,sys
from fractions import Fraction as F
import numpy as np
import sympy as s
import gr_full_subcell_reference as reference

g=reference.g;OUT=g.OUT/'gr-shared-heat-reference'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def run():
    assert not OUT.exists();OUT.mkdir();reference.bindings()
    assert json.loads((reference.OUT/'pilot.json').read_text())['all_passed']
    paths=[g.ROOT/'verification/gr_shared_heat_reference.py',reference.OUT/'pilot-manifest.json',
        reference.OUT/'plan.json',g.OUT/'gr-increment-structure/path-4.npz',g.OUT/'initial-state-17-4.npz',
        g.OUT/'gr-transport/diagnostics.npz',g.OUT/'gr-microphysics/auxiliaries.npz']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='b84cb63',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},arithmetic_relative_tolerance=1e-12,
        reconstruction='One shared saved initial L_infinity on each material face, zero heat at centre and outer photosphere. Interpolate L_infinity linearly in enclosed baryon mass inside each cell, retaining the same face values. Q=L_infinity/(4*pi*r^2*N^2*c). No constant nonzero heat density is extended to r=0.',
        lapse='Use the historical saved face lapse with the increment-corrected face radius/mass. At the six pilot cells, reconstruct the node lapse by the isentropic enthalpy ratio relative to the saved midpoint. Finite inconsistencies between these separately computed geometries/lapses are not enclosed.',
        scope='Full shared-face initial flux inventory plus six native subcell initial heat profiles. Prescribed initial data, not a solved transport law, physical atmosphere or GR time trajectory.'))
    r,rho,L,B,N,c=s.symbols('r rho L B N c',positive=True)
    enclosed=4*s.pi*rho*r**3/3
    Q=L*enclosed/B/(4*s.pi*r*r*N*N*c)
    assert s.simplify(Q/r-L*rho/(3*B*N*N*c))==0
    assert s.limit(Q,r,0,dir='+')==0
    a,P,v=s.symbols('a P v',real=True)
    energy_flux=4*s.pi*r*r*N*c*(P*v+Q)/a
    assert s.simplify(energy_flux.subs(v,0)-L*enclosed/(B*N*a))==0
    assert s.simplify((4*s.pi*r*r*N*c*P*v/a).subs(v,0))==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        centre='For a smooth finite positive central baryon density and lapse, B_enclosed=(4*pi/3)*rho_c*r^3+O(r^5), so linear L_infinity(B) with L_infinity(0)=0 gives Q=L_out*rho_c*r/(3*B_cell*N_c^2*c)+O(r^3). The flat constant-density leading coefficient and Q(0)=0 were checked symbolically.',
        shared_flux='At v=0, the material shell-energy flux is F_E=L_infinity/(N*a), not L_infinity everywhere inside matter. For moving material add 4*pi*r^2*N*c*P*v/a. The same interface value enters the two adjoining shell updates with opposite signs.',
        centre_limit='Regular-centre conclusion assumes a smooth central solution; the seed/native EOS approximation is not certified by this leading-order identity.'))
    state=dict(np.load(g.OUT/'initial-state-17-4.npz'));grid=dict(np.load(g.OUT/'gr-increment-structure/path-4.npz'))
    diag=dict(np.load(g.OUT/'gr-transport/diagnostics.npz'));mid=np.load(g.OUT/'gr-microphysics/auxiliaries.npz')['eos']
    ld=np.longdouble;cc=ld(g.c.gr.C)*100;G=ld(g.c.gr.G)*1000
    radius=grid['radius_m']*100;mass=grid['mass_geom_m']*100;lapse=np.exp(state['nu_faces']).astype(ld)
    metric=np.ones(len(radius),dtype=ld);metric[:-1]=1/np.sqrt(1-2*mass[:-1]/radius[:-1])
    luminosity=np.r_[0.,diag['interior_Linf'],0.].astype(ld)
    assert len(luminosity)==len(radius) and luminosity[0]==luminosity[-1]==0
    flux=luminosity/(lapse*metric)
    exact_flux=[F(*value.as_integer_ratio()) for value in flux]
    exact_rates=[b-a for a,b in zip(exact_flux,exact_flux[1:])]
    assert sum(exact_rates,F())==0
    rates=flux[1:]-flux[:-1]
    Qfaces=np.zeros_like(flux);Qfaces[:-1]=luminosity[:-1]/(4*np.pi*radius[:-1]**2*lapse[:-1]**2*cc)
    control_velocity=ld('1e-8');surface_P=np.exp(grid['faces'][0,2])
    surface_work=4*np.pi*radius[0]**2*lapse[0]/metric[0]*cc*surface_P*control_velocity
    assert surface_work>0
    np.savez_compressed(OUT/'shared-faces.npz',radius_cm=radius,mass_geom_cm=mass,lapse=lapse,metric_a=metric,
        Linfinity_erg_s=luminosity,Q_erg_cm3=Qfaces,material_energy_flux_erg_s=flux,
        shell_energy_rate_erg_s=rates)
    data=dict(np.load(reference.OUT/'pilot-nodes-16.npz'));knots,_=np.polynomial.legendre.leggauss(16)
    split=len(grid['outer'])-1;n=len(state['dm']);records=[];fields=[]
    for index,i in enumerate(data['cells']):
        outside=i<split;j=int(i) if outside else n-1-int(i)
        branch=grid['outer'] if outside else grid['inner'];low=ld(0) if j==0 else branch[j,0];high=branch[j+1,0]
        if outside:fraction=(1-knots.astype(ld))/2
        else:
            left=np.cbrt(low);right=np.cbrt(high)
            q=((left+right)/2+(right-left)*knots.astype(ld)/2)**3;fraction=(q-low)/(high-low)
        assert np.all((fraction>0)&(fraction<1))
        native=data['eos'][index];r=data['radius_cm'][index];aa=data['metric_a'][index]
        rest=ld(data['C_X'][index])*cc*cc
        Hmid=ld(mid[i,2])+ld(mid[i,1])/ld(mid[i,0]);H=native[:,2].astype(ld)+native[:,1]/native[:,0]
        N=np.exp(ld(state['nu'][i])+np.log1p((Hmid-H)/(rest+H)))
        L=luminosity[i+1]+fraction*(luminosity[i]-luminosity[i+1]);heat=L/(4*np.pi*r*r*N*N*cc)
        qgeom=G*heat/cc**4;extrinsic=4*np.pi*r*aa*qgeom
        momentum_defect=2*extrinsic/r-8*np.pi*aa*qgeom
        scale=8*np.pi*aa*abs(qgeom);relative=np.divide(abs(momentum_defect),scale,out=np.zeros_like(scale),where=scale>0)
        eps=native[:,0]*(rest+native[:,2]);P=native[:,1];w=eps+P
        discriminant=w*w-4*heat*heat;assert np.all(discriminant>0)
        root=np.sqrt(discriminant);landau_energy=(eps-P+root)/2;landau_pressure=(-eps+P+root)/2
        margin=np.minimum(landau_energy-abs(landau_pressure),landau_energy-abs(P))
        record=dict(cell=int(i),nodes=16,maximum_relative_momentum_defect=float(relative.max()),
            minimum_dominant_energy_margin_over_enthalpy=float(np.min(margin/w)),
            maximum_abs_heat_over_enthalpy=float(np.max(abs(heat/w))),regular_centre_cell=bool(i==n-1))
        record['passed']=bool(relative.max()<=1e-12 and np.all(margin>=0));records.append(record)
        fields.append(dict(lapse=N,Q_erg_cm3=heat,extrinsic_per_cm=extrinsic,
            dlnmetric_dt=-cc*N*extrinsic,dmass_geom_cm_dt=-4*np.pi*r*r*N/aa*qgeom*cc,
            enclosed_cell_baryon_fraction=fraction))
    np.savez_compressed(OUT/'pilot-heat-fields.npz',cells=data['cells'],
        **{key:np.array([row[key] for row in fields]) for key in fields[0]})
    save('result.json',dict(classification='Counterexample candidate',completed=True,shared_material_faces=len(flux),
        shell_cells=len(rates),exact_saved_flux_rate_sum=str(sum(exact_rates,F())),
        binary_longdouble_summed_rate_erg_s=float(rates.sum(dtype=ld)),rows=records,
        all_pilot_point_gates_passed=all(row['passed'] for row in records),
        surface_pressure_work_control_erg_s=float(surface_work),control_velocity_over_c=float(control_velocity),
        surface_control_is_actual_evolution=False,full_node_grid_completed=False,full_GR_evolution=False,
        physical_transport_certified=False,continuous_lapse_or_Hamiltonian_error_certified=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify();print('SHARED HEAT',records,flush=True)


def verify():
    for name,key in [('plan.json','bindings'),('manifest.json','sha256')]:
        for rel,digest in json.loads((OUT/name).read_text())[key].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'symbolic.json').read_text())['passed']
    assert json.loads((OUT/'result.json').read_text())['completed']
    print('PASS shared material heat faces and regular-centre initial reference bindings',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
