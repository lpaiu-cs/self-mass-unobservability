"""Conditional Vaidya readout and an exact obstruction at the saved closed boundary."""
import json,sys,urllib.request
import numpy as np
import sympy as sp
from mpmath import iv
import gr_heat_nonlinear_balance_runner as balance
import gr_shared_heat_reference as shared
import gr_outer_product_pilot as intervals

g=shared.g;OUT=g.OUT/'gr-radiating-surface-connection'
URL='https://arxiv.org/pdf/1301.1417'


def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2)+'\n')


def run():
    assert not OUT.exists();balance.verify();shared.verify();OUT.mkdir()
    with urllib.request.urlopen(URL,timeout=60) as response:(OUT/'maharaj-govender-govender2013.pdf').write_bytes(response.read())
    assert (OUT/'maharaj-govender-govender2013.pdf').read_bytes().startswith(b'%PDF')
    files=[g.ROOT/'verification/gr_radiating_surface_connection.py',balance.OUT/'manifest.json',balance.OUT/'result.json',
        shared.OUT/'manifest.json',shared.OUT/'shared-faces.npz',g.OUT/'gr-increment-structure/path-4.npz',OUT/'maharaj-govender-govender2013.pdf']
    save('plan.json',dict(classification='Proven',checkpoint='4bfd9c95',bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in files},
        assumptions=['G=c=1; spherical, smooth timelike material surface outside trapped regions',
            'Interior baryon-frame isotropic P and proper radial heat Q; exterior pure outgoing null dust, not a general material atmosphere',
            'No surface layer: Darmois matching imposes continuous induced metric/extrinsic curvature and the corresponding normal stress/energy flux projections',
            'For the luminosity relation, interior Misner-Sharp mass equals the asymptotically flat outgoing Vaidya mass M(u) on the matched surface'],
        source=dict(classification='Imported from prior work',title='Radiating stars with generalised Vaidya atmospheres',
            authors=['S. D. Maharaj','G. Govender','M. Govender'],url=URL,arxiv='1301.1417',
            inspected_pdf_pages_zero_based=[2,3,4,5,6],equations=['1','2','6','12b','20','21'],
            scope='Their shear-free coordinates give p=qB for the ordinary Vaidya limit. Here Q is the proper comoving heat flux, so qB=Q. The frame-projection necessity and moving-surface readout are rederived below; no shear-free assumption is added to the interior balance theorem.'),
        boundary='A necessary junction condition and conditional bolometric readout, not a matched atmosphere, a finite GR trajectory, data inference, or a new dynamic-chi observable.'))
    eps,P,Q,psi=sp.symbols('epsilon P Q psi',real=True);eta=sp.diag(-1,1);u=sp.Matrix([1,0]);n=sp.Matrix([0,1]);k=u+n
    interior=sp.Matrix([[eps,-Q],[-Q,P]]);kc=eta*k;exterior=psi*kc*kc.T
    normal=(n.T*(interior-exterior)*n)[0];flux=(u.T*(interior-exterior)*n)[0]
    assert normal==P-psi and flux==-Q+psi
    assert sp.solve([normal,flux],[P,psi],dict=True)==[{P:Q,psi:Q}]
    Gamma,U,R=sp.symbols('Gamma U R',positive=True);Ufree=sp.symbols('Ufree',real=True)
    f=Gamma**2-Ufree**2;udot=1/(Gamma+Ufree)
    assert sp.cancel(f*udot**2+2*Ufree*udot-1)==0
    # Gamma>|U| is a declared untrapped condition, so this is the future root.
    a,W,v=sp.symbols('a W v',real=True)
    gamma=W/a;ur=W*v/a
    assert sp.cancel((gamma**2-ur**2).subs(W**2,1/(1-v*v))-1/a**2)==0
    proper_mass_rate=-4*sp.pi*R*R*Q*(Gamma+Ufree)
    luminosity=-proper_mass_rate/udot
    assert sp.cancel(luminosity-4*sp.pi*R*R*Q*(Gamma+Ufree)**2)==0
    ratio=((gamma+ur)**2).subs(W**2,1/(1-v*v))
    assert sp.cancel(ratio-(1+v)/(a*a*(1-v)))==0
    assert sp.simplify(ratio.subs(v,0)-1/a**2)==0
    # The stored grid is ordered from surface (0) to centre (-1).
    faces=dict(np.load(shared.OUT/'shared-faces.npz'));grid=dict(np.load(g.OUT/'gr-increment-structure/path-4.npz'))
    assert faces['radius_cm'][0]>faces['radius_cm'][-1]==0
    assert faces['Linfinity_erg_s'][0]==faces['Q_erg_cm3'][0]==0
    logP=grid['faces'][0,2];assert np.isfinite(logP);iv.prec=256
    log_exact=intervals.F(*logP.as_integer_ratio());pressure=iv.exp(intervals.I(log_exact));assert intervals.low(pressure)>0
    save('result.json',dict(classification='Proven',passed=True,
        necessary_junction='In the comoving orthonormal surface frame, T_in(n,n)=P and T_in(u,n)=-Q. Pure outgoing null dust gives psi and -psi. Continuity therefore requires P=Q=psi. This is necessary, not a claim that this single equation suffices for complete matching.',
        outgoing_metric='ds^2=-(1-2M(u)/R)du^2-2 du dR+R^2 dOmega^2; u is retarded Bondi time at infinity.',
        proper_time_relation='Let sigma be surface proper time, U=dR/dsigma and Gamma=sqrt(U^2+1-2M/R)>|U|. Timelike normalization gives du/dsigma=1/(Gamma+U). On the interior material worldtube U=W*v/a and Gamma=W/a.',
        mass_and_readout='With P=Q and full mass matching, dM/dsigma=-4*pi*R^2*Q*(Gamma+U). Hence L_infinity=-dM/du=4*pi*R^2*Q*(Gamma+U)^2. This is bolometric mass/energy flux; it does not determine a spectrum or instrument likelihood.',
        redshift_factor='L_infinity/L_comoving=(Gamma+U)^2=(1-2M/R)*(1+v)/(1-v). At v=0 this reduces to 1/a^2, equal to N^2 only with the appropriate static exterior lapse normalization.',
        cgs_restoration='For physical local heat flux q in erg/(cm^2 s), L_comoving=4*pi*R_cm^2*q. Use the same dimensionless redshift/Doppler factor. If M_phys is in grams and u in physical seconds, L_infinity=-c^2*dM_phys/du. Q=q/c in the interior stress tensor; it is not q itself.',
        saved_boundary=dict(surface_index=0,log_pressure_exact=str(log_exact),pressure_cgs_enclosure=intervals.cusp.interval_text(pressure),
            pressure_cgs_display=float(pressure.mid),proper_Q_erg_cm3_exact='0',outward_Linfinity_erg_s_exact='0',
            normalized_P_minus_Q_over_P='1',necessary_Vaidya_junction_satisfied=False,
            verdict='The current positive-pressure, zero-heat material endpoint cannot be identified with a no-layer vacuum or pure outgoing Vaidya surface. It is a closed interior boundary. This does not invalidate its earlier closed-boundary control calculations.'),
        missing_boundary='Retain or solve the exterior pressure-bearing atmosphere/radiation stress, or explicitly add a physically justified surface layer. Simply turning on a flux or relabeling the photospheric endpoint does not solve the matching problem. A general atmosphere need not obey the pure-null-dust condition P=Q.',
        actual_atmosphere_solved=False,actual_time_trajectory_computed=False,observational_inference_completed=False,new_dynamic_chi_observable_proven=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    for name,key in [('plan.json','bindings'),('manifest.json','sha256')]:
        for rel,digest in json.loads((OUT/name).read_text())[key].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());assert r['passed'] and not r['saved_boundary']['necessary_Vaidya_junction_satisfied']
    assert r['saved_boundary']['normalized_P_minus_Q_over_P']=='1' and not r['observational_inference_completed']
    print('PASS conditional radiating-surface readout; saved closed surface fails the necessary pure-Vaidya junction with exact relative defect 1',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
