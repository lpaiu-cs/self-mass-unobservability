"""Differentiate the subcell moment map and audit saved nonlinear recoveries."""
import json,sys
from fractions import Fraction as F
import numpy as np
import sympy as s
import gr_subcell_inverse as inverse

g=inverse.g;OUT=g.OUT/'gr-subcell-audit'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def run():
    assert not OUT.exists();OUT.mkdir();inverse.verify()
    paths=[g.ROOT/'verification/gr_subcell_audit.py',inverse.OUT/'manifest.json',inverse.reference.OUT/'manifest.json']
    save('plan.json',dict(classification='Proven',checkpoint='6338252',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        scope='Direct node-map derivative, positive frozen cell factors, finite quadrature and stored nonlinear recovery checks. No fresh native or continuum certificate.'))
    eta,theta,v=s.symbols('eta theta v',real=True)
    rho0,C,u0,ur,cvT,P0,Pr,Pt,Q=s.symbols('rho0 C u0 ur cvT P0 Pr Pt Q',real=True)
    rho=rho0*s.exp(eta);u=u0+ur*eta+cvT*theta;P=P0+Pr*eta+Pt*theta;eps=rho*(C+u)
    W=1/s.sqrt(1-v*v);D=rho*W;E=(eps+P*v*v+2*v*Q)*W*W;J=((eps+P)*v+Q*(1+v*v))*W*W
    matrix=s.Matrix([D,E,J]).jacobian([eta,theta,v]).subs({eta:0,theta:0,v:0})
    expected=s.Matrix([[rho0,0,0],[rho0*(C+u0+ur),rho0*cvT,2*Q],[0,0,rho0*(C+u0)+P0]])
    assert all(s.simplify(x)==0 for x in matrix-expected)
    B,V,H,A,Z=s.symbols('B V H A Z',real=True)
    assert s.Matrix([[B,0,0],[A,V,2*Z],[0,0,H]]).det()==B*V*H
    save('symbolic.json',dict(classification='Proven',passed=True,
        construction='Directly differentiate the full Lorentz map with rho=rho0*exp(eta) and arbitrary first-order EOS jet u=u0+u_lnrho*eta+cvT*theta, P=P0+P_lnrho*eta+P_lnT*theta, fixed Q. The initial node Jacobian has rows [rho0,0,0], [rho0*(C+u0+u_lnrho),rho0*cvT,2Q], [0,0,rho0*(C+u0)+P0].',
        integration='Sum row 1 with positive proper-volume weights, row 2 with positive coordinate-volume weights, row 3 with proper-volume weights. Subtract C*kappa times row 1 from row 2. Linearity gives determinant B0*Cv_coordinate*enthalpy_proper>0 for any finite positive quadrature with positive rho,cvT,enthalpy. The continuum version additionally requires differentiation under the integral.',
        qualification='Local inverse theorem still requires a differentiable EOS. This exact frozen-coefficient result does not enclose native evaluation, branch or global inverse errors.'))
    plan=json.loads((inverse.OUT/'plan.json').read_text());report=json.loads((inverse.OUT/'result.json').read_text())
    quads=json.loads((inverse.reference.OUT/'subcell-quadrature.json').read_text())['rows']
    factors=[];recovery=[]
    for i in plan['cells']:
        data=dict(np.load(inverse.reference.OUT/f'cell-{i}-tolerance-1-nodes-16.npz'))
        C=F(float(data['C_X']))*(F(float(g.c.gr.C))*100)**2;B=F(0);Cv=F(0);H=F(0)
        for a,cw,pw in zip(data['eos'],data['coordinate_weights_cm3'],data['proper_weights_cm3']):
            rho,P,u,cv=map(lambda x:F(float(x)),[a[0],a[1],a[2],a[10]])
            cw,pw=F(float(cw)),F(float(pw));assert rho>0 and cw>0 and pw>0 and cv>0
            enthalpy=rho*(C+u)+P;assert enthalpy>0
            B+=rho*pw;Cv+=rho*cv*cw;H+=enthalpy*pw
        assert B>0 and Cv>0 and H>0
        factors.append(dict(cell=i,exact_positive_factors=True,B=str(B),Cv_coordinate=str(Cv),enthalpy_proper=str(H)))
        for case in range(3):
            d=dict(np.load(inverse.OUT/f'cell-{i}-case-{case}.npz'))
            errors=abs(d['true_parameters']-d['recovered_parameters'])
            assert errors[0]<=plan['root_log_density_tolerance'] and errors[1]<=plan['root_logT_tolerance'] and errors[2]<=plan['root_velocity_tolerance']
            assert max(abs(d['residual']))<=plan['scaled_residual_tolerance']
            recovery.append([*map(float,errors),float(abs(d['recovered_moments'][0]/d['target'][0]-1)),float(max(abs(d['residual'])))])
    assert report['all_passed'] and len(recovery)==12 and all(r['passed'] for r in quads)
    summary=dict(classification='Counterexample candidate',subcell_cases=len(recovery),all_passed=True,
        maximum_errors_eta_theta_v_B_residual=list(map(float,np.max(recovery,axis=0))),
        maximum_quadrature_errors={k:max(abs(r[k]) for r in quads) for k in [
            'baryon_relative_direct_difference','baryon_relative_inventory_difference','mass_relative_direct_difference','internal_energy_relative_direct_difference']},
        reference_EOS_statistics=json.loads((inverse.reference.OUT/'result.json').read_text())['entropy_root_statistics'],
        inverse_EOS_evaluations=report['EOS_evaluations'],
        interpretation='All four nonuniform reference cells, both node counts and DOP tolerances passed the registered finite integral comparisons without normalization. All twelve native moment inverses passed. Fixed geometry, composition and Q/profile family; no full-star trajectory or continuum/physical certification.')
    save('positive-factors.json',dict(classification='Proven',rows=factors));save('summary.json',summary)
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    print('PASS direct subcell Jacobian, exact positive factors and 12 nonlinear moment inverses',flush=True);verify()


def verify():
    for rel,digest in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'symbolic.json').read_text())['passed'] and json.loads((OUT/'summary.json').read_text())['all_passed']
    print('PASS independent subcell audit bindings',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
