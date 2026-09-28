import json
from pathlib import Path
import sympy as s
h,a,a0,c,D,K,e,v,p,r,u,F,dF=s.symbols('h a a0 c D K e v p r u F dF')
U=a*K+(a-a0)*c*D
changed=a*(1+h*e)*K+(a*(1+h*e)-a0)*c*D
assert s.expand(s.diff(changed,h).subs(h,0)-e*(U+a0*c*D))==0
assert s.expand(s.diff(v*(changed+a*(1+h*e)*p),h).subs(h,0)-v*e*(U+a0*c*D+a*p))==0
assert s.expand(s.diff((1+h*e)*(1+h*u)**2*(F+h*dF),h).subs(h,0)-(dF+(e+2*u)*F))==0
Path('native-full-incident181-work/geometric-face-symbolic.json').write_text(json.dumps(dict(
    classification='Proven',passed=True,
    scope='At fixed primitive state, dU_energy=ell*(U_energy+a_surface*cx*D), dF_energy=v*(dU_energy+a*ell*p); lapse times area contributes(ell+2*u)*F. Exact finite-owner product rules, not native EOS or full-GR certification.'),indent=2)+'\n')
