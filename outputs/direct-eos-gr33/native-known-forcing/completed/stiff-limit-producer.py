"""Proven: scalar stiff-limit boundary of the attempted source quadratures."""
from pathlib import Path
import sympy as s
import integrate_native_known_forcing as r
z,t=s.symbols('z t',real=True)
A=s.Matrix([[s.Rational(5,12),-s.Rational(1,12)],[s.Rational(3,4),s.Rational(1,4)]])
c=[s.Rational(1,3),s.Integer(1)];S=s.Matrix([v**3 for v in c]);H=s.Matrix([v**4/4 for v in c])
J=s.Matrix([v**5/20 for v in c]);Heff=A.inv()*J
raw=(s.eye(2)-z*A).inv()*A*S
lifted=(s.eye(2)-z*A).inv()*H
fitted=(s.eye(2)-z*A).inv()*(z*A*Heff)+H
limit=lambda v:v.applyfunc(lambda x:s.limit(x,z,-s.oo))
assert limit(z*raw)==-S
assert limit(z*lifted)==-A.inv()*H
assert limit(fitted)==H-Heff and H-Heff!=s.zeros(2,1)
result=dict(classification='Proven',passed=True,equation='xprime=lambda*x+t^3, step h=1, x(0)=0, lambda->-infinity',
    direct_Radau_scaled_limit=[str(v) for v in limit(z*raw)],
    exact_primitive_scaled_limit=[str(v) for v in limit(z*lifted)],
    fitted_known_stream_unscaled_limit=[str(v) for v in limit(fitted)],
    interpretation='Direct Radau preserves lambda*x_j=-S_j at leading stiff order. Exact primitive stage forcing gives -A_inverse*H instead. Moment-fitting L*H creates a nonzero H-Heff limit, although the actual stiff solution tends to zero.',
    scope='Exact scalar counterexample to a generic stiff-accuracy claim for these quadratures; it is not proof of the unique cause in the full time-dependent transport system.',
    bindings={str(Path(__file__)):r.sha(__file__)})
r.write(r.OUT/'stiff-limit.json',result);print(r.original.reuse.json_text(result),flush=True)
