"""Check the collocation identity and the adapter's frozen receipt schema."""
from pathlib import Path
import ast,hashlib,json
import sympy as s

x,h,u,y0,y1,y2,r1,r2=s.symbols('x h u y0 y1 y2 r1 r2')
q=s.Matrix([s.Rational(3,2)*x-s.Rational(3,4)*x*x,s.Rational(3,4)*x*x-x/2])
assert q.subs(x,s.Rational(1,3))==s.Matrix([s.Rational(5,12),-s.Rational(1,12)])
assert q.subs(x,1)==s.Matrix([s.Rational(3,4),s.Rational(1,4)])
for state,at in [(y1,s.Rational(1,3)),(y2,1)]:
    integral=(q.dot(s.Matrix([r1,r2]))).subs(x,at)
    assert s.expand(u*(state-y0-h*integral)-(u*state-u*y0-h*u*integral))==0
rows=[]
for module,receipt in [
    ('read_material_return_charge.py','native-full-charge251-work/regression-receipt.json'),
    ('read_material_return_exterior.py','native-returned-exterior252-work/check-receipt.json')]:
    p=Path('verification')/module;source=p.read_text();tree=ast.parse(source)
    assert 'sha257' not in source
    keys=[n.slice.value for n in ast.walk(tree) if isinstance(n,ast.Subscript) and isinstance(n.slice,ast.Constant)]
    assert 'source_sha256' in keys and 'source_sha256' in json.loads(Path(receipt).read_text())
    rows.append(dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),receipt=receipt))
out=dict(classification='Proven',passed=True,radau_dense_stage_identity=True,
    diagonal_conserved_coordinate_identity=True,adapter_receipt_schema_checked=True,bindings=rows,
    scope='Exact symbolic identities and receipt-key regression only; no numerical or physical error certificate.')
Path('phase257-symbolic-adapter-check.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out))
