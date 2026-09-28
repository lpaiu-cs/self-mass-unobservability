"""Independent exact-rational replay of exported electron thermodynamics."""
from fractions import Fraction as F
import json, sys
import electron_thermo_certificate as e


class Box:
    """Only the four rational interval operations needed by this audit."""
    def __init__(self,a,b=None):
        self.a=F(a);self.b=F(a if b is None else b);assert self.a<=self.b
    def __add__(self,other):
        y=other if isinstance(other,Box) else Box(other)
        return Box(self.a+y.a,self.b+y.b)
    __radd__=__add__
    def __neg__(self): return Box(-self.b,-self.a)
    def __sub__(self,other): return self+-(other if isinstance(other,Box) else Box(other))
    def __mul__(self,other):
        y=other if isinstance(other,Box) else Box(other)
        ends=[a*b for a in [self.a,self.b] for b in [y.a,y.b]]
        return Box(min(ends),max(ends))
    __rmul__=__mul__
    def __truediv__(self,other):
        y=other if isinstance(other,Box) else Box(other)
        assert y.a*y.b>0,'Division by an interval containing zero'
        return self*Box(1/y.b,1/y.a)


def read(path): return json.loads(path.read_text())
def box(text): return Box(*e.audit.interval(text))


def run():
    e.verify();e.audit.verify()
    result=read(e.OUT/'result.json');plan=read(e.OUT/'plan.json')
    for rel,digest in plan['bindings'].items(): assert e.g.c.sha(e.g.ROOT/rel)==digest,rel
    cert=read(e.audit.OUT/'uniform-certificate.json');errors={}
    for row in cert['rows']:
        errors[(row['nu'],row['eta_order'],row['beta_order'])]=F(row['panels'],8192)**8*box(row['uniform_error_upper']).b+box(row['tail_upper']).b
    coefficients=[('D','0.5','1.5',F(1),F(1)),('K','1.5','2.5',F(2,3),F(1,3)),('H','1.5','2.5',F(1),F(1))]
    recorded={(r['quantity'],r['eta_order'],r['beta_order']):box(r['error_upper']).b for r in result['uniform_error_components']}
    for name,first,second,a,b in coefficients:
        for k,l in e.u.ORDERS:
            value=a*errors[(first,k,l)]+b*F('0.006')*errors[(second,k,l)]
            if l: value+=l*b*errors[(second,k,l-1)]
            assert recorded[(name,k,l)]>=value
    slopes=result['monotonicity_slabs'];assert [r['eta'] for r in slopes]==[[a,a+1] for a in range(-17,24)]
    m=min(box(r['density_slope_lower']).a for r in slopes);inverse=result['inverse_transfer']
    E0=recorded[('D',0,0)];E1=recorded[('D',1,0)]
    # Recorded common error upper bounds may include outward decimal padding.
    assert F(inverse['density_slope_lower_rational'])==m and m>E1>0
    eta_bound=F(inverse['eta_error_from_zero_quadrature_residual_upper_rational'])
    # Replay from the exact primitive bounds, without the later display padding.
    original_E0=errors[('0.5',0,0)]+F('0.006')*errors[('1.5',0,0)]
    original_E1=errors[('0.5',1,0)]+F('0.006')*errors[('1.5',1,0)]
    assert eta_bound>=original_E0/m
    assert F(inverse['quadrature_density_slope_lower_rational'])<=m-original_E1
    assert F(inverse['reciprocal_slope_error_at_same_eta_upper_rational'])>=original_E1/(m*(m-original_E1))
    checks=[];count=0
    for i,control in enumerate(result['point_controls']):
        point=read(e.audit.OUT/f'point-{i}.json');beta=F(str(point['point'][1]));eta=F(point['point'][0]);fields={}
        for r in point['rows']: fields[(r['nu'],r['eta_order'],r['beta_order'])]=box(r['enclosure'])
        values={}
        for name,first,second,a,b in coefficients:
            for k,l in e.u.ORDERS:
                value=a*fields[(first,k,l)]+b*beta*fields[(second,k,l)]
                if l: value+=l*b*fields[(second,k,l-1)]
                values[(name,k,l)]=value
        D,De,Db=[values[('D',*order)] for order in [(0,0),(1,0),(0,1)]]
        K,Kb=[values[('K',*order)] for order in [(0,0),(0,1)]]
        H,He,Hb=[values[('H',*order)] for order in [(0,0),(1,0),(0,1)]]
        A=(F(3,2)*D+beta*Db)/De
        quantities=dict(reduced_density=D,reduced_pressure=K,reduced_energy=H,
            energy_per_particle_over_kT=H/D,entropy_per_particle_over_k=(H+K)/D-eta,
            deta_dlnnumber_at_T=D/De,deta_dlnT_at_number=-A,
            electron_pressure_chi_number=D*D/(K*De),
            electron_pressure_chi_temperature=Box(F(5,2))+beta*Kb/K-D*A/K,
            electron_cv_per_particle_over_k=(F(5,2)*H+beta*Hb-He*A)/D)
        for key,value in quantities.items():
            exported=box(control['thermodynamics'][key])
            assert exported.a<=value.a and exported.b>=value.b,(i,key)
            count+=1
        native={(r['nu'],r['eta_order'],r['beta_order']):F.from_float(r['native_value']) for r in point['rows']}
        target=native[('0.5',0,0)]+beta*native[('1.5',0,0)]
        assert target==F(control['target_reduced_density_rational'])
        residual=max(abs(D.a-target),abs(D.b-target))
        stored=F(control['residual_upper_rational']);assert stored>=residual
        local_m=F(control['local_density_slope_lower_rational']);delta=F(control['eta_radius_upper_rational'])
        assert local_m>0 and delta==stored/local_m<plan['inverse_control_half_width']
        a,b=map(F,control['eta_enclosure_rational']);assert a==eta-delta and b==eta+delta
        assert control['inverse_exists_and_is_unique'] and quantities['electron_cv_per_particle_over_k'].a>0
        checks.append(dict(point=i,exact_residual_upper_rational=str(residual),
            exported_residual_covers_original_primitive_box=True,inverse_bracket_passed=True))
    path=e.OUT/'audit.json';assert not path.exists()
    e.save('audit.json',dict(classification='Proven',passed=True,exact_rational_thermodynamic_enclosures=count,
        exact_uniform_error_replays=27,inverse_controls=checks,
        point_scope='Thirty derived thermodynamic intervals independently enclose exact rational operations on the original primitive boxes. Root residuals are audited directly against those primitive boxes; widened display intervals are not mistaken for new input uncertainty.',
        physical_EOS_certified=False))
    paths=[e.OUT/name for name in ['plan.json','result.json','symbolic.json','manifest.json','audit.json']]
    paths += [e.g.ROOT/'verification/verify_electron_thermo.py']
    e.save('audit-manifest.json',dict(classification='Proven',sha256={p.relative_to(e.g.ROOT).as_posix():e.g.c.sha(p) for p in paths}))
    print('PASS EXACT ELECTRON AUDIT',count,'thermodynamic intervals, 27 error transfers, 3 unique inverse roots',flush=True)


def verify():
    for rel,digest in read(e.OUT/'audit-manifest.json')['sha256'].items(): assert e.g.c.sha(e.g.ROOT/rel)==digest,rel
    e.verify();e.audit.verify();print('PASS ELECTRON AUDIT BINDINGS',flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
