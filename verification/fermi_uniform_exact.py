"""Repeat the frozen Fermi protocol with outward decimal serialization."""
import builtins, inspect, json, shutil, sys
from pathlib import Path
import mpmath.libmp.libmpi as mpi
import mpmath.libmp.libmpf as mpf
import fermi_uniform as u
from interval_records import interval_text, control

OLD=u.OUT
u.OUT=u.g.OUT/'fermi-uniform-exact';u.CACHE=u.g.CACHE/'fermi-uniform-exact'
# Reuse the entire frozen calculation. Only interval serialization is replaced.
u.str=lambda value: interval_text(value) if hasattr(value,'_mpi_') else builtins.str(value)


def prepare():
    u.prepare();plan=json.loads((u.OUT/'plan.json').read_text())
    plan['checkpoint']='f896be7'
    plan['serialization']='Exact binary endpoint -> rational -> directed Decimal division. Interval decimal literals enclose the live interval, including when read at a higher precision.'
    plan['preserved_original']=dict(plan_sha256=u.g.c.sha(OLD/'plan.json'),
        certificate_sha256=u.g.c.sha(OLD/'uniform-certificate.json'),
        point_sha256={p.name:u.g.c.sha(p) for p in OLD.glob('point-*.json')})
    for name in ['fermi_uniform_exact.py','interval_records.py']:
        path=u.g.ROOT/'verification'/name;plan['bindings'][path.relative_to(u.g.ROOT).as_posix()]=u.g.c.sha(path)
    plan['native_point_argument_convention']='Reference eta/beta are the exact declared decimal values. The native call receives their binary64 conversions; the reported discrepancy includes that representation effect. No native uniform-in-argument error is inferred.'
    u.save('plan.json',plan)
    u.save('serialization-control.json',control())
    (u.OUT/'mpmath-serialization-source.txt').write_text(inspect.getsource(mpi.mpi_to_str)+'\n'+inspect.getsource(mpf.to_str))


if __name__=='__main__':
    name=sys.argv[1]
    prepare() if name=='prepare' else getattr(u,name)()
