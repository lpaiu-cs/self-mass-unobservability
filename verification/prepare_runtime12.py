"""Create isolated source/run copies and build the exponential-drive extension in WSL."""
import difflib
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

root = Path(__file__).resolve().parents[1]
runtime = Path.home()/'work/nutimo_pilot'
source = runtime/'nutimo_sepdyn/src'
target = runtime/'nutimo_request12/src'
assert not target.exists(), 'Do not replace a prior Request12 build'
shutil.copytree(source, target)
shutil.copytree(runtime/'run_planetGR', runtime/'run_request12')
p = target/'AllTheories3Bodies.cpp'
old = p.read_text()
changes = {
    'value_type SEPdyn_A, SEPdyn_w, SEPdyn_ph;': 'value_type SEPdyn_A, SEPdyn_w, SEPdyn_ph, SEPdyn_tau;',
    'SEPdyn_ph = (e_ = getenv("SEPDYN_PH")) ? atof(e_) : 0.0;':
    'SEPdyn_ph = (e_ = getenv("SEPDYN_PH")) ? atof(e_) : 0.0;\n      SEPdyn_tau = (e_ = getenv("SEPDYN_TAU")) ? atof(e_) : 0.0;',
    'SEPdyn_A * cos( SEPdyn_w * t + SEPdyn_ph )':
    '(SEPdyn_A == zero ? zero : SEPdyn_A * (SEPdyn_tau > zero ? exp(-t/SEPdyn_tau) : cos( SEPdyn_w * t + SEPdyn_ph )))',
}
new = old
for a, b in changes.items():
    assert new.count(a) == 1, a
    new = new.replace(a, b)
p.write_text(new)
out = root/'outputs/research-completion/runtime12'
out.mkdir(parents=True, exist_ok=True)
(out/'transient.patch').write_text(''.join(difflib.unified_diff(old.splitlines(True), new.splitlines(True), fromfile='archived/AllTheories3Bodies.cpp', tofile='request12/AllTheories3Bodies.cpp')))
third = runtime/'install/third_party'
sources = 'AllTheories3Bodies.cpp Delay_brut.cpp Fittriple-compute.cpp Fittriple-init.cpp Fittriple-IO.cpp Fittriple-diagnostics.cpp IO.cpp Orbital_elements.cpp Spline.cpp Utilities.cpp Diagnostics.cpp Parameters.cpp'.split()
command = ['g++-9', '-shared', '-fPIC', '-O2', '-fopenmp', '-std=gnu++11',
           '-I'+str(third/'boost_1_55_0'), '-I'+str(third/'tempo2/include'),
           '-I'+str(third/'Minuit2-5.34.14/include'), *sources,
           '-L'+str(third/'libstatictempo2/with_fpic'), '-ltempo2', '-lsofa',
           '-llapack', '-lcblas', '-lgfortran', '-lf77blas', '-latlas',
           str(third/'Minuit2-5.34.14/lib/libMinuit2.a'), '-o', 'libFittriplecpp.so']
with (out/'build.log').open('w') as log:
    result = subprocess.run(command, cwd=target, stdout=log, stderr=subprocess.STDOUT)
manifest = dict(command=command, returncode=result.returncode,
                source_before_sha256=hashlib.sha256(old.encode()).hexdigest(),
                source_after_sha256=hashlib.sha256(new.encode()).hexdigest())
(out/'build.json').write_text(json.dumps(manifest, indent=2)+'\n')
assert result.returncode == 0, 'See build.log'
print('Built isolated runtime', target)
