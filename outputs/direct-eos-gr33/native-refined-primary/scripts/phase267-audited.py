"""Phase267 diagnostic wrapper: run a script while recording opened files and existence checks.

Existence checks (Path.exists/is_file/is_dir, os.path.exists/isfile/isdir) raise no audit event but can switch
code paths, so both are written at exit. Usage: python3 .phase267-audited.py <record json> <script> [args...]
"""
import atexit, json, os, pathlib, runpy, sys
record = sys.argv[1]; opened = {}; checked = {}
def hook(event, args):
    if event == 'open' and args and isinstance(args[0], (str, bytes, os.PathLike)):
        p = os.path.abspath(os.fsdecode(args[0])); m = args[1] if len(args) > 1 else None
        e = opened.setdefault(p, [0, 0]); e[1 if isinstance(m, str) and any(c in m for c in 'wax+') else 0] += 1
sys.addaudithook(hook)
def recorder(original, is_method):
    def wrapped(*a, **k):
        r = original(*a, **k); checked.setdefault(os.path.abspath(str(a[0]) if is_method else os.fsdecode(a[0])), bool(r)); return r
    return wrapped
for name in ['exists', 'is_file', 'is_dir']: setattr(pathlib.Path, name, recorder(getattr(pathlib.Path, name), True))
for name in ['exists', 'isfile', 'isdir']: setattr(os.path, name, recorder(getattr(os.path, name), False))
atexit.register(lambda: pathlib.Path(record).write_text(json.dumps(dict(opened=opened, checked=checked), indent=1)))
sys.argv = sys.argv[2:]; runpy.run_path(sys.argv[0], run_name='__main__')
