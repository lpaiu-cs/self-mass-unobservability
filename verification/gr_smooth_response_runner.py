"""Resolve CAPD's interval abs namespace; retain the original pre-run build failure."""
import inspect,json,shutil,sys
import gr_smooth_response_certificate as g

ORIGINAL=g.OUT;OUT=ORIGINAL.parent/'gr-smooth-response-defined';CACHE=g.native.g.CACHE/'smooth-response-defined'


def prepare():
    p=g.bindings();assert not OUT.exists() and not CACHE.exists();OUT.mkdir();CACHE.mkdir()
    error=(ORIGINAL/'build.log').read_text();assert 'error: no matching function for call to' in error and not list(ORIGINAL.glob('block-*.jsonl*'))
    source=(g.ROOT/'verification/gr_smooth_response_certificate.cpp').read_text();anchor='using Six=std::array<MI,6>;';assert source.count(anchor)==1
    (OUT/'effective-source.cpp').write_text(source.replace(anchor,'using capd::abs;\n'+anchor))
    for name in ['inputs.tsv','rule.json','rule.tsv','rule-manifest.json']:shutil.copy2(ORIGINAL/name,OUT/name)
    shutil.copy2(ORIGINAL/'build.log',OUT/'original-build-error.log')
    p['bindings'].update({path.relative_to(g.ROOT).as_posix():g.sha(path) for path in [g.ROOT/'verification/gr_smooth_response_runner.py',
        ORIGINAL/'plan.json',ORIGINAL/'rule-manifest.json',OUT/'inputs.tsv',OUT/'effective-source.cpp',OUT/'original-build-error.log']})
    p['execution_overlay']='Only import capd::abs for the interval absolute-value calls. Original compile failure preceded all numerical cases. Quadrature rule, formulas, input states, exact partition checks, error budgets and controls are unchanged.'
    (OUT/'plan.json').write_text(json.dumps(p,ensure_ascii=False,indent=2)+'\n')


def activate():
    g.OUT=OUT;g.CACHE=CACHE
    source=inspect.getsource(g.build);old="ROOT/'verification/gr_smooth_response_certificate.cpp'";assert source.count(old)==1
    exec(compile(source.replace(old,"OUT/'effective-source.cpp'"),str(OUT/'build-overlay.py'),'exec'),g.__dict__)


def build():activate();g.build()
def run():activate();g.run()
def controls():activate();g.controls()
def finalize():activate();g.finalize()
def verify():activate();g.verify()


if __name__=='__main__':globals()[sys.argv[1]]()
