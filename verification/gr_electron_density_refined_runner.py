"""Give the isolated frozen evaluator an importable process-pool entry point."""
import json,shutil,sys
import gr_electron_density_refinement as original

g=original.g;OUT=original.OUT;raw_evaluate=original.ns['evaluate']


def evaluate(first,count,label):return raw_evaluate(first,count,label)
original.ns['evaluate']=evaluate


def prepare():
    assert not list(OUT.glob('block-*.jsonl'))
    assert not (OUT/'worker-binding.json').exists()
    log=g.ROOT/'outputs/gr-electron-density-refined33-run.log'
    assert "import of module '<run_path>' failed" in log.read_text()
    shutil.copy2(log,OUT/'worker-dispatch-failure.log')
    files=[g.ROOT/'verification/gr_electron_density_refined_runner.py',OUT/'worker-dispatch-failure.log']
    original.ns['save']('worker-binding.json',dict(classification='Proven',
        reason='The frozen runpy evaluator was not importable by ProcessPoolExecutor. A module-level wrapper dispatches exactly the same evaluator and arguments; zero numerical blocks existed at the failure.',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in files},
        algorithm_or_input_changed=False))


def worker_bindings():
    for rel,digest in json.loads((OUT/'worker-binding.json').read_text())['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel


def run():worker_bindings();original.run()
def verify():worker_bindings();original.verify()
export_namespace=original.export_namespace
def export_run():verify();original.export_run()
def export_verify():verify();original.export_verify()
if __name__=='__main__':globals()[sys.argv[1]]()
