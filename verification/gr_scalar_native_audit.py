"""Preserve the failed bitwise scalar-entry assertion and audit its EOS call boundary."""
import inspect,json,sys
import numpy as np
import gr_scalar_reference as original

g=original.g;OUT=g.OUT/'gr-scalar-native-audit'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def diagnose():
    assert not OUT.exists();OUT.mkdir()
    state=dict(np.load(g.OUT/'initial-state-17-4.npz'));native=np.load(g.OUT/'gr-microphysics/auxiliaries.npz')['eos']
    delta=native[:,2].astype(np.longdouble)-state['u_W'].astype(np.longdouble)
    budget=np.maximum(2.,32*np.spacing(abs(state['u_W'])))
    np.savez_compressed(OUT/'difference.npz',pressure_root_energy=state['u_W'],density_replay_energy=native[:,2],
        difference_erg_g=delta,budget_erg_g=budget,EOS=native)
    record=dict(classification='Counterexample candidate',completed=True,bitwise_internal_energy_equal=False,
        different_states=int(np.count_nonzero(delta)),maximum_absolute_difference_erg_g=float(abs(delta).max()),
        maximum_difference_over_cvT=float(np.max(abs(delta)/native[:,10])),
        maximum_existing_energy_budget_score=float(np.max(abs(delta)/budget)),
        existing_energy_budget_passed=bool(np.all(abs(delta)<=budget)),
        pressure_root_source='common_eos.structure calls the native pressure/entropy inverse and stores its energy. gr_microphysics reevaluates at saved log-density/log-temperature/composition; these are distinct rounded calls.',
        scope='The original bitwise scalar-entry assertion failed before any scalar IVP. This does not establish physical scalar failure or authorize silently replacing the original GR state.')
    save('result.json',record)
    log=g.ROOT/'outputs/gr-scalar-reference33.log'
    original.save('failure.json',dict(classification='Counterexample candidate',completed=False,
        failed_precondition='np.array_equal(native density replay u, stored pressure-root u_W)',
        scalar_IVP_started=False,log_sha256=g.c.sha(log)))
    original.save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p)
        for p in original.OUT.iterdir() if p.is_file() and p.name!='manifest.json'}))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in [
        g.ROOT/'verification/gr_scalar_native_audit.py',original.OUT/'manifest.json',log,
        OUT/'result.json',OUT/'difference.npz']}))
    print('SCALAR ENTRY ENERGY AUDIT',record,flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
