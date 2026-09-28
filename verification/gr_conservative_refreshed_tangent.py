"""Counterexample candidate: refresh the conservative native tangent each iteration.

Same equations, BDF history, residual and conservation gates as the frozen
parent. Replay its four accepted coarse states and warm-start step five from
its preserved failed iterate. This does not overwrite the original failure.
"""
from types import FunctionType
import gr_conservative_two_carrier as prior
from gr_conservative_two_carrier import *

OUT = e.g.OUT/'gr-conservative-refreshed-tangent'
REFINEMENT = None


def initialize(pool):
    star = prior.initialize(pool)
    star.operator_step = 0
    star.operator_time = ld(0)
    return star


def stage(star, previous, older, h, coefficients, log):
    star.operator_step += 1
    star.operator_time += h
    step = star.operator_step
    delta, history, factor = previous[0].copy(), [], None
    if REFINEMENT == 1 and step <= 5:
        path = prior.OUT/'path-1'/('failed-iterate.npz' if step == 5 else f'step-{step:04d}.npz')
        cp = np.load(path)
        assert cp['time_seconds'] == star.operator_time
        delta = cp['delta'].copy()
        y = star.base+delta
        star.material_cache = {e.material_key(row):aux for row,aux in zip(zip(y[:,0],y[:,1],y[:,5:]),cp['aux'])}
        value,z = residual(star,delta,previous,older,h,coefficients)
        if step <= 4:
            for key in ['m','mf','a','N','Q','aux','dU']:
                assert np.array_equal(z[key],cp[key]),(step,key)
            norm = float(np.max(abs(value)/ATOL))
            assert norm <= 1
            log(dict(iteration=0,residual_norm=norm,maximum_absolute_residual=np.max(abs(value),axis=0).astype(float).tolist(),
                     imported_accepted_native_state=str(path.relative_to(e.ROOT))))
            return delta,z
    value,z = residual(star,delta,previous,older,h,coefficients)
    for iteration in range(24):
        star.last_delta, star.last_state = delta.copy(), z
        norm = float(np.max(abs(value)/ATOL))
        log(dict(iteration=iteration,residual_norm=norm,maximum_absolute_residual=np.max(abs(value),axis=0).astype(float).tolist()))
        if norm <= 1:
            return delta,z
        if iteration == 23:
            break
        # Refresh at the actual iterate; the prior step-five stale directions
        # increased the native residual for every tested line-search fraction.
        factor = splu(jacobian(tangent(star,delta,z),delta,previous,older,h,coefficients))
        x,f = (delta[:, :5]/SCALE).ravel(),(value[:, :5]/SCALE).ravel()
        mixed,raw = parent.parent.direction(x,f,factor,history)
        accepted = None
        trials = []
        for proposal in ([mixed,raw] if history else [raw]):
            correction = proposal.reshape(star.n,5).astype(ld)*SCALE
            fraction = min(1.,.1/max(float(abs(correction[:, :2]).max()),1e-300),.01/max(float(abs(correction[:, 2]).max()),1e-300))
            for backtrack in range(8):
                candidate = delta.copy()
                candidate[:, :5] += fraction*correction
                try:
                    candidate[:, 5:] = species(star,candidate,previous,older,h,coefficients,z)
                    trial,state = residual(star,candidate,previous,older,h,coefficients)
                    score = float(np.max(abs(trial)/ATOL))
                    trials.append(dict(fraction=fraction,norm=score))
                    if score <= 1 or score < norm*(1-1e-4*fraction):
                        accepted = candidate,trial,state
                        break
                except (AssertionError,np.linalg.LinAlgError,FloatingPointError) as error:
                    trials.append(dict(fraction=fraction,invalid_trial=repr(error)))
                fraction /= 2
            if accepted is not None:
                break
        if accepted is None:
            raise RuntimeError(('Conservative native line search failed',norm,trials))
        history = (history+[(x.copy(),f.copy())])[-6:]
        delta,value,z = accepted
    raise RuntimeError(('Conservative native BDF stage did not converge',norm))


def prepare():
    assert not OUT.exists()
    plan = prior.bindings()
    assert prior.symbolic()['passed'] and plan['implementation_check']['passed']
    folder = prior.OUT/'path-1'
    failure = json.loads((folder/'failure.json').read_text())
    assert failure['step'] == 5 and 'line search failed' in failure['reason']
    files = [folder/f'step-{step:04d}.npz' for step in range(5)]
    files += [folder/'failed-iterate.npz',folder/'failure.json',prior.OUT/'plan.json',Path(__file__)]
    for path in files:
        plan['bindings'][path.relative_to(e.ROOT).as_posix()] = e.digest(path)
    plan.update(checkpoint=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        candidate='Same native conservative two-carrier equations with a refreshed tangent at each nonlinear iterate.',
        nonlinear_operator='Rebuild the five-field approximate native tangent every nonlinear iteration, retain the same Anderson six-vector history, 24-iteration limit, eight backtracks, native acceptance tolerances, time grid, physics and conservation gates.',
        imported_prefix='Only coarse path 1 replays parent accepted steps 0-4 with native residual and bitwise state/metric checks. Step 5 starts from the preserved failed nonlinear iterate. Paths 2 and 4 start from the original initial state without imported iterates. A successful new path never changes the parent failure.',
        parent_plan_sha256=e.digest(prior.OUT/'plan.json'))
    OUT.mkdir()
    e.write(OUT/'plan.json',plan)
    print('FROZEN REFRESHED TANGENT; UNCHANGED CONSERVATIVE EQUATIONS AND GATES',flush=True)


bindings = FunctionType(prior.bindings.__code__,globals())
completed = FunctionType(prior.completed.__code__,globals())
compare = FunctionType(prior.compare.__code__,globals())
_native_run = FunctionType(prior.run.__code__,globals())


def run(refinement,workers):
    global REFINEMENT
    REFINEMENT = refinement
    return _native_run(refinement,workers)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command',choices=['prepare','run','completed','compare'])
    parser.add_argument('--refinement',type=int,choices=[1,2,4],default=1)
    parser.add_argument('--workers',type=int,default=15)
    args = parser.parse_args()
    if args.command == 'run':
        run(args.refinement,args.workers)
    elif args.command == 'completed':
        print(json.dumps(completed(args.refinement)[1],indent=2))
    else:
        globals()[args.command]()
