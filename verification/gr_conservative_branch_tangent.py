"""Counterexample candidate: differentiate the current material upwind branch.

The actual native residual, EOS, fluxes, BDF scheme and scientific gates stay
unchanged. Only the iteration tangent keeps its anchor donor during finite
differences. Accepted parent coarse states are replayed with full residuals.
"""
import gr_conservative_refreshed_tangent as previous
from gr_conservative_refreshed_tangent import *

OUT = e.g.OUT/'gr-conservative-branch-tangent'
REFINEMENT, PREFIX = None, 0


class BranchTangent(prior.LocalTangent):
    def finish_moments(self,delta,z,dE):
        prior.ConservativeStar.finish_moments(self,delta,z,dE)
        anchor = self.linearization
        direction = self.faces(anchor['N']*anchor['v']/anchor['a'],odd=True)
        velocity = self.faces(z['N']*z['v']/z['a'],odd=True)
        # Differentiate the selected branch even when the finite perturbation
        # crosses zero velocity. The actual native flux still selects freely.
        fB = velocity*prior.upwind(z['B'],direction)
        fE = velocity*prior.upwind(z['E'],direction)+self.faces(z['N']/z['a']*(z['P']*z['v']+z['Q']),odd=True)
        fS = velocity*prior.upwind(z['AS'],direction)+self.faces(z['N']*(z['P']+z['Q']*z['v']))
        fS[-1] = z['N'][-1]*z['R'][-1]
        z['fluxes'] = (fB,fE,fS)
        rate = -e.GRAV*e.C*4*np.pi*self.rf**2*fE
        z['at'] = z['a']**3*(rate[:-1]+self.fraction*np.diff(rate))/self.r


def tangent(star,delta,z):
    model = prior.tangent(star,delta,z)
    model.__class__ = BranchTangent
    return model


def stage(star, previous_state, older, h, coefficients, log):
    star.operator_step += 1
    star.operator_time += h
    step = star.operator_step
    delta, history, factor = previous_state[0].copy(), [], None
    if REFINEMENT == 1 and step <= PREFIX:
        path = previous.OUT/'path-1'/f'step-{step:04d}.npz'
        cp = np.load(path)
        assert cp['time_seconds'] == star.operator_time
        delta = cp['delta'].copy()
        y = star.base+delta
        star.material_cache = {e.material_key(row):aux for row,aux in zip(zip(y[:,0],y[:,1],y[:,5:]),cp['aux'])}
        value,z = residual(star,delta,previous_state,older,h,coefficients)
        for key in ['m','mf','a','N','Q','aux','dU']:
            assert np.array_equal(z[key],cp[key]),(step,key)
        norm = float(np.max(abs(value)/ATOL))
        assert norm <= 1
        log(dict(iteration=0,residual_norm=norm,maximum_absolute_residual=np.max(abs(value),axis=0).astype(float).tolist(),
                 imported_accepted_native_state=str(path.relative_to(e.ROOT))))
        return delta,z
    value,z = residual(star,delta,previous_state,older,h,coefficients)
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
        factor = splu(jacobian(tangent(star,delta,z),delta,previous_state,older,h,coefficients))
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
                    candidate[:, 5:] = species(star,candidate,previous_state,older,h,coefficients,z)
                    trial,state = residual(star,candidate,previous_state,older,h,coefficients)
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


def check():
    assert prior.symbolic()['passed']
    plan = prior.bindings()
    star = prior.initialize(None)
    times = parent.wall.prior.time_nodes(plan,1)
    states = []
    for name in ['step-0003.npz','step-0004.npz','failed-iterate.npz']:
        cp = np.load(prior.OUT/'path-1'/name)
        delta = cp['delta'].copy()
        y = star.base+delta
        star.material_cache = {e.material_key(row):aux for row,aux in zip(zip(y[:,0],y[:,1],y[:,5:]),cp['aux'])}
        states.append((delta,star.evaluate(delta)))
    older,last,(delta,z) = states
    h = times[5]-times[4]
    coefficients = weights(h,times[4]-times[3])
    value,z = residual(star,delta,last,older,h,coefficients)
    physical = prior.tangent(star,delta,z)
    errors = {}
    for name,model in [('old',physical),('selected_branch',tangent(star,delta,z))]:
        factor = splu(jacobian(model,delta,last,older,h,coefficients))
        change = np.zeros_like(delta)
        change[:,:5] = -factor.solve(np.asarray(value[:,:5]/SCALE,float).ravel()).reshape(star.n,5).astype(ld)*SCALE
        step = ld('.0001')
        plus = residual(physical,delta+step*change,last,older,h,coefficients)[0]
        minus = residual(physical,delta-step*change,last,older,h,coefficients)[0]
        defect = (plus-minus)/(2*step)+value
        errors[name] = (np.max(abs(defect[:,:5]),axis=0)/np.maximum(np.max(abs(value[:,:5]),axis=0),1e-99)).astype(float).tolist()
    assert max(errors['old'][:2]) > 1e-4,errors
    assert max(errors['selected_branch'][:2]) < 1e-6,errors
    assert max(errors['selected_branch']) < 1e-4,errors
    result = dict(classification='Counterexample candidate',passed=True,
        relative_direction_errors=errors,actual_new_native_step=False,
        scope='Actual preserved failure-state tangent regression. No EOS, global nonlinear convergence or physical certification.')
    print('SELECTED BRANCH TANGENT CHECK',json.dumps(result),flush=True)
    return result


def prepare():
    assert not OUT.exists()
    plan = previous.bindings()
    checked = check()
    prefix = json.loads((previous.OUT/'path-1/progress.json').read_text())['step']
    assert prefix >= 5
    files = [previous.OUT/'path-1'/f'step-{step:04d}.npz' for step in range(prefix+1)]
    files += [previous.OUT/'plan.json',Path(__file__)]
    for path in files:
        plan['bindings'][path.relative_to(e.ROOT).as_posix()] = e.digest(path)
    plan.update(checkpoint=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        candidate='Native conservative evolution with a tangent of the selected material upwind branch.',
        nonlinear_operator='Refresh the local five-field tangent every iterate, holding the current donor choice only while differentiating. Actual native residuals choose their own donors. Same Anderson history, 24 iterations, eight backtracks and all original scientific gates. Metric and native EOS composition derivatives remain approximate in this iteration operator.',
        imported_prefix_steps=prefix,
        imported_prefix='Coarse path 1 replays only the bound accepted refreshed-tangent prefix, with bitwise state/metric and actual BDF residual checks; then starts each next step from its last accepted state. No failed-iterate warm start. Refined paths start from the original initial state.',
        parent_plan_sha256=e.digest(previous.OUT/'plan.json'),
        operator_check=checked)
    OUT.mkdir()
    e.write(OUT/'plan.json',plan)
    print('FROZEN SELECTED BRANCH TANGENT; IMPORTED ACCEPTED STEPS',prefix,flush=True)


bindings = FunctionType(prior.bindings.__code__,globals())
completed = FunctionType(prior.completed.__code__,globals())
compare = FunctionType(prior.compare.__code__,globals())
_native_run = FunctionType(prior.run.__code__,globals())


def run(refinement,workers):
    global REFINEMENT,PREFIX
    REFINEMENT = refinement
    PREFIX = bindings()['imported_prefix_steps'] if refinement == 1 else 0
    return _native_run(refinement,workers)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command',choices=['check','prepare','run','completed','compare'])
    parser.add_argument('--refinement',type=int,choices=[1,2,4],default=1)
    parser.add_argument('--workers',type=int,default=15)
    args = parser.parse_args()
    if args.command == 'run':
        run(args.refinement,args.workers)
    elif args.command == 'completed':
        print(json.dumps(completed(args.refinement)[1],indent=2))
    else:
        globals()[args.command]()
