"""Counterexample candidate: include native EOS composition response in the solve.

Two neighboring composition directions provide a local iteration basis. They
do not restrict the actual 26-species state or change the native residual.
The physical EOS, fluxes, BDF histories and acceptance gates stay unchanged.
"""
from gr_conservative_branch_tangent import *
import gr_conservative_branch_tangent as previous

HISTORY = previous.previous
OUT = e.g.OUT/'gr-conservative-composition-tangent'
REFINEMENT,PREFIX = None,0


def composition_response(star,delta,z):
    y = star.base+delta
    X = y[:,5:]
    basis = np.stack([np.vstack([X[:1],X[:-1]])-X,np.vstack([X[1:],X[-1:]])-X],axis=1)
    inverse = np.linalg.pinv(np.asarray(basis.transpose(0,2,1),float),rcond=1e-10)
    keys,jobs = [],[]
    for i,k in np.argwhere(np.any(basis != 0,axis=2)):
        amount = min(ld(1),ld('.00001')/np.max(abs(basis[i,k])))
        neighbor = X[max(0,i-1)] if k == 0 else X[min(star.n-1,i+1)]
        mixture = (1-amount)*X[i]+amount*neighbor
        keys.append((int(i),int(k),amount))
        jobs.append((y[i,0],y[i,1],mixture))
    coefficients = np.zeros((star.n,2,4),dtype=ld)
    native = parent.parent.material
    answers = map(native,jobs) if star.pool is None else star.pool.map(native,jobs,chunksize=16)
    for (i,k,amount),aux in zip(keys,answers):
        for q,j in enumerate([1,2,24,27]):
            difference = aux[j]-z['aux'][i,j] if j == 2 else np.log(aux[j]/z['aux'][i,j])
            coefficients[i,k,q] = difference/amount
    assert np.all(np.isfinite(coefficients))
    print('NATIVE COMPOSITION RESPONSES',len(jobs),flush=True)
    # ponytail: two local directions and a lagged stage response are iteration
    # approximations; use a wider basis if native residual convergence needs it.
    return dict(inverse=inverse,coefficients=coefficients,evaluations=len(jobs))


class CompositionTangent(previous.BranchTangent):
    def evaluate(self,delta):
        last,older,h,coefficients = self.composition_context
        delta = delta.copy()
        projected = species(self,delta,last,older,h,coefficients,self.linearization)
        delta[:,5:] = self.anchor[:,5:]+(projected-self.projected_anchor)
        change = delta-self.anchor
        aux = self.linearization['aux'].copy()
        factors = np.einsum('nki,ni->nk',self.composition_data['inverse'],change[:,5:])
        extra = np.einsum('nk,nkj->nj',factors,self.composition_data['coefficients'])
        aux[:,1] *= np.exp(aux[:,5]*change[:,0]+aux[:,6]*change[:,1]+extra[:,0])
        aux[:,2] += aux[:,9]*change[:,0]+aux[:,10]*change[:,1]+extra[:,1]
        for q,j in [(2,24),(3,27)]:
            aux[:,j] *= np.exp(aux[:,j+1]*change[:,0]+aux[:,j+2]*change[:,1]+extra[:,q])
        aux[:,21] = aux[:,24]*aux[:,27]/(aux[:,24]+aux[:,27])
        y = self.base+delta
        self.material_cache = {e.material_key(row):value for row,value in zip(zip(y[:,0],y[:,1],y[:,5:]),aux)}
        z = prior.ConservativeStar.evaluate(self,delta)
        for key in ['m','mf','a','da','N','nur']:
            z[key] = self.linearization[key]
        self.finish_moments(delta,z,z['dU'][:,1])
        return z


def tangent(star,delta,z):
    model = previous.tangent(star,delta,z)
    model.__class__ = CompositionTangent
    last,older,h,coefficients = star.composition_context
    model.projected_anchor = species(star,delta,last,older,h,coefficients,z)
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
    star.last_delta,star.last_state = delta.copy(),z
    star.composition_context = (previous_state,older,h,coefficients)
    star.composition_data = composition_response(star,delta,z)
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
    plan = HISTORY.bindings()
    times = parent.wall.prior.time_nodes(plan,1)
    with ProcessPoolExecutor(max_workers=15,initializer=e.worker_init) as pool:
        star = prior.initialize(pool)
        states = []
        for name in ['step-0004.npz','step-0005.npz','failed-iterate.npz']:
            cp = np.load(HISTORY.OUT/'path-1'/name)
            delta = cp['delta'].copy()
            y = star.base+delta
            star.material_cache = {e.material_key(row):aux for row,aux in zip(zip(y[:,0],y[:,1],y[:,5:]),cp['aux'])}
            states.append((delta,star.evaluate(delta)))
        older,last,(delta,z) = states
        h = times[6]-times[5]
        coefficients = weights(h,times[5]-times[4])
        value,z = residual(star,delta,last,older,h,coefficients)
        initial = float(np.max(abs(value)/ATOL))
        assert initial == 72.5041790539357
        star.composition_context = (last,older,h,coefficients)
        star.composition_data = composition_response(star,delta,z)
        scores = {}
        for name,make in [('without_composition',previous.tangent),('native_composition',tangent)]:
            factor = splu(jacobian(make(star,delta,z),delta,last,older,h,coefficients))
            change = np.zeros_like(delta)
            change[:,:5] = -factor.solve(np.asarray(value[:,:5]/SCALE,float).ravel()).reshape(star.n,5).astype(ld)*SCALE
            assert np.max(abs(change[:,:2])) < .1 and np.max(abs(change[:,2])) < .01
            trial = delta+change
            trial[:,5:] = species(star,trial,last,older,h,coefficients,z)
            native,state = residual(star,trial,last,older,h,coefficients)
            scores[name] = float(np.max(abs(native)/ATOL))
        assert scores['without_composition'] > 1 and scores['native_composition'] <= 1,scores
        cp4,cp5 = [np.load(HISTORY.OUT/f'path-1/step-{step:04d}.npz') for step in [4,5]]
        c0,c1,c2 = coefficients
        area = 4*np.pi*star.rf**2
        energy = (-c1*cp5['integrated_energy_flux']-c2*cp4['integrated_energy_flux']+h*e.C*area*state['fluxes'][1])/c0
        baryon = (-c1*cp5['integrated_baryon_flux']-c2*cp4['integrated_baryon_flux']+h*e.C*area*state['fluxes'][0])/c0
        conservation,_,_ = budget(star,state,energy,baryon)
        cone = parent.cones(state)
        assert budget_passed(conservation,plan) and cone['sampled_cone_inside_light_cone']
        result = dict(classification='Counterexample candidate',passed=True,initial_residual=initial,
            corrected_native_residuals=scores,native_composition_evaluations=star.composition_data['evaluations'],
            conservation=conservation,cone=cone,actual_native_stage_check=True,full_trajectory_completed=False,
            scope='One native correction from the frozen step-six failure, including all 31 residuals and the actual BDF budgets. This is not a full-time, physical EOS or continuum certificate.')
    print('NATIVE COMPOSITION COUPLED CHECK',json.dumps(result),flush=True)
    return result


def prepare():
    assert not OUT.exists()
    plan = previous.bindings()
    assert plan['operator_check']['passed']
    checked = check()
    prefix = json.loads((previous.OUT/'path-1/progress.json').read_text())['step']
    assert prefix >= 5
    files = [previous.OUT/f'path-1/step-{step:04d}.npz' for step in range(prefix+1)]
    files += [previous.OUT/'plan.json',HISTORY.OUT/'path-1/failed-iterate.npz',
              HISTORY.OUT/'path-1/failure.json',Path(__file__)]
    for path in files:
        plan['bindings'][path.relative_to(e.ROOT).as_posix()] = e.digest(path)
    plan.update(checkpoint=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        candidate='Same native conservative evolution with transport-induced native EOS composition response in its iteration operator.',
        nonlinear_operator='At each time-stage initial iterate, evaluate the actual EOS and both opacity channels along the two neighboring composition directions. Retain those local response coefficients during that stage. Each iteration differentiates the species projection and those native responses while retaining the selected material flux branch. Actual acceptance still uses all 31 native equations. Same 24 iterations, eight backtracks, Anderson history, physical model, time grid and scientific gates.',
        composition_probe='Convex mixtures with at most 1e-5 maximum absolute species change, mixture fraction at most one; pressure/opacity log differences and absolute specific-energy differences. Batched local two-direction pseudoinverse rcond=1e-10. Exact-zero directions require no EOS calls.',
        iteration_limitations='Lagged local two-direction composition basis, five-field radius-two sparse approximation to the projected transport response, and frozen metric derivatives are only iteration approximations. No change or reduction of the actual 26-species native state.',
        imported_prefix_steps=prefix,
        imported_prefix='Replay the bound accepted parent coarse states with bitwise state/metric and native BDF residual checks; then start each new time stage from its preceding accepted state. Refined paths begin at the original initial state.',
        parent_plan_sha256=e.digest(previous.OUT/'plan.json'),operator_check=checked)
    OUT.mkdir()
    e.write(OUT/'plan.json',plan)
    print('FROZEN NATIVE COMPOSITION ITERATION; IMPORTED ACCEPTED STEPS',prefix,flush=True)


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
