"""Test a target-only mechanical predictor on two frozen inverse failures."""
import json,subprocess,sys
import numpy as np
import gr_metric_coupled_subcell_newton as prior

ROOT=prior.ROOT;OUT=prior.OUT.parent/'gr-metric-predictor-diagnostic';sha=prior.sha


def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2)+'\n')


def mechanical_predictor(cell,target,inner_shift):
    # Temperature is fixed only during this search predictor, not in the equations.
    x=np.array([0.,.003,0.]);history=[];indices=[0,2]
    for iteration in range(8):
        value,J,_,_=cell.evaluate(x,inner_shift)
        f=np.asarray(value-target,float)[indices];A=np.asarray(J,float)[np.ix_(indices,indices)]
        merit=float(np.max(abs(f)))
        row=dict(iteration=iteration,x=x.tolist(),mechanical_merit=merit,
            full_merit=float(np.max(abs(value-target))),condition=float(np.linalg.cond(A)),trials=[])
        history.append(row)
        if merit<=1e-12:break
        step=np.linalg.solve(A,-f)
        step*=min(1.,.1/max(abs(step[0]),1e-300),.01/max(abs(step[1]),1e-300))
        for backtrack in range(24):
            proposed=x.copy();proposed[indices]+=step*2.**(-backtrack)
            if abs(proposed[2])>=.5:continue
            score=float(np.max(abs(cell.evaluate(proposed,inner_shift)[0]-target)[indices]))
            row['trials'].append(dict(backtrack=backtrack,mechanical_merit=score))
            if score<merit:x=proposed;row['accepted_backtrack']=backtrack;break
        else:break
    cell.predictor_history=history
    assert all(b['mechanical_merit']<a['mechanical_merit'] for a,b in zip(history,history[1:]))
    return x


def transformed():
    source,_=prior.source();start=source.index('def newton(');end=source.index('\ndef one_cell(',start)
    before=source[start:end]
    old='x=np.array([0.,.003,0.]);history=[];success=False'
    new='x=mechanical_predictor(cell,target,inner_shift);history=[];success=False'
    assert before.count(old)==1;after=before.replace(old,new)
    assert after.replace(new,old)==before
    return before,after


def prepare():
    assert not OUT.exists();OUT.mkdir();before,after=transformed()
    (OUT/'original-newton.py').write_text(before);(OUT/'predicted-newton.py').write_text(after)
    files=[ROOT/'verification/gr_metric_predictor_diagnostic.py',
        prior.OUT/'plan.json',prior.OUT/'candidate.py',prior.OUT/'preflight.json',
        prior.OUT/'block-0064/result.json',prior.OUT/'block-0080/result.json',
        OUT/'original-newton.py',OUT/'predicted-newton.py']
    save('plan.json',dict(classification='Counterexample candidate',
        checkpoint=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        bindings={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        cases=[dict(cell=76,nodes=8,case=0),dict(cell=83,nodes=8,case=2),dict(cell=0,nodes=8,case=1)],
        predictor_iterations=8,predictor_mechanical_gate=1e-12,
        change='Begin at the same [0,.003,0]. At fixed initial temperature shift, first solve only the original baryon and momentum rows for density and velocity, with the full metric chain. Then run the original 30-iteration full Newton solve unchanged. The predictor receives only the target and inner mass shift, not the manufactured true primitive. Its full residual need not decrease; only its two mechanical rows decrease. The original full-phase monotonic residual rule and all root/physical acceptance gates remain unchanged.',
        controls='For both failures, check the analytic Jacobian at the initial and last archived iterate using both original finite-difference steps and the original 1e-4 row-relative gate. Include the original zero-target cell as a positive control. Compare exact target moment dyadics with the archived original in every case.',
        boundary='Post-failure numerical diagnosis on three explicitly selected cases; no frozen population rescue, all-cell convergence, global uniqueness, native/continuous error or finite GR evolution claim.'))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert sha(ROOT/rel)==digest,rel
    before,after=transformed()
    assert before==(OUT/'original-newton.py').read_text() and after==(OUT/'predicted-newton.py').read_text()
    return p,after


def archived(cell,n,case):
    path=prior.OUT/'preflight.json' if cell==0 else prior.OUT/f'block-{cell//16*16:04d}/result.json'
    row=next(c for c in json.loads(path.read_text())['rows'] if c['cell']==cell)
    return next(r for r in row['rows'] if r['nodes']==n and r['case']==case)


def run():
    p,source=bindings();assert not (OUT/'result.json').exists();engine=prior.engine();engine.initialize()
    namespace=dict(engine.__dict__,mechanical_predictor=mechanical_predictor)
    exec(compile(source,str(OUT/'predicted-newton.py'),'exec'),namespace);solver=namespace['newton'];rows=[]
    for task in p['cases']:
        i,n,k=task['cell'],task['nodes'],task['case'];cell=engine.Cell(i,n)
        truth=engine.PLAN['probes'][k];shift=engine.PLAN['inner_log_mass_shifts'][k]
        target,_,moments,_=cell.evaluate(truth,shift);old=archived(i,n,k)
        exact=[[str(a),str(b)] for a,b in (q.as_integer_ratio() for q in moments)]
        assert exact==old['target_changes_exact'];controls=[]
        if not old['passed']:
            for label,x in [('initial',[0.,.003,0.]),('archived-last',old['history'][-1]['x'])]:
                x=np.array(x);_,J,_,_=cell.evaluate(x,shift)
                for h in engine.PLAN['finite_Jacobian_steps']:
                    fd=[]
                    for j in range(3):
                        dx=np.zeros(3);dx[j]=h
                        fd.append((cell.evaluate(x+dx,shift)[0]-cell.evaluate(x-dx,shift)[0])/(2*engine.ld(h)))
                    fd=np.array(fd).T;den=np.maximum(np.max(abs(J),axis=1),engine.ld('1e-30'))
                    score=float(np.max(abs(fd-J)/den[:,None]))
                    controls.append(dict(location=label,step=h,score=score,passed=score<=engine.PLAN['finite_Jacobian_relative_gate'],
                        full_Jacobian_condition=float(np.linalg.cond(np.asarray(J,float)))))
        recovered,success,history=solver(cell,target,shift)
        residual=float(np.max(abs(cell.evaluate(recovered,shift)[0]-target)));error=abs(recovered-np.array(truth))
        passed=bool(success and error[0]<=engine.PLAN['root_log_density_tolerance']
            and error[1]<=engine.PLAN['root_logT_tolerance'] and error[2]<=engine.PLAN['root_velocity_tolerance']
            and residual<=engine.PLAN['scaled_residual_tolerance'])
        rows.append(dict(**task,original_passed=old['passed'],original_scaled_residual=old['scaled_residual'],
            target_changes_exact=exact,errors=error.tolist(),scaled_residual=residual,solver_success=success,passed=passed,
            Jacobian_controls=controls,predictor_history=cell.predictor_history,history=history))
        save('progress.json',dict(rows=rows));print('MECHANICAL PREDICTOR',task,passed,residual,flush=True)
    save('result.json',dict(classification='Counterexample candidate',completed=True,rows=rows,
        all_passed=all(r['passed'] and all(c['passed'] for c in r['Jacobian_controls']) for r in rows),
        original_failures_preserved=True,full_GR_evolution=False,physical_EOS_certified=False))
    save('manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    p,_=bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    result=json.loads((OUT/'result.json').read_text());assert result['completed'] and len(result['rows'])==len(p['cases'])
    for row,task in zip(result['rows'],p['cases']):
        assert all(row[k]==v for k,v in task.items())
        old=archived(task['cell'],task['nodes'],task['case'])
        assert old['target_changes_exact']==row['target_changes_exact'] and old['passed']==row['original_passed']
        assert row['predictor_history'][0]['x']==[0.,.003,0.]
        for key,merit in [('predictor_history','mechanical_merit'),('history','merit')]:
            assert all(b[merit]<a[merit] for a,b in zip(row[key],row[key][1:]))
    print('PASS frozen predictor diagnosis bindings and phase-specific descent; inspect all_passed:',result['all_passed'],flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
