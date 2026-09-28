def newton(cell,target,inner_shift):
    x=np.array([0.,.003,0.]);history=[];success=False
    for iteration in range(PLAN['max_iterations']):
        value,J,_,_=cell.evaluate(x,inner_shift);residual=np.asarray(value-target,float);J=np.asarray(J,float)
        merit=float(np.max(abs(residual)));row=dict(iteration=iteration,x=x.tolist(),merit=merit,
            condition=float(np.linalg.cond(J)),trials=[]);history.append(row)
        if merit<=PLAN['Newton_scaled_residual']:success=True;break
        step=np.linalg.solve(J,-residual);step*=min(1.,.1/max(abs(step[0]),1e-300),.1/max(abs(step[1]),1e-300),.01/max(abs(step[2]),1e-300))
        for backtrack in range(PLAN['max_backtracks']):
            proposed=x+step*2.**(-backtrack)
            if abs(proposed[2])>=.5:continue
            new_merit=float(np.max(abs(cell.evaluate(proposed,inner_shift)[0]-target)))
            row['trials'].append(dict(backtrack=backtrack,merit=new_merit))
            if new_merit<merit:x=proposed;row['accepted_backtrack']=backtrack;break
        else:break
    assert all(b['merit']<a['merit'] for a,b in zip(history,history[1:]))
    return x,success,history

