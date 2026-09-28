def neighbours(m,op,rhs,sol,log,triplets=False):
    scale=max(np.linalg.norm(rhs),LD('1e-290'));goal=LD('1e-14')*scale;dim=len(rhs)//2
    blocks=inspect.getclosurevars(op._CustomLinearOperator__matvec_impl).nonlocals['blocks_by_component']
    for attempt in range(4):
        residual=rhs-op.matvec(sol);norm=np.linalg.norm(residual)
        if norm<goal:break
        target=[]
        for stage,row in enumerate(residual.reshape(2,dim)):
            _,gas=m.unpack(row)
            target += [(stage,int(k),c) for c in [2,3] for k in np.flatnonzero(abs(gas[:,c])>goal/8)]
        columns=sorted({col for stage,cell,c in target for col,value in blocks[c][stage][cell] if value})
        data=[];images=[]
        for col in columns:
            index=(col//(4*m.n))*dim+m.size+col%(4*m.n);quantum=abs(np.spacing(sol[index])) if sol[index] else LD(0)
            unit=np.zeros_like(sol);unit[index]=1;image=op.matvec(unit);cn=np.linalg.norm(image)
            if cn==0 or not np.isfinite(cn):continue
            data.append((index,quantum,cn));images.append(sparse.csc_matrix(np.asarray(image/cn,float)[:,None]))
        if not images:break
        matrix=sparse.hstack(images,format='csc');rows=np.unique(matrix.nonzero()[0]);M=matrix[rows].toarray();r=np.asarray(residual[rows],float)
        fine=[j for j,(_,q,cn) in enumerate(data) if q*cn<goal/4];hard=[j for j in range(len(data)) if j not in fine]
        if not fine or not hard:break
        F=M[:,fine];u,s,vt=np.linalg.svd(F,full_matrices=False);keep=s>s[0]*1e-14;Q=u[:,keep];pinv=(vt[keep].T/s[keep])@Q.T
        pr=r-Q@(Q.T@r);H=M[:,hard]*np.array([float(data[j][1]*data[j][2]) for j in hard]);PH=H-Q@(Q.T@H)
        sizes=np.linalg.norm(PH,axis=0);order=[int(i) for i in np.argsort(sizes) if sizes[i]>float(goal)/128][:8]
        order=list(range(len(hard)));guesses=integer_guesses(PH,pr,goal,integer_log)
        best=sol;best_norm=norm;tested=[]
        for key,score in sorted(guesses.items(),key=lambda item:item[1])[:24]:
            candidate=sol.copy();shifts=np.array(key)
            for j,k in enumerate(hard):
                if shifts[j]:index,q,_=data[k];candidate[index]+=LD(int(shifts[j]))*q
            delta=pinv@(r-H@shifts)
            for j,d in zip(fine,delta):index,_,cn=data[j];candidate[index]+=LD(d)/cn
            actual=rhs-op.matvec(candidate);newnorm=np.linalg.norm(actual);moments=joint.physical_norm(m,actual)/joint.scales(m,rhs,candidate)
            passed=bool(newnorm<best_norm and max(moments)<1e-13 and np.all(np.isfinite(candidate)))
            tested.append(dict(indices=[data[hard[j]][0] for j in np.flatnonzero(shifts)],ULPs=shifts[shifts!=0].tolist(),relative=float(newnorm/scale),physical=moments.astype(float).tolist(),improved=passed))
            if passed:best,best_norm=candidate,newnorm
            if best_norm<goal:break
        log.append(dict(pass_index=attempt,before=float(norm/scale),after=float(best_norm/scale),fine_columns=len(fine),hard_columns=len(hard),searched_hard_columns=len(order),projected_hard_norms=sizes.tolist(),tested=tested))
        if best_norm>=norm:break
        sol=best
    return sol
