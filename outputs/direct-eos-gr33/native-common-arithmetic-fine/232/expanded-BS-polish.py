def polish(m,op,rhs,sol):
    scale=max(np.linalg.norm(rhs),LD('1e-290'));goal=LD('1e-14')*scale
    cv=inspect.getclosurevars(op._CustomLinearOperator__matvec_impl).nonlocals
    blocks=cv['blocks_by_component'];dim=len(rhs)//2;records=[]
    for attempt in range(4):
        residual=rhs-op.matvec(sol);norm=np.linalg.norm(residual)
        if norm<goal:break
        if norm>goal*10000:break
        target=[]
        for stage,row in enumerate(residual.reshape(2,dim)):
            _,gas=m.unpack(row)
            target += [(stage,int(k),component) for component in [2,3] for k in np.flatnonzero(abs(gas[:,component])>goal/8)]
        if not target:break
        columns=sorted({col for stage,cell,component in target for col,value in blocks[component][stage][cell] if value and col%4!=2})
        candidates=[]
        for col in columns:
            index=(col//(4*m.n))*dim+m.size+col%(4*m.n)
            quantum=abs(np.spacing(sol[index])) if sol[index] else LD(0)
            score=max((abs(LD(str(dict(blocks[component][stage][cell]).get(col,0))))*quantum for stage,cell,component in target),default=LD(0))
            if score<goal/4:candidates.append((float(score),index,quantum))
        candidates.sort();selected=[];images=[];norms=[]
        for _,index,quantum in candidates[:24]:
            unit=np.zeros_like(sol);unit[index]=1;image=op.matvec(unit);cn=np.linalg.norm(image)
            if not np.isfinite(cn) or cn==0 or cn*quantum>=goal/4:continue
            selected.append(index);norms.append(cn);images.append(sparse.csc_matrix(np.asarray(image/cn,float)[:,None]))
        if not images:break
        matrix=sparse.hstack(images,format='csc');rows=np.unique(matrix.nonzero()[0])
        delta,_,rank,_=np.linalg.lstsq(matrix[rows].toarray(),np.asarray(residual[rows],float),rcond=1e-14)
        candidate=sol.copy();candidate[selected]+=np.asarray(delta,LD)/np.asarray(norms,LD)
        actual=rhs-op.matvec(candidate);newnorm=np.linalg.norm(actual)
        moments=joint.physical_norm(m,actual)/joint.scales(m,rhs,candidate)
        accepted=bool(newnorm<norm and max(moments)<1e-13 and np.all(np.isfinite(candidate)))
        records.append(dict(components=sorted({v[2] for v in target}),before=float(norm/scale),after=float(newnorm/scale),rank=int(rank),columns=selected,
            rows=len(rows),physical_relative=moments.astype(float).tolist(),accepted=accepted))
        if not accepted:break
        sol=candidate
    path=OUT/'polish.json';old=read(path) if path.exists() else []
    write(path,old+[dict(classification='Counterexample candidate',attempts=records)])
    return sol
