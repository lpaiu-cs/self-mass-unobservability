def constrain(self,lr,lt,target,guess=None,target_molecules=None,tolerance=1e-10):
    if guess is not None:self.fields[:]=guess
    # Chemical fields enforce stage inventories, not physical rate laws.
    # Electron degeneracy/nonideal terms are iterated by the native EOS.
    active=target>target.sum(1)[:,None]*1e-18
    for iteration in range(16):
        a=self.snapshot(lr,lt);current=a['number_fractions']
        err=float(np.max(abs(current-target))/target.sum())
        if err<tolerance:return a,iteration+1,err
        for e,z in enumerate(self.Z):
            if target[e].sum()==0:continue
            anchor=int(np.argmax(target[e]))
            desired=np.where(active[e,:z+1],target[e,:z+1],target[e].sum()*1e-18)
            ratios=np.log(np.maximum(desired,1e-290)/np.maximum(current[e,:z+1],1e-290))
            ratios-=ratios[anchor]
            # Fix neutral affinity to zero. All ionic affinities use the
            # same neutral reference, including a tiny neutral inventory.
            ratios-=ratios[0]
            # Do not invert native underflow zeros into finite ions. Only
            # enforce initially populated coordinates; the remaining
            # coordinates equilibrate and their total is checked below.
            take=active[e,1:z+1]|(current[e,1:z+1]>target[e].sum()*1e-18)
            self.fields[self.starts[e]:self.starts[e+1]]+=np.where(take,ratios[1:],0.)
        if target_molecules is not None:
            now=a['molecular_H_fractions']
            correction=np.log(np.maximum(target_molecules,1e-290)/np.maximum(now,1e-290))
            correction-=2*np.log(target[0,0]/max(current[0,0],1e-290))
            self.fields[316:]+=np.where(now>1e-18,correction,0.)
    raise AssertionError(('Constrained ion solve',lr,lt,err))
