def solve(self,n,label):
    start=time.monotonic();m=self.model;z=np.clongdouble(-1j*n*self.omega);w=float(-z.imag)
    G=(self.K+complex(z*z)*self.M).tocsc();Gx=self.Kx+z*z*self.Mx
    rhs=-(self.kl+z*z*self.ml).astype(np.clongdouble)
    sc=np.sqrt(abs(G.diagonal()));D=diags(1/sc);lu=splu((D@G@D).astype(complex).tocsc(),permc_spec='NATURAL')
    qa=(lu.solve(np.asarray(rhs/sc,complex))/sc).astype(np.clongdouble)
    for _ in range(3):qa+=(lu.solve(np.asarray((rhs-self.Kx@qa-z*z*(self.Mx@qa))/sc,complex))/sc).astype(np.clongdouble)
    qa_error=float(np.max(abs(rhs-self.Kx@qa-z*z*(self.Mx@qa))/(abs(rhs)+self.absK@abs(qa)+abs(z*z)*(self.absM@abs(qa))+1e-100)))
    gain=np.sum(self.conductance*m.heat.geometry.tc*self.lam/(z*(z+self.lam)),axis=1)
    B=self.original_force-z*z*(self.Mx@self.H)
    C=diags(1/gain)-self.GEraw
    block=bmat([[Gx,-B@diags(self.energy_scale)],[-self.Gq,C@diags(self.energy_scale)]],format='csc')
    forcing=self.Gq@qa+self.gL
    bvec=np.r_[np.zeros(m.size,np.clongdouble),forcing]
    perm=self.permutation;A=block[perm,:][:,perm].astype(complex).tocsc();bp=bvec[perm]
    row=np.asarray(abs(A).max(axis=1).toarray()).ravel();Ar=diags(1/row)@A
    col=np.asarray(abs(Ar).max(axis=0).toarray()).ravel();scaled=(Ar@diags(1/col)).tocsc()
    factor=splu(scaled,permc_spec='NATURAL')
    def invert(b):
        yp=factor.solve(np.asarray(b[perm]/row,complex))/col
        y=np.empty_like(yp);y[perm]=yp
        return y.astype(np.clongdouble)
    answer=invert(bvec)
    for _ in range(4):answer+=invert(bvec-block@answer)
    residual=bvec-block@answer
    error=float(np.max(abs(residual)/(abs(bvec)+abs(block)@abs(answer)+1e-100)))
    dq=answer[:m.size];E=answer[m.size:]*self.energy_scale
    heat_drive=self.Gq@(qa+dq)+self.gL+self.GEraw@E
    heat_defect=E-gain*heat_drive
    heat_error=float(abs(heat_defect).max()/max(abs(E).max(),1e-100))
    # Full physical-coordinate and heat equations, independent of the
    # assembled correction block and its row/column scaling.
    gr_defect=self.Kx@dq+z*z*(self.Mx@dq)-B@E
    gr_error=float(np.max(abs(gr_defect)/(self.absK@abs(dq)+abs(z*z)*(self.absM@abs(dq))+abs(B)@abs(E)+1e-100)))
    assert max(error,qa_error,gr_error)<1e-9 and heat_error<1e-9,(error,qa_error,gr_error,heat_error)
    Ra=-rhs@qa+self.kll+z*z*self.mll
    dR=-rhs@dq-self.fL@E
    Za=Ra/(2*self.F);dZ=dR/(2*self.F)
    wave=orbit.exterior.outgoing(self.mu,self.flux,2*w);h=wave['h'];Zo=wave['impedance']
    drive=np.exp(-1j*w*(1+self.delay))/(self.F*h);Da=Za-Zo
    aa=drive/Da;af=drive/(Da+dZ);da=-drive*dZ/(Da*(Da+dZ))
    tail=2*da/h*np.exp(-1j*w*self.delay);charge=-self.R/self.mass*tail
    fullE=np.zeros(len(m.heat.edges),np.clongdouble);fullE[self.active]=af*E
    temp=af*(self.Tq@(qa+dq)+self.TL)+self.TE@fullE
    contrast_temp=af*(self.Tq@dq)+da*(self.Tq@qa+self.TL)+self.TE@fullE
    balance=float(abs(np.sum(-np.diff(fullE),dtype=np.clongdouble))/max(abs(fullE).max(),1e-100))
    # Real dissipated quadratic form of the same positive thermal poles.
    flux=z*E/m.heat.geometry.tc
    dissipation=float(np.real(np.vdot(heat_drive,flux)))
    assert dissipation>=0 and balance<2e-13
    old=json.loads((orbit.OUT/'result.json').read_text());incident=old['rows'][n]['drive_amplitude']
    pair=lambda v:[float(v.real),float(v.imag)]
    np.savez_compressed(OUT/f'{label}-{n}.npz',qa=qa,dq=dq,E=E,temperature=temp,temperature_correction=contrast_temp,
        active_faces=self.active,native_radius=m.original.native,grid=m.grid,indices=m.indices)
    result=dict(classification='Counterexample candidate',harmonic=n,degree=m.degree,seconds=time.monotonic()-start,
        scalar_boundary_value=pair(af),adiabatic_boundary_value=pair(aa),adiabatic_impedance=pair(Za),thermal_impedance_difference=pair(dZ),
        outgoing_conduction_contribution=pair(tail),radiative_charge_gain=pair(charge),
        actual_drive_amplitude=incident,delta_alpha_radiative_over_phi0=float(abs(charge)*incident/.001),
        maximum_actual_Eulerian_delta_lnT=float(abs(temp).max()*incident),
        maximum_actual_thermal_temperature_correction=float(abs(contrast_temp).max()*incident),
        linear_residual=error,adiabatic_residual=qa_error,physical_GR_residual=gr_error,heat_law_residual=heat_error,
        heat_balance=balance,positive_pole_dissipation=dissipation,
        active_faces=len(self.active),dofs=block.shape[0],matrix_nonzeros=block.nnz,
        factor_nonzeros=factor.L.nnz+factor.U.nnz,memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
        full_radial_photons=False,full_nonlinear_or_thermal_background=False,full_goal_complete=False)
    write(OUT/f'{label}-{n}.json',result);print('ORBITAL HEAT',label,n,json.dumps(result),flush=True)
    return result
