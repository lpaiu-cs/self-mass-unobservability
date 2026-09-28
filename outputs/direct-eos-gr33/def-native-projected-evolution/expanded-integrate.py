def integrate(d,model,order=8,kind='linear',stride=1):
    """Actual finite-volume source at retarded times, with no future values."""
    xg,wg=np.polynomial.legendre.leggauss(order);edges=d['edges'];r=(edges[:-1,None]+np.diff(edges)[:,None]*(xg+1)/2).ravel();wq=np.tile(wg/2,len(edges)-1)
    ids=np.repeat(np.arange(len(edges)-1),order);geom=green.Geometry(model.m);w,delay,_,B,_=geom(r-model.m.RJ)
    measure=(wq*r*r*B).reshape(len(edges)-1,order)
    wq=(measure/measure.sum(1)[:,None]).ravel()
    # The shared readout clock is fixed by the outer FACE, not the last
    # quadrature point (which changes with fluid resolution and Gauss order).
    outer_delay=float(geom(np.array([edges[-1]-model.m.RJ]))[1][0])
    end=prior.END-max(outer_delay,0.);times=np.linspace(0,end,129)
    trace=d['nonrest_trace_erg'];mass=np.asarray(d['baryon_g'],float);n=int(d['deep_cells']);t=d['t'][::stride]
    inputs=[trace[:,:n],trace[:,n:],mass[:,n:]];paths=[]
    # Integrate source histories first; each source cell is then evaluated at
    # its actual radial quadrature delays. No fitted relaxation parameter.
    for k,values in enumerate(inputs):
        subset=ids<n if k==0 else ids>=n;source_ids=ids[subset] if k==0 else ids[subset]-n
        antiderivative=green.polynomial(t,values[::stride],kind).antiderivative();rows=[]
        for tt in times:
            at=tt+delay[subset];cut=np.clip(at,0,t[-1]);interval=np.clip(np.searchsorted(antiderivative.x,cut,side='right')-1,0,len(antiderivative.x)-2);dt=cut-antiderivative.x[interval];v=np.zeros_like(dt)
            for coeff in antiderivative.c:v=v*dt+coeff[interval,source_ids]
            v[at<=0]=0;assert at.max()<=t[-1]+2e-15
            rows.append(np.sum(v*w[subset]*wq[subset],dtype=np.longdouble))
        factor=-G/(2*C**3) if k<2 else -G*float(d['cx'])/(2*C)
        paths.append(np.array(rows,float)*factor)
    # Complete represented baryon conservation at the *actual* inner fluid
    # face only as an explicit conditional point debit. Its physical profile
    # and thermal/stress response are not asserted to have been simulated.
    wi,di,*_=geom(np.array([float(d['inner_material_face_cm'])-model.m.RJ]));debit=-mass.sum(1)
    H=green.polynomial(t,debit[::stride],kind).antiderivative();debit_wave=-G*float(d['cx'])/(2*C)*wi[0]*green.paired(H,times+di[0])
    pieces=np.vstack([*paths,debit_wave]);Q=pieces.sum(0);M=float(d['M_cm'])
    return times,-Q/M,-pieces/M,dict(order=order,history=kind,stride=stride,endpoint=float(-Q[-1]/M),maximum_future_source_seconds=float(times[-1]+delay.max()-prior.END))
