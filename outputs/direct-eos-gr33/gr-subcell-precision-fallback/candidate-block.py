def block(job):
    label,cells=job;plan=bindings();ld=np.longdouble
    data=dict(np.load(g.OUT/'reference-state.npz'));saved=dict(np.load(structure.OUT/'path-4.npz'))
    solver=structure.Structure(data,4);m=solver.mat;n=len(m.lp);B=ld(m.B)
    stats=dict(calls=0,evaluations=0,maximum_score=0.,label=label)
    inverse=make_inverse(stats)
    arrays={number:[] for number in plan['nodes']};rows=[]
    c=ld(g.c.gr.C)*100;G=ld(g.c.gr.G)*1000
    for i in cells:
        outside=i<m.split;j=i if outside else n-1-i
        branch=saved['outer'] if outside else saved['inner'];start=branch[j,1:].astype(ld)
        low,high=map(float,branch[j:j+2,0]);assert 0<low<high
        lower=ld(0) if j==0 else ld(low);upper=ld(high)
        finish=branch[j+1,1:].astype(ld)
        scale=np.maximum(abs(finish-start),[ld('1e-20'),ld('1e-25'),ld('1e-14')])
        def rhs(x,delta):
            absolute=np.asarray(start+scale*delta,float)
            return np.asarray(solver.rhs(x,absolute,i,float(B),outside)/scale,float)
        result=solve_ivp(rhs,(np.log(low),np.log(high)),np.zeros(3),method='DOP853',
            rtol=plan['rtol'],atol=plan['atol'],max_step=(np.log(high)-np.log(low))/4,dense_output=True)
        assert result.success,result.message
        endpoint=start+scale*result.y[:,-1]
        endpoint_error=float(abs(endpoint-finish).max())
        cellrows=[]
        for number in plan['nodes']:
            knots,gauss=np.polynomial.legendre.leggauss(number);knots=knots.astype(ld);gauss=gauss.astype(ld)
            if outside:
                coordinate=(lower+upper)/2+(upper-lower)*knots/2
                dB=ld(m.baryon_g)*(upper-lower)*gauss/2
            else:
                left=np.cbrt(lower);right=np.cbrt(upper)
                z=(left+right)/2+(right-left)*knots/2;coordinate=z**3
                dB=ld(m.baryon_g)*3*z*z*(right-left)*gauss/2
            normal=coordinate>=low;absolute=np.zeros((number,3),dtype=ld)
            absolute[normal]=start+result.sol(np.asarray(np.log(coordinate[normal]),float)).T*scale
            if np.any(~normal):
                # Original analytic seed, never extrapolate the numerical dense path.
                if outside:
                    surface=saved['faces'][0];slope=(start-surface)/ld(low)
                    absolute[~normal]=surface+coordinate[~normal,None]*slope
                else:
                    pc=ld(saved['faces'][-1,2]);ratio=coordinate[~normal]/ld(low)
                    absolute[~normal,0]=start[0]*np.cbrt(ratio)
                    absolute[~normal,1]=start[1]*ratio
                    absolute[~normal,2]=pc+np.log1p(np.expm1(start[2]-pc)*ratio**(ld(2)/3))
            radius=absolute[:,0]*ld(m.R)*100;mass=absolute[:,1]*B*100
            aa=1/np.sqrt(1-2*mass/radius);assert np.all(aa>=1)
            native=[];temperatures=[]
            for lp in absolute[:,2]:
                value,lt,_=inverse(solver.eos,float(lp),solver.ref[i,3],m.eps[i],m.lt[i])
                native.append(value);temperatures.append(lt)
            native=np.array(native);rho=native[:,0].astype(ld)
            proper=dB/rho;weights=proper/aa
            assert np.all(weights>0) and np.all(native[:,10]>0)
            C=ld(m.cx[i])*c*c
            mass_integral=np.sum((C+native[:,2])*dB/aa)*G/c**4
            # Independent geometric volume: use a factored radius difference.
            radii=saved['radius_m'][i:i+2].astype(ld)*100;ro,ri=radii
            volume=4*ld(np.pi)/3*(ro-ri)*(ro*ro+ro*ri+ri*ri)
            coordinate_volume=np.sum(weights)
            inventory=np.sum(dB)/ld(data['dm'][i])-1
            mass_gap=mass_integral/(saved['shell_mass_geom_m'][i]*100)-1
            volume_gap=coordinate_volume/volume-1
            row=dict(cell=i,nodes=number,endpoint_difference=endpoint_error,
                coordinate_inventory_relative_difference=float(inventory),
                native_shell_mass_relative_difference=float(mass_gap),
                native_coordinate_volume_relative_difference=float(volume_gap),
                analytic_seed_nodes=int(np.count_nonzero(~normal)),
                coordinate_volume_cm3=float(coordinate_volume),
                proper_internal_energy_erg=float(np.sum(native[:,2]*dB)),
                native_shell_mass_geom_cm=float(mass_integral))
            row['passed']=bool(endpoint_error<=plan['endpoint_tolerance'] and
                abs(inventory)<=plan['inventory_relative_tolerance'] and
                abs(mass_gap)<=plan['shell_mass_relative_tolerance'] and
                abs(volume_gap)<=plan['volume_relative_tolerance'])
            cellrows.append(row)
            arrays[number].append(dict(eos=native,lnT=np.array(temperatures),radius_cm=radius,mass_geom_cm=mass,
                metric_a=aa,coordinate_weights_cm3=weights,proper_weights_cm3=proper,baryon_weights_g=dB,
                logP=absolute[:,2],C_X=m.cx[i]))
        gap=max(abs(cellrows[1][key]/cellrows[0][key]-1) for key in
            ['coordinate_volume_cm3','proper_internal_energy_erg','native_shell_mass_geom_cm'])
        rows.append(dict(cell=i,quadratures=cellrows,finite_quadrature_relative_difference=gap,
            passed=all(r['passed'] for r in cellrows) and gap<=plan['finite_quadrature_relative_tolerance']))
    for number,items in arrays.items():
        np.savez_compressed(OUT/f'{label}-nodes-{number}.npz',cells=np.array(cells),
            **{key:np.array([a[key] for a in items]) for key in items[0]})
    record=dict(classification='Counterexample candidate',rows=rows,all_passed=all(r['passed'] for r in rows),
        entropy_roots=stats,table_or_fallback_calls=solver.calls,geometry_direct_calls=solver.direct_calls)
    save(label+'.json',record);print('FULL SUBCELL',label,len(rows),record['all_passed'],stats,flush=True)
    return record
