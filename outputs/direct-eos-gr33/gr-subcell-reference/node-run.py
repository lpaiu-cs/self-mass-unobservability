def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    state,aux=preflight.heat.old.char.prior.micro.inputs();ref=dict(np.load(g.OUT/'reference-state.npz'))
    eos=g.EOS();stats=dict(calls=0,evaluations=0,maximum_score=0.,label='direct-cell-entropy')
    # Reuse the exact checked inversion code, changing only its failure-log
    # destination and counters in a private globals mapping. No live module,
    # EOS implementation, stopping rule or historical output is changed.
    source=g.audit.strict_invert
    inverse=FunctionType(source.__code__,dict(source.__globals__,ROOT_STATS=stats,e=SimpleNamespace(OUT=OUT,save=save)))
    c=g.c.gr.C*100;G=g.c.gr.G*1000;records=[]
    for i in plan['cells']:
        X=state['X'][i];entropy=ref['s_B'][i];lp0=float(state['logP'][i]);lt0=float(state['lnT'][i])
        amid,tmid,_=inverse(eos,lp0,entropy,X,lt0)
        Hmid=amid[2]+amid[1]/amid[0];CX=state['CX'][i]
        cpT=amid[10]-amid[1]/amid[0]*amid[8]
        grad_ad=amid[1]/amid[0]*(-amid[8])/cpT
        calls=[]
        def at(lp):
            a,t,used=inverse(eos,float(lp),entropy,X,tmid+grad_ad*(float(lp)-lp0))
            return a,t
        boundaries=[]
        for face in [i,i+1]:
            increment=(CX*c*c+Hmid)*np.expm1(state['nu'][i]-state['nu_faces'][face])
            target=Hmid+increment;budget=max(2.,32*np.spacing(abs(target)))
            lp=lp0+np.clip(increment/(amid[1]/amid[0]),-1.,1.);best=None
            for iteration in range(32):
                a,t=at(lp);defect=a[2]+a[1]/a[0]-target
                if best is None or abs(defect)<best[0]:best=(abs(defect),lp,a,t,defect)
                if abs(defect)<budget*.25:break
                proposed=lp-np.clip(defect/(a[1]/a[0]),-.3,.3)
                if proposed==lp:break
                lp=proposed
            if best[0]>budget*.25:
                for direction in [-np.inf,np.inf]:
                    lp=best[1]
                    for _ in range(3):
                        lp=np.nextafter(lp,direction);a,t=at(lp);defect=a[2]+a[1]/a[0]-target
                        if abs(defect)<best[0]:best=(abs(defect),lp,a,t,defect)
            record=dict(face=face,logP=float(best[1]),lnT=float(best[3]),enthalpy_residual=float(best[4]),
                budget_erg_g=float(budget),passed=bool(best[0]<=budget))
            boundaries.append(record);assert record['passed'],record
        save(f'cell-{i}-boundaries.json',dict(classification='Counterexample candidate',cell=i,rows=boundaries))
        left,right=[r['logP'] for r in boundaries];assert left<right
        rout=float(state['radius_faces_m'][i]*100);rin=float(state['radius_faces_m'][i+1]*100)
        mout=float(state['mass_faces_geom'][i]*100);minner=float(state['mass_faces_geom'][i+1]*100)
        width=rout-rin;dm=float(state['dm'][i]);mass_scale=G*dm/c**2
        volume_scale=dm/amid[0];u_scale=max(abs(amid[2]),amid[10]);Pmid=amid[1]
        def rhs(lp,y):
            a,t=at(lp);rho,P,u=a[:3]
            radius=rout-width*y[0];mass=mout-mass_scale*y[1];f=1-2*mass/radius
            assert radius>0 and mass>0 and f>0
            pgeom=G*P/c**4;epsgeom=G*rho*(CX*c*c+u)/c**4
            dr=-pgeom*radius*(radius-2*mass)/((epsgeom+pgeom)*(mass+4*np.pi*radius**3*pgeom))
            dvolume=-4*np.pi*radius**2*dr/np.sqrt(f);db=rho*dvolume
            return np.array([-dr/width,-4*np.pi*radius**2*epsgeom*dr/mass_scale,
                db/dm,dvolume/volume_scale,u*db/(dm*u_scale),P*dvolume/(Pmid*volume_scale)])
        paths=[]
        for j,tolerance in enumerate(plan['relative_tolerances']):
            result=solve_ivp(rhs,(left,right),np.zeros(6),method='DOP853',rtol=tolerance,
                atol=plan['absolute_scaled_tolerance'],max_step=(right-left)/4,dense_output=True)
            assert result.success,result.message;y=result.y[:,-1]
            integrated_baryon=dm*y[2];volume=volume_scale*y[3];uint=dm*u_scale*y[4]
            pvolume=Pmid*volume_scale*y[5];rho_avg=integrated_baryon/volume;u_avg=uint/integrated_baryon
            # Evaluate a uniform state with the same mean density and initial
            # entropy. Its energy need not equal the nonuniform cell average.
            t=tmid;best=None
            for _ in range(15):
                a=eos(2,float(np.log(rho_avg)),float(t),X)
                defect=np.exp(t)*(a[3]-entropy);budget=max(2.,32*np.spacing(abs(a[2])))
                if best is None or abs(defect)<best[0]:best=(abs(defect),a.copy(),t)
                if abs(defect)<=budget:break
                proposed=t-np.clip(defect/a[10],-.15,.15)
                if proposed==t:break
                t=proposed
            assert best[0]<=budget,(i,j,'mean-density entropy root',best[0],budget)
            row=dict(classification='Counterexample candidate',cell=i,tolerance=tolerance,nfev=result.nfev,
                integrated_baryon_relative_difference=float(y[2]-1),
                radius_endpoint_difference_over_cell_width=float(1-y[0]),
                integrated_geometric_mass_cm=float(y[1]*mass_scale),
                stored_face_geometric_mass_difference_cm=float(mout-minner),
                mean_baryon_density_cgs=float(rho_avg),mean_internal_energy_erg_g=float(u_avg),
                mean_pressure_dyn_cm2=float(pvolume/volume),
                uniform_same_mean_density_entropy_energy_difference_erg_g=float(best[1][2]-u_avg),
                uniform_same_mean_density_entropy_energy_difference_over_cvT=float((best[1][2]-u_avg)/best[1][10]),
                uniform_same_mean_density_entropy_pressure_relative_difference=float(best[1][1]/(pvolume/volume)-1),
                same_mean_density_entropy_lnT=float(best[2]))
            np.savez_compressed(OUT/f'cell-{i}-path-{j}.npz',logP=result.t,scaled_integrals=result.y)
            retain(i,j,result,left,right,at,rout,width,mout,mass_scale,CX,dm,volume_scale,u_scale,Pmid);paths.append(row);print('DIRECT CELL INTEGRAL',row,flush=True)
        difference=float(abs(np.load(OUT/f'cell-{i}-path-1.npz')['scaled_integrals'][:,-1]-
            np.load(OUT/f'cell-{i}-path-0.npz')['scaled_integrals'][:,-1]).max())
        records.append(dict(cell=i,paths=paths,maximum_finite_tolerance_difference_scaled=difference))
        save('progress.json',dict(classification='Counterexample candidate',records=records,entropy_root_statistics=stats))
    save('result.json',dict(classification='Counterexample candidate',completed=True,records=records,
        entropy_root_statistics=stats,full_star_recomputed=False,physical_EOS_certified=False,
        continuous_radial_or_native_error_certified=False,full_GR_evolution=False))
    save('manifest.json',dict(classification='Counterexample candidate',
        sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify()
