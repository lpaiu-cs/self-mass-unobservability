"""Retain direct-EOS TOV quadrature nodes for conservative subcell recovery."""
import inspect,json,sys
import numpy as np
import sympy as s
import gr_direct_cell_integrals as original

g=original.g;OUT=g.OUT/'gr-subcell-reference';records=[]


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def transformed_source():
    source=inspect.getsource(original.run)
    replacements={
        "atol=plan['absolute_scaled_tolerance'],max_step=(right-left)/4)":
        "atol=plan['absolute_scaled_tolerance'],max_step=(right-left)/4,dense_output=True)",
        "paths.append(row);print('DIRECT CELL INTEGRAL',row,flush=True)":
        "retain(i,j,result,left,right,at,rout,width,mout,mass_scale,CX,dm,volume_scale,u_scale,Pmid);paths.append(row);print('DIRECT CELL INTEGRAL',row,flush=True)"}
    changed=source
    for old,new in replacements.items():assert changed.count(old)==1;changed=changed.replace(old,new)
    reverse=changed
    for old,new in replacements.items():reverse=reverse.replace(new,old)
    assert reverse==source
    return source,changed,replacements


def prepare():
    assert not OUT.exists();OUT.mkdir();original.verify()
    source,changed,replacements=transformed_source()
    (OUT/'original-run.py').write_text(source);(OUT/'node-run.py').write_text(changed)
    plan=json.loads((original.OUT/'plan.json').read_text())
    plan.update(checkpoint='6338252',node_counts=[8,16],baryon_mass_relative_tolerance=1e-7,
        geometric_mass_relative_tolerance=1e-7,internal_energy_relative_tolerance=1e-6,
        change='Retain the exact direct EOS TOV RHS and boundary equations. Request dense DOP853 output and export positive Gauss quadrature nodes with fresh native EOS at each node. Two exact source substitutions are stored and reversible.',
        scope='Four reconstructed nonuniform reference cells and two finite quadratures. No geometry or baryon normalization is applied to force a pass. Dense-output, native and continuum errors are not enclosed. Subcell-node export is not a full-star evolution.',
        substitutions=replacements)
    plan['bindings'].update({p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in [
        g.ROOT/'verification/gr_subcell_reference.py',original.OUT/'manifest.json',OUT/'node-run.py']})
    save('plan.json',plan)
    B,Cv,W,H,Z=s.symbols('B Cv W H Z',positive=True);A=s.symbols('A',real=True)
    matrix=s.Matrix([[B,0,0],[A,Cv,2*Z],[0,0,W]])
    assert s.expand(matrix.det()-B*Cv*W)==0
    save('cell-inverse-symbolic.json',dict(classification='Proven',passed=True,
        ansatz='Within one fixed-metric cell retain positive nonuniform reference rho0(x),T0(x), fixed composition and Q(x). Set rho=exp(eta)*rho0, T=exp(theta)*T0, and a common cell velocity v. This is a specified reconstruction family, not an arbitrary unresolved fluid state.',
        conserved_moments='B=integral[a*D dV_coordinate], U=integral[E dV_coordinate], Pi=integral[a*J dV_coordinate], dV_coordinate=4*pi*r^2 dr. Set kappa=integral[rho0 dV_coordinate]/integral[a*rho0 dV_coordinate]. Since D=exp(eta)*W*rho0, integral[D dV_coordinate]=kappa*B exactly in this family.',
        stable_energy='At fixed geometry use Ubar=U-C*kappa*B=integral[(E-C*D) dV_coordinate], C=C_X*c^2. This preserves the gravitational versus proper-volume weights rather than subtracting C*B from U. When geometry changes, kappa changes and its energy contribution must be included.',
        local_inverse='At eta=theta=v=0, the Jacobian of (B,Ubar,Pi) has rows [B0,0,0], [A,Cv_coordinate,2*integral Q dV_coordinate], [0,0,W_proper], where Cv_coordinate=integral rho0*cv*T0 dV_coordinate>0 and W_proper=integral a*(epsilon+P) dV_coordinate>0. Hence determinant=B0*Cv_coordinate*W_proper>0 for this fixed-metric reconstruction family and a differentiable EOS.',
        limits='No single uniform EOS state replaces the reference. Local inverse theorem requires differentiability and does not certify native EOS, a finite parameter neighbourhood, positive quadrature error, full GR constraints, shocks, arbitrary subcell moments or nonlinear evolution.'))


def retain(i,j,result,left,right,at,rout,width,mout,mass_scale,CX,dm,volume_scale,u_scale,Pmid):
    plan=json.loads((OUT/'plan.json').read_text());c=g.c.gr.C*100;G=g.c.gr.G*1000
    for number in plan['node_counts']:
        nodes,weights=np.polynomial.legendre.leggauss(number);lp=(left+right)/2+(right-left)*nodes/2
        weights=weights*(right-left)/2;geometry=result.sol(lp)
        radius=rout-width*geometry[0];mass=mout-mass_scale*geometry[1]
        values=[];temperatures=[]
        for p in lp:
            a,t=at(float(p));values.append(a);temperatures.append(t)
        values=np.array(values);rho,P,u=values[:,:3].T;f=1-2*mass/radius;aa=1/np.sqrt(f)
        pg=G*P/c**4;eg=G*rho*(CX*c*c+u)/c**4
        dr=-pg*radius*(radius-2*mass)/((eg+pg)*(mass+4*np.pi*radius**3*pg))
        coordinate_weights=-4*np.pi*radius**2*dr*weights;proper_weights=aa*coordinate_weights
        assert np.all(coordinate_weights>0) and np.all(proper_weights>0)
        baryon=np.sum(rho.astype(np.longdouble)*proper_weights)
        mquad=np.sum((rho*(CX*c*c+u)).astype(np.longdouble)*coordinate_weights)*np.longdouble(G/c**4)
        energy=np.sum((rho*u).astype(np.longdouble)*proper_weights)
        final=result.y[:,-1];bm=dm*final[2];mm=mass_scale*final[1];um=dm*u_scale*final[4]
        row=dict(cell=i,tolerance_index=j,nodes=number,baryon_relative_direct_difference=float(baryon/bm-1),
            baryon_relative_inventory_difference=float(baryon/dm-1),mass_relative_direct_difference=float(mquad/mm-1),
            internal_energy_relative_direct_difference=float(energy/um-1))
        row['passed']=bool(abs(row['baryon_relative_direct_difference'])<plan['baryon_mass_relative_tolerance'] and
            abs(row['mass_relative_direct_difference'])<plan['geometric_mass_relative_tolerance'] and
            abs(row['internal_energy_relative_direct_difference'])<plan['internal_energy_relative_tolerance'])
        np.savez_compressed(OUT/f'cell-{i}-tolerance-{j}-nodes-{number}.npz',logP=lp,lnT=np.array(temperatures),
            eos=values,radius_cm=radius,mass_geom_cm=mass,metric_a=aa,coordinate_weights_cm3=coordinate_weights,
            proper_weights_cm3=proper_weights,baryon_g=baryon,mass_geom_integral_cm=mquad,proper_internal_energy_erg=energy,
            C_X=CX,kappa_rest_weight=np.sum(rho.astype(np.longdouble)*coordinate_weights)/baryon)
        records.append(row);save('subcell-quadrature.json',dict(classification='Counterexample candidate',rows=records))
        print('SUBCELL QUADRATURE',row,flush=True)


def run():
    source,changed,_=transformed_source();assert changed==(OUT/'node-run.py').read_text()
    namespace=dict(original.run.__globals__,OUT=OUT,save=save,retain=retain,verify=verify)
    exec(compile(changed,str(OUT/'node-run.py'),'exec'),namespace);namespace['run']()


def verify():
    for rel,digest in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'result.json').read_text())['completed']
    assert json.loads((OUT/'cell-inverse-symbolic.json').read_text())['passed']
    rows=json.loads((OUT/'subcell-quadrature.json').read_text())['rows'];assert len(rows)==16
    print('PASS retained native subcell reference bindings; finite quadrature results only',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
