"""Positive material-coordinate quadratures of every nonuniform GR reference cell."""
import json,sys
from concurrent.futures import ProcessPoolExecutor
from types import FunctionType,SimpleNamespace
import numpy as np
from scipy.integrate import solve_ivp
import sympy as s
import gr_increment_structure as structure

g=structure.g;OUT=g.OUT/'gr-full-subcell-reference'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    paths=[g.ROOT/'verification/gr_full_subcell_reference.py',g.ROOT/'verification/gr_increment_structure.py',
        g.ROOT/'verification/common_eos.py',g.ROOT/'verification/baryon_entropy.py',
        g.ROOT/'verification/audit_structured_enthalpy.py',g.OUT/'initial-adiabats-17.npz',
        g.OUT/'reference-state.npz',structure.OUT/'path-4.npz',structure.OUT/'path-4.json',
        structure.OUT/'plan.json',structure.OUT/'seed-unit-correction.json']
    assert json.loads((structure.OUT/'path-4.json').read_text())['interface_passed']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='b84cb63',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},nodes=[8,16],
        cells=5735,workers=4,block_size=128,pilot_cells=[0,1,2,2971,2972,5734],
        rtol=2e-10,atol=1e-12,endpoint_tolerance=1e-8,volume_relative_tolerance=1e-5,
        shell_mass_relative_tolerance=1e-5,inventory_relative_tolerance=1e-10,
        finite_quadrature_relative_tolerance=1e-5,
        geometry='Reintegrate each saved cell branch with DOP853 on local scaled coordinate increments, using the unchanged 17-point entropy table with strict direct EOS fallback. No full-star refit. Evaluate the native entropy inverse independently at every exported node.',
        quadrature='Use enclosed baryon fraction to the one-third power on the inner branch, and exterior baryon fraction on the outer branch. The full material interval includes the central and surface seeds. Nodes below a seed use that original analytic seed expansion; all others use dense DOP853 output.',
        weights='dB is the original coordinate measure, not a fitted normalization. Positive proper weights are dB/rho_native and coordinate weights dB/(a*rho_native). Baryon recovery is therefore an algebraic coordinate identity, not independent validation. Compare their coordinate volume against the separate radius-face volume and their native energy against independently retained shell mass.',
        seed='Original 1 metre centre seed and original surface first-order seed. Their approximation error is not enclosed. No central cell is omitted.',
        scope='Full-grid finite quadrature/reference test. Native-at-nodes does not make the intervening tabulated geometry a direct-EOS curve. No continuum, physical EOS, complete primitive-root, transport, metric-time or observational certificate.'))
    x,a,b,rho,metric=s.symbols('x a b rho metric',positive=True)
    q=a**3+(b**3-a**3)*x
    assert s.integrate(s.diff(q,x),(x,0,1))==b**3-a**3
    assert s.simplify(rho*metric/(rho*metric)-1)==0
    z=s.symbols('z',positive=True)
    assert s.integrate(3*z*z,(z,a,b))==b**3-a**3
    # Uniform-density regular-centre control: r~q^(1/3), coordinate volume is exact.
    nodes,weights=np.polynomial.legendre.leggauss(8);z=(nodes+1)/2
    assert abs(np.dot(weights/2,3*z*z)-1)<1e-14
    save('symbolic.json',dict(classification='Proven',passed=True,
        identity='For positive rho,a, dVproper=dB/rho and dVcoord=dB/(a*rho). The cube-root inner coordinate removes the uniform-density regular-centre radius singularity. Eight-point Gauss integrates the manufactured flat uniform-density volume exactly to the recorded floating check.',
        limit='The baryon identity alone cannot certify radial geometry or native EOS.'))


def bindings():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    return plan


def block(job):
    label,cells=job;plan=bindings();ld=np.longdouble
    data=dict(np.load(g.OUT/'reference-state.npz'));saved=dict(np.load(structure.OUT/'path-4.npz'))
    solver=structure.Structure(data,4);m=solver.mat;n=len(m.lp);B=ld(m.B)
    stats=dict(calls=0,evaluations=0,maximum_score=0.,label=label)
    inverse=FunctionType(g.audit.strict_invert.__code__,dict(g.audit.strict_invert.__globals__,
        ROOT_STATS=stats,e=SimpleNamespace(OUT=OUT,save=save)))
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


def pilot():
    plan=bindings();record=block(('pilot',plan['pilot_cells']))
    save('pilot-manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p)
        for p in OUT.glob('pilot*') if p.name!='pilot-manifest.json'}))
    assert record['all_passed'],'pilot failure retained; full-grid start is gated'


def run():
    plan=bindings();assert json.loads((OUT/'pilot.json').read_text())['all_passed']
    for rel,digest in json.loads((OUT/'pilot-manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest
    jobs=[(f'block-{start:04}',list(range(start,min(start+plan['block_size'],plan['cells']))))
        for start in range(0,plan['cells'],plan['block_size'])]
    assert not any(OUT.glob('block-*.json')),'do not overwrite a started full-grid run'
    records=[]
    with ProcessPoolExecutor(max_workers=plan['workers']) as pool:
        for record in pool.map(block,jobs):
            records.append(record)
            save('progress.json',dict(classification='Counterexample candidate',completed_cells=sum(len(r['rows']) for r in records)))
    rows=[r for part in records for r in part['rows']]
    save('result.json',dict(classification='Counterexample candidate',completed=True,cells=len(rows),
        all_passed=all(r['passed'] for r in rows),failed_cells=[r for r in rows if not r['passed']],
        maximum_finite_quadrature_difference=max(r['finite_quadrature_relative_difference'] for r in rows),
        full_GR_evolution=False,physical_EOS_certified=False,continuous_errors_certified=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p)
        for p in OUT.iterdir() if p.is_file() and p.name!='manifest.json'}))
    verify()


def verify():
    plan=bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    result=json.loads((OUT/'result.json').read_text());assert result['completed'] and result['cells']==plan['cells']
    cells=[]
    for path in sorted(OUT.glob('block-*.json')):cells.extend(r['cell'] for r in json.loads(path.read_text())['rows'])
    assert cells==list(range(plan['cells']))
    print('PASS complete full-grid subcell bindings; consult failed_cells for scientific gates',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
