"""Measure radial-chart consistency without changing any reference node or weight."""
from pathlib import Path
import json,sys
import numpy as np
import gr_metric_coupled_subcell as metric

ROOT=metric.ROOT;OUT=metric.OUT.parent/'gr-subcell-coordinate-chart';sha=metric.sha;ld=np.longdouble


def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2)+'\n')


def prepare():
    assert not OUT.exists();prior=metric.bindings();OUT.mkdir()
    files=[Path(__file__),metric.OUT/'plan.json']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='0beea1d2',
        bindings=dict(prior['bindings'],**{p.relative_to(ROOT).as_posix():sha(p) for p in files}),
        cells=5735,nodes=[8,16],diagnostic_threshold=1e-5,
        definition='Use the unchanged positive coordinate weights, original radius nodes and material-face radii. Compare A*(weights/w) with 4*pi/3*(r^3-r_inner^3), normalized by the exact endpoint shell volume. Also compare the sum of weights with the endpoint shell volume. Evaluate cube differences in factored form to avoid cancellation.',
        scope='A finite 64-significand-bit diagnostic of the supplied chart and collocation rule. The threshold is a diagnostic scale, not a certified continuous error. No rescaling, EOS calls, state evolution, physical calibration or reference-profile repair.'))


def bindings():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert sha(ROOT/rel)==digest,rel
    return plan


def volume(r,inner):return ld(4)*ld(np.pi)/3*(r-inner)*(r*r+r*inner+inner*inner)


def calculate():
    plan=bindings();reference=json.loads((metric.OUT/'plan.json').read_text())
    grid=np.load(metric.g.OUT/'gr-increment-structure/path-4.npz');faces=grid['radius_m'].astype(ld)*100
    rows=[];controls=[]
    assert np.finfo(ld).nmant+1==64
    for number in plan['nodes']:
        s,w,A,errors=metric.integration_matrix(number)
        # The central affine-radius chart has a quadratic volume Jacobian.
        V=4*ld(np.pi)*s*s;truth=volume(s,ld(0))
        control=float(np.max(abs(A@V-truth))/(4*np.pi/3))
        negative=float(np.max(abs(A@(2*V)-truth))/(4*np.pi/3))
        assert control<1e-13 and negative>plan['diagnostic_threshold']
        controls.append(dict(nodes=number,affine_radius_error=control,altered_weight_defect=negative,
            polynomial_errors=errors,passed=True))
        for rel in sorted(set(reference['cell_sources'].values())):
            source=ROOT/rel;data=np.load(source.with_name(source.stem+f'-nodes-{number}.npz'))
            for j,index in enumerate(data['cells']):
                index=int(index);r=data['radius_cm'][j].astype(ld);order=np.argsort(r);r=r[order]
                weights=data['coordinate_weights_cm3'][j,order].astype(ld)
                inner,outer=faces[index+1],faces[index]
                assert 0<=inner<outer and np.all((r>inner)&(r<outer)) and np.all(weights>0)
                shell=volume(outer,inner);cumulative=A@(weights/w);geometric=volume(r,inner)
                defect=cumulative-geometric
                rows.append(dict(cell=index,nodes=number,
                    maximum_node_volume_defect=float(np.max(abs(defect))/shell),
                    endpoint_volume_defect=float(abs(np.sum(weights)-shell)/shell),
                    minimum_cumulative_volume_fraction=float(np.min(cumulative)/shell),
                    maximum_cumulative_volume_fraction=float(np.max(cumulative)/shell)))
    rows.sort(key=lambda r:(r['cell'],r['nodes']))
    assert [(r['cell'],r['nodes']) for r in rows]==[(i,n) for i in range(5735) for n in [8,16]]
    return dict(classification='Counterexample candidate',completed=True,cells=5735,representations=len(rows),
        controls=controls,rows=rows,summary=[dict(nodes=n,
            maximum_node_volume_defect=max(r['maximum_node_volume_defect'] for r in rows if r['nodes']==n),
            maximum_endpoint_volume_defect=max(r['endpoint_volume_defect'] for r in rows if r['nodes']==n),
            node_diagnostic_exceedances=sum(r['maximum_node_volume_defect']>plan['diagnostic_threshold'] for r in rows if r['nodes']==n),
            endpoint_diagnostic_exceedances=sum(r['endpoint_volume_defect']>plan['diagnostic_threshold'] for r in rows if r['nodes']==n)) for n in [8,16]],
        continuous_chart_error_certified=False,physical_EOS_certified=False,full_GR_evolution=False)


def run():
    assert not (OUT/'result.json').exists();result=calculate();save('result.json',result)
    save('manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in OUT.iterdir() if p.is_file()}));verify()
    print(json.dumps(result['summary']),flush=True)


def verify():
    bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    result=json.loads((OUT/'result.json').read_text());assert result==calculate()
    print('PASS reproducible finite chart diagnostic for 5735 cells at both node counts; inspect defects',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
