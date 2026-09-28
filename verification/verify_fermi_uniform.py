"""Independent rational audit of the exported Fermi enclosures and bindings."""
from fractions import Fraction as F
from pathlib import Path
import json, sys
import direct_eos_gr as g

OUT=g.OUT/'fermi-uniform-exact';OLD=g.OUT/'fermi-uniform'


def read(path): return json.loads(path.read_text())


def interval(text):
    assert text.startswith('[') and text.endswith(']')
    a,b=map(F,(s.strip() for s in text[1:-1].split(',')));assert a<=b
    return a,b


def export():
    """Include endpoint serialization in each exported discrepancy bound."""
    target=OUT/'exported-records.json';assert not target.exists()
    rows=[];failed=0;maximum_padding=F(0)
    for i in range(3):
        for j,row in enumerate(read(OUT/f'point-{i}.json')['rows']):
            lower,upper=interval(row['enclosure']);native=F.from_float(row['native_value'])
            original=interval(row['native_absolute_error_upper'])[1]
            exact=max(abs(native-lower),abs(native-upper));bound=max(original,exact)
            padding=bound-original;failed+=int(padding>0);maximum_padding=max(maximum_padding,padding)
            rows.append(dict(point=i,component=j,enclosure=row['enclosure'],native_value=row['native_value'],
                native_error_upper_rational=str(bound),serialization_padding_rational=str(padding)))
    assert failed>0,'The preserved failed audit must reproduce'
    result=dict(classification='Proven',rows=rows,raw_record_consistency_passed=False,
        raw_inconsistent_components=failed,maximum_padding_rational=str(maximum_padding),
        display_only_maximum_padding=float(maximum_padding),
        sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in
            [OUT/f'point-{i}.json' for i in range(3)]+[OUT/'before-export-audit.py',OUT/'before-export-audit.log']},
        scope='The raw live-interval discrepancy omitted outward serialization enlargement. These exact rational bounds cover the saved enclosure endpoints and retain every original bound and scientific threshold. Raw records and failed audit remain unchanged.')
    target.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print('EXPORTED FERMI RECORDS',len(rows),'raw inconsistent components',failed,
        'maximum serialization padding',float(maximum_padding),flush=True)


def audit():
    plan=read(OUT/'plan.json');old_plan=read(OLD/'plan.json')
    for rel,digest in plan['bindings'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    for key in ['eta_interval','beta_interval','cases','transformed_interval','coarse_intervals',
            'fine_panels','maximum_uniform_quadrature_error','point_controls',
            'native_point_error_relative_with_floor_one']:
        assert plan[key]==old_plan[key],key
    for name in ['plan','certificate']:
        path=OLD/('plan.json' if name=='plan' else 'uniform-certificate.json')
        assert g.c.sha(path)==plan['preserved_original'][name+'_sha256']
    for name,digest in plan['preserved_original']['point_sha256'].items(): assert g.c.sha(OLD/name)==digest
    bridge=read(OUT/'bridge.json');library=g.d.CACHE/'full-integral-build/src/libfree_eos_direct24_integral_full.so.1.0.0'
    assert g.c.sha(library)==bridge['library_sha256']
    assert g.c.sha(g.CACHE/'fermi-uniform-exact/primitive.so')==bridge['bridge_sha256']
    cert=read(OUT/'uniform-certificate.json');assert cert['passed']
    assert cert['plan_sha256']==g.c.sha(OUT/'plan.json') and len(cert['rows'])==28
    for case,row in zip(plan['cases'],cert['rows']):
        assert all(row[k]==v for k,v in case.items())
        assert interval(row['uniform_error_upper'])[1]<=F(plan['maximum_uniform_quadrature_error'])
    exported=read(OUT/'exported-records.json');assert len(exported['rows'])==84
    for rel,digest in exported['sha256'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    count=0;max_error=F(0);max_scaled=F(0);max_radius=F(0);legacy_contained=0;raw_failed=0
    for i,point in enumerate(plan['point_controls']):
        result=read(OUT/f'point-{i}.json');old=read(OLD/f'point-{i}.json')
        assert result['passed'] and result['point']==point and len(result['rows'])==28
        for j,(case,row,previous) in enumerate(zip(plan['cases'],result['rows'],old['rows'])):
            assert all(row[k]==v for k,v in case.items())
            assert row['native_value']==previous['native_value']
            lower,upper=interval(row['enclosure']);native=F.from_float(row['native_value'])
            saved=exported['rows'][count]
            assert saved['point']==i and saved['component']==j
            assert saved['enclosure']==row['enclosure'] and saved['native_value']==row['native_value']
            original=interval(row['native_absolute_error_upper'])[1]
            bound=F(saved['native_error_upper_rational']);assert bound>=original
            exact=max(abs(native-lower),abs(native-upper));assert exact<=bound
            assert bound==max(exact,original)
            assert F(saved['serialization_padding_rational'])==bound-original
            raw_failed+=int(exact>original)
            scale=max(F(1),abs(native));assert bound<=F(plan['native_point_error_relative_with_floor_one'])*scale
            radius=(upper-lower)/2;assert radius<=F('1e-11')+F('1e-40')
            a,b=interval(previous['enclosure']);legacy_contained+=int(lower<=a and upper>=b)
            max_error=max(max_error,bound);max_scaled=max(max_scaled,bound/scale);max_radius=max(max_radius,radius);count+=1
    assert raw_failed==exported['raw_inconsistent_components'] and raw_failed>0
    result=dict(classification='Proven',passed=True,primitive_components=28,native_point_components=count,
        raw_record_consistency_passed=False,raw_inconsistent_components=raw_failed,
        corrected_exported_records_passed=True,scientific_thresholds_unchanged=True,
        exported_records_sha256=g.c.sha(OUT/'exported-records.json'),
        prior_native_values_bitwise_equal=True,legacy_decimal_intervals_contained=legacy_contained,
        maximum_point_enclosure_radius_rational=str(max_radius),
        maximum_native_error_upper_rational=str(max_error),maximum_native_scaled_error_upper_rational=str(max_scaled),
        display_only=dict(maximum_radius=float(max_radius),maximum_native_absolute_error=float(max_error),maximum_native_scaled_error=float(max_scaled)),
        old_native_uniform_error_certified=False,full_EOS_certified=False,
        scope='Exact rational validation of exported decimals, independent native value identity, protocol/runtime hashes and preserved original records. The three native comparisons include decimal-to-binary64 input conversion and do not establish a uniform native error bound.')
    (OUT/'audit.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print('PASS EXACT FERMI RECORDS',count,'native components; scaled bound',float(max_scaled),'radius',float(max_radius),flush=True)


def freeze():
    assert read(OUT/'audit.json')['passed'] and not (OUT/'manifest.json').exists()
    paths=[p for folder in [OUT,OLD] for p in folder.rglob('*') if p.is_file()]
    paths += [g.ROOT/'verification'/(name+'.py') for name in ['fermi_uniform','fermi_uniform_exact','interval_records','verify_fermi_uniform']]
    result=dict(classification='Proven',scope='Continuous-domain declared Fermi quadrature and finite native controls only.',
        sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in sorted(paths)},
        full_EOS_certified=False,full_GR_evolution=False)
    (OUT/'manifest.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    verify()


def verify():
    manifest=read(OUT/'manifest.json')
    for rel,digest in manifest['sha256'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    print('PASS FROZEN FERMI CERTIFICATE',len(manifest['sha256']),'SHA bindings',flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
