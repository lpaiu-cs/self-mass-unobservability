"""Read actual OP absorption means and audit their new-EOS stellar coverage.

Counterexample candidate. This is a declared OP/FreeEOS hybrid data connection;
unsupported elements, ionization inconsistency and physical errors stay explicit.
"""
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
import hashlib, json, shutil, struct, subprocess, sys
import numpy as np
from scipy.io import FortranFile, FortranEOFError
import gr_molecular_reference_runner as reference

g=reference.g;OUT=g.OUT/'op-planck-connection';CACHE=g.CACHE/'opacity-project'
SOURCE=CACHE/'source';ARCHIVE=CACHE/'OP4STARS_1.3.tar.xz';NATIVE=CACHE/'readop-control'
ELEMENTS=[1,2,6,7,8,10,11,12,13,14,16,18,20,24,25,26,28]


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists() and not NATIVE.exists();OUT.mkdir();reference.verify()
    metadata=json.loads((g.ROOT/'outputs/op-monochromatic-4390522-metadata.json').read_text())
    file=metadata['files'][0];assert file['key']==ARCHIVE.name and ARCHIVE.stat().st_size==file['size']
    assert hashlib.md5(ARCHIVE.read_bytes()).hexdigest()==file['checksum'].split(':')[1]
    assert g.c.sha(ARCHIVE)=='aeae2b31e62c7cebc100be2813e9b976de0681a31faa5fa2716c766cc1a6e809'
    for src,name in [(g.ROOT/'outputs/op-monochromatic-4390522-metadata.json','zenodo-metadata.json'),
                     (CACHE/'extracted.json','source-files.json')]:shutil.copy2(src,OUT/name)
    selected=['README','README_MESA','README_2.1','CHANGES_HU','manuals/OPCD.pdf',
        'manuals/OPserver.pdf','src/readop.f','src/monop.f','src/unform_bash.f',
        'opserver/opserver_dp.f','opserver/opax.f']
    selected += [f'mono/m{z:02}.{suffix}' for z in ELEMENTS for suffix in ['index','smry']]
    for rel in selected:
        target=OUT/'sources'/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(SOURCE/rel,target)
    constants=reference.original.model.OUT.parent/'gr-radiation-eos-split/mod_free_eos_constants.f90'
    assert 'avogadro = 6.02214076e23_fp_kind' in constants.read_text() and g.c.NA==6.02214076e23
    paths=[g.ROOT/'verification/op_planck_connection.py',reference.OUT/'manifest.json',constants,
           OUT/'zenodo-metadata.json',OUT/'source-files.json']+list((OUT/'sources').rglob('*'))
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='9115c50',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths if p.is_file()},
        archive=dict(path=str(ARCHIVE),sha256=g.c.sha(ARCHIVE),metadata_checksum=file['checksum']),
        imported_source=dict(classification='Imported from prior work',
            dataset='https://zenodo.org/records/4390522',doi='10.5281/zenodo.4390522',
            manual='OPCD manual sections 4.4.2, 4.7 and 5.2: OPLNCK is the Planck mean absorption cross section excluding scattering; EPATOM is free electrons per atom. Source ion index NE counts bound electrons minus one, so charge=Z-NE-1.'),
        controls='Compile the unchanged archived readop.f in a private directory, symlink input index/numeric files only, and regenerate its summary outputs. Independently parse every native binary record and require each header to lie within the exact last-decimal printing interval of the independent Fortran output. Check every record length, density index, ion support and frequency array; preserve negative-data diagnostics without clipping.',
        target='All 5735 saved new-molecular-EOS rho_B,T,X states, with Ne=N_A*rmue from their actual native output. Do not fit composition, density or an effective electron abundance.',
        interpolation='Positive bilinear interpolation of log(Planck absorption cross section) on each four-corner (logT,logNe) rectangle. Interpolate EPATOM linearly with the same nonnegative weights. Reject missing corners; do not extrapolate.',
        normalization='Partial supported-element contribution per baryon gram: N_A*(10^-16.55280 cm^2)*sum_z Y_z*sigma_P,z. The area conversion follows the archived OP code. Y_z=sum_{isotopes of z} X_i/A_i is not renormalized after unsupported elements are removed.',
        boundary='OP uses its own element-group populations at T,Ne. Compare its predicted supported free-electron density with FreeEOS; this difference is not automatically attributable to unsupported elements. Missing Li,Be,B,F and molecular/isotope effects are not certified small. A Planck mean can supply an LTE emission integral for its declared model; it is not the non-LTE radiation-weighted absorption mean or a frequency-resolved scattering kernel. No total stellar opacity or radiation evolution is claimed.'))


def bindings():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert g.c.sha(ARCHIVE)==plan['archive']['sha256']
    return plan


def native_control():
    NATIVE.mkdir()
    # ponytail: use the original data-reader program, not a new native bridge.
    for z in ELEMENTS:
        for path in (SOURCE/'mono').glob(f'm{z:02}.*'):
            if path.suffix[1:].isdigit() or path.suffix=='.index':(NATIVE/path.name).symlink_to(path)
    exe=NATIVE/'readop';cmd=['gfortran','-std=legacy','-O2',str(SOURCE/'src/readop.f'),'-o',str(exe)]
    build=subprocess.run(cmd,capture_output=True,text=True);(OUT/'build.log').write_text(build.stdout+build.stderr)
    assert build.returncode==0,build.stderr
    run=subprocess.run([str(exe)],cwd=NATIVE,capture_output=True,text=True)
    (OUT/'native-readop.log').write_text(run.stdout+run.stderr);assert run.returncode==0,run.stderr
    save('native-build.json',dict(classification='Counterexample candidate',command=cmd,
        executable=str(exe),sha256=g.c.sha(exe),source_sha256=g.c.sha(SOURCE/'src/readop.f')))
    for z in ELEMENTS:shutil.copy2(NATIVE/f'm{z:02}.smry',OUT/f'native-m{z:02}.smry')


def printed_intervals(path):
    lines=iter(path.read_text().splitlines());next(lines);lo,hi,step=map(int,next(lines).split());result={}
    for expected_t in range(lo,hi+1,step):
        it,jlo,jhi,jstep=map(int,next(lines).split());assert it==expected_t
        for expected_j in range(jlo,jhi+1,jstep):
            row=next(lines).split();j=int(row[0]);assert j==expected_j
            result[it,j]=[Decimal(x) for x in row[1:]]
    assert next(lines,None) is None
    return result


def printing_score(value,printed):
    half=Fraction(10)**printed.as_tuple().exponent/2
    return abs(Fraction(float(value))-Fraction(printed))/half


def records(path):
    with FortranFile(path,'r',header_dtype='<u4') as f:
        h=f.read_record(np.uint8).tobytes();assert len(h)==44
        z,it,mass,umin,umax,ncoarse,ntot,dpack,jlo,jhi,jstep=struct.unpack('<iifffiifiii',h)
        assert ntot==10000 and ncoarse==1 and jstep==2
        for expected_j in range(jlo,jhi+1,jstep):
            row=f.read_record(np.uint8).tobytes();j,epa,planck,ross,nlo,nhi=struct.unpack('<ifffii',row[:24])
            assert j==expected_j and -1<=nlo<=nhi<=z and len(row)==24+4*(nhi-nlo+1)
            ions=np.frombuffer(row[24:],dtype='<f4').astype(float)
            n=int(f.read_ints('<i4')[0]);assert 0<=n<=ntot
            spectrum=f.read_record(np.uint8).tobytes()
            if n:
                packed=np.frombuffer(spectrum,dtype=[('index','<i4'),('sigma','<f4')])
                assert len(packed)==n and packed['index'][0]==1 and packed['index'][-1]==ntot
                assert np.all(np.diff(packed['index'])>0);values=packed['sigma']
            else:
                values=np.frombuffer(spectrum,dtype='<f4');assert len(values)==ntot
            assert np.all(np.isfinite(values)) and np.all(np.isfinite(ions))
            assert np.all(np.isfinite([epa,planck,ross]))
            yield [z,it,j,epa,planck,ross,ions.sum(),ions@np.arange(z-nlo-1,z-nhi-2,-1),
                   float(values.min()),n,mass,umin,umax,dpack]
        try:f.read_record(np.uint8)
        except FortranEOFError:pass
        else:raise AssertionError(('Unexpected extra record',path))


def headers():
    rows=[];worst=Fraction(0);negative=[]
    for z in ELEMENTS:
        printed=printed_intervals(OUT/f'native-m{z:02}.smry');element=[]
        for path in sorted((SOURCE/'mono').glob(f'm{z:02}.[0-9][0-9][0-9]')):
            assert str(path.relative_to(SOURCE)) in source_digests
            for row in records(path):
                assert row[0]==z;key=int(row[1]),int(row[2]);tokens=printed.pop(key)
                for value,token in zip(row[3:6],tokens):
                    score=printing_score(value,token);assert score<=1,(path,key,value,token,float(score));worst=max(worst,score)
                if min(row[3:6])<0 or row[8]<0:negative.append(row[:6]+[row[8]])
                element.append(row)
        assert not printed;rows.extend(element)
        print('OP HEADERS',z,len(element),flush=True)
    values=np.array(rows);np.savez_compressed(OUT/'binary-headers.npz',values=values)
    save('header-control.json',dict(classification='Counterexample candidate',passed=True,rows=len(rows),
        maximum_exact_printing_score=float(worst),negative_rows=negative,
        columns=['Z','ITE','JNE','free_electrons_per_atom','absorption_Planck_cross_section_a0_squared',
            'Rosseland_cross_section_a0_squared','saved_ion_fraction_sum','saved_mean_ion_charge',
            'minimum_saved_monochromatic_cross_section','packed_points','native_mean_atomic_mass',
            'u_min','u_max','packing_tolerance'],
        scope='Independent native Fortran/Python reading agreement. Printed header accuracy, spectral packing accuracy and physical opacity accuracy are distinct.'))
    return values


def connect(values):
    data=dict(np.load(reference.OUT/'reference-state.npz'));n=len(data['X'])
    eos=np.concatenate([np.load(reference.OUT/f'block-{i:04}.npz')['new'] for i in range(0,n,128)])
    ne=g.c.NA*eos[:,13];logT=data['lnT']/np.log(10);logne=np.log10(ne)
    abundance=data['X']/g.c.A;by_z={int(z):abundance[:,g.c.Z==z].sum(axis=1) for z in np.unique(g.c.Z)}
    unsupported=sorted(set(by_z)-set(ELEMENTS));supported=sorted(set(by_z)&set(ELEMENTS))
    grid={(int(row[0]),int(row[1]),int(row[2])):row[3:6] for row in values}
    results=np.full((n,2),np.nan);missing=[]
    for i in range(n):
        it=2*int(np.floor(logT[i]*20));jn=2*int(np.floor(logne[i]*2))
        u=(40*logT[i]-it)/2;v=(4*logne[i]-jn)/2
        weights=np.array([(1-u)*(1-v),(1-u)*v,u*(1-v),u*v]);assert min(weights)>=0 and abs(weights.sum()-1)<1e-14
        sums=np.zeros(2);bad=[]
        for z in supported:
            keys=[(z,it+a,jn+b) for a,b in [(0,0),(0,2),(2,0),(2,2)]]
            if any(key not in grid for key in keys):bad.append(z);continue
            corners=np.array([grid[key] for key in keys]);assert np.all(corners[:,1]>0)
            sums += by_z[z][i]*np.array([np.exp(weights@np.log(corners[:,1])),weights@corners[:,0]])
        if bad:missing.append(dict(cell=i,elements=bad,logT=float(logT[i]),logNe=float(logne[i]),ITE=it,JNE=jn))
        else:results[i]=[g.c.NA*10**-16.55280*sums[0],np.exp(data['lnd'][i])*g.c.NA*sums[1]]
    good=np.isfinite(results[:,0]);assert good.any()
    missing_X=data['X'][:,np.isin(g.c.Z,unsupported)].sum(axis=1)
    electron_gap=results[:,1]/ne-1
    np.savez_compressed(OUT/'stellar-connection.npz',partial_Planck_opacity_per_baryon_gram=results[:,0],
        OP_supported_predicted_electron_density=results[:,1],FreeEOS_electron_density=ne,
        supported_grid_covered=good,unsupported_mass_fraction=missing_X,electron_relative_difference=electron_gap,
        dm=data['dm'],X=data['X'],logT=logT,logNe=logne)
    save('result.json',dict(classification='Counterexample candidate',completed=True,stellar_cells=n,
        supported_elements=supported,unsupported_elements=unsupported,
        unsupported_species=[name for name,z in zip(g.c.NAMES,g.c.Z) if z in unsupported],
        supported_grid_covered_cells=int(good.sum()),missing_grid=missing,
        mass_fraction_outside_supported_grid=float(data['dm'][~good].sum()/data['dm'].sum()),
        maximum_unsupported_mass_fraction=float(missing_X.max()),
        global_unsupported_mass_fraction=float(data['dm']@missing_X/data['dm'].sum()),
        partial_Planck_opacity_range=[float(np.nanmin(results[:,0])),float(np.nanmax(results[:,0]))],
        supported_OP_vs_FreeEOS_electron_difference_range=[float(np.nanmin(electron_gap)),float(np.nanmax(electron_gap))],
        recovered_absorption_Planck_data=True,complete_mixture_absorption=False,
        physical_EOS_or_opacity_certified=False,nonLTE_absorption_known=False,
        scattering_kernel_identified=False,radiation_evolved=False,full_GR_evolution=False))


def run():
    bindings();global source_digests
    source_digests={row['path']:row['sha256'] for row in json.loads((OUT/'source-files.json').read_text())}
    for rel,digest in source_digests.items():assert g.c.sha(SOURCE/rel)==digest,rel
    native_control();values=headers();connect(values)
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.rglob('*') if p.is_file()}))
    verify()


def verify():
    bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    native=json.loads((OUT/'native-build.json').read_text());assert g.c.sha(native['executable'])==native['sha256']
    assert json.loads((OUT/'header-control.json').read_text())['passed']
    assert json.loads((OUT/'result.json').read_text())['completed']
    print('PASS sourced OP Planck data connection; inspect domain/composition/EOS gaps',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
