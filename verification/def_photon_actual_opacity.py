"""Direct-temperature SYNSPEC opacity bridge; population matching is a separate gate.

The pinned MIT upstream stays in the external cache. Every edit and acquired file
is hashed; one-state continuum runs do not certify missing lines or EOS closure.
"""
from pathlib import Path
import argparse
import difflib
import hashlib
import json
import shutil
import subprocess
import time
import urllib.request
import numpy as np
import def_photon_collective as previous

OUT=previous.OUT.parent/'def-photon-actual-opacity'
CACHE=Path('/home/lpaiu/work/direct-eos-gr33/photon-actual-opacity')
COMMIT='b9149f7208eeca9b4fdd38dd11d9f736c7a050d7'
BASE=f'https://raw.githubusercontent.com/callendeprieto/synple/{COMMIT}/'
BUILD=CACHE/'build-balanced'
write=previous.old.ex.write
digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()


def fetch(paths,limit):
    tree={r['path']:r for r in json.loads((OUT/'upstream-tree.json').read_text())['tree']}
    assert sum(tree[p]['size'] for p in paths)<=limit
    start=time.monotonic();rows=[]
    for name in paths:
        path=CACHE/name;path.parent.mkdir(parents=True,exist_ok=True)
        if not path.exists():
            with urllib.request.urlopen(BASE+name,timeout=20) as response:raw=response.read(limit+1)
            assert len(raw)==tree[name]['size']
            assert hashlib.sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest()==tree[name]['sha']
            path.write_bytes(raw)
        raw=path.read_bytes()
        assert hashlib.sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest()==tree[name]['sha']
        rows.append(dict(path=name,bytes=len(raw),sha256=digest(path),git_blob=tree[name]['sha']))
    return dict(files=rows,seconds=time.monotonic()-start)


def build():
    assert not (OUT/(BUILD.name+'.json')).exists()
    sources=['data/'+n for n in ['h1.dat','he1.dat','he2.dat','hydprf.dat','he1prf.dat','he2prf.dat']]
    acquisition=fetch(sources,800000)
    target=BUILD;target.mkdir(exist_ok=True)
    for path in (CACHE/'synspec').glob('*.FOR'):shutil.copyfile(path,target/path.name)
    before=(CACHE/'synspec/PARAMS.FOR').read_text()
    after=before.replace('MTTAB   =     100','MTTAB   =       3').replace('MRTAB   =     100','MRTAB   =       1')
    (target/'PARAMS.FOR').write_text(after)
    patches=list(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='upstream/PARAMS.FOR',tofile='build/PARAMS.FOR'))
    before=(CACHE/'synspec/synspec54.f').read_text();after=before
    # Root guard covers both INGRID density-grid branches, including one density.
    assert after.count('if(ndens.gt.0) dr=(at2-at1)/(ndens-1)')==2
    after=after.replace('if(ndens.gt.0) dr=(at2-at1)/(ndens-1)','if(ndens.gt.1) dr=(at2-at1)/(ndens-1)')
    # The ninth threshold is outside ENEV's eight columns; its INPOT class
    # is above He II. Avoid the upstream out-of-bounds read in all callers.
    after=after.replace('            if(enev(i,j).ge.enhe2) then',
        '            if(j.eq.9) then\n               inpot(i,j)=3\n             else if(enev(i,j).ge.enhe2) then')
    assert after.count('      DO JJ=1,NATOMS')==1
    after=after.replace('      DO JJ=1,NATOMS','      DO JJ=1,MATEX')
    # INTERP's F77 length-one dummy arrays are used with the passed lengths.
    after=after.replace('DIMENSION X(1),Y(1),XX(1),YY(1)',
        'DIMENSION X(NX),Y(NX),XX(NXX),YY(NXX)')
    after=after.replace('      ILLAST=INDLIN(NLIN)',
        '      ILLAST=0\n      IF(NLIN.GT.0) ILLAST=INDLIN(NLIN)')
    # HYDLIN writes levels 1..50; use the occupation-probability capacity.
    after=after.replace('DIMENSION PJ(40),PRF0(54)','DIMENSION PJ(NLMX),PRF0(54)')
    after=after.replace('DIMENSION ABLIN(1),EMLIN(1),OSCHE2',
        'DIMENSION ABLIN(*),EMLIN(*),OSCHE2')
    after=after.replace('dimension wltab(1),absop(1),wlgrid(1),abgrd(1)',
        'dimension wltab(nfr),absop(nfr),wlgrid(nfgrid),abgrd(nfgrid)')
    anchor='''      if(fr0.lt.freqc(ijcon)) then
         ijcon=ijcon+1
         absta=0.5*(absoc(ijcon)+scatc(ijcon)+
     *             absoc(ijcon-1)+scatc(ijcon-1))
      end if'''
    assert after.count(anchor)==1
    after=after.replace(anchor,'''      do while(ijcon.lt.nfreqc)
         if(fr0.ge.freqc(ijcon)) exit
         ijcon=ijcon+1
      end do
      absta=0.5*(absoc(ijcon)+scatc(ijcon)+
     *          absoc(ijcon-1)+scatc(ijcon-1))''')
    assert after.count('IL0=INDLIP(IPRSET-1)+1')==1
    after=after.replace('IL0=INDLIP(IPRSET-1)+1',
        'IF(IPRSET.GT.1) IL0=INDLIP(IPRSET-1)+1')
    # Table-mode continuum interpolation independently interpolates j and chi,
    # breaking Kirchhoff on broad batches. Evaluate CROSS and OPAC at all
    # table frequencies, retaining the existing separate atomic line routines.
    split=after.index('      SUBROUTINE OPACW(')
    head=after[:split];tail=after[split:]
    assert head.count('IF(IMODE.EQ.2) IJ0=NFREQ')==2
    head=head.replace('IF(IMODE.EQ.2) IJ0=NFREQ',
        'IF(IMODE.EQ.2.OR.IMODE0.EQ.-3) IJ0=NFREQ')
    anchor='''      DO IJ=3,NFREQ
         ABSO(IJ)=FRX1(IJ)*ABSO(2)+FRX2(IJ)*ABSO(1)
         EMIS(IJ)=FRX1(IJ)*EMIS(2)+FRX2(IJ)*EMIS(1)
         SCAT(IJ)=FRX1(IJ)*SCAT(2)+FRX2(IJ)*SCAT(1)
      END DO'''
    assert head.count(anchor)==1
    head=head.replace(anchor,'      IF(IMODE0.NE.-3) THEN\n'+anchor+'\n      END IF')
    after=head+tail
    start=after.index('      SUBROUTINE LINOP(')
    end=after.index('      SUBROUTINE LINOPW(')
    line_source=after[start:end]
    anchor='         EMLIN(IJ)=EMLIN(IJ)+ABLIN(IJ)*PLAN(ID)'
    assert line_source.count(anchor)==1
    line_source=line_source.replace(anchor,'''         if(imode0.eq.-3) then
            fr=freq(ij)
            xx=exp(-HK*fr/TEMP(ID))
            ablin(ij)=ablin(ij)*(1.d0-xx)/stim(id)
            bplan=BN*(fr*1.d-15)**3*xx/(1.d0-xx)
            emlin(ij)=emlin(ij)+ablin(ij)*bplan
          else
            EMLIN(IJ)=EMLIN(IJ)+ABLIN(IJ)*PLAN(ID)
         end if''')
    after=after[:start]+line_source+after[end:]
    after=after.replace('637 format(i10,f14.5,0pf12.5)','637 format(i10,1p2e26.17)')
    anchor='         SCAT(IJ)=SCAD+SCLY+sce'
    assert after.count(anchor)==1
    after=after.replace(anchor,anchor+'\nC     Direct audit before spectral interpolation; cgs volume absorption.\n      if(imode0.le.-3) write(29,\'(1p8e26.17)\')\n     *   FR,ABSO(IJ),EMIS(IJ)*X1/(BNU*X),\n     *   ABF,ANE*X*EBF,ANE*X1*AFF,ABAD,SCAD')
    anchor='         POPUL(J,ID)=POPLTE(J)\n      END DO'
    assert after.count(anchor)==1
    after=after.replace(anchor,anchor+'''\n      write(28,'(a,1p3e26.17)') 'STATE',temp(id),dens(id),elec(id)
      do i=1,natoms
         do j=1,mion0
            write(28,'(a,2i5,1p3e26.17)') 'ION',i,j-1,
     *          rrr(id,j,i),pfstd(j,i),rrr(id,j,i)*pfstd(j,i)
         end do
      end do
      do j=1,nlevel
         write(28,'(a,3i5,1p4e26.17)') 'LEVEL',j,iatm(j),iel(j),
     *        popul(j,id),g(j),enion(j),wop(j,id)
      end do
      do i=1,nion
         sn=0.d0
         do j=nfirst(i),nlast(i)
            sn=sn+popul(j,id)
         end do
         write(28,'(a,2i5,1pe26.17)') 'STAGE',
     *        numat(iatm(nfirst(i))),iz(i)-1,sn
      end do
      do i=1,natom
         j=nka(i)
         write(28,'(a,2i5,1pe26.17)') 'STAGE',
     *        numat(i),iz(iel(j)),popul(j,id)
      end do''')
    anchor='      SCAT(2)=SCAT(2)-SCLY\n      RETURN'
    assert after.count(anchor)==1
    after=after.replace(anchor,'''      SCAT(2)=SCAT(2)-SCLY
      if(imode0.le.-3) then
         do ij=3,nfreq-1
            fr=freq(ij)
            bplan=BN*(fr*1.d-15)**3/(exp(HK*fr/T)-1.d0)
            write(88,'(1p3e26.17)') fr,abso(ij),emis(ij)/bplan
         end do
      end if
      RETURN''')
    (target/'synspec54.f').write_text(after)
    patches+=list(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='upstream/synspec54.f',tofile='build/synspec54.f'))
    patch=OUT/(BUILD.name+'.patch');patch.write_text(''.join(patches))
    command=['gfortran','-std=legacy','-fno-automatic','-mcmodel=medium','-g','-fcheck=all','-fbacktrace','synspec54.f','-o','synspec54']
    start=time.monotonic()
    result=subprocess.run(command,cwd=target,capture_output=True,text=True,timeout=90)
    (OUT/(BUILD.name+'.log')).write_text(result.stdout+result.stderr)
    write(OUT/(BUILD.name+'.json'),dict(command=command,returncode=result.returncode,seconds=time.monotonic()-start,
        acquisition=acquisition,patch_sha256=digest(patch),
        binary_sha256=digest(target/'synspec54') if result.returncode==0 else None))
    assert result.returncode==0,result.stderr[-4000:]


def run(name='pilot-final',temperature=None,step=2.,wstart=200.,wend=20000.,metals=False,lines=False,fixed_state=False):
    destination=OUT/name;assert not destination.exists();destination.mkdir()
    work=CACHE/name;assert not work.exists();work.mkdir()
    (work/'data').symlink_to(CACHE/'data',target_is_directory=True)
    if metals:(work/'atom').symlink_to(CACHE/'models/atom',target_is_directory=True)
    if lines:(work/'fort.19').symlink_to(CACHE/'selected-full.19')
    bank=np.load(previous.old.OUT/'bank.npz');inventory=np.load(previous.OUT/'inventory.npz')
    T=float(bank['T']) if temperature is None else temperature
    ne=float(bank['ne'])/1e6
    Z=previous.inventory_reader.g.d.CHARGES
    abundance={int(z):float(n) for z,n in zip(Z,inventory['ni'].sum(axis=1)/inventory['ni'][0].sum()) if n>0}
    # Declare absent elements explicitly: -3 otherwise enables unlisted solar
    # elements above NATOMS, which would silently change the composition.
    atoms=[f'{T:.17g} 4.0','T F',"'tas'",'50','99']
    explicit=[1,2,6,7,8,10,12,20] if metals else [1,2]
    atoms += [f'{2 if z in explicit else 1 if z in abundance else 0} {abundance.get(z,0.):.17g} 0' for z in range(1,100)]
    atoms += ["1 0 9 0 0 0 'H 1' 'data/h1.dat'","1 1 1 1 0 0 'H 2' ' '",
        "2 0 24 0 0 0 'He 1' 'data/he1.dat'","2 1 20 0 0 0 'He 2' 'data/he2.dat'","2 2 1 1 0 0 'He 3' ' '"]
    if metals:
        groups={6:[(0,40,'c1_28+12lev.dat'),(1,22,'c2_17+5lev.dat'),(2,46,'c3_34+12lev.dat')],
            7:[(0,34,'n1_27+7lev.dat'),(1,42,'n2_32+10lev.dat'),(2,32,'n3_25+7lev.dat')],
            8:[(0,33,'o1_23+10lev.dat'),(1,48,'o2_36+12lev.dat'),(2,41,'o3_28+13lev.dat')],
            10:[(0,35,'ne1_23+12lev.dat'),(1,32,'ne2_23+9lev.dat'),(2,34,'ne3_22+12lev.dat')],
            12:[(1,25,'mg2_21+4lev.dat')],20:[(0,66,'Ca1kas_F_zat.sy'),(1,24,'Ca2kas_F_zat.sy')]}
        for z,records in groups.items():
            for q,n,file in records:
                parent='data' if z==20 else 'atom'
                atoms.append(f"{z} {q} {n} 0 0 0 'Z{z:02}' '{parent}/{file}'")
            atoms.append(f"{z} {records[-1][0]+1} 1 1 0 0 'top' ' '")
    atoms.append("0 0 0 -1 0 0 ' ' ' '")
    inputs={'fort.5':'\n'.join(atoms)+'\n','tas':'ND=1\nIFMOL=0\nTMOLIM=10000.\n',
        'fort.2':f"1 {T:.17g} {T:.17g}\n0\n1 {ne:.17g} {ne:.17g}\n10000 1 {wstart} {wend}\n'opacity.txt' 0\n",
        'fort.55':f'{-3 if lines else -4} 0 4\n1 0 0 0\n0 0 0 0 0\n1 0 0 0 0\n0 0 0\n{wstart} {-wend} 200. 0 1e-4 {step}\n0 0\n0.0\n'}
    if fixed_state:
        inputs['fort.2']=f"0 0 0\n0\n1 1 1\n10000 1 {wstart} {wend}\n'opacity.txt' 0\n"
        inputs['fort.8']=f"1 3\n1\n{T:.17g} {ne:.17g} {float(bank['rho']):.17g}\n"
    for filename,content in inputs.items():
        (work/filename).write_text(content);(destination/filename).write_text(content)
    start=time.monotonic()
    with (work/'fort.5').open() as inp,(destination/'stdout.log').open('w') as log:
        result=subprocess.run([str(BUILD/'synspec54')],stdin=inp,stdout=log,stderr=subprocess.STDOUT,cwd=work,timeout=120)
    for filename in ['fort.27','fort.28','fort.29','fort.88','opacity.txt']:
        if (work/filename).exists():shutil.copyfile(work/filename,destination/filename)
    write(destination/'run.json',dict(classification='Counterexample candidate',returncode=result.returncode,
        seconds=time.monotonic()-start,T=T,ne_cm3=ne,EOS_rho=float(bank['rho']),
        explicit_elements=explicit,lines=lines,fixed_rho_ne=fixed_state,
        other_elements='F continuum absent; Mg I/Mg III and Ca III continua absent; higher metal stages are terminal reference ions.',
        mode='Native independent atomic populations; NOT injected EOS levels. Mode -4 still includes native H/He II lines, as established from OPAC source.',
        binary_sha256=digest(BUILD/'synspec54')))
    assert result.returncode==0,(destination/'stdout.log').read_text()[-3000:]
    assert (destination/'fort.27').stat().st_size>0
    raw=np.loadtxt(destination/'fort.88');assert np.all(np.isfinite(raw)) and np.all(raw[:,1]>0)
    print(name, 'seconds',time.monotonic()-start,'samples',len(raw),'native Kirchhoff max',np.max(abs(raw[:,2]/raw[:,1]-1)),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['build','pilot']);args=parser.parse_args()
    if args.action=='build':build()
    else:run()
