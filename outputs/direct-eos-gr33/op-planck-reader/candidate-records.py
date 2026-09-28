def records(path):
    with FortranFile(path,'r',header_dtype='<u4') as f:
        h=f.read_record(np.uint8).tobytes();assert len(h)==44
        z,it,mass,umin,umax,ncoarse,ntot,dpack,jlo,jhi,jstep=struct.unpack('<iifffiifiii',h)
        assert ntot==10000 and ncoarse==int((SOURCE/'mono'/f'm{z:02}.index').read_text().splitlines()[3].split()[0]) and jstep==2
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
