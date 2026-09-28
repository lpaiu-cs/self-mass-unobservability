c Code monop.f
c
c Writes file of monochromatic opacity, SIGMA_k(u)
c
c Prompts for zz (*2), ttt (*3), jne and name of output file 
c Reads files mzz.ttt, mzz.mesh
c
      parameter(nptot=1000000)
      DIMENSION NX(nptot),Y(nptot),FION(-1:30),umesh(nptot)
      character num(0:9)*1
      data num/'0','1','2','3','4','5','6','7','8','9'/
C
      CHARACTER LAB1*20,LAB2*20,zlab*3,tlab*3
C
c  Read zlab 
      PRINT*,' Read zlab (character*3)'
      READ(5,555)zlab
  555 format(a3)
c
c  Open file .mesh. unit 1
      print*,' Opening '//zlab,'.mesh'
      open(1,file=zlab//'.mesh',status='old',form='unformatted')
      read(1)dv,ntot,(umesh(n),n=1,ntot)
      close(1)
c
c  Read tlab, open file zz.ttt, unit 1
      print*,' Read tlab (character*3)'
      read(5,555)tlab
      print*,' Opening '//zlab//'.'//tlab
      OPEN(1,FILE=zlab//'.'//tlab,
     + STATUS='OLD',FORM='UNFORMATTED')  
      READ(1)IZT,ITT,AMAMUT,UMINT,UMAXT,NCRSET,Ntot2,DPACKT,JN1,JN2,JN3
      if(ntot2.ne.ntot)then
	print*,' ntot=',ntot,', ntot2=',ntot2
	stop
      endif
c
c  Print range of jne, get selected jne
      WRITE(6,620)JN1,JN2,JN3
      READ*,JNE
c
c  Get output file name, open output file, unit 2
      PRINT*,' READ OUTPUT FILENAME, 9 CHARACTERS'
      READ(5,500)LAB2
  500 FORMAT(A20)
      OPEN(2,FILE=LAB2,STATUS='unknown')
      WRITE(2,2000)IZT,ITT,AMAMUT,UMINT,UMAXT,NCRSET,NTOT,DPACKT
C
      DO 10 J=JN1,JN2,JN3
         READ(1)JN,EPATOM,OPLNCK,OROSS,NE1,NE2,(FION(NE),NE=NE1,NE2)
	 read(1)np
	 if(np.gt.0)then
           read(1)(NX(N),Y(N),N=1,NP)
	 else
	   read(1)(y(n),n=1,ntot)
	 endif
         IF(JN.EQ.JNE)THEN
            WRITE(2,2020)EPATOM,OPLNCK,OROSS
            WRITE(2,2040)NE1,NE2
c
c  write fion
	    write(2,2050)(n,fion(n),n=ne1,ne2)
c
c  write umesh and y
	    if(np.gt.0)then
               WRITE(2,2030)(umesh(nx(n)),Y(N),n=1,np)
	    else
	       write(2,2030)(umesh(n),y(n),n=1,ntot)
	    endif
            WRITE(6,600)
            CLOSE(2,STATUS='KEEP')
	    close(1)
	    stop
         ENDIF
   10    CONTINUE
c
c Required file not found
         WRITE(6,610)JNE
         CLOSE(2,STATUS='DELETE')
	 close(1)
	 stop
C
  600 FORMAT(/5X,'FILE WRITTEN'/)
  610 FORMAT(/5X,'JNE=',I3,' NOT FOUND ON INPUT FILE'/)
  620 FORMAT(5X,'JN1,JN2,JN3=',3I5//5X,'SELECT JNE'/)
C
 2000 FORMAT(2I5,1P,3E11.4,0P,2I8,1P,E10.2)
 2010 FORMAT(I3,5X,3I5)
 2020 FORMAT(1P,3E12.3)
 2030 FORMAT(1P,E14.5,E10.2)
 2040 FORMAT(2I5)
 2050 FORMAT(I5,1PE10.2)
C
      END
C
C**********************************************************************
