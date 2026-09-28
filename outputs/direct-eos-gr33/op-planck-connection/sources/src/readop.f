c  Code readop.f
c
c  For all elements zz, reads files mzz.index
c  For all zz and all ttt, reads files mzz.ttt
c  Write files giving summary information:
c    mzz.smry for electros oer atom, Planck, Ross cross-sections;
c    mzz.ion for ionization equilibrium
c 
      parameter(nptot=1000000)
      DIMENSION M(nptot),Y(nptot),FION(-1:30)
C
      CHARACTER NUM(0:9)*1,ZLAB*4,TELAB*3,lab1*3
      DATA NUM/'0','1','2','3','4','5','6','7','8','9'/
      CHARACTER ZZ(17)*2
      DATA ZZ/'01','02','06','07','08','10','11','12',
     + '13','14','16','18','20','24','25','26','28'/      
C
      DO L=1,17
      ZLAB='m'//ZZ(L)//'.'
C
C
c  Read index file
      PRINT*,' OPENING 8 ',  ZLAB,'index' 
      OPEN(8, FILE= ZLAB//'index',STATUS='OLD')
      READ(8,*)IZ,AMAMU
      READ(8,*)ITE1,ITE2,ITE3
      READ(8,*)UMIN,UMAX
      READ(8,*)NCRSE,NTOT
      READ(8,*)DPACK
      CLOSE(8)
c
c  Initialise .smry file, unit 2
      OPEN(2,FILE=ZLAB//'smry',STATUS='UNKNOWN')
      WRITE(2,2000)IZ,AMAMU,UMIN,UMAX,NCRSE,NTOT,DPACK
      WRITE(2,2010)ITE1,ITE2,ITE3
c
c  Initialise .ion file, unit 10
      OPEN(10,FILE=ZLAB//'ion',STATUS='UNKNOWN')
      WRITE(10,100)ITE1,ITE2,ITE3
C
c  Start I loop
      DO 20 I=ITE1,ITE2,ITE3
         TELAB=NUM(I/100)//NUM(I/10-10*(I/100))//NUM(I-10*(I/10))
c        Open .ttt file, unit 1
         OPEN(1,FILE=ZLAB//TELAB,STATUS='OLD',
     +   FORM='UNFORMATTED')
         READ(1)IZZ,ITE,AMAMU,UMIN,UMAX,NCRSE,NTOT,DPACK,JN1,JN2,JN3
         WRITE(2,2010)ITE,JN1,JN2,JN3
         WRITE(10,100)ITE,JN1,JN2,JN3
         DO 10 J=JN1,JN2,JN3
            READ(1)JNE,EPATOM,OPLNCK,OROSS,NE1,NE2,
     +      (FION(NE),NE=NE1,NE2)
	    read(1)np
	    if(np.gt.0)then
              read(1)(M(N),Y(N),N=1,NP)
	    else
	      read(1)(y(n),n=1,ntot)
	    endif
            WRITE(2,2020)JNE,EPATOM,OPLNCK,OROSS
            WRITE(10,109)JNE,NE1,NE2,(NE,FION(NE),NE=NE1,MIN(NE2,NE1+3))
            IF(NE1+4.LE.NE2)THEN
               WRITE(10,110)(NE,FION(NE),NE=MIN(NE2,NE1+3)+1,NE2)
            ENDIF
   10    CONTINUE
         CLOSE(1)
   20 CONTINUE
      CLOSE(2,STATUS='KEEP')
C
C
      ENDDO
C
C
      STOP
C
  620 FORMAT(//5X,' *** ERROR OPENING FILE ',A11,' ***'//)
C
 2000 FORMAT(I5,1P,3E11.4,0P,2I8,1P,E10.2)
 2010 FORMAT(I3,5X,3I5)
 2020 FORMAT(I5,1P,2E12.3,e13.5)
  100 FORMAT(4I5)
  109 FORMAT(3I4,2X,4(0P,I4,',',1P,E9.3))
  110 FORMAT(14X,
     + 0P,I4,',',1P,E9.3,
     + 0P,I4,',',1P,E9.3,
     + 0P,I4,',',1P,E9.3,
     + 0P,I4,',',1P,E9.3)
C
      END
