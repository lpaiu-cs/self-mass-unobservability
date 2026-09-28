c  Code unform_bash.f
c
c  Used by bash script BASHOP
c  Reads fomratted file fmzz.ttt
c  Write unformatted file mzz.ttt
c
      parameter(nptot=1000000)
      DIMENSION M(NPTOT),Y(NPTOT),FION(-1:30),UMESH(NPTOT)
c      
c Modified by CM on 20-11-04
C
        OPEN(2,FILE='LABO',STATUS='UNKNOWN',FORM='UNFORMATTED')
C
C       Get files mzz.ttt	
           READ(5,100)
     +	   IZZ,ITE,AMAMU,UMIN,UMAX,NCRSE,NTOT,DPACK,JN1,JN2,JN3
           WRITE(2)IZZ,ITE,AMAMU,UMIN,UMAX,NCRSE,NTOT,DPACK,JN1,JN2,JN3
           DO 10 J=JN1,JN2,JN3
             READ(5,101)JNE,EPATOM,OPLNCK,OROSS,NE1,NE2,
     +       (FION(NE),NE=NE1,NE2)
             WRITE(2)JNE,EPATOM,OPLNCK,OROSS,NE1,NE2,
     +       (FION(NE),NE=NE1,NE2)
	     READ(5,102)NP
             WRITE(2)NP
	     IF(NP.GT.0)THEN
	       READ(5,103)(M(N),Y(N),N=1,NP)
	       WRITE(2)(M(N),Y(N),N=1,NP)
	     ELSE
	       READ(5,104)(Y(N),N=1,NTOT)
	       WRITE(2)(Y(N),N=1,NTOT)
	     ENDIF  
   10      CONTINUE
           close(2)
      STOP
C
  100 FORMAT(I3,I4,3E12.4,2I6,E12.4,3I4)
  101 FORMAT(I4,3E12.4,2I3/28E12.4)
  102 FORMAT(I10)
  103 FORMAT(I10,E12.4)
  104 FORMAT(1P,E12.4)
  105 FORMAT(E12.4,i10/(E12.4))
  500 format(A50)
C  
      END
c***********************************************************      
