      FUNCTION GAUNT(I,FR)
C     ====================
C
C     Hydrogenic bound-free Gaunt factor for the principal quantum
C     number I and frequency FR
C
      INCLUDE 'PARAMS.FOR'
      X=FR/2.99793E14
      GAUNT=1.
      IF(I.EQ.1) THEN
      GAUNT=1.2302628+X*(-2.9094219E-3+X*(7.3993579E-6-8.7356966E-9*X))
     *+(12.803223/X-5.5759888)/X
       ELSE IF(I.EQ.2) THEN
      GAUNT=1.1595421+X*(-2.0735860E-3+2.7033384E-6*X)+(-1.2709045+
     *(-2.0244141/X+2.1325684)/X)/X
       ELSE IF(I.EQ.3) THEN
      GAUNT=1.1450949+X*(-1.9366592E-3+2.3572356E-6*X)+(-0.55936432+
     *(-0.23387146/X+0.52471924)/X)/X
       ELSE IF(I.EQ.4) THEN
      GAUNT=1.1306695+X*(-1.3482273E-3+X*(-4.6949424E-6+2.3548636E-8*X))
     *+(-0.31190730+(0.19683564-5.4418565E-2/X)/X)/X
       ELSE IF(I.EQ.5) THEN
      GAUNT=1.1190904+X*(-1.0401085E-3+X*(-6.9943488E-6+2.8496742E-8*X))
     *+(-0.16051018+(5.5545091E-2-8.9182854E-3/X)/X)/X
       ELSE IF(I.EQ.6) THEN
      GAUNT=1.1168376+X*(-8.9466573E-4+X*(-8.8393133E-6+3.4696768E-8*X))
     *+(-0.13075417+(4.1921183E-2-5.5303574E-3/X)/X)/X
       ELSE IF(I.EQ.7) THEN
      GAUNT=1.1128632+X*(-7.4833260E-4+X*(-1.0244504E-5+3.8595771E-8*X))
     *+(-9.5441161E-2+(2.3350812E-2-2.2752881E-3/X)/X)/X
       ELSE IF(I.EQ.8) THEN
      GAUNT=1.1093137+X*(-6.2619148E-4+X*(-1.1342068E-5+4.1477731E-8*X))
     *+(-7.1010560E-2+(1.3298411E-2 -9.7200274E-4/X)/X)/X
       ELSE IF(I.EQ.9) THEN
      GAUNT=1.1078717+X*(-5.4837392E-4+X*(-1.2157943E-5+4.3796716E-8*X))
     *+(-5.6046560E-2+(8.5139736E-3-4.9576163E-4/X)/X)/X
       ELSE IF(I.EQ.10) THEN
      GAUNT=1.1052734+X*(-4.4341570E-4+X*(-1.3235905E-5+4.7003140E-8*X))
     *+(-4.7326370E-2+(6.1516856E-3-2.9467046E-4/X)/X)/X
      END IF
      RETURN
      END
C
C
C     ****************************************************************
C
C
 
      FUNCTION GFREE(T,FR)
C     ====================
C
C     Hydrogenic free-free Gaunt factor, for temperature T and
C     frequency FR
C
      INCLUDE 'PARAMS.FOR'
      THET=5040.4/T
      IF(THET.LT.4.E-2) THET=4.E-2
      X=FR/2.99793E14
      IF(X.GT.1) GO TO 10
      IF(X.LT.0.2) X=0.2
      GFREE=(1.0823+2.98E-2/THET)+(6.7E-3+1.12E-2/THET)/X
      RETURN
   10 C1=(3.9999187E-3-7.8622889E-5/THET)/THET+1.070192
      C2=(6.4628601E-2-6.1953813E-4/THET)/THET+2.6061249E-1
      C3=(1.3983474E-5/THET+3.7542343E-2)/THET+5.7917786E-1
      C4=3.4169006E-1+1.1852264E-2/THET
      GFREE=((C4/X-C3)/X+C2)/X+C1
      RETURN
      END
C
C ********************************************************************
C ********************************************************************
C
      SUBROUTINE STARK0(I,J,IZZ,XKIJ,WL0,FIJ,FIJ0)
C
C     Auxiliary procedure for evaluating the approximate Stark profile
C     of hydrogen lines - sets up necessary frequency independent
C     parameters
C
C     Input:  I     - principal quantum number of the lower level
C             J     - principal quantum number of the upper level
C             IZZ   - ionic charge (IZZ=1 for hydrogen, etc.)
C     Output: XKIJ  - coefficients K(i,j) for the Hotzmark profile;
C                     exact up to j=6, asymptotic for higher j
C             WL0   - wavelength of the line i-j
C             FIJ   - Stark f-value for the line i-j
C             FIJ0  - f-value for the undisplaced component of the line
C
C
      INCLUDE 'PARAMS.FOR'
      PARAMETER (RYD1=911.763811,RYD2=911.495745,CXKIJ=5.5E-5)
      PARAMETER (WI1=911.753578, WI2=227.837832)
      PARAMETER (UN=1.,TEN=10.,TWEN=20.,HUND=100.)
      DIMENSION FSTARK(10,4),XKIJT(5,4),FOSC0(10,4),FADD(5,5)
      DATA XKIJT/3.56E-4,5.23E-4,1.09E-3,1.49E-3,2.25E-3,.0125,.0177,
     * .028,.0348,.0493,.124,.171,.223,.261,.342,.683,.866,1.02,1.19,
     * 1.46/
      DATA FSTARK/  .1387,    .0791,   .02126,   .01394,   .00642,
     *           4.814E-3, 2.779E-3, 2.216E-3, 1.443E-3, 1.201E-3,
     *              .3921,    .1193,   .03766,   .02209,   .01139,
     *           8.036E-3, 5.007E-3,  3.85E-3, 2.658E-3, 2.151E-3,
     *              .6103,    .1506,   .04931,   .02768,   .01485,
     *             .01023, 6.588E-3, 4.996E-3, 3.524E-3, 2.838E-3,
     *              .8163,    .1788,   .05985,   .03189,   .01762,
     *             .01196, 7.825E-3, 5.882E-3, 4.233E-3, 3.375E-3/
      DATA FOSC0 / 0.27746,  0., 0.00773,  0., 0.00134, 0., 
     *             0.000404, 0., 0.000162, 0.,
     *             0.24869,  0., 0.00701,  0., 0.00131, 0.,
     *             0.000422, 0., 0.000177, 0.,
     *             0.23175,  0., 0.00653,  0., 0.00118, 0.,
     *             0.000392, 0., 0.000169, 0.,
     *             0.22148,  0.0005, 0.00563, 0.0004, 0.00108, 0.,
     *             0.000362, 0., 0.000159, 0./ 
      DATA FADD /  1.231, 0.2069, 7.448E-2, 3.645E-2, 2.104E-2,
     *             1.424, 0.2340, 8.315E-2, 4.038E-2, 2.320E-2,
     *             1.616, 0.2609, 9.163E-2, 4.416E-2, 2.525E-2,
     *             1.807, 0.2876, 1.000E-1, 4.787E-2, 2.724E-2,
     *             1.999, 0.3143, 1.083E-1, 5.152E-2, 2.918E-2/
C
      II=I*I
      JJ=J*J
      JMIN=J-I
      IF(JMIN.LE.5.and.i.le.4) THEN
         XKIJ=XKIJT(JMIN,I)
       ELSE
         XKIJ=CXKIJ*(II*JJ)*(II*JJ)/(JJ-II)
      END IF
      IF(I.LE.4) THEN
         IF(JMIN.LE.10) THEN
            FIJ=FSTARK(JMIN,I)
            FIJ0=FOSC0(JMIN,I)
          ELSE 
            CFIJ=((TWEN*I+HUND)*J/(I+TEN)/(JJ-II))
            FIJ=FSTARK(10,I)*CFIJ*CFIJ*CFIJ
            FIJ0=0.
         END IF
       ELSE IF(I.LE.9) THEN
         IF(JMIN.LE.5) THEN
            FIJ=FADD(JMIN,I-4)
            FIJ0=0.
          ELSE 
            CFIJ=((TEN*I+25.)*J/(I+5.)/(JJ-II))
            FIJ=FADD(5,I-4)*CFIJ*CFIJ*CFIJ
            FIJ0=0.
         END IF
       ELSE
         CFIJ=UN*J/(JJ-II)
         FIJ=1.96*I*CFIJ*CFIJ*CFIJ
         FIJ0=0.
      END IF
C
C     wavelength with an explicit correction to the air wavalength
C
      w0=wi1
      if(izz.eq.2) w0=wi2
      WL0=W0/(UN/II-UN/JJ)
      IF(WL0.GT.vaclim) THEN
         ALM=1.E8/(WL0*WL0) 
         XN1=64.328+29498.1/(146.-ALM)+255.4/(41.-ALM)
         WL0=WL0/(XN1*1.D-6+UN)
      END IF        
      RETURN
      END
C
C ********************************************************************
C
