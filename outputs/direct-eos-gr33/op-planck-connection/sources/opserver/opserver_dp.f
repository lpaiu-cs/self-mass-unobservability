c FORTRAN 90 module for calculation of radiative accelerations,
c based on the Opacity Project (OP) code "OPserver".
c See CHANGES_HU for changes made to the original code.
c
c Haili Hu 2010
c
      module opserver_dp
c
c     integer, parameter :: iprec = kind(0.0)
      integer, parameter, private :: iprec = kind(0.0d0)
      integer, parameter, private :: nptot = 10000
      integer, parameter, private :: ipe = 17
      logical, private :: screening = .true.   ! if true, use screening corrections
c      
      contains
C******************************************************************
      subroutine op_dload(ierr)
      parameter(ipz=28)
      common /mesh/ ntotv,dv,dv1,umesh(nptot)
      common/atomdata/ ite1,ite2,ite3,jn1(91),jn2(91),jne3,umin,umax,ntot,
     +  nc,nf,int(17),epatom(17,91,25),oplnck(17,91,25),ne1(17,91,25),
     +  ne2(17,91,25),fion(-1:28,28,91,25),np(17,91,25),kp1(17,91,25),
     +  kp2(17,91,25),kp3(17,91,25),npp(17,91,25),mx(33417000),
     +  yy1(33417000),yy2(120000000),nx(19305000),yx(19305000)
      dimension ifl(ipe),iflp(ipe)
      character num(0:9)*1,zlab(ipe)*3,tlab*6,zlabp(ipe)*3, path*500
      DATA NUM/'0','1','2','3','4','5','6','7','8','9'/
c HH: module opm not needed, few 
      integer :: kz(17)
      data kz/1, 2, 6, 7, 8, 10, 11, 12, 13, 14, 16, 18, 20, 24, 25, 26, 28/      
      save /atomdata/, /mesh/    !HH: put common block in static memory
c      
      ierr=0
      ite1=140
      ite2=320
      ite3=2
      do n=1,ipe
        ifl(n)=50+n
        zlab(n)='m'//num(kz(n)/10)//num(kz(n)-10*(kz(n)/10))
        iflp(n)=70+n
        zlabp(n)='a'//num(kz(n)/10)//num(kz(n)-10*(kz(n)/10))
      enddo
      
      call getenv("oppath", path)
      if (len(trim(path)) == 0) then
         write(6,*) 'Define environmental variable oppath (directory of OP data)' 
         stop
      endif

C  READ INDEX FILES
C     FIRST FILE
      NN=1
c     print*,' Opening '//'./mono/'//zlab(1)//'.index'
      OPEN(1,FILE=trim(path)//'/mono/'//ZLAB(1)//'.index',STATUS='OLD',
     + ERR=1010)
      READ(1,*)IZZ,AMM
      READ(1,*)ITTE1,ITTE2,ITTE3
      READ(1,*)UMIN,UMAX
      READ(1,*)NC,NF
      READ(1,*)DPACK
      CLOSE(1)
      IF(IZZ.NE.KZ(1))then
         write(6,6001)zlab(1),izz,nn,kz(1)
         ierr=1
         return
      endif
      NTOT=NF
      IF(NTOT.GT.nptot)then
         write(6,6002)ntot,nptot
         ierr=2         
         return
      endif
      INT(1)=1
      IF(ITTE3.NE.ITE3)then
         write(6,6077)ite3,itte3,nn
         ierr=3
         return
      endif
c      ITE1=MAX(ITE1,ITTE1)
c      ITE2=MIN(ITE2,ITTE2)
C
c  READ MESH FILES                 
      OPEN(1,FILE=trim(path)//'/mono/'//ZLAB(1)//'.mesh',status='old',
     + form='unformatted',err=1011)
      READ(1)DV,NTOTV,(UMESH(N),N=1,NTOTV)   
      umin=umesh(1)
      umax=umesh(ntotv)
      DV1=DV
      CLOSE(1)               
C      
C  GET MESH FOR SCREEN   
!      CALL IMESH(UMESH,NTOTV)   !HH: Not needed
C
C     SUBSEQUENT FILES
      DO 40 N=2,ipe
         NN=N
         OPEN(1,FILE=trim(path)//'/mono/'//ZLAB(N)//'.index',
     +   STATUS='OLD')
         READ(1,*)IZZ,AMM
         READ(1,*)ITE11,ITE22,ITE33
         READ(1,*)UMINN,UMAXX
         READ(1,*)NC,NF
         READ(1,*)DPACK
         CLOSE(1)
         IF(ITE33.NE.ITE3)then
            write(6,6077)ite3,ite33,nn
            ierr=4
            return
         endif
c         ITE1=MAX(ITE1,ITE11)
c         ITE2=MIN(ITE2,ITE22)
         IF(IZZ.NE.KZ(N))then
            write(6,6001)zlab(n),izz,nn,kz(nn)
            ierr=5
            return
         endif
         NTOTT=NF
         IF(NTOTT.GT.NTOT)then
            write(6,6006)nn,ntott,ntot
            ierr=6
            return
         endif
c!!         IF(UMIN.NE.UMINN.OR.UMAX.NE.UMAXX) GOTO 1003  !!
         INT(N)=NTOT/NTOTT
         IF(INT(N)*NTOTT.NE.NTOT)then
            WRITE(6,6009)NN,NTOTT,NTOT
            ierr=7
            return
         endif
         IF(INT(N).NE.1)WRITE(6,6007)N,INT(N)
c
c        READ MESH FILES
c
         OPEN(1,FILE=trim(path)//'/mono/'//ZLAB(N)//'.mesh',
     +   status='old',form='unformatted',err=1011)
         READ(1)DV
       IF(DV.NE.DV1)THEN
          PRINT*,' OP: N=',N,', DV=',DV,' NOT EQUAL TO DV1=',DV1
            ierr=8
          return
       ENDIF
         CLOSE(1)
   40 CONTINUE
C
C  START TEMPERATURE LOOP
C
      ncount1=0
      ncount2=0
      ncount3=0
      do it=ite1,ite2,ite3
c
C        OPEN FILES
c
         TLAB='.'//NUM(IT/100)//NUM(IT/10-10*(IT/100))//
     +    NUM(IT-10*(IT/10))
         do n=1,ipe
c            IF(SKIP(N))GOTO 70
            NN=N
            OPEN(IFL(N),FILE=trim(path)//'/mono/'//ZLAB(N)//TLAB,
     +      FORM='UNFORMATTED',STATUS='OLD',ERR=1004)
            if(n.gt.2) then
            OPEN(IFLP(N),FILE=trim(path)//'/mono/'//ZLABP(N)//TLAB,
     +      FORM='UNFORMATTED',STATUS='OLD',ERR=1004)
            endif
         enddo
C        READ HEADINGS
         NN=1
         READ(IFL(1))IZZ,ITE,AM,UM,UX,NCCC,NFFF,DP,JNE1,JNE2,JNE3
         do n=2,ipe
c            IF(SKIP(N))GOTO 80
            NN=N
            READ(IFL(N))IZZ,ITE,AM,UM,UX,NC,NF,DP,JNE11,JNE22,JNE33
            if(n.gt.2) read(iflp(n))
            IF(JNE33.NE.JNE3)then
              write(6,6099)jne3,jne33,nn
              ierr=9
              return
            endif
            JNE1=MAX(JNE1,JNE11)
            JNE2=MIN(JNE2,JNE22)
         enddo
         itt=(it-ite1)/2+1
         jn1(itt)=jne1
         jn2(itt)=jne2
C
c         WRITE(98,9802)ITE,JNE1,JNE2,JNE3
C
C        START DENSITY LOOP
C
         do n=1,ipe
           do jn=jne1,jne2,jne3
             jnn=(jn-jne1)/2+1
C
C           START LOOP ON ELEMENTS
C
   95        READ(IFL(N))JNE,EPATOM(n,itt,jnn),OPLNCK(n,itt,jnn),ORSS,
     +         NE1(n,itt,jnn),NE2(n,itt,jnn),
     +         (FION(NE,n,itt,jnn),NE=NE1(n,itt,jnn),NE2(n,itt,jnn))
             read(ifl(n))np(n,itt,jnn)
             if(np(n,itt,jnn).gt.0)then
          read(ifl(n))(mx(k+ncount1),yy1(k+ncount1),k=1,np(n,itt,jnn))
                 kp1(n,itt,jnn)=ncount1
                 ncount1=ncount1+np(n,itt,jnn)
             else
               read(ifl(n))(yy2(k+ncount2),k=1,ntot)
                 kp2(n,itt,jnn)=ncount2
                 ncount2=ncount2+ntot
             endif
               if(n.gt.2) then
                 read(iflp(n))ja,npp(n,itt,jnn)
                 if(npp(n,itt,jnn).gt.0) then
           read(iflp(n))(nx(k+ncount3),yx(k+ncount3),k=1,npp(n,itt,jnn))
                   kp3(n,itt,jnn)=ncount3
                   ncount3=ncount3+npp(n,itt,jnn)
                 endif
               endif
            enddo
          enddo
c          
c     write(6,610)it
c     write(6,*)'ncount1 = ',ncount1
c     write(6,*)'ncount2 = ',ncount2
c     write(6,*)'ncount3 = ',ncount3
c
C        CLOSE FILES
c
         DO 150 N=1,ipe
            CLOSE(IFL(N))
            close(iflp(n))
  150    CONTINUE
c
      enddo      
c     write(6,*)'ncount1 = ',ncount1
c     write(6,*)'ncount2 = ',ncount2
c     write(6,*)'ncount3 = ',ncount3
      return
610   format(10x,'Done IT= ',i3)
1004  WRITE(6,6004)ZLAB(NN),TLAB
      STOP
6001  FORMAT(//5X,'*** OP: FILE ',A3,' GIVES IZZ=',I3,
     + 'NOT EQUAL TO IZ(',I2,')=',I2,' ***')
6002  FORMAT(//5X,'*** OP: NTOT=',I7,' GREATER THAN nptot=',
     + I7,' ***')
6003  FORMAT(//5X,'*** OP: DISCREPANCY BETWEEN DATA ON FILES ',
     + A3,' AND ',A3,' ***')
6004  FORMAT(//5X,'*** OP: ERROR OPENING FILE ',A3,A6,'  ***')
6006  FORMAT(//5X,'OP: N=',I2,', NTOTT=',I7,', GREATER THAN NTOT=',I7)
6007  FORMAT(/5X,'OP: N=',I2,', INT(N)=',I4)
6009  FORMAT(' OP: N=',I5,', NTOTT=',I10,', NTOT=',I10/
     + '   NTOT NOT MULTIPLE OF NTOTT')
c6012  FORMAT(/10X,'ERROR, SEE WRITE(6,6012)'/
c     + 10X,'IT=',I3,', JN=',I3,', N=',I3,', JNE=',I3/)
6077  FORMAT(//5X,'OP: DISCREPANCY IN ITE3'/10X,I5,' READ FROM UNIT 5'/
     + 10X,I5,' FROM INDEX FILE ELEMENT',I5)
6099  FORMAT(//5X,'OP: DISCREPANCY IN JNE3'/10X,I5,' READ FOR N=1'/
     + 10X,I5,' READ FOR N=',I5)

c8000  FORMAT(5X,I5,F10.4/5X,3I5/2E10.2/2I10/10X,E10.2)
1010  print*,' OP: ERROR OPENING FILE '//'./mono/'//ZLAB(1)//'.index'
      stop
1011  print*,' OP: ERROR OPENING FILE '//'./mono/'//ZLAB(1)//'.mesh'
      stop
      end subroutine op_dload
C******************************************************************
c HH: Based on "op_ax.f"
c Input:   kk = number of elements to calculate g_rad for
c          iz1(kk) = charge of element to calculate g_rad for
c          nel = number of elements in mixture
c          izzp(nel) = charge of elements
c          fap(nel) = number fractions of elements
c          flux = local radiative flux (Lrad/4*pi*r^2)
c          fltp = log T
c          flrhop = log rho
c Output: g1 = log kappa
c         gx1 = d(log kappa)/d(log T)
c         gy1 = d(log kappa)/d(log rho)
c         gp1(kk) = d(log kappa)/d(log xi) 
c         grl1(kk) = log grad
c         fx1(kk) = d(log grad)/d(log T) 
c         fy1(kk) = d(log grad)/d(log rho)
c         gr1p1(kk) = d(log grad)/d(log xi)
c         meanZ(nel) = average ionic charge of elements
c         zetx1(nel) = d(meanZ)/d(log T) 
c         zety1(nel) = d(meanZ)/d(log rho)
c         ierr = 0 for correct use 
      subroutine op_radacc(kk, izk, nel, izzp, fap, flux, fltp, flrhop, 
     : g1, gx1, gy1, gp1, grl1, fx1, fy1, grlp1, meanZ, zetx1, zety1, ierr)
      use opax
      implicit none
      integer, intent(in) :: kk, nel
      integer, intent(in) :: izk(kk), izzp(nel)
      real(kind=iprec), intent(in) :: fap(nel)
      real(kind=iprec), intent(in) :: flux, fltp
      real(kind=iprec), intent(inout) :: flrhop
      real(kind=iprec), intent(out) :: g1, gx1, gy1
      real(kind=iprec), intent(out) :: grl1(kk), meanZ(nel), grlp1(kk), gp1(kk), fx1(kk), fy1(kk),
     :  zetx1(nel), zety1(nel) 
      integer,intent(out) :: ierr
c local variables      
      integer :: n, i, k2, i3, ntot, jhmin, jhmax
      integer :: ih(4), jh(4), ilab(4), kzz(nrad), izz(ipe), iz1(nrad)
      real :: const, gx, gy, flt, flrho, flmu, dscat, dv, xi, flne,
     : epa, eta, ux, uy, g
      real :: umesh(nptot), uf(0:100), rion(28,4,4), rossl(4,4), flr(4,4),
     : ff(nptot, ipe, 4, 4), rr(28, ipe, 4, 4), ta(nptot, nrad, 4, 4), fa(ipe), 
     : gaml(4, 4, nrad), f(nrad), zetal(ipe, 4, 4), zetb(ipe), am1(nrad), 
     : rs(nptot, 4, 4), s(nptot, nrad, 4, 4), gamlp(4, 4, nrad), fp(nrad),
     : fmu1(nrad), rosslp(4, 4, nrad), gp(nrad), fx(nrad), fy(nrad),
     : zetx(ipe), zety(ipe)
c  Save variables to avoid blow up in certain compilers. CM, 06/06
      save ff, ta, rs, s       
c
c  Initialisations
      ierr = 0      
      if(nel.le.0.or.nel.gt.ipe) then
         write(6,*)'OP - NUMBER OF ELEMENTS OUT OF RANGE:', nel
         ierr = 1
         return
      endif
c  Get i3 for mesh type q='m'
      i3=2      
c      
c HH: k2 loops over elements for which to calculate grad.
      do k2 = 1, kk
         do n = 1, ipe
            if(izk(k2).eq.kz(n)) then
               iz1(k2) = izk(k2)
               exit
            endif   
            if(n.eq.ipe) then
               write(6,*)'OP - SELECTED ELEMENT CANNOT BE TREATED: IZ1 = ', izk(k2)
               ierr = 5
               return
            endif
         enddo
      enddo
c      
      outer: do i = 1, nel
         inner: do n = 1, ipe
            if(izzp(i).eq.kz(n)) then
               izz(i) = izzp(i)
               fa(i) = fap(i)
               if(fa(i).lt.0.0) then
                  write(6,*)'OP - NEGATIVE FRACTIONAL ABUNDANCE:',fa(i)
                  ierr = 7
                  return
               endif
               cycle outer
            endif
         enddo inner
         write(6,*)'OP - CHEM. ELEMENT CANNOT BE INCLUDED: Z = ', izzp(i)
         ierr = 8
         return
      enddo outer
c
c Calculate mean atomic weight (flmu) and 
c array kzz indicating elements for which to calculate g_rad
      call abund(nel, izz, kk, iz1, fa,     !input variables 
     :   kzz, flmu, am1, fmu1)              !output variables
c           
c  Other initialisations
c       dv = interval in frequency variable v
c       ntot=number of frequency points
c       umesh, values of u=(h*nu/k*T) on mesh points
c       uf, dscat used in scattering correction
      call msh(dv, ntot, umesh, uf, dscat)  !output variables
c
c  Start loop on temperature-density points
c  flt=log10(T, K)
c  flrho=log10(rho, cgs)
c
      flt = fltp
      flrho = flrhop
c     Get temperature indices
c       Let ite(i) be temperature index used in mono files
c       Put ite(i)=2*ih(i)
c       Use ih(i), i=1 to 4
c       xi=interpolation variable
c       log10(T)=flt=0.025*(ite(1)+xi+3)
c       ilab(i) is temperature label
      call xindex(i3, flt,             
     :  ih, ilab, xi) 
c      
c     Get density indices
c       Let jne(j) be density index used  in mono files
c       Put jne(j)=2*jh(j)
c       Use jh(j), j=1 to 4
c       Get extreme range for jh
      call jrange(i3, ih,
     : jhmin, jhmax)
c
c     Get electron density flne=log10(Ne) for specified mass density flrho
c       Also:  UY=0.25*[d log10(rho)]/[d log10(Ne)] 
c              epa=electrons per atom
      call findne(i3, ih, ilab, jhmin, jhmax, nel, fa, flrho, flt, xi, flmu,
     : flne, flr, epa, uy, ierr)
      if (ierr .ne. 0 ) return
c
c     Get density indices jh(j), j=1 to 4,
c       Interpolation variable eta
c       log10(Ne)=flne=0.25*(jne(1)+eta+3)
      call yindex(i3, jhmin, jhmax, 
     : flne, jh, eta)
c      	
c     Get ux=0.025*[d log10(rho)]/[d log10(T)]
      call findux(flr, xi, eta, 
     : ux)
            
c    rossl(i,j)=log10(Rosseland mean) on mesh points (i,j)
c     Get new mono opacities, ff(n,k,i,j)
      call rd(i3, kk, kzz, nel, izz, ilab, jh, ntot, umesh,
     : ff, rr, ta, zetal)
c
c     Get rs = weighted sum of monochromatic opacity cross sections
      call mix(kk, kzz, ntot, nel, fa, ff, rr,
     : rs, rion, s)
c
c     Screening corrections      
      if (screening) then 
c        data in /COMMON/CIMESH/ used for screening correction
         call imesh(umesh, ntot)
c     
c        Get Boercker scattering correction
         call scatt(ih, jh, rion, uf, rs, umesh, dscat, ntot, epa, ierr)
         if (ierr .ne. 0) return 
c
c        Get correction for Debye screening
         call screen1(ih, jh, rion, umesh, ntot, epa, rs)
      endif               
c      
c     Get rossl, array of log10(Rosseland mean in cgs)    
      call ross(kk, flmu, fmu1, dv, ntot, rs, s, 
     : rossl, gaml, ta, rosslp, gamlp)
c
c     Interpolate to required flt, flrho
c     g=log10(ross, cgs)
      call interp(nel, kk, rossl, gaml, xi, eta, g, i3, f, zetal, 
     : zetb, zetx, zety, ux, uy, gx, gy, rosslp, gp, gamlp, fp, fx, fy)
c      
c Write grad in terms of local radiative flux instead of (Teff, r/R*):
      const = 13.30295 + log10(flux) ! -log10(c) - log10(amu) + log10(flux)
      do k2 = 1, kk 
         gp1(k2) = gp(k2)            
         grl1(k2) = const + flmu - log10(am1(k2)) + f(k2) + g    ! log g_rad 
         grlp1(k2) = fmu1(k2)/10.d0**flmu + fp(k2) + gp(k2)      ! d(log g_rad)/d(log xi)
         fx1(k2) = fx(k2) + gx
         fy1(k2) = fy(k2) + gy
      enddo   !k2
      zetx1(1:nel) = zetx(1:nel)
      zety1(1:nel) = zety(1:nel)
c      
      g1 = g                ! log kappa
      gx1 = gx              ! dlogkappa/dlogt
      gy1 = gy              ! dlogkappa/dlogrho   
      meanZ(1:nel) = zetb(1:nel)  ! average ionic chanrge
      flrhop = flrho        ! take min/max allowed value of log rho
c
      return
c
      end subroutine op_radacc
c***********************************************************************
c HH: Based on "op_mx.f", opacity calculations to be used for stellar evolution calculations 
c Input:   nel = number of elements in mixture
c          izzp(nel) = charge of elements
c          fap(nel) = number fractions of elements
c          fltp = log (temperature)
c          flrhop = log (mass density) 
c Output: g1 = log kappa
c         gx1 = d(log kappa)/d(log T)
c         gy1 = d(log kappa)/d(log rho)
c         ierr = 0 for correct use 
      subroutine op_ev(nel, izzp, fap, fltp, flrhop, g1, gx1, gy1, gp1, ierr)
      use opac_ev
      implicit none
      integer, intent(in) :: nel
      integer, intent(in) :: izzp(nel)
      real(kind=iprec), intent(in) :: fap(nel)
      real(kind=iprec), intent(in) :: fltp, flrhop
      real(kind=iprec), intent(out) :: g1, gx1, gy1, gp1(nel)
      integer,intent(out) :: ierr
c local variables      
      integer :: n, i, i3, jhmin, jhmax, ntot
      integer :: ih(4), jh(4), ilab(4), izz(ipe) 
      real :: flt, flrho, flmu, flne, dv, dscat, const, gx, gy, g,
     : eta, epa, xi, ux, uy
      real :: umesh(nptot), uf(0:100), rion(28, 1:4, 1:4), rossl(4, 4), flr(4, 4), 
     : ff(nptot, ipe, 4, 4), rs(nptot, 4, 4), fa(ipe), rr(28, ipe, 4, 4),
     : s(nptot, nrad, 4, 4), rosslp(4, 4, nrad), gp(nrad),  fmu1(nrad)
c  Save variables to avoid blow up in certain compilers. CM, 06/06
      save ff, rs, s
c
c  Initialisations
      ierr=0      
      if(nel.le.0.or.nel.gt.ipe) then
         write(6,*)'OP - NUMBER OF ELEMENTS OUT OF RANGE:',nel
         ierr=1
         return
      endif
c      
c  Get i3 for mesh type q='m'
      i3=2
c      
        outer: do i=1,nel
          inner: do n=1,ipe
            if(izzp(i).eq.kz(n)) then
              izz(i)=izzp(i)
              fa(i)=fap(i)
              if(fa(i).lt.0.0) then
                write(6,*)'OP - NEGATIVE FRACTIONAL ABUNDANCE:',fa(i)
                ierr=7
                return
              endif
              cycle outer
            endif
          enddo inner
          write(6,*)'OP - CHEM. ELEMENT CANNOT BE INCLUDED: Z = ',
     +      izzp(i)
          ierr=8
          return
       enddo outer

c Calculate mean atomic weight (flmu) 
      call abund(nel, izz, fa, flmu, fmu1)
c          
c  Other initialisations
c       dv = interval in frequency variable v
c       ntot=number of frequency points
c       umesh, values of u=(h*nu/k*T) on mesh points
c       uf, dscat used in scattering correction
      call msh(dv, ntot, umesh, uf, dscat)
c
c  Start loop on temperature-density points
c  flt=log10(T, K)
c  flrho=log10(rho, cgs)
      flt = fltp
      flrho = flrhop
c     Get temperature indices
c       Let ite(i) be temperature index used in mono files
c       Put ite(i)=2*ih(i)
c       Use ih(i), i=1 to 4
c       xi=interpolation variable
c       log10(T)=flt=0.025*(ite(1)+xi+3)
c       ilab(i) is temperature label
      call xindex(flt, ilab, xi, ih, i3)
c      
c     Get density indices
c       Let jne(j) be density index used  in mono files
c       Put jne(j)=2*jh(j)
c       Use jh(j), j=1 to 4
c       Get extreme range for jh
      call jrange(ih, jhmin, jhmax, i3)   
c
c     Get electron density flne=log10(Ne) for specified mass density flrho
c       Also:  UY=0.25*[d log10(rho)]/[d log10(Ne)] 
c              epa=electrons per atom
      call findne(ilab, fa, nel, jhmin, jhmax, ih, flrho, flt,
     + xi, flne, flmu, flr, epa, uy, i3)
c
c     Get density indices jh(j), j=1 to 4,
c       Interpolation variable eta
c       log10(Ne)=flne=0.25*(jne(1)+eta+3)
      call yindex(jhmin, jhmax, flne, jh, i3, eta)
c      	
c     Get ux=0.025*[d log10(rho)]/[d log10(T)]
      call findux(flr, xi, eta, ux)
            
c    rossl(i,j)=log10(Rosseland mean) on mesh points (i,j)
c     Get new mono opacities, ff(n,k,i,j)
      call rd(nel, izz, ilab, jh, ntot, ff, rr, i3, umesh)
c
c     Up-date mixture
      call mix(ntot, nel, fa, ff, rs, rr, rion, s)  
c
      if(screening) then 
c        data in /COMMON/CIMESH/ used for screening correction
         call imesh(umesh, ntot)      
c        Get Boercker scattering correction
         call scatt(ih, jh, rion, uf, rs, umesh, dscat, ntot, epa)
c        Get correction for Debye screening
         call screen1(ih, jh, rion, umesh, ntot, epa, rs)
      endif     
c      
c     Get rossl, array of log10(Rosseland mean in cgs)
      call ross(flmu, fmu1, dv, ntot, rs, s, rossl, rosslp)
c      
c     Interpolate to required flt, flrho
c     g=log10(ross, cgs)
      call interp(nel, rossl, rosslp, xi, eta, 
     : g, i3, ux, uy, gx, gy, gp)
c            
      gp1(1:nel) = gp(1:nel)
      g1 = g                ! log kappa
      gx1 = gx              ! dlogkappa/dt
      gy1 = gy              ! dlogkappa/drho   
c        
      return
c
      end subroutine op_ev
***********************************************************************
c HH: Based on "op_mx.f", opacity calculations to be used for non-adiabatic pulsation calculations
c Special care is taken to ensure smoothness of opacity derivatives
c Input:   nel = number of elements in mixture
c          izzp(nel) = charge of elements
c          fap(nel) = number fractions of elements
c          fltp = log (temperature)
c          flrhop = log (mass density) 
c Output: g1 = log kappa
c         gx1 = d(log kappa)/d(log T)
c         gy1 = d(log kappa)/d(log rho)
c         ierr = 0 for correct use 
      subroutine op_osc(nel, izzp, fap, fltp, flrhop, g1, gx1, gy1, ierr)
      use opac_osc
      implicit none
      integer, intent(in) :: nel
      integer, intent(in) :: izzp(nel)
      real(kind=iprec), intent(in) :: fap(nel)
      real(kind=iprec), intent(in) :: fltp, flrhop
      real(kind=iprec), intent(out) :: g1, gx1, gy1
!      real(kind=iprec), intent(out) :: meanZ(nel)
      integer,intent(out) :: ierr
c local variables      
      integer :: n, i, i3, jhmin, jhmax, ntot
      integer :: ih(0:5), jh(0:5), ilab(0:5), izz(ipe) 
      real :: flt, flrho, flmu, flne, dv, dscat, const, gx, gy, g,
     : eta, epa, xi, ux, uy
      real :: umesh(nptot), uf(0:100), rion(28, 0:5, 0:5), rossl(0:5, 0:5), flr(4, 4), 
     : ff(nptot, ipe, 0:5, 0:5), rs(nptot, 0:5, 0:5), fa(ipe), rr(28, ipe, 0:5, 0:5)
c  Save variables to avoid blow up in certain compilers. CM, 06/06
      save ff, rs
c
c  Initialisations
      ierr=0      
      if(nel.le.0.or.nel.gt.ipe) then
         write(6,*)'OP - NUMBER OF ELEMENTS OUT OF RANGE:',nel
         ierr=1
         return
      endif
c      
c  Get i3 for mesh type q='m'
      i3=2
c      
        outer: do i=1,nel
          inner: do n=1,ipe
            if(izzp(i).eq.kz(n)) then
              izz(i)=izzp(i)
              fa(i)=fap(i)
              if(fa(i).lt.0.0) then
                write(6,*)'OP - NEGATIVE FRACTIONAL ABUNDANCE:',fa(i)
                ierr=7
                return
              endif
              cycle outer
            endif
          enddo inner
          write(6,*)'OP - CHEM. ELEMENT CANNOT BE INCLUDED: Z = ',
     +      izzp(i)
          ierr=8
          return
       enddo outer

c Calculate mean atomic weight (flmu) 
      call abund(nel, izz, fa, flmu)
c          
c  Other initialisations
c       dv = interval in frequency variable v
c       ntot=number of frequency points
c       umesh, values of u=(h*nu/k*T) on mesh points
c       uf, dscat used in scattering correction
      call msh(dv, ntot, umesh, uf, dscat)
c
c  Start loop on temperature-density points
c  flt=log10(T, K)
c  flrho=log10(rho, cgs)
      flt = fltp
      flrho = flrhop
c     Get temperature indices
c       Let ite(i) be temperature index used in mono files
c       Put ite(i)=2*ih(i)
c       Use ih(i), i=1 to 4
c       xi=interpolation variable
c       log10(T)=flt=0.025*(ite(1)+xi+3)
c       ilab(i) is temperature label
      call xindex(flt, ilab, xi, ih, i3)
c      
c     Get density indices
c       Let jne(j) be density index used  in mono files
c       Put jne(j)=2*jh(j)
c       Use jh(j), j=1 to 4
c       Get extreme range for jh
      call jrange(ih, jhmin, jhmax, i3)   
c
c     Get electron density flne=log10(Ne) for specified mass density flrho
c       Also:  UY=0.25*[d log10(rho)]/[d log10(Ne)] 
c              epa=electrons per atom
      call findne(ilab, fa, nel, jhmin, jhmax, ih, flrho, flt,
     + xi, flne, flmu, flr, epa, uy, i3)
c
c     Get density indices jh(j), j=1 to 4,
c       Interpolation variable eta
c       log10(Ne)=flne=0.25*(jne(1)+eta+3)
      call yindex(jhmin, jhmax, flne, jh, i3, eta)
c      	
c     Get ux=0.025*[d log10(rho)]/[d log10(T)]
      call findux(flr, xi, eta, ux)
            
c    rossl(i,j)=log10(Rosseland mean) on mesh points (i,j)
c     Get new mono opacities, ff(n,k,i,j)
      call rd(nel, izz, ilab, jh, ntot, ff, rr, i3, umesh)
c
c     Up-date mixture
      call mix(ntot, nel, fa, ff, rs, rr, rion)  
c
      if(screening) then 
c        data in /COMMON/CIMESH/ used for screening correction
         call imesh(umesh, ntot)      
c        Get Boercker scattering correction
         call scatt(ih, jh, rion, uf, rs, umesh, dscat, ntot, epa)
c        Get correction for Debye screening
         call screen1(ih, jh, rion, umesh, ntot, epa, rs)
      endif     
c      
c     Get rossl, array of log10(Rosseland mean in cgs)
      call ross(flmu, dv, ntot, rs, rossl)
c      
c     Interpolate to required flt, flrho
c     g=log10(ross, cgs)
      call interp(nel, rossl, xi, eta, g, i3, ux, uy, gx, gy)
c      
      g1 = g                ! log kappa
      gx1 = gx              ! dlogkappa/dt
      gy1 = gy              ! dlogkappa/drho   
c        
      return
c
      end subroutine op_osc
c***********************************************************************
      end module opserver_dp
