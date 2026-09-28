      module opax

      integer, parameter, private :: nptot = 10000
      integer, parameter, public :: nrad = 17
      integer, parameter, private :: ipe = 17      
      integer, public :: kz(17)
      integer, private :: i, nkz(17)
      integer, private :: JS(140:320), JE(140:320)
      real, private :: amass(17)
      character(len=2), dimension(17), private ::  name
      
      DATA (KZ(I),NAME(I),AMASS(I),I=1,17)/
     +   1,     'H ',     1.0080,
     +   2,     'He',     4.0026,
     +   6,     'C ',    12.0111,
     +   7,     'N ',    14.0067,
     +   8,     'O ',    15.9994,
     +  10,     'Ne',    20.179,
     +  11,     'Na',    22.9898,
     +  12,     'Mg',    24.305,
     +  13,     'Al',    26.9815,
     +  14,     'Si',    28.086,
     +  16,     'S ',    32.06,
     +  18,     'Ar',    39.948,
     +  20,     'Ca',    40.08,
     +  24,     'Cr',    51.996,
     +  25,     'Mn',    54.9380,
     +  26,     'Fe',    55.847,
     +  28,     'Ni',    58.71/

      data js/
     +  14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 15, 18, 19, 22, 23,
     +  26, 27, 30, 31, 34, 34, 34, 34, 36, 36, 36, 36, 36, 36, 38, 38,
     +  38, 38, 38, 39, 40, 40, 40, 40, 40, 41, 42, 42, 42, 42, 42, 43,
     +  44, 44, 44, 44, 44, 45, 46, 46, 46, 46, 46, 47, 48, 48, 48, 48,
     +  48, 49, 50, 50, 50, 50, 52, 52, 52, 52, 52, 52, 54, 54, 54, 54,
     +  54, 54, 56, 56, 56, 56, 56, 56, 58, 58, 58, 58, 58, 58, 60, 60,
     +  60, 60, 60, 60, 62, 62, 62, 62, 62, 63, 64, 64, 64, 64, 64, 65,
     +  66, 66, 66, 66, 66, 67, 68, 68, 68, 68, 68, 69, 70, 70, 70, 70,
     +  70, 71, 72, 72, 72, 72, 72, 73, 74, 74, 74, 74, 76, 76, 76, 76,
     +  76, 76, 78, 78, 78, 78, 78, 78, 80, 80, 80, 80, 80, 80, 82, 82,
     +  82, 82, 82, 82, 84, 84, 84, 84, 84, 84, 86, 86, 86, 86, 86, 87,
     +  88, 88, 88, 88, 88/

      data je/
     + 52, 56, 56, 58, 58, 60, 60, 60, 60, 62, 62, 64, 64, 66, 66, 68,
     + 68, 70, 70, 72, 72, 74, 74, 74, 74, 76, 76, 78, 78, 80, 80, 80,
     + 80, 81, 82, 82, 82, 82, 82, 83, 84, 84, 84, 84, 84, 86, 86, 86,
     + 86, 86, 86, 88, 88, 88, 88, 89, 90, 90, 90, 88, 88, 88, 88, 92,
     + 92, 92, 92, 94, 94, 94, 94, 94, 94, 89, 90, 90, 90, 90, 90, 91,
     + 92, 92, 92, 92, 92, 94, 94, 90, 90, 92, 92, 92, 92, 92, 92, 93,
     + 94, 94, 94, 94, 94, 94, 94, 96, 96, 96, 96, 96, 96, 98, 98, 98,
     + 98, 98, 98, 99,100,100,100,100,100,100,100,102,102,102,102,102,
     +102,103,104,104,104,104,104,104,104,106,106,106,106,106,106,107,
     +108,108,108,108,108,108,108,110,110,110,110,110,110,111,112,112,
     +112,112,112,112,112,114,114,114,114,114,114,115,116,116,116,116,
     +116,116,116,118,118/

      contains
c***********************************************************************1
      subroutine abund(nel, izz, kk, iz1, fa, kzz, flmu, am1, fmu1)
      implicit none
      integer, intent(in) :: nel, izz(ipe), kk, iz1(nrad)
      real, intent(in) :: fa(ipe)
      integer,  intent(out) :: kzz(nrad)
      real, intent(out) :: flmu, am1(nrad), fmu1(nrad)
c local variables     
      integer :: k, k1, k2, m
      real :: fmu, a1, c1, fmu0, amamu(ipe)
c
c Get k1,get amamu(k)
      do k = 1, nel              
         do m = 1, ipe
            if(izz(k).eq.kz(m))then
               amamu(k) = amass(m)
               nkz(k) = m
            goto 1
            endif
         enddo  
         print*,' k=',k,', izz(k)=',izz(k)
         print*,' kz(m) not found'
         stop
    1    continue           
      enddo
c
      do k2 = 1, kk
         do k = 1, nel
            if(izz(k).eq.iz1(k2))then
               k1 = k
            endif
         enddo    
         kzz(k2) = k1
         am1(k2) = amamu(k1)  
      enddo       
c   
c  Mean atomic weight = fmu
      fmu = 0.
      do k = 1, nel
         fmu = fmu + fa(k)*amamu(k)
      enddo
c
      do k2 = 1, kk
         a1 = fa(kzz(k2))
         c1 = 1./(1.-a1)
         if ( a1 < 1.d-10 ) then                       ! for very small f, get derivative to f instead of log xi
            fmu0 = fmu - fa(kzz(k2))*amamu(kzz(k2))  
            fmu1(k2) = (am1(k2)-fmu0)*1.660531e-24     !dmu/df         
         else
            fmu0 = c1*(fmu - fa(kzz(k2))*amamu(kzz(k2)))
            fmu1(k2) = a1*(am1(k2)-fmu0)*1.660531e-24  !dmu/dlog xi
         endif   
      enddo   
c      
      fmu = fmu*1.660531e-24 ! Convert to cgs
      flmu = log10(fmu)
c  
      return
      end subroutine abund
c**********************************************************************
      subroutine msh(dv, ntot, umesh, uf, dscat)
      implicit none
      integer, intent(out) :: ntot
      real, intent(out) :: dv, umesh(nptot), uf(0:100), dscat
      integer :: i, k, ntotv
      real :: dvp, dv1, umin, umax, umeshp(nptot)
      common /mesh/ ntotv, dvp, dv1, umeshp    
      save /mesh/    
c
      ntot = ntotv
      dv = dvp
      umesh(1:ntot) = umeshp(1:ntot)
c
      umin = umesh(1)
      umax = umesh(ntot)
      dscat = (umax - umin)*0.01
      do i = 0, 100
         uf(i) = umin + i*dscat
      enddo  
c
      return
c
      end subroutine msh
c**********************************************************************
      subroutine xindex(i3, flt, ih, ilab, xi)
      implicit none
      integer, intent(in) :: i3
      real, intent(in) :: flt
      integer, intent(out) :: ih(4), ilab(4)
      real, intent(out) :: xi
      integer :: i, ih2
      real :: x
c
      if(flt.lt.3.5) then
        write(6,600) flt
        stop
      elseif(flt.gt.8.) then
        write(6,601) flt
        stop
      endif
c      
      x = 40.*flt/real(i3)
      ih2 = x
      ih2 = max(ih2, 140/i3+1)
      ih2 = min(ih2, 320/i3-2)
      do i = 1, 4
         ih(i) = ih2 + i - 2
         ilab(i) = i3*ih(i)
      enddo
         
      xi = 2.*(x-ih2) - 1
c
      return
c
  600 format(' flt=',1p,e11.3,' less than lower limit of 3.5')
  601 format(' flt=',1p,e11.3,' larger than upper limit of 8')
      end subroutine xindex
c**********************************************************************
      subroutine jrange(i3, ih, jhmin, jhmax) 
      implicit none
      integer, intent(in) :: ih(4), i3
      integer, intent(out) :: jhmin, jhmax  
      integer :: i
c      
      jhmin = 0
      jhmax = 1000
      do i = 1, 4
         jhmin = max(jhmin, js(ih(i)*i3)/i3)
         jhmax = min(jhmax, je(ih(i)*i3)/i3)
      enddo
c      
      return
      end subroutine jrange
c**********************************************************************
      subroutine findne(i3, ih, ilab, jhmin, jhmax, nel, fa, flrho, flt, xi, flmu, 
     : flne, flr, epa, uy, ierr)
      implicit none
      integer, intent(in) :: i3, ih(4), ilab(4), nel, jhmin
      integer, intent(inout) :: jhmax, ierr
      real, intent(in) :: fa(ipe), flt, xi, flmu
      real,intent(out) :: flne, uy, epa, flr(4,4)
      real, intent(inout) :: flrho
c local variables
      integer :: i, j, n, jh, jm, itt, jne, jnn
      real :: flrmin, flrmax, uyi(4), efa(4, 7:118), 
     : flrh(4, 7:118), u(4), flnei(4), y, zeta, efa_temp
c declare variables in common block, by default: real (a-h, o-z), integer (i-n)      
      integer :: ite1, ite2, ite3, jn1, jn2, jne3, ntot, nc, nf, int,
     : ne1, ne2, np, kp1, kp2, kp3, npp, mx, nx   
      real :: umin, umax, epatom, oplnck, fion, yy1, yy2, yx
      common/atomdata/ ite1,ite2,ite3,jn1(91),jn2(91),jne3,umin,umax,ntot,
     + nc,nf,int(17),epatom(17,91,25),oplnck(17,91,25),ne1(17,91,25),
     + ne2(17,91,25),fion(-1:28,28,91,25),np(17,91,25),kp1(17,91,25),
     + kp2(17,91,25),kp3(17,91,25),npp(17,91,25),mx(33417000),
     + yy1(33417000),yy2(120000000),nx(19305000),yx(19305000)
      save /atomdata/
c
c  efa(i,jh)=sum_n epa(i,jh,n)*fa(n)
c  flrh(i,jh)=log10(rho(i,jh))
c
c  Get efa        
      do i = 1, 4
         itt = (ilab(i)-ite1)/2 + 1
         do jne = jn1(itt), jn2(itt), i3
            jnn = (jne-jn1(itt))/2 + 1
            jh = jne/i3
            efa_temp = 0.
            do n = 1, nel               
               efa_temp = efa_temp + epatom(nkz(n), itt, jnn)*fa(n)
            enddo !n    
            efa(i, jh) = efa_temp             
         enddo !jne
      enddo !i 
c        
c Get range for efa.gt.0
      do i = 1, 4
         do jh = jhmin, jhmax
            if(efa(i, jh) .le. 0.)then
               jm = jh - 1
               goto 3
            endif
         enddo
         goto 4
    3    jhmax = MIN(jhmax, jm)
    4    continue
      enddo
c    
c Get flrh
!      do i = 1, 4
!        do jh = jhmin,jhmax
      forall(i=1:4, jh=jhmin:jhmax)
         flrh(i, jh) = flmu + 0.25*i3*jh - log10(efa(i, jh))
      end forall    
!        enddo
!      enddo  
c
c  Find flrmin and flrmax
      flrmin = -1000
      flrmax = 1000
      do i = 1, 4
         flrmin = max(flrmin, flrh(i,jhmin))
         flrmax = min(flrmax, flrh(i,jhmax))
      enddo
c      
c  Check range of flrho
      if(flrho .lt. flrmin .or. flrho .gt. flrmax)then
         write(6,601) flt, flrho, flrmin, flrmax
         if (flrho .lt. flrmin) flrho = flrmin
         if (flrho .gt. flrmax) flrho = flrmax
!         stop
      endif
c
c  Interpolations in j for flne
      do jh = jhmin, jhmax
         if(flrh(2,jh) .gt. flrho)then
            jm = jh - 1
            goto 5
         endif
      enddo
      print*,' Interpolations in j for flne'
      print*,' Not found, i=',i
      ierr = 9
      return
!      stop
    5 jm=max(jm,jhmin+1)
      jm=min(jm,jhmax-2)
c      
      do i = 1, 4
         do j = 1, 4
            u(j) = flrh(i, jm+j-2)
            flr(i,j) = flrh(i, jm+j-2)
         enddo
         call solve(u, flrho, zeta, uyi(i), ierr)
         y = jm + 0.5*(zeta+1)
         flnei(i) = .25*i3*y
      enddo  
c
c  Interpolations in i
      flne = fint(flnei, xi)
      uy = fint(uyi, xi)
c Get epa      
      epa = 10**(flne + flmu - flrho)     
c      
      return
c   
  601 format(' For flt=',1p,e11.3,', flrho=',e11.3,' is out of range'/
     +  ' Allowed range for flrho is ',e11.3,' to ',e11.3) 
      end subroutine findne
c***********************************************************************

      subroutine yindex(i3, jhmin, jhmax, flne, jh, eta)
      implicit none
      integer, intent(in) :: i3, jhmin, jhmax
      real, intent(in) :: flne
      integer, intent(out) :: jh(4)
      real, intent(out) :: eta
c local variables      
      integer :: j, k
      real :: y
c
      y = 4.*flne/real(i3)
      j = y
      j = max(j,jhmin+1)
      j = min(j,jhmax-2)
      do k = 1, 4
         jh(k) = j + k - 2
      enddo
      eta = 2.*(y-j)-1
c      
      return
      end subroutine yindex
c***********************************************************************
      subroutine findux(flr, xi, eta, ux)
      implicit none
      real, intent(in) :: flr(4, 4), xi, eta
      real, intent(out) :: ux
c local variables      
      integer :: i, j      
      real :: uxj(4), u(4)
c
      do j = 1, 4
         do i = 1, 4
            u(i) = flr(i, j)
         enddo
         uxj(j) = fintp(u, xi)
      enddo
      ux = fint(uxj, eta)
c
      return
      end subroutine findux   
c**********************************************************************
      subroutine rd(i3, kk, kzz, nel, izz, ilab, jh, ntot, umesh, 
     : ff, rr, ta, zetal)
      implicit none
      integer, intent(in) :: i3, kk, kzz(nrad), nel, izz(ipe), ilab(4), jh(4),
     : ntot
      real, intent(in) :: umesh(nptot)
      real, intent(out) :: ff(nptot, ipe, 4, 4), rr(28, ipe, 4, 4), 
     : ta(nptot, nrad, 4, 4), zetal(ipe, 4, 4)
c local variables   
      integer :: i, j, k, k2, l, n, m, itt, jnn, izp, ne1, ne2, ne, ib, ia
      real :: ya, yb, d, se, u, fion(-1:28), zet!, ff_temp(nptot), ta_temp(nptot)
c declare variables in common block, by default: real (a-h, o-z), integer (i-n)   
      integer :: ite1, ite2, ite3, jn1, jn2, jne3, ntotp, nc, nf, int,
     : ne1p, ne2p, np, kp1, kp2, kp3, npp, mx, nx   
      real :: umin, umax, epatom, oplnck, fionp, yy1, yy2, yx
      common/atomdata/ ite1,ite2,ite3,jn1(91),jn2(91),jne3,umin,umax,ntotp,
     + nc,nf,int(17),epatom(17,91,25),oplnck(17,91,25),ne1p(17,91,25),
     + ne2p(17,91,25),fionp(-1:28,28,91,25),np(17,91,25),kp1(17,91,25),
     + kp2(17,91,25),kp3(17,91,25),npp(17,91,25),mx(33417000),
     + yy1(33417000),yy2(120000000),nx(19305000),yx(19305000)
      save /atomdata/
c
c  i=temperature inex
c  j=density index
c  k=frequency index
c  n=element index
c  Get:
c    mono opacity cross-section ff(k,n,i,j)
c    modified cross-section for selected element, ta(k,i,j)
c    zet(i,j) for diffusion coefficient
c
c     Initialisations    
!      ff = 0.
      rr = 0.
!      ta = 0.
c     
c  Start loop on i (temperature index)
      do i = 1, 4
         itt = (ilab(i) - ite1)/2 + 1                 
c        Read mono opacities  
         do j = 1, 4
            jnn = (jh(j)*i3 - jn1(itt))/2 + 1            
            do n = 1, nel                       
               izp = izz(n)
               ne1 = ne1p(nkz(n), itt, jnn)
               ne2 = ne2p(nkz(n), itt, jnn)
               do ne = ne1, ne2
                  fion(ne) = fionp(ne, nkz(n), itt, jnn)
                  if (ne .le. min(ne2, izp-2)) rr(izp-1-ne, n, i, j) = fion(ne)
               enddo                           
               call zetbarp(izp, ne1, ne2, fion, zet, i3)
               zetal(n, i, j) = zet 
                         
               if(np(nkz(n), itt, jnn).eq.0) then
                  do k = 1, ntot
                     ff(k, n, i, j) = yy2(k+kp2(nkz(n), itt, jnn))
                  enddo  
               else  
                  ib = 1
                  yb = yy1(1+kp1(nkz(n),itt,jnn))
                  ff(1, n, i, j) = yb
                  do m = 2, np(nkz(n), itt, jnn)
                     ia = ib
                     ya = yb
                     ib = mx(m+kp1(nkz(n), itt, jnn))
                     yb = yy1(m+kp1(nkz(n), itt, jnn))
                     d = (yb-ya)/float(ib-ia)
                     do l = ia+1, ib-1
                        ff(l, n, i, j) = ya + (l-ia)*d
                     enddo
                     ff(ib, n, i, j) = yb
                  enddo
               endif 
               do k2 = 1, kk             
                  if (kzz(k2)== n) then                                    
                     ib = 1
                     yb = yx(1+kp3(nkz(kzz(k2)), itt, jnn))
                     u = umesh(ib)
                     if(u.lt.0.01) then
                        se = u*(1.-0.5*u)
                     else
                        se = 1. - exp(-u)
                     endif            
                     ta(1, k2, i, j) = se*ff(1, n, i, j) - yb
                     do m = 2, npp(nkz(kzz(k2)), itt, jnn)
                        ia = ib
                        ya = yb
                        ib = nx(m+kp3(nkz(kzz(k2)), itt, jnn))
                        yb = yx(m+kp3(nkz(kzz(k2)), itt, jnn))
                        d = (yb-ya)/float(ib-ia)
                        do l = ia+1, ib-1
                           u = umesh(l)
                           if(u.lt.0.01) then
                           se = u*(1.-0.5*u)
                           else
                              se = 1. - exp(-u)
                           endif   
                           ta(l, k2, i, j) = se*ff(l, n, i, j) -(ya + (l-ia)*d)
                        enddo
                        u = umesh(ib)
                        if(u.lt.0.01) then
                           se = u*(1.-0.5*u)
                        else
                           se = 1. - exp(-u)
                        endif            
                        ta(ib, k2, i, j) = se*ff(ib, n, i, j) - yb                        
                     enddo
                     goto 101  ! get out of k2-loop  
                  endif                  
               enddo !k2
 101           continue                                           
            enddo !n
         enddo !j  
      enddo !i
c
      return
c
      end subroutine rd
c***********************************************************************
      subroutine zetbarp(iz, ne1, ne2, fion, zet, i3)
      implicit none
      integer, intent(in) :: iz, i3, ne1, ne2
      real, intent(in) :: fion(-1:28)
      real, intent(out) :: zet
      integer :: ne
      real :: sz   !, fne, t, b     
c      
!      fne = 10**(0.25*i3*jhj)
!      t = 10**(0.025*i3*ihi)
!      b = 2.7285e8*t**3/fne
      zet = 0.
      do ne = ne1, ne2
         sz = iz - ne - 1
         zet = zet + fion(ne)*sz
      enddo
c        
      return
      end subroutine zetbarp
      
c***********************************************************************
      subroutine mix(kk, kzz, ntot, nel, fa, ff, rr, rs, rion, s)
      implicit none
      integer, intent(in) :: kk, kzz(nrad), ntot, nel
      real, intent(in) :: fa(ipe), ff(nptot, ipe, 4, 4), rr(28, 17, 4, 4)
      real, intent(out) :: rs(nptot, 4, 4), rion(28, 4, 4), s(nptot, nrad, 4, 4)
c local variables      
      integer :: i, j, k, n, m, k2
      real :: rs_temp, rion_temp, a1, c1
c
      do i = 1, 4
         do j = 1, 4
            do n = 1, ntot
               rs_temp = ff(n,1,i,j)*fa(1)
               do k = 2, nel
                  rs_temp = rs_temp + ff(n,k,i,j)*fa(k)
               enddo
               rs(n, i, j) = rs_temp  
            enddo
            
            do m = 1, 28
               rion_temp = rr(m, 1, i, j)*fa(1)
               do k = 2, nel
                  rion_temp = rion_temp + rr(m,k,i,j)*fa(k)
               enddo
               rion(m, i, j) = rion_temp
            enddo
             
            do k2 = 1, kk  
               a1 = fa(kzz(k2))
               c1 = 1./(1.-a1)                      
               do n = 1 , ntot
                  if ( a1 < 1.d-10) then
                     s(n, k2, i, j) = (1.+a1)*ff(n, kzz(k2), i, j) - rs(n, i, j)   ! d/d(fa) 
                  else
                     s(n, k2, i, j) =   rs(n, i, j) - c1*(rs(n, i, j) -            ! d/d(log xi)
     :                             ff(n, kzz(k2), i, j)*a1 )               
!                   a1*(ff(n, kzz(k2), i, j) - 
!     :             c1*(rs(n, i, j) - ff(n, kzz(k2), i, j)*a1)) 
                  endif                    
               enddo   
            enddo                   
         enddo
      enddo
c
      return
      end subroutine mix
c***********************************************************************
      subroutine ross(kk, flmu, fmu1, dv, ntot, rs, s, 
     : rossl, gaml, ta, rosslp, gamlp) 
      implicit none
      integer, intent(in) :: kk, ntot
      real, intent(in) :: flmu, dv, rs(nptot, 4, 4), ta(nptot, nrad, 4, 4), 
     : s(nptot, nrad, 4, 4), fmu1(nrad)
      real, intent(out) :: rossl(4, 4), gaml(4, 4, nrad), gamlp(4, 4, nrad),
     : rosslp(4, 4, nrad)
c local variables      
      integer :: k2, i, j, n
      double precision :: drs, dd,  dgm(nrad)
      real :: fmu, oross, tt, dd2, ss, drsp(nrad), dgmp(nrad)
c
c  oross=cross-section in a.u.
c  rossl=log10(ROSS in cgs)
         do i = 1, 4
            do j = 1, 4
               drs = 0.d0
               drsp(:) = 0.
               dgm(:) = 0.d0
               dgmp(:) = 0.
               do n = 1, ntot   !10000
                  dd = 1.d0/rs(n, i, j)  
                  dd2 = dd**2                                                
                  drs = drs + dd        
                  do k2 = 1, kk
                     ss = s(n, k2, i, j)   
                     drsp(k2) = drsp(k2) + ss*dd2                                      
                     tt = ta(n, k2, i, j)                
                     dgm(k2) = dgm(k2) + tt*dd
                     dgmp(k2) = dgmp(k2) + tt*ss*dd2 
                  enddo   
               enddo              
               oross = 1./(drs*dv)               
               rossl(i, j) = log10(oross) - 16.55280 - flmu
               do k2 = 1, kk
                  drsp(k2) = drsp(k2)*dv 
                  rosslp(i, j, k2) = oross*drsp(k2)-fmu1(k2)/10.**flmu     
                  if(dgm(k2).gt.0) then
                     dgm(k2) = dgm(k2)*dv                     
                     gaml(i, j, k2) = log10(dgm(k2))
                     dgmp(k2) = dgmp(k2)*dv
                     gamlp(i, j, k2) = - dgmp(k2)/dgm(k2) 
                  else
                     gaml(i, j, k2) = -30. 
                     gamlp(i, j, k2) = 0.!-30.            
                  endif               
               enddo                           
            enddo !j
         enddo !i
c         
      return
      end subroutine ross
c***********************************************************************
      subroutine interp(nel, kk, rossl, gaml, xi, eta, g, i3, f, 
     : zet, zetb, zetx, zety, ux, uy, gx, gy, rosslp, gp, gamlp, fp, fx1, fy1)
      implicit none
      integer, intent(in) :: nel, kk, i3
      real, intent(in) :: ux, uy, gaml(4, 4, nrad), zet(ipe, 4, 4), rossl(4, 4), eta, 
     : rosslp(4, 4, nrad), gamlp(4, 4, nrad)
      real, intent(out) :: gx, gy, f(nrad), zetb(ipe), g, gp(nrad), fp(nrad), fx1(nrad), fy1(nrad),
     : zetx(ipe), zety(ipe)
      integer :: l, i, j, k2
      real ::  V(4), U(4), vyi(4), xi, x, fx(4, 4), fy(4, 4), fxy(4, 4)
c      
c     Interpolation of zet (=mean ionic charge)
      do l = 1, nel
         do i = 1, 4
            do j = 1, 4
               u(j) = zet(l, i, j)
            enddo
            v(i) = fint(u, eta)
            vyi(i) = fintp(u, eta)
         enddo
         zetb(l) = fint(v, xi)
         zety(l) = fint(vyi, xi)
         zetx(l) = fintp(v, xi)
         zety(l) = zety(l)/uy
         zetx(l) = (80./real(i3))*(zetx(l)-zety(l)*ux) 
      enddo
c        
c     interpolation of g (=rosseland mean opacity)
      DO I = 1, 4
         DO J = 1, 4
            U(J) = rossl(I, J)
         ENDDO
         V(I) = FINT(U, ETA)
         vyi(i) = fintp(u, eta)
      ENDDO
      g = FINT(V, XI)
      gy = fint(vyi, xi)
      gx = fintp(v, xi)
      gy = gy/uy
      gx = (80./real(i3))*(gx-gy*ux)     
c
      do k2 = 1, kk
c     interpolation of gp
c     gp=[d g]/[d log10(chi)]
         do I = 1, 4
	         do J = 1, 4
	            U(J)=rosslp(I, J, k2)
	         ENDDO
	         V(I) = FINT(U, ETA)
	      ENDDO
	      gp(k2)=FINT(V, XI)     
c
c     Interpolation of gam (=dimensionless parameter in g_rad)
c     f=log10(gamma)
         DO I = 1, 4
            DO J = 1, 4
! HH: This gives irregularities, perhaps it is preferable to assign nonnegative values of 
! neighbouring interpolation points (in subroutine "ross") ?   
!               if(gaml(i, j, k2).eq.-30.)then
!                  f(k2) = -30.
!                  fp(k2) = -30.
!                  fx1(k2) = -30.
!                  fy1(k2) = -30.
!                  write(6,*) 'Warning: negative g_rad for element ', k2 
!                  goto 1
!               endif      
               U(J) = gaml(I, J, k2)
            ENDDO
            V(I) = FINT(U, ETA)
            vyi(i) = fintp(u, eta)
         ENDDO  
         F(k2) = FINT(V, XI)
         fy1(k2) = fint(vyi, xi)
         fx1(k2) = fintp(v, xi)
         fy1(k2) = fy1(k2)/uy
         fx1(k2) = (80./real(i3))*(fx1(k2)-fy1(k2)*ux)
c      
c     Interpolation of gamp  
c     fp=[d f]/[d log10(chi)]
         do i = 1, 4
            do j = 1, 4
! HH: This gives irregularities, it is preferable to assign nonnegative values of 
! neighbouring interpolation points, see subroutine "ross".           
               U(J) = gamlp(i, j, k2)
            enddo
            V(i) = FINT(U, ETA)
         enddo  
         fp(k2) = FINT(V, XI)          
    1 continue                     
      enddo !k2           
c
      return
      end subroutine interp
C**************************************
      function fint(u,r)
      dimension u(4)
c
c  If  P(R) =   u(1)  u(2)  u(3)  u(4)
c  for   R  =    -3    -1     1     3
c  then a cubic fit is:
      P(R)=( 
     +  27*(u(3)+u(2))-3*(u(1)+u(4)) +R*(
     +  27*(u(3)-u(2))-(u(4)-u(1))   +R*(
     +  -3*(u(2)+u(3))+3*(u(4)+u(1)) +R*(
     +  -3*(u(3)-u(2))+(u(4)-u(1)) ))))/48.
c
        fint=p(r)
c
      return
      end function fint
c***********************************************************************
      function fintp(u,r)
      dimension u(4)
c
c  If  P(R) =   u(1)  u(2)  u(3)  u(4)
c  for   R  =    -3    -1     1     3
c  then a cubic fit to the derivative is:
      PP(R)=( 
     +  27*(u(3)-u(2))-(u(4)-u(1))   +2.*R*(
     +  -3*(u(2)+u(3))+3*(u(4)+u(1)) +3.*R*(
     +  -3*(u(3)-u(2))+(u(4)-u(1)) )))/48.
c
        fintp=pp(r)
c
      return
      end function fintp
c***********************************************************************
      subroutine solve(u,v,z,uz,ierr)
      integer, intent(inout) :: ierr
      dimension u(4)
c
c  If  P(R) =   u(1)  u(2)  u(3)  u(4)
c  for   R  =    -3    -1    1     3
c  then a cubic fit is:
      P(R)=( 
     +  27*(u(3)+u(2))-3*(u(1)+u(4)) +R*(
     +  27*(u(3)-u(2))-(u(4)-u(1))   +R*(
     +  -3*(u(2)+u(3))+3*(u(4)+u(1)) +R*(
     +  -3*(u(3)-u(2))+(u(4)-u(1)) ))))/48.
c  First derivative is:
      PP(R)=( 
     +  27*(u(3)-u(2))-(u(4)-u(1))+ 2*R*(
     +  -3*(u(2)+u(3))+3*(u(4)+u(1)) +3*R*(
     +  -3*(u(3)-u(2))+(u(4)-u(1)) )))/48.
c
!      ierr = 0
c  Find value of z giving P(R)=v
c  First estimate
      z=(2.*v-u(3)-u(2))/(u(3)-u(2))
c  Newton-Raphson iterations
      do k=1,10
         uz=pp(z)
         d=(v-p(z))/uz
         z=z+d
         if(abs(d).lt.1.e-4)return
      enddo
c      
      print*,' Not converged after 10 iterations in SOLVE'
      print*,' v=',v
      DO N=1,4
         PRINT*,' N, U(N)=',N,U(N)
      ENDDO  
      ierr = 10
      return
!      stop
c      
      end subroutine solve
c***********************************************************************
      subroutine scatt(ih,jh,rion,uf,f,umesh,dscat,ntot,epa, ierr)
      integer, intent(inout) :: ierr
      dimension rion(28,4,4),uf(0:100),f(nptot,4,4),umesh(nptot),
     +  fscat(0:100),p(nptot),rr(28),ih(4),jh(4)      
        integer i,j,k,n
c HH: always use meshtype q='m'
      ite3=2
      umin=umesh(1)
        CSCAT=EPA*2.37567E-8
c
        do i=1,4
        do j=1,4
!          if(.not.x(i,j))then
            ft=10**(ITE3*real(ih(i))/40.)
            fne=10**(ITE3*real(jh(j))/4.)
            do k=1,ntot
              p(k)=f(k,i,j)
            enddo
            do m=1,28
              rr(m)=rion(m,i,j)
            enddo      
            CALL BRCKR(FT,FNE,RR,28,UF,100,FSCAT, ierr)
            do n=0,100
              u=uf(n)
              fscat(n)=cscat*(fscat(n)-1)
            enddo
            do n=2,ntot-1
              u=umesh(n)
              if(u.lt.0.01)then      
                  se=u*(1.-.5*u)
              else
                se=1.-exp(-u)
              endif
              m=(u-umin)/dscat
              ua=umin+dscat*m
              ub=ua+dscat
              p(n)=p(n)+((ub-u)*fscat(m)+(u-ua)*fscat(m+1))/(dscat*se)
            enddo
            u=umesh(ntot)
            p(ntot)=p(ntot)+fscat(100)/(1-exp(-u))
            p(1)=p(1)+fscat(1)/(1.-.5*umin)
            do k=1,ntot
              f(k,i,j)=p(k)
            enddo  
!          endif
        enddo
      enddo  
C
      return
      end subroutine scatt
c***********************************************************************
      SUBROUTINE BRCKR(T,FNE,RION,NION,U,NFREQ,SF, ierr)
      integer, intent(inout) :: ierr
C
C  CODE FOR COLLECTIVE EFFECTS ON THOMSON SCATTERING.
C  METHOD OF D.B. BOERCKER, AP. J., 316, L98, 1987.
C
C  INPUT:-
C     T=TEMPERATURTE IN K
C     FNE=ELECTRON DENSITY IN CM**(-3)
C     ARRAY RION (DIMENSIONED FOR 30 IONS).
C        RION(IZ) IS NUMBER OF IONS WITH NET CHARGE IZ.
C        NORMALISATION OF RION IS OF NO CONSEQUENCE.
C     NION=NUMBER OF IONS INCLUDED.
C     ARRAY U (DIMENSIONED FOR 1000). VALUES OF (H*NU/K*T).
C     NFREQ=NUMBER OF FREQUENCY POINTS.
C
C  OUTPUT:-
C     ARRAY SF, GIVING FACTORS BY WHICH THOMSON CROSS SECTION
C     SHOULD BE MULTIPLIED TO ALLOW FOR COLLECTIVE EFFECTS.
C
C  MODIFFICATIONS:-
C     (1) REPLACE (1.-Y) BY EXP(-Y) TO AVOID NEGATIVE FACTORS FOR
C         HIGHLY-DEGENERATE CASES.
C     (2) INCLUDE RELATIVISTIC CORRECTION.
C
      PARAMETER (IPZ=28,IPNC=100)
      DIMENSION RION(IPZ),U(0:IPNC),SF(0:IPNC)
C
      AUNE=1.48185E-25*FNE
      AUT=3.16668E-6*T
      C1=-1.0650E-4*AUT
      C2=+1.4746E-8*AUT**2
      C3=-2.0084E-12*AUT**3
      V=7.8748*AUNE/(AUT*SQRT(AUT))
      CALL FDETA(V,ETA, ierr)  ! 23.10.93
      W=EXP(ETA)         ! 23.10.93
   11 R=FMH(W)/V
      A=0.
      B=0.
      DO 20 I=1,NION
         A=A+I*RION(I)
         B=B+I**2*RION(I)
   20 CONTINUE
      X=R+B/A
C
      Y=.353553*W
      C=1.1799E5*X*(AUNE/AUT**3)
      DO 30 N=0,NFREQ
         D=C/U(N)**2
         IF(D.GT.5.)THEN
            D=-2./D
            F=2.666667*(1.+D*(.7+D*(.55+.341*D)))
         ELSE
            G=2.*D*(1+D)
            F=D*((G+D**3)*LOG(D/(2.+D))+G+2.6666667)
         ENDIF
         DELTA=.375*R*F/X
         SF(N)=(1.-R*DELTA-Y*FUNS(W))*
     +   (1.+U(N)*(C1+U(N)*(C2+U(N)*C3)))   !SAMPSON CORRECTION
   30 CONTINUE
C
      RETURN
C
  600 FORMAT(5X,'NOT CONVERGED IN LOOP 10 OF BRCKR'/
     +       5X,'T=',1P,E10.2,', FNE=',E10.2)
C
      END SUBROUTINE BRCKR
C***********************************************************************
      FUNCTION FUNS(A)
C
      IF(A.LE.0.001)THEN
         FUNS=1.
      ELSEIF(A.LE.0.01)THEN
         FUNS=(1.+A*(-1.0886+A*(1.06066+A*1.101193)))/
     +     (1.+A*(0.35355+A*(0.19245+A+0.125)))
      ELSE
         FUNS=(  1./(1.+0.81230*A)**2+
     +        0.92007/(1.+0.31754*A)**2+
     +        0.05683/(1.+0.04307*A)**2 )/
     +     (  1./(1.+0.65983*A)+
     +        0.92007/(1.+0.10083*A)+
     +        0.05683/(1.+0.00186*A)    )
      ENDIF
      RETURN
      END FUNCTION FUNS
C***********************************************************************
      FUNCTION FMH(W)
C
C  CALCULATES FD INTERGAL I_(-1/2)(ETA). INCLUDES FACTOR 1/GAMMA(1/2).
C  ETA=LOG(W)
C
      IF(W.LE.2.718282)THEN
         FMH=W*(1+W*(-.7070545+W*(-.3394862-W*6.923481E-4))
     +   /(1.+W*(1.2958546+W*.35469431)))
      ELSEIF(W.LE.54.59815)THEN
         X=LOG(W)
         FMH=(.6652309+X*(.7528360+X*.6494319))
     +   /(1.+X*(.8975007+X*.1153824))
      ELSE
         X=LOG(W)
         Y=1./X**2
         FMH=SQRT(X)*(1.1283792+(Y*(-.4597911+Y*(2.286168-Y*183.6074)))
     +   /(1.+Y*(-10.867628+Y*384.61501)))
      ENDIF
C
      RETURN
      END FUNCTION FMH
C***********************************************************************
      SUBROUTINE FDETA(X,ETA, ierr)
C
C  GIVEN X=N_e/P_e, CALCULATES FERMI-DIRAC ETA
C  USE CHEBYSHEV FITS OF W.J. CODY AND H.C. THACHER,
C  MATHS. OF COMP., 21, 30, 1967.
C
      integer, intent(inout) :: ierr
      DIMENSION D(2:12)
      DATA D/
     +  3.5355339E-01, 5.7549910E-02, 5.7639604E-03, 4.0194942E-04,
     +  2.0981899E-05, 8.6021311E-07, 2.8647149E-08, 7.9528315E-10,
     +  1.8774422E-11, 3.8247505E-13, 6.8427624E-15/
C
      integer n,k
c
!      ierr = 0 
      a=x*0.88622693
c
      IF(X.LT.1)THEN
         v=x
         S=V
         U=V
         DO 10 N=2,12
             S=S*V
            SS=S*D(N)
            U=U+SS
            IF(ABS(SS).LT.1.E-6*U)GOTO 11
   10    CONTINUE
         PRINT*,' COMPLETED LOOP 10 IN FDETA'
         ierr = 11
         return
!         STOP
   11    ETA=LOG(U)
c
      ELSE
         if(a.lt.2)then      
            E=LOG(X)
         else
            e=(1.5*a)**0.667
         endif
         do 20 k=1,10
            CALL FDF1F2(E,F1,F2)
            DE=(A-F2)*2./F1
            E=E+DE
            if(abs(dE).lt.1.e-4*abs(E))goto 21
   20    continue
         print*,' completed loop 20 IN FDETA'
         ierr = 12
         return
!         stop
   21    ETA=E
c
      ENDIF
C
      RETURN
      END SUBROUTINE FDETA
C***********************************************************************
      SUBROUTINE FDF1F2(ETA,F1,F2)
C
C  CALCULATES FD INTEGRALS F1, F2=F(-1/2), F(+1/2)
C  USE CHEBYSHEV FITS OF W.J. CODY AND H.C. THACHER,
C  MATHS. OF COMP., 21, 30, 1967.
C
      IF(ETA.LE.1)THEN
         X=EXP(ETA)
         F1=X*(1.772454+X*(-1.2532215+X*(-0.60172359-X*0.0012271551))/
     +      (1.+X*(1.2958546+X*0.35469431)))
         F2=X*(0.88622693+X*(-0.31329180+X*(-0.14275695-
     +      X*0.0010090890))/
     +      (1.+X*(0.99882853+X*0.19716967)))
      ELSEIF(ETA.LE.4)THEN
         X=ETA
         F1=(1.17909+X*(1.334367+X*1.151088))/
     +      (1.+X*(0.8975007+X*0.1153824))
         F2=(0.6943274+X*(0.4918855+X*0.214556))/
     +      (1.+X*(-0.0005456214+X*0.003648789))
      ELSE
         X=SQRT(ETA)
         Y=1./ETA**2
         F1=X*(2.+Y*(-0.81495847+Y*(4.0521266-Y*325.43565))/
     +      (1.+Y*(-10.867628+Y*384.61501)))
         F2=ETA*X*(0.666666667+Y*(0.822713535+Y*(5.27498049+
     +      Y*290.433403))/
     +      (1.+Y*(5.69335697+Y*322.149800)))
      ENDIF
C
      RETURN
      END SUBROUTINE FDF1F2
c***********************************************************************
        SUBROUTINE IMESH(UMESH,NTOT)
C      
      DIMENSION UMESH(nptot)
      COMMON/CIMESH/U(100),AA(nptot),BB(nptot),IN(nptot),ITOT,NN
      save /cimesh/
      
      UMIN=UMESH(1)
      UMAX=UMESH(NTOT)
c      
      II=100
      A=(II*UMIN-UMAX)/REAL(II-1)
      B=(UMAX-UMIN)/REAL(II-1)
      DO I=1,II
        U(I)=A+B*I
      ENDDO  
c
      ib=2
      ub=u(ib)
      ua=u(ib-1)
      d=ub-ua
      do n=2,ntot
        if(umesh(n).gt.ub)then
          ua=ub
          ib=ib+1
          ub=u(ib)
          d=ub-ua
          if(umesh(n).gt.ub)then
            nn=n-1
            ibb=ib-1
            goto 1
          endif 
        endif  
        in(n)=ib
        aa(n)=(ub-umesh(n))/d
        bb(n)=(umesh(n)-ua)/d
      enddo
c
    1      ib=ibb
      do n=nn+1,ntot
        ib=ib+1
        in(n)=ib
        u(ib)=umesh(n)
      enddo  
      itot=ib
c      
        return
      end SUBROUTINE IMESH
c***********************************************************************      
      subroutine screen1(ih,jh,rion,umesh,ntot,epa,f)
      dimension uf(0:100),umesh(nptot),
     +  fscat(0:100),ih(4),jh(4)    
      real, target :: f(nptot,4,4), rion(28,4,4)   
      integer i,j,k,m
      real, pointer :: p(:), rr(:)
c
      ite3=2
      umin=umesh(1)
      umax=umesh(ntot)
c
        do i=1,4
        do j=1,4
            ft=10**(ITE3*real(ih(i))/40.)
            fne=10**(ITE3*real(jh(j))/4.)
!            do k=1,ntot
              p => f(1:ntot,i,j)
!            enddo
!            do m=1,28
              rr => rion(1:28,i,j)
!            enddo      
            call screen2(ft,fne,rr,epa,ntot,umin,umax,umesh,p)
!            do k=1,ntot
!              f(k,i,j)=p(k)
!            enddo            
        enddo
      enddo  
C
      return
      end subroutine screen1
c************************************************************************
      subroutine screen2(ft,fne,rion,epa,ntot,umin,umax,umesh,p)
      parameter(ipz=28)
      dimension rion(ipz),p(nptot),f(100),umesh(nptot)
      dimension x(3),wt(3)
      data x/0.415775,2.294280,6.289945/
      data wt/0.711093,0.278518,0.0103893/
      data twopi/6.283185/
      COMMON/CIMESH/U(100),AA(nptot),BB(nptot),IN(nptot),ITOT,NN
      save /cimesh/
c      
      rydt=ft/157894.
      aune=1.48185e-25*fne
c      
c       get alp2=1/(Debye)**2      
      b=0
      do i=1,ipz
        b=b+rion(i)*i**2
      enddo
        alp2=(5.8804e-19)*fne*b/(epa*ft)
      if(alp2/ft.lt.5e-8)return !!!!!!!!!!!
c      
      c=1.7337*aune/sqrt(rydt)
c      
      do i=1,itot
        w=u(i)*rydt
        f(i)=0.
        do 1 k=1,ipz
          if(rion(k).le.0.01)goto 1
          crz=c*rion(k)*k**2
          ff=0
          do j=1,3
            e=x(j)*rydt
            fk=sqrt(e)
            fkp=sqrt(e+w)
            x1=1.+alp2/(fkp+fk)**2
            x2=1.+alp2/(fkp-fk)**2
            q=(1./x2-1./x1+log(x1/x2))*
     +        (fkp*(1.-exp(-twopi*k/fkp)))/(fk*(1.-exp(-twopi*k/fk)))
              ff=ff+wt(j)*q
          enddo
          f(i)=f(i)+crz*ff
    1     continue
        enddo
c  
      p(1)=f(1)
      do n=2,nn
        w=umesh(n)*rydt
        p(n)=p(n)+(aa(n)*f(in(n)-1)+bb(n)*f(in(n)))/w**3
      enddo
      do n=nn+1,ntot
        w=umesh(n)*rydt
        p(n)=p(n)+f(in(n))/w**3
      enddo  
c
      return
      end subroutine screen2        
c***********************************************************************
      end module opax
