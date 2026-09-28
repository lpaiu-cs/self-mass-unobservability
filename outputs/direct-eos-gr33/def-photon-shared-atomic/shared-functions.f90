logical function shared_has(ion)
  integer,intent(in) :: ion
  shared_has=shared_start(ion)>0
end function

subroutine shared_atomic(ion,tl,x,q,qt,qx,qtt,qtx,qxx)
  use mod_ionization_data, only: nion
  use mod_statistical_weight_data, only: iqion
  integer,intent(in) :: ion
  real(fp_kind),intent(in) :: tl,x(5)
  real(fp_kind),intent(out) :: q,qt,qx(5),qtt,qtx(5),qxx(5,5)
  integer :: lo,hi
  lo=shared_start(ion);hi=shared_end(ion)
  if(lo<1.or.hi<lo) error stop 'shared atomic catalog missing'
  qx=0; qtx=0; qxx=0
  ! Finite imported levels only. Native radius approximation ignores ell.
  call qmhd_calc(1,1._fp_kind,.true.,.false.,1._fp_kind,11,10,tl,nion(ion),1._fp_kind, &
       shared_weight(lo:hi),shared_neff(lo:hi),0._fp_kind*shared_neff(lo:hi),x,q,qt,qx,qtt,qtx,qxx)
  q=q/iqion(ion);qt=qt/iqion(ion);qx=qx/iqion(ion)
  qtt=qtt/iqion(ion);qtx=qtx/iqion(ion);qxx=qxx/iqion(ion)
end subroutine

subroutine shared_terms(ion,tl,x,values) bind(C)
  use iso_c_binding
  use mod_ionization_data, only: nion
  integer(c_int),value :: ion
  real(c_double),value :: tl
  real(c_double),intent(in) :: x(5)
  real(c_double),intent(out) :: values(38,700)
  real(fp_kind) :: q,qt,qtt,qx(5),qtx(5),qxx(5,5)
  integer :: j,k
  values=0
  if(.not.shared_has(ion)) error stop 'shared optical catalog missing'
  do j=shared_start(ion),shared_end(ion)
     qx=0;qtx=0;qxx=0
     call qmhd_calc(1,1._fp_kind,.true.,.false.,1._fp_kind,11,10,tl,nion(ion),1._fp_kind, &
          shared_weight(j:j),shared_neff(j:j),0._fp_kind*shared_neff(j:j),x,q,qt,qx,qtt,qtx,qxx)
     k=j-shared_start(ion)+1
     values(:,k)=[q,qt,qtt,qx,qtx,reshape(qxx,[25])]
  end do
end subroutine
