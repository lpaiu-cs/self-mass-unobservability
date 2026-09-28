subroutine hydrogen_rr(t,r) bind(C)
 use iso_c_binding
 implicit none
 real(c_double),value::t
 real(c_double),intent(out)::r
 real::a,b
 a=real(t)
 call rrfit(1,1,a,b)
 r=real(b,c_double)
end subroutine
