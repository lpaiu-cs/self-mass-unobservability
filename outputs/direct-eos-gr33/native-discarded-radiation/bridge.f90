subroutine hydrogen_cross(n,z,level,t,fr,bf,ff) bind(C)
 use iso_c_binding
 implicit none
 integer(c_int),value::n,z,level
 real(c_double),value::t
 real(c_double),intent(in)::fr(n)
 real(c_double),intent(out)::bf(n),ff(n)
 real(c_double),external::gaunt,gfree
 integer::i
 do i=1,n
  bf(i)=2.815d29*real(z,c_double)**4/real(level,c_double)**5/fr(i)**3*gaunt(level,fr(i)/real(z*z,c_double))
  ff(i)=gfree(t,fr(i)/real(z*z,c_double))
 enddo
end subroutine
subroutine hydrogen_line(i,j,z,f) bind(C)
 use iso_c_binding
 implicit none
 integer(c_int),value::i,j,z
 real(c_double),intent(out)::f
 real(c_double)::x,w,f1,f0
 call stark0(i,j,z,x,w,f1,f0)
 f=f1+f0
end subroutine
