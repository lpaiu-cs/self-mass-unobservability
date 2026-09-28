subroutine exchange_cached(res) bind(C)
  use iso_c_binding
  use mod_free_eos_constants, only: boltzmann,cpe
  use mod_master_exchange_data, only: iforder,fexprime,n_e,t,p_e,pstarprime, &
       psiprime,psi,muex,muex2,muexf,muex2f
  implicit none
  real(c_double), intent(out) :: res(5)
  res=[fexprime+boltzmann*t*(psiprime-psi)*n_e-(cpe*pstarprime(1)-p_e), &
       muex+muex2,muexf+muex2f,n_e,real(iforder,c_double)]
end subroutine exchange_cached

subroutine exchange_probe(fl,tl,res) bind(C)
  use iso_c_binding
  use mod_free_eos_constants, only: boltzmann,c_e
  use mod_exchange, only: master_exchange,exchange_free
  use mod_master_exchange_data, only: iforder,psi,psiprime,dpsiprimedf,flprimef
  implicit none
  real(c_double), value :: fl,tl
  real(c_double), intent(out) :: res(16)
  real(c_double) :: r(9),p(9),ss(3),uu(3),dve,dvef,dvet,free,freef,ne,kt,psif
  call master_exchange(0,fl,tl,r,p,ss,uu,21,14,dve,dvef,dvet)
  call exchange_free(r,p,free,freef)
  ne=c_e*r(1);kt=boltzmann*exp(tl);psif=sqrt(1.d0+exp(fl))
  res=[free,-dve,-dvef,-dvet,ne,r(2),kt,freef, &
       -dvef/r(2),(psif-dvef)/r(2),psi,psiprime,dpsiprimedf*flprimef, &
       -dve-(psiprime-psi),psif-dvef-dpsiprimedf*flprimef,real(iforder,c_double)]
end subroutine exchange_probe
