program native_constants
  use mod_free_eos_constants, only: c2,rydberg,electron_mass,h_mass
  implicit none
  write(*,'(4es26.17)') c2,rydberg,electron_mass,h_mass
end program
