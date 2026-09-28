program export_constants
  use mod_free_eos_constants, only: ergspercmm1,rydberg,c2,clight,boltzmann
  use mod_statistical_weight_data, only: iqneutral,iqion
  use mod_ionization_data, only: bi
  implicit none
  print *,ergspercmm1,rydberg,c2,clight,boltzmann
  print *,iqneutral
  print *,iqion
  print *,bi
end program
