module datatype
  use omp_lib
 implicit none
 integer, parameter, public :: r1=selected_real_kind(1)
 integer, parameter, public :: r2=selected_real_kind(9)
end module datatype