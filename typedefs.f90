  module typedefs
  implicit none
  integer, parameter :: i1b=selected_int_kind(2)
  integer, parameter :: i2b=selected_int_kind(4)
  integer, parameter :: i4b=selected_int_kind(8)
  integer, parameter :: sp=kind(1.0)
  integer, parameter :: spc=kind((1.0,1.0))
  integer, parameter :: dp=selected_real_kind(15)
  integer, parameter :: qp=selected_real_kind(2*precision(1.0_dp))
  integer, parameter :: dpc=kind((1.0_dp,1.0_dp))
  integer, parameter :: qpc=kind((1.0_qp,1.0_qp))
#ifdef DP
  integer, parameter :: r_type = dp 
  integer, parameter :: q_type = qp
  integer, parameter :: c_type = dpc 
  integer, parameter :: qc_type = qpc 
#else
  integer, parameter :: r_type = sp 
  integer, parameter :: c_type = spc 
#endif
end module typedefs
