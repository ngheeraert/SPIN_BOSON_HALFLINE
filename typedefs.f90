! ==============================================================================
!  MODULE: typedefs.f90
!
!  PURPOSE & CONTEXT
!    Foundation layer defining the numerical precision (Kinds)
!
! ==============================================================================

MODULE typedefs

  implicit none

  ! ----------------------------------------------------------------------------
  ! 1. INTEGER KINDS
  ! ----------------------------------------------------------------------------
  integer, parameter :: i1b=selected_int_kind(2)	! 1-byte integer (up to 10^2)
  integer, parameter :: i2b=selected_int_kind(4)	! 2-byte integer (up to 10^4)
  integer, parameter :: i4b=selected_int_kind(8)	! 4-byte integer (up to 10^8)

  ! ----------------------------------------------------------------------------
  ! 2. FLOATING POINT & COMPLEX KINDS (Base Definitions)
  ! ----------------------------------------------------------------------------
  integer, parameter :: sp=kind(1.0)			! 32-bit Single Precision
  integer, parameter :: spc=kind((1.0,1.0))		! 32-bit Single Complex
  integer, parameter :: dp=selected_real_kind(15)	! 64-bit Double Precision (15 decimal digits)
  integer, parameter :: qp=selected_real_kind(2*precision(1.0_dp)) ! 64-bit Double Complex
  integer, parameter :: dpc=kind((1.0_dp,1.0_dp))	! 128-bit Quadruple Precision
  integer, parameter :: qpc=kind((1.0_qp,1.0_qp))	! 128-bit Quadruple Complex

  ! ----------------------------------------------------------------------------
  ! 3. GLOBAL PRECISION ALIASES (Compiler Preprocessor)
  ! ----------------------------------------------------------------------------
  ! If compiled with the flag `-DDP` (e.g., gfortran -cpp -DDP ...), the 
  ! simulation defaults to highly stable Double Precision. Otherwise, it 
  ! compiles in Single Precision (faster, but highly prone to round-off errors 
  ! during the TDVP matrix inversions).

#ifdef DP
! -- Double Precision Configuration --
  integer, parameter :: r_type = dp 
  integer, parameter :: q_type = qp
  integer, parameter :: c_type = dpc 
  integer, parameter :: qc_type = qpc 
#else
! -- Single Precision Configuration --
  integer, parameter :: r_type = sp 
  integer, parameter :: c_type = spc 
#endif
END MODULE typedefs
