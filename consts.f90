! ==============================================================================
!  MODULE: consts.f90
!
!  PURPOSE & CONTEXT
!    Central repository for shared mathematical constants, physical constants, 
!    and low-level safety utilities for the waveguide QED simulation.
!
!  PHYSICS CONTEXT
!    All physics within the simulation is expressed in natural units where 
!    Planck's constant \hbar = 1. Consequently, energy and frequency share 
!    the same dimensions. 
!
!  CORE RESPONSIBILITIES
!    1. Constants     : Defines precision-independent versions of Pi, 0, 1, 
!                       and the imaginary unit `i` to prevent type-mismatch 
!                       errors during complex tensor algebra.
!    2. Memory Safety : Provides `alloc_check` to gracefully catch and report 
!                       RAM allocation failures before causing a segmentation fault.
! ==============================================================================
MODULE consts

  use typedefs 
  implicit none

  ! ----------------------------------------------------------------------------
  ! 1. REAL CONSTANTS (Standard / Double Precision via r_type)
  ! ----------------------------------------------------------------------------
  real(r_type),   parameter :: c_light = 1_r_type
  real(r_type),   parameter :: pi=(4.0_r_type)*atan(1.0_r_type)
  real(r_type),   parameter :: one_r=1.0_r_type
  real(r_type),   parameter :: zero_r=0.0_r_type

  ! ----------------------------------------------------------------------------
  ! 2. QUADRUPLE PRECISION CONSTANTS (Extended Precision via q_type)
  ! ----------------------------------------------------------------------------
  real(q_type),   parameter :: pi_q=(4.0_q_type)*atan(1.0_q_type)
  real(q_type),   parameter :: one_q=1.0_q_type
  real(q_type),   parameter :: zero_q=0.0_q_type

  ! ----------------------------------------------------------------------------
  ! 3. COMPLEX CONSTANTS (Double Complex via c_type)
  ! ----------------------------------------------------------------------------
  complex(c_type),parameter :: im=(0.0d0,1.0d0)
  complex(c_type),parameter :: zero=(0.0d0,0.0d0)
  complex(c_type),parameter :: one=(1.0d0,0.0d0),two=(2.0d0,0.0d0)

contains

  !> -------------------------------------------------------------------------
  !> SUBROUTINE: alloc_check
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Defensive programming utility for safe dynamic memory allocation. 
  !>   It checks the status flag returned by Fortran's `allocate(..., stat=X)` 
  !>   statements. If the operating system fails to grant memory (e.g., due 
  !>   to the exponentially growing basis size), it catches the error, prints 
  !>   the specific array name, and safely halts execution to prevent a hard 
  !>   segmentation fault.
  !> Arguments:
  !>   - stat : The integer status flag returned by the `allocate` call.
  !>   - str  : A string identifier indicating which array failed to allocate.
  !>
 subroutine alloc_check(stat,str)
    implicit none 
    integer           :: stat
    character(len=50) :: str

    if(stat.ne.0)then 
       print*, ': Allocation faliure of the matrix... quitting', &
            ' - ', str
       stop  
    end if
  end subroutine alloc_check

END MODULE consts
