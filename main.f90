PROGRAM main
!==============================================================================
! MAIN PROGRAM
!
! Reads parameters from the command line, allocates the initial state, and runs
! the waveguide spontaneous-emission simulation via output:printTrajectory_HL.
!
! Typical run:
!   ./mpol -prep 2 -del 0.1 -al 1.0 -nm 4000 -dt 0.02 -tmax 200 -npini 2 -npadd 4
!==============================================================================


  USE system
  USE output

  IMPLICIT NONE

  TYPE(param) 	 					::  sys,sys2
  type(state)						::  st, ini_st,wp_st
  real(rl)							::  t1,t2

  print*, "-- Starting parameter initialisation - MAIN"
  CALL getParameters(sys)
  print*, "-- Parameters initialised - MAIN"

  CALL printTrajectory_HL(sys,st)

  IF ( (sys%prep > 99) .and. ( sys%prep < 110 ) ) then

	 sys2=sys
	 sys2%npini = sys%npini + sys%npadd
	 sys%prep = sys%prep - 2600
	 CALL allocate_state(sys2,st)

  END IF

END PROGRAM main







