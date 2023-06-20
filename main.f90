PROGRAM main

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







