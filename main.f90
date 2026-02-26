! ==============================================================================
!  FILE: main.f90
!
!  PURPOSE & CONTEXT
!    Primary entry point for the Time-Dependent Variational Principle (TDVP) 
!    simulation of spontaneous emission in ultrastrong-coupling waveguide QED.
!    Generates the data underlying the time-domain results in:
!      N. Gheeraert et al., "Spontaneous emission of Schrödinger cats in a 
!      waveguide at ultrastrong coupling", New J. Phys. 19 (2017) 023036.
!
!  CORE RESPONSIBILITIES
!    1. Initialization : Parses command-line arguments to configure the 
!                        Hamiltonian and the multi-polaron basis.
!    2. Orchestration  : Dispatches the constructed parameters to the `output` 
!                        module to execute the target trajectory.
!    3. Post-Processing: Handles conditional flag shifting for multi-stage 
!                        runs (e.g., reloading states for phase-space tomography).
!
!  TYPICAL RUN
!    ./mpol -prep 2 -del 0.1 -al 1.0 -nm 4000 -dt 0.02 -tmax 200 -npini 2 -npadd 4
!
!  EXECUTION HIERARCHY
!      main.f90 
!       ├── system:getParameters (Parses CLI and builds spatial/momentum grids)
!       └── output:printTrajectory_HL (Manages time-stepping and basis expansion)
! ==============================================================================
PROGRAM main

  USE system
  USE output

  IMPLICIT NONE

  TYPE(param) 	 					::  sys,sys2
  type(state)						::  st, ini_st,wp_st
  real(rl)							::  t1,t2

  ! ----------------------------------------------------------------------------
  ! 1. INITIALIZATION
  ! ----------------------------------------------------------------------------
  ! Parse the command-line arguments to build the `sys` (param) object.
  ! This establishes the grid sizes (N_modes), physical couplings (alpha, delta),
  ! and the thresholds for the adaptive basis expansion.
  print*, "-- Starting parameter initialisation - MAIN"
  CALL getParameters(sys)
  print*, "-- Parameters initialised - MAIN"

  ! ----------------------------------------------------------------------------
  ! 2. PRIMARY TIME EVOLUTION TRAJECTORY
  ! ----------------------------------------------------------------------------
  ! Hand control over to the output module's 'Half-Line' (HL) orchestrator.
  ! This routine will handle memory allocation, initial state preparation 
  ! (e.g., pumping the atom to the excited state), RK4 integration, and data I/O.
  CALL printTrajectory_HL(sys,st)

  ! ----------------------------------------------------------------------------
  ! 3. CONDITIONAL POST-PROCESSING / RESTART LOGIC
  ! ----------------------------------------------------------------------------
  ! If a specific `prep` flag block is triggered (e.g., prep values between 
  ! 100 and 109), the program prepares for a secondary processing stage. 
  ! It clones the system parameters, expands the initial basis size to account 
  ! for polarons added dynamically during the primary run, and shifts the 
  ! prep flag (e.g., subtracting 2600) to trigger specific diagnostic routines 
  ! (like high-resolution Wigner tomography) upon reallocation.
  IF ( (sys%prep > 99) .and. ( sys%prep < 110 ) ) then

	 sys2=sys
	 sys2%npini = sys%npini + sys%npadd
	 sys%prep = sys%prep - 2600
	 CALL allocate_state(sys2,st)

  END IF

END PROGRAM main







