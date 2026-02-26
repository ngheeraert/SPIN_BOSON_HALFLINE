! ==============================================================================
!  MODULE: system.f90 (SYSTEM)
!
!  PURPOSE & CONTEXT
!    Core physics engine and data structures for the time-dependent variational 
!    simulation of spontaneous emission in waveguide QED at ultrastrong coupling.
!    Generates the data underlying the time-domain results in:
!      N. Gheeraert et al., "Spontaneous emission of Schrödinger cats in a 
!      waveguide at ultrastrong coupling", New J. Phys. 19 (2017) 023036.
!
!  PHYSICAL MODEL
!    Spin-boson model in the continuum limit. A two-level system (qubit) with 
!    tunneling splitting `del` is coupled to a 1D bosonic waveguide.
!      H = (Δ/2)σ_z + Σ_k ω_k a_k^† a_k + (σ_x) Σ_k g_k (a_k + a_k^†)
!    The coupling g_k enforces an Ohmic spectral density with an exponential cutoff.
!
!  CORE RESPONSIBILITIES
!    1. Data Structures : Defines the Multi-Polaron Ansatz parameters (`state`).
!    2. TDVP Evaluator  : Implements `CalcDerivatives`, which translates the 
!                         Dirac-Frenkel variational principle into a solvable 
!                         set of explicit Ordinary Differential Equations (ODEs).
!    3. Tensor Algebra  : Implements the $O(N_{cs}^6)$ flattening trick (detailed 
!                         in the paper's Appendix) to resolve the implicit 
!                         dependence of $\dot{f}$ on $\kappa$.
!    4. Memory Caching  : Precomputes non-orthogonal overlaps to drastically 
!                         reduce integration bottleneck times.
!
!  EXECUTION HIERARCHY (Integration Context)
!      main.f90 
!       ├── system:getParameters (initializes grid and Hamiltonian parameters)
!       └── output:printTrajectory_HL 
!            ├── system:allocate_state / initialise_vaccum (setup)
!            └── output:evolveState_HL 
!                 ├── output:Evolve_RK4 (integrator)
!                 │    ├── system:update_sums (caches macroscopic overlaps)
!                 │    └── system:CalcDerivatives (solves explicit TDVP EOMs)
!                 └── system:Energy, error, sigmaZ, n_up_x (evaluates observables)
! ==============================================================================

MODULE SYSTEM 

	USE consts
	USE lapackmodule
	USE inverse
	USE typedefs, only : cx => c_type, rl => r_type

	IMPLICIT NONE 

!------------------------------------------------------------------------------
! TYPE: param
!   Central configuration object. Stores the discretized grid (k-space), the 
!   Hamiltonian parameters (α, Δ, ω_c), and the dynamic run controls (time steps, 
!   error thresholds for adaptive basis expansion).
!------------------------------------------------------------------------------
	TYPE param
		integer    				 ::  npini
		integer    				 ::  npadd
		real(rl), dimension(60) ::  merr_arr   	!-- time for adding additional polaron
		real(rl), dimension(60) ::  tref_arr   	!-- time for adding additional polaron
		real(rl)    				 ::  p0 			!-- initial value for added polarons 

		real(rl)					 ::  alpha     
		real(rl)   				 ::  del    	!-- delta
		integer    				 ::  nmode  	!-- number of modes
		character(len=4)			 ::  file_path
		real(rl)    				 ::  p1i 		!-- initial value for p(1) 
		real(rl)					 ::  length    !-- length on which to plot the field distribution
		real(rl)					 ::  dx     	!-- distance between 2 sites
		real(rl)   				 ::  dt     	!-- time step
		real(rl)   				 ::  tmax   	!-- final time
		real(rl)   				 ::  tref   	!-- time to add polaron
		real(rl)   				 ::  merr   	!-- time for adding additional polaron
		real(rl)   				 ::  wigxmin  	!-- final time
		real(rl)   				 ::  wigxmax   !-- final time
		real(rl)   				 ::  wc     	!-- frequency cutoff
		real(rl)					 ::  wmax      !-- maximum frequency (cutoff at wc is exponential)
		real(rl), allocatable   ::  w(:)      !-- frequency spectrum
		real(rl)					 ::  dk1
		real(rl), allocatable   ::  dk(:)     !-- spacing between k values
		real(rl), allocatable   ::  g(:)      !-- coupling
		integer						 ::  getWigner
		real(rl)		 			 ::  max_deltaE
		integer						 ::  prep
		integer 					 :: verbose

	END TYPE param

!------------------------------------------------------------------------------
! TYPE: STATE
!   The dynamic Variational Wavefunction |Psi(t)>. 
!   Constructed as a superposition of N_{cs} multimode coherent states.
!     |Psi(t)> = Σ [ p_m(t)|↑>⊗|f_m(t)>  +  q_m(t)|↓>⊗|h_m(t)> ]
!   To avoid O(N_modes * N_cs^3) scaling in the tight RK4 inner loops, macroscopic 
!   overlaps (ov_ff) and energy sums (bigW, bigL) are aggressively cached.
!------------------------------------------------------------------------------
	TYPE STATE
		real(rl)					 	::  t 								!-- Time
		integer    				 	::  np
		complex(cx), allocatable  ::  f(:,:),h(:,:)   				!-- value of the fiel displacements and momenta
		complex(cx), allocatable  ::  p(:),q(:)   					!-- probability amplitudes of the polarons
		complex(cx), allocatable  ::  fdot(:,:),hdot(:,:)   		!-- the time derivatives
		complex(cx), allocatable  ::  pdot(:),qdot(:)   			!-- the time derivatives
		complex(cx), allocatable  ::  ov_ff(:,:), ov_hh(:,:)  	!-- matrices of overlaps
		complex(cx), allocatable  ::  ov_fh(:,:), ov_hf(:,:)  	!-- matrices of overlaps
		complex(cx), allocatable  ::  bigW_f(:,:), bigW_h(:,:)  !-- Ws (as in the notes)
		complex(cx), allocatable  ::  bigL_f(:,:), bigL_h(:,:)  !-- Ls (as in the notes)
	END TYPE STATE

!------------------------------------------------------------------------------
! TYPE: TRAJ
!   Trajectory container. Buffers the time-series arrays for macroscopic 
!   observables before they are flushed to disk to prevent I/O bottlenecking.
!------------------------------------------------------------------------------
	TYPE TRAJ
		integer							::  i_time
		real(rl), allocatable  	::  time_ar(:), energy_ar(:),  &
			norm_ar(:), error_ar(:), spinXYZ_ar(:,:), ps_ar(:,:) !,spinXYZ(:,:)
		integer, allocatable  		::  np_ar(:)
	END TYPE TRAJ

	COMPLEX(cx), PARAMETER :: Ic = ( 0._rl , 1._rl )
	REAL(rl),	   PARAMETER :: IP_p = 0.000005_rl 		!-- initialisation parameter


CONTAINS

  !> -------------------------------------------------------------------------
  !> SUBROUTINE: getParameters
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Parses the command line to populate the `sys` object. Crucially, it 
  !>   initializes the discretized Ohmic bath parameters. The coupling array `g_k` 
  !>   is strictly scaled so that the macroscopic spectral density $J(\omega)$ 
  !>   reproduces the continuum limit: $J(\omega) \propto \alpha \omega e^{-\omega/\omega_c}$.
  !> Arguments:
  !>   - sys : Parameter structure to be populated.
  !>
	SUBROUTINE getParameters(sys)
		type(param)           			  ::  sys
		type(state)           	     	  ::  st
		character(len=100)    			  ::  buffer
		integer 				  			  ::  i, nargs
		character							  ::  path_input

		!== these are the parameters provided in the command line 
		nargs = iargc()
		call get_command_argument(1, buffer)

		sys%prep = 2
		do i=1,nargs
			call get_command_argument(i, buffer)
			if(buffer=='-prep')then 
				call get_command_argument(i+1, buffer)
				read(buffer,*) sys%prep
			end if
		end do

		sys%npini= 1
		sys%npadd= 0
		sys%p0 = IP_p
		sys%merr = 0.000001_rl
		sys%tref = 0.2_rl
		sys%dt = 0.05_rl
		sys%file_path = "data"
		sys%verbose = 1

		!-- for half-chain
		sys%alpha = 0.1_rl
		sys%max_deltaE= 1.0e-7_rl
		sys%del = 0.1_rl
		sys%nmode = 400
		sys%tmax = 200._rl
		sys%wc = 1000._rl
		sys%wmax = 1._rl
		sys%getWigner = 0
		sys%wigxmin = sys%tmax - 200._rl
		sys%wigxmax = sys%tmax + 100._rl


		if(buffer == '--help' .or. buffer == '-h' .or. buffer=='-help')then
			print*, '========================================'
			print*, '-npini [1]','     	', 'Initial number of coherent states'
			print*, '-nm [400]','      	', 'Number of modes'
			print*, '-del [0.1]','     	', 'Delta'
			print*, '-al [0.1]','      	', 'Alpha'
			print*, '-dt [0.05]','      	', 'Time step'
			print*, '-tmax [00]','    		', 'Final time'
			print*, '-wc [1000]','     	', 'Frequency cutoff'
			print*, '-wmax [1]','      	', 'Maximum frequency (smooth cutoff)'
			print*, ' '
			print*, '-- Polaron adding parameters ------------------------------'
			print*, '-npadd [0]','     	', 'Number of coherent states to be added'
			print*, '-tref [0.2]','	     	', 'Minimum time before adding additional CS'
			print*, '-merr [10^-7]','	   ', 'Error threshold at which to add a CS'
			print*, ' '
			print*, '-- State preparation parameters ----------------------------'
			print*, '-prep [2]','			', '1=|down>+|up>, 2=|up>, 3=|down>'
			print*, '========================================'
			stop

		elseif(nargs==0)then
			stop '# Use input parameters parameters'
		else
			do i=1,nargs
				call get_command_argument(i, buffer)
				if(buffer=='-p1i')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%p1i
				else if(buffer=='-p0')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%p0
				else if(buffer=='-merr')then
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%merr
				elseif(buffer=='-npini')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%npini
				elseif(buffer=='-npadd')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%npadd
				else if(buffer=='-nm')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%nmode
				else if(buffer=='-del')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%del
				else if(buffer=='-al')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%alpha
				else if(buffer=='-wc')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%wc
				else if(buffer=='-wmax')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%wmax
				else if(buffer=='-dt')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%dt
				else if(buffer=='-tmax')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%tmax
				else if(buffer=='-tref')then 
				    call get_command_argument(i+1, buffer)
					read(buffer,*) sys%tref
				else if(buffer=='-xmin')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%wigxmin
				else if(buffer=='-xmax')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%wigxmax
				else if(buffer=='-getwig')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%getWigner
				else if(buffer=='-verbose')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%verbose
				else if(buffer=='-path')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) path_input
					if ( path_input == "." ) then
						sys%file_path = "   ."
					else if (path_input == "d") then
						sys%file_path = "data"
					end if

				end if

			end do

		end if

		sys%merr_arr = sys%merr
		sys%tref_arr = sys%tref

		sys%merr_arr(1) = sys%merr/3._rl
		sys%merr_arr(2:) = sys%merr	 
		sys%tref_arr(1) = 0.05_rl
		sys%tref_arr(2:) = sys%tref

		allocate(sys%w(sys%nmode))
		allocate(sys%dk(sys%nmode))
		allocate(sys%g(sys%nmode))

		sys%dk(:) = sys%wmax / dble( sys%nmode )
		sys%dk1 = sys%dk(1)

		do i=1, sys%nmode
			sys%w(i) = sys%dk(1) * (i-0.5_rl)
		end do

		!-- define de coupling g(k)
!       ! Coupling amplitudes g_k (discretized): chosen such that the bath spectral density
!       ! matches an Ohmic form with exponential cutoff:
!       !   J(omega) = pi Sum_k g_k^2 delta(omega-omega_k)  ->  2pialpha omega e^{-omega/omegac}.
!       ! With a linear grid omega_k and spacing dk, this gives roughly g_k^2 proportional to alpha omega_k dk e^{-omega_k/omegac}.
		sys%g(:) = sqrt( 2._rl*sys%alpha * sys%w(:) * sys%dk(:) * exp(-sys%w(:)/sys%wc) ) 

		!-- system variables that depend on a possible input
!       ! Real-space grid length for Fourier transforms (k-space <-> x-space):
!       ! with constant dk, a periodic box length L approx pi/dk is used (convention of this code).
		sys%length = pi/sys%dk1
		sys%dx  = sys%length/sys%nmode
	END SUBROUTINE  

  !> -------------------------------------------------------------------------
  !> SUBROUTINE: allocate_state
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Dynamically allocates memory for the complex arrays constituting the 
  !>   multi-polaron ansatz. It is called both at initialization and dynamically 
  !>   during RK4 execution whenever the Dirac-Frenkel error mandates basis 
  !>   expansion (increasing `np`).
  !> Arguments:
  !>   - sys     : Parameter structure providing grid dimensions.
  !>   - st      : State object to be allocated.
  !>   - npvalue : Optional override for the target basis dimension.
  !>
	SUBROUTINE allocate_state(sys,st,npvalue)
		type(param), intent(in)      	::  sys
		type(state), intent(out)  		::  st
		type(state)							::  tmpst !-- if st has to be reallocated
		integer, intent(in),optional		::  npvalue
		integer									::  np, spinXYZ_dim, step_numb

		if (sys%verbose==1) then
			print*, "-- BEGINNING ARRAYS ALLOCATION"
		end if

		if (present(npvalue)) then
			tmpst%np = npvalue
		else
			tmpst%np = sys%npini
		end if
		np = tmpst%np

		!spinXYZ_dim = (sys%tmax/sys%dt)*2._rl
		!allocate(tmpst%spinXYZ(spinXYZ_dim,7))

		allocate(tmpst%f(np,sys%nmode))
		allocate(tmpst%h(np,sys%nmode))
		allocate(tmpst%p(np))
		allocate(tmpst%q(np))
		allocate(tmpst%fdot(np,sys%nmode))
		allocate(tmpst%hdot(np,sys%nmode))
		allocate(tmpst%pdot(np))
		allocate(tmpst%qdot(np))
		allocate(tmpst%ov_ff(np,np))
		allocate(tmpst%ov_hh(np,np)) 
		allocate(tmpst%ov_fh(np,np)) 
		allocate(tmpst%ov_hf(np,np)) 
		allocate(tmpst%bigW_f(np,np)) 
		allocate(tmpst%bigW_h(np,np)) 
		allocate(tmpst%bigL_f(np,np)) 
		allocate(tmpst%bigL_h(np,np))

		tmpst%t = 0._rl
		!tmpst%spinXYZ(:,:) = 0._rl
		tmpst%f(:,:) = 0._cx
		tmpst%h(:,:) = 0._cx
		tmpst%p(:) = 0._cx
		tmpst%q(:) = 0._cx
		tmpst%fdot(:,:) = 0.0_rl
		tmpst%hdot(:,:) = 0.0_rl
		tmpst%pdot(:) = 0.0_rl
		tmpst%qdot(:) = 0.0_rl
		tmpst%ov_ff(:,:) = 0._cx
		tmpst%ov_hh(:,:) = 0._cx
		tmpst%ov_fh(:,:) = 0._cx
		tmpst%ov_hf(:,:) = 0._cx

		st = tmpst

		if (sys%verbose==1) then
			print*, "-- ARRAYS ALLOCATED"
		end if

	END SUBROUTINE

  !> -------------------------------------------------------------------------
  !> SUBROUTINE: allocate_trajectory
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Allocates the contiguous memory buffers for the logging trajectory. 
  !>   Sized slightly larger (1.5x) than the raw `tmax / dt` limit to safely 
  !>   accommodate potential variable-timestep rewinds triggered by `checkTimestep`.
  !> Arguments:
  !>   - sys : Parameter structure.
  !>   - tr  : Trajectory object to be allocated.
  !>
	SUBROUTINE allocate_trajectory(sys,tr)
		type(param), intent(in)      	::  sys
		type(traj), intent(out)      	::  tr
		integer									::  step_numb

		if (sys%verbose==1) then
			print*, "-- BEGINNING TRANJECTORY ALLOCATION"
		end if

		!spinXYZ_dim = (sys%tmax/sys%dt)*2._rl
		!allocate(tmpst%spinXYZ(spinXYZ_dim,7))
		step_numb = int( (sys%tmax / sys%dt)*1.5 ) 

		allocate( tr%time_ar(step_numb) )
		allocate( tr%error_ar(step_numb) )
		allocate( tr%norm_ar(step_numb) )
		allocate( tr%energy_ar(step_numb) )
		allocate( tr%np_ar( step_numb ) )
		allocate( tr%spinXYZ_ar(step_numb,3) )
		allocate( tr%ps_ar(step_numb,sys%npini+sys%npadd) )

		tr%i_time=1
		tr%time_ar = 0._rl
		tr%norm_ar = 0._rl
		tr%error_ar = 0._rl
		tr%energy_ar = 0._rl
		tr%np_ar = 0._rl
		tr%spinXYZ_ar = 0._rl
		tr%ps_ar = 0._rl

		if (sys%verbose==1) then
			print*, "-- TRAJECTORY ALLOCATED"
		end if
	END SUBROUTINE

  !> -------------------------------------------------------------------------
  !> SUBROUTINE: initialise_vaccum
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Constructs the initial $t=0$ boundary condition. The bosonic field is 
  !>   forced strictly to the global vacuum (all displacement arrays zeroed), 
  !>   while the atomic system is pumped into the specified initial configuration 
  !>   (e.g., fully excited `prep=2`) to begin spontaneous emission.
  !> Arguments:
  !>   - sys     : Parameter structure.
  !>   - st      : State object to initialize.
  !>   - time_in : The initialization time (usually 0.0).
  !>
	SUBROUTINE initialise_vaccum(sys,st,time_in)
		type(param), intent(in)      	::  sys
		type(state), intent(in out) 	 	::  st
		real(rl), intent(in)				::  time_in
		integer									::  i,np

		st%t = time_in
		st%f = 0._rl
		st%h = 0._rl

		if ( sys%prep == 1 ) then
			!-- preapre 1st excited state: up + down
			st%p(1) = 1._rl/sqrt(2._rl)
			st%q(1) = 1._rl/sqrt(2._rl)
			if (sys%verbose==1) then
				print*, "=="
				print*, "STATE PREPARED IN UP + DOWN"
				print*, "=="
			end if
		else if ( (sys%prep == 2) ) then
			!-- prepare in the up state
			st%p(1) = 0.99999999999999_rl
			st%q(1) = sqrt( 1._rl - conjg(st%p(1))*st%p(1) )
			if (sys%verbose==1) then
				print*, "=="
				print*, "STATE PREPARED IN UP"
				print*, "=="
			end if
		else if ( sys%prep == 3 ) then
			!-- prepare in the up state
			st%q(1) = 0.99999999999999_rl
			st%p(1) = sqrt( 1._rl - conjg(st%q(1))*st%q(1) )
			if (sys%verbose==1) then
				print*, "=="
				print*, "STATE PREPARED IN DOWN"
				print*, "=="
			end if
		else 
			if (sys%verbose==1) then
				print*, "Error in the value of prep"
			end if
		end if

		!-- updating the sums over k
		CALL update_sums(sys,st)
		CALL normalise(st)
	END SUBROUTINE

  !> -------------------------------------------------------------------------
  !> SUBROUTINE: initialise_from_file
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Post-processing and restart utility. Directly loads a previously evolved 
  !>   multi-polaron state ($\mathbf{f}, \mathbf{h}, p, q$) from the exported 
  !>   ASCII files. This allows the user to re-evaluate expensive spatial 
  !>   diagnostics (like the 2D Wigner tomography) without re-integrating the EOMs.
  !> Arguments:
  !>   - sys : Parameter structure governing filenames.
  !>   - st  : Target state to be overwritten with disk data.
  !>
	SUBROUTINE initialise_from_file(sys,st)
		type(param), intent(in)			  	::  sys
		type(state), intent(in out)	  		::  st
		integer  									::  i,j,m,k
		character(len=200)						::  fks_file,ps_file
		real(rl)									::  f_r,f_i,h_r,h_i,p_r,q_r,p_i,q_i,a

		if (sys%verbose==1) then
			print*, "Initialising from: ", parameterchar(sys)
		end if
		fks_file="data/fks_fst_"//trim(adjustl(parameterchar(sys)))//".d"
		ps_file="data/ps_fst_"//trim(adjustl(parameterchar(sys)))//".d"
		open (unit=101,file=ps_file,action="read",status="old")
		open (unit=100,file=fks_file,action="read",status="old")

		do  k=1,sys%nmode
			read(100,'(f25.15)',advance='no') a
			do i=1,st%np
				read(100,'(2f25.15)',advance='no') f_r, h_r
				st%f(i,k) = f_r
				st%h(i,k) = h_r
			end do
			do i=1,st%np
				read(100,'(2f25.15)',advance='no') f_i, h_i
				st%f(i,k) = st%f(i,k) + Ic*f_i
				st%h(i,k) = st%h(i,k) + Ic*h_i
			end do
			read(100,*)
		end do

		do i=1,st%np
			read(101,'(2f25.15)',advance='no') p_r, q_r
			st%p(i) = p_r
			st%q(i) = q_r
		end do
		do i=1,st%np
			read(101,'(2f25.15)',advance='no') p_i, q_i
			st%p(i) = st%p(i) + Ic*p_i
			st%q(i) = st%q(i) + Ic*q_i
		end do

		st%t = sys%tmax
		close(100)
		close(101)

		!-- updating the sums over k
		CALL update_sums(sys,st)
		CALL normalise(st)
	END SUBROUTINE

  !> -------------------------------------------------------------------------
  !> FUNCTION: parameterchar
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Standardized file-tagging utility. Concatenates the critical physical 
  !>   and numerical parameters (α, Δ, error thresholds, grid cuts) into a 
  !>   unique string. Prevents collision and overwriting during cluster parameter 
  !>   sweeps.
  !> Arguments:
  !>   - sys : Parameter structure.
  !> Return:
  !>   - character(len=100) : The formatted signature string.
  !>
	FUNCTION parameterchar(sys)
		type(param), intent(in)		::   sys
		character(len=100)      		:: delchar,alChar,npiniChar,nmChar,npaddChar,&
			tmaxchar,dtchar,p1ichar, merrchar,trefchar, &
			p0char,wigxminchar,wigxmaxchar, &
			wmaxchar,wcchar,prepchar
		character(len=200)				:: parameterChar, addchar

		write(delchar, '(f6.1)') sys%del
		write(alChar, '(f8.4)') sys%alpha
		write(p1iChar, '(f6.2)') sys%p1i
		write(npinichar, '(i2)') sys%npini
		write(npaddchar, '(i2)') sys%npadd
		write(p0char, '(f10.8)') sys%p0
		write(merrchar, '(f11.8)') sys%merr
		write(nmChar, '(I5)') sys%nmode
		write(tmaxchar, '(I10)') int(sys%tmax)
		write(trefchar, '(f7.2)') sys%tref
		write(wigxminchar, '(I5)') int(sys%wigxmin)
		write(wigxmaxchar, '(I5)') int(sys%wigxmax)
		write(dtchar, '(f6.4)') sys%dt
		write(wmaxchar, '(I4)') int(sys%wmax)
		write(wcchar, '(I4)') int(sys%wc)
		write(prepChar, '(I2)') sys%prep

		addchar="_"
		if (sys%npadd .ne. 0) then
			addchar="_tr"//trim(adjustl(trefchar))//&
				!"_"//trim(adjustl(tac1char))//&
			"_me"//trim(adjustl(merrchar))//&
				"_p"//trim(adjustl(p0char))//"_"
		end if

		parameterchar=trim(adjustl(nmchar))//"m"//&
			"np"//trim(adjustl(npinichar))//"_"//trim(adjustl(npaddchar))//&
			"_al"//trim(adjustl(alchar))//&
			"_del"//trim(adjustl(delchar))//&
			"_dt"//trim(adjustl(dtchar))//&
			trim(adjustl(addchar))//&
			"tmax"//trim(adjustl(tmaxchar))//&
			"_wc"//trim(adjustl(wcchar))//&
			"_wmax"//trim(adjustl(wmaxchar))//&
			"_p"//trim(adjustl(prepchar))
	END FUNCTION

  !> -------------------------------------------------------------------------
  !> FUNCTION: gs_filename
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Returns the expected filename for the pre-calculated dressed ground state. 
  !>   Because the USC vacuum contains virtual photons, several scattering 
  !>   preparations require loading this state to properly isolate the scattered 
  !>   wavepacket from the bound dressing cloud.
  !> Arguments:
  !>   - sys     : Parameter structure.
  !>   - st      : Current state (determines the required polaron count).
  !>   - path_in : Optional base directory override.
  !>
	FUNCTION gs_filename(sys,st,path_in)
		type(param), intent(in)   		::  sys
		type(state), intent(in)			::  st
		character, intent(in),optional 	::  path_in
		character(len=200)				 	::  gs_path
		character(len=10)					::  alchar,delchar,npchar,nmchar,wmchar,wcchar
		character(len=200)					::  gs_filename

		gs_path = "gs_data"
		if (present(path_in)) then
			gs_path = path_in
		end if

		write(alchar,'(f6.3)') sys%alpha
		write(delchar,'(f5.2)') sys%del
		write(npchar,'(I2)') st%np
		write(nmchar,'(I7)') sys%nmode
		write(wcchar,'(I4)') int(sys%wc+0.0001)		
		write(wmchar,'(I4)') int(sys%wmax+0.0001)

		gs_filename= trim(adjustl(gs_path))//"/FinalFksP_al"//trim(adjustl(alchar))//&
			"_del"//trim(adjustl(delchar))//&
			"_pols"//trim(adjustl(npchar))//"_nk"//trim(adjustl(nmchar))//&
			"_wm"//trim(adjustl(wmchar))//"_wc"//trim(adjustl(wcchar))//".d"
	END FUNCTION

  !> -------------------------------------------------------------------------
  !> SUBROUTINE: CalcDerivatives
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   THE CORE PHYSICS ENGINE. This routine solves the highly non-linear, 
  !>   implicit Euler-Lagrange equations derived from the Time-Dependent 
  !>   Variational Principle (TDVP).
  !>   
  !>   Methodology (Ref: Appendix of 2017 Paper):
  !>   The exact functional derivative yields implicit equations where $\dot{f}$ 
  !>   appears on both the LHS and within the $\kappa$ overlap matrices on the RHS. 
  !>   To resolve this, the routine algebraically maps the tensor overlaps into a 
  !>   linear system $\mathbf{A}_{2D} \cdot \vec{X} = \vec{B}$. 
  !>   By flattening the indices and solving this system via LAPACK `directInverse`, 
  !>   it extracts the explicit derivatives $\dot{p}, \dot{q}, \dot{f}, \dot{h}$. 
  !>   This mapping reduces the computational complexity from an intractable 
  !>   $O(N_{modes} \times N_{cs}^3)$ to $O(N_{cs}^6)$, enabling continuum simulation.
  !> Arguments:
  !>   - sys : Parameter structure.
  !>   - st  : Current state (in), returns updated time derivatives (out).
  !>
	SUBROUTINE CalcDerivatives(sys,st)
		type(param), intent(in)                         :: sys
		type(state), intent(in out)                     :: st
		logical                                         :: superInverseFlag
		complex(cx), dimension(st%np)                   :: bigP,bigQ
		complex(cx), dimension(st%np, sys%nmode)        :: bigF, bigH
		complex(cx), dimension(st%np,st%np)             :: inv_ov_ff,inv_ov_hh
		complex(cx), dimension(st%np,st%np)             :: a_f, a_h, b_f, b_h
		complex(cx), dimension(st%np, sys%nmode)        :: rhsA_f, rhsA_h
		complex(cx), dimension(st%np**2)                :: packed_dRHS_f, packed_dRHS_h, packedSol_f, packedSol_h, solsave_f, solsave_h
		complex(cx), dimension(st%np,st%np)             :: tempMatrix_f, tempMatrix_h, dRHS_f, dRHS_h
		complex(cx), dimension(st%np)                   :: tempMatrixTer_f, tempMatrixTer_h
		complex(cx), dimension(st%np, sys%nmode)        :: tempfDot, temphDot
		complex(cx), dimension(st%np,st%np,st%np)       :: alphaT_f,alphaT_h
		integer                                         :: info,i,j,k,n,m,l

		integer,save                                    :: counter
		real(rl)                                        :: startTime, endTime


		!-- initialisations
		!print*, 'initialisations'
		info=0
		do i=1,st%np
			bigP(i) = P_j(sys,st,i)
			bigQ(i) = Q_j(sys,st,i)
			bigF(i,:) =  F_j(sys,st,i)
			bigH(i,:) =  H_j(sys,st,i)
		end do

		!-- invert overlap matrices
		!print*, 'invert overlap matrix'
		inv_ov_ff=st%ov_ff
		inv_ov_hh=st%ov_hh
		CALL invertH(inv_ov_ff,info)
		CALL invertH(inv_ov_hh,info)

		!-- build b matrices
		!print*, 'build b matrices'
		b_f=matmul(CONJG(st%f),TRANSPOSE(st%f))
		b_h=matmul(CONJG(st%h),TRANSPOSE(st%h))
		b_f=matmul(CONJG(st%f),TRANSPOSE(st%f))
		b_h=matmul(CONJG(st%h),TRANSPOSE(st%h))

		!-- build RHS
		!print*, 'build rhs'
		do k=1, sys%nmode
			do n=1, st%np
				rhsA_f(n,k)=sum(inv_ov_ff(n,:)*(bigF(:,k)-st%f(n,k)*bigP(:)))
				rhsA_h(n,k)=sum(inv_ov_hh(n,:)*(bigH(:,k)-st%h(n,k)*bigQ(:)))
			end do
		end do

		dRHS_f=matmul(CONJG(st%f),TRANSPOSE(rhsA_f))
		dRHS_f=st%ov_ff*dRHS_f
		dRHS_f=matmul(inv_ov_ff,dRHS_f)
		dRHS_h=matmul(CONJG(st%h),TRANSPOSE(rhsA_h))
		dRHS_h=st%ov_hh*dRHS_h
		dRHS_h=matmul(inv_ov_hh,dRHS_h)

		do i=1, st%np
			do n=1, st%np
				packed_dRHS_f((n-1)*st%np+i)=dRHS_f(i,n)
				packed_dRHS_h((n-1)*st%np+i)=dRHS_h(i,n)
			end do
		end do


		!-- build alphaTensor
		do i=1, st%np
			do n=1, st%np
				do m=1, st%np
					alphaT_f(i,n,m)=sum(inv_ov_ff(i,:)*st%ov_ff(:,n)*(b_f(:,m)-b_f(:,n)) )
					alphaT_h(i,n,m)=sum(inv_ov_hh(i,:)*st%ov_hh(:,n)*(b_h(:,m)-b_h(:,n)) )
				end do
			end do
		end do

		! -- SuperInverse detection
		CALL directInverse(st%np, alphaT_f, packedSol_f, packed_dRHS_f)
		CALL directInverse(st%np, alphaT_h, packedSol_h, packed_dRHS_h)

		!-- system unpack
		! print*, 'unpack'
		do i=1, st%np
			do n=1, st%np
				a_f(n,i)=packedSol_f((n-1)*st%np+i)
				a_h(n,i)=packedSol_h((n-1)*st%np+i)
			end do
		end do

		a_f=matmul(st%ov_ff, a_f)
		a_f=TRANSPOSE(a_f/st%ov_ff)  !! -- This TRANSPOSE is a mystery.

		a_h=matmul(st%ov_hh, a_h)
		a_h=TRANSPOSE(a_h/st%ov_hh)  !! -- This TRANSPOSE is a mystery.

		!-- fDot and hdot extraction
		! print*, 'fDot and hdot'

		tempMatrix_f=matmul(inv_ov_ff,st%ov_ff*TRANSPOSE(a_f))
		do i=1, st%np
			tempfDot(i,:)=SUM(tempMatrix_f(i,:))*st%f(i,:)
		end do
		st%fDot= rhsA_f-matmul(tempMatrix_f,st%f)+tempfDot
		do i=1, st%np
			st%fDot(i,:)=st%fDot(i,:)/st%p(i)
		end do


		tempMatrix_h=matmul(inv_ov_hh,st%ov_hh*TRANSPOSE(a_h))
		do i=1, st%np
			temphDot(i,:)=SUM(tempMatrix_h(i,:))*st%h(i,:)
		end do
		st%hDot= rhsA_h-matmul(tempMatrix_h,st%h)+temphDot
		do i=1, st%np
			st%hDot(i,:)=st%hDot(i,:)/st%q(i)
		end do

		!-- evaluate pDot and qDot
		! print*, 'pdot and qdot'
		tempMatrixTer_f= MATMUL(inv_ov_ff,bigP)
		st%pDot= 0.5_rl*( (/ (a_f(n,n), n=1, st%np) /) + st%p*(/ (conjg(a_f(m,m)), m=1, st%np) /) /CONJG(st%p) )
		st%pdot=st%pDot + tempMatrixTer_f
		st%pDot=st%pDot - SUM(tempMatrix_f, dim=2)

		tempMatrixTer_h= MATMUL(inv_ov_hh,bigQ)
		st%qDot= 0.5_rl*( (/ (a_h(m,m), m=1, st%np) /) + st%q*(/ (conjg(a_h(m,m)), m=1, st%np) /) /CONJG(st%q) )
		st%qdot=st%qDot + tempMatrixTer_h
		st%qDot=st%qDot - SUM(tempMatrix_h, dim=2)
	END SUBROUTINE


  !> -------------------------------------------------------------------------
  !> FUNCTION: P_j, Q_j, F_j, H_j
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   These four functions evaluate the raw local gradients of the total system 
  !>   energy with respect to the complex conjugate of the variational parameters 
  !>   ($\partial E / \partial p^*_j$, etc.). They form the foundational Right-Hand 
  !>   Side driving terms for the EOMs prior to the $\kappa$-tensor inversion.
  !> Arguments:
  !>   - sys : Parameter structure.
  !>   - st  : Current state.
  !>   - m/j : The specific polaron index being evaluated.
  !> Return:
  !>   - complex(cx) : Evaluated derivative vector/scalar.
  !>
	FUNCTION P_j(sys,st,m)

		type(param),intent(in)         ::  sys
		type(state),intent(in)         ::  st
		complex(cx) 					     ::  P_j
		integer,intent(in)             ::  m     !-- index of the p
		integer 							  ::  j

		P_j = 0._rl
		do j=1, size(st%p,1)
			P_j = P_j  &
				+ 0.5_cx * sys%del * st%q(j) * st%ov_fh(m,j) & 
				+ st%p(j)*st%ov_ff(m,j) * ( st%bigW_f(m,j) - 0.5_cx * st%bigL_f(m,j) )
		end do
		P_j = -Ic*P_j

	END FUNCTION 
	FUNCTION Q_j(sys,st,m)

		type(param),intent(in)  		   ::  sys
		type(state),intent(in)          ::  st
		complex(cx) 				         ::  Q_j
		integer,intent(in)              ::  m     !-- index of the q
		integer 								::  j

		Q_j = 0._rl
		do j=1, size(st%q,1)
			Q_j = Q_j &
				+ 0.5_rl*sys%del * st%p(j) * st%ov_hf(m,j) & 
				+ st%q(j)*st%ov_hh(m,j) * ( st%bigW_h(m,j) + 0.5_rl * st%bigL_h(m,j) )
		end do
		Q_j = -Ic*Q_j

	END FUNCTION 
	FUNCTION F_j(sys,st,j)

		type(param),intent(in)  		   			 ::  sys
		type(state),intent(in)          			 ::  st
		complex(cx) 				         			 ::  F_j(sys%nmode)
		integer,intent(in)              			 ::  j     !-- index of the f
		integer 											 ::  m,s

		F_j  = 0._cx
		do m=1, size(st%f,1)
			F_j(:) = F_j(:) &
				+ 0.5_rl*sys%del*st%q(m)*st%h(m,:)*st%ov_fh(j,m) &
				+ st%p(m)*( st%f(m,:)*st%bigW_f(j,m) + sys%w(:)*st%f(m,:) )*st%ov_ff(j,m) &
				- 0.5_rl*st%p(m)*st%ov_ff(j,m)*( st%f(m,:)*st%bigL_f(j,m) + sys%g(:) )
		end do
		F_j(:) = - Ic*F_j(:) 

	END FUNCTION
	FUNCTION H_j(sys,st,j)

		type(param),intent(in)  		   			 ::  sys
		type(state),intent(in)          			 ::  st
		complex(cx) 				         			 ::  H_j(sys%nmode)
		integer,intent(in)              			 ::  j     !-- index of the f
		integer 											 ::  m,s

		H_j  = 0._cx
		do m=1, size(st%f,1)
			H_j(:) = H_j(:) &
				+ 0.5_rl*sys%del*st%p(m)*st%f(m,:)*st%ov_hf(j,m) &
				+ st%q(m)*( st%h(m,:)*st%bigW_h(j,m) + sys%w(:)*st%h(m,:) )*st%ov_hh(j,m) &
				+ 0.5_rl*st%q(m)*st%ov_hh(j,m)*( st%h(m,:)*st%bigL_h(j,m) + sys%g(:) )
		end do
		H_j(:) = - Ic*H_j(:) 

	END FUNCTION

  !> -------------------------------------------------------------------------
  !> SUBROUTINE: update_sums
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Critical performance optimizer. The coherent states are non-orthogonal, 
  !>   meaning their inner products $\langle f_i | f_j \rangle$ appear everywhere 
  !>   in the energy landscape. By pre-computing these macroscopic overlaps and 
  !>   the free/interaction energy projections (the `bigW` and `bigL` matrices) 
  !>   once per RK4 step, it prevents the $O(N_{cs}^3)$ derivative loops from 
  !>   having to iterate over the massive momentum grid ($N_{modes}$).
  !> Arguments:
  !>   - sys : Parameter structure.
  !>   - st  : Target state (caches updated in-place).
  !>
	SUBROUTINE update_sums(sys,st)

		type(param),intent(in)			::  sys
		type(state),intent(in out)   ::  st
		integer								::  i,j,s

		do i=1, size(st%f,1)
			do j=1, size(st%f,1)

				st%ov_ff(i,j) = ov(st%f(i,:),st%f(j,:))
				st%ov_hh(i,j) = ov(st%h(i,:),st%h(j,:))
				st%ov_fh(i,j) = ov(st%f(i,:),st%h(j,:))
				st%ov_hf(i,j) = ov(st%h(i,:),st%f(j,:))

				if (sys%nmode .ne. 1) then
					st%bigW_f(i,j) = dot_product( sys%w , conjg(st%f(i,:)) * st%f(j,:) )
					st%bigW_h(i,j) = dot_product( sys%w , conjg(st%h(i,:)) * st%h(j,:) )
					st%bigL_f(i,j) = dot_product( sys%g , conjg(st%f(i,:)) + st%f(j,:) )
					st%bigL_h(i,j) = dot_product( sys%g , conjg(st%h(i,:)) + st%h(j,:) )
				else if (sys%nmode == 1) then
					st%bigW_f(i,j) = sys%w(1)*conjg(st%f(i,1)) * st%f(j,1)
					st%bigW_h(i,j) = sys%w(1)*conjg(st%h(i,1)) * st%h(j,1)
					st%bigL_f(i,j) = sys%g(1)*( conjg(st%f(i,1)) + st%f(j,1) )
					st%bigL_h(i,j) = sys%g(1)*( conjg(st%h(i,1)) + st%h(j,1) )
				end if

			end do
		end do

	END SUBROUTINE

	!=======================================
	!== State functions
	!======================================

  !> -------------------------------------------------------------------------
  !> FUNCTION: Energy
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Evaluates the total variational energy functional $E[\Psi] = \langle \Psi | H | \Psi \rangle$. 
  !>   This is the central quantity minimized by the Dirac-Frenkel
  !>   Time-Dependent Variational Principle (TDVP). It sums the bare atomic energy, 
  !>   the free-field bosonic energy (via `bigW`), and the spin-boson dipole 
  !>   interaction energy (via `bigL`). Because the coherent states are non-orthogonal, 
  !>   every term is weighted by the complex cross-overlaps.
  !> Arguments:
  !>   - sys : Parameter structure (provides $\Delta$).
  !>   - st  : Current multi-polaron state object.
  !> Return:
  !>   - real(rl) : The exact expectation value of the Hamiltonian. Throws a 
  !>                warning if numerical instability produces an imaginary component.
  !>
	FUNCTION Energy(sys,st)
		type(param),intent(in)         				::  sys
		type(state),intent(in)         				::  st
		complex(cx)						  				::  tmp
		real(rl)							  				::  energy
		integer 							  				::  i,j

		tmp=0._cx
		energy=0._rl

		do i=1,size(st%p,1)
			do j=1,size(st%p,1)

				tmp = tmp &
					+ 0.5_rl * sys%del*( conjg(st%p(i))*st%q(j)*st%ov_fh(i,j) + conjg(st%q(j))*st%p(i)*st%ov_hf(j,i) ) &
					+ ( conjg(st%p(i))*st%p(j)*st%ov_ff(i,j) * st%bigW_f(i,j) + conjg(st%q(i))*st%q(j)*st%ov_hh(i,j) * st%bigW_h(i,j) ) &
					- 0.5_rl * ( conjg(st%p(i))*st%p(j)*st%ov_ff(i,j) * st%bigL_f(i,j) - conjg(st%q(i))*st%q(j)*st%ov_hh(i,j) * st%bigL_h(i,j) )

			end do
		end do

		if (aimag(tmp) > 1e-10_rl) then
			print*, "Error: Energy is complex"
		end if

		energy = real(tmp,KIND=8)

	END FUNCTION

  !> -------------------------------------------------------------------------
  !> FUNCTION: norm
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes the overall $L^2$ norm of the multi-polaron state, $\langle \Psi | \Psi \rangle$, 
  !>   to verify total probability conservation. It evaluates the double sum over 
  !>   all coherent state cross-terms, factoring in both the atomic spin amplitudes 
  !>   (`p`, `q`) and the bosonic inner products (`ov_ff`, `ov_hh`).
  !> Arguments:
  !>   - st : Current state object to evaluate.
  !> Return:
  !>   - real(rl) : The scalar norm of the variational wavefunction.
  !>
	FUNCTION norm(st)

		real(rl) 			       		::  norm
		complex(cx)						::  tmp
		type(state),intent(in)  		::  st
		integer						 		::  m,n

		tmp = 0._cx
		do m=1,size(st%p,1)
			do n=1,size(st%p,1)
				tmp = tmp  &
					+ conjg(st%p(m))*st%p(n)*st%ov_ff(m,n)  &
					+ conjg(st%q(m))*st%q(n)*st%ov_hh(m,n)
			end do
		end do

		norm = sqrt(tmp)

	END FUNCTION

  !> -------------------------------------------------------------------------
  !> SUBROUTINE: normalise
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Enforces strict probability conservation. Evaluates the instantaneous 
  !>   state norm and scales the complex atomic probability amplitudes (`p` and `q`) 
  !>   to restore $\langle \Psi | \Psi \rangle = 1$. The bosonic displacement arrays 
  !>   (`f`, `h`) are left untouched, as scaling them would unphysically alter 
  !>   the photon number rather than the overall state weight.
  !> Arguments:
  !>   - st : Target state to be normalized in-place.
  !>
	SUBROUTINE normalise(st)

		type(state), intent(in out)  ::  st
		real(rl)							:: normval

		normval=norm(st)

		st%p = st%p/normval
		st%q = st%q/normval

	END SUBROUTINE

  !> -------------------------------------------------------------------------
  !> FUNCTION: error
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   The adaptive trigger for basis expansion. Calculates the McLachlan 
  !>   variational error metric, defined as the squared distance between the 
  !>   exact Schrödinger evolution and the restricted TDVP trajectory: 
  !>   $Err = || (i \partial_t - H)|\Psi(t)\rangle ||^2$. 
  !>   
  !>   Physics Context:
  !>   In the ultrastrong coupling regime, inelastic scattering produces highly 
  !>   entangled, non-Gaussian multi-photon states (Schrödinger cats). When the 
  !>   current basis dimension ($N_{cs}$) can no longer capture this complexity, 
  !>   this error diverges. The main loop monitors this function and injects 
  !>   a new polaron when the threshold is breached.
  !> Arguments:
  !>   - sys  : Parameter structure.
  !>   - rost : Unused/legacy variable (historically $t-2dt$).
  !>   - st   : Current state object at time $t$.
  !> Return:
  !>   - complex(cx) : The complex error magnitude (imaginary part should be 0).
  !>
	FUNCTION error(sys,rost,st)

		type(param),intent(in)		::  sys
		type(state),intent(in) 		::  rost,st
		type(state) 				::  ost
		complex(cx)					::  tmp1,tmp2, tmp3
		complex(cx)	   				::  error
		complex(cx), dimension(st%np,sys%nmode)  ::  off, ohh, ofh, ohf,ofdd,ohdd
		integer					   	::  i,j

		tmp1 = 0._cx
		tmp2 = 0._cx
		tmp3 = 0._cx
		error = 0._cx

		ost = st

		do i=1, size(st%f,1)
			do j=1, size(ost%f,1)
				off(i,j) = dot_product( ost%fdot(i,:) , ost%f(j,:) )
				ohh(i,j) = dot_product( ost%hdot(i,:) , ost%h(j,:) )
				ofh(i,j) = dot_product( ost%fdot(i,:) , ost%h(j,:) )
				ohf(i,j) = dot_product( ost%hdot(i,:) , ost%f(j,:) )
			end do
		end do

		tmp3 = tmp3 + 0.25_rl*sys%del**2 + 0.25_rl*sum(sys%g(:)*sys%g(:))

		do i=1,ost%np
			do j=1,ost%np

				!- - dot p  * dot p
				tmp1 = tmp1 + ost%ov_ff(i,j)*( &
					+ conjg(ost%pdot(i))*ost%pdot(j) - 0.5_rl*conjg(ost%pdot(i))*ost%p(j)*( conjg(off(j,j)) + off(j,j) - 2_rl*conjg(off(j,i)) )&
					- 0.5_rl*conjg(ost%p(i))*ost%pdot(j)*( conjg(off(i,i)) + off(i,i) - 2_rl*off(i,j) ) &
					+ 0.25_rl*conjg(ost%p(i))*ost%p(j)*( (conjg(off(i,i))+off(i,i))*( conjg(off(j,j)) + off(j,j) - 2_rl*conjg(off(j,i))) &
					- 2_rl*off(i,j)*( conjg(off(j,j)) + off(j,j) ) &
					+ 4_rl*( off(i,j)*conjg(off(j,i)) + sum(conjg(ost%fdot(i,:))*ost%fdot(j,:)) ) ) )

				tmp2 = tmp2 + ost%ov_ff(i,j)*ost%p(j)* ( &
					( conjg(ost%pdot(i)) - 0.5_rl*conjg(ost%p(i))*(conjg(off(i,i)) + off(i,i)) )*( ost%bigW_f(i,j) &
					- 0.5_rl*ost%bigL_f(i,j) ) &
					+ conjg(ost%p(i))*off(i,j)*( ost%bigW_f(i,j) - 0.5_rl*ost%bigL_f(i,j) ) &
					+ conjg(ost%p(i))*sum( sys%w(:)*conjg(ost%fdot(i,:))*ost%f(j,:) - 0.5_rl*conjg(ost%fdot(i,:))*sys%g(:)) ) &
					+ 0.5_rl*sys%del*ost%q(j)*ost%ov_fh(i,j)*( conjg(ost%pdot(i)) & 
					- 0.5_rl*conjg(ost%p(i))*( conjg(off(i,i)) + off(i,i) - 2_rl*ofh(i,j) ) )


				tmp3 = tmp3 + conjg(ost%p(i))*ost%p(j)*ost%ov_ff(i,j)*( sum( sys%w(:)*sys%w(:)*conjg(ost%f(i,:))*ost%f(j,:) ) &
					+ (ost%bigW_f(i,j))**2 &
					+ 0.25_rl*(ost%bigL_f(i,j))**2 &
					- 0.5_rl*sum( sys%g(:)*sys%w(:)*(conjg(ost%f(i,:)) + ost%f(j,:)) ) &
					- ost%bigL_f(i,j)*ost%bigW_f(i,j) ) &
					+ sys%del*conjg(ost%p(i))*ost%q(j)*ost%ov_fh(i,j)*sum( sys%w(:)*conjg(ost%f(i,:))*ost%h(j,:) )

				!tmp4 = tmp4 + conjg(ost%p(j))*ost%ov_ff(j,i)*( &
				!	+ opdd(i) &
				!	- ost%pdot(i)*( off(i,i) + conjg(off(i,i)) - 2_rl*conjg(off(i,j)) ) &
				!	- 0.5_rl*ost%p(i)*( sum( 2_rl*real(conjg(ofdd(i,:))*ost%f(i,:)) &
				!	+ 2_rl*conjg(ost%fdot(i,:))*ost%fdot(i,:) &
				!	- 2_rl*ofdd(i,:)*conjg(ost%f(j,:)) ) ) &
				!	+ 0.25_rl*ost%p(i)*( off(i,i) + conjg(off(i,i)) - 2_rl*conjg(off(i,j)) )**2 &
				!	)


				!-- the q and h parts 

				tmp1 = tmp1 + ost%ov_hh(i,j)*( &
					+ conjg(ost%qdot(i))*ost%qdot(j) - 0.5_rl*conjg(ost%qdot(i))*ost%q(j)*( conjg(ohh(j,j)) + ohh(j,j) - 2_rl*conjg(ohh(j,i)) )&
					- 0.5_rl*conjg(ost%q(i))*ost%qdot(j)*( conjg(ohh(i,i)) + ohh(i,i) - 2_rl*ohh(i,j) ) &
					+ 0.25_rl*conjg(ost%q(i))*ost%q(j)*( (conjg(ohh(i,i))+ohh(i,i))*( conjg(ohh(j,j)) + ohh(j,j) - 2_rl*conjg(ohh(j,i))) &
					- 2_rl*ohh(i,j)*( conjg(ohh(j,j)) + ohh(j,j) ) &
					+ 4_rl*( ohh(i,j)*conjg(ohh(j,i)) + sum(conjg(ost%hdot(i,:))*ost%hdot(j,:)) ) ) )


				tmp2 = tmp2 + ost%ov_hh(i,j)*ost%q(j)* ( &
					( conjg(ost%qdot(i)) - 0.5_rl*conjg(ost%q(i))*(conjg(ohh(i,i)) + ohh(i,i)) )*( ost%bigW_h(i,j) &
					- 0.5_rl*(-ost%bigL_h(i,j)) ) &
					+ conjg(ost%q(i))*ohh(i,j)*( ost%bigW_h(i,j) - 0.5_rl*(-ost%bigL_h(i,j)) ) &
					+ conjg(ost%q(i))*sum( sys%w(:)*conjg(ost%hdot(i,:))*ost%h(j,:) - 0.5_rl*conjg(ost%hdot(i,:))*(-sys%g(:)) ) ) &
					+ 0.5_rl*sys%del*ost%p(j)*ost%ov_hf(i,j)*( conjg(ost%qdot(i)) & 
					- 0.5_rl*conjg(ost%q(i))*( conjg(ohh(i,i)) + ohh(i,i) - 2_rl*ohf(i,j) ) )

				tmp3 = tmp3 + conjg(ost%q(i))*ost%q(j)*ost%ov_hh(i,j)*( sum( sys%w(:)*sys%w(:)*conjg(ost%h(i,:))*ost%h(j,:) ) &
					+ (ost%bigW_h(i,j))**2 &
					+ 0.25_rl*(-ost%bigL_h(i,j))**2 &
					- 0.5_rl*sum( (-sys%g(:))*sys%w(:)*(conjg(ost%h(i,:)) + ost%h(j,:)) ) &
					- (-ost%bigL_h(i,j))*ost%bigW_h(i,j) ) &
					+ sys%del*conjg(ost%q(i))*ost%p(j)*ost%ov_hf(i,j)*sum( sys%w(:)*conjg(ost%h(i,:))*ost%f(j,:) )

				!tmp4 = tmp4 + conjg(ost%q(j))*ost%ov_hh(j,i)*( &
				!	+ oqdd(i) &
				!	- ost%qdot(i)*( ohh(i,i) + conjg(ohh(i,i)) - 2_rl*conjg(ohh(i,j)) ) &
				!	- 0.5_rl*ost%q(i)*( sum( 2_rl*real(conjg(ohdd(i,:))*ost%h(i,:)) &
				!	+ 2_rl*conjg(ost%hdot(i,:))*ost%hdot(i,:) &
				!	- 2_rl*ohdd(i,:)*conjg(ost%h(j,:)) ) ) &
				!	+ 0.25_rl*ost%q(i)*( ohh(i,i) + conjg(ohh(i,i)) - 2_rl*conjg(ohh(i,j)) )**2 &
				!	)

			end do
		end do

		error =  tmp1 - 2._rl*aimag(tmp2) + tmp3
		!error =  -0.5_rl*real(tmp4) + 0.5_rl*tmp1 - 2._rl*aimag(tmp2) + tmp3

	END FUNCTION	

  !> -------------------------------------------------------------------------
  !> FUNCTION: ov_states
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes the full macroscopic quantum overlap between two entirely 
  !>   separate multi-polaron states, $\langle \Psi_1 | \Psi_2 \rangle$. Useful 
  !>   for calculating fidelities, transition probabilities, or verifying 
  !>   orthogonality between scattered wavepackets and the bound ground state.
  !> Arguments:
  !>   - st1 : First (bra) state object.
  !>   - st2 : Second (ket) state object.
  !> Return:
  !>   - complex(cx) : The complex scalar inner product.
  !>
	FUNCTION ov_states(st1,st2)

		type(state), intent(in)  ::  st1,st2
		complex(cx) 				  ::  ov_states 					
		complex(cx)				  ::  tmp
		integer						  ::  i,j

		tmp = 0._cx

		do i=1,st1%np
			do j=1,st1%np
				tmp = tmp + conjg(st1%p(i))*st2%p(j)*ov( st1%f(i,:),st2%f(j,:) ) &
					+ conjg(st1%q(i))*st2%q(j)*ov( st1%h(i,:),st2%h(j,:) )
			end do
		end do

		ov_states = tmp

	END FUNCTION

  !> -------------------------------------------------------------------------
  !> FUNCTION: ov_scalar
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Evaluates the fundamental inner product between two scalar (single-mode) 
  !>   coherent states $|\alpha\rangle$ and $|\beta\rangle$. 
  !>   Formula: $\exp( -0.5|\alpha|^2 - 0.5|\beta|^2 + \alpha^* \beta )$.
  !> Arguments:
  !>   - f1 : First scalar complex amplitude ($\alpha$).
  !>   - f2 : Second scalar complex amplitude ($\beta$).
  !> Return:
  !>   - complex(cx) : The non-orthogonal overlap scalar.
  !>
	FUNCTION ov_scalar(f1,f2)

		complex(cx), intent(in) :: f1, f2
		complex(cx)             :: ov_scalar
		complex(cx)				  :: tmp1, tmp2, tmp3

		tmp1 = conjg(f1)*f1
		tmp2 = conjg(f2)*f2
		tmp3 = conjg(f1)*f2

		ov_scalar = exp( -0.5_rl*tmp1 - 0.5_rl*tmp2 + tmp3 ) 

	END FUNCTION	ov_scalar

  !> -------------------------------------------------------------------------
  !> FUNCTION: ov
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Generalizes the coherent state overlap to the continuous multi-mode 
  !>   waveguide. Computes the macroscopic overlap between two displacement 
  !>   fields by taking the dot product across the entire momentum grid $k$. 
  !>   This is the core mathematical object that dictates the metric tensor 
  !>   of the variational manifold.
  !> Arguments:
  !>   - f1 : First continuous momentum array (bra).
  !>   - f2 : Second continuous momentum array (ket).
  !> Return:
  !>   - complex(cx) : The total macroscopic overlap.
  !>
	FUNCTION ov(f1,f2)

		complex(cx), intent(in) :: f1( : ), f2( : )
		complex(cx)             :: ov
		!= internal variables
		complex(cx)    :: tmp1, tmp2, tmp3

		!= initialize
		tmp1 = 0._rl
		tmp2 = 0._rl
		tmp3 = 0._rl

		if (size(f1,1) .ne. 1) then
			tmp1 = dot_product(f1, f1)
			tmp2 = dot_product(f2, f2)
			tmp3 = dot_product(f1, f2)
		else if (size(f1,1) == 1) then
			tmp1 = conjg(f1(1))*f1(1)
			tmp2 = conjg(f2(1))*f2(1)
			tmp3 = conjg(f1(1))*f2(1)
		end if

		ov = exp( -0.5_rl*tmp1 - 0.5_rl*tmp2 + tmp3 ) 

	END FUNCTION	ov

  !> -------------------------------------------------------------------------
  !> FUNCTION: upProb
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Calculates the reduced density matrix probability of finding the artificial 
  !>   atom in the excited $|\uparrow\rangle$ state. It traces out the bosonic 
  !>   environment by summing the squared atomic amplitudes (`p`), weighted by 
  !>   their associated multi-mode bosonic overlaps (`ov_ff`).
  !> Arguments:
  !>   - st : Current state object to evaluate.
  !> Return:
  !>   - real(rl) : The spin-up probability (0.0 to 1.0).
  !>
	FUNCTION upProb(st)

		real(rl) 			       		::  upProb
		type(state),intent(in)  		::  st
		integer						 		::  m,n

		upProb = 0_rl
		do m=1,size(st%p,1)
			do n=1,size(st%p,1)
				upProb =  upProb  &
					+ conjg(st%p(m))*st%p(n)*st%ov_ff(m,n)
			end do
		end do

	END FUNCTION

  !> -------------------------------------------------------------------------
  !> FUNCTION: downProb
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Calculates the reduced density matrix probability of finding the artificial 
  !>   atom in the ground $|\downarrow\rangle$ state. Traces out the bosonic 
  !>   environment by evaluating the `q` amplitudes and the `ov_hh` cross-overlaps.
  !> Arguments:
  !>   - st : Current state object to evaluate.
  !> Return:
  !>   - real(rl) : The spin-down probability (0.0 to 1.0).
  !>
	FUNCTION downProb(st)

		real(rl) 			       		::  downProb
		type(state),intent(in)  		::  st
		integer						 		::  m,n

		downProb = 0_rl
		do m=1,size(st%p,1)
			do n=1,size(st%p,1)
				downProb  =  downProb  &
					+ conjg(st%q(m))*st%q(n)*st%ov_hh(m,n)
			end do
		end do

	END FUNCTION

  !> -------------------------------------------------------------------------
  !> FUNCTION: sigmaX
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes the expectation value of the Pauli X operator, $\langle \sigma_x \rangle$. 
  !>   Because $\sigma_x$ flips the atomic spin, this evaluates the off-diagonal 
  !>   coherences between the "up" (`p`) and "down" (`q`) manifolds, weighted 
  !>   by the corresponding cross-branch bosonic overlaps (`ov_fh` and `ov_hf`). 
  !>   Crucial for tracking tunneling dynamics driven by $\Delta$.
  !> Arguments:
  !>   - st : Current state object to evaluate.
  !> Return:
  !>   - real(rl) : The expectation value of $\sigma_x$.
  !>
	FUNCTION sigmaX(st)

		type(state), intent(in)  :: st
		real(rl) 					  :: sigmaX
		complex(cx) 				  :: tmp
		integer						  :: n,m

		tmp = 0_cx
		sigmaX = 0_rl

		do n=1,size(st%p,1)
			do m=1,size(st%p,1)
				tmp = tmp + conjg(st%p(n))*st%q(m)*st%ov_fh(n,m) &
					+ conjg(st%q(n))*st%p(m)*st%ov_hf(n,m)
			end do	
		end do
		sigmaX = real(tmp)

	END FUNCTION

  !> -------------------------------------------------------------------------
  !> FUNCTION: sigmaZ
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes the expectation value of the Pauli Z operator, $\langle \sigma_z \rangle$, 
  !>   representing the atomic population inversion. Evaluated simply as the 
  !>   difference between the spin-up and spin-down reduced probabilities.
  !>   In the USC regime, the true ground state is dressed, meaning $\langle \sigma_z \rangle$ 
  !>   will not perfectly reach -1 even after full spontaneous emission.
  !> Arguments:
  !>   - st : Current state object to evaluate.
  !> Return:
  !>   - real(rl) : The expectation value of $\sigma_z$ (ranges from -1 to 1).
  !>
	FUNCTION sigmaZ(st)

		type(state), intent(in)  :: st
		real(rl) 					  :: sigmaZ
		integer						  :: n,m

		sigmaZ = upProb(st) - downProb(st)

	END FUNCTION

  !> -------------------------------------------------------------------------
  !> FUNCTION: sigmaY
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes the expectation value of the Pauli Y operator, $\langle \sigma_y \rangle$. 
  !>   Similar to $\sigma_x$, it evaluates the off-diagonal coherences between 
  !>   the atomic manifolds, but applies the appropriate complex $\pm i$ phase 
  !>   factor during the bosonic trace.
  !> Arguments:
  !>   - st : Current state object to evaluate.
  !> Return:
  !>   - real(rl) : The expectation value of $\sigma_y$.
  !>
	FUNCTION sigmaY(st)

		type(state), intent(in)  :: st
		real(rl) 					  :: sigmaY
		real(rl) 					  :: tmp
		integer						  :: n,m

		tmp = 0_cx
		sigmaY = 0_rl

		do n=1,size(st%p,1)
			do m=1,size(st%p,1)
				tmp = tmp + 2._rl*real( Ic * conjg(st%q(n))*st%p(m)*st%ov_hf(n,m) )
			end do	
		end do
		sigmaY = tmp

	END FUNCTION

  !> -------------------------------------------------------------------------
  !> FUNCTION: f_nk_FT
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   A specialized spatial windowing and Fourier Transform utility. 
  !>   To analyze the emitted Schrödinger cat state without interference from 
  !>   the virtual photons dressing the atom, this function applies a spatial 
  !>   cutoff (`xmin`, `xmax`) to isolate the freely propagating wavepacket. 
  !>   It then transforms the filtered real-space data back into the momentum 
  !>   domain to extract the true emission spectrum.
  !> Arguments:
  !>   - sys        : Parameter structure defining the grids.
  !>   - fnx        : The real-space amplitude array to be filtered.
  !>   - xmin, xmax : Optional spatial boundaries to isolate the propagating cat.
  !> Return:
  !>   - complex(cx) : The filtered, momentum-space amplitude array.
  !>
	FUNCTION f_nk_FT(sys,fnx,xmin,xmax)

		type(param),intent(in)   						::  sys
		complex(cx), intent(in)  						::  fnx(:,:)
		real(rl),intent(in),optional					::  xmin,xmax
		complex(cx), dimension( size(fnx,1),size(fnx,2) )  	::  f_nk_FT, new_fnx
		complex(cx), dimension(sys%nmode)  		::  x_arr,fnx_imin
		integer												::  n,i, imin, imax

		if (present(xmin)) then
			imin = int(xmin/sys%dx) + 1
			imax = int(xmax/sys%dx) + 1
		else
			imin = 1
			imax = sys%nmode
		end if

		do i=1, sys%nmode
			x_arr(i) = sys%dx * (i-0.5_rl)
		end do
		new_fnx = 0._cx
		new_fnx(:,imin:) = fnx(:,imin:) 

		do n=1,size(fnx,1)
			do i=1,sys%nmode
				f_nk_FT(n,i) = ( 2._rl/sqrt(2._rl*pi) )*sqrt(sys%dx*sys%dk1)*sum( cos(sys%w(i)*x_arr(:))*new_fnx(n,:) )
			end do
		end do

	END FUNCTION

  !> -------------------------------------------------------------------------
  !> FUNCTION: n_up_k
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Calculates the momentum-space photon number spectrum $n(k)$ strictly 
  !>   projected onto the artificial atom's excited $|\uparrow\rangle$ state. 
  !>   It collapses the "down" amplitudes to zero, normalizes the subspace, 
  !>   and evaluates the expectation value of the photon number operator.
  !> Arguments:
  !>   - sys : Parameter structure.
  !>   - st  : Current multi-polaron state.
  !>   - i   : The specific momentum grid index ($k_i$) to evaluate.
  !> Return:
  !>   - real(rl) : The spin-projected photon density at mode $i$.
  !>
	FUNCTION n_up_k(sys,st,i)

		type(state), intent(in)   						::  st
		type(param), intent(in)							::  sys
		integer, intent(in)									::  i
		real(rl)												::  n_up_k

		type(state)											::  upst
		complex(cx)											::  tmp,sum_k
		integer													::  n,m

		upst = st
		upst%q(:) = 0._cx
		CALL normalise(upst)

		tmp = 0._cx
		do n=1,st%np
			do m=1,st%np

				tmp = tmp + conjg(upst%p(n))*upst%p(m)*(conjg(upst%f(n,i))*upst%f(m,i))*upst%ov_ff(n,m)

			end do
		end do

		n_up_k = real( tmp )

	END FUNCTION

  !> -------------------------------------------------------------------------
  !> FUNCTION: n_up_x
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Calculates the spatially integrated, real-space photon density conditioned 
  !>   on the atom being in the $|\uparrow\rangle$ state. It accepts spatial 
  !>   windowing arguments to specifically quantify the macroscopic photon content 
  !>   of the propagating wavepacket branch associated with the excited atom.
  !> Arguments:
  !>   - sys        : Parameter structure.
  !>   - st         : Current state.
  !>   - xmin, xmax : Optional spatial bounds.
  !> Return:
  !>   - real(rl) : Total projected macroscopic photon count in the spatial window.
  !>
	FUNCTION n_up_x(sys,st,xmin,xmax)

		type(state), intent(in)   						::  st
		type(param), intent(in)							::  sys
		real(rl),intent(in),	optional					::  xmax,xmin
		real(rl)												::  n_up_x

		type(state)											::  upst
		complex(cx)											::  tmp,sum_i
		integer													::  imin, imax,n,m,i
		complex(cx)										   ::  fnx(st%np,sys%nmode)
		real(rl)												::  x

		upst = st
		upst%q(:) = 0._cx
		CALL normalise(upst)
		fnx = f_nx(sys,st)

		if (present(xmin)) then
			imin = int(xmin/sys%dx) + 1
			imax = int(xmax/sys%dx) + 1
		else
			imin = 1
			imax = sys%nmode
		end if

		tmp = 0._cx
		do n=1,st%np
			do m=1,st%np

				sum_i = 0._cx
				do i=imin,imax
					sum_i = sum_i + conjg( fnx(n,i) ) * fnx(m,i) 
				end do
				tmp = tmp + conjg(upst%p(n))*upst%p(m)*sum_i*upst%ov_ff(n,m)

			end do
		end do
		n_up_x = real( tmp )

	END FUNCTION

  !> -------------------------------------------------------------------------
  !> FUNCTION: n_up
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Calculates the total macroscopic number of photons associated exclusively 
  !>   with the qubit "up" state across the entire unwindowed waveguide.
  !> Arguments:
  !>   - sys : Parameter structure.
  !>   - st  : Current state.
  !> Return:
  !>   - real(rl) : Total integrated photon count in the "up" subspace.
  !>
	FUNCTION n_up(sys,st)

		type(state), intent(in)   						::  st
		type(param), intent(in)							::  sys
		real(rl)												::  n_up

		type(state)											::  upst
		complex(cx)											::  tmp,sum_k
		integer													::  n,m,k

		upst = st
		upst%q(:) = 0._cx
		CALL normalise(upst)

		tmp = 0._cx
		do n=1,st%np
			do m=1,st%np

				sum_k = sum( conjg( upst%f(n,:) ) * upst%f(m,:) )
				tmp = tmp + conjg(upst%p(n))*upst%p(m)*sum_k*upst%ov_ff(n,m)

			end do
		end do

		n_up = real( tmp )

	END FUNCTION

  !> -------------------------------------------------------------------------
  !> FUNCTION: f_k / h_k
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Extracts the effective global momentum-space probability amplitude for 
  !>   a specific mode `k` within the atomic "up" (`f_k`) or "down" (`h_k`) manifold. 
  !>   It properly traces over the non-orthogonal multi-polaron superposition, 
  !>   weighting the individual coherent displacements by their atomic amplitudes 
  !>   and bosonic overlaps.
  !> Arguments:
  !>   - st : Current state.
  !>   - k  : Momentum grid index.
  !> Return:
  !>   - complex(cx) : The total coherent amplitude at mode `k`.
  !>
	FUNCTION f_k(st,k)

		type(state), intent(in)  ::  st
		complex(cx)				  ::  f_k
		integer, intent(in)		  ::  k
		integer						  ::  n,m

		f_k=0._cx
		do n=1, size(st%p,1)
			do m=1, size(st%p,1)
				f_k = f_k + conjg(st%p(n))*st%p(m) * st%f(m,k) * st%ov_ff(n,m)
			end do
		end do

	END FUNCTION
	FUNCTION h_k(st,k)

		type(state), intent(in)  ::  st
		complex(cx)				  ::  tmp
		complex(cx)				  ::  h_k
		integer, intent(in)		  ::  k
		integer						  ::  n,m

		h_k = 0_cx
		do n=1, size(st%q,1)
			do m=1, size(st%q,1)
				h_k = h_k + conjg(st%q(n))*st%q(m) * st%h(m,k) * st%ov_hh(n,m)
			end do
		end do

	END FUNCTION

  !> -------------------------------------------------------------------------
  !> FUNCTION: f_nx / h_nx
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   High-performance bulk real-space transforms. Rather than evaluating the 
  !>   global field at a single point, these functions take the raw momentum 
  !>   matrices (`f`, `h`) and transform the entirety of all constituent polarons 
  !>   onto the discretized spatial grid simultaneously. Essential for calculating 
  !>   the Wigner function and the spatially-resolved photon density $n(x)$.
  !> Arguments:
  !>   - sys : Parameter structure (defines the spatial bounds and grid `dx`).
  !>   - st  : Current state.
  !> Return:
  !>   - complex(cx) : 2D array of real-space fields (N_cs $\times$ N_modes).
  !>
	FUNCTION f_x(sys,st,x)

		type(param), intent(in)  ::  sys
		type(state), intent(in)  ::  st
		complex(cx)				  ::  f_x
		real(rl), intent(in)	  ::  x
		integer						  ::  k

		f_x=0._cx
		do k=1,size(st%f,2)
			f_x = f_x + ( 2._rl/sqrt(2._rl*pi) )*sqrt(sys%dx*sys%dk(k)) * cos(sys%w(k)*x) * f_k(st,k)
		end do

	END FUNCTION
	FUNCTION h_x(sys,st,x)

		type(param), intent(in)  ::  sys
		type(state), intent(in)  ::  st
		complex(cx)				  ::  h_x
		real(rl), intent(in)	  ::  x
		integer						  ::  k

		h_x=0._cx
		do k=1,size(st%f,2)
			h_x = h_x + ( 2._rl/sqrt(2._rl*pi) )*sqrt(sys%dx*sys%dk(k)) * cos(sys%w(k)*x) * h_k(st,k)
		end do

	END FUNCTION

  !> -------------------------------------------------------------------------
  !> FUNCTION: f_nx / h_nx
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   High-performance bulk real-space transforms. Rather than evaluating the 
  !>   global field at a single point, these functions take the raw momentum 
  !>   matrices (`f`, `h`) and transform the entirety of all constituent polarons 
  !>   onto the discretized spatial grid simultaneously. Essential for calculating 
  !>   the Wigner function and the spatially-resolved photon density $n(x)$.
  !> Arguments:
  !>   - sys : Parameter structure (defines the spatial bounds and grid `dx`).
  !>   - st  : Current state.
  !> Return:
  !>   - complex(cx) : 2D array of real-space fields (N_cs $\times$ N_modes).
  !>
	FUNCTION f_nx(sys,st)

		type(param),intent(in)   						::  sys
		type(state),intent(in)    					::  st
		real(rl)				 	   					::  x
		complex(cx), dimension(st%np,sys%nmode)  ::  f_nx
		integer												::  n,i


		do n=1,st%np
			do i=1,sys%nmode
				x = (pi/sys%wmax) * (i-0.5_rl)
				f_nx(n,i) = ( 2._rl/sqrt(2._rl*pi) )*sqrt(sys%dx)*sum( sys%dk1 * cos(sys%w(:)*x)*st%f(n,:) )
			end do
		end do


	END FUNCTION
	FUNCTION h_nx(sys,st)

		type(param),intent(in)   						::  sys
		type(state),intent(in)    					::  st
		real(rl)				 	   					::  x
		complex(cx), dimension(st%np,sys%nmode)  ::  h_nx
		integer												::  n,i


		do n=1,st%np
			do i=1,sys%nmode
				x = (pi/sys%wmax) * (i-0.5)
				h_nx(n,i) = ( 2._rl/sqrt(2._rl*pi) )*sqrt(sys%dx)*sum( sys%dk(:)*cos(sys%w(:)*x)*st%h(n,:) )
			end do
		end do


	END FUNCTION

  !> -------------------------------------------------------------------------
  !> FUNCTION: factorial
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Standard mathematical utility calculating the integer factorial $n!$.
  !> Arguments:
  !>   - n : Input integer.
  !> Return:
  !>   - integer : The factorial result.
  !>
	FUNCTION factorial(n)

		integer, intent(in)		::  n
		integer						::  factorial
		integer						::  tmp, i 

		factorial=1
		if (n > 1) then
			do i=1,n
				factorial = factorial*i
			end do
		end if

	END FUNCTION


END MODULE SYSTEM
