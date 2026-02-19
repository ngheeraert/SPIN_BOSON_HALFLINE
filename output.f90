MODULE output
!==============================================================================
! MODULE OUTPUT
!
! High-level time evolution and diagnostics for the multipolaron simulation.
!
! Main responsibilities:
!   - Time integration (RK4) of the variational equations (calls system:CalcDerivatives).
!   - Adaptive control: monitor energy drift / error estimator and reduce timestep if needed.
!   - Optional adaptive basis growth: add a new coherent state (polaron) when Err(t) exceeds
!     a threshold (cf. paper's convergence discussion).
!   - Write observables (spin components, energy, photon distributions, Wigner function).
!==============================================================================


	USE system
	USE consts

	implicit none

CONTAINS

!------------------------------------------------------------------------------
! printTrajectory_HL
!   Entry point for a 'half-line' waveguide simulation:
!     - allocate the initial state & trajectory arrays,
!     - evolve until tmax (possibly with polaron additions),
!     - write summary observables to disk.
!------------------------------------------------------------------------------
	SUBROUTINE printTrajectory_HL(sys,st) 
		type(param), intent(in out)			::  sys
		type(state), intent(in out)			::  st
		type(traj)									::  tr
		real(rl), allocatable					::  Qf_x_var(:), Qf_p_var(:), Qf(:,:),Qf_x_normalised(:)
		real(rl)									::  factor, projected_value, rand_n,x,xmin,xmax,pmin,pmax,dx
		real(rl), allocatable					::  array_elements(:)
		integer										::  element_numb, index_i, i, index_f, numb_points, xnum, rand_times_numb
		integer :: date_values(1:8)

		!======================================================== 
		!-- Allocation and initialisation
		!=======================================================

		CALL allocate_state(sys,st)
		CALL allocate_trajectory(sys,tr)
		CALL initialise_vaccum(sys,st,0._rl)

		if (sys%verbose==1) then
			print*,"-- state initialisation complete"
		end if

		print*, parameterchar(sys)

		!========================================================
		!-- Evolving the state until sys%tmax
		!=======================================================

		!-- Evolution including the adding of polarons
		ADDING_IF: IF (sys%npadd .ne. 0) THEN

			ADDING_LOOP: DO

				CALL evolveState_HL(sys,st,tr,sys%tmax,1) !-- 1 specifies that routine should stop for adding pol

				if ( st%t >= sys%tmax ) exit adding_loop

				if (sys%npadd .ne. 0) then

					!CALL addpolaron(sys,st,sys%rtadd_arr(st%np - sys%npini + 1)) Spontaneous emisison
					CALL addpolaron(sys,st)
					tr%np_ar(tr%i_time) = st%np

				end if

				if ( st%np == (sys%npini + sys%npadd) ) exit adding_loop

			END DO ADDING_LOOP

		END IF ADDING_IF

		!-- Evolution without adding of polarons
		CALL evolveState_HL(sys,st,tr,sys%tmax, 0)   !-- 0 specifies no adding of polarons


		!========================================================
		!== HALF CHAIN FINAL COMPUTATIONS
		!========================================================

		!-- print out f(x) and f(k) at t=tmax
		CALL print_fks(sys,st,"  fst")
		CALL print_ps(sys,st,"  fst")

		CALL print_nk(sys,st)
		CALL print_nx(sys,st,"  fst")

		CALL print_evolution_data(sys,tr)
	END SUBROUTINE

!------------------------------------------------------------------------------
! evolveState_HL
!   Drive the time evolution loop up to tmax, storing observables in `tr` and calling
!   `addpolaron` when the TDVP error estimator exceeds `sys%merr_arr` after `sys%tref_arr`.
!------------------------------------------------------------------------------
	SUBROUTINE evolveState_HL(sys,st,tr,tmax, add) 
		type(param), intent(in)				::  sys
		real(rl), intent(in)					::  tmax
		type(state), intent(in out)	  		::  st
		type(traj), intent(in out)	  		::  tr
		integer,intent(in),optional			::  add
		type(state)								::  ost,oost,rot_st
		integer  									::  i,m,k,n,counter_n_plot,counter
		real(rl)									::  tini,phasediff,tcounter,dtprint,tprinted,slowfactor
		real(rl)									::  err_instant,err_reference, d_sigma_val, angle
		logical										::  err_ref_defined

		!-- initialise the states corresponding to previous time-steps
		ost = st
		oost = st

		tini = st%t

		!-- error monitoring variables
		err_instant = 0._rl
		err_reference = 1.0e8_rl
		err_ref_defined = .false. 

		!-- variables for printing interval
		tprinted = 0._rl
		counter_n_plot = 0
		counter=0

		EVOLUTION: DO

			if ( tmax <= st%t ) exit evolution

			if ( st%t >= tini + sys%tref_arr(st%np) &
				.and.  ( add == 1 ) &
				.and.  (err_instant >= sys%merr_arr( st%np - sys%npini + 1 ) ) ) then 
				!.and.  (err_instant -err_reference >= sys%merr_arr( st%np - sys%npini + 1 ) ) ) then 

				exit evolution

			end if

			!!-- ERROR REFERNCE VALUE
			!IF ( (st%t - tini > sys%tref_arr(st%np) ) .and. (.not. err_ref_defined) ) then

			!	CALL calcderivatives( sys,st )
			!	!CALL evolve_rk4( sys,st,ost,oost,1.0e4_rl )
			!	err_instant = real( error(sys,ost,st) )
			!	err_reference = err_instant

			!	err_ref_defined = .true.
			!	write( *,'(a30,f9.2)' ) "ERROR REFERENCE DEFINED AT t= ", st%t 

			!END IF

			slowfactor = 1._rl
			dtprint = 4*sys%dt/dble(slowfactor)

			if ((sys%npadd > 0) .and. (st%np>sys%npini)) then

				IF ( (st%t < 0.00001_rl+tini) ) then
					slowfactor = 1.0e5_rl
					dtprint = 30*sys%dt/dble(slowfactor)
				elseif ( (st%t < 0.0001_rl+tini) .and. (st%t > 0.00001_rl+tini) ) then
					slowfactor = 1.0e4_rl
					dtprint = 30*sys%dt/dble(slowfactor)
				elseif ( (st%t < 0.001_rl+tini) .and. (st%t > 0.0001_rl+tini) ) then
					slowfactor = 1.0e3_rl
					dtprint = 30*sys%dt/dble(slowfactor)
				elseif ( (st%t < 0.01_rl+tini) .and. (st%t > 0.001_rl+tini) ) then
					slowfactor = 1.0e2_rl
					dtprint = 30*sys%dt/dble(slowfactor)
				elseif ( (st%t < 0.05_rl+tini) .and. (st%t > 0.01_rl+tini) ) then
					slowfactor = 50
					dtprint = 30*sys%dt/dble(slowfactor)
				elseif ( (st%t < 6._rl+tini) .and. (st%t > 0.05_rl+tini) ) then
					slowfactor = 5
					dtprint = 30*sys%dt/dble(slowfactor)
				END IF

			end if

			!-- PRINTOUT
			PRINTOUT_IF: IF ( ((st%t-tprinted) > dtprint) .or. (st%t<1e-8) ) THEN

				!-- calculate 3 points with dt/10000, to obtain smooth curve for error
				!CALL evolve_rk4( sys,st,ost,oost,1.0e4_rl )
				CALL calcderivatives( sys,st )
				err_instant = real( error(sys,ost,st) )

				tr%time_ar(tr%i_time) = st%t
				tr%error_ar(tr%i_time) = err_instant
				tr%energy_ar(tr%i_time) = energy(sys,st)
				tr%norm_ar(tr%i_time) = norm(st)
				tr%spinXYZ_ar(tr%i_time,1) = sigmaX(st)
				tr%spinXYZ_ar(tr%i_time,2) = sigmaY(st)
				tr%spinXYZ_ar(tr%i_time,3) = sigmaZ(st)
				tr%i_time = tr%i_time + 1

				tprinted=st%t

			ENDIF PRINTOUT_IF

			!========================================================
			!== Evolve the system by one time-step

			CALL evolve_RK4(sys,st,ost,oost,slowfactor)
			if ( (sys%max_deltaE > 0) &
				.and. (ost%np == st%np ) &
				.and. (oost%np == st%np ) ) then
				CALL checkTimestep(sys,st,ost,oost,slowfactor,sys%max_deltaE)
			end if

			!========================================================

		END DO EVOLUTION
	END SUBROUTINE

!------------------------------------------------------------------------------
! print_evolution_data
!   Write diagnostic data to file for later plotting (field displacements, spin amplitudes,
!   photon-number distributions, etc.). Filenames are constructed from the parameter set.
!------------------------------------------------------------------------------
	SUBROUTINE print_evolution_data(sys,tr)
		type(param), intent(in)  		::  sys
		type(traj), intent(in)			::  tr
		character(len=200)		  		::  name_eren,name_spinXYZ,name_np,name_deltaR,name_ps
		real(rl)								::  t
		integer								::  i,n

		name_spinXYZ=trim(adjustl(sys%file_path))//"/spinXYZ_"//trim(adjustl(parameterchar(sys)))//".d"
		name_np=trim(adjustl(sys%file_path))//"/np_"//trim(adjustl(parameterchar(sys)))//".d"
		name_ErEN=trim(adjustl(sys%file_path))//"/ErEN_"//trim(adjustl(parameterchar(sys)))//".d"
		name_ps=trim(adjustl(sys%file_path))//"/ps_"//trim(adjustl(parameterchar(sys)))//".d"

		open (unit=10,file=name_ErEN,action="write",status="replace")
		open (unit=11,file= name_spinXYZ,action="write",status="replace")
		open (unit=12,file= name_np,action="write",status="replace")
		open (unit=13,file= name_ps,action="write",status="replace")

		do i=1, tr%i_time	-1

			t = tr%time_ar(i)
			write(10,'(4f25.10)') t, tr%error_ar(i), tr%energy_ar(i), tr%norm_ar(i)
			write(11,'(f25.15,3f25.10)') t, tr%spinXYZ_ar(i,1), tr%spinXYZ_ar(i,2), tr%spinXYZ_ar(i,3)
			write(12,*) t, tr%np_ar(i)

			write(13,'(f25.15)',advance='no') t
			do n=1,sys%npini + sys%npadd
				write(13,'(f25.15)',advance='no') tr%ps_ar(i,n)
			end do
			write(13,*)

		end do
		close(10)
		close(11)
		close(12)
		close(13)

	END SUBROUTINE

!------------------------------------------------------------------------------
! evolve_RK4
!   One integration step using classic 4th-order Runge-Kutta on the variational ODEs.
!   The parameter `slowfactor` can be used to effectively reduce dt (dt/slowfactor) when
!   the adaptive timestep control detects too large an energy drift.
!------------------------------------------------------------------------------
	SUBROUTINE Evolve_RK4(sys,st,ost,oost,slowfactor)
		type(state),intent(in out)   ::  st, ost, oost
		type(param),intent(in)   		::  sys
		real(rl), optional				::  slowfactor
		type(state)   					::  midSt
		real(rl)							::  dt

		oost = ost
		ost = st 
		midst = st

		if (present(slowfactor)) then
			dt = sys%dt/dble(slowfactor)
		else
			dt = sys%dt
		end if

		!== FIRST STEP OF RK4 ==========! 
		CALL calcderivatives(sys,st)

		st%f = st%f + st%fdot*dt/6._rl
		st%h = st%h + st%hdot*dt/6._rl
		st%p = st%p + st%pdot*dt/6._rl
		st%q = st%q + st%qdot*dt/6._rl
		CALL update_sums(sys,st)

		!== SECOND STEP OF RK4 ==========! 
		midSt%f = ost%f + 0.5_rl*dt*st%fdot
		midSt%h = ost%h + 0.5_rl*dt*st%hdot
		midSt%p = ost%p + 0.5_rl*dt*st%pdot
		midSt%q = ost%q + 0.5_rl*dt*st%qdot
		CALL update_sums(sys,midst)

		!-- for 1 mode only one cannot use the intrinsic fortran functions
		!-- hence we use here the manual summations in the case of 1 mode
		CALL calcderivatives(sys,midst)

		st%f = st%f + midst%fdot*dt/3._rl
		st%h = st%h + midst%hdot*dt/3._rl
		st%p = st%p + midst%pdot*dt/3._rl
		st%q = st%q + midst%qdot*dt/3._rl
		CALL update_sums(sys,st)

		!== THIRD STEP OF RK4 ==========! 
		midSt%f = ost%f + 0.5_rl*dt*midst%fdot 
		midSt%h = ost%h + 0.5_rl*dt*midst%hdot 
		midSt%p = ost%p + 0.5_rl*dt*midst%pdot 
		midSt%q = ost%q + 0.5_rl*dt*midst%qdot 
		CALL update_sums(sys,midst)

		CALL calcderivatives(sys,midst)

		st%f = st%f + midst%fdot*dt/3._rl
		st%h = st%h + midst%hdot*dt/3._rl
		st%p = st%p + midst%pdot*dt/3._rl
		st%q = st%q + midst%qdot*dt/3._rl
		CALL update_sums(sys,st)

		!== FOURTH STEP OF RK4 ==========! 
		midSt%f = ost%f + dt*midst%fdot 
		midSt%h = ost%h + dt*midst%hdot 
		midSt%p = ost%p + dt*midst%pdot 
		midSt%q = ost%q + dt*midst%qdot 
		CALL update_sums(sys,midst)

		CALL calcderivatives(sys,midst)

		st%f = st%f + midst%fdot*dt/6._rl
		st%h = st%h + midst%hdot*dt/6._rl
		st%p = st%p + midst%pdot*dt/6._rl
		st%q = st%q + midst%qdot*dt/6._rl
		CALL update_sums(sys,st)

		st%t = st%t + dt
	END SUBROUTINE evolve_RK4

	!-- routine to add polarons
!------------------------------------------------------------------------------
! addpolaron
!   Increase Ncs by 1 by reallocating the state and appending an extra coherent state
!   initialized near vacuum (small amplitude `sys%p0`). This is the adaptive basis growth
!   mechanism used to keep the TDVP error small while maintaining Ncs << Nmodes.
!------------------------------------------------------------------------------
	SUBROUTINE addpolaron(sys,st)
		type(param), intent(in)			::   sys
		type(state), intent(in out)     ::   st
		type(state)					 	::   ba_st
		real(rl)						::   a,b
		integer							::   n,k
		real(rl)						::   E1,E2,SZ1,SZ2,SX1,SX2

		CALL update_sums(sys,st)

		E1 = energy(sys,st)
		SZ1 = sigmaZ(st)
		SX1 = sigmaX(st)
		ba_st = st 

		!-- Re-allocate st with more polarons
		CALL allocate_state(sys,st,ba_st%np+1)

		st%t = ba_st%t

		do n=1, ba_st%np
			st%f(n,:) = ba_st%f(n,:)
			st%h(n,:) = ba_st%h(n,:)
			st%p(n) = ba_st%p(n)
			st%q(n) = ba_st%q(n)
		end do
		st%f(ba_st%np+1,:) = 0._cx
		st%h(ba_st%np+1,:) = 0._cx
		st%p(ba_st%np+1) = sys%p0
		st%q(ba_st%np+1) = sys%p0

		CALL update_sums(sys,st)
		!CALL normalise(st)

		if (sys%verbose==1) then
			print*,"================================="
			write(*,'(a21,I5,a4,I3,a6,f9.2)') "POLARON ADDED: from ",ba_st%np," to ",st%np, "at t=", st%t 
			E2 = energy(sys,st)
			SZ2 = sigmaZ(st)
			SX2 = sigmaX(st)
			print*, "delta( Sz ) = ",SZ2-SZ1
			print*, "delta( Sx ) = ",SX2-SX1
			print*, "delta( E ) = ",E2-E1
			print*,"================================="
			print*, " "
		end if

	END SUBROUTINE

!------------------------------------------------------------------------------
! checkTimestep
!   Practical stabilizer: if the variational integration produces an energy drift larger
!   than `errorlimit`, the routine rewinds two steps and re-evolves with a reduced dt
!   (increasing `slowfactor`) until the drift is acceptable or a max number of retries is hit.
!------------------------------------------------------------------------------
	SUBROUTINE checkTimestep(sys,st,ost,oost,slowfactor,errorlimit)
		type(param), intent(in)			:: sys
		type(state), intent(in out)   :: st, ost, oost
		real(rl), intent(in) 	      :: slowfactor
		integer, save       				:: cnt
		real(rl), intent(in)          :: errorlimit
		type(state)							:: st_tofix,ost_tofix,oost_tofix,best_st, best_ost
		integer                       :: fixtry
		real(rl)								:: best_deltaE, deltaE, ini_deltaE,slowfactorFix
		character(len=200)		  		::  name_fixes
		logical 								::  od


		st_tofix = st
		ost_tofix = ost
		oost_tofix = oost
		deltaE = energy(sys,st) - energy(sys,ost)
		ini_deltaE = deltaE
		best_deltaE = deltaE
		best_st = st
		best_ost = ost
		slowfactorFix = 1._rl
		fixtry = 0

		FIXING: do

			if ( (abs(deltaE) < errorlimit) .or. (fixtry > 20) ) exit fixing

			if (sys%verbose==1) then
				print*, "FIXING: t,deltaE=",st%t,deltaE
			end if

			fixtry = fixtry + 1 	!-- keep track of the number of tries
			st = oost_tofix				 	!-- rewind evolution by 2 steps

			!if (slowfactor < 100) then
			slowfactorFix = slowfactor*( fixtry+fixtry/dble(20) )  !-- first try to decrease dt
			!else 
			! slowfactorFix = dble(slowfactor)/dble(fixtry+1)  !-- first try to decrease dt
			!end if

			RE_EVOLVE: do

				if ( st_tofix%t+2*sys%dt  < st%t ) exit RE_EVOLVE
				if (sys%verbose==1) then
					write(*,'(I10)',advance='no') int( slowfactorFix )
				end if
				CALL evolve_RK4(sys,st,ost,oost,slowfactorFix)

			end do RE_EVOLVE

			!-- recalculate the energy error
			deltaE = (energy(sys,st) - energy(sys,ost_tofix))

			if (abs(deltaE) < abs(best_deltaE)) then
				best_deltaE = deltaE
				best_st = st
				best_ost = ost
			end if

			if (sys%verbose==1) then
				print*, "END FIXING"
			end if
		end do FIXING

		if ( (fixtry .ne. 0) .and. (abs(deltaE) < errorlimit) ) then

			write(*,*) "No fix for ",(cnt)," steps"
			write(*,*) "SUCCESSFUL, fix attempts:"&
				,(fixtry)," time=",(st%t) ," New_deltaE=",(deltaE)," Ini_deltaE=", (ini_deltaE)
			write(*,*)
			cnt=0

		else if ( (fixtry > 20) .and. (abs(deltaE) > errorlimit) ) then

			if (abs(deltaE) > 1.0e-5_rl) then

				if (sys%verbose==1) then
					print*,"ABORT -- error in the energy is too high"
				end if

				!print*,"DELETING ALL evolution files: ",trim(adjustl(parameterchar(sys)))
				!CALL system('rm data/*'//trim(adjustl(parameterchar(sys)))//'*')

				name_fixes=trim(adjustl(sys%file_path))//"/FAILED_"//trim(adjustl(parameterchar(sys)))//".d"
				open (unit=21,file= name_fixes,action="write",status="replace")
				write(21,*) "No fix for ",cnt," steps"
				write(21,*) "FAILED, new state set to best outcome: deltaE=",best_deltaE
				write(21,*)

				close(21)

				!== CLOSING THE PROGRAM
				!=======================
				stop
				!=======================

			end if

			write(*,*) "No fix for ",cnt," steps"
			write(*,*) "FAILED, new state set to best outcome: deltaE=",best_deltaE
			write(*,*)
			st = best_st
			ost = best_ost
			cnt=0
		else if ( (fixtry == 0) .and. (abs(deltaE) < errorlimit) ) then
			!-- count the number of steps betweeen 2 fixings
			cnt = cnt + 1
		end if
	END SUBROUTINE

	!== printing the state vector
!------------------------------------------------------------------------------
! print_fks
!   Write diagnostic data to file for later plotting (field displacements, spin amplitudes,
!   photon-number distributions, etc.). Filenames are constructed from the parameter set.
!------------------------------------------------------------------------------
	SUBROUTINE print_fks(sys,st,label, makefilename)
		type(param), intent(in)			::  sys
		type(state), intent(in)			::  st
		character(len=5), intent(in)    ::  label
		integer, intent(in), optional   ::  makefilename
		integer									::  k,i,min_index,max_index
		complex(cx), allocatable  		::  fnk(:,:), hnk(:,:)
		real(rl), allocatable				::   w(:)
		character(len=200)					::  name_fks

		if ( .not. present(makefilename) ) then
			name_fks=trim(adjustl(sys%file_path))//"/fks_"//trim(adjustl(label))//"_"//&
				trim(adjustl(parameterchar(sys)))//".d"
		else
			name_fks = " "//trim(adjustl(label))//" "
		end if
		open (unit=100,file= name_fks,action="write",status="replace")

		allocate(fnk(st%np,sys%nmode))
		allocate(hnk(st%np,sys%nmode))
		allocate(w(sys%nmode))
		fnk = st%f
		hnk = st%h
		w = sys%w
		min_index = 1
		max_index = sys%nmode

		do k = min_index, max_index
			write(100,'(f25.15)',advance='no') w(k)
			do i=1,st%np
				write(100,'(2f25.15)',advance='no') real( fnk(i,k) ), real( hnk(i,k) )
			end do
			do i=1,st%np
				write(100,'(2f25.15)',advance='no') aimag( fnk(i,k) ), aimag( hnk(i,k) )
			end do
			write(100,*)
		end do
		write(100,*)
		close(100)
	END SUBROUTINE

!------------------------------------------------------------------------------
! print_ps
!   Write diagnostic data to file for later plotting (field displacements, spin amplitudes,
!   photon-number distributions, etc.). Filenames are constructed from the parameter set.
!------------------------------------------------------------------------------
	SUBROUTINE print_ps(sys, st, label, makefilename)
		type(param), intent(in)			::  sys
		type(state), intent(in)			::  st
		character(len=5), intent(in)    ::  label
		integer, intent(in), optional   ::  makefilename
		integer									::  i
		character(len=200)					::  name_ps

		if ( .not. present(makefilename) ) then
			name_ps=trim(adjustl(sys%file_path))//"/ps_"//trim(adjustl(label))//"_"//&
				trim(adjustl(parameterchar(sys)))//".d"
		else
			name_ps = " "//trim(adjustl(label))//" "
		end if
		open (unit=100,file= name_ps,action="write",status="replace")

		do i=1,st%np
			write(100,'(2f25.15)',advance='no') real( st%p(i) ), real( st%q(i) )
		end do
		do i=1,st%np
			write(100,'(2f25.15)',advance='no') aimag( st%p(i) ), aimag( st%q(i) ) 
		end do
		write(100,*)

		close(100)
	END SUBROUTINE

!------------------------------------------------------------------------------
! print_fxs
!   Write diagnostic data to file for later plotting (field displacements, spin amplitudes,
!   photon-number distributions, etc.). Filenames are constructed from the parameter set.
!------------------------------------------------------------------------------
	SUBROUTINE print_fxs(sys,st,label,makefilename)
		type(param), intent(in)		 ::  sys
		type(state), intent(in)		 ::  st
		character(len=5), intent(in)	 ::  label
		integer, intent(in), optional ::  makefilename
		integer								 ::  n,i,min_index,max_index
		complex(cx), allocatable  	 ::  fnx(:,:), hnx(:,:)
		character(len=200)				 ::  name_fxs,tchar

		name_fxs=trim(adjustl(sys%file_path))//"/fxs_"//trim(adjustl(label))//"_"//&
			trim(adjustl(parameterchar(sys)))//".d"
		open (unit=100,file= name_fxs,action="write",status="replace")

		allocate(fnx(st%np,sys%nmode))
		allocate(hnx(st%np,sys%nmode))
		fnx = f_nx(sys,st)
		hnx = h_nx(sys,st)
		min_index = 1
		max_index = sys%nmode

		do i= min_index, max_index

			write(100,'(f25.15)',advance='no') sys%dx*dble(i-0.5_rl)
			do n=1, st%np
				write(100,'(2f25.15)',advance='no') real( fnx(n,i) ), real( hnx(n,i) )
			end do
			do n=1, st%np
				write(100,'(2f25.15)',advance='no') aimag( fnx(n,i) ), aimag( hnx(n,i) )
			end do
			write(100,*)

		end do
		write(100,*)
		close(100)
	END SUBROUTINE

	!== printing MEAN PHOTON NUMBERS
!------------------------------------------------------------------------------
! print_nk
!   Write diagnostic data to file for later plotting (field displacements, spin amplitudes,
!   photon-number distributions, etc.). Filenames are constructed from the parameter set.
!------------------------------------------------------------------------------
	SUBROUTINE print_nk(sys,st_in) 
		type(param), intent(in)			  	::  sys
		type(state), intent(in), optional	::  st_in
		type(state)								::	 st, st_filtered
		character(len=200)						::  wigxminchar,wigxmaxchar,filename,name1,name2
		integer										::  i

		if (present(st_in)) then
			st = st_in	
		else
			CALL allocate_state(sys,st,sys%npini+sys%npadd)
			CALL initialise_from_file(sys,st)
		end if

		write(wigxminchar, '(I5)') int(sys%wigxmin)
		write(wigxmaxchar, '(I5)') int(sys%wigxmax)
		filename=trim(adjustl(sys%file_path))//"/nk_end_"&
			//trim(adjustl(parameterchar(sys)))//"_"&
			//trim(adjustl(wigxminchar))//"_"//trim(adjustl(wigxmaxchar))//".d"
		open (unit=105,file= filename,action="write",status="replace")

		st_filtered = st
		st_filtered%f = f_nk_FT(sys,f_nx(sys,st),sys%wigxmin,sys%wigxmax)

		CALL print_fxs(sys,st_filtered,"  fil")
		CALL print_fks(sys,st_filtered,"  fil")

		do i=1,5 !sys%nmode
			write(105,'(f10.5,f14.8)') sys%w(i), n_up_k(sys,st_filtered,i)/sys%dk1
		end do
		close(105)
	END SUBROUTINE

!------------------------------------------------------------------------------
! print_nx
!   Write diagnostic data to file for later plotting (field displacements, spin amplitudes,
!   photon-number distributions, etc.). Filenames are constructed from the parameter set.
!------------------------------------------------------------------------------
	SUBROUTINE print_nx(sys,st,label) 
		type(param), intent(in)		 ::  sys
		type(state), intent(in)		 ::  st
		character(len=5), intent(in)	 ::  label
		type(state)						 ::  upst
		character(len=200)				 ::  filename
		integer								 ::  n,m,i
		complex(cx)						 ::  fnx( st%np,1:sys%nmode )
		complex(cx)						 ::  tmp

		filename=trim(adjustl(sys%file_path))//"/nx_"//trim(adjustl(label))//"_"&
			//trim(adjustl(parameterchar(sys)))//".d"
		open (unit=105,file= filename,action="write",status="replace")

		upst = st
		upst%q(:) = 0._cx
		CALL normalise(upst)

		fnx = f_nx(sys,st)

		tmp = 0._cx
		do i=1,sys%nmode
			do n=1,st%np
				do m=1,st%np

					tmp = tmp + conjg(upst%p(n))*upst%p(m)*(conjg(fnx(n,i))*fnx(m,i))*ov( fnx(n,:) , fnx(m,:) )

				end do
			end do
			write(105,'(f25.15,f25.15)') sys%dx*(dble(i)-0.5_rl), real(tmp)/sys%dx
			tmp=0._cx
		end do
		close(105)
	END SUBROUTINE

	!-- wigner function
!------------------------------------------------------------------------------
! wignerise / printwigner
!   Compute and export a phase-space Wigner function for the emitted wavepacket.
!   This is used to visualize Schrodinger cat structure (interference fringes).
!------------------------------------------------------------------------------
	SUBROUTINE wignerise(sys,st) 
		type(param), intent(in)			  	::  sys
		type(state), intent(in)		  		::  st
		character(len=200)						::  wigxminchar,wigxmaxchar,name_wigner

		write(wigxminchar, '(I5)') int(sys%wigxmin)
		write(wigxmaxchar, '(I5)') int(sys%wigxmax)
		name_wigner=trim(adjustl(sys%file_path))//"/wigner_"&
			//trim(adjustl(parameterchar(sys)))//"_"&
			//trim(adjustl(wigxminchar))//"_"//trim(adjustl(wigxmaxchar))//".d"
		open (unit=100,file= name_wigner,action="write",status="replace")

		CALL printwigner(100,sys,st,50)
		close(100)
	END SUBROUTINE
	
	SUBROUTINE printwigner(filenb,sys,st,res,xxminIN,xxmaxIN)
		type(param), intent(in)						 	::  sys
		type(state), intent(in) 		    				::  st
		integer, intent(in)								   ::  res
		integer,intent(in)									::  filenb
		real(rl), intent(in), optional					::  xxminIN, xxmaxIN
		type(state)					  						::  upst
		complex(cx), allocatable  						::  fnx(:,:)
		complex(cx), allocatable  						::  z(:)
		real(rl)												::  normz
		complex(cx)					  						::  wigner_f,tmp,wig_up,wig_down,zz
		integer							  						::  i,n,m,xnum,j
		real(rl)												::  xmin,xmax,pmin,pmax,x,p
		complex(cx), allocatable							::  ovmat(:,:)
		complex(cx), allocatable							::  zTfnx(:)
		complex(cx), allocatable							::  nullarray(:)
		real(rl)												::  xxmin, xxmax

		allocate( ovmat(st%np,st%np), zTfnx(st%np) )
		allocate( fnx(st%np,sys%nmode), z(sys%nmode), nullarray(sys%nmode) )
		nullarray = 0._cx

		if ( present(xxminIN) ) then
			xxmin = xxminIN
			xxmax = xxmaxIN
		else
			xxmin = sys%wigxmin
			xxmax = sys%wigxmax
		end if

		upst = st
		upst%q(:) = 0._cx
		CALL normalise(upst)

		fnx = f_nx(sys,st)

		z = 0._cx
		do i=int(xxmin/sys%dx)+1,int(xxmax/sys%dx)
			tmp = 0._cx
			do n=1,st%np
				tmp = tmp + upst%p(n)*fnx(n,i)
			end do
			z(i) = tmp
		end do

		normz = 0._rl
		normz = sqrt( sum( abs(z(:))**2 ) )
		if (normz > 1.0e-8) then
			z =(z/normz)
		end if

		do i=0,sys%nmode
			write(153,*) i*sys%dx, real(z(i)), aimag(z(i)) 
		end do

		ovmat = 0._cx
		zTfnx = 0._cx

		ovmat = st%ov_ff
		do n=1,size(upst%p,1)
			zTfnx(n) = sum(conjg(z(:))*fnx(n,:))
		end do

		xmin = -3.5
		xmax = 3.5
		pmin = -3.5
		pmax = 3.5
		xnum = (xmax - xmin) * res

		do i=0,xnum
			do j=0,xnum

				wigner_f = 0._rl
				wig_up = 0._cx
				wig_down = 0._cx
				x = xmin + (xmax-xmin)*i/xnum
				p = pmin + (pmax-pmin)*j/xnum

				zz =  ( x + Ic*p )

				do n=1,size(upst%p,1)
					do m=1,size(upst%p,1)
						zz = ( x + Ic*p )
						wig_up = wig_up + conjg(upst%p(n))*upst%p(m) &
							* exp( -2._rl * ( real(zz)+Ic*aimag(zz) - zTfnx(m) ) * ( real(zz)-Ic*aimag(zz) - conjg(zTfnx(n)) ) ) &
							* ovmat(n,m)
						zz = - zz
						wig_down = wig_down + conjg(upst%p(n))*upst%p(m) &
							* exp( -2._rl * ( real(zz)+Ic*aimag(zz) - zTfnx(m) ) * ( real(zz)-Ic*aimag(zz) - conjg(zTfnx(n)) ) ) &
							* ovmat(n,m)
					end do
				end do

				wigner_f = (2._rl/pi) * 0.5_rl * (wig_up + wig_down)

				write(filenb,*) x,p,real(wigner_f)

			end do
			write(filenb,*)
		end do
	END SUBROUTINE

	SUBROUTINE re_wignerise(sys) 
		type(param), intent(in)			  	::  sys
		type(state)						  		::  st
		character(len=200)						::  wigxminchar,wigxmaxchar,name_wigner

		CALL allocate_state(sys,st,sys%npini+sys%npadd)
		CALL initialise_from_file(sys,st)

		write(wigxminchar, '(I5)') int(sys%wigxmin)
		write(wigxmaxchar, '(I5)') int(sys%wigxmax)
		name_wigner=trim(adjustl(sys%file_path))//"/wigner_2_"&
			//trim(adjustl(parameterchar(sys)))//"_"&
			//trim(adjustl(wigxminchar))//"_"//trim(adjustl(wigxmaxchar))//".d"
		open (unit=100,file= name_wigner,action="write",status="replace")

		CALL update_sums(sys,st)
		CALL normalise(st)
		close(100)
	END SUBROUTINE

END MODULE output


