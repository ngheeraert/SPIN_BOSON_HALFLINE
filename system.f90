MODULE SYSTEM 

	USE consts
	USE lapackmodule
	USE inverse
	USE typedefs, only : cx => c_type, rl => r_type

	IMPLICIT NONE 

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

	END TYPE param

	TYPE STATE
		real(rl)					 	::  t 								!-- Time
		integer    				 	::  np
		complex(cx), allocatable  ::  f(:,:),h(:,:)   				!-- value of the fiel displacements and momenta
		complex(cx), allocatable  ::  p(:),q(:)   					!-- probability amplitudes of the polarons
		complex(cx), allocatable  ::  fdot(:,:),hdot(:,:)   		!-- the time derivatives
		complex(cx), allocatable  ::  pdot(:),qdot(:)   			!-- the time derivatives
		!-- for storing the large sums over k
		complex(cx), allocatable  ::  ov_ff(:,:), ov_hh(:,:)  	!-- matrices of overlaps
		complex(cx), allocatable  ::  ov_fh(:,:), ov_hf(:,:)  	!-- matrices of overlaps
		complex(cx), allocatable  ::  bigW_f(:,:), bigW_h(:,:)  !-- Ws (as in the notes)
		complex(cx), allocatable  ::  bigL_f(:,:), bigL_h(:,:)  !-- Ls (as in the notes)
	END TYPE STATE

	TYPE TRAJ
		integer							::  i_time
		real(rl), allocatable  	::  time_ar(:), energy_ar(:),  &
			norm_ar(:), error_ar(:), spinXYZ_ar(:,:), ps_ar(:,:) !,spinXYZ(:,:)
		integer, allocatable  		::  np_ar(:)
	END TYPE TRAJ

	COMPLEX(cx), PARAMETER :: Ic = ( 0._rl , 1._rl )
	REAL(rl),	   PARAMETER :: IP_p = 0.000005_rl 		!-- initialisation parameter


CONTAINS

	!== initialize the parameters and get those in the command line
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
		sys%merr = 0.0000001_rl
		sys%tref = 0.2_rl
		sys%dt = 0.05_rl
		sys%file_path = "data"

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
				else if(buffer=='-me')then
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
		sys%g(:) = sqrt( 2._rl*sys%alpha * sys%w(:) * sys%dk(:) * exp(-sys%w(:)/sys%wc) ) 

		!-- system variables that depend on a possible input
		sys%length = pi/sys%dk1
		sys%dx  = sys%length/sys%nmode
	END SUBROUTINE  

	!== Initialisation routines
	SUBROUTINE allocate_state(sys,st,npvalue)
		type(param), intent(in)      	::  sys
		type(state), intent(out)  		::  st
		type(state)							::  tmpst !-- if st has to be reallocated
		integer, intent(in),optional		::  npvalue
		integer									::  np, spinXYZ_dim, step_numb

		print*, "-- BEGINNING ARRAYS ALLOCATION"

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

		print*, "-- ARRAYS ALLOCATED"
	END SUBROUTINE

	SUBROUTINE allocate_trajectory(sys,tr)
		type(param), intent(in)      	::  sys
		type(traj), intent(out)      	::  tr
		integer									::  step_numb

		print*, "-- BEGINNING TRANJECTORY ALLOCATION"

		!spinXYZ_dim = (sys%tmax/sys%dt)*2._rl
		!allocate(tmpst%spinXYZ(spinXYZ_dim,7))
		step_numb = int( (sys%tmax / sys%dt)*1.1 ) 

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

		print*, "-- TRAJECTORY ALLOCATED"
	END SUBROUTINE

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
			print*, "=="
			print*, "STATE PREPARED IN UP + DOWN"
			print*, "=="
		else if ( (sys%prep == 2) ) then
			!-- prepare in the up state
			st%p(1) = 0.99999999999999_rl
			st%q(1) = sqrt( 1._rl - conjg(st%p(1))*st%p(1) )
			print*, "=="
			print*, "STATE PREPARED IN UP"
			print*, "=="
		else if ( sys%prep == 3 ) then
			!-- prepare in the up state
			st%q(1) = 0.99999999999999_rl
			st%p(1) = sqrt( 1._rl - conjg(st%q(1))*st%q(1) )
			print*, "=="
			print*, "STATE PREPARED IN DOWN"
			print*, "=="
		else 
			print*, "Error in the value of prep"
		end if

		!-- updating the sums over k
		CALL update_sums(sys,st)
		CALL normalise(st)
	END SUBROUTINE

	SUBROUTINE initialise_from_file(sys,st)
		type(param), intent(in)			  	::  sys
		type(state), intent(in out)	  		::  st
		integer  									::  i,j,m,k
		character(len=200)						::  fks_file,ps_file
		real(rl)									::  f_r,f_i,h_r,h_i,p_r,q_r,p_i,q_i,a

		print*, "Initialising from: ", parameterchar(sys)
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
		write(p0char, '(f6.3)') sys%p0*1000._rl
		write(merrchar, '(f6.3)') sys%merr*1000000._rl
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
			"_p"//trim(adjustl(prepchar))
	END FUNCTION

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

	!=======================================
	!== Calculation of the derivatives
	!======================================

	!-- Calculate the state derivatives form the st variables
	SUBROUTINE calcDerivatives_OLD(sys,st) 
		type(param), intent(in)                            		::  sys
		type(state), intent(in out)	                      		::  st
		complex(cx), dimension(st%np)							 		::  bigP, bigQ, v_p, v_q
		complex(cx), dimension(st%np,sys%nmode)				 		::  bigF, bigH, v_f, v_h
		complex(cx), dimension(st%np,st%np)              		::  M_p, M_q, M_pi, M_qi,N_p, N_q, N_pi, N_qi
		complex(cx), dimension(st%np,st%np,st%np)       		::  A_p, A_q
		complex(cx), dimension(st%np,st%np,st%np,sys%nmode)  ::  A_f, A_h
		real(rl), dimension(2*st%np**2,2*st%np**2)   	 		::  EqMat_p, EqMat_q
		real(rl), dimension(2*st%np**2)    		 			 		::  EqRhs_p, EqRhs_q, kappaVec_p, kappaVec_q
		real(rl), dimension(st%np,st%np)	   				 		::  Kappa_p_r,Kappa_q_r,kappa_p_c,kappa_q_c
		complex(cx), dimension(st%np,st%np)					 		::  Kappa_p,Kappa_q  !-- FINAL KAPPA MATRICES
		complex(cx),dimension(st%np,sys%nmode)  			 		::  f,h,fc,hc  		
		complex(cx),dimension(st%np)			   			 		::  p,q,pc,qc
		complex(cx),dimension(st%np)			   			 		::  der_E_pc_save,der_E_qc_save
		complex(cx),dimension(st%np,sys%nmode)			   	::  der_E_fc_save,der_E_hc_save
		integer												          		::  i,j,k,l,m,cj,ii,jj,np,info,s

		complex(cx), dimension(st%np) 				          		::  tmpV1_p,tmpV1_q,tmpV2_p,tmpV2_q
		complex(cx), dimension(st%np,st%np,st%np,sys%nmode) 	::  tmpA_f1,tmpA_h1,tmpA_f2,tmpA_h2
		complex(cx), dimension(st%np,sys%nmode) 				   ::  tmpV1_f,tmpV1_h
		complex(cx), dimension(st%np,st%np,st%np,st%np)		::  D_f,D_h
		complex(cx), dimension(st%np,st%np,st%np)				::  E_f,E_h
		complex(cx), dimension(st%np,st%np)						::  K_f,K_h


		!-- A few shortcuts
		f=st%f;h=st%h;p=st%p;q=st%q
		fc=conjg(f);hc=conjg(h);pc=conjg(p);qc=conjg(q)
		np = st%np

		!- set all variables to zero
		bigP = 0._cx;bigQ = 0._cx;bigF = 0._cx;bigH = 0._cx
		M_p=0._cx;M_q=0._cx;N_p=0._cx;N_q=0._cx
		M_pi=0._cx;M_qi=0._cx;N_pi=0._cx;N_qi=0._cx
		A_p=0._cx;A_q=0._cx;A_f=0._cx;A_h=0._cx
		eqMat_p=0._cx;eqMat_q=0._cx;eqrhs_p=0._cx;eqrhs_q=0._cx
		eqRhs_p = 0._cx;eqRhs_q = 0._cx
		v_p=0._cx;v_q=0._cx;v_f=0._cx;v_h=0._cx
		Kappa_p_r=0._cx;Kappa_q_r=0._cx;kappa_p_c=0._cx;kappa_q_c=0._cx
		kappaVec_p=0._cx; kappaVec_q=0._cx;Kappa_p=0._cx;Kappa_q=0._cx
		der_E_fc_save(:,:) = 0._cx
		der_E_hc_save(:,:) = 0._cx
		der_E_pc_save(:) = 0._cx
		der_E_qc_save(:) = 0._cx
		D_f = 0._cx;D_h = 0._cx
		E_f = 0._cx;E_h = 0._cx
		K_f = 0._cx;K_h = 0._cx

		!==================================================
		!-- STEP 1: Computation of A_p,A_q and v_p,v_q
		!-- pDot(i) = sum_(j,k) [ A_p(i,j,k)*Kappa(k,j) + v_p(i) ] 
		!-- qDot(i) = sum_(j,k) [ A_q(i,j,k)*Kappa(k,j) + v_q(i) ] 

		do i=1,np
			bigP(i) = P_j(sys,st,i)
			bigQ(i) = Q_j(sys,st,i)
			bigF(i,:) =  F_j(sys,st,i) 
			bigH(i,:) =  H_j(sys,st,i) 
		end do

		!-- Matrix M: to be inverted
		M_p = st%ov_ff
		M_q = st%ov_hh

		!-- perform the inversion
		M_pi = M_p
		M_qi = M_q
		CALL invertH(M_pi,info)
		CALL invertH(M_qi,info)

		!-- Define v_p, v_q and A_p, A_q
		v_p = M_pi .matprod. bigP
		v_q = M_qi .matprod. bigQ

		do i=1,size(A_p,1)
			do j=1,size(A_p,2)
				do k=1,size(A_p,3)
					A_p(i,j,k) = 0.5_rl* M_pi(i,j) * st%ov_ff(j,k) * p(k)
					A_q(i,j,k) = 0.5_rl* M_qi(i,j) * st%ov_hh(j,k) * q(k)
				end do
			end do
		end do


		!==================================================
		!-- STEP 2: Computation of A_f,A_h and v_f,v_h
		!-- fDot(i) = sum_(j,k) [ A_f(i,j,k)*Kappa(k,j) ] + v_f(i)  
		!-- hDot(i) = sum_(j,k) [ A_h(i,j,k)*Kappa(k,j) ] + v_h(i)  

		do j=1,np
			do m=1,np
				N_p(j,m) = p(m)*st%ov_ff(j,m)
				N_q(j,m) = q(m)*st%ov_hh(j,m) 
			end do
		end do

		!-- Calculate Nij and is inverse
		N_pi = N_p
		N_qi = N_q
		CALL invertGeneral(N_pi,info)
		CALL invertGeneral(N_qi,info)


		!-- Computation of v_f and v_h (H_i in the notes)
		tmpV1_f=0._cx
		tmpV1_h=0._cx

		!-- First term of v_f (and v_h)
		tmpV1_f = bigF
		tmpV1_h = bigH
		!-- Second term of v_f
		do s=1,sys%nmode
			do j=1,np
				tmpV1_f(j,s) = tmpV1_f(j,s) &
					- sum( st%ov_ff(j,:)*v_p(:)*f(:,s) )
				tmpV1_h(j,s) = tmpV1_h(j,s) &
					- sum( st%ov_hh(j,:)*v_q(:)*h(:,s) )
			end do
			v_f(:,s) = N_pi .matprod. tmpV1_f(:,s)
			v_h(:,s) = N_qi .matprod. tmpV1_h(:,s)
		end do

		tmpA_f1=0._cx; tmpA_h1=0._cx
		tmpA_f2=0._cx; tmpA_h2=0._cx

		!-- Defining the Betref matrices, here designated by A_f and A_h

		do s=1,sys%nmode
			do k=1,size(A_f,2)
				do j=1,size(A_f,3)

					do i=1,size(A_f,1)
						tmpA_f1(i,j,k,s) = 0.5_rl* N_pi(i,j) * p(k)*f(k,s) * st%ov_ff(j,k)
						tmpA_h1(i,j,k,s) = 0.5_rl* N_qi(i,j) * q(k)*h(k,s) * st%ov_hh(j,k)
					end do

					do l=1,np
						tmpA_f2(l,j,k,s) = tmpA_f2(l,j,k,s) & 
							- sum( A_p(:,j,k)*f(:,s)*st%ov_ff(l,:) )
						tmpA_h2(l,j,k,s) = tmpA_h2(l,j,k,s) &
							- sum( A_q(:,j,k)*h(:,s)*st%ov_hh(l,:) )
					end do

				end do


				A_f(:,:,k,s) = tmpA_f1(:,:,k,s) + ( N_pi .matprod. tmpA_f2(:,:,k,s) )
				A_h(:,:,k,s) = tmpA_h1(:,:,k,s) + ( N_qi .matprod. tmpA_h2(:,:,k,s) )

			end do
		end do

		!==================================================
		!-- STEP 3: Solve the system of equations for Kappa
		!-- Define Q_ij, a 2*np**2 matrix: j in [1,np**2] correpsonds to Re(kappa)
		!--    while j in [1+np**2,2*np**2] corresponds to Im(Kappa)

		cj = st%np**2   !-- a shortcut to get to the imagainary part of kappaVec

		do i=1,st%np
			do j=1,st%np
				K_f(i,j) = sum( conjg(f(i,:))*v_f(i,:) + f(i,:)*conjg(v_f(i,:)) - 2._rl*conjg(f(j,:))*v_f(i,:) )
				K_h(i,j) = sum( conjg(h(i,:))*v_h(i,:) + h(i,:)*conjg(v_h(i,:)) - 2._rl*conjg(h(j,:))*v_h(i,:) )
				do k=1,st%np
					E_f(i,j,k) = sum( f(i,:)*conjg(A_f(i,j,k,:))  )
					E_h(i,j,k) = sum( h(i,:)*conjg(A_h(i,j,k,:))  )
					do m=1,st%np
						D_f(i,m,j,k) = sum( (conjg(f(i,:)) - 2._rl*conjg(f(m,:)))*A_f(i,j,k,:) )
						D_h(i,m,j,k) = sum( (conjg(h(i,:)) - 2._rl*conjg(h(m,:)))*A_h(i,j,k,:) )
					end do
				end do
			end do
		end do


		do i=1, np
			do m=1, np

				ii = np*(i-1)+m

				eqRhs_p(ii) 	  =  real( K_f(i,m) )
				eqRhs_q(ii) 	  =  real( K_h(i,m) )
				eqRhs_p(ii+cj)  =  aimag( K_f(i,m) )
				eqRhs_q(ii+cj)  =  aimag( K_h(i,m) )

				do j=1, np
					do k=1, np

						jj = np*(k-1) + j

						eqMat_p(ii,jj) = kroneckerDelta(ii,jj) - real(  D_f(i,m,j,k) + E_f(i,j,k) )
						eqMat_p(ii,jj+cj) = kroneckerDelta(ii,jj+cj) + aimag(  D_f(i,m,j,k) - E_f(i,j,k) )
						eqMat_p(ii+cj,jj) = kroneckerDelta(ii+cj,jj) - aimag(  D_f(i,m,j,k) + E_f(i,j,k) )
						eqMat_p(ii+cj,jj+cj) = kroneckerDelta(ii+cj,jj+cj) - real(  D_f(i,m,j,k) - E_f(i,j,k) )

						eqMat_q(ii,jj) = kroneckerDelta(ii,jj) - real(  D_h(i,m,j,k) + E_h(i,j,k) )
						eqMat_q(ii,jj+cj) = kroneckerDelta(ii,jj+cj) + aimag(  D_h(i,m,j,k) - E_h(i,j,k) )
						eqMat_q(ii+cj,jj) = kroneckerDelta(ii+cj,jj) - aimag(  D_h(i,m,j,k) + E_h(i,j,k) )
						eqMat_q(ii+cj,jj+cj) = kroneckerDelta(ii+cj,jj+cj) - real(  D_h(i,m,j,k) - E_h(i,j,k) )

					end do
				end do
			end do
		end do


		!-- Apply the Lapack algorithm to solve the real system of equations
		!-- the solution replaces the 2nd argument
		kappaVec_p = eqRhs_p
		CALL solveEq_r(eqMat_p,kappaVec_p)
		kappaVec_q = eqRhs_q
		CALL solveEq_r(eqMat_q,kappaVec_q)


		! Convert the matrices to np*np matrix
		kappa_p_r = transpose(reshape(kappaVec_p(1:np**2),(/np,np/)) )
		kappa_p_c = transpose(reshape(kappaVec_p(1+np**2:2*np**2),(/np,np/)) )
		kappa_q_r = transpose(reshape(kappaVec_q(1:np**2),(/np,np/)) )
		kappa_q_c = transpose(reshape(kappaVec_q(1+np**2:2*np**2),(/np,np/)) )

		kappa_p = kappa_p_r + Ic * kappa_p_c
		kappa_q = kappa_q_r + Ic * kappa_q_c

		!-- Compute the pdots and qdots
		!-- pdot(i) = sum_j,k 0.5*M_pi(i,j)*ov(f(j),f(k))*p(k)*Kappa_p(k,j) + v_p(i)

		tmpV1_p=0._rl;tmpV2_p=0._rl
		tmpV1_q=0._rl;tmpV2_q=0._rl
		do j=1,np
			tmpV1_p(j)=tmpV1_p(j) + 0.5_rl*sum( p(:)*st%ov_ff(j,:)*Kappa_p(:,j) )
			tmpV1_q(j)=tmpV1_q(j) + 0.5_rl*sum( q(:)*st%ov_hh(j,:)*Kappa_q(:,j) )
		end do

		tmpV2_p = M_pi .matprod. tmpV1_p
		tmpV2_q = M_qi .matprod. tmpV1_q

		st%pdot = tmpV2_p + v_p
		st%qdot = tmpV2_q + v_q


		!-- Compute the fdots and hdots
		!-- fdot(i) = sum_j,k [ -N_pi(i,j)*p(k)*Kappa_p(k,j)*pc(j)*f(k)*ov(f(j),f(k)) 
		!								+ sum_l,m [ -N_pi(i,l)*A_p(m,j,k)*Kappa_p(k,j)*pc(l)f(m)*ov(f(l),f(m)) ] ] + v_f(i)
		tmpV1_f=0._cx
		tmpV1_h=0._cx

		do s=1, sys%nmode
			do i=1,np
				do j=1,np
					tmpV1_f(i,s) = tmpV1_f(i,s) + sum( A_f(i,j,:,s) * Kappa_p(:,j) )
					tmpV1_h(i,s) = tmpV1_h(i,s) + sum( A_h(i,j,:,s) * Kappa_q(:,j) )
				end do
			end do
		end do

		st%fdot = tmpV1_f + v_f
		st%hdot = tmpV1_h + v_h
	END SUBROUTINE calcDerivatives_OLD

	!-- Calculate the state derivatives form the st variables with KRYLOV
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


		!print*, '-----------------------------------'
		!print*, "fastCalcDerivatives used ! number : "
		! print*, counter
		! if (counter > 100000) then
		!   stop "counter>100000 !"
		! end if

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


	!-- Derivatives of E with respect to sepcified variable STARRED
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

		energy = real(tmp)

	END FUNCTION
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
	SUBROUTINE normalise(st)

		type(state), intent(in out)  ::  st
		real(rl)							:: normval

		normval=norm(st)

		st%p = st%p/normval
		st%q = st%q/normval

	END SUBROUTINE
	FUNCTION error(sys,oost,ost,st)

		type(param),intent(in)		::  sys
		type(state),intent(in) 	::  oost,ost,st
		complex(cx)					::  tmp1,tmp2, tmp3, tmp4
		complex(cx)	   			::  error
		complex(cx), dimension(ost%np,sys%nmode)  ::  off, ohh, ofh, ohf,ofdd,ohdd
		complex(cx), dimension(ost%np)  ::  opdd,oqdd
		integer					   	::  i,j

		tmp1 = 0._cx
		tmp2 = 0._cx
		tmp3 = 0._cx
		tmp4 = 0._cx
		error = 0._cx

		ofdd(:,:) = (st%fdot(:,:) - oost%fdot(:,:))/(st%t-oost%t)
		ohdd(:,:) = (st%hdot(:,:) - oost%hdot(:,:))/(st%t-oost%t)
		opdd(:) = (st%pdot(:) - oost%pdot(:))/(st%t-oost%t)
		oqdd(:) = (st%qdot(:) - oost%qdot(:))/(st%t-oost%t)

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

				tmp4 = tmp4 + conjg(ost%p(j))*ost%ov_ff(j,i)*( &
					+ opdd(i) &
					- ost%pdot(i)*( off(i,i) + conjg(off(i,i)) - 2_rl*conjg(off(i,j)) ) &
					- 0.5_rl*ost%p(i)*( sum( 2_rl*real(conjg(ofdd(i,:))*ost%f(i,:)) &
					+ 2_rl*conjg(ost%fdot(i,:))*ost%fdot(i,:) &
					- 2_rl*ofdd(i,:)*conjg(ost%f(j,:)) ) ) &
					+ 0.25_rl*ost%p(i)*( off(i,i) + conjg(off(i,i)) - 2_rl*conjg(off(i,j)) )**2 &
					)


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

				tmp4 = tmp4 + conjg(ost%q(j))*ost%ov_hh(j,i)*( &
					+ oqdd(i) &
					- ost%qdot(i)*( ohh(i,i) + conjg(ohh(i,i)) - 2_rl*conjg(ohh(i,j)) ) &
					- 0.5_rl*ost%q(i)*( sum( 2_rl*real(conjg(ohdd(i,:))*ost%h(i,:)) &
					+ 2_rl*conjg(ost%hdot(i,:))*ost%hdot(i,:) &
					- 2_rl*ohdd(i,:)*conjg(ost%h(j,:)) ) ) &
					+ 0.25_rl*ost%q(i)*( ohh(i,i) + conjg(ohh(i,i)) - 2_rl*conjg(ohh(i,j)) )**2 &
					)

			end do
		end do

		error =  -0.5_rl*real(tmp4) + 0.5_rl*tmp1 - 2._rl*aimag(tmp2) + tmp3

	END FUNCTION	

	!-- Calculate the overlap between two coherent states
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
	FUNCTION ov_scalar(f1,f2)

		complex(cx), intent(in) :: f1, f2
		complex(cx)             :: ov_scalar
		complex(cx)				  :: tmp1, tmp2, tmp3

		tmp1 = conjg(f1)*f1
		tmp2 = conjg(f2)*f2
		tmp3 = conjg(f1)*f2

		ov_scalar = exp( -0.5_rl*tmp1 - 0.5_rl*tmp2 + tmp3 ) 

	END FUNCTION	ov_scalar
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

	!-- Probability of being in the up or down state
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

	!-- Expectation value of sigmaX
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
	FUNCTION sigmaZ(st)

		type(state), intent(in)  :: st
		real(rl) 					  :: sigmaZ
		integer						  :: n,m

		sigmaZ = upProb(st) - downProb(st)

	END FUNCTION
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

	!======================================================
	!== INTERPOLATION ROUTINES
	!======================================================

	!======================================================
	!== HALF LINE FUNCITONS
	!======================================================
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

	!======================================================
	!== MATH FUNCTIONS
	!======================================================

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


!
!  SUBROUTINE de_interpolate(sys,st)
!
!	 type(param),intent(in out)   ::  sys
!	 type(state),intent(in out)   ::  st
!	 type(param)				   	::  sys_dint
!	 type(state)  				   	::  st_dint
!	 integer								::  i, m, i_int, nmode1, nmode_rest
!
!	 m = sys%dk_ratio
!	 sys_dint = sys
!
!	 sys_dint%nmode = sys%nmode1 + (sys%nmode-sys%nmode1)/m
!	 nmode_rest = sys%nmode - sys%nmode1 - ((sys%nmode-sys%nmode1)/m)*m
!	 print*,"NMODE_REST=",nmode_rest
!
!	 deallocate(sys_dint%w)
!	 deallocate(sys_dint%dk)
!	 deallocate(sys_dint%g)
!	 allocate(sys_dint%w(sys_dint%nmode))
!	 allocate(sys_dint%dk(sys_dint%nmode))
!	 allocate(sys_dint%g(sys_dint%nmode))
!
!	 CALL allocate_state(sys_dint,st_dint,st%np)
!	 st_dint%p = st%p
!	 st_dint%q = st%q
!	 st_dint%t = st%t
!
!	 !-- the first nmode1 terms are indentical
!	 sys_dint%w(1:sys%nmode1) = sys%w(1:sys%nmode1)
!	 sys_dint%dk(1:sys%nmode1) = sys%dk(1:sys%nmode1)
!	 sys_dint%g(1:sys%nmode1) = sys%g(1:sys%nmode1)
!	 st_dint%f(:,1:sys%nmode1) = st%f(:,1:sys%nmode1)
!	 st_dint%fo(:,1:sys%nmode1) = st%fo(:,1:sys%nmode1)
!	 st_dint%h(:,1:sys%nmode1) = st%h(:,1:sys%nmode1)
!	 st_dint%ho(:,1:sys%nmode1) = st%ho(:,1:sys%nmode1)
!
!	 !-- for the first de_interpollation, just move in the array by 1/2 of the
!	 !			new dk
!	 sys_dint%w(sys%nmode1+1) = sys%w(sys%nmode1+(m+1)/2)
!	 sys_dint%dk(sys%nmode1+1) = sys%dk(sys%nmode1+(m+1)/2)*m
!	 sys_dint%g(sys%nmode1+1) = sys%g(sys%nmode1+(m+1)/2)*sqrt(dble(m))
!	 st_dint%f(:,sys%nmode1+1) = st%f(:,sys%nmode1+(m+1)/2)*sqrt(dble(m))
!	 st_dint%fo(:,sys%nmode1+1) = st%fo(:,sys%nmode1+(m+1)/2)*sqrt(dble(m))
!	 st_dint%h(:,sys%nmode1+1) = st%h(:,sys%nmode1+(m+1)/2)*sqrt(dble(m))
!	 st_dint%ho(:,sys%nmode1+1) = st%ho(:,sys%nmode1+(m+1)/2)*sqrt(dble(m))
!
!	 !-- now move by dk, except the last
!	 do i=1, sys_dint%nmode  -sys%nmode1 -1 -1
!	 	i_int = sys%nmode1+(m+1)/2+m*i
!	   sys_dint%w(sys%nmode1+1+i) = sys%w(i_int)
!	   sys_dint%dk(sys%nmode1+1+i) = sys%dk(i_int)*m
!	   sys_dint%g(sys%nmode1+1+i) = sys%g(i_int)*sqrt(dble(m))
!	   st_dint%f(:,sys%nmode1+1+i) = st%f(:,i_int)*sqrt(dble(m))
!	   st_dint%fo(:,sys%nmode1+1+i) = st%fo(:,i_int)*sqrt(dble(m))
!	   st_dint%h(:,sys%nmode1+1+i) = st%h(:,i_int)*sqrt(dble(m))
!	   st_dint%ho(:,sys%nmode1+1+i) = st%ho(:,i_int)*sqrt(dble(m))
!	 end do
!
!	 !-- int the last we need to incorporate the modes that are left: nmode_rest
!	 i=sys_dint%nmode  -sys%nmode1 -1
!	 i_int = sys%nmode1+(m+1)/2+m*i
!	 sys_dint%w(sys_dint%nmode) = sys%w(i_int)
!	 sys_dint%dk(sys_dint%nmode) = sys%dk(i_int)*(m+nmode_rest)
!	 sys_dint%g(sys_dint%nmode) = sys%g(i_int)*sqrt(dble(m+nmode_rest))
!	 st_dint%f(:,sys_dint%nmode) = st%f(:,i_int)*sqrt(dble(m+nmode_rest))
!	 st_dint%fo(:,sys_dint%nmode) = st%fo(:,i_int)*sqrt(dble(m+nmode_rest))
!	 st_dint%h(:,sys_dint%nmode) = st%h(:,i_int)*sqrt(dble(m+nmode_rest))
!	 st_dint%ho(:,sys_dint%nmode) = st%ho(:,i_int)*sqrt(dble(m+nmode_rest))
!
!
!	 CALL update_sums(sys_dint,st_dint)
!	 CALL normalise(st_dint)
!
!	 print*,"ARRAYS DE-INTERPOLATED"
!	 print*,"-- number of modes =", sys%nmode1 ," / ", sys_dint%nmode - sys%nmode1
!	 print*,"-- dk =", sys%dk(1) ," / ", sys_dint%dk(sys%nmode1+1)
!
!	 sys = sys_dint
!	 st = st_dint
!
!  END SUBROUTINE
!  SUBROUTINE interpolate(sys,st,st_int)
!
!	 type(param),intent(in)   		::  sys
!	 type(state),intent(in)    	::  st
!	 type(state),intent(in out)   ::  st_int
!	 integer								::  i, j, m, i_int
!
!	 m = sys%dk_ratio
!
!	 CALL allocate_state(sys,st_int,st%np)
!
!	 st_int%p = st%p
!	 st_int%q = st%q
!	 st_int%t = st%t
!
!	 !-- first part of the arrays is the same
!	 st_int%f(:,1:sys%nmode1) = st%f(:,1:sys%nmode1)
!	 st_int%fo(:,1:sys%nmode1) = st%fo(:,1:sys%nmode1)
!	 st_int%h(:,1:sys%nmode1) = st%h(:,1:sys%nmode1)
!	 st_int%ho(:,1:sys%nmode1) = st%ho(:,1:sys%nmode1)
!
!	 !-- first point to be interpollated
!	 st_int%f(:,sys%nmode1+(m+1)/2) = st%f(:,sys%nmode1+1)
!	 st_int%fo(:,sys%nmode1+(m+1)/2) = st%fo(:,sys%nmode1+1)
!	 st_int%h(:,sys%nmode1+(m+1)/2) = st%h(:,sys%nmode1+1)
!	 st_int%ho(:,sys%nmode1+(m+1)/2) = st%ho(:,sys%nmode1+1)
!
!	 i_int = sys%nmode1+(m+1)/2
!	 do j=1,(m-1)/2
!
!		st_int%f(:,i_int-j) = st%f(:,sys%nmode1+1) &
!			 + (st%f(:,sys%nmode1)*sqrt(dble(m)) - st%f(:,sys%nmode1+1))*dble(j)/dble( (m-1)/2+1 )
!		st_int%fo(:,i_int-j) = st%fo(:,sys%nmode1+1) &
!			 + (st%fo(:,sys%nmode1)*sqrt(dble(m)) - st%fo(:,sys%nmode1+1))*dble(j)/dble( (m-1)/2+1 )
!		st_int%h(:,i_int-j) = st%h(:,sys%nmode1+1) &
!			 + (st%h(:,sys%nmode1)*sqrt(dble(m)) - st%h(:,sys%nmode1+1))*dble(j)/dble( (m-1)/2+1 )
!		st_int%ho(:,i_int-j) = st%ho(:,sys%nmode1+1) &
!			 + (st%ho(:,sys%nmode1)*sqrt(dble(m)) - st%ho(:,sys%nmode1+1))*dble(j)/dble( (m-1)/2+1 )
!
!	 end do
!
!	 !-- interpollate all the middle terms except the lst of the input array
!	 do i=1, size(st%f,2) - sys%nmode1-1 -1
!
!	 	i_int = sys%nmode1+(m+1)/2+m*i
!	   st_int%f(:,i_int) = st%f(:,sys%nmode1+1+i) 
!	   st_int%fo(:,i_int) = st%fo(:,sys%nmode1+1+i)
!	   st_int%h(:,i_int) = st%h(:,sys%nmode1+1+i)
!	   st_int%ho(:,i_int) = st%ho(:,sys%nmode1+1+i)
!
!		do j=1,m-1
!
!		  st_int%f(:,i_int-j) = st%f(:,sys%nmode1+1+i) &
!						  + (st%f(:,sys%nmode1+i) - st%f(:,sys%nmode1+i+1))*dble(j)/dble(m)
!		  st_int%fo(:,i_int-j) = st%fo(:,sys%nmode1+1+i) &
!						  + (st%fo(:,sys%nmode1+i) - st%fo(:,sys%nmode1+i+1))*dble(j)/dble(m)
!		  st_int%h(:,i_int-j) = st%h(:,sys%nmode1+1+i) &
!						  + (st%h(:,sys%nmode1+i) - st%h(:,sys%nmode1+i+1))*dble(j)/dble(m)
!		  st_int%ho(:,i_int-j) = st%ho(:,sys%nmode1+1+i) &
!						  + (st%ho(:,sys%nmode1+i) - st%ho(:,sys%nmode1+i+1))*dble(j)/dble(m)
!
!		end do
!
!	 end do
!
!	 !-- last interpollation
!	 i = size(st%f,2) - sys%nmode1 -1-1
!	 i_int = sys%nmode1+(m+1)/2+m*i
!	 do j=1,size(st_int%f,2)-i_int
!
!		st_int%f(:,i_int+j) = st_int%f(:,i_int) + (st_int%f(:,i_int) - st_int%f(:,i_int-1))*j
!		st_int%fo(:,i_int+j) = st_int%fo(:,i_int) + (st_int%fo(:,i_int) - st_int%fo(:,i_int+1))*j
!		st_int%h(:,i_int+j) = st_int%h(:,i_int) + (st_int%h(:,i_int) - st_int%h(:,i_int-1))*j
!		st_int%ho(:,i_int+j) = st_int%ho(:,i_int) + (st_int%ho(:,i_int) - st_int%ho(:,i_int-1))*j
!
!	 end do
!
!	 st_int%f(:,sys%nmode1+1:) = st_int%f(:,sys%nmode1+1:)/sqrt(dble(m))
!	 st_int%fo(:,sys%nmode1+1:) = st_int%fo(:,sys%nmode1+1:)/sqrt(dble(m))
!	 st_int%h(:,sys%nmode1+1:) = st_int%h(:,sys%nmode1+1:)/sqrt(dble(m))
!	 st_int%ho(:,sys%nmode1+1:) = st_int%ho(:,sys%nmode1+1:)/sqrt(dble(m))
!
!	 CALL update_sums(sys,st_int)
!	 CALL normalise(st_int)
!
!	 print*,"ARRAYS INTERPOLATED"
!
!  END SUBROUTINE
!










