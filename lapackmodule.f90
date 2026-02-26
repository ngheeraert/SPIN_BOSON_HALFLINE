! ==============================================================================
!  MODULE: lapackmodule.f90
!
!  PURPOSE & CONTEXT
!    High-performance linear algebra interface for the waveguide QED simulation.
!    Provides Fortran wrappers around standard BLAS and LAPACK routines to handle 
!    the dense, complex tensor operations required by the TDVP equations of motion.
!    Associated paper: N. Gheeraert et al., New J. Phys. 19 (2017) 023036.
!
!  PHYSICS & NUMERICS CONTEXT
!    The multi-polaron ansatz yields a non-orthogonal basis. The resulting metric 
!    tensors (overlap matrices) and explicit derivative matrices must be inverted s
!    or solved continuously during the RK4 time integration. Because these matrices 
!    can become ill-conditioned as the basis grows, robust LAPACK solvers with 
!    pivoting (LU, Bunch-Kaufman) are essential to maintain numerical stability.
!
!  CORE RESPONSIBILITIES
!    1. Inversions      : Cholesky (`ZPOTRF`), LU (`ZGETRF`), and Bunch-Kaufman (`ZHETRF`).
!    2. Linear Solvers  : Direct linear system solvers (`DGESV`, `ZGESV`).
!    3. Operator Overloads: Maps custom operators (`.matadd.`, `.matprod.`) to 
!                           highly optimized Level-3 BLAS (`ZGEMM`, `DGEMM`).
! ==============================================================================
MODULE lapackmodule

USE typedefs, only : cx => c_type, rl=> r_type
USE consts

implicit none

  ! ----------------------------------------------------------------------------
  ! OPERATOR INTERFACE
  !   These overloads allow the physics code in `system.f90` to remain readable 
  !   (e.g., `A .matprod. B`) while guaranteeing the underlying execution utilizes 
  !   cache-optimized BLAS routines rather than unoptimized intrinsic operations.
  ! ----------------------------------------------------------------------------
  INTERFACE OPERATOR(.matprod.)
    module procedure matmultiply_c_rank12, matmultiply_c_rank21, matmultiply_c_rank22,&
           matmultiply_r_rank12, matmultiply_r_rank21, matmultiply_r_rank22
  END INTERFACE OPERATOR(.matprod.)
  
  
CONTAINS

  !> -------------------------------------------------------------------------
  !> SUBROUTINE: InvertHPD
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Inverts a Hermitian Positive Definite (HPD) matrix using Cholesky 
  !>   factorization (`ZPOTRF` / `ZPOTRI`). Highly efficient, but requires 
  !>   strict positive-definiteness. Often used for perfectly conditioned 
  !>   overlap metrics. Explicitly reconstructs the lower triangular half.
  !> Arguments:
  !>   - A : The HPD matrix to be inverted (modified in-place).
  !>
  SUBROUTINE InvertHPD(A)

    COMPLEX(8), intent(in out)                      ::  A(:,:)
	 INTEGER                	   	            ::  INFO, LDA,N,i,j

	 LDA = size(A,1)
	 N = size(A,2)
	 info=0
	 
	 CALL ZPOTRF('U',N,A,LDA,INFO)  !-- Performs the Choelesky factorisation

	 if (info==0) then

		CALL ZPOTRI('U',N,A,LDA,INFO) !-- CAREFUL: returns only triangular part
		do i=1,N
		  do j=1,N
			 if (i>j) then
				a(i,j) = a(j,i) 
			 end if
		  end do
		end do

		if (info /= 0) then
		  print*, "Failure in the inversion step, ZPOTRI"
		end if
	 else 
	 	print*, "Failure in ZPOTRF, choelesky factorisation"
	 end if

  END SUBROUTINE InvertHPD

  !> -------------------------------------------------------------------------
  !> SUBROUTINE: InvertGeneral
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Inverts a general dense complex matrix via LU factorization with partial 
  !>   pivoting (`ZGETRF` / `ZGETRI`). Used as a robust fallback when matrices 
  !>   lose symmetry or positive-definiteness during highly non-linear dynamics.
  !> Arguments:
  !>   - A    : The general complex matrix to be inverted (modified in-place).
  !>   - info : Output error flag (0 = success).
  !>
  SUBROUTINE InvertGeneral(A,info)
    COMPLEX(cx), intent(in out)                  ::  A(:,:)
    integer, intent(out)								 ::  info
	 INTEGER                	   		          ::  LDA,N,LWORK,i,j
	 INTEGER, dimension(size(A,1))		          ::  IPIV 	
	 COMPLEX(cx), allocatable 				          ::  WORK(:)
	 
	 LDA = size(A,1)
	 N = size(A,1)
	 info=0
	 LWORK = N
	 allocate(WORK(LWORK))
	 
	 CALL ZGETRF(N,N,A,LDA,IPIV,INFO)  

	 if (info==0) then

		CALL ZGETRI(N,A,LDA,IPIV,WORK,LWORK,INFO) 

		if (info /= 0) then
		  print*, "Failure in the inversion step, ZGETRI"
		  print*,"info=", info
		end if

	 else 
		print*, "Failure in ZGETRF"
		print*,"info=", info
	 end if

  END SUBROUTINE InvertGeneral

  !> -------------------------------------------------------------------------
  !> SUBROUTINE: InvertH
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Inverts a Hermitian Indefinite matrix using Bunch-Kaufman diagonal 
  !>   pivoting (`ZHETRF` / `ZHETRI`). Safe for Hermitian overlap matrices 
  !>   that may have developed near-zero eigenvalues due to basis over-completeness.
  !> Arguments:
  !>   - A    : The complex Hermitian matrix (modified in-place).
  !>   - info : Output error flag (0 = success).
  !>
  SUBROUTINE InvertH(A,info)
    COMPLEX(cx), intent(in out)                  ::  A(:,:)
    integer, intent(out)								 ::  info
	 INTEGER                	   		          ::  LDA,N,LWORK,i,j
	 INTEGER, dimension(size(A,1))		          ::  IPIV 	
	 COMPLEX(cx), allocatable 				          ::  WORK(:)
	 
	 LDA = size(A,1)
	 N = size(A,2)
	 info=0
	 LWORK = N
	 allocate(WORK(LWORK))
	 
	 CALL ZHETRF('U',N,A,LDA,IPIV,WORK,LWORK,INFO)  !-- Performs the Bunch-Kaufman factorisation

	 if (info==0) then

		CALL ZHETRI('U',N,A,LDA,IPIV,WORK,INFO) !-- CAREFUL: returns only triangular part
		do i=1,N
		  do j=1,N
			 if (i>j) then
				a(i,j) = conjg(a(j,i))
			 end if
		  end do
		end do

		if (info /= 0) then
		  print*, "Failure in the inversion step, ZHETRI"
		  print*,"info=", info
		end if

	 else 
		print*, "Failure in ZHETRF, Bunch-Kaufman factorisation"
		print*,"info=", info
	 end if

  END SUBROUTINE InvertH

  !> -------------------------------------------------------------------------
  !> SUBROUTINE: SolveEq_r
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Solves the real linear system $A \cdot x = B$ using standard LU 
  !>   decomposition (`DGESV`). Overwrites $B$ with the solution vector.
  !>
  SUBROUTINE SolveEq_r(A,B)
    REAL(rl), intent(in out)               ::  A(:,:)
	 REAL(rl), intent(in out)					 ::  B(:)
	 INTEGER                	   		    ::  INFO,LDA,LDB,N,NRHS
	 INTEGER, dimension(size(A,1))		    ::  IPIV   !-- pivot indices	
	 
	 NRHS = 1  						!-- number of right hand sides
	 N = size(A,1)					!-- the number of linear equations
	 LDA = size(A,1)				!-- the leading dimension of A 
	 LDB = size(B,1)				!-- the leading dimension of B
	 info=0							!-- 0 is successful

	 CALL DGESV(N,NRHS,A,LDA,IPIV,B,LDB,INFO)  !-- Solve by performing the Bunch-Kaufman factorisation

	 if (info /= 0) then
	 	print*, "Failure in DGESV - solving real system of equations"
	 	print*,"INFO = ",info
	 end if

  END SUBROUTINE SolveEq_r

  !> -------------------------------------------------------------------------
  !> SUBROUTINE: SolveEq_c
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Solves the complex linear system $A \cdot x = B$ using LU decomposition 
  !>   (`ZGESV`). Used predominantly by `inverse.f90` to resolve the flattened 
  !>   variational $\kappa$ tensor system.
  !>
    SUBROUTINE SolveEq_c(A,B)

    COMPLEX(cx), intent(in out)               ::  A(:,:)
	 COMPLEX(cx), intent(in out)					 ::  B(:)
	 INTEGER                	   		    	 ::  INFO,LDA,LDB,N,NRHS
	 INTEGER, dimension(size(A,1))		    	 ::  IPIV   !-- pivot indices

	 NRHS = 1  						!-- number of right hand sides
	 N = size(A,1)					!-- the number of linear equations
	 LDA = size(A,1)				!-- the leading dimension of A multiply_c
	 LDB = size(B,1)				!-- the leading dimension of B
	 info=0							!-- 0 is successful

	 CALL ZGESV(N,NRHS,A,LDA,IPIV,B,LDB,INFO)  !-- Solve by performing the Bunch-Kaufman factorisation

	 if (info /= 0) then
	 	print*, "Failure in DGESV - solving real system of equations"
	 	print*,"INFO = ",info
	 end if

 END SUBROUTINE

  !> -------------------------------------------------------------------------
  !> FUNCTIONS: matmultiply_c_rank** / matmultiply_r_rank**
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   A suite of operator overloads mapping `.matprod.` to Level-3 BLAS (`ZGEMM` 
  !>   and `DGEMM`). 
  !>   - rank22: Matrix * Matrix.
  !>   - rank21: Matrix * Column Vector.
  !>   - rank12: Row Vector * Matrix (transposes the 1D input array automatically).
  !>
  ! -- Complex Matrix Multiplication --
  FUNCTION matmultiply_c_rank22(A,B) RESULT(C)

    complex(cx), dimension(:,:), intent(in)            ::  A
    complex(cx), dimension(:,:), intent(in)            ::  B
    complex(cx), dimension(size(A,1),size(B,2))        ::  C
	 integer													       ::  M,N,K,LDA,LDB,LDC

	 M = size(A,1)    !-- number of rows of A
	 N = size(B,2)	   !-- number of columns of B
	 K = size(A,2)    !-- number of cols of A/rows of B
	 LDA = size(A,1)  !-- leading dimension of A
	 LDB = size(B,1)  !-- leading dimension of B
	 LDC = size(C,1)	!-- leading dimension of C


    CALL zgemm ('N','N',M,N,K,one_r,A,LDA,B,LDB,zero_r,C,LDC)

  END FUNCTION matmultiply_c_rank22
  FUNCTION matmultiply_c_rank21(A,B) RESULT(C)

    complex(cx), dimension(:,:), intent(in)            :: A
    complex(cx), dimension(:), intent(in)              :: B
    complex(cx), dimension(size(A,1))                  :: C
	 integer													       :: M,N,K,LDA,LDB,LDC

	 M = size(A,1)    !-- number of rows of A
	 N = 1			   !-  number of columns of B
	 K = size(A,2)    !-- number of cols of A/rows of B
	 LDA = size(A,1)  !-- leading dimension of A
	 LDB = size(B,1)  !-- leading dimension of B
	 LDC = size(C,1)	!-- leading dimension of C

    CALL zgemm ('N','N',M,N,K,one_r,A,LDA,B,LDB,zero_r,C,LDC)

  END FUNCTION matmultiply_c_rank21
  FUNCTION matmultiply_c_rank12(A,B) RESULT(C)

    complex(cx), dimension(:), intent(in)        :: A
    complex(cx), dimension(:,:), intent(in)      :: B
    complex(cx), dimension(size(B,2))            :: C
    complex(cx), dimension(1,size(B,2))			 :: Cbis
	 integer												    :: M,N,K,LDA,LDB,LDCbis

	 M = 1		            !-- number of rows of A
	 N = size(B,2)		      !-- number of columns of B
	 K = size(A,1)          !-- number of cols of A/rows of B
	 LDA = size(A,1)        !-- leading dimension of A
	 LDB = size(B,1)        !-- leading dimension of B
	 LDCbis = size(Cbis,1)	!-- leading dimension of C

    CALL zgemm ('T','N',M,N,K,one_r,A,LDA,B,LDB,zero_r,C,LDCbis)
    C(:)=Cbis(1,:)

  END FUNCTION matmultiply_c_rank12

  !-- Real Matrix Multiplication
  FUNCTION matmultiply_r_rank22(A,B) RESULT(C)

    real(rl), dimension(:,:), intent(in)            :: A
    real(rl), dimension(:,:), intent(in)            :: B
    real(rl), dimension(size(A,1),size(B,2))        :: C
	 integer													    :: M,N,K,LDA,LDB,LDC

	 M = size(A,1)    !-- number of rows of A
	 N = size(B,2)	   !-- number of columns of B
	 K = size(A,2)    !-- number of cols of A/rows of B
	 LDA = size(A,1)  !-- leading dimension of A
	 LDB = size(B,1)  !-- leading dimension of B
	 LDC = size(C,1)	!-- leading dimension of C

    CALL dgemm ('N','N',M,N,K,one_r,A,LDA,B,LDB,zero_r,C,LDC)

  END FUNCTION matmultiply_r_rank22
  FUNCTION matmultiply_r_rank21(A,B) RESULT(C)

    real(rl), dimension(:,:), intent(in)            :: A
    real(rl), dimension(:), intent(in)              :: B
    real(rl), dimension(size(A,1))                :: C
	 integer													 :: M,N,K,LDA,LDB,LDC

	 M = size(A,1)    !-- number of rows of A
	 N = 1			   !-  number of columns of B
	 K = size(A,2)    !-- number of cols of A/rows of B
	 LDA = size(A,1)  !-- leading dimension of A
	 LDB = size(B,1)  !-- leading dimension of B
	 LDC = size(C,1)	!-- leading dimension of C

    CALL dgemm ('N','N',M,N,K,one_r,A,LDA,B,LDB,zero_r,C,LDC)

  END FUNCTION matmultiply_r_rank21
  FUNCTION matmultiply_r_rank12(A,B) RESULT(C)

    real(rl), dimension(:), intent(in)        :: A
    real(rl), dimension(:,:), intent(in)      :: B
    real(rl), dimension(size(B,2))            :: C
    real(rl), dimension(1,size(B,2))			 :: Cbis
	 integer												 :: M,N,K,LDA,LDB,LDCbis

	 M = 1		      !-- number of rows of A
	 N = size(B,2)		!-- number of columns of B
	 K = size(A,1)    !-- number of cols of A/rows of B
	 LDA = size(A,1)  !-- leading dimension of A
	 LDB = size(B,1)  !-- leading dimension of B
	 LDCbis = size(Cbis,1)	!-- leading dimension of C

    CALL dgemm ('T','N',M,N,K,one_r,A,LDA,B,LDB,zero_r,Cbis,LDCbis)
    C(:)=Cbis(1,:)

  END FUNCTION matmultiply_r_rank12


END MODULE lapackmodule
