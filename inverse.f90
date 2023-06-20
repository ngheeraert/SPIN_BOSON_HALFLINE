MODULE INVERSE
  USE consts
  USE lapackmodule
  USE typedefs, only : cx => c_type, rl=> r_type
  IMPLICIT none


  CONTAINS

  subroutine directInverse(n,A,sol,rhs)
    complex(cx),intent(in)     :: A(n,n,n)
    complex(cx),intent(in)     :: rhs(n**2)
    complex(cx),intent(in out) :: sol(n**2)
    complex(cx)           :: A2D(n**2, n**2)
    integer               :: n,nn, i,j,k  !-- n=block size, nn=array size


    !-- Build Matrix
    A2D=0._cx
    do i=1, n
      do j=1, n
        do k=1, n
          A2D((j-1)*n+i,(j-1)*n+k)=A(i,j,k)
        end do
        A2D((j-1)*n+i,(i-1)*n+j)=1.0_cx
      end do
    end do
    ! -- resolve
    sol=rhs
    Call solveEq_c(A2D, sol)
  end subroutine directInverse

  FUNCTION kroneckerDelta(a,b)

 		integer,intent(in)   ::  a,b
 		real(rl)				   ::  kroneckerDelta

 		if (a==b) then
 		  kroneckerDelta = 1
 		else
 		  kroneckerDelta = 0
 		end if
  END FUNCTION


END MODULE INVERSE
