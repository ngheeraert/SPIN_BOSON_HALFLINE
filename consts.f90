module consts

  use typedefs 
  implicit none

  real(r_type),   parameter :: c_light = 1_r_type
  real(r_type),   parameter :: pi=(4.0_r_type)*atan(1.0_r_type)
  real(r_type),   parameter :: one_r=1.0_r_type
  real(r_type),   parameter :: zero_r=0.0_r_type

  real(q_type),   parameter :: pi_q=(4.0_q_type)*atan(1.0_q_type)
  real(q_type),   parameter :: one_q=1.0_q_type
  real(q_type),   parameter :: zero_q=0.0_q_type

  complex(c_type),parameter :: im=(0.0d0,1.0d0)
  complex(c_type),parameter :: zero=(0.0d0,0.0d0)
  complex(c_type),parameter :: one=(1.0d0,0.0d0),two=(2.0d0,0.0d0)

contains

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

end module consts
