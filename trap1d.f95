program trap1d
implicit none

integer :: i
integer, parameter :: N=40
real :: s, a, b, h, f, func1d, I_num, x, i_an, abs_err

s=0.150443474 !student number
a=-1.0-s      ! limits
b=1.0+s
h=(b-a)/N     !step size

I_num=0

do i=0,N
   x=a+i*h
   f=func1d(x)
   if (i/=0 .and. i/=N) then
      I_num=I_num+ f
   else
   I_num=I_num+ 0.5*f
   end if
end do

I_num=I_num*h

I_an=1.367007239  !analytical I

abs_err=abs(I_an-I_num)

print'(e15.7)', I_num

end program trap1d

real function func1d(x)

func1d=exp(-abs(x))

end function func1d

!For N=10, I=0.1373033E+01, absolute error=0.6025553E-02

!For N=20, I=0.1368514E+01, absolute error=0.1507163E-02

!For N=40, I=0.1367384E+01, absolute error=0.3769398E-03



