program rk4
implicit none
real ::  h, u, unum, f1, f2, f3, f4, ti, ti1, f, uanalytic, relerror
integer :: N, i, counter

N=1000
u=1.0+0.150443474
ti=0
counter=0


open(1, file="rk2.txt", status="replace", position="append")!opens text file

do i=0,N-1!1 less than N as loop starts at 0
ti1=ti+(1./N)
h=ti1-ti
f1=h*f(u,ti)                 !runge-kutta 4th order method
f2=h*f(u+(f1/2.), ti+(h/2.))
f3=h*f(u+(f2/2.), ti+(h/2.))
f4=h*f(u+f3, ti+h)

unum=u+ (1./6.)*(f1+2*f2+2*f3+f4)
uanalytic=1.0 + 0.150443474*exp(-ti**2)
if (counter==N/10.) then
    counter=0!resets counter so next 0.1 step can be calculated
    write(1, '(3e15.7)') ti, unum, uanalytic
    write(6, '(3e15.7)') ti, unum, uanalytic
end if

ti=ti1
u=unum
counter=counter+1
end do

relerror=abs((uanalytic-unum)/uanalytic)

print'(e15.7)',relerror

print'(e15.7)', unum

end program


real function f(u,t)
implicit none
real :: u, t
f=-2*u*t +2*t
end function f
