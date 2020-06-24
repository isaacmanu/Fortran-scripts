program bessel
implicit none
integer :: i
real :: v, xinit, xfinal, h, xchange, step, f
real, allocatable, dimension(:) :: x, y1, y2 !creates three allocatable 1d arrays

v=0.0     !initial conditions
xinit=0
xfinal=25
h=0.001
xchange=0.05
step= 25000/500

allocate(x(int(25000))) !creates array dimensions
allocate(y2(int(25000)))
allocate(y1(int(25000)))

y1(1)=1.150443474
y2(1)=0.0
x(1)=0.0

y1(2)=(h*y2(1))+y1(1)!euler step
y2(2)=0.0
x(2)=h

open(1, file='bessel.txt')

do i=3,25001!euler step did i=1,2
    x(i)= h + x(i-1)
    if (x(i-2)==0) then  !solution is not singular at x=0
    y1(i)=y1(i-1) + 0.5*h*(3.0*y2(i-1)-y2(i-2))
    y2(i)=0.0
    else
    y1(i)= y1(i-1) + 0.5*h*(3.0*y2(i-1)-y2(i-2))
    y2(i)= y2(i-1) + 0.5*h*(3.0*f(x(i-1), v, y1(i-1), y2(i-1)) -f(x(i-2), v, y1(i-2),y2(i-2)))
    if (mod(real(i-1), real(50)) == 0) then 
    print'(f15.7, f15.7)', x(i), y1(i) !prints multiples of 0.05
    write(1, fmt='(f10.7, a, f10.7)') x(i), ',', y1(i)!writes to text file
end if
end if
end do
close(1)
print'(e15.7)', y1(25001)
end program  
function f(x,v,y1,y2)
implicit none
real :: f, x, v, y1, y2

f=((v**2)/(x**2)-1.0)*y1 -y2/x !ODE from question

end function
