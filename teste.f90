program testando
real*8 integral, func
integer, parameter :: n=100
real*8 w(n), y(n)
integer i

open(10,file='quadrature_x.txt')
open(11,file='quadrature_w.txt')

read(10,*) y
read(11,*) w

        !integral=0.; dx=1d-6
        !x=-1.
        !do while (x<=9.)
        !integral = integral + exp(-x)*dx
        !x = x+dx
        !enddo

        integral=0.
        do i=1,n
        func = exp(-y(i))*cos(y(i)**2.)/log(sqrt(y(i)-1)/(y(i)+1.))
        integral = integral + w(i)*func 
        enddo

        !print*, integral


end program
