program correlation
implicit none
complex*16 jj,rho,n_k
complex*16, allocatable :: rho_quad(:)
real*8 pi,x,k,dk,w,dx
real*8, allocatable :: quad_w(:), quad_x(:)
integer i,j,quad_order,io,rho_order

open(1,file='./data/rho_x_exact_quadrature_a1_b2.dat')
rho_order=0
do
read(1,*,iostat=io)
if(io.ne.0) exit
rho_order = rho_order + 1
enddo
close(1)

print*, rho_order

allocate(rho_quad(rho_order))
allocate(quad_w(rho_order)); allocate(quad_x(rho_order))

pi=acos(-1.d0); jj = cmplx(0.d0,1.d0)


open(20,file='./data/n_k_exact_a1_b2.dat')
open(22,file='./data/rho_x_exact_quadrature_a1_b2.dat')

do i=1,rho_order
read(22,*) quad_x(i), quad_w(i), rho_quad(i)
enddo
close(22)



k = -10.; dk = 1d-2
do while (k.le.10.)

n_k = 0.
 do i=1,rho_order
 x = quad_x(i) 
 w = quad_w(i) 
 rho = rho_quad(i)
 !if(x.ge.0.) n_k = n_k + exp(jj*k*x)*rho*w
 n_k = n_k + exp(jj*k*x)*rho*w
enddo 

n_k = n_k/pi
print*, k, n_k
write(20,*) k, real(n_k) !, aimag(n_k), abs(n_k)
k=k+dk
enddo

!close(20);close(1);close(10)
!close(2)

end program
