program correlation
implicit none
complex*16 jj,rho,n_k,norm, norm_ex
complex*16, allocatable :: rho_quad(:)
real*8 pi,x,k,dk,w,dx
real*8, allocatable :: quad_w(:), quad_x(:)
integer i,j,quad_order,io,rho_order

open(1,file='./data/rho_x_effective_quadrature_a1_b1.dat')
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


open(20,file='./data/n_k_effective_a1_b1.dat')
open(22,file='./data/rho_x_effective_quadrature_a1_b1.dat')

do i=1,rho_order
read(22,*) quad_x(i), quad_w(i), rho_quad(i)
enddo
close(22)

call ex(norm_ex)

call eff(norm)

k = -10.; dk = 1d-2
do while (k.le.10.)

n_k = 0.
 do i=1,rho_order
 x = quad_x(i)
 if(x.ge.0) then
 w = quad_w(i)
 rho = rho_quad(i)
 n_k = n_k + exp(jj*k*x)*rho*w
 endif
enddo



n_k = n_k/pi

print*, k, n_k
write(20,*) k, real(n_k*norm_ex/norm)!, aimag(n_k), abs(n_k)
k=k+dk
enddo


end program





subroutine eff(norm)
implicit none
complex*16 jj,rho,n_k,norm
complex*16, allocatable :: rho_quad(:)
real*8 pi,x,k,dk,w,dx
real*8, allocatable :: quad_w(:), quad_x(:)
integer i,j,quad_order,io,rho_order

open(1,file='./data/rho_x_effective_quadrature_a1_b1.dat')
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

open(22,file='./data/rho_x_effective_quadrature_a1_b1.dat')

do i=1,rho_order
read(22,*) quad_x(i), quad_w(i), rho_quad(i)
enddo
close(22)


n_k = 0.
 do i=1,rho_order
 x = quad_x(i)
 if(x.ge.0) then
 w = quad_w(i)
 rho = rho_quad(i)
 n_k = n_k + exp(jj*k*x)*rho*w
 endif
enddo

n_k = n_k/pi
norm = n_k
end subroutine



subroutine ex(norm_ex)
implicit none
complex*16 jj,rho,n_k, norm_ex
complex*16, allocatable :: rho_quad(:)
real*8 pi,x,k,dk,w,dx
real*8, allocatable :: quad_w(:), quad_x(:)
integer i,j,quad_order,io,rho_order

open(1,file='./data/rho_x_exact_quadrature_a1_b1.dat')
rho_order=0
do
read(1,*,iostat=io)
if(io.ne.0) exit
rho_order = rho_order + 1
enddo
close(1)

allocate(rho_quad(rho_order))
allocate(quad_w(rho_order)); allocate(quad_x(rho_order))

pi=acos(-1.d0); jj = cmplx(0.d0,1.d0)


open(22,file='./data/rho_x_exact_quadrature_a1_b1.dat')

do i=1,rho_order
read(22,*) quad_x(i), quad_w(i), rho_quad(i)
enddo
close(22)



k = 0.
n_k = 0.
 do i=1,rho_order
 x = quad_x(i) 
 w = quad_w(i) 
 rho = rho_quad(i)
 !if(x.ge.0.) n_k = n_k + exp(jj*k*x)*rho*w
 n_k = n_k + exp(jj*k*x)*rho*w
enddo 

n_k = n_k/pi

norm_ex = n_k


end subroutine
