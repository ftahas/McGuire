program correlation
implicit none
complex*16 jj,rho,nu_p,nu_q,J_x,arg_int,gprime,rho_l
complex*16 intJ_x,rhoint1, rhoint2, rhoint3, rhoint3_2
real*8 pi,delta_p,delta_q,beta,n_p, n_q
real*8 p,q,lambda,g,mu,w_p,w_q
real*8 x,dx,ec,k,fd,w,alpha, w_l
real*8, allocatable :: quad_w(:), quad_x(:)
integer i,j,l,m,quad_order,io

quad_order=0
open(1,file='quadrature_x.txt')
do
read(1,*,iostat=io)
if(io.ne.0) exit
quad_order = quad_order + 1
enddo
close(1)

allocate(quad_w(quad_order)); allocate(quad_x(quad_order))

pi=acos(-1.d0); jj = cmplx(0.d0,1.d0)


alpha = 1d-8
g=2./alpha
beta=1. 
mu=1.

open(20,file='./data/rho_x_77_a0_b1.dat')


open(10, file='quadrature_w.txt')
open(11, file='quadrature_x.txt')

do i=1,quad_order
read(10,*) quad_w(i)
read(11,*) quad_x(i) 
enddo 

!print*, quad_w, quad_x

close(10); close(11)

x = 0.; dx = 1d-1
do while(x.le.5.)
if (x .lt. 1.d0) dx = 1d-1
if (x .ge. 1.d0) dx = 1d-1

rho = 0.

do m=1, quad_order
lambda = quad_x(m); w_l = quad_w(m)

J_x = 0.; rhoint1 = 0.; rhoint2 = 0.
rhoint3_2 = 0.

do i=1, quad_order
q = quad_x(i); w_q = quad_w(i) !!!VARIABLE HERE IS Q!!!

delta_q = pi/2. -atan(lambda -2.*q/g)
n_q = 1./(1. +exp(beta*(q**2. -mu)))
nu_q = (1./(2.*pi*jj)) * log(1. +n_q*(exp(2.*jj*delta_q) -1.))


gprime = 2.*exp(beta*(q**2. -mu))*(4.* q**3.*beta -2.*g*(jj +2.*q**2.*beta*lambda) +g**2.*q*beta*(1.+lambda**2.)) / &
        ((exp(beta*(q**2.-mu))-exp(2.*jj*atan(2.*q/g-lambda)))*(4.*q**2. +g*(g-4.*q*lambda+g*lambda**2.)))

!print*, gprime

intJ_x = 0.; rhoint3 = 0.
do j=1, quad_order
p = quad_x(j); w_p = quad_w(j) !!!VARIABLE HERE IS P!!!

delta_p = pi/2. -atan(lambda -2.*p/g)
n_p = 1./(1. +exp(beta*(p**2. -mu)))
nu_p = (1./(2.*pi*jj)) * log(1. +n_p*(exp(2.*jj*delta_p) -1.))

if(abs(p-q).le.1d-3) then
        arg_int = (2.*exp(beta*mu)*g*(exp(beta*mu)+exp(beta* p**2.)*(1. +p*beta*(2.*p-g*(lambda-jj)))))/ &
               ((exp(beta* p**2.)+exp(beta*mu))*pi *(-2.*p +g*(lambda-jj))* (exp(beta* p**2.)*(-2.*p +g*(lambda-jj)) &
              +exp(beta*mu)*(-2.*p +g*(lambda+jj)))) 
else 
        arg_int = (nu_p - nu_q)/(p-q)
endif

rhoint3 = rhoint3 + arg_int**2. *w_p !DOUBLE INTEGRAL

intJ_x = intJ_x + 2.*arg_int*w_p
enddo

rhoint1 = rhoint1 + nu_q*w_q
rhoint2 = rhoint2 + nu_q*gprime*w_q

rhoint3_2 = rhoint3_2 + rhoint3*w_q !DOUBLE INTEGRAL

J_x = J_x + (1./pi) * n_q*(sin(delta_q)**2.)*exp(-jj*x*q +2.*intJ_x)*w_q
enddo

!print*, rho

rho_l = J_x*exp(jj*x*rhoint1 -rhoint2 -rhoint3_2/2.)

!print*, rhoint1, rhoint2, rhoint3

!print*, rho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!calc of ec(lambda)!!!!!!!!!!!!!!!!!!!!!!!!!!
ec = 0.
do l=1,quad_order
k = quad_x(l); w = quad_w(l)
fd = 1./(1.+exp(beta*(k**2. -mu)))
ec = ec -2.*(k/pi)*fd*(pi/2. - atan(lambda -2.*k/g))*w
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

rho = rho + exp(-beta*ec)*rho_l*w_l
enddo

print*, x, rho

write(20,*) x, real(rho)
x = x+dx
enddo


close(20)
end program
