program correlation
implicit none
complex*16 jj,rho1,nu_p,nu_q,J_x1,arg_int,gprime,rho2,J_x2
complex*16 intJ_x, intJ_x_2, rhoint1, rhoint2, rhoint3, rhoint3_2
real*8 pi,delta_p,delta_q,beta,n_p, n_q
real*8 p,q,lambda,g,mu,w_p,w_q,A,Gam,Gam2
real*8 x,dx,x0,alpha,xi,dlambda,x1,x2 
real*8, allocatable :: quad_w(:), quad_x(:)
integer i,j,l,quad_order,io

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

g=1d0; beta=0.1d0; mu=1.22713
alpha = 2./g

open(22,file='./data/xi_lambda_77.dat') !asymptotics

open(10, file='quadrature_w.txt')
open(11, file='quadrature_x.txt')

do i=1,quad_order
read(10,*) quad_w(i)
read(11,*) quad_x(i) 
enddo 

!print*, quad_w, quad_x

close(10); close(11)

x1 = 10.; x2 = 20.
lambda = -10.; dlambda = 0.1d0

do while(lambda.le.10.)

rho1 = 0.; rho2 = 0.
J_x1 = 0.; J_x2 = 0.
rhoint1 = 0.; rhoint2 = 0.
rhoint3_2 = 0.


do i=1, quad_order
q = quad_x(i); w_q = quad_w(i) !!!VARIABLE HERE IS Q!!!

delta_q = pi/2. -atan(lambda -2.*q/g)
n_q = 1./(1. +exp(beta*(q**2. -mu)))
nu_q = (1./(2.*pi*jj)) * log(1. +n_q*(exp(2.*jj*delta_q) -1.))


gprime = 2.*exp(beta*(q**2. -mu))*(4.* q**3.*beta -2.*g*(jj +2.*q**2.*beta*lambda) +g**2.*q*beta*(1.+lambda**2.)) / &
        ((exp(beta*(q**2.-mu))-exp(2.*jj*atan(2*q/g-lambda)))*(4.*q**2. +g*(g-4.*q*lambda+g*lambda**2.)))

!print*, gprime

intJ_x = 0.; rhoint3 = 0.

do j=1, quad_order
p = quad_x(j); w_p = quad_w(j) !!!VARIABLE HERE IS P!!!

delta_p = pi/2. -atan(lambda -2.*p/g)
n_p = 1./(1. +exp(beta*(p**2. -mu)))
nu_p = (1./(2.*pi*jj)) * log(1. +n_p*(exp(2.*jj*delta_p) -1.))

if(p.eq.q) then
        arg_int = (2.*exp(beta*mu)*g*(exp(beta*mu)+exp(beta* p**2.)*(1. +p*beta*(2.*p-g*(lambda-jj)))))/ &
               ((exp(beta* p**2.)+exp(beta*mu))*pi *(-2.*p +g*(lambda-jj))* (exp(beta* p**2.)*(-2.*p +g*(lambda-jj)) &
              +exp(beta*mu)*(-2.*p +g*(lambda+jj)))) 
else 
        arg_int = (nu_p - nu_q)/(p-q)
endif

rhoint3 = rhoint3 + arg_int**2. *w_p !DOUBLE INTEGRAL

intJ_x = intJ_x + arg_int*w_p
enddo

rhoint1 = rhoint1 + nu_q*w_q
rhoint2 = rhoint2 + nu_q*gprime*w_q

rhoint3_2 = rhoint3_2 + rhoint3*w_q !DOUBLE INTEGRAL

J_x1 = J_x1 + (1./pi) * n_q*(sin(delta_q)**2.)*exp(-jj*x1*q -2.*intJ_x)*w_q
J_x2 = J_x2 + (1./pi) * n_q*(sin(delta_q)**2.)*exp(-jj*x2*q -2.*intJ_x)*w_q
enddo

rho1 = J_x1*exp(jj*x1*rhoint1 -rhoint2 -rhoint3_2/2.)
rho2 = J_x2*exp(jj*x2*rhoint1 -rhoint2 -rhoint3_2/2.)


xi = log(abs(rho2/rho1))/(x1-x2)
A = abs(rho1)*exp(x1*xi)
Gam = acos(real(rho2)*exp(x2*xi)/A)/x2
Gam2 = asin(-aimag(rho2)*exp(x2*xi)/A)/x2



if(lambda .gt. 1.9999 .and. lambda .lt. 2.00001) print*, xi

write(22,*) lambda, xi !xi !x, xi !asymptotics

lambda = lambda+dlambda
enddo

close(22)

end program
