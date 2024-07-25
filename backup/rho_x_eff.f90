program correlation
implicit none
complex*16 jj,rho,nu_p,nu_q,J_x,arg_int,gprime,norm,dnu_q
complex*16 intJ_x,rhoint1,rhoint2,rhoint3,gfun,dnu_p,arg_int_J
real*8 pi,delta_p,delta_q,beta,n_p, n_q
real*8 p,q,lambda,g,mu,w_p,w_q,w_l
real*8 x,dx,ec,k,fd,w,x0,xf,alpha,dbeta
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

print*, 'quadrature = ', quad_order

allocate(quad_w(quad_order)); allocate(quad_x(quad_order))

pi=acos(-1.d0); jj = cmplx(0.d0,1.d0)

mu=1.d0

beta = 1.d0
g = 2.d0

alpha = 2.d0/g

open(20,file='./data/rho_x_eff.dat')

open(10, file='quadrature_w.txt')
open(11, file='quadrature_x.txt')

do i=1,quad_order
read(10,*) quad_w(i)
read(11,*) quad_x(i) 
enddo 

!print*, quad_w, quad_x

close(10); close(11)

!write(*,*) 'Enter x0 and xf'
!read(*,*) x0, xf

!x = x0
!do while(x.le.xf)
!if (x .lt. 1.d0) dx = 1d-1
!if (x .ge. 1.d0) dx = 1d-1

x = 0.

beta = 0.5; dbeta = 0.1
do while (beta.le.2.)
rho = 0. 
do l=1, quad_order
lambda = quad_x(l)
w_l = quad_w(l)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!GOOD NU!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (abs(lambda).ge.alpha*sqrt(mu)) then

J_x = 0.; rhoint1 = 0.; rhoint2 = 0.; rhoint3 = 0.
do i=1, quad_order
q = quad_x(i); w_q = quad_w(i) !!!VARIABLE HERE IS Q!!!

delta_q = pi/2. -atan(lambda -2.*q/g)
n_q = 1./(1. +exp(beta*(q**2. -mu)))
nu_q = (1./(2.*pi*jj)) * log(1. +n_q*(exp(2.*jj*delta_q) -1.))


gprime = 2.*exp(beta*(q**2. -mu))*(4.* q**3.*beta -2.*g*(jj +2.*q**2.*beta*lambda) +g**2.*q*beta*(1.+lambda**2.)) / &
        ((exp(beta*(q**2.-mu))-exp(2.*jj*atan(2*q/g-lambda)))*(4.*q**2. +g*(g-4.*q*lambda+g*lambda**2.)))

rhoint1 = rhoint1 + nu_q*w_q
rhoint2 = rhoint2 + nu_q*gprime*w_q


intJ_x=0.
do j=1, quad_order
p = quad_x(j); w_p = quad_w(j) !!!VARIABLE HERE IS P!!!

delta_p = pi/2. -atan(lambda -2.*p/g)
n_p = 1./(1. +exp(beta*(p**2. -mu)))
nu_p = (1./(2.*pi*jj)) * log(1. +n_p*(exp(2.*jj*delta_p) -1.))

if(p.eq.q) then
      arg_int = (2.*exp(beta*(p**2. +mu))*p*beta/((exp(q**2. *beta) +exp(beta*mu))*(exp(p**2. *beta) &
              -exp(-beta*mu +2.*jj*atan(p*alpha - lambda)))) + alpha/(-jj+p*alpha+exp(beta*(p**2. -mu))*(jj+p*alpha-lambda) &
              -lambda)) /(pi*(jj+p*alpha-lambda))
      arg_int_J = (2.*exp(beta*(p**2. +mu))*p*beta/((exp(q**2. *beta) +exp(beta*mu))*(exp(p**2. *beta) &
              -exp(-beta*mu +2.*jj*atan(p*alpha - lambda)))) + alpha/(-jj+p*alpha+exp(beta*(p**2. -mu))*(jj+p*alpha-lambda) &
              -lambda)) /(pi*(jj+p*alpha-lambda))
else
        arg_int = (nu_p-nu_q)/(p-q)
        arg_int_J = nu_p/(p-q)
endif

rhoint3 = rhoint3 + arg_int**2. *w_p*w_q !DOUBLE INTEGRAL
intJ_x = intJ_x + arg_int_J*w_p
enddo


J_x = J_x + (1./pi) * n_q*(sin(delta_q)**2.)*exp(-jj*x*q -2.*intJ_x)*w_q
enddo

rho = J_x*exp(jj*x*rhoint1 -rhoint2 -rhoint3/2.)
endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!JUMPY NU!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (abs(lambda).lt.alpha*sqrt(mu)) then 

J_x = 0.; rhoint1 = 0.; rhoint2 = 0.; rhoint3 = 0.
do i=1, quad_order
q = quad_x(i); w_q = quad_w(i) !!!VARIABLE HERE IS Q!!!

delta_q = pi/2. -atan(lambda -2.*q/g)
n_q = 1./(1. +exp(beta*(q**2. -mu)))
nu_q = (1./(2.*pi*jj)) * log(1. +n_q*(exp(2.*jj*delta_q) -1.))

gfun = 2.*pi*jj*nu_q-2.*jj*delta_q-log(n_q)

dnu_q = -jj*(exp(beta*(q**2. -mu))*(1.+exp(2.*jj*atan(q*alpha-lambda)))*q*beta &
        +jj*alpha*(1.+exp(beta*(q**2. -mu)))/((jj-q*alpha-lambda)**2.))/ &
        ((1.+exp(beta*q**2. -mu))**2. * (pi -2.*pi*jj/((1.+exp(beta*(q**2.-mu)))*(jj+q*alpha-lambda))))

rhoint1 = rhoint1 + q*dnu_q*w_q
rhoint2 = rhoint2 + gfun*dnu_q*w_q

do j=1, quad_order
p = quad_x(j); w_p = quad_w(j) !!!VARIABLE HERE IS P!!!

delta_p = pi/2. -atan(lambda -2.*p/g)
n_p = 1./(1. +exp(beta*(p**2. -mu)))
nu_p = (1./(2.*pi*jj)) * log(1. +n_p*(exp(2.*jj*delta_p) -1.))


dnu_p = -jj*(exp(beta*(p**2. -mu))*(1.+exp(2.*jj*atan(p*alpha-lambda)))*p*beta &
        +jj*alpha*(1.+exp(beta*(p**2. -mu)))/((jj-p*alpha-lambda)**2.))/ &
        ((1.+exp(beta*p**2. -mu))**2. * (pi -2.*pi*jj/((1.+exp(beta*(p**2.-mu)))*(jj+p*alpha-lambda))))

if(p.eq.q) then
        arg_int = 0.
else
        arg_int = dnu_p*dnu_q*log(abs(p-q))
endif

rhoint3 = rhoint3 + arg_int*w_p*w_q !DOUBLE INTEGRAL
enddo
enddo 

rho = 0.5*exp(-jj*x*rhoint1 +rhoint2 +rhoint3)
endif




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!calc of ec(lambda)!!!!!!!!!!!!!!!!!!!!!!!!!!
ec = 0.
do m=1,quad_order
k = quad_x(m); w = quad_w(m)
fd = 1./(1.+exp(beta*(k**2. -mu)))
ec = ec + (k/pi)*fd*(pi/2. -atan(lambda -2.*k/g))*w
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

rho = rho + exp(-beta*ec)*rho*w_l
enddo

!if(x.eq.0) norm = rho
norm =   (2.44950023737424728E-002,-1.24036027420767100E-020) ! from exact formula, quadrature 101, -5, 5
!print*, x,rho/norm
print*, beta, rho/norm
write(20,*) beta, real(rho/norm)
!write(20,*) x, real(rho/norm), aimag(rho/norm), abs(rho/norm)
!x = x+dx
beta = beta+dbeta
enddo


close(20)

end program
