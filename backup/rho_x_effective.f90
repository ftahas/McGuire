program effective
implicit none
character (len = 6) answer
complex*16 jj,rho,nu_p,nu_q,J_x,arg_int,gprime,norm,dnu_q,arg_int_J
complex*16 intJ_x,rhoint1,rhoint2,rhoint3,rhoint3_2,dnu_p,gfun,rho_l
real*8 pi,delta_p,delta_q,beta,n_p, n_q
real*8 p,q,lambda,g,mu,w_p,w_q,w_l,w_x
real*8 x,dx,ec,k,fd,w,x0,xf,alpha
real*8, allocatable :: quad_w(:), quad_x(:), quad_w_l(:), quad_x_l(:), quad2_w(:), quad2_x(:)
real*8, allocatable ::  quad_ec_x(:), quad_ec_w(:), quad_nk_x(:), quad_nk_w(:), quad1_w(:), quad1_x(:)
integer i,j,l,m,quad_order,io,jx, quad_order1
integer quad_order_ec, quad_order_lambda, quad_order_nk,quad_order2

quad_order=0
open(1,file='quadrature_x.txt')
do
read(1,*,iostat=io)
if(io.ne.0) exit
quad_order = quad_order + 1
enddo
close(1)

quad_order_lambda=0
open(1,file='quad_lambda_x.txt')
do
read(1,*,iostat=io)
if(io.ne.0) exit
quad_order_lambda = quad_order_lambda + 1
enddo
close(1)

quad_order_ec=0
open(1,file='quad_ec_x.txt')
do
read(1,*,iostat=io)
if(io.ne.0) exit
quad_order_ec = quad_order_ec + 1
enddo
close(1)


quad_order1=0
open(1,file='quadrature1_x.txt')
do
read(1,*,iostat=io)
if(io.ne.0) exit
quad_order1 = quad_order1 + 1
enddo
close(1)

quad_order2=0
open(1,file='quadrature2_x.txt')
do
read(1,*,iostat=io)
if(io.ne.0) exit
quad_order2 = quad_order2 + 1
enddo
close(1)

allocate(quad2_w(quad_order2)); allocate(quad2_x(quad_order2))
allocate(quad1_w(quad_order1)); allocate(quad1_x(quad_order1))


allocate(quad_ec_w(quad_order)); allocate(quad_ec_x(quad_order_ec))
allocate(quad_w(quad_order)); allocate(quad_x(quad_order))
allocate(quad_w_l(quad_order_lambda)); allocate(quad_x_l(quad_order_lambda))


pi=acos(-1.d0); jj = cmplx(0.d0,1.d0)

mu=1.
alpha = 1.
beta = 1.
g = 2./alpha

open(20,file='./data/rho_x_effective_a1_b1.dat')
open(22,file='./data/rho_x_effective_quadrature_a1_b1.dat')
open(33,file='./data/rho_x_effective_nonorm_a1_b1.dat')


open(10, file='quadrature_w.txt')
open(11, file='quadrature_x.txt')
open(12, file='quad_lambda_w.txt')
open(13, file='quad_lambda_x.txt')
open(14, file='quad_ec_w.txt')
open(15, file='quad_ec_x.txt')
open(16, file='quadrature2_w.txt')
open(17, file='quadrature2_x.txt')
open(18, file='quadrature1_w.txt')
open(19, file='quadrature1_x.txt')

do i=1,quad_order
read(10,*) quad_w(i)
read(11,*) quad_x(i) 
enddo

do i=1,quad_order_lambda
read(12,*) quad_w_l(i)
read(13,*) quad_x_l(i)
enddo 

do i=1,quad_order_ec
read(14,*) quad_ec_w(i)
read(15,*) quad_ec_x(i)
enddo

do i=1,quad_order2
read(16,*) quad2_w(i)
read(17,*) quad2_x(i)
enddo

do i=1,quad_order1
read(18,*) quad1_w(i)
read(19,*) quad1_x(i)
enddo

!print*, quad_w, quad_x

close(10); close(11)
close(12); close(13)
close(14); close(15)
close(16); close(17)
close(18); close(19)

call calculanorma(quad_order,quad_order_lambda,quad_order_ec,alpha,beta,norm)
print*, norm


print*, "DO U WANT QUADRATURE OR NOT YOUR MOTHERFUCKER?"
read*, answer

if (answer.eq.'no') then
write(*,*) 'Enter x0 and xf'
read(*,*) x0, xf

x = x0
do while(x.le.xf)
if (x .lt. 1.d0) dx = 1d-1
if (x .ge. 1.d0) dx = 1d-1
!do jx=1,quad_order
!x = quad_x(jx)

rho = 0.
do l=1, quad_order_lambda
lambda = quad_x_l(l)
w_l = quad_w_l(l)





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!GOOD NU!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (abs(lambda).ge.alpha*sqrt(mu)) then

J_x = 0.; rhoint1 = 0.; rhoint2 = 0.;rhoint3 = 0.
do i=1, quad_order
q = quad_x(i); w_q = quad_w(i) !!!VARIABLE HERE IS Q!!!

delta_q = pi/2. -atan(lambda -2.*q/g)
n_q = 1./(1. +exp(beta*(q**2. -mu)))
nu_q = (1./(2.*pi*jj)) * log(1. +n_q*(exp(2.*jj*delta_q) -1.))


gprime = 2.*exp(beta*(q**2. -mu))*(4.* q**3.*beta -2.*g*(jj +2.*q**2.*beta*lambda) +g**2.*q*beta*(1.+lambda**2.)) / &
        ((exp(beta*(q**2.-mu))-exp(2.*jj*atan(2*q/g-lambda)))*(4.*q**2. +g*(g-4.*q*lambda+g*lambda**2.)))

!print*, gprime

rhoint1 = rhoint1 + nu_q*w_q
rhoint2 = rhoint2 + nu_q*gprime*w_q

intJ_x=0.
do j=1, quad_order
p = quad_x(j); w_p = quad_w(j) !!!VARIABLE HERE IS P!!!

delta_p = pi/2. -atan(lambda -2.*p/g)
n_p = 1./(1. +exp(beta*(p**2. -mu)))
nu_p = (1./(2.*pi*jj)) * log(1. +n_p*(exp(2.*jj*delta_p) -1.))

if(p.eq.q) then
        arg_int = (2.*exp(beta*mu)*g*(exp(beta*mu)+exp(beta* p**2.)*(1. +p*beta*(2.*p-g*(lambda-jj)))))/ &
               ((exp(beta* p**2.)+exp(beta*mu))*pi *(-2.*p +g*(lambda-jj))* (exp(beta* p**2.)*(-2.*p +g*(lambda-jj)) &
              +exp(beta*mu)*(-2.*p +g*(lambda+jj))))  
       arg_int_J = (2.*exp(beta*(p**2. +mu))*p*beta/((exp(q**2. *beta) +exp(beta*mu))*(exp(p**2. *beta) &
              -exp(-beta*mu +2.*jj*atan(p*alpha - lambda)))) + alpha/(-jj+p*alpha+exp(beta*(p**2. -mu))*(jj+p*alpha-lambda) &
              -lambda)) /(pi*(jj+p*alpha-lambda))
else 
        arg_int = (nu_p - nu_q)/(p-q)
        arg_int_J = nu_p/(p-q)
endif

rhoint3 = rhoint3 + arg_int**2. *w_p*w_q !DOUBLE INTEGRAL
intJ_x = intJ_x + arg_int_J*w_p
enddo

J_x = J_x + (1./pi) * n_q*(sin(delta_q)**2.)*exp(-jj*x*q -2.*intJ_x)*w_q
enddo

rho_l = J_x*exp(jj*x*rhoint1 -rhoint2 -rhoint3/2.)
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!











!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!JUMPY NU!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(abs(lambda).lt.alpha*sqrt(mu)) then

J_x = 0.; rhoint1 = 0.; rhoint2 = 0.; rhoint3 = 0.
do i=1, quad_order
q = quad_x(i); w_q = quad_w(i) !!!VARIABLE HERE IS Q!!!

delta_q = pi/2. -atan(lambda -2.*q/g)
n_q = 1./(1. +exp(beta*(q**2. -mu)))
nu_q = (1./(2.*pi*jj)) * log(1. +n_q*(exp(2.*jj*delta_q) -1.))
!phi_q = 1. + (2.*jj*n_q)/(lambda-alpha*q-jj)

!gfun = 2.*pi*jj*nu_q -2.*jj*delta_q -log(n_q)
gfun = -log((n_q*exp(2.*jj*delta_q))/(exp(2.*pi*jj*nu_q)))

dnu_q = -jj*(exp(beta*(q**2. -mu))*(1.+exp(2.*jj*atan(q*alpha-lambda)))*q*beta &
        +jj*alpha*(1.+exp(beta*(q**2. -mu)))/((jj+q*alpha-lambda)**2.))/ &
        ((1.+exp(beta*q**2. -mu))**2. * (pi -2.*pi*jj/((1.+exp(beta*(q**2.-mu)))*(jj+q*alpha-lambda))))


rhoint1 = rhoint1 + q*dnu_q*w_q
rhoint2 = rhoint2 + gfun*dnu_q*w_q

do j=1, quad_order2
p = quad2_x(j); w_p = quad2_w(j) !!!VARIABLE HERE IS P!!!

dnu_p = -jj*(exp(beta*(p**2. -mu))*(1.+exp(2.*jj*atan(p*alpha-lambda)))*p*beta &
        +jj*alpha*(1.+exp(beta*(p**2. -mu)))/((jj+p*alpha-lambda)**2.))/ &
        ((1.+exp(beta*p**2. -mu))**2. * (pi -2.*pi*jj/((1.+exp(beta*(p**2.-mu)))*(jj+p*alpha-lambda))))

!if(p.eq.q) then
!        arg_int = 0.
!else
        arg_int = dnu_p*dnu_q*log(abs(p-q))
!endif



rhoint3 = rhoint3 + arg_int*w_p*w_q !DOUBLE INTEGRAL
enddo
enddo

rho_l = 0.5*exp(-jj*x*rhoint1 +rhoint2 +rhoint3)
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!calc of ec(lambda)!!!!!!!!!!!!!!!!!!!!!!!!!!
ec = 0.
do m=1,quad_order_ec
k = quad_ec_x(m); w = quad_ec_w(m)
fd = 1./(1.+exp(beta*(k**2. -mu)))
ec = ec -2.*(k/pi)*fd*(pi/2. - atan(lambda -2.*k/g))*w
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


rho = rho + exp(-beta*ec)*rho_l*w_l
enddo !!!!!!!!!!!!!lambda loop sum!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




print*, x, rho/norm
write(20,*) x, real(rho/norm), aimag(rho/norm), abs(rho/norm)
write(33,*) x, real(rho), aimag(rho), abs(rho)
x = x+dx
enddo

endif!!!!!!!!!!!!!!!!!!!!!!!answer yes or no if statement


















if (answer.eq.'yes') then
quad_order_nk=0
open(1,file='quad_nk_x.txt')
do
read(1,*,iostat=io)
if(io.ne.0) exit
quad_order_nk = quad_order_nk + 1
enddo
close(1)

allocate(quad_nk_w(quad_order_nk)); allocate(quad_nk_x(quad_order_nk))



open(10, file='quad_nk_w.txt')
open(11, file='quad_nk_x.txt')

do i=1,quad_order_nk
read(10,*) quad_nk_w(i)
read(11,*) quad_nk_x(i)
enddo

do jx=1,quad_order_nk
x = quad_nk_x(jx); w_x = quad_nk_w(jx)

rho = 0.
do l=1, quad_order_lambda
lambda = quad_x_l(l)
w_l = quad_w_l(l)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!GOOD NU!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (abs(lambda).ge.alpha*sqrt(mu)) then

J_x = 0.; rhoint1 = 0.; rhoint2 = 0.;rhoint3 = 0.
do i=1, quad_order
q = quad_x(i); w_q = quad_w(i) !!!VARIABLE HERE IS Q!!!

delta_q = pi/2. -atan(lambda -2.*q/g)
n_q = 1./(1. +exp(beta*(q**2. -mu)))
nu_q = (1./(2.*pi*jj)) * log(1. +n_q*(exp(2.*jj*delta_q) -1.))


gprime = 2.*exp(beta*(q**2. -mu))*(4.* q**3.*beta -2.*g*(jj +2.*q**2.*beta*lambda) +g**2.*q*beta*(1.+lambda**2.)) / &
        ((exp(beta*(q**2.-mu))-exp(2.*jj*atan(2*q/g-lambda)))*(4.*q**2. +g*(g-4.*q*lambda+g*lambda**2.)))

!print*, gprime

rhoint1 = rhoint1 + nu_q*w_q
rhoint2 = rhoint2 + nu_q*gprime*w_q

intJ_x=0.
do j=1, quad_order
p = quad_x(j); w_p = quad_w(j) !!!VARIABLE HERE IS P!!!

delta_p = pi/2. -atan(lambda -2.*p/g)
n_p = 1./(1. +exp(beta*(p**2. -mu)))
nu_p = (1./(2.*pi*jj)) * log(1. +n_p*(exp(2.*jj*delta_p) -1.))

if(p.eq.q) then
        arg_int = (2.*exp(beta*mu)*g*(exp(beta*mu)+exp(beta* p**2.)*(1. +p*beta*(2.*p-g*(lambda-jj)))))/ &
               ((exp(beta* p**2.)+exp(beta*mu))*pi *(-2.*p +g*(lambda-jj))* (exp(beta* p**2.)*(-2.*p +g*(lambda-jj)) &
              +exp(beta*mu)*(-2.*p +g*(lambda+jj))))  
       arg_int_J = (2.*exp(beta*(p**2. +mu))*p*beta/((exp(q**2. *beta) +exp(beta*mu))*(exp(p**2. *beta) &
              -exp(-beta*mu +2.*jj*atan(p*alpha - lambda)))) + alpha/(-jj+p*alpha+exp(beta*(p**2. -mu))*(jj+p*alpha-lambda) &
              -lambda)) /(pi*(jj+p*alpha-lambda))
else 
        arg_int = (nu_p - nu_q)/(p-q)
        arg_int_J = nu_p/(p-q)
endif

rhoint3 = rhoint3 + arg_int**2. *w_p*w_q !DOUBLE INTEGRAL
intJ_x = intJ_x + arg_int_J*w_p
enddo

J_x = J_x + (1./pi) * n_q*(sin(delta_q)**2.)*exp(-jj*x*q -2.*intJ_x)*w_q
enddo

rho_l = J_x*exp(jj*x*rhoint1 -rhoint2 -rhoint3/2.)
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!











!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!JUMPY NU!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(abs(lambda).lt.alpha*sqrt(mu)) then

J_x = 0.; rhoint1 = 0.; rhoint2 = 0.; rhoint3 = 0.
do i=1, quad_order
q = quad_x(i); w_q = quad_w(i) !!!VARIABLE HERE IS Q!!!

delta_q = pi/2. -atan(lambda -2.*q/g)
n_q = 1./(1. +exp(beta*(q**2. -mu)))
nu_q = (1./(2.*pi*jj)) * log(1. +n_q*(exp(2.*jj*delta_q) -1.))
!phi_q = 1. + (2.*jj*n_q)/(lambda-alpha*q-jj)

!gfun = 2.*pi*jj*nu_q -2.*jj*delta_q -log(n_q)
gfun = -log((n_q*exp(2.*jj*delta_q))/(exp(2.*pi*jj*nu_q)))

dnu_q = -jj*(exp(beta*(q**2. -mu))*(1.+exp(2.*jj*atan(q*alpha-lambda)))*q*beta &
        +jj*alpha*(1.+exp(beta*(q**2. -mu)))/((jj+q*alpha-lambda)**2.))/ &
        ((1.+exp(beta*q**2. -mu))**2. * (pi -2.*pi*jj/((1.+exp(beta*(q**2.-mu)))*(jj+q*alpha-lambda))))


rhoint1 = rhoint1 + q*dnu_q*w_q
rhoint2 = rhoint2 + gfun*dnu_q*w_q

do j=1, quad_order2
p = quad2_x(j); w_p = quad2_w(j) !!!VARIABLE HERE IS P!!!

dnu_p = -jj*(exp(beta*(p**2. -mu))*(1.+exp(2.*jj*atan(p*alpha-lambda)))*p*beta &
        +jj*alpha*(1.+exp(beta*(p**2. -mu)))/((jj+p*alpha-lambda)**2.))/ &
        ((1.+exp(beta*p**2. -mu))**2. * (pi -2.*pi*jj/((1.+exp(beta*(p**2.-mu)))*(jj+p*alpha-lambda))))

!if(p.eq.q) then
!        arg_int = 0.
!else
        arg_int = dnu_p*dnu_q*log(abs(p-q))
!endif



rhoint3 = rhoint3 + arg_int*w_p*w_q !DOUBLE INTEGRAL
enddo
enddo

rho_l = 0.5*exp(-jj*x*rhoint1 +rhoint2 +rhoint3)
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!calc of ec(lambda)!!!!!!!!!!!!!!!!!!!!!!!!!!
ec = 0.
do m=1,quad_order_ec
k = quad_ec_x(m); w = quad_ec_w(m)
fd = 1./(1.+exp(beta*(k**2. -mu)))
ec = ec -2.*(k/pi)*fd*(pi/2. - atan(lambda -2.*k/g))*w
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


rho = rho + exp(-beta*ec)*rho_l*w_l
enddo !!!!!!!!!!!!!lambda loop sum!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print*, x, rho/norm
write(20,*) x, real(rho/norm), aimag(rho/norm), abs(rho/norm)
write(22,*) x, w_x, rho/norm !, real(rho/norm), aimag(rho/norm), abs(rho/norm)
enddo
endif

close(20); close(22)
end program





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calculanorma(quad_order,quad_order_lambda,quad_order_ec,alpha,beta,norm)
implicit none
character (len = 6) answer
complex*16 jj,ep1,ep2,em1,em2,KK,WW,rho,norm,rho_l 
complex*16, allocatable :: KpW(:,:), IK(:,:)
real*8 pi,delta1,delta2,beta,fd1,fd2,x0,xf
real*8 w1,w2,q1,q2,lambda,dlambda,g,mu,Z,alpha 
real*8 x,dx,k,dk,ec,fd,w,w_l,w_x
real*8, allocatable :: quad_w(:), quad_x(:), quad_w_2(:), quad_x_2(:)
real*8, allocatable :: II(:,:),  quad_ec_x(:), quad_ec_w(:), quad_w_l(:), quad_x_l(:), quad_nk_x(:), quad_nk_w(:)
integer i,j,l,m,quad_order,io,jx
integer quad_order_ec, quad_order_lambda, quad_order_nk


pi=acos(-1.d0); jj = cmplx(0.d0,1.d0)
mu=1.d0
g = 2./alpha

allocate(quad_ec_w(quad_order_ec)); allocate(quad_ec_x(quad_order_ec))
allocate(quad_w_l(quad_order_lambda)); allocate(quad_x_l(quad_order_lambda))
allocate(quad_w(quad_order)); allocate(quad_x(quad_order))
allocate(II(quad_order,quad_order)); allocate(KpW(quad_order,quad_order))
allocate(IK(quad_order,quad_order))

open(10, file='quadrature_w.txt')
open(11, file='quadrature_x.txt')

do i=1,quad_order
read(10,*) quad_w(i)
read(11,*) quad_x(i)
do j=1,quad_order
if (i.eq.j) II(i,j)=1.
if (i.ne.j) II(i,j)=0.
enddo 
enddo 
close(10); close(11)

open(12, file='quad_lambda_w.txt')
open(13, file='quad_lambda_x.txt')
open(14, file='quad_ec_w.txt')
open(15, file='quad_ec_x.txt')

do i=1,quad_order_lambda
read(12,*) quad_w_l(i)
read(13,*) quad_x_l(i)
enddo

do i=1,quad_order_ec
read(14,*) quad_ec_w(i)
read(15,*) quad_ec_x(i)
enddo

close(12); close(13)
close(14); close(15)

x = 0.
rho=0.
do l=1, quad_order_lambda
lambda = quad_x_l(l)
w_l = quad_w_l(l)

do i=1, quad_order
do j=1, quad_order

q1 = quad_x(i); q2 = quad_x(j)
w1 = quad_w(i); w2 = quad_w(j)

delta1 = pi/2. - atan(lambda-2.*q1/g); delta2 = pi/2. - atan(lambda-2.*q2/g)

ep1 = (1./pi) * exp(jj*q1*x/2. +jj*delta1); ep2 = (1./pi) * exp(jj*q2*x/2. +jj*delta2)
em1 = exp(-jj*q1*x/2.)*sin(delta1); em2 = exp(-jj*q2*x/2.)*sin(delta2)

if (q1.ne.q2) then
                KK = 1./(q1-q2) * (ep1*em2-em1*ep2)
        else
                KK = exp(jj*atan(q2*2./g -lambda))*(jj*2./g + x*(jj -q2*2./g + lambda)) / &
                        (pi*(-jj +q2*2./g -lambda)*sqrt(1. +(lambda -q2*2./g)**2.))
        end if
!Z = 1./(alpha*pi) *(atan(alpha-lambda)-atan(alpha+lambda))
Z = 1.
WW = 1./(Z*pi) * em1*em2

fd1 = 1./(1.+exp(beta*(q1**2. -mu))); fd2 = 1./(1.+exp(beta*(q2**2. -mu)))

KpW(i,j) = II(i,j) + sqrt(w1*fd1)*(KK+WW)*sqrt(w2*fd2) 
IK(i,j) = II(i,j) + sqrt(w1*fd1)*KK*sqrt(w2*fd2)
enddo
enddo
rho_l = det(KpW)-det(IK)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!calc of ec(lambda)!!!!!!!!!!!!!!!!!!!!!!!!!!
ec = 0.
do m=1,quad_order_ec
k = quad_ec_x(m); w = quad_ec_w(m)
fd = 1./(1.+exp(beta*(k**2. -mu)))
ec = ec -2.*(k/pi)*fd*(pi/2. - atan(lambda -2.*k/g))*w
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


rho = rho + exp(-beta*ec)*rho_l*w_l
!print*, rho
enddo
norm = rho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!DETERMINANT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
complex*16 function det(aa)
complex*16 aa(:,:),tmp,c(size(aa,dim=1),size(aa,dim=2)),maxx
integer i,j,k,l,n,m,num(size(aa,dim=1))
n=size(aa,dim=1)
det=1.
do k=1,n
maxx=aa(k,k);num(k)=k;
do i=k+1,n
if(abs(maxx)<abs(aa(i,k))) then
maxx=aa(i,k)
num(k)=i
endif
enddo
if (num(k)/=k) then
do l=k,n
tmp=aa(k,l)
aa(k,l)=aa(num(k),l)
aa(num(k),l)=tmp
enddo
det=-1.*det
endif
do m=k+1,n
c(m,k)=aa(m,k)/aa(k,k)
do l=k,n
aa(m,l)=aa(m,l)-c(m,k)*aa(k,l)
enddo
enddo !There we made matrix triangular!
enddo

do i=1,n
det=det*aa(i,i)
enddo
return
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

endsubroutine
