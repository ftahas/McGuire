program effective
implicit none
complex*16 jj,rho,nu_p,nu_q,J_x,arg_int,gprime,dnu_q,arg_int_J,rho_l
complex*16 intJ_x,rhoint1,rhoint2,rhoint3,dnu_p,gfun,DI,K1,K2
complex*16, allocatable :: KR1(:,:), KR2(:,:)
real*8 pi,delta_p,delta_q,beta,n_p, n_q,dlambda
real*8 p,q,lambda,g,mu,w_p,w_q,w_l
real*8 ec,k,fd,w,alpha,x,dq
real*8, allocatable :: quad_w(:), quad_x(:), quad_w_l(:), quad_x_l(:), quad2_w(:), quad2_x(:), II(:,:)
real*8, allocatable ::  quad_ec_x(:), quad_ec_w(:), quad_nk_x(:), quad_nk_w(:), quad1_w(:), quad1_x(:)
integer i,j,l,m,quad_order,io,jx, quad_order1
integer quad_order_ec, quad_order_lambda, quad_order_nk,quad_order2

pi=acos(-1.d0); jj = cmplx(0.d0,1.d0)

mu = 1d0 
alpha = 1d-20
beta = 10d0
g = 2.d0/alpha

open(22,file='./data/rho_lambda_effective.dat')


quad_order=0
open(1,file='quadrature_x.txt')
do
read(1,*,iostat=io)
if(io.ne.0) exit
quad_order = quad_order + 1
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


quad_order_lambda=0
open(1,file='quad_lambda_x.txt')
do
read(1,*,iostat=io)
if(io.ne.0) exit
quad_order_lambda = quad_order_lambda + 1
enddo
close(1)

allocate(quad_w_l(quad_order_lambda)); allocate(quad_x_l(quad_order_lambda))

allocate(II(quad_order,quad_order)); allocate(KR2(quad_order,quad_order))
allocate(KR1(quad_order,quad_order))

allocate(quad2_w(quad_order2)); allocate(quad2_x(quad_order2))
allocate(quad1_w(quad_order1)); allocate(quad1_x(quad_order1))

allocate(quad_w(quad_order)); allocate(quad_x(quad_order))

open(10, file='quadrature_w.txt')
open(11, file='quadrature_x.txt')
open(12, file='quad_lambda_w.txt')
open(13, file='quad_lambda_x.txt')
open(16, file='quadrature2_w.txt')
open(17, file='quadrature2_x.txt')
open(18, file='quadrature1_w.txt')
open(19, file='quadrature1_x.txt')

do i=1,quad_order
read(10,*) quad_w(i)
read(11,*) quad_x(i) 
do j=1,quad_order
if (i.eq.j) II(i,j)=1.
if (i.ne.j) II(i,j)=0.
enddo
enddo


do i=1,quad_order_lambda
read(12,*) quad_w_l(i)
read(13,*) quad_x_l(i)
enddo 

do i=1,quad_order2
read(16,*) quad2_w(i)
read(17,*) quad2_x(i)
enddo

do i=1,quad_order1
read(18,*) quad1_w(i)
read(19,*) quad1_x(i)
enddo


close(10); close(11)
close(12); close(13)
close(16); close(17)
close(18); close(19)

x = 0.
do l=1, quad_order_lambda
lambda = quad_x_l(l)

rho_l=0.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!GOOD NU!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (abs(lambda).ge.alpha*sqrt(mu)) then

!print*, 'GOOD NU'

rhoint1 = 0.; rhoint2 = 0.;rhoint3 = 0.
do i=1, quad_order
q = quad_x(i); w_q = quad_w(i) !!!VARIABLE HERE IS Q!!!

delta_q = pi/2. -datan(lambda -2.*q/g)
n_q = 1./(1. +dexp(beta*(q**2. -mu)))
nu_q = (1./(2.*pi*jj)) * cdlog(1. +n_q*(exp(2.*jj*delta_q) -1.))


gprime = 2.*dexp(beta*(q**2. -mu))*(4.* q**3.*beta -2.*g*(jj +2.*q**2.*beta*lambda) +g**2.*q*beta*(1.+lambda**2.)) / &
        ((dexp(beta*(q**2.-mu))-exp(2.*jj*datan(2*q/g-lambda)))*(4.*q**2. +g*(g-4.*q*lambda+g*lambda**2.)))

!print*, gprime

rhoint1 = rhoint1 + nu_q*w_q
rhoint2 = rhoint2 + nu_q*gprime*w_q

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
        arg_int = (nu_p-nu_q)/(p-q)
endif

rhoint3 = rhoint3 + arg_int**2. *w_p*w_q !DOUBLE INTEGRAL
enddo
enddo




J_x = 0.
do i=1, quad_order
q = quad1_x(i); w_q = quad1_w(i) !!!VARIABLE HERE IS Q!!!

delta_q = pi/2. -atan(lambda -2.*q/g)
n_q = 1./(1. +exp(beta*(q**2. -mu)))
nu_q = (1./(2.*pi*jj)) * log(1. +n_q*(exp(2.*jj*delta_q) -1.))


intJ_x=0.
do j=1, quad_order2
p = quad2_x(j); w_p = quad2_w(j) !!!VARIABLE HERE IS P!!!

delta_p = pi/2. -atan(lambda -2.*p/g)
n_p = 1./(1. +exp(beta*(p**2. -mu)))
nu_p = (1./(2.*pi*jj)) * log(1. +n_p*(exp(2.*jj*delta_p) -1.))

if(p.eq.q) then              
        arg_int_J = (2.*exp(beta*(p**2. +mu))*p*beta/((exp(q**2. *beta) +exp(beta*mu))*(exp(p**2. *beta) &
              -exp(-beta*mu +2.*jj*atan(p*alpha - lambda)))) + alpha/(-jj+p*alpha+exp(beta*(p**2. -mu))*(jj+p*alpha-lambda) &
              -lambda)) /(pi*(jj+p*alpha-lambda))
else 
        arg_int_J = (nu_p-nu_q)/(p-q)
endif

intJ_x = intJ_x + arg_int_J*w_p
enddo
intJ_x = intJ_x + nu_q*log(abs((p-q)/(p+q))) !correction in PV integral

J_x = J_x + (1./pi) * n_q*(sin(delta_q)**2.)*exp(-jj*x*q -2.*intJ_x)*w_q
enddo







rho_l = J_x*exp(jj*x*rhoint1 -rhoint2 -rhoint3/2.)

print*, 'x=', x, 'lambda=', lambda, 'J(x)=', J_x, 'rhoint1=', rhoint1, 'rhoint2=', rhoint2, 'rhoint3=', rhoint3
print*, '------------------------------------------------------------------------------------------------------'
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!











!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!JUMPY NU!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(abs(lambda).lt.alpha*sqrt(mu)) then

       !print*, 'JUMPY NU'


rhoint1 = 0.; rhoint2 = 0.; rhoint3 = 0.
do i=1, quad_order1
q = quad1_x(i); w_q = quad1_w(i) !!!VARIABLE HERE IS Q!!!
!q=-3.; dq=1d-6
!do while (q.le.3.)


delta_q = pi/2. -atan(lambda -2.*q/g)
n_q = 1./(1. +exp(beta*(q**2. -mu)))
nu_q = (1./(2.*pi*jj)) * log(1. +n_q*(exp(2.*jj*delta_q) -1.))

!if(q.gt.lambda/alpha) nu_q = nu_q +1.

gfun = -log((n_q*exp(2.*jj*delta_q))/(exp(2.*pi*jj*nu_q)))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!BAD DNU!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!dnu_q = -jj*(exp(beta*(q**2. -mu))*(1.+exp(2.*jj*atan(q*alpha-lambda)))*q*beta &
!        +jj*alpha*(1.+exp(beta*(q**2. -mu)))/((jj+q*alpha-lambda)**2.))/ &
!        ((1.+exp(beta*q**2. -mu))**2. * (pi -2.*pi*jj/((1.+exp(beta*(q**2.-mu)))*(jj+q*alpha-lambda))))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!BAD DNU!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
dnu_q = (alpha +exp(beta*(q**2. -mu))*(alpha+2.*alpha*beta*q**2. -2.*q*beta*(lambda-jj))) &
        /((1.+exp(beta*(q**2. -mu)))*pi*(jj+q*alpha-lambda)*(-jj+q*alpha+exp(beta*(q**2.-mu))*(jj+q*alpha-lambda)-lambda))





rhoint1 = rhoint1 + q*dnu_q*w_q
rhoint2 = rhoint2 + gfun*dnu_q*w_q

!q=q+dq

!do j=1, quad_order2
!p = quad2_x(j); w_p = quad2_w(j) !!!VARIABLE HERE IS P!!!

!dnu_p = -jj*(exp(beta*(p**2. -mu))*(1.+exp(2.*jj*atan(p*alpha-lambda)))*p*beta &
 !       +jj*alpha*(1.+exp(beta*(p**2. -mu)))/((jj+p*alpha-lambda)**2.))/ &
 !       ((1.+exp(beta*p**2. -mu))**2. * (pi -2.*pi*jj/((1.+exp(beta*(p**2.-mu)))*(jj+p*alpha-lambda))))

!if(p.eq.q) then
!        arg_int = 0.
!else
 !       arg_int = dnu_p*dnu_q*log(abs(p-q))
!endif



!rhoint3 = rhoint3 + arg_int*w_p*w_q !DOUBLE INTEGRAL
!enddo
enddo




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PREFACTOR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
do i=1, quad_order
do j=1, quad_order
p = quad_x(i); q = quad_x(j)
w_p = quad_w(i); w_q = quad_w(j)

delta_p = pi/2. -atan(lambda -alpha*p)
n_p = 1./(1. +exp(beta*(p**2. -mu)))
nu_p = 1./(2.*pi*jj) * log(1. +n_p*(exp(2.*jj*delta_p) -1.))

!if(p.gt.lambda/alpha) nu_p = nu_p +1.


delta_q = pi/2. -atan(lambda -2.*q/g)
n_q = 1./(1. +exp(beta*(q**2. -mu)))
nu_q = (1./(2.*pi*jj)) * log(1. +n_q*(exp(2.*jj*delta_q) -1.))

!if(q.gt.lambda/alpha) nu_q = nu_q +1.


if (p.eq.q) then
                K1 = (-(1./pi) *(alpha + beta*(exp(beta*(q**2.-mu))*(2.*jj*q +2.*alpha*q**2. +alpha/beta -2.*q*lambda)))) &
                        *sin(pi*nu_p)*sin(pi*nu_q)
                K2 = (-(1./pi) *(alpha + beta*(exp(beta*(q**2.-mu))*(2.*jj*q +2.*alpha*q**2. +alpha/beta -2.*q*lambda)))) &
                        *sin(pi*nu_p)*sin(pi*nu_q) +sin(pi*nu_p)*sin(pi*nu_q)/pi
        else
                K1 = (-(1./pi) *((p-q)*alpha +exp(beta*(p**2.-mu))*(jj+p*alpha -lambda) &
                        +exp(beta*(q**2.-mu))*(lambda-jj-q*alpha))/(p-q)) *sin(pi*nu_p)*sin(pi*nu_q)
                K2 = (-(1./pi) *((p-q)*alpha +exp(beta*(p**2.-mu))*(jj+p*alpha -lambda) &
                        +exp(beta*(q**2.-mu))*(lambda-jj-q*alpha))/(p-q))*sin(pi*nu_p)*sin(pi*nu_q) +sin(pi*nu_p)*sin(pi*nu_q)/pi
        end if

KR1(i,j) = II(i,j) + sqrt(w_p)*K1*sqrt(w_q)
KR2(i,j) = II(i,j) + sqrt(w_p)*K2*sqrt(w_q)
enddo
enddo
DI = det(KR2)-det(KR1)

!rho_l = 0.5*exp(-jj*x*rhoint1 +rhoint2 +rhoint3)
rho_l = DI*exp(-jj*x*rhoint1 +rhoint2)

print*, 'lambda=', lambda, 'DI=', DI, 'rhoint1=', rhoint1, 'rhoint2=', rhoint2 !, 'rhoint3=', rhoint3
print*, '-------------------------------------------------------------------------'
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!REGION THAT FUCKED MY LIFE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!if(abs(lambda).ge.alpha*sqrt(mu)-1. .AND. abs(lambda).le.1.1*alpha*sqrt(mu)+1.) then
!       call calculaex(quad_order,alpha,beta,lambda,x,rho_l)
!endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!calc of ec(lambda)!!!!!!!!!!!!!!!!!!!!!!!!!!
ec = 0.
do m=1,quad_order
k = quad_x(m); w = quad_w(m)
fd = 1./(1.+exp(beta*(k**2. -mu)))
ec = ec -2.* (k/pi)*fd*(pi/2. - atan(lambda -alpha*k))*w
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

rho = exp(-beta*ec)*rho_l
!print*, lambda, rho
write(22,*) lambda, real(rho)
enddo



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



end program 








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calculaex(quad_order,alpha,beta,lambda,x,rho_l)
implicit none
complex*16 jj,ep1,ep2,em1,em2,KK,WW,rho,norm,rho_l 
complex*16, allocatable :: KpW(:,:), IK(:,:)
real*8 pi,delta1,delta2,beta,fd1,fd2,x0,xf
real*8 w1,w2,q1,q2,lambda,dlambda,g,mu,Z,alpha 
real*8 x,dx,k,dk,ec,fd,w,w_l,w_x
real*8, allocatable :: quad_w(:), quad_x(:), II(:,:)
integer i,j,l,m,quad_order,io,jx


pi=acos(-1.d0); jj = cmplx(0.d0,1.d0)
mu=1.d0
g = 2./alpha

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

rho_l=0.


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
