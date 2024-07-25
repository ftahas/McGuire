program correlation
implicit none
complex*16 jj,rho,nu_p,nu_q,J_x,arg_int,gprime,dnu_q,dnu_p,arg_int_J
complex*16 intJ_x,rhoint1,rhoint2,rhoint3,gfun,rhoint3_2,phi_q,DI,K1,K2
complex*16, allocatable :: KR1(:,:), KR2(:,:)
real*8 pi,delta_p,delta_q,beta,n_p, n_q
real*8 p,q,lambda,g,mu,w_p,w_q
real*8 x,dx,x0,alpha,xi,xf,dp,dq
real*8, allocatable :: quad_w(:), quad_x(:), quad2_w(:), quad2_x(:),II(:,:)
real*8, allocatable :: quad1_w(:), quad1_x(:)
integer i,j,l,quad_order,io,quad_order2, quad_order1

pi=acos(-1.d0); jj = cmplx(0.d0,1.d0)

lambda = 8.
beta = 1.
alpha = 10.
g = 2./alpha
mu=1.

open(20,file='./data/rho_x_lambda8_a10_b1_effective.dat')

print*, 'alpha=', alpha
print*, 'beta=', beta
print*, 'lambda=', lambda

quad_order=0
open(1,file='quadrature1_x.txt')
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
 
allocate(quad_w(quad_order)); allocate(quad_x(quad_order))
allocate(quad2_w(quad_order2)); allocate(quad2_x(quad_order2))
allocate(quad1_w(quad_order1)); allocate(quad1_x(quad_order1))


allocate(II(quad_order,quad_order)); allocate(KR2(quad_order,quad_order))
allocate(KR1(quad_order,quad_order))

open(10, file='quadrature_w.txt')
open(11, file='quadrature_x.txt')

open(12, file='quadrature1_w.txt')
open(13, file='quadrature1_x.txt')

open(14, file='quadrature2_w.txt')
open(15, file='quadrature2_x.txt')

do i=1,quad_order
read(10,*) quad_w(i)
read(11,*) quad_x(i) 
do j=1,quad_order
if (i.eq.j) II(i,j)=1.
if (i.ne.j) II(i,j)=0.
enddo
enddo 

do i=1,quad_order1
read(12,*) quad1_w(i)
read(13,*) quad1_x(i)
enddo

do i=1,quad_order2
read(14,*) quad2_w(i)
read(15,*) quad2_x(i)
enddo

!print*, quad_w, quad_x

close(10); close(11)
close(12); close(13)
close(14); close(15)

!x = 0.; x0 = 20.; dx = 1d-1

write(*,*) 'enter x0 and xf'
read(*,*) x0, xf

x=x0; dx = 1d-1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!GOOD NU!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
if (abs(lambda).gt.alpha*sqrt(mu)) then
do while(x.le.xf)
if (x.le.0.5) dx = 1d-2
if (x.gt.1d0) dx = 1d-1

rho = 0.


J_x = 0.; rhoint1 = 0.; rhoint2 = 0.;rhoint3 = 0.
do i=1, quad_order1
q = quad1_x(i); w_q = quad1_w(i) !!!VARIABLE HERE IS Q!!!

delta_q = pi/2. -atan(lambda -2.*q/g)
n_q = 1./(1. +exp(beta*(q**2. -mu)))
nu_q = (1./(2.*pi*jj)) * log(1. +n_q*(exp(2.*jj*delta_q) -1.))


gprime = 2.*exp(beta*(q**2. -mu))*(4.* q**3.*beta -2.*g*(jj +2.*q**2.*beta*lambda) +g**2.*q*beta*(1.+lambda**2.)) / &
        ((exp(beta*(q**2.-mu))-exp(2.*jj*atan(2*q/g-lambda)))*(4.*q**2. +g*(g-4.*q*lambda+g*lambda**2.)))

!print*, gprime

rhoint1 = rhoint1 + nu_q*w_q
rhoint2 = rhoint2 + nu_q*gprime*w_q

intJ_x=0.
do j=1, quad_order2
p = quad2_x(j); w_p = quad2_w(j) !!!VARIABLE HERE IS P!!!

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
        arg_int_J = (nu_p - nu_q)/(p-q)
endif

rhoint3 = rhoint3 + arg_int**2. *w_p*w_q !DOUBLE INTEGRAL
intJ_x = intJ_x + arg_int_J*w_p
enddo
intJ_x = intJ_x + nu_q*log(abs((p-q)/(p+q))) !CORRECTION IN PV INTEGRAL

J_x = J_x + (1./pi) * n_q*(sin(delta_q)**2.)*exp(-jj*x*q -2.*intJ_x)*w_q
enddo

!print*, rho

rho = J_x*exp(jj*x*rhoint1 -rhoint2 -rhoint3/2.)

print*, x, rho
write(20,*) x, real(rho), aimag(rho), abs(rho)
x = x+dx
enddo 
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!JUMPY NU!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
if (abs(lambda).le.alpha*sqrt(mu)) then  
do while(x.le.xf)
if (x.le.0.5) dx = 1d-2
if (x.gt.1.) dx = 1d-1

rho = 0.

rhoint1 = 0.; rhoint2 = 0.; rhoint3 = 0.

do i=1, quad_order
q = quad_x(i); w_q = quad_w(i) !!!VARIABLE HERE IS Q!!!

!q=-5.;dq=1d-2
!do while(q.le.5)

delta_q = pi/2. -atan(lambda -2.*q/g)
n_q = 1./(1. +exp(beta*(q**2. -mu)))
nu_q = (1./(2.*pi*jj)) * log(1. +n_q*(exp(2.*jj*delta_q) -1.))
phi_q = 1. + (2.*jj*n_q)/(lambda-alpha*q-jj)

gfun = -log((n_q*exp(2.*jj*delta_q))/(exp(2.*pi*jj*nu_q)))

dnu_q = -jj*(exp(beta*(q**2. -mu))*(1.+exp(2.*jj*atan(q*alpha-lambda)))*q*beta &
        +jj*alpha*(1.+exp(beta*(q**2. -mu)))/((jj+q*alpha-lambda)**2.))/ &
        ((1.+exp(beta*q**2. -mu))**2. * (pi -2.*pi*jj/((1.+exp(beta*(q**2.-mu)))*(jj+q*alpha-lambda)))) 


rhoint1 = rhoint1 + q*dnu_q*w_q
rhoint2 = rhoint2 + gfun*dnu_q*w_q

!p=-5.1;dp=1d-2
!do while(p.le.5.1)
!do j=1, quad_order2
!p = quad2_x(j); w_p = quad2_w(j) !!!VARIABLE HERE IS P!!!

!dnu_p = -jj*(exp(beta*(p**2. -mu))*(1.+exp(2.*jj*atan(p*alpha-lambda)))*p*beta &
 !       +jj*alpha*(1.+exp(beta*(p**2. -mu)))/((jj+p*alpha-lambda)**2.))/ &
 !       ((1.+exp(beta*p**2. -mu))**2. * (pi -2.*pi*jj/((1.+exp(beta*(p**2.-mu)))*(jj+p*alpha-lambda))))

!if(p.eq.q) then
!        arg_int = 0.
!else
  !      arg_int = dnu_p*dnu_q*log(abs(p-q))
!endif

!rhoint3 = rhoint3 + arg_int*dp*dq!w_p*w_q !DOUBLE INTEGRAL
!p=p+dp
!enddo
!q=q+dq
enddo


do i=1, quad_order
do j=1, quad_order
p = quad_x(i); q = quad_x(j)
w_p = quad_w(i); w_q = quad_w(j)

delta_p = pi/2. -atan(lambda -alpha*p)
n_p = 1./(1. +exp(beta*(p**2. -mu)))
nu_p = 1./(2.*pi*jj) * log(1. +n_p*(exp(2.*jj*delta_p) -1.))


delta_q = pi/2. -atan(lambda -2.*q/g)
n_q = 1./(1. +exp(beta*(q**2. -mu)))
nu_q = (1./(2.*pi*jj)) * log(1. +n_q*(exp(2.*jj*delta_q) -1.))


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


!rho = 0.5*exp(-jj*x*rhoint1 +rhoint2 +rhoint3)

rho = DI*exp(-jj*x*rhoint1 +rhoint2)
!print*, x, log(2.*DI), rhoint2 
print*, x, rho
write(20,*) x, real(rho), aimag(rho), abs(rho)
x = x+dx
enddo 
endif

close(20)

 
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
