program exact
implicit none
character (len = 6) answer
complex*16 jj,ep1,ep2,em1,em2,KK,WW,rho,norm,rho_l
complex*16, allocatable :: KpW(:,:), IK(:,:)
real*8 pi,delta1,delta2,beta,fd1,fd2,x0,xf
real*8 w1,w2,q1,q2,lambda,dlambda,g,mu,Z,alpha 
real*8 x,dx,k,dk,ec,fd,w,w_l,w_x
real*8, allocatable :: quad_w(:), quad_x(:), quad_w_l(:), quad_x_l(:), II(:,:)
real*8, allocatable ::  quad_ec_x(:), quad_ec_w(:), quad_nk_x(:), quad_nk_w(:)
integer i,j,l,m,quad_order,io,jx
integer quad_order_ec, quad_order_lambda, quad_order_nk

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

allocate(quad_ec_w(quad_order)); allocate(quad_ec_x(quad_order_ec))
allocate(quad_w(quad_order)); allocate(quad_x(quad_order))
allocate(quad_w_l(quad_order_lambda)); allocate(quad_x_l(quad_order_lambda))


allocate(II(quad_order,quad_order)); allocate(KpW(quad_order,quad_order))
allocate(IK(quad_order,quad_order))

pi = acos(-1.d0); jj = cmplx(0.d0,1.d0)

mu = 1.d0

!write(*,*) 'Enter alpha and beta'
!read(*,*) alpha, beta

alpha = 1d-16
beta = 1.

open(20,file='./data/rho_x_exact_a0_b1.dat')
open(22,file='./data/rho_x_exact_quadrature_a0_b1.dat')

g = 2.d0/alpha

!print*, x_min, x_max


open(10, file='quadrature_w.txt')
open(11, file='quadrature_x.txt')
open(12, file='quad_lambda_w.txt')
open(13, file='quad_lambda_x.txt')
open(14, file='quad_ec_w.txt')
open(15, file='quad_ec_x.txt')

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

do i=1,quad_order_ec
read(14,*) quad_ec_w(i)
read(15,*) quad_ec_x(i)
enddo


close(10); close(11)
close(12); close(13)
close(14); close(15)

print*, "DO U WANT QUADRATURE OR NOT YOUR MOTHERFUCKER?"
read*, answer

if (answer.eq.'no') then 
write(*,*) 'Enter x0 and xf'
read(*,*) x0, xf
x = x0
do while (x.le.xf)
if (x .lt. 1.d0) dx = 1d-1
if (x .ge. 1.d0) dx = 1d-1

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
if (x .eq. 0.) norm = rho
!norm = (2.44198660696074744E-002,-6.82034317324653596E-022) !for quadrature 201, -5, 5
!norm =   (2.44950023737424728E-002,-1.24036027420767100E-020) !for quadrature 101, -5, 5
!norm = (5.91572988137469783E-003,2.16613318451298365E-022) ! for quadrature 101, -10, 10
!print*, norm
!norm = 1.
print*,   x, rho/norm
          write(20,*) x, real(rho/norm), aimag(rho/norm), abs(rho/norm)
          x=x+dx
enddo
endif









if(answer.eq.'yes') then
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



do jx=0,quad_order_nk

if(jx.eq.0) x = 0.
if(jx.ne.0) x = quad_nk_x(jx); w_x = quad_nk_w(jx)

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
!norm = (2.44198660696074744E-002,-6.82034317324653596E-022) !for quadrature 201, -5, 5
!norm =   (2.44950023737424728E-002,-1.24036027420767100E-020) !for quadrature 101, -5, 5
!norm = (5.91572988137469783E-003,2.16613318451298365E-022) ! for quadrature 101, -10, 10
!print*, norm
!norm = 1.
if(x.eq.0.) norm = rho
print*, x, rho/norm
if(jx.ne.0) write(22,*) x, w_x, rho/norm !!, real(rho/norm), aimag(rho/norm), abs(rho/norm)
enddo
endif



close(20);close(22)

!100 continue 

!print*, rho0

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
