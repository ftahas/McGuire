program exact
implicit none
character (len = 6) answer
complex*16 jj,ep1,ep2,em1,em2,KK,WW,rho,rho_l
complex*16, allocatable :: KpW(:,:), IK(:,:)
real*8 pi,delta1,delta2,beta,fd1,fd2,x0,xf
real*8 w1,w2,q1,q2,lambda,dlambda,mu,Z,alpha 
real*8 x,dx,k,dk,ec,fd,w,w_l,w_x,dq1,dq2
real*8, allocatable :: quad_w(:), quad_x(:), II(:,:)
real*8, allocatable ::quad_w_l(:), quad_x_l(:), quad_w_2(:), quad_x_2(:)
integer i,j,l,m,quad_order,io,jx, quad_order_2,quad_order_lambda

pi=acos(-1.d0); jj = cmplx(0.d0,1.d0)

mu = 1.

alpha = 0.1
beta = 10.

open(20,file='./data/rho_lambda_exact.dat')

!print*, x_min, x_max


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

allocate(quad_w_l(quad_order_lambda)); allocate(quad_x_l(quad_order_lambda))

allocate(quad_w(quad_order)); allocate(quad_x(quad_order))
allocate(II(quad_order,quad_order)); allocate(KpW(quad_order,quad_order))
allocate(IK(quad_order,quad_order))

open(10, file='quadrature_w.txt')
open(11, file='quadrature_x.txt')

open(12, file='quad_lambda_w.txt')
open(13, file='quad_lambda_x.txt')

do i=1,quad_order_lambda
read(12,*) quad_w_l(i)
read(13,*) quad_x_l(i)
enddo

do i=1,quad_order
read(10,*) quad_w(i)
read(11,*) quad_x(i)
do j=1,quad_order
if (i.eq.j) II(i,j)=1.
if (i.ne.j) II(i,j)=0.
enddo 
enddo 

!print*, quad_w, quad_x

close(10); close(11)
close(12);close(13)

!print*, II

!go to 100

x = 0.

do l=1, quad_order_lambda
lambda = quad_x_l(l)

rho_l=0.

do i=1, quad_order
do j=1, quad_order

q1 = quad_x(i); q2 = quad_x(j)
w1 = quad_w(i); w2 = quad_w(j)

delta1 = pi/2. - atan(lambda-alpha*q1); delta2 = pi/2. - atan(lambda-alpha*q2)

ep1 = (1./pi) * exp(jj*q1*x/2. +jj*delta1); ep2 = (1./pi) * exp(jj*q2*x/2. +jj*delta2)
em1 = exp(-jj*q1*x/2.)*sin(delta1); em2 = exp(-jj*q2*x/2.)*sin(delta2)

if (q1.ne.q2) then
                KK = 1./(q1-q2) * (ep1*em2-em1*ep2)
        else
                KK = exp(jj*atan(alpha*q2 -lambda))*(jj*alpha + x*(jj -q2*alpha + lambda)) / &
                        (pi*(-jj +q2*alpha -lambda)*sqrt(1. +(lambda -q2*alpha)**2.))
        end if
WW = 1./(pi) * em1*em2

fd1 = 1./(1.+exp(beta*(q1**2. -mu))); fd2 = 1./(1.+exp(beta*(q2**2. -mu)))

KpW(i,j) = II(i,j) + sqrt(w1*fd1)*(KK+WW)*sqrt(w2*fd2) 
IK(i,j) = II(i,j) + sqrt(w1*fd1)*KK*sqrt(w2*fd2)
enddo
enddo
rho_l = det(KpW)-det(IK)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!calc of ec(lambda)!!!!!!!!!!!!!!!!!!!!!!!!!!
ec = 0.d0
do m=1,quad_order
k = quad_x(m); w = quad_w(m)
fd = 1./(1.+exp(beta*(k**2. -mu)))
ec = ec -2.*(k/pi)*fd*(pi/2. - atan(lambda -alpha*k))*w
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


rho = exp(-beta*ec)*rho_l
print*,   lambda, rho
          write(20,*) lambda, real(rho)
enddo

close(20)

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