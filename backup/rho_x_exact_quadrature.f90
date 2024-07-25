program correlation
implicit none
complex*16 jj,ep1,ep2,em1,em2,KK,WW,rho,norm 
complex*16, allocatable :: KpW(:,:), IK(:,:)
real*8 pi,delta1,delta2,beta,fd1,fd2,x0,xf
real*8 w1,w2,q1,q2,lambda,dlambda,g,mu,Z,alpha 
real*8 x,dx,k,dk,ec,fd,w,w_l,w_x
real*8, allocatable :: quad_w(:), quad_x(:), II(:,:)
integer i,j,l,m,quad_order,io,jx

quad_order=0
open(1,file='quadrature_x.txt')
do 
read(1,*,iostat=io)
if(io.ne.0) exit
quad_order = quad_order + 1
enddo
close(1)

print*, quad_order

allocate(quad_w(quad_order)); allocate(quad_x(quad_order))
allocate(II(quad_order,quad_order)); allocate(KpW(quad_order,quad_order))
allocate(IK(quad_order,quad_order))

pi=acos(-1.d0); jj = cmplx(0.d0,1.d0)

mu = 1.d0

!write(*,*) 'Enter alpha and beta'
!read(*,*) alpha, beta

alpha = 1.d0
beta = 1.d0

open(22,file='./data/rho_x_exact_quadrature.dat')

g = 2.d0/alpha

!print*, x_min, x_max

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

!print*, quad_w, quad_x

close(10); close(11)

!print*, II

!go to 100

do jx=1,quad_order
x = quad_x(jx); w_x = quad_w(jx)

rho=0.
do m=1, quad_order
lambda = quad_x(m)
w_l = quad_w(m)

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
rho = (det(KpW)-det(IK))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!calc of ec(lambda)!!!!!!!!!!!!!!!!!!!!!!!!!!
ec = 0.d0
do l=1,quad_order
k = quad_x(l); w = quad_w(l)
fd = 1./(1.+exp(beta*(k**2. -mu)))
ec = ec + (k/pi)*fd*(pi/2. - atan(lambda -2.*k/g))*w
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


rho = rho + exp(-beta*ec)*rho*w_l
!print*, rho
enddo
!norm = (2.44198660696074744E-002,-6.82034317324653596E-022) !for quadrature 201, -5, 5
norm =   (2.44950023737424728E-002,-1.24036027420767100E-020) !for quadrature 101, -5, 5
!norm = (5.91572988137469783E-003,2.16613318451298365E-022) ! for quadrature 101, -10, 10
!print*, norm
!norm = 1.
print*, x, rho/norm
          write(22,*) x, w_x, rho/norm, real(rho/norm), aimag(rho/norm), abs(rho/norm)
enddo
close(22)

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
