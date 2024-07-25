program correlation
implicit none
complex*16 jj,ep1,ep2,em1,em2,KK,WW,rho, KK1, WW1
complex*16, allocatable :: KpW(:,:), IK(:,:)
real*8 pi,delta1,delta2,beta,fd1,fd2,b,a
real*8 w1,w2,q1,q2,lambda,dlambda,g,mu 
real*8 x,dx,k,dk,ec,fd,Z,alpha,xi
real*8, allocatable :: quad_w(:), quad_x(:), II(:,:)
integer i,j,l,quad_order,io, n_xi

quad_order=0
open(1,file='quadrature_x.txt')
do
read(1,*,iostat=io)
if(io.ne.0) exit
quad_order = quad_order + 1
enddo
close(1)

!print*, quad_order

allocate(quad_w(quad_order)); allocate(quad_x(quad_order))
allocate(II(quad_order,quad_order)); allocate(KpW(quad_order,quad_order))
allocate(IK(quad_order,quad_order))

pi=acos(-1.); jj = cmplx(0.,1.)

g=10.d0; beta=1.d0; mu = 1.22713d0
alpha = 2.d0/g

open(20,file='./data/A_lambda_48.dat')

!Z = 1.d0/(alpha*pi) *(atan(alpha-lambda)+atan(alpha-lambda))
Z = 1.d0

open(9,file='quadrature_r.txt')
read(9,*) a, b
close(9)


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

n_xi=0
open(2,file='./data/xi_lambda_48.dat')
do
read(2,*,iostat=io)
if(io.ne.0) exit
n_xi = n_xi + 1
enddo
close(2)


open(2,file='./data/xi_lambda_48.dat')

x = 0.
do l=1, n_xi 

read(2,*) lambda, xi

do i=1, quad_order
do j=1, quad_order

q1 = quad_x(i); q2 = quad_x(j)
w1 = quad_w(i); w2 = quad_w(j)

delta1 = pi/2. - atan(lambda-2.*q1/g); delta2 = pi/2. - atan(lambda-2.*q2/g)

ep1 =  zexp(jj*q1*x/2.+jj*delta1)/pi; ep2 =  zexp(jj*q2*x/2.+jj*delta2)/pi
em1 = zexp(-jj*q1*x/2.)*sin(delta1); em2 = zexp(-jj*q2*x/2.)*sin(delta2)


if (q1 .ne. q2) then
                KK = (ep1*em2-em1*ep2)/(q1-q2)
                KK1 = zexp(-jj*x*(q1+q2)/2.)/(sqrt(1.+(alpha*q1-lambda)**2.)*sqrt(1.+(alpha*q2-lambda)**2.)*pi) &
                        *(zexp(jj*q1*x)*(lambda+jj-alpha*q1) -zexp(jj*q2*x)*(lambda+jj-alpha*q2))/(q1-q2)
        elseif (q1 .eq. q2) then
                KK = zexp(jj*atan(q2*2./g -lambda))*(jj*2./g + x*(jj -q2*2./g + lambda)) / &
                        (pi*(-jj +q2*2./g -lambda)*sqrt(1. +(lambda -q2*2./g)**2.))
                 KK1 = ((lambda+jj-alpha*q1)*jj*x-alpha)/(pi*(1.+(alpha*q1-lambda)**2.))
endif


WW = 1./(Z*pi) * em1*em2
WW1 = zexp(-jj*x*(q1+q2)/2.)/(sqrt(1.+(alpha*q1-lambda)**2.)*sqrt(1.+(alpha*q2-lambda)**2.)*pi)

if(abs(KK1-KK).ge.1d-6) print*, abs(KK1-KK), q1-q2

!WW = WW1; KK = KK1

fd1 = 1./(1.+exp(beta*(q1**2. -mu))); fd2 = 1./(1.+exp(beta*(q2**2. -mu)))

KpW(i,j) = II(i,j) + sqrt(w1*fd1)*(KK+WW)*sqrt(w2*fd2) 
IK(i,j) = II(i,j) + sqrt(w1*fd1)*KK*sqrt(w2*fd2)
enddo
enddo
rho = det(KpW)-det(IK)


!print*, rho
write(20,*) lambda, abs(rho)
enddo

close(2)
close(20)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!DETERMINANT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
complex*16 function det(aa)
complex*16 aa(:,:), tmp, c(size(aa,dim=1),size(aa,dim=2)), maxx
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
enddo
enddo

do i=1,n
det=det*aa(i,i)
enddo
return
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program
