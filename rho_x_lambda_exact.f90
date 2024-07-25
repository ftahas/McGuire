program correlation
implicit none
complex*16 jj,ep1,ep2,em1,em2,KK,WW,rho,KK1,WW1
complex*16, allocatable :: KpW(:,:), IK(:,:)
real*8 pi,delta1,delta2,beta,fd1,fd2,b,a
real*8 w1,w2,q1,q2,lambda,dlambda,mu 
real*8 x,dx,k,dk,ec,fd,Z,alpha,norm,x0, xf
real*8, allocatable :: quad_w(:), quad_x(:), II(:,:)
integer i,j,l,quad_order,io

pi=dacos(-1.d0); jj = cmplx(0.d0,1.d0)

alpha = 1d0 
beta = 1.d0
lambda = 1.2d0

mu = 1.d0

open(20,file='./data/rho_x_lambda12_a1_b1_exact.dat')


print*, 'alpha=', alpha
print*, 'beta=', beta
print*, 'lambda=', lambda


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


!x = 0.d0; dx = 1d-1

write(*,*) 'Enter x0 and xf'
read(*,*) x0, xf

x = x0; dx = 1d-1
do while (x.le.xf)
!if (x<0.5) dx = 0.01d0
!if (x>=1.d0) dx = 0.1d0
do i=1, quad_order
do j=1, quad_order

q1 = quad_x(i); q2 = quad_x(j)
w1 = quad_w(i); w2 = quad_w(j)

delta1 = pi/2.d0 - datan(lambda-alpha*q1); delta2 = pi/2.d0 - datan(lambda-alpha*q2)

ep1 = cdexp(jj*q1*x/2.d0 +jj*delta1)/pi; ep2 = cdexp(jj*q2*x/2.d0 +jj*delta2)/pi
em1 = cdexp(-jj*q1*x/2.d0)*dsin(delta1); em2 = cdexp(-jj*q2*x/2.d0)*dsin(delta2)



if (q1.ne.q2) then !dabs(q1-q2).gt.1d-6) then
        
        KK = (ep1*em2-em1*ep2)/(q1-q2)

                KK1 = zexp(-jj*x*(q1+q2)/2.d0)/(dsqrt(1.+(alpha*q1-lambda)**2.)*dsqrt(1.+(alpha*q2-lambda)**2.d0)*pi) &
                        *(zexp(jj*q1*x)*(lambda+jj-alpha*q1) -zexp(jj*q2*x)*(lambda+jj-alpha*q2))/(q1-q2)
        
        elseif (q1.eq.q2) then !dabs(q1-q2).le.1d-6) then
        
                KK = cdexp(jj*datan(q2*alpha -lambda))*(jj*alpha + x*(jj -q2*alpha + lambda))/ &
                        (pi*(-jj +q2*alpha -lambda)*dsqrt(1.d0 +(lambda-q2*alpha)**2.d0))

                 KK1 = ((lambda+jj-alpha*q1)*jj*x-alpha)/(pi*(1.d0+(alpha*q1-lambda)**2.d0))
endif

WW = em1*em2/(Z*pi)
WW1 = cdexp(-jj*x*(q1+q2)/2.)/(dsqrt(1.+(alpha*q1-lambda)**2.)*dsqrt(1.+(alpha*q2-lambda)**2.)*pi)

if(abs(KK-KK1).ge.1d-6) print*, abs(KK-KK1), abs(q1-q2)

!KK = KK1; WW = WW1

fd1 = 1./(1.+dexp(beta*(q1**2. -mu))); fd2 = 1./(1.+dexp(beta*(q2**2. -mu)))

KpW(i,j) = II(i,j) + dsqrt(w1*fd1)*(KK+WW)*dsqrt(w2*fd2) 
IK(i,j) = II(i,j) + dsqrt(w1*fd1)*KK*dsqrt(w2*fd2)

enddo
enddo
rho = det(KpW)-det(IK)

print*, x, rho
write(20,*) x, real(rho), aimag(rho), abs(rho)
x=x+dx
enddo

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
