program exact
implicit none 
real*8 contact,pi,beta,k,w,dbeta,lambda
real*8 mu,alpha,ec,fd,w_l, s0,s1,s2,norm
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

print*, quad_order

allocate(quad_w(quad_order)); allocate(quad_x(quad_order))

pi=acos(-1.d0)

mu = 1.
alpha = 10.

open(20,file='./data/contact_a10.dat')


open(10, file='quadrature_w.txt')
open(11, file='quadrature_x.txt')

do i=1,quad_order
read(10,*) quad_w(i)
read(11,*) quad_x(i)
enddo 

close(10); close(11)

beta=0. 
dbeta = 1d-2
do while(beta.le.10.)

contact=0.; norm=0.
do j=1,quad_order
lambda = quad_x(j)
w_l = quad_w(j)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!calc of ec(lambda)!!!!!!!!!!!!!!!!!!!!!!!!!!
ec = 0.d0
do l=1,quad_order
k = quad_x(l); w = quad_w(l)
fd = 1./(1.+exp(beta*(k**2. -mu)))
ec = ec + (k/pi)*fd*(pi/2. - atan(lambda -alpha*k))*w
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call calculasn(quad_order,0,alpha,beta,lambda,s0)
call calculasn(quad_order,1,alpha,beta,lambda,s1)
call calculasn(quad_order,2,alpha,beta,lambda,s2)
contact = contact +exp(-beta*ec)*(s0*s2 -s1**2.)*w_l
norm = norm +pi*exp(-beta*ec)*s0*w_l
enddo

print*, beta, contact/norm
write(20,*) beta, contact/norm
beta = beta+dbeta
enddo

end program










subroutine calculasn(quad_order,n,alpha,beta,lambda,sn)
        implicit none
        real*8 mu,pi,alpha,beta,lambda,sn,k,w,fd
        real*8, allocatable :: quad_w(:), quad_x(:)
        integer i, n, quad_order

allocate(quad_w(quad_order)); allocate(quad_x(quad_order))
pi=acos(-1.d0)

mu = 1.

open(10, file='quadrature_w.txt')
open(11, file='quadrature_x.txt')

do i=1,quad_order
read(10,*) quad_w(i)
read(11,*) quad_x(i)
enddo

close(10);close(11)

        sn = 0.
        do i = 1, quad_order
        k = quad_x(i); w = quad_w(i)
        fd = 1./(1.+exp(beta*(k**2. -mu)))
        sn = sn + k**n *fd/((alpha*k -lambda)**2. +1.)*w
        enddo
        sn = sn/pi
endsubroutine
