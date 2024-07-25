program exact
implicit none
complex*16 jj, phi, nu
real*8 pi,delta, q, dq, alpha, beta, mu, lambda, fd
integer i,j

pi=acos(-1.d0); jj = cmplx(0.d0,1.d0)

open(22,file='./data/fig1_winding.dat')
open(33,file='./data/fig1_nonwinding.dat')

alpha = 1.
beta = 1.
mu = 1.

!winding region
lambda = 0.8

q = -10.; dq = 1d-2
do while(q.le.10.)
delta = pi/2. - atan(lambda -alpha*q)
fd = 1./(1.+exp(beta*(q**2. -mu)))
phi = 1. + fd*(exp(2.*jj*delta)-1.)
nu = log(phi)/(2.*pi*jj)

if(q.ge.0.8) nu = nu + 1.

write(22,*) real(phi), aimag(phi), q, real(nu), aimag(nu)
q = q+dq
enddo

!non-winding region
lambda = 1.2

q = -10.; dq = 1d-2
do while(q.le.10.)

delta = pi/2. - atan(lambda- alpha*q)
fd = 1./(1.+exp(beta*(q**2. -mu)))
phi = 1. + fd*(exp(2.*jj*delta)-1.)
nu = log(phi)/(2.*pi*jj)

write(33,*) real(phi), aimag(phi), q, real(nu), aimag(nu)
q = q+dq
enddo

end program
