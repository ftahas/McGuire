program teste
real*8  m(3,3)
integer i,j

m(1,1)=1.; m(1,2) = 2.; m(1,3) = 1.
m(2,1)=0.3; m(2,2) = 0.2; m(2,3) = 0.1
m(3,1)=1.; m(3,2) = 2.; m(3,3) = 2.

print*, det(m)


contains 
        real*8 function DET(aa)
real*8 aa(:,:)
real*8 tmp,c(size(aa,dim=1),size(aa,dim=2))
real*8 max
        integer i,j,k,l,m,num(size(aa,dim=1))
        n=size(aa,dim=1)
        det=1.
        do k=1,n
                max=aa(k,k);num(k)=k;
                do i=k+1,n
                        if(abs(max)<abs(aa(i,k))) then
                                max=aa(i,k)
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
        end program
