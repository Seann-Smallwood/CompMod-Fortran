program HW_1

    use numtype
    implicit none
            real(qp) :: u, umin, umax, du, w
            integer :: imax, i, counter, qmin, qmax, k, q, n
            real(qp), dimension(2001,0:10) :: H, Psi
            
            counter = 1
            imax = 10
            umin = -10
            umax = 10
            du = .01
            qmin = -1000
            qmax = 1000
            
            u = umin
            do q = qmin, qmax
                
                H(counter,0) = 1
                H(counter,1) = 2._qp*u
                do i = 1, imax-1
                    H(counter,i+1) = 2._qp*u*H(counter,i)-2*i*H(counter,i-1)
                end do
                
                w = exp(-((u**2._dp)/2))
            
                do n = 0,10
                    Psi(counter,n) = (pi**(-(1._dp/4._dp))) * (((2**n)*(product((/(k,k=1,n)/))))**(-(1._dp/2._dp))) * H(counter,n)*w
                end do
                
                write(14,*) u, Psi(counter,0:10)
                counter = counter + 1
                u = u + du
            end do
         
    print *,'Problem 1 result :', prob1(10)
    
    
    
    contains
    
        function prob1(imax) result(C)
            implicit none
            integer :: imax, i
            real(dp), dimension(imax) :: coeff, y
            real(dp) :: C 
        
            coeff = 0               
            coeff(1) = 0._dp
            coeff(2) = 1._dp
        
            do i = 3, imax
                coeff(i) = coeff(i-1) + coeff(i-2)
            end do
            
            y = 0
            y(imax) = coeff(imax)
            do i = imax-1, 1, -1
                y(i) = 2**coeff(i) + (1/y(i+1))
            end do

            C = 1/y(1)
            

        end function prob1

end program HW_1