program prob1
    use numtype
    implicit none
    real(dp), external :: func1
    real(dp) :: a, b, res, eps
    integer :: nint

    b = 1
    a = 0
    nint = 10
    eps = 1.e-10_dp

    call rombint(a,b, func1, res, nint, eps)
    print *, res

end program prob1
    
subroutine rombint( a, b, func, res, n, eps )
    ! Romberg integration
    !   res = int_a^b f(x) dx
    !   eps  required accuracy
    !   n    n_max in approx (input)
    !        accuaracy reached after n steps (output)
        use numtype
        implicit none
        integer, parameter :: maxint = 300
        real(dp) :: a, b, eps, res
        real(dp), external :: func
        integer ::  np, i, j, k, m, n
        real(dp) :: h, sumt, r(maxint,maxint)

        h = b - a 
        np = 1 
        r(1,1) = h/2 * ( func(a) + func(b) )
        res = r(1,1)

        do i=2,n 
            h = h/2
            np = 2*np 
            sumt = 0.0_dp
            do k=1,(np-1),2
                sumt = sumt + func( a + k*h)
            end do
            r(i,1) = 0.5_dp * r(i-1,1) + h * sumt
            m = 1
            do j=2,i 
                m = 4*m
                r(i,j) = r(i,j-1) + (r(i,j-1)-r(i-1,j-1))/(m-1)
            end do
            if ( abs(res-r(i,i)) < eps ) then
                n = i
                res = r(i,i)
                return
            end if
            res = r(i,i)
        end do
        !print *,' romint :',eps,r(i-1,i-1),res

end subroutine rombint
        
        
       
    

        
function func1(x)
    use numtype
    implicit none
    real(dp) :: x, func1
        
    func1 = sqrt(x)
        
end function func1
