module setup

    use numtype
    implicit none

    integer, parameter :: n_eq = 2
     

end module setup


program prob4

    use setup
    implicit none

    real(dp), dimension(n_eq) :: phi
    real(dp) :: r, dr, rmax

    r = 0
    dr = 0.01_dp
    rmax = 6._dp 
    

    phi(1) = 0
    phi(2) = 1

    do while ( r < rmax )

        write(17,*) r, phi(1)
        call rk4step ( r, dr, phi )

    end do

end program prob4


subroutine rk4step ( x, h, y )

    use setup
    implicit none
    real(dp), intent(inout) :: x
    real(dp), intent(in) :: h
    real(dp), intent(inout), dimension(n_eq) :: y
    real(dp), dimension(n_eq) :: k1, k2, k3, k4, dy

    k1 = kv( x, h, y )
    k2 = kv( x+h/2, h, y + k1/2 )
    k3 = kv( x+h/2, h, y + k2/2 ) 
    k4 = kv( x+h , h,  y + k3 )  
    dy = ( k1 + 2*k2 + 2*k3 + k4 ) / 6

    x = x + h 
    y = y + dy

    contains

        function kv( x, dx, phi ) result(k)

            use setup
            implicit none
            real(dp), intent(in) :: x, dx
            real(dp), intent(in), dimension(n_eq) :: phi
            real(dp), dimension(n_eq) :: f, k
!----------------------------------------------------------
            f(1) =  phi(2)
            f(2) =  -2/x*phi(2) - 4*pi*rho(x)
!____________________________________________________________
            k = dx * f 

        end function kv 

        function rho(r)

            use setup
            real(dp) :: r, rho

            rho = 1._dp/(8._dp*pi)*exp(-r)

        end function rho

end subroutine rk4step 