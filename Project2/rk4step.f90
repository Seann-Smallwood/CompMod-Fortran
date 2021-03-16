

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

        function kv(t, dt, y ) result(k)

            use setup 
	        implicit none
	        real(dp), intent(in) :: t, dt
	        real(dp), dimension(n_eq), intent(in) :: y
	        real(dp), dimension(n_eq) :: f, k

            !---------------------
			

            f(1) =  (-gravity*(2*mass_1 + mass_2)*sin(y(2))-mass_2*gravity*sin(y(2)-2*y(4))-2*sin(y(2) &
                -y(4))*mass_2*((y(3))**2*l_2 + (y(1))**2*l_1*cos(y(2) - y(4))))/ &
                    (l_1*(2*mass_1 + mass_2 - mass_2*cos(2 * y(2) - 2 * y(4))))
            
            f(2) = y(1)


            f(3) = (2*sin(y(2)-y(4))*((y(1))**2*l_1*(mass_1+mass_2)*cos(y(2)) + &
                (y(3))**2*l_2*mass_2*cos(y(2)-y(4))))/(l_2*(2*mass_1+mass_2-mass_2*cos(2*y(2)-2*y(4))))
            
            f(4) = y(3)
		

              
            !----------------

            k = dt * f

        end function kv 

end subroutine rk4step 





