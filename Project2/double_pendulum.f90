
module setup

    use numtype
    implicit none

    integer, parameter :: n_eq = 4
    real(dp), parameter :: gravity = 9.8_dp, l_1 = 0.5_dp, &
        mass_1 = 3.0_dp, l_2 = 0.4_dp, mass_2 = 2, &
            alpha_1 = pi/4, alpha_2 = pi/4

end module setup


program double_pendulum 

    use setup
    implicit none

    real(dp) :: t, dt,  tmax, eps, delt 
    real(dp), dimension(n_eq) :: y

    eps = .01
    delt = .01
    t = 0
    dt = .01                                    !time step
    tmax = 200                                  !200 second observation

    y(1:2)   =  (/ alpha_1, 0._dp /)            !initial (theta 1, omega 1) 
    y(3:4)   =  (/ alpha_1, 0._dp /)            !initial (theta 2, omega 2))              

    do while ( t < tmax )

        write(39,*) l_1 * sin(y(1)), -l_1 * cos(y(1))
        write(40,*) l_1 * sin(y(1)) + l_2 * sin(y(3)), -l_1 * cos(y(1)) - l_2 * cos(y(3))
        write(71,*) t, dt

        if ((y(1) - eps) < 0 .AND. y(2) > 0) then
             write(41,*) y(3), y(4)
        endif

        call rk4step ( t, dt, y )

    end do


end program double_pendulum


