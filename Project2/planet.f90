
module setup

    use numtype
    implicit none

    integer, parameter :: n_eq = 12
    real(dp), parameter :: gravity = 6.673e-11_dp, &
        mass_sun = 1.981e+30_dp, mass_earth = 6.1e+24_dp, mass_comet = 5e+14_dp

end module setup


program planet 

    use setup
    implicit none

    real(dp) :: t, dt,  tmax 
    real(dp), dimension(n_eq) :: y

    t = 0
    dt = 7 * 24 * 60*60                                 !weekly time steps
    tmax =  10 * 60*60*24*365                           !10 year observational period

    y(1:3)   =  (/ 1.496e+11_dp, 0._dp, 0._dp /)        !initial (x,y,z) of earth
    y(4:6)   =  (/ 0._dp, 29.726e+3_dp, 0._dp /)  !initial compontent velocities of earth                 

    ! addition of satelitte
    
    y(7:9)   = (/1.496e+11_dp * cos(pi/3), 1.496e+11_dp * sin(pi/3), 0._dp /)           !initial (x,y,z) of satelitte
    y(10:12)   =  (/ -29.726e+3 * sin(pi/3),29.726e+3 * cos(pi/3) , 0._dp /)        !initial compontent velocities of satelitte

    do while ( t < tmax )

        write(39,*) y(1), y(2)
        write(40,*) y(7), y(8)
        write(71,*) t, dt

        call rk4step ( t, dt, y )

    end do


end program planet


