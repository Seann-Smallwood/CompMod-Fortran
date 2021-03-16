
module setup

    use numtype
    implicit none

    real(dp), parameter :: mass = 1.0_dp, g = 9.8_dp, &
        spring_k = 10_dp, length_0 = 2.0_dp, tmax = 30_dp, &
        r_sus = 0.5_dp


end module setup


program mass_spring 

    use setup
    implicit none

    real(dp) :: t, dt, dx, dy, mag_l, stretch
    real(dp), dimension(3) :: bob_pos, bob_vel, bob_mom, F_spring, F_g, &
        F, delta, l, l_hat

    t = 0
    dt = .01 

    bob_pos = (/3.5_dp, 0._dp, 0._dp/)              !initial position of bob
    bob_vel = (/0._dp, 0._dp, 0._dp/)               !initial veloctiy of bob
    bob_mom = mass * bob_vel

    !F_g = (/0._dp, -g*mass, 0._dp/)



    do while (t < tmax)

        dx = r_sus * cos(7*t)                       !time parameterization
        dy = r_sus * sin(7*t)                       !of suspension point
        
        delta = (/dx, dy, 0._dp/)                   !position of sus point at time t
        
        l = bob_pos - delta
        mag_l = sqrt(l(1)**2 + l(2)**2 + l(3)**2)
        l_hat = l/mag_l
        stretch = mag_l - length_0

        F_spring = -spring_k * stretch * l_hat
        F = F_spring! + F_g 

        bob_mom = bob_mom + F * dt
        bob_vel = bob_mom/mass
        bob_pos = bob_pos + bob_vel * dt

        write(13,*) bob_pos(1), bob_pos(2)
        write(14,*) t

        t = t + dt
    end do
end program mass_spring


