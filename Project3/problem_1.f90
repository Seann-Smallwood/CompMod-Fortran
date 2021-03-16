module setup

    use numtype
    implicit none
 
    integer, parameter :: npmax = 400, npar = 3
    integer, parameter :: nspmin = 1, nspmax = 250
    real(dp) :: xx(1:npmax), yy(1:npmax)
    integer ::  nsp, ical, iprint 

end module setup 

program problem_1

    use setup
    implicit none

    integer :: stat, i, itmin, itmax, in
    real(dp), external :: chi2
    real(dp) :: xstart(npar), fstart, stepi, epsf, dx
    
    
    xx(1) = 1.0_dp                                                      !----------------------------------------
    dx = .02_dp

    do in = 1, 250                                                      !This loop populates x and y vectors from 
                                                                        !given function
        yy(in) = -2.0_dp + (2 * ((1 - (0.5*(exp((-xx(in)/2)+2))))**2))

        write(11,*) xx(in), yy(in)                                      !used to graph the minimum of the given function

        xx(in+1) = xx(in) + dx
    end do                                                              !----------------------------------------

   
    xstart(1:npar) = (/ 2.0_dp , 2.5_dp, -2.0_dp /)
    xstart(1:npar) = (/ 0.5573_dp , 3.2226_dp, -2.0384_dp /)

    ical = 0
    iprint = 7
    fstart = chi2 (xstart)
    stepi= 0.05_dp
    epsf = 0.001_dp
    itmin = 100
    itmax = 1000

    iprint = 0
    call downhill(npar,chi2,xstart,fstart,stepi,epsf,itmin,itmax)

    iprint = 17
    fstart = chi2 (xstart)
    print *, xstart(1:npar)

end program problem_1

function chi2( par ) 

    use setup
    implicit none
    real(dp) :: chi2
    real(dp) :: par(npar)
    real(dp) :: K, mid, x, fi, height 
    integer :: i 

    ical = ical + 1
    K = par(1); mid = par(2); height = par(3)          !using harmonic function V(x) = 1/2 K (x-mid)^2 + height
    
    chi2 = 0        ! chi^2

    do i = nspmin, nspmax

        x = xx(i) 
        fi = 0.5_dp * K * (x - mid)**2 + height         !harmonic function V(x) = 1/2 K (x-mid)^2 + height

        chi2 = chi2 +  ( yy(i) - fi )**2  * 1/sqrt(  2._dp + yy(i) )
 
    end do 
    chi2 = chi2 / abs(nspmax-nspmin)
    print '(i4,2x,3f12.2,3x,f20.4)',ical, par(1:npar), chi2

    ! printing 
    if (  iprint /= 0 ) then

        do i = nspmin, nspmax

            x = xx(i) 
            fi =  0.5_dp * K * (x-mid)**2 + height
            write( unit=iprint, fmt='(3f15.4)') xx(i), yy(i) 
            write( unit=iprint + 1, fmt='(3f15.4)') xx(i),  fi 

        end do 

    end if

end function chi2
