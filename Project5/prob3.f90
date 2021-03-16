module setup

    use numtype
    implicit none

    !real(dp), parameter :: 
    

    integer, parameter :: npmax = 100, npar = 2
    integer, parameter :: nspmin = 1, nspmax = 43
    real(dp) :: xx(1:npmax), yy(1:npmax)  
    integer ::  nsp, ical, iprint 

end module setup

program prob3

    use setup
    implicit none
    
    integer :: stat, i, itmin, itmax
    real(dp), external :: scatpot
    real(dp) :: xstart(npar), fstart, stepi, epsf
    
    
    open(unit=2, file='reflectiondata.txt')
    stat = 0

    i = 1
    do 
        read(2,*, iostat= stat) xx(i),yy(i) 
        if( stat /= 0 ) exit
        print *,i,xx(i),yy(i)
        i = i + 1
    end do
    nsp = i-1
    close(2)

    xstart(1:npar) = (/ 1 , 1 /)
       

    ical = 0
    iprint = 7
    fstart = scatpot (xstart)
    stepi= 0.05_dp
    epsf = 0.001_dp
    itmin = 100
    itmax = 1000

    iprint = 0
    call downhill(npar,scatpot,xstart,fstart,stepi,epsf,itmin,itmax)

    iprint = 17
    fstart = scatpot (xstart)
    print *, xstart(1:npar)

end program prob3



function scatpot (par) 

    use setup
    implicit none
    real(dp) :: scatpot
    real(dp) :: par(npar)
    real(dp) :: scale, temp, x, fi 
    integer :: i 

    ical = ical + 1
    scale = par(1); temp = par(2)
    
    scatpot = 0        

    do i = nspmin, nspmax

        x = xx(i) 
        !fi = scale * x**3 * 1/( exp( hc * x / ( k_b * temp) ) - 1 )

        !scatpot = scatpot +  ( yy(i) - fi )**2  * 1/sqrt(  2._dp + yy(i) )
 
    end do 
    scatpot = scatpot / abs(nspmax-nspmin)
    print '(i4,2x,2f12.2,3x,f20.4)',ical, par(1:npar), scatpot

    ! printing 
    if (  iprint /= 0 ) then

        do i = nspmin, nspmax

            x = xx(i) 
            !fi =  scale * x**3 * 1/( exp( hc * x / ( k_b * temp) ) - 1 )
            write( unit=iprint, fmt='(3f15.4)') xx(i), yy(i) 
            write( unit=iprint + 1, fmt='(3f15.4)') xx(i),  fi 

        end do 

    end if

end function scatpot
