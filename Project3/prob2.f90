program prob_2

    use numtype
    use thiele_approx
    implicit none
    
    integer, parameter :: np = 11
    integer :: inmax, in, i
    real(dp), dimension(1:np) :: sums, xx, yy
    real(dp) :: x

   
    inmax = 10**7
    sums = 0
    

    
    do in = 1,inmax                             !itereate sums for function chi(s)
        sums(1) = sums(1) + 1/(in**2.0_dp)      !chi(s) = sum from 1 to inf. of n^-s
        sums(2) = sums(2) + 1/(in**2.2_dp)      !at various values of s
        sums(3) = sums(3) + 1/(in**2.4_dp)
        sums(4) = sums(4) + 1/(in**2.6_dp)
        sums(5) = sums(5) + 1/(in**2.8_dp)
        sums(6) = sums(6) + 1/(in**3.0_dp)
        sums(7) = sums(7) + 1/(in**3.2_dp)
        sums(8) = sums(8) + 1/(in**3.4_dp)
        sums(9) = sums(9) + 1/(in**3.6_dp)
        sums(10) = sums(10) + 1/(in**3.8_dp)
        sums(11) = sums(11) + 1/(in**4.0_dp)
    end do
    
    xx(1) = 2.0_dp                              !chosen values for s to analyze sum
    xx(2) = 2.2_dp                              !to input into thiele cf approximation
    xx(3) = 2.4_dp                              
    xx(4) = 2.6_dp                              
    xx(5) = 2.8_dp                              
    xx(6) = 3.0_dp                              
    xx(7) = 3.2_dp                              
    xx(8) = 3.4_dp                              
    xx(9) = 3.6_dp                              
    xx(10) = 3.8_dp                              
    xx(11) = 4.0_dp       
    
    yy(1) = sums(1)                             !resulting sums of corresponding s values
    yy(2) = sums(2)                             !yy is a function of xx
    yy(3) = sums(3)
    yy(4) = sums(4)
    yy(5) = sums(5)
    yy(6) = sums(6)
    yy(7) = sums(7)
    yy(8) = sums(8)
    yy(9) = sums(9)
    yy(10) = sums(10)
    yy(11) = sums(11)

    write(10,*) xx, yy                          !make sure s values and partial sums are good

    call thiele_coef( np, xx, yy, an )          !generate continued fraction coefficients

    x = 1.0_dp  
        
    print *, x, thiele_cf (x, np, xx, an)       !use cf coeffs to evalutae the function at x

   

end program