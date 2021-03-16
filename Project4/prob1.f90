

program prob1

    use numtype
    implicit none

    integer, parameter :: ndim = 6, lwork = 4*ndim          !2*j + 1 (spin states)
    complex(dp), dimension(ndim,ndim) :: j_plus, j_minus, j_x, j_y, j_z, j, j_2, sum_j_x, sum_j_y, sum_j_z,sum_j_xn, sum_j_yn, lhs
    complex(dp) ::  i, work(lwork)
    real(dp), dimension(ndim) :: w_x, w_y, w_2
    real(dp) :: jj, mm, kplus, kminus, rwork(lwork), alpha, in
    integer :: n, info, ii

    i = (0._dp,1._dp)
    jj = 5._dp/2._dp

    mm = -jj
   
    do n = 1, ndim                          !loop to generate J_z
        J_z(n,n) = cmplx(mm,0._dp)
        mm = mm + 1._dp
        
    end do

    j_plus = 0._dp                              !initialize matrices
    j_minus = 0._dp
    j_x = 0._dp
    j_y = 0._dp
    j_2 = 0._dp
    sum_j_x = 0._dp
    sum_j_y = 0._dp
    sum_j_z = 0._dp
    sum_j_xn = 0._dp
    sum_j_yn = 0._dp

    mm = -jj
    
    if ((jj*(jj+1._dp))-(mm*(mm+1._dp)) >= 0._dp) then  !first j_plus
            
        kplus = sqrt((jj*(jj+1._dp))-(mm*(mm+1._dp)))
        j_plus(2,1) = cmplx(kplus, 0._dp) 
    
    else
       
        kplus = sqrt(-1._dp*((jj*(jj+1._dp))-(mm*(mm+1._dp))))
        j_plus(2,1) = cmplx(0._dp, kplus)
    end if 

    
    mm = mm + 1._dp
    
    Do n = 2, ndim - 1                          !loop to generate middle jplus/minus
        
        if ((jj*(jj+1._dp))-(mm*(mm+1._dp)) >= 0._dp) then
            
            kplus = sqrt((jj*(jj+1._dp))-(mm*(mm+1._dp)))
            j_plus(n+1,n) = cmplx(kplus, 0._dp) 
        
        else
           
            kplus = sqrt(-((jj*(jj+1._dp))-(mm*(mm+1._dp))))
            j_plus(n+1,n) = cmplx(0._dp, kplus)
        end if     
        
        if ((jj*(jj+1._dp))-(mm*(mm-1._dp)) >= 0._dp) then
            
            kminus = sqrt((jj*(jj+1._dp))-(mm*(mm-1._dp)))
            j_minus(n-1,n) = cmplx(kminus, 0._dp) 
        
        else
            
            kminus = sqrt(-((jj*(jj+1._dp))-(mm*(mm-1._dp))))
            j_minus(n-1,n) = cmplx(0._dp, kminus)
        end if

        mm = mm + 1._dp
    end do

    if ((jj*(jj+1._dp))-(mm*(mm-1._dp)) >= 0._dp) then !last j_min
            
        kminus = sqrt((jj*(jj+1._dp))-(mm*(mm-1._dp)))
        j_minus(ndim-1,ndim) = cmplx(kminus, 0._dp) 
    
    else
        
        kminus = sqrt(-1._dp*((jj*(jj+1._dp))-(mm*(mm-1._dp))))
        j_minus(ndim-1,ndim) = cmplx(0._dp, kminus)
    end if
   

    j_x = (j_plus + j_minus)/2.0_dp
    j_y = (j_plus - j_minus)/(2.0_dp*i)

    j_2 = (j_x ** 2) + (j_y ** 2) + (j_z ** 2)

!--------------------------------------------------------
   
    alpha = pi/20
    in = 1
    do n = 1, 20

    sum_j_x = sum_j_x + ((i*j_x*alpha)**in)/product((/(ii,ii=1,n)/))
    sum_j_y = sum_j_y + ((i*j_y*alpha)**in)/product((/(ii,ii=1,n)/))

    sum_j_xn = sum_j_xn + ((-i*j_x*alpha)**in)/product((/(ii,ii=1,n)/))
    sum_j_yn = sum_j_yn + ((-i*j_y*alpha)**in)/product((/(ii,ii=1,n)/))

    sum_j_z = sum_j_z + ((i*j_y*(alpha**2))**in)/product((/(ii,ii=1,n)/))

    in = in + 1
    end do

    lhs = sum_j_y*sum_j_x*sum_j_yn*sum_j_xn
 
    do n = 1, ndim
        !print '(36f7.2)',j_plus(n,1:ndim)
        !print '(36f7.2)',j_minus(n,1:ndim)
        !print '(36f7.2)', j_x(1:ndim,n)
        !print '(36f7.2)', j_y(1:ndim,n)
        !print '(36f7.2)', j_z(1:ndim,n)
        print '(12f7.2)', j_2(1:ndim,n)
        !print '(12f7.2)', lhs(1:ndim,n) 
    end do

    print *, '-----------------------'
    
    info = 0
    call zheev	( 'v', 'u', ndim , j_x, ndim, w_x, work, lwork, rwork, info )
    if( info /= 0 ) stop ' info /= 0 '
    !print *,w_x(1:ndim)

    info = 0
    call zheev	( 'v', 'u', ndim , j_y, ndim, w_y, work, lwork, rwork, info )
    if( info /= 0 ) stop ' info /= 0 '
    !print *,w_y(1:ndim)

    info = 0
    call zheev	( 'v', 'u', ndim , j_2, ndim, w_2, work, lwork, rwork, info )
    if( info /= 0 ) stop ' info /= 0 '
    !print *,w_2(1:ndim)

    do n = 1, ndim
        !print '(10f7.2)', w_x(n)      !eigenvalues/vectors for j_x
        !print '(36f7.2)', j_x(1:ndim,n) 
        
        !print '(10f7.2)', w_y(n)      !eigenvalues/vectors for j_y
        !print '(36f7.2)', j_y(1:ndim,n)
    
        print '(10f7.2)', w_2(n)      !eigenvalues/vectors for j_y
        print '(12f7.2)', j_2(1:ndim,n)
    end do
    
    do n = 1, ndim
        !print '(36f7.2)',j_plus(n,1:ndim)
        !print '(36f7.2)',j_minus(n,1:ndim)
        !print '(36f7.2)', j_x(1:ndim,n)
        !print '(36f7.2)', j_y(1:ndim,n)
        !print '(36f7.2)', j_z(1:ndim,n)
        !print '(36f7.2)', j_2(1:ndim,n)
        !print '(12f7.2)', sum_j_z(1:ndim,n) 
    end do

end program prob1

