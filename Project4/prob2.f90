

program prob2

    use numtype
    implicit none

    integer, parameter :: ndim = 30, lwork = 4*ndim, ldvr = 1, ldvl = 1      !dimensions + 1 for t^0 term
    real(dp), dimension(ndim,ndim) :: comp_p
    real(dp), dimension(ndim) :: wr, wi
    real(dp) :: nn, nmax, work(lwork), vl(ldvl,ndim), vr(ldvr,ndim)
    integer :: n, info

    nmax = 31
    nn = 0
    
    comp_p = 0._dp
    
    do n = 1, ndim - 1

        comp_p(n+1,n) = 1._dp

    end do
    
    do n = 1, ndim
        
        comp_p(n,ndim) = -(nmax - nn)

        nn = nn + 1
        !print '(21f7.2)',comp_p(n,1:ndim)

    end do
    print *, '--------------------------' 
    info = 0
    call dgeev ('n','n',ndim,comp_p,ndim,wr, wi, vl, ldvl, vr, ldvr, work,lwork,info)

    print '(30f7.2)', wr(1:ndim)
    print '(30f7.2)', wi(1:ndim)
end program prob2

