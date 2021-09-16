program tridiag_test

  use set_precision, only : prec
  use set_constants, only : set_derived_constants, zero, one
  use tridiag_operations, only : trivec, tprec

  implicit none

  integer :: imax = 2
  integer :: neq  = 2

  integer :: j  = -99
  integer :: k  = -99
  integer :: m  = -99
  integer :: jj = -99
  real(prec) :: percerr
  real(prec), allocatable :: da(:,:,:), dd(:,:,:), dc(:,:,:), &
                                c(:,:), x(:,:), xexact(:,:)

  allocate( da(imax,neq,neq), &
            dd(imax,neq,neq), &
            dc(imax,neq,neq), &
            c(imax,neq),      &
            x(imax,neq),      &
            xexact(imax,neq) )

  jj = imax

  dd(1,1,1) = 10.0_prec
  dd(1,1,2) = 2.0_prec
  dd(1,2,1) = 3.0_prec
  dd(1,2,2) = 8.0_prec

  dd(2,1,1) = 7.0_prec
  dd(2,1,2) = 1.0_prec
  dd(2,2,1) = 4.0_prec
  dd(2,2,2) = 5.0_prec

  da(1,1,1) = 0.0_prec
  da(1,1,2) = 0.0_prec
  da(1,2,1) = 0.0_prec
  da(1,2,2) = 0.0_prec

  da(2,1,1) = 3.0_prec
  da(2,1,2) = 1.0_prec
  da(2,2,1) = 2.0_prec
  da(2,2,2) = 5.0_prec

  dc(1,1,1) = 1.0_prec
  dc(1,1,2) = 2.0_prec
  dc(1,2,1) = 3.0_prec
  dc(1,2,2) = 4.0_prec

  dc(2,1,1) = 0.0_prec
  dc(2,1,2) = 0.0_prec
  dc(2,2,1) = 0.0_prec
  dc(2,2,2) = 0.0_prec

  do j = 1, jj
    do k = 1, neq
      do m = 1, neq
        write(*,*) 'j,k,m: ',j,k,m
        write(*,*) '    da,dd,dc: ',da(j,k,m),dd(j,k,m),dc(j,k,m)
      enddo
    enddo
  enddo

  c(1,1) = 16.0_prec
  c(1,2) = 14.0_prec
  c(2,1) = 10.0_prec
  c(2,2) = 6.0_prec

! Exact Solution (via Matlab)
  xexact(1,1) = 1.52901023890785_prec
  xexact(1,2) = 1.80375426621160_prec
  xexact(2,1) = 0.77815699658703_prec
  xexact(2,2) = -1.83788395904437_prec

  do j = 1, jj
    do k = 1, neq
      write(*,*) '         j,k,c: ',j,k,c(j,k)
      x(j,k) = c(j,k)
    enddo
  enddo

  call trivec(jj,neq,da,dd,dc)
  call tprec(jj,neq,da,dd,dc,x)
  write(*,*)
  write(*,*) 'SOLUTION'
  write(*,*) '  j, k, x(j,k),           xexact(j,k),     %error'
  do j = 1, jj
    do k = 1, neq
      percerr = (x(j,k) - xexact(j,k))/xexact(j,k)*100.0_prec
      write(*,*) j,k,x(j,k),xexact(j,k),percerr
    enddo
  enddo

  deallocate( da, dd, dc, c, x, xexact )

end program tridiag_test
