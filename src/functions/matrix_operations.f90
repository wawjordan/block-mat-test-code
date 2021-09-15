module matrix_operations

  use set_precision, only : prec, dp
  use set_constants, only : zero, one

  implicit none

  public :: ludcmp, lubksb, banmul, bandec, banbks, mprove

contains

  subroutine ludcmp( A, N, NP, indx, d )
    use set_constants, only : near_zero
    integer, intent(in) :: N, NP
    integer, intent(out) :: indx(N)
    integer, parameter :: NMAX = 500
    real(prec) :: d, A(NP,NP)
    integer :: i, imax, j, k
    real(prec) :: Aamax, dum, runsum, vv(NMAX)
    d = one                   ! no row interchanges yet
    do i = 1,N                ! loop over rows to get the implicit scaling
      Aamax = zero
      do j = 1,N
        if (abs(A(i,j))>Aamax) then
          Aamax = abs(A(i,j))
        end if
      end do
      if (Aamax<near_zero) then  ! no nonzero largest element
        write(*,*) 'singular matrix in ludcmp'
        stop
      end if
      vv(i) = one/Aamax     ! save the scaling
    end do
    do j = 1,N              ! loop over columns of Crout's method
      do i = 1,j-1
        runsum = A(i,j)
        do k = 1,i-1
          runsum = runsum - A(i,k)*A(k,j)
        end do
        A(i,j) = runsum
      end do
      Aamax = zero         ! initialize search for largest pivot element
      do i = j,N
        runsum = A(i,j)
        do k = 1,j-1
          runsum = runsum - A(i,k)*A(k,j)
        end do
        A(i,j) = runsum
        dum = vv(i)*abs(runsum) ! figure of merit of the pivot
        if (dum>=Aamax) then    ! better than the best so far?
          imax = i
          Aamax = dum
        end if
      end do
      if (j/=imax) then   ! do we need to interchange rows
        do k = 1,N
          dum = A(imax,k)
          A(imax,k) = A(j,k)
          A(j,k) = dum
        end do
        d = -d            ! change parity of d
        vv(imax) = vv(j)  ! also interchange the scale factor
      end if
      indx(j) = imax
      if (abs(A(j,j))<near_zero) then ! if the matrix is nearly singular
        A(j,j) = near_zero
      end if
      if (j/=N) then     ! divide by the pivot element
        dum = one/A(j,j)
        do i = j+1,N
          A(i,j) = A(i,j)*dum
        end do
      end if
    end do              ! Next column in the reduction
  end subroutine ludcmp

  subroutine lubksb( A, N, NP, indx, b )
    use set_constants, only : small
    integer, intent(in) :: N, NP, indx(N)
    real(prec), intent(in) :: A(NP,NP)
    real(prec), intent(inout) :: b(N)
    integer :: i, ii, j, LL
    real(prec) :: runsum
    ii = 0            ! when ii is positive, it is the index of the first
    do i = 1,N        ! non-vanishing element of b. Loop does forward
      LL = indx(i)    ! substitution, unscrambling along the way
      runsum = b(LL)
      b(LL) = b(i)
      if (ii/=0) then
        do j = ii,i-1
          runsum = runsum - A(i,j)*b(j)
        end do
      elseif (abs(runsum)>small) then
        ii = i     ! a nonzero element was encountered
      end if
      b(i) = runsum
    end do
    do i = n,1,-1  ! do backsubstitution
      runsum = b(i)
      do j = i+1,N
        runsum = runsum - A(i,j)*b(j)
      end do
      b(i) = runsum/A(i,i)  ! store a component of the solution vector
    end do
  end subroutine lubksb

  subroutine banmul( A, N, M1, M2, NP, MP, x, b )
    integer, intent(in) :: M1, M2, MP, N, NP
    real(prec), intent(in)  :: A(NP,MP), x(N)
    real(prec), intent(out) :: b(N)
    integer :: i, j, k

    do i = 1,N
      b(i) = zero
      k = i - M1 - 1
      do j = max(1,1-k),min(M1+M2+1,N-k)
        b(i) = b(i) + A(i,j)*x(j+k)
      end do
    end do
  end subroutine banmul

  subroutine bandec( A, N, M1, M2, NP, MP, AL, MPL, indx, d )
    use set_constants, only : small
    integer, intent(in)       :: M1, M2, MP, MPL, N, NP
    integer, intent(inout)    :: indx(N)
    real(prec), intent(inout) :: d, A(NP,MP), AL(NP,MPL)
    integer :: i, j, k, L, mm
    real(prec) :: dum

    MM = M1 + M2 + 1
    if ( (MM>MP).or.(M1>MPL).or.(N>NP) ) then
      write(*,*) 'Bad args in bandec'
      stop
    end if
    L = M1
    do i = 1,M1            ! rearrange storage
      do j = M1+2-i,MM
        A(i,j-L) = A(i,j)
      end do
      L = L-1
      do j = MM-L,MM
        A(i,j) = zero
      end do
    end do
    d = one
    L = M1
    do k = 1,n             ! for each row
      dum = A(k,1)
      i = k
      if (L<N) then
        L = L + 1
      end if
      do j = k+1,L        ! find the pivot element
        if (abs(A(j,1))>abs(dum)) then
          dum = A(j,1)
          i = j
        end if
      end do
      indx(k) = i
      if (abs(dum)<small) then
        A(k,1) = small    ! matrix is algorithmically singular
      end if
      if (i/=k) then      ! interchange rows
        d = -d
        do j = 1,MM
          dum = A(k,j)
          A(k,j) = A(i,j)
          A(i,j) = dum
        end do
      end if
      do i = k+1,L       ! do the elimination
        dum = A(i,1)/A(k,1)
        AL(k,i-k) = dum
        do j = 2,MM
          A(i,j-1) = A(i,j) - dum*A(k,j)
        end do
        A(i,MM) = zero
      end do
    end do
  end subroutine bandec

  subroutine banbks( A, N, M1, M2, NP, MP, AL, MPL, indx, b )
    integer, intent(in)       :: M1, M2, MP, MPL, N, NP
    integer, intent(inout)    :: indx(N)
    real(prec), intent(inout) :: A(NP,MP), AL(NP,MPL), b(N)
    integer    :: i, k, L, MM
    real(prec) :: dum

    MM = M1 + M2 + 1
    if ( (MM>MP).or.(M1>MPL).or.(N>NP) ) then
      write(*,*) 'Bad args in banbks'
      stop
    end if
    L = M1
    do k = 1,N        ! forward substitution, unscrambling along the way
      i = indx(k)
      if (i/=k) then
        dum = b(k)
        b(k) = b(i)
        b(i) = dum
      end if
      if (L<N) then
        L = L + 1
      end if
      do i = k+1,L
        b(i) = b(i) - AL(k,i-k)*b(k)
      end do
    end do
    L = 1
    do i = N,1,-1    ! backsubstitution
      dum = b(i)
      do k = 2,L
        dum = dum - A(i,k)*b(k+i-1)
      end do
      b(i) = dum/A(i,1)
      if (L<MM) then
        L = L + 1
      end if
    end do
  end subroutine banbks

  subroutine mprove( A, Alud, N, Np, indx, b, x )
    integer, intent(in) :: N, NP, indx(N)
    integer, parameter :: NMAX = 500
    real(prec), intent(in) :: A(NP,NP), Alud(NP,NP), b(N)
    real(prec), intent(inout) :: x(N)
    integer :: i, j
    real(prec) :: r(NMAX)
    real(prec) :: sdp

    do i = 1,N                 ! calculate RHS, accumulating residual
      sdp = -b(i)
      do j = 1,N
        sdp = sdp + A(i,j)*x(j)
      end do
      r(i) = sdp
    end do
    call lubksb( Alud, N, NP, indx, r ) ! solve for error term
    do i = 1,N                 ! subtract from the old solution
      x(i) = x(i) - r(i)
    end do
  end subroutine mprove

end module matrix_operations
