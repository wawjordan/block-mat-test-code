program band_test

  use set_precision, only : prec
  use set_constants, only : set_derived_constants, zero, one
  use matrix_operations, only : banmul, bandec, banbks

  implicit none

  integer :: M1, M2, MP, MPL, N, NP
  integer, allocatable :: indx(:)
  real(prec), allocatable :: Aold(:,:), A(:,:), AL(:,:), b(:), x(:)
  real(prec) :: d

  character(LEN=30) :: rowfmt
  integer :: i, j


  M1  = 2
  M2  = 1
  MP  = 7
  MPL = 7
  N   = 7
  NP  = 7



  write(rowfmt,'(A,I4,A)') '(',N,'(1X,F10.4))'

  allocate( indx(N), Aold(N,M1+M2+1), A(N,M1+M2+1), AL(N,M1), b(N), x(n) )

  A = reshape(                                                           &
  (/  0.0_prec,0.0_prec,9.0_prec,3.0_prec,7.0_prec,3.0_prec,2.0_prec,    &
      0.0_prec,4.0_prec,2.0_prec,5.0_prec,9.0_prec,8.0_prec,4.0_prec,    &
      3.0_prec,1.0_prec,6.0_prec,8.0_prec,3.0_prec,4.0_prec,4.0_prec,    &
      1.0_prec,5.0_prec,5.0_prec,9.0_prec,2.0_prec,6.0_prec,0.0_prec /), &
      shape(A) )

  Aold = A

  b = one

  write(*,*)
  write(*,'(A3)') 'A ='
  do i = 1,N
    write(*,FMT=rowfmt) (A(i,j), j = 1,M1+M2+1)
  end do
  write(*,*)


  write(*,'(A3)') 'b ='
  do i = 1,N
    write(*,'(F10.4)') b(i)
  end do
  write(*,*)

  call bandec( A, N, M1, M2, NP, MP, AL, MPL, indx, d )

  write(*,'(A3)') 'U ='
  do i = 1,N
    write(*,FMT=rowfmt) (A(i,j), j = 1,M1+M2+1)
  end do
  write(*,*)
  write(*,'(A3)') 'L ='
  do i = 1,N
    write(*,FMT=rowfmt) (AL(i,j), j = 1,M1)
  end do
  write(*,*)

  call banbks( A, N, M1, M2, NP, MP, AL, MPL, indx, b )

  write(*,'(A3)') 'x ='
  do i = 1,N
    write(*,'(F10.4)') b(i)
  end do
  write(*,*)

  x = b

  call banmul( Aold, N, M1, M2, NP, MP, x, b )

  write(*,'(A5)') 'A*x ='
  do i = 1,N
    write(*,'(F10.4)') b(i)
  end do
  write(*,*)

  deallocate( indx, A, b, x )

end program band_test
