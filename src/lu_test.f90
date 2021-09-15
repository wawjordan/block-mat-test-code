program lu_test

  use set_precision, only : prec
  use set_constants, only : set_derived_constants, zero, one
  use matrix_operations, only : ludcmp, lubksb, mprove

  implicit none

  integer :: N, NP
  integer, allocatable :: indx(:)
  real(prec), allocatable :: A(:,:), Alud(:,:), b(:), x(:), xold(:)
  real(prec) :: d

  character(LEN=30) :: rowfmt
  integer :: i, j

  N  = 5
  NP = N

  write(rowfmt,'(A,I4,A)') '(',N,'(1X,F12.4))'

  allocate( indx(N), A(N,N), Alud(N,N), b(N), x(N), xold(N) )

  A = transpose( reshape(                                     &
  (/ 17.0_prec, 24.0_prec,  1.0_prec,  8.0_prec, 15.0_prec,   &
     23.0_prec,  5.0_prec,  7.0_prec, 14.0_prec, 16.0_prec,   &
      4.0_prec,  6.0_prec, 13.0_prec, 20.0_prec, 22.0_prec,   &
     10.0_prec, 12.0_prec, 19.0_prec, 21.0_prec,  3.0_prec,   &
     11.0_prec, 18.0_prec, 25.0_prec,  2.0_prec,  9.0_prec /),&
     shape(A) ) )

  Alud = A

  b = one

  write(*,'(A3)') 'A ='
  do i = 1, N
    write(*,FMT=rowfmt) (A(i,j), j = 1,N)
  end do
  write(*,*)

  write(*,'(A6)') 'A_LU ='
  do i = 1, N
    write(*,FMT=rowfmt) (Alud(i,j), j = 1,N)
  end do
  write(*,*)

  write(*,'(A3)') 'b ='
  do i = 1, N
    write(*,'(F10.4)') b(i)
  end do
  write(*,*)

  call ludcmp( Alud, N, NP, indx, d )

  write(*,'(A6)') 'A_LU ='
  do i = 1, N
    write(*,FMT=rowfmt) (Alud(i,j), j = 1,N)
  end do
  write(*,*)

  x = b
  call lubksb( Alud, N, NP, indx, x)
  write(*,'(A3)') 'x ='
  do i = 1, N
    write(*,'(F10.4)') b(i)
  end do
  write(*,*)

  xold = x
  call mprove( A, Alud, N, NP, indx, b, x )
  write(*,'(A4)') 'DX ='
  do i = 1, N
    write(*,*) abs(x(i) - xold(i))
  end do
  write(*,*)

  deallocate( indx, A, Alud, b, x, xold )

end program lu_test
