program main_program

  use set_precision, only : dp, prec
  use set_constants, only : set_derived_constants, zero, one

  implicit none

  external DGETRF, DGETRS
  character(1) :: TRANS
  integer :: N, NRHS, LDA, LDB, INFO
  integer, allocatable :: IPIV(:)
  real(prec), allocatable :: A(:,:), B(:,:)

  character(LEN=30) :: rowfmt
  integer :: i, j, k

  TRANS = 'N'
  N    = 5
  NRHS = 5
  LDA  = 5
  LDB  = 5
  INFO = 0

  write(rowfmt,'(A,I4,A)') '(',N,'(1X,F10.4))'

  allocate( IPIV(N), A(LDA,N), B(LDB,NRHS) )

  IPIV = 0
  A = transpose( reshape(                                     &
  (/ 17.0_prec, 24.0_prec,  1.0_prec,  8.0_prec, 15.0_prec,   &
     23.0_prec,  5.0_prec,  7.0_prec, 14.0_prec, 16.0_prec,   &
      4.0_prec,  6.0_prec, 13.0_prec, 20.0_prec, 22.0_prec,   &
     10.0_prec, 12.0_prec, 19.0_prec, 21.0_prec,  3.0_prec,   &
     11.0_prec, 18.0_prec, 25.0_prec,  2.0_prec,  9.0_prec /),&
     shape(A) ) )

  B = 0.0_prec
  do i = 1, LDB
    B(i,i) = 1.0_prec
  end do

  write(*,*)
  do i = 1, LDB
    write(*,FMT=rowfmt) (B(i,j), j = 1,NRHS)
  end do
  write(*,*)

  write(*,*)
  do i = 1, LDA
    write(*,FMT=rowfmt) (A(i,j), j = 1,N)
  end do
  write(*,*)

  call DGETRF( N, N, A, LDA, IPIV, INFO )
  call DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )

  write(*,*)
  do i = 1, LDB
    write(*,FMT=rowfmt) (B(i,j), j = 1,NRHS)
  end do
  write(*,*)

  deallocate( IPIV, A, B )

end program main_program
