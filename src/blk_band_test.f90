program blk_band_test

  use set_precision, only : prec
  use set_constants, only : set_derived_constants, zero, one, two
  use block_matrix_operations, only : blk_bandec, blk_banbks

  implicit none

  integer :: M1, M2, MP, MPL, N, NP, q
  integer, allocatable :: indx(:)
  real(prec), allocatable :: Aold(:,:,:,:), A(:,:,:,:), &
                               AL(:,:,:,:), b(:,:), x(:,:), &
                               tmpA(:,:), tmpB(:,:), tmpC(:,:)
  real(prec) :: d
  character(LEN=30) :: rowfmt
  character(LEN=30) :: rowfmt2
  character(LEN=30) :: rowfmt3
  integer :: i, j, k, L

  M1 = 1
  M2 = 1
  MP  = M1 + M2 + 1
  MPL = M1
  N   = 3
  NP  = N
  q   = 5

  write(rowfmt,'(A,I4,A)') '(',N,'(1X,G14.6))'
  write(rowfmt2,'(A,I4,A)') '(',q,'(1X,F8.2))'
  write(rowfmt3,'(A,I4,A)') '(',q,'(1X,G14.6))'

  allocate( indx(N), &
          Aold(q,q,N,M1+M2+1), &
          A(q,q,N,M1+M2+1), &
          AL(q,q,N,M1), &
          b(q,N), &
          x(q,N), &
          tmpA(q,q), tmpB(q,q), tmpC(q,q) )
  Aold = zero
  A = zero
  AL = zero
  b = real(reshape( (/(j, j = 1,N*q)/), shape(b) ),prec )
  x = b

  tmpA = transpose( reshape(                                  &
  (/ 17.0_prec, 24.0_prec,  1.0_prec,  8.0_prec, 15.0_prec,   &
     23.0_prec,  5.0_prec,  7.0_prec, 14.0_prec, 16.0_prec,   &
      4.0_prec,  6.0_prec, 13.0_prec, 20.0_prec, 22.0_prec,   &
     10.0_prec, 12.0_prec, 19.0_prec, 21.0_prec,  3.0_prec,   &
     11.0_prec, 18.0_prec, 25.0_prec,  2.0_prec,  9.0_prec /),&
     shape(tmpA) ) )
  tmpB = -tmpA
  tmpC = two*tmpA

  do i = 2,N
    A(:,:,i,1) = tmpA
  end do
  do i = 1,N
    A(:,:,i,2) = tmpC
  end do
  do i = 1,N-1
    A(:,:,i,3) = tmpB
  end do

  write(*,'(A3)') 'A ='
  do L = 1,3
  do k = 1,N
  write(*,'(A6,I1,A1,I1,A3)') 'A(:,:,',k,',',L,') ='
  do i = 1,q
    write(*,FMT=rowfmt2) (A(i,j,k,L), j = 1,q)
  end do
  write(*,*)
  end do
  end do

  Aold = A

  call blk_bandec( A, N, q, M1, M2, NP, MP, AL, MPL, indx, d )

  write(*,'(A5)') 'A_U ='
  do L = 1,3
  do k = 1,N
  write(*,'(A7,I1,A1,I1,A3)') 'AU(:,:,',k,',',L,') ='
  do i = 1,q
    write(*,FMT=rowfmt3) (A(i,j,k,L), j = 1,q)
  end do
  write(*,*)
  end do
  end do
  call blk_banbks( A, N, q, M1, M2, NP, MP, AL, MPL, indx, x )
  write(*,'(A5)') 'A_L ='
  do L = 1,1
  do k = 1,N
  write(*,'(A7,I1,A1,I1,A3)') 'AU(:,:,',k,',',L,') ='
  do i = 1,q
    write(*,FMT=rowfmt3) (A(i,j,k,L), j = 1,q)
  end do
  write(*,*)
  end do
  end do

  write(*,'(A3)') 'x ='
  do i = 1,q
    write(*,FMT=rowfmt) (x(i,j), j = 1,N)
  end do
  write(*,*)

  deallocate( indx, Aold, A, AL, b, x, tmpA, tmpB, tmpC )

end program blk_band_test
