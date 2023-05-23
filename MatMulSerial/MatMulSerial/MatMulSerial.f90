!  MatMulSerial.f90 
!
!  FUNCTIONS:
!  MatMulSerial - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: MatMulSerial
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program MatMulSerial
    !use :: random_numbers
    
    implicit none

    ! Variables
    integer, parameter :: n = 3
    integer :: i, j, k
    real :: a(n,n), b(n,n), c(n,n)
    
    ! Body of MatMulSerial
    call random_seed()

    do i = 1, n
        do j = 1, n
            call RANDOM_NUMBER(a(i,j))
            call RANDOM_Number(b(i,j))
        end do
    end do
    
     write(*, '(a)') "Matrix A:"
    do i = 1, n
        write(*, '(3f8.4)') (a(i,j),j = 1, n)
    end do
   
    write(*, '(a)') "Matrix B:"
    do i = 1, n
        write(*, '(3f8.4)') (b(i,j),j = 1, n)
    end do
    
    do i = 1, n
        do j = 1, n
            c(i,j) = 0.0
            do k = 1, n
                c(i, j) = c(i, j) + a(i, k)*b(k, j)
            end do
        end do
    end do
    
    write(*, '(a)') "Result Matrix C:"
    do i = 1, n
        write(*, '(3f8.4)') (c(i,j),j = 1, n)
    end do
    Read(*,*)
    
    end program MatMulSerial

