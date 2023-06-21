!  StupidWay.f90 
!
!  FUNCTIONS:
!  StupidWay - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: StupidWay
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program StupidWay

    implicit none

    ! Variables
    integer :: i,j,k,l
    integer, parameter :: n=4
    real :: a(n,n), b(n,n), c(n,n), c_act(n,n)
    
    ! Body of StupidWay
    do i = 1, n
        do j = 1, n
            a(i,j) = i+j
            b(i,j) = i-j
        end do
    end do
    
    !do i = 1, n
    !    do j = 1, n
    !        write(*, '(f5.1)') a(i, j)
    !    end do
    !    write(*,*)
    !end do
    
    !do i = 1, n
    !    do j = 1, n
    !        write(*, '(f5.1)') b(i, j)
    !    end do
    !    write(*,*)
    !end do
 
    !multiplication
    do i = 1, n
        do j = 1, n
            c(i,j) = 0
            c_act(i,j) = 0
        end do
    end do
    
    !actual solution
    do i = 1,n
        do j = 1,n
            do k = 1,n
                c_act(i,j) = c_act(i,j)+a(i,k)*b(k,j)
            end do
        end do
    end do
    
    !do i = 1, n
    !    do j = 1, n
    !        write(*, '(f15.1)') c_act(i, j)
    !    end do
    !    write(*,*)
    !end do
    
    
    !stupid way
    do i = 1,n
        do j = 1,n
            do k = 1,n
                c(i,k) = c(i,k)+a(i,j)*b(j,k)
            end do
        end do
    end do
    
    !do i = 1, n
    !    do j = 1, n
    !        write(*, '(f15.1)') c(i, j)
    !    end do
    !    write(*,*)
    !end do
    
    do i = 1, n
        write(*, '(4f5.0)') (c(i,j),j = 1, n)
    end do


    read(*,*)
    end program StupidWay

