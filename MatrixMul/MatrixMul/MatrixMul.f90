!  MatrixMul.f90 
!
!  FUNCTIONS:
!  MatrixMul - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: MatrixMul
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program MatrixMul

    use mpi
    implicit none

    ! Variables
    integer, parameter :: n = 3
    integer :: i,j,k
    real :: a(n,n), b(n,n), c(n,n)
    
    integer :: rank, num_procs, local_n, local_start, local_end, i, ierr, remaining_elements, offset
    integer, allocatable :: local_array(:)
    integer :: global_sum, local_sum, recv_sum

    ! Body of MatrixMul
    call RANDOM_SEED()
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, num_procs, ierr)
    
    !local_n = n / num_procs
    !local_start = rank * local_n + 1
    !local_end = (rank + 1) * local_n
    
    do i = 1,n
        do j = 1,n
            call RANDOM_NUMBER(a(i,j))
            call RANDOM_NUMBER(b(i,j))
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

    !divide tasks
    !scatter matrix A
    temp = min(num_procs - 1, n)
    do j = 1,temp
        a_col(1:n) = a(1:n,j)
        destination = j
        tag = j
        
        call MPI_Send(a_col, n, MPI_Real, destination, tag, MPI_COMM_WORLD, ierr)
        
        num_sent = num_sent + 1
    end do
    
    !entire matrix-B copied to all processes
    call MPI_Bcast(b, n*n, MPI_Real, 0, MPI_COMM_WORLD, ierr)
    
    !local matrix mul
    
    
    
    end program MatrixMul

