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
    integer, parameter :: n = 4
    integer :: i,j,k
    real :: a(n,n), b(n,n), c(n,n), local_A(n,n), local_C(n,n)
    
    integer :: rank, num_procs, local_n, local_start, local_end, ierr, remaining_elements, offset
    integer, allocatable :: local_array(:)
    integer :: global_sum, local_sum, recv_sum

    ! Body of MatrixMul
    call RANDOM_SEED()
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, num_procs, ierr)
    
    local_n = n / num_procs
    local_start = rank * local_n + 1
    local_end = (rank + 1) * local_n
    
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
    
    !temp = min(num_procs - 1, n)
    !do j = 1,temp
     !   a_col(1:n) = a(1:n,j)
      !  destination = j
       ! tag = j
        
       ! call MPI_Send(a_col, n, MPI_Real, destination, tag, MPI_COMM_WORLD, ierr)
        
       ! num_sent = num_sent + 1
    !end do
    
    call MPI_Scatter(a, local_n * n, MPI_REAL, local_A, local_n * n, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
    
    !entire matrix-B copied to all processes
    call MPI_Bcast(b, n*n, MPI_Real, 0, MPI_COMM_WORLD, ierr)
    
    !local matrix mul
     do i = local_start, local_end
        do j = 1, n
            local_C(i, j) = 0.0
            do k = 1, n
                local_C(i, j) = local_C(i, j) + local_A(i, k) * B(k, j)
            end do
        end do
     end do
     
     ! Gather local results to root process
    call MPI_Gather(local_C, local_n * n, MPI_REAL, c, local_n * n, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
    
    ! Print the final matrix multiplication result from the root process
    if (rank == 0) then
        write(*, '(a)') "Matrix C (Result):"
        do i = 1, n
            do j = 1, n
                write(*, '(f5.1)', advance="no") c(i, j)
            end do
            write(*,*)
        end do
    end if
    
    !write(*,'(F6.2)') c(i, j)
    write(*, '(a)') "Matrix C:"
    do i = 1, n
        write(*, '(3f8.4)') (c(i,j),j = 1, n)
    end do
    
    ! Finalize MPI
    call MPI_Finalize(ierr)
    
    read(*,*)
    end program MatrixMul

