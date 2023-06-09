!  MatMul.f90 
!
!  FUNCTIONS:
!  MatMul - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: MatMul
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program MatMul

    use mpi
    implicit none

    ! Variables
    integer, parameter :: n = 5
    integer :: i,j,k
    real :: a(n,n), b(n,n), c(n,n), local_A(n,n), local_C(n,n)
    
    integer :: rank, num_procs, local_n, local_start, local_end, ierr, remaining_elements, offset
    integer, allocatable :: local_array(:)
    integer :: global_sum, local_sum, recv_sum
    real :: start_time, end_time, elapsed_time
    !call SYSTEM_CLOCK(start_time)
    

    ! Body of MatrixMul
    call RANDOM_SEED()
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, num_procs, ierr)
    
    start_time = MPI_Wtime()
    
    local_n = n / num_procs
    local_start = rank * local_n + 1
    local_end = (rank + 1) * local_n
    
    do i = 1,n
        do j = 1,n
            !call RANDOM_NUMBER(a(i,j))
            a(i,j) = i+j
            !call RANDOM_NUMBER(b(i,j))
            b(i,j) = abs(i-j)
        end do
    end do
    
    do i = 1,n
        do j = 1,n
            c(i,j) = 0.0
            do k = 1,n
                c(i,j) = c(i,j) + a(i,k) * b(k,j)
            end do
        end do
    end do
    
    write(*, '(a)') "Matrix C (Exact Result):"
    do i = 1, n
        do j = 1, n
            write(*, '(f5.2)', advance="no") c(i, j)
        end do
        write(*,*)
    end do
    
    !exact value of C:
    !if ( rank == 0 ) then

    !c(1:n,1:n) = matmul( a(1:n,1:n), b(1:n,1:n) )

    !write ( *, '(a)' ) ' '
    !write ( *, '(a)' ) 'MATMAT - Master process:'
    !write ( *, '(a)' ) '  Initial 5 x 5 block of exact product matrix C:'
    !write ( *, '(a)' ) ' '

    !do i = 1, n
     ! write ( *, '(50g14.6)' ) c(i,1:n)
    !end do
    !end if
    
    
    !write(*, '(a)') "Matrix A:"
    !do i = 1, n
    !    write(*, '(3f8.4)') (a(i,j),j = 1, n)
    !end do
   
    !write(*, '(a)') "Matrix B:"
    !do i = 1, n
    !    write(*, '(3f8.4)') (b(i,j),j = 1, n)
    !end do

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
                write(*, '(f5.2)', advance="no") c(i, j)
            end do
            write(*,*)
        end do
    end if
    
    !write(*,'(F6.2)') c(i, j)
    !write(*, '(a)') "Matrix C:"
    !do i = 1, n
    !    write(*, '(3f8.4)') (c(i,j),j = 1, n)
    !end do
    
    ! Finalize MPI
    !call MPI_Finalize(ierr)
    
    !call SYSTEM_CLOCK(end_time)
    end_time = MPI_Wtime()
    elapsed_time = end_time - start_time
    
    write(*, '(a, f6.2,"seconds")') "Run time: ", elapsed_time
    !, "for process: ",rank 
    call MPI_Finalize(ierr)
    read(*,*)
 



    end program MatMul

