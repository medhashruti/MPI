!  MatMul_eachRC.f90 
!
!  FUNCTIONS:
!  MatMul_eachRC - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: MatMul_eachRC
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program MatMul_eachRC
    use mpi
    implicit none
    
    integer, parameter :: n = 4
    integer :: i,j,k,r,cb,cin,destin,tag,id
    real :: a(n,n), b(n,n), c(n,n), local_A(n,n), local_B(n,n), local_C(n,n)
    real :: b_col(n)
    
    integer :: rank, num_procs, local_n, local_start, local_end, ierr, remaining_elements, offset
    integer, allocatable :: local_array(:)
    integer :: global_sum, local_sum, recv_sum
    integer :: num_sent
    !integer :: start_time, end_time, elapsed_time
    !call SYSTEM_CLOCK(start_time)

    ! Body of MatrixMul
    call RANDOM_SEED()
    
    do i = 1,n
        do j = 1,n
            call RANDOM_NUMBER(a(i,j))
            call RANDOM_NUMBER(b(i,j))
        end do
    end do
    !write(*, '(a)') "Matrix A:"
    !do i = 1, n
    !    write(*, '(3f8.4)') (a(i,j),j = 1, n)
    !end do


    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, num_procs, ierr)

    !call MPI_Scatter(a, local_n * n, MPI_REAL, local_A, local_n * n, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
    !Exact solution of C = A.B   
    if (rank == 0) then
        c(1:n,1:n) = matmul(a(1:n,1:n), b(1:n,1:n))
        write(*,'(a)') ' Exact result of C: '
        write(*,'(a)') '  '
        
        do i = 1,n
            write(*,'(4g6.2)') c(i,1:n)
        end do
    end if
    
    !distribute each column of A to first 4 processes
    if(id == 0) then
        num_sent = 0
        
        !send column c of B to process c
        cin = min(num_procs-1, n)
        do cb = 1, cin
            b_col(1:n) = b(1:n,cb)
            
            destin = cb
            tag = cb
            
            call MPI_send(b_col,n,MPI_REAL,destin,tag,MPI_COMM_WORLD,ierr)
            num_sent = num_sent + 1
        end do
        
        !similarly send rows of A to left over processes
        !Perform MAtrix multiplication
    !    .
     !   .
      !  .
      !  .
    end if
    
    call MPI_Finalize(ierr)
    read(*,*)
    end program MatMul_eachRC

