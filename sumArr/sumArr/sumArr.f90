!  sumArr.f90 
!
!  FUNCTIONS:
!  sumArr - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: sumArr
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program sumArr
    use mpi
    implicit none

    !read(*,*)
    ! Variables
    integer, parameter :: n = 24
    integer :: rank, num_procs, local_n, i, ierr, remaining_elements, offset
    integer, allocatable :: local_array(:)
    integer :: global_sum, local_sum, recv_sum
    
    ! Body of sumArr
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, num_procs, ierr)
    !write(*,*) "rank= ", rank, "number of processors: ", num_procs
    
    local_n = n/num_procs
!    write(*,*), "Local n : ", local_n
    !allocate(local_array(local_n))
    if (n < num_procs) then
        remaining_elements = n - (local_n * num_procs)
    else
        remaining_elements = mod(n, num_procs)
    end if
    
    if (rank < remaining_elements) then
        local_n = local_n + 1
        offset = rank * local_n
    else
        offset = remaining_elements * (local_n + 1) + (rank - remaining_elements) * local_n
    end if
    
    allocate(local_array(local_n))
    
    do i = 1, local_n
        local_array(i) = offset + i
    end do
    
    !do i = 1, local_n
     !   local_array(i) = rank*local_n + i
        !write(*, local_array(i))
    !end do
    
    local_sum = sum(local_array)
    write(*,*) "local sum: ", local_sum, "for rank: ", rank
    !call MPI_Allreduce(local_sum, global_sum, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD)
    if (rank == 0) then
        global_sum = local_sum
        do i = 1, num_procs-1
            call MPI_Recv(recv_sum, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
            global_sum = global_sum + recv_sum
        end do
    else
        call MPI_Send(local_sum, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr)
    end if
    
    if (rank == 0) then
        write(*, '(a,i)') "global sum: ", global_sum
    end if
    
    read(*,*)
    
    deallocate(local_array)

    end program sumArr

