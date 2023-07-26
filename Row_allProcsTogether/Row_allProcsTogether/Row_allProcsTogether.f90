!  Row_allProcsTogether.f90 
!
!  FUNCTIONS:
!  Row_allProcsTogether - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Row_allProcsTogether
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program Row_allProcsTogether

    use mpi
    implicit none

    ! Variables

    integer :: i,j,k,l
    !, extra_rows, num_rows_per_process
    integer, parameter :: n=28
    real :: a(n,n), b(n,n), c(n,n), c_act(n,n), a_row(n), c_row(1:n), b_col(n), c_col(1:n)
    real :: val
    
    integer :: rank, num_procs, new_rank, new_num_procs, ierr, choice
    integer :: num_sent, num_recv, dest, tag, source, status(mpi_status_size)
    integer :: procs_av, subset_color, subset_key, new_comm
    
    !logical :: terminate
    
    ! Body of RowMajor
    do i = 1, n
        do j = 1, n
            a(i,j) = i+j
            b(i,j) = i-j
            c_act(i,j) = 0
            c(i,j) = 0
        end do
    end do

    !exact solution
    
    c_act(1:n,1:n) = matmul(a(1:n,1:n),b(1:n,1:n))
    
    !do i = 1, n
     !   write(*, '(4f5.0)') (c_act(i,j),j = 1, n)
    !end do
    
   
    !using MPI
    
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, procs_av, ierr)
    
    
    !send copy of entire B matrix to all other processes:
    write(*,*) "Processes available: ", procs_av
    write(*,*) "Enter number of processes you want to use: "
    read(*,*) num_procs
    !write(*,*) "num procs:", num_procs
    
    if(rank<=num_procs) then
        subset_color = 0
    else
        subset_color = MPI_UNDEFINED
    endif
    
    subset_key = rank
    
    ! Create a new communicator with the desired subset of processes
    call MPI_Comm_split(MPI_COMM_WORLD, subset_color, subset_key, new_comm, ierr)
    call MPI_Comm_rank(new_comm, new_rank, ierr)
    call MPI_Comm_size(new_comm, new_num_procs, ierr)
    
    if (subset_color == 0) then
        call MPI_Bcast(b, n*n, mpi_real, 0, new_comm, ierr)
    
        !send rows of A to each process and collects
        ! product C(I,1:n) = A(I,1:n)*B(1:n,1:n)
    !extra_rows = MOD(n, num_procs)
    !write(*,*) "Reminder: ", extra_rows
    
    !num_rows_per_process = n / (num_procs - 1)
    !write(*,*) "no. of rows per process: ",num_rows_per_process
    
    
        num_sent = 0
        if (new_rank <n ) then
    !call MPI_Bcast(b, n*n, mpi_real, 0, MPI_COMM_WORLD, ierr)
            !l = min(num_procs - 1, n)
            do j = new_rank+1, n, new_num_procs
            
                a_row(1:n) = a(j,1:n)
                
                dest = j
                tag = j
            
            !write(*,*) "Dest and tag: ", dest, tag
                
                call MPI_Send(a_row, n, MPI_Real, dest, tag, new_comm, ierr)
            
                num_sent = num_sent+1
            end do
            
        
            num_recv = 0
            do while (num_recv < n)
            !write(*,*) "Entered do loop"
                call MPI_Recv(c_row, n, MPI_Real, MPI_ANY_SOURCE, MPI_ANY_TAG, new_comm, status, ierr)
            !write(*,*) "MPI Receive successful"
                num_recv = num_recv+1
                source = status(MPI_Source)
                tag = status(MPI_Tag)

                c(tag,1:n) = c_row(1:n)
            
                if(num_sent < n) then
                    num_sent = num_sent+1
                    a_row(1:n) = a(num_sent,1:n)
                    dest = source
                    tag = num_sent

                    call MPI_Send(a_row, n, MPI_REAL, dest, tag, new_comm, ierr)
                
                else
                    val = 1.0E+00
                    dest = source
                    tag = 0
                    call MPI_Send(val, 1, MPI_REAL, dest, tag, new_comm, ierr)
                
                end if
            end do
        
        else if (new_rank<=n) then
            do
                call MPI_Recv(a_row,n,MPI_Real,0,MPI_ANY_TAG,new_comm, status, ierr)
            
                tag = status(MPI_Tag)
            
                if(tag == 0) then
                    exit
                end if
                
                !c_row(1:n)=0.0
            
                do i = 1,n
                    do j = 1,n
                        c_row(j) = c_row(j)+a_row(i)*b(i,j)
                    end do
                end do
            
            
                call MPI_Send(c_row, n, MPI_Real, 0, tag, new_comm, ierr)
            
            end do
        end if
        
        ! Finalize MPI for the subset of processes
        call MPI_Finalize(ierr)
    else
        ! Finalize MPI for the excluded processes
        call MPI_Finalize(ierr)
    endif
    
    if(rank == 0) then
        do i = 1,n
            write(*, '(28f7.0)') (c(i,j),j = 1, n)
        end do
    end if
    
    read(*,*)
        
    
    call MPI_Finalize(ierr)


    end program Row_allProcsTogether

