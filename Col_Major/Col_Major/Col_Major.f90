!  Col_Major.f90 
!
!  FUNCTIONS:
!  Col_Major - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Col_Major
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program Col_Major

    use mpi
    implicit none

    ! Variables

    integer :: i,j,k,l
    integer, parameter :: n=140
    real :: a(n,n), b(n,n), c(n,n), c_act(n,n), b_col(n), c_col(1:n)
    real :: val
    
    integer :: rank, num_procs, ierr
    integer :: num_sent, num_recv, dest, tag, source, status(mpi_status_size)
    
    ! Body of RowMajor
    do i = 1, n
        do j = 1, n
            a(i,j) = 1
            !i+j
            b(i,j) = 1
            !i-j
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
    call MPI_Comm_size(MPI_COMM_WORLD, num_procs, ierr)
    
    !write(*,*) "number of processes: ", num_procs
    
    !send copy of entire A matrix to all other processes:
    
    call MPI_Bcast(a, n*n, mpi_real, 0, MPI_COMM_WORLD, ierr)
    
    !send rows of A to each process and collects
    ! product C(I,1:n) = A(I,1:n)*B(1:n,1:n)
    
    if (rank == 0) then
        num_sent = 0
        
        l = min(num_procs - 1, n)
        do j = 1, l
            b_col(1:n) = b(1:n,j)
            
            dest = j
            tag = j
            
            call MPI_Send(b_col, n, MPI_Real, dest, tag, MPI_COMM_WORLD, ierr)
            
            num_sent = num_sent+1
        end do
        
        num_recv = 0
        do while (num_recv<n)
            
            call MPI_Recv(c_col, n, MPI_Real, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
            num_recv = num_recv+1
            source = status(MPI_Source)
            tag = status(MPI_Tag)
            
            c(1:n, tag) = c_col(1:n)
            
            if(num_sent<n) then
                num_sent = num_sent+1
                b_col(1:n) = b(1:n,num_sent)
                dest = source
                tag = num_sent
                
                call MPI_Send(b_col, n, MPI_Real, dest, tag, MPI_COMM_WORLD, ierr)
                
            else
                val = 1.0E+00
                dest = source
                tag = 0
                call MPI_Send(val, 1, MPI_Real, dest, tag, MPI_COMM_WORLD, ierr)
                
            end if
        end do
        
    else if (rank<=n) then
        do
            call MPI_Recv(b_col,n,MPI_Real,0,MPI_ANY_TAG,MPI_COMM_WORLD, status, ierr)
            
            tag = status(mpi_tag)
            
            if(tag == 0) then
                exit
            end if
            
            c_col(1:n) = matmul(a(1:n,1:n),b_col(1:n))
            
            
            call MPI_Send(c_col, n, MPI_Real, 0, tag, MPI_COMM_WORLD, ierr)
            
        end do
    end if
    
    if(rank == 0) then
        do i = 1,n
            write(*, '(14f7.0)') (c(i,j),j = 1, n)
        end do
    end if
    

    
    call MPI_Finalize(ierr)
    
    read(*,*)
    end program Col_Major

