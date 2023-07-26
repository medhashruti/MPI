!  RowMajor.f90 
!
!  FUNCTIONS:
!  RowMajor - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: RowMajor
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program RowMajor

    use mpi
    implicit none

    ! Variables

    integer :: i,j,k,l
    integer, parameter :: n=4
    real :: a(n,n), b(n,n), c(n,n), c_act(n,n), a_row(n), c_row(1:n), b_col(n), c_col(1:n)
    real :: val
    
    integer :: rank, num_procs, ierr, choice
    integer :: num_sent, num_recv, dest, tag, source, status(mpi_status_size)
    integer :: procs_av
    
    logical :: terminate
    
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
    
    !if (rank == 0) then
    write(*,*) "total number of processes available: ", procs_av
    write(*,*) "Enter number of processes you want to use: "
    write(*,*) " "
    read(*,*) num_procs
    write(*,*) " "

    
    
    if ( num_procs < 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MATMAT - Fatal error!'
        write ( *, '(a)' ) '  Must have at least 2 processes!'
        stop
    end if



    !if ( num_procs < 2 ) then
        !call MPI_Abort (  MPI_COMM_WORLD, 1, ierr )
        !stop
    !end if
    
    if(rank == 0) then
        write(*,*) "Enter 1 for Row Major operation"
        write(*,*) "Enter 2 for Column Major operation"
        write(*,*) " "
        read(*,*) choice
        
        write(*,*) "number of processes being used: ", num_procs
    
        if (choice /= 1 .and. choice /= 2) then
            write(*, *) "Choice not valid!... Exiting"
            stop
        end if
    end if
    
    call MPI_Bcast(choice, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    if(choice == 1) then
        !send copy of entire B matrix to all other processes:
        write(*,*) "Starting Row operation"
        !write(*,*) "num procs:", num_procs
        !call MPI_Bcast(b, n*n, mpi_real, 0, MPI_COMM_WORLD, ierr)
    
        !send rows of A to each process and collects
        ! product C(I,1:n) = A(I,1:n)*B(1:n,1:n)
    
        if (rank == 0) then
            call MPI_Bcast(b, n*n, mpi_real, 0, MPI_COMM_WORLD, ierr)
            num_sent = 0
            !write(*,*) "Current rank:", rank
            l = min(num_procs - 1, n)
            !write(*,*) "L = ",l
            do j = 1, l
                a_row(1:n) = a(j,1:n)
                
                dest = j
                tag = j
            
                !write(*,*) "Dest and tag: ", dest, tag
                
                call MPI_Send(a_row, n, MPI_Real, dest, tag, MPI_COMM_WORLD, ierr)
            
                num_sent = num_sent+1
            end do
            
        
            num_recv = 0
            do while (num_recv < n)
                write(*,*) "Entered do loop"
                call MPI_Recv(c_row, n, MPI_Real, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
                write(*,*) "MPI Receive successful"
                num_recv = num_recv+1
                source = status(MPI_Source)
                tag = status(MPI_Tag)
            
                write(*,*) "Tag1", tag
                c(tag,1:n) = c_row(1:n)
            
                if(num_sent < n) then
                    num_sent = num_sent+1
                    a_row(1:n) = a(num_sent,1:n)
                    dest = source
                    tag = num_sent
                
                    write(*,*) "Tag2", tag
                    call MPI_Send(a_row, n, MPI_REAL, dest, tag, MPI_COMM_WORLD, ierr)
                
                else
                    val = 1.0E+00
                    dest = source
                    tag = 0
                    call MPI_Send(val, 1, MPI_REAL, dest, tag, MPI_COMM_WORLD, ierr)
                
                end if
            end do
        
        else if (rank<=n) then
            do
                call MPI_Recv(a_row,n,MPI_Real,0,MPI_ANY_TAG,MPI_COMM_WORLD, status, ierr)
            
                tag = status(MPI_Tag)
                write(*,*) "Tag3", tag
            
                if(tag == 0) then
                    exit
                end if
                
                !c_row(1:n)=0.0
            
                do i = 1,n
                    do j = 1,n
                        c_row(j) = c_row(j)+a_row(i)*b(i,j)
                    end do
                end do
            
            
                call MPI_Send(c_row, n, MPI_Real, 0, tag, MPI_COMM_WORLD, ierr)
            
            end do
        end if
    
        if(rank == 0) then
            do i = 1,n
                write(*, '(4f7.0)') (c(i,j),j = 1, n)
            end do
        end if
        !read(*,*)
        
    else if (choice == 2) then
        !send copy of entire A matrix to all other processes:
        write(*,*) "Starting Column operation"
        !call MPI_Bcast(a, n*n, mpi_real, 0, MPI_COMM_WORLD, ierr)
    
        !send rows of A to each process and collects
        ! product C(I,1:n) = A(I,1:n)*B(1:n,1:n)
    
        if (rank == 0) then
            call MPI_Bcast(a, n*n, mpi_real, 0, MPI_COMM_WORLD, ierr)
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
    end if
    
        
    
    call MPI_Finalize(ierr)
    
    !read(*,*)
    end program RowMajor