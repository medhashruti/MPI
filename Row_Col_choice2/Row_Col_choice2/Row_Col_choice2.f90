!  Row_Col_choice2.f90 
!
!  FUNCTIONS:
!  Row_Col_choice2 - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Row_Col_choice2
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program Row_Col_choice2
    
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
    logical :: valid_choice
    
    ! Declare valid_choice variable
    valid_choice = .false.
    
    ! Body of RowMajor
    do i = 1, n
        do j = 1, n
            a(i,j) = i+j
            b(i,j) = i-j
            c_act(i,j) = 0
            c(i,j) = 0
        end do
    end do
    
    !using MPI
    
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, num_procs, ierr)
    
    !if (rank == 0) then
    write(*,*) "total number of processes available: ", num_procs
    write(*,*) "Enter number of processes you want to use: "
    write(*,*) " "
    read(*,*) procs_av
    write(*,*) " "

    
    
    if ( num_procs < 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MATMAT - Fatal error!'
        write ( *, '(a)' ) '  Must have at least 2 processes!'
        stop
    end if

    
    if(rank == 0) then
        write(*,*) "Enter 1 for Row Major operation"
        write(*,*) "Enter 2 for Column Major operation"
        write(*,*) " "
        valid_choice = .false.
        do while (.not. valid_choice)
            read(*,*) choice

            if (choice == 1 .or. choice == 2) then
                valid_choice = .true.
            else
                write(*, *) "Choice not valid! Please enter 1 or 2."
            end if
        end do
        !read(*,*) choice
        
        write(*,*) "number of processes being used: ", num_procs
    
        !if (choice /= 1 .and. choice /= 2) then
         !   write(*, *) "Choice not valid!... Exiting"
          !  stop
        !end if
    end if
    
    call MPI_Bcast(choice, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    if(choice == 1) then
        write(*,*) "Starting Row operation"
        call MPI_Bcast(b, n*n, mpi_real, 0, MPI_COMM_WORLD, ierr)
    else
        write(*,*) "Starting Column operation"
        call MPI_Bcast(a, n*n, mpi_real, 0, MPI_COMM_WORLD, ierr)
    end if
    
    if (rank == 0) then
        num_sent = 0
        
        l = min(num_procs - 1, n)
        if(choice==1) then
            write(*,*) "entered row operation for separating rows of A"
            do j = 1, l
                a_row(1:n) = a(j,1:n)
                
                dest = j
                tag = j
            
                !write(*,*) "Dest and tag: ", dest, tag
                
                call MPI_Send(a_row, n, MPI_Real, dest, tag, MPI_COMM_WORLD, ierr)
            
                num_sent = num_sent+1
            end do
        else
            write(*,*) "entered column operation for separating columns of B"
            do j = 1, l
                b_col(1:n) = b(1:n,j)
            
                dest = j
                tag = j
            
                call MPI_Send(b_col, n, MPI_Real, dest, tag, MPI_COMM_WORLD, ierr)
            
                num_sent = num_sent+1
            end do
        end if
        
        num_recv = 0
        if (choice==1) then
            write(*,*) "Receiving rows of C"
            do while (num_recv < n)
                call MPI_Recv(c_row, n, MPI_Real, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
                num_recv = num_recv+1
                source = status(MPI_Source)
                tag = status(MPI_Tag)
            
                c(tag,1:n) = c_row(1:n)
            
                if(num_sent < n) then
                    num_sent = num_sent+1
                    a_row(1:n) = a(num_sent,1:n)
                    dest = source
                    tag = num_sent
                
                    call MPI_Send(a_row, n, MPI_REAL, dest, tag, MPI_COMM_WORLD, ierr)
                
                else
                    val = 1.0E+00
                    dest = source
                    tag = 0
                    call MPI_Send(val, 1, MPI_REAL, dest, tag, MPI_COMM_WORLD, ierr)
                
                end if
            end do
        else
            write(*,*) "Receiving columns of C"
            do while (num_recv<n) 
                !call MPI_Recv(c_col, n, MPI_Real, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
                num_recv = num_recv+1
                source = status(MPI_Source)
                tag = status(MPI_Tag)
            
                c(1:n, tag) = c_col(1:n)
            
                if(num_sent<n) then
                    num_sent = num_sent+1
                    b_col(1:n) = b(1:n,num_sent)
                    dest = source
                    tag = num_sent
                
                    !call MPI_Send(b_col, n, MPI_Real, dest, tag, MPI_COMM_WORLD, ierr)
                
                else
                    val = 1.0E+00
                    dest = source
                    tag = 0
                    !call MPI_Send(val, 1, MPI_Real, dest, tag, MPI_COMM_WORLD, ierr)
                
                end if
            end do
        end if
        
    else if (rank<= n) then
        if (choice == 1) then
            do
                write(*,*) "Receiving rows of A for matrix multiplication"
                call MPI_Recv(a_row,n,MPI_Real,0,MPI_ANY_TAG,MPI_COMM_WORLD, status, ierr)
            
                tag = status(MPI_Tag)
            
                if(tag == 0) then
                    exit
                end if
                
                do i = 1,n
                    do j = 1,n
                        c_row(j) = c_row(j)+a_row(i)*b(i,j)
                    end do
                end do
                write(*,*) "Row operation performed on row of C"
            
                call MPI_Send(c_row, n, MPI_Real, 0, tag, MPI_COMM_WORLD, ierr)
            
            end do
        else
            do
                !call MPI_Recv(b_col,n,MPI_Real,0,MPI_ANY_TAG,MPI_COMM_WORLD, status, ierr)
            
                tag = status(mpi_tag)
            
                if(tag == 0) then
                    exit
                end if
            
                c_col(1:n) = matmul(a(1:n,1:n),b_col(1:n))
            
            
                !call MPI_Send(c_col, n, MPI_Real, 0, tag, MPI_COMM_WORLD, ierr)
            
            end do
        end if
    end if
    
    if(rank == 0) then
        do i = 1,n
            write(*, '(14f7.0)') (c(i,j),j = 1, n)
        end do
    end if
    
    read(*,*)
        
    
    call MPI_Finalize(ierr)
    

    end program Row_Col_choice2

