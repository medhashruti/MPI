!  EqualPatchMatMul.f90 
!
!  FUNCTIONS:
!  EqualPatchMatMul - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: EqualPatchMatMul
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program EqualPatchMatMul
    use mpi
    use ieee_arithmetic
    implicit none

    ! Variables
    integer, parameter :: n = 9
    integer :: i,j,k
    integer :: proc = 1
    integer :: p, q
    real :: a(n,n),b(n,n),c(n,n)
    real, allocatable :: a_sub(:,:), b_sub(:,:), c_sub(:,:)
    
     
    !initialize mpi
    integer :: rank, size, ierr !,p ,q
    call MPI_Init(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    
    !grid dimension
    p = int(sqrt(real(size)))
    q = size/p
    
    
    do i = 1, n
        do j = 1, n
            a(i,j) = i+j
            b(i,j) = i-j
        end do
    end do
    
    !a=reshape([(i+j, i=1,n), j=1,n], shape(a))
    !b=reshape([(i+j, i=1,n), j=1,n], shape(b))
    
    ! Body of EqualPatchMatMul
    allocate(a_sub(n/p, n/p), b_sub(n/p, n/p), c_sub(n/p, n/p))
    
    do i = 1,p
        do j = 1,p
            a_sub = a((i-1)*(n/p)+1:i*(n/p), (j-1)*(n/p)+1:j*(n/p))
            b_sub = b((i-1)*(n/p)+1:i*(n/p), (j-1)*(n/p)+1:j*(n/p))
            !sending each submatrix to corresponding process
            call MPI_Send(a_sub,(n/p)*(n/p),MPI_Real,proc,0,MPI_COMM_WORLD,ierr)
            call MPI_Send(b_sub,(n/p)*(n/p),MPI_Real,proc,0,MPI_COMM_WORLD,ierr)
            proc = proc+1
        end do
    end do
    
    !receiving each submatrix into the correct process
    call MPI_Recv(a_sub,(n/p)*(n/p),MPI_Real,0,0,MPI_COMM_WORLD,ierr)
    call MPI_Recv(b_sub,(n/p)*(n/p),MPI_Real,0,0,MPI_COMM_WORLD,ierr)
    
    !matrix multiplication
    c_sub = matmul(a_sub,b_sub)
    
    !Gather
    if(rank==0) then
        print *, "Matrix C:"
        do i=1, n
            do j=1,n
                print *, c(i,j)
            end do
        end do
    end if
    
    deallocate(a_sub, b_sub, c_sub)
    
    !Finalize MPI
    call MPI_Finalize(ierr)

    read(*,*)
    end program EqualPatchMatMul

