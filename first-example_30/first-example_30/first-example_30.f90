PROGRAM first_example

!==============================================================!
!                                                              !
! This file has been written as a sample solution to an        !
! exercise in a course given at the High Performance           !
! Computing Centre Stuttgart (HLRS).                           !
! The examples are based on the examples in the MPI course of  !
! the Edinburgh Parallel Computing Centre (EPCC).              !
! It is made freely available with the understanding that      !
! every copy of this file must include this header and that    !
! HLRS and EPCC take no responsibility for the use of the      !
! enclosed teaching material.                                  !
!                                                              !
! Authors: Joel Malard, Alan Simpson,            (EPCC)        !
!          Rolf Rabenseifner, Traugott Streicher (HLRS)        !
!                                                              !
! Contact: rabenseifner@hlrs.de                                !
!                                                              !
! Purpose: A first MPI example for this course                 !
!                                                              !
! Contents: F-Source                                           !
!                                                              !
!==============================================================!

  USE mpi_f08

  IMPLICIT NONE

  INTEGER :: n               ! application-related data
  DOUBLE PRECISION :: result ! application-related data
  INTEGER :: my_rank, num_procs, rank  ! MPI-related data

  CALL MPI_Init()

  CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank)
  CALL MPI_Comm_size(MPI_COMM_WORLD, num_procs)

  IF (my_rank == 0) THEN  
    !  reading the application data "n" from stdin only by process 0:
    WRITE(*,*) "Enter the number of elements (n):"
    READ(*,*) n
  ENDIF
  !  broadcasting the content of variable "n" in process 0 
  !  into variables "n" in all other processes:
  CALL MPI_Bcast(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD)

  !  doing some application work in each process, e.g.:
  result = 1.0 * my_rank * n
  WRITE(*,'(A,I3,A,I3,A,I2,A,I5,A,F9.2)') &
   &        'I am process ', my_rank, ' out of ', num_procs, &
   &        ' handling the ', my_rank, 'th part of n=', n, ' elements, result=', result

  IF (my_rank /= 0) THEN  
    !  sending some results from all processes (except 0) to process 0:
    CALL MPI_Send(result, 1, MPI_DOUBLE_PRECISION, 0, 99, MPI_COMM_WORLD)
  ELSE
    !  receiving all these messages and, e.g., printing them 
    WRITE(*,'(A,F9.2)') &
     &      'I''m proc 0: My own result is ', result
    DO rank=1, num_procs-1
      CALL MPI_Recv(result, 1, MPI_DOUBLE_PRECISION, rank, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE)
      WRITE(*,'(A,I3,A,F9.2)') &
       &      'I''m proc 0: received result of process ', rank, ' is ', result 
    END DO
  ENDIF

  CALL MPI_Finalize()

  read(*,*)
  
END PROGRAM


