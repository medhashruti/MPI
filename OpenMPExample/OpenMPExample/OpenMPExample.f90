program OpenMPMatrixMultiplication
    use omp_lib
    
    implicit none
    integer, parameter :: N = 1000
    real, dimension(N,N) :: A, B, C
    integer :: i, j, k
    
    ! Initialize matrices A and B
    do i = 1, N
        do j = 1, N
            A(i,j) = real(i + j)
            B(i,j) = real(i - j)
        end do
    end do
    
    ! Perform matrix multiplication in parallel using OpenMP
    !$OMP PARALLEL DO PRIVATE(i, j, k)
    do i = 1, N
        do j = 1, N
            C(i,j) = 0.0
            do k = 1, N
                C(i,j) = C(i,j) + A(i,k) * B(k,j)
            end do
        end do
    end do
    !$OMP END PARALLEL DO
    
    ! Print a sample element of matrix C to verify the result
    print *, "Result:", C(1, 1)
    
    read(*,*)
    
end program OpenMPMatrixMultiplication
