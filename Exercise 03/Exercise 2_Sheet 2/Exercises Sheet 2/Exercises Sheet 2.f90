!   ExercisesSheet2.f90 
!   a) Write a function which multiplies an N-element vector with a NxN matrix
!   b) Write a function that solves the system of equations Mx=b
!   M is an upper-diagonal matrix
!   c) Write a code for solving linear systems of equations with the Gauss method with pivoting
    

program matrix_vector_prod
    implicit none
    
    integer, parameter :: kk = SELECTED_REAL_KIND(7)
    real(kind=kk), allocatable :: matrix(:,:), vector(:), outcome(:)
    integer :: n
    
    real(kind=kk), allocatable :: matrix2(:,:), vector2(:), row(:), temp(:), temp2(:), solution(:)
    real(kind=kk) :: t1, sigma, tolerance = 1e-6 
    integer :: i, k, j, index_max
    
! Part 1: Matrix-Vector Multiplication
    print*, "Part 1: Matrix Multiplication"
    n = 3
    allocate(matrix(n, n))
    allocate(vector(n))
    allocate(outcome(n))
    
    matrix = RESHAPE((/1, 4, 7, 2, 5, 8, 3, 6, 9/), (/n, n/))
    vector = (/10, 11, 12/)
    !print*, matrix(1,1), matrix(1,2), matrix(1,3)
    
    print*, "Matrix:"
    print*, matrix
    print*, "Vector:"
    print*, vector
    
    ! perform matrix-vector multiplication
    outcome = prod(matrix, vector, n)
    print*, "Outcome =", outcome
    print*, ""
    
    deallocate(matrix, vector, outcome)
    
    
! ------------------------------------------------------------------------------------
! Part 2: Gauss method
    print*, "Part 2: Gaussian Elimination with pivoting"
    print*, "Test"
    n = 4
    allocate(matrix2(n,n))
    allocate(vector2(n))
    allocate(row(n))
    allocate(temp2(n))
    
    matrix2 = RESHAPE((/1, 0, 1, 2, 0, 1, 2, 1, 1, -2, -1, 3, 2, 0, 0, -2/), (/n, n/))   ! test
    vector2 = (/6, -3, -2, 0/)
    !print*, matrix2(2,1), matrix2(2,2), matrix2(2,3), matrix2(2,4)
    
    print*, "Initial Matrix:"
    print*, matrix2
    print*, "Initial Vector:"
    print*, vector2
    
    ! turn a matrix into an upper-triangular matrix
    main_row: do i = 1, n-1
        row = matrix2(i,:)
        
        ! if the pivot is zero, search for a non-zero pivot
        if (abs(matrix2(i, i)) < tolerance) then
            allocate(temp(n - i))
            
            temp = matrix2(i+1:n, i)   ! NOTE: the index of the vector slice starts from 1
            
            ! find the index of the biggest absolute value in the column
            index_max = maxloc(abs(temp), 1)   
            index_max = i + index_max
            
            if (matrix2(index_max, i) == 0.0_kk) then
                if (i == n-1) then   ! in this case we are done
                    exit
                else
                    cycle   ! otherwise go directly to the next step
                end if
            end if
            
            ! swap the rows of the matrix
            temp2 = matrix2(index_max,:)
            matrix2(index_max,:) = matrix2(i,:)
            matrix2(i,:) = temp2
            
            ! swap the elements of the vector
            t1 = vector2(i)
            vector2(i) = vector2(index_max)
            vector2(index_max) = t1
            
            row = matrix2(i,:)
            
            deallocate(temp)
        end if
        
        ! Gauss method
        rows: do k = i+1, n
            t1 = matrix2(k, i) / matrix2(i, i)
            matrix2(k,:) = t1 * row - matrix2(k,:)
            vector2(k) = t1 * vector2(i) - vector2(k)
        end do rows
    end do main_row
    
    print*, "Upper Triangular Matrix:"
    print*, matrix2
    print*, "Modified Vector:"
    print*, vector2
    
    ! find the solution of the linear system
    allocate(solution(n))
    solution(n) = vector2(n) / matrix2(n,n)
    sigma = 0
    do i = n-1, 1, -1
        do j = i+1, n
            sigma = sigma + matrix2(i,j) * solution(j)
        end do
        solution(i) = (vector2(i) - sigma) / matrix2(i,i)
        sigma = 0
    end do
    
    print*, "Solution of the linear system:"
    print*, solution
    print*, ""
    
    deallocate(matrix2, vector2, row, temp2, solution)
    
! First linear system (see sheet 2)    
! Gauss method without pivoting --------------------------------------------------------------
    print*, "First linear system"
    print*, "Gaussian Elimination without pivoting"
    n = 3
    allocate(matrix2(n,n))
    allocate(vector2(n))
    allocate(row(n))
    allocate(temp2(n))
    
    matrix2 = RESHAPE((/2.0, 0.05, 0.12, 0.1, 4.2, -0.07, -0.2, 0.032, 5.0/), (/n, n/))
    vector2 = (/10, 11, 12/)
    
    print*, "Initial Matrix:"
    print*, matrix2
    print*, "Initial Vector:"
    print*, vector2
    
    ! turn a matrix into an upper-triangular matrix
    main_row2: do i = 1, n-1
        row = matrix2(i,:)
        
        ! Gauss method
        rows2: do k = i+1, n
            t1 = matrix2(k, i) / matrix2(i, i)
            matrix2(k,:) = t1 * row - matrix2(k,:)
            vector2(k) = t1 * vector2(i) - vector2(k)
        end do rows2
    end do main_row2
    
    print*, "Upper Triangular Matrix:"
    print*, matrix2
    print*, "Modified Vector:"
    print*, vector2
    
    ! find the solution of the linear system
    allocate(solution(n))
    solution(n) = vector2(n) / matrix2(n,n)
    sigma = 0
    do i = n-1, 1, -1
        do j = i+1, n
            sigma = sigma + matrix2(i,j) * solution(j)
        end do
        solution(i) = (vector2(i) - sigma) / matrix2(i,i)
        sigma = 0
    end do
    
    print*, "Solution of the linear system:"
    print*, solution
    print*, ""

    deallocate(matrix2, vector2, row, temp2, solution)
    
! Comparison between Gauss method with pivoting and without pivoting --------------------------------------
    print*, "Solution of the linear system solved in the last section with pivoting"
    
    allocate(matrix2(n,n))
    allocate(vector2(n))
    allocate(row(n))
    allocate(temp2(n))
    
    matrix2 = RESHAPE((/2.0, 0.05, 0.12, 0.1, 4.2, -0.07, -0.2, 0.032, 5.0/), (/n, n/))
    vector2 = (/10, 11, 12/)
    
    ! turn a matrix into an upper-triangular matrix
    main_row3: do i = 1, n-1
        row = matrix2(i,:)
        
        ! if the pivot is zero, search for a non-zero pivot
        if (abs(matrix2(i, i)) < tolerance) then
            allocate(temp(n - i))
            
            temp = matrix2(i+1:n, i)   ! NOTE: the index of the vector slice starts from 1
            
            ! find the index of the biggest absolute value in the column
            index_max = maxloc(abs(temp), 1)   
            index_max = i + index_max
            
            if (matrix2(index_max, i) == 0.0_kk) then
                if (i == n-1) then   ! in this case we are done
                    exit
                else
                    cycle   ! otherwise go directly to the next step
                end if
            end if
            
            ! swap the rows of the matrix
            temp2 = matrix2(index_max,:)
            matrix2(index_max,:) = matrix2(i,:)
            matrix2(i,:) = temp2
            
            ! swap the elements of the vector
            t1 = vector2(i)
            vector2(i) = vector2(index_max)
            vector2(index_max) = t1
            
            row = matrix2(i,:)
            
            deallocate(temp)
        end if
        
        ! Gauss method
        rows3: do k = i+1, n
            t1 = matrix2(k, i) / matrix2(i, i)
            matrix2(k,:) = t1 * row - matrix2(k,:)
            vector2(k) = t1 * vector2(i) - vector2(k)
        end do rows3
    end do main_row3
    
    print*, "Upper Triangular Matrix:"
    print*, matrix2
    print*, "Modified Vector:"
    print*, vector2
    
    ! find the solution of the linear system
    allocate(solution(n))
    solution(n) = vector2(n) / matrix2(n,n)
    sigma = 0
    do i = n-1, 1, -1
        do j = i+1, n
            sigma = sigma + matrix2(i,j) * solution(j)
        end do
        solution(i) = (vector2(i) - sigma) / matrix2(i,i)
        sigma = 0
    end do
    
    print*, "Solution of the linear system:"
    print*, solution
    print*, ""

    deallocate(matrix2, vector2, row, temp2, solution)
    
    
! Second linear system (see sheet 2)
! Gauss method without pivoting -----------------------------------------------------------------
    print*, "Second linear system"
    print*, "Gaussian Elimination without pivoting"
    n = 3
    allocate(matrix2(n,n))
    allocate(vector2(n))
    allocate(row(n))
    allocate(temp2(n))
    
    matrix2 = RESHAPE((/1.0, 2.0, 0.0, 1.0, 2.0, 3.0, 0.0, -2.0, 15.0/), (/n, n/))
    vector2 = (/1.0, -2.0, 33.0/)
    
    print*, "Initial Matrix:"
    print*, matrix2
    print*, "Initial Vector:"
    print*, vector2
    
    ! turn a matrix into an upper-triangular matrix
    main_row4: do i = 1, n-1
        row = matrix2(i,:)
        
        ! Gauss method
        rows4: do k = i+1, n
            if (abs(matrix2(i,i)) < tolerance) then   ! if I find a vanishing pivot, skip this section
                print*, "Vanishing pivot"
                allocate(solution(n))
                go to 10
            end if
            
            t1 = matrix2(k, i) / matrix2(i, i)
            matrix2(k,:) = t1 * row - matrix2(k,:)
            vector2(k) = t1 * vector2(i) - vector2(k)
        end do rows4
    end do main_row4
    
    print*, "Upper Triangular Matrix:"
    print*, matrix2
    print*, "Modified Vector:"
    print*, vector2
    
    ! find the solution of the linear system
    allocate(solution(n))
    solution(n) = vector2(n) / matrix2(n,n)
    sigma = 0
    do i = n-1, 1, -1
        do j = i+1, n
            sigma = sigma + matrix2(i,j) * solution(j)
        end do
        solution(i) = (vector2(i) - sigma) / matrix2(i,i)
        sigma = 0
    end do
    
    print*, "Solution of the linear system:"
    print*, solution
10  print*, ""

    deallocate(matrix2, vector2, row, temp2, solution)
    
! Comparison between Gauss method with pivoting and without pivoting --------------------------------------
    print*, "Solution of the linear system solved in the last section with pivoting"
    
    allocate(matrix2(n,n))
    allocate(vector2(n))
    allocate(row(n))
    allocate(temp2(n))
    
    matrix2 = RESHAPE((/1.0, 2.0, 0.0, 1.0, 2.0, 3.0, 0.0, -2.0, 15.0/), (/n, n/))
    vector2 = (/1.0, -2.0, 33.0/)
    
    ! turn a matrix into an upper-triangular matrix
    main_row5: do i = 1, n-1
        row = matrix2(i,:)
        
        ! if the pivot is zero, search for a non-zero pivot
        if (abs(matrix2(i, i)) < tolerance) then
            allocate(temp(n - i))
            
            temp = matrix2(i+1:n, i)   ! NOTE: the index of the vector slice starts from 1
            
            ! find the index of the biggest absolute value in the column
            index_max = maxloc(abs(temp), 1)   
            index_max = i + index_max
            
            if (matrix2(index_max, i) == 0.0_kk) then
                if (i == n-1) then   ! in this case we are done
                    exit
                else
                    cycle   ! otherwise go directly to the next step
                end if
            end if
            
            ! swap the rows of the matrix
            temp2 = matrix2(index_max,:)
            matrix2(index_max,:) = matrix2(i,:)
            matrix2(i,:) = temp2
            
            ! swap the elements of the vector
            t1 = vector2(i)
            vector2(i) = vector2(index_max)
            vector2(index_max) = t1
            
            row = matrix2(i,:)
            
            deallocate(temp)
        end if
        
        ! Gauss method
        rows5: do k = i+1, n
            t1 = matrix2(k, i) / matrix2(i, i)
            matrix2(k,:) = t1 * row - matrix2(k,:)
            vector2(k) = t1 * vector2(i) - vector2(k)
        end do rows5
    end do main_row5
    
    print*, "Upper Triangular Matrix:"
    print*, matrix2
    print*, "Modified Vector:"
    print*, vector2
    
    ! find the solution of the linear system
    allocate(solution(n))
    solution(n) = vector2(n) / matrix2(n,n)
    sigma = 0
    do i = n-1, 1, -1
        do j = i+1, n
            sigma = sigma + matrix2(i,j) * solution(j)
        end do
        solution(i) = (vector2(i) - sigma) / matrix2(i,i)
        sigma = 0
    end do
    
    print*, "Solution of the linear system:"
    print*, solution

    deallocate(matrix2, vector2, row, temp2, solution)
    
    
contains   ! **********************************************************************************

    ! function to calculate the product
    function prod(m, v, n) result(res)
        implicit none
        
        integer, parameter :: kk = SELECTED_REAL_KIND(7)
        integer, intent(in) :: n
        integer :: i, j
        real(kind=kk), intent(in) :: m(n, n), v(n)
        real(kind=kk), dimension(n) :: res
        
        res = 0.0_kk   ! initialization of the vector
        
        ! perform multiplication
        row: do i = 1, n
            column: do j = 1, n
                res(i) = res(i) + m(i, j) * v(j)
            end do column
        end do row
    end function prod

end program matrix_vector_prod