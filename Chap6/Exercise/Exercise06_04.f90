! ピボット選択付きガウスの消去法
!
program main
    use solver
    implicit none
    integer, parameter :: n = 3
    real(8) a(n, n), b(n), x(n), r(n)

    a(1, 1) = 3.5d0
    a(1, 2) = 1.0d0
    a(1, 3) = 0.2d0

    a(2, 1) = 1.3d0
    a(2, 2) = 0.0d0
    a(2, 3) = 1.0d0

    a(3, 1) = 1.9d0
    a(3, 2) = 5.3d0
    a(3, 3) = 0.7d0

    b(1) = 0.2d0
    b(2) = 1.3d0
    b(3) = 1.1d0
    
    x(:) = Gaussian_Elimination(a, b, n)
    r(:) = b(:) - matmul(a, x)
    write(*, *) 'x =', x
    write(*, *) 'r**2 =', DOT_PRODUCT(r, r)
    
end program main

module solver
    implicit none
    contains
    function Gaussian_Elimination(a0, b0, n) result(x)
        integer, intent(in) :: n
        real(8), intent(in) :: a0(n, n), b0(n)
        integer i, j, m
        real(8) adiagmax, t, ar, a(n, n), b(n), x(n), w(n)
        a(:, :) = a0(:, :)
        b(:) = b0(:)
        do i = 1, n
            m = i
            adiagmax = abs(a(i, i))
            do j = i+1, n
                if (abs(a( j, i)) > adiagmax) then
                    m = j
                    adiagmax = a(j, i)
                end if
            end do
            if ( adiagmax < 1.0e-15) stop 'A is singular'
            ! stopping when A is a singular matrix
            if (i /= m) then
            ! swapping k the row and m the row
                w(   i:n) = a(i, i:n)
                a(i, i:n) = a(m, i:n)
                a(m, i:n) = w(   i:n)
                t    = b(i)
                b(i) = b(m)
                b(m) = t  
            end if
            ! Gauussian Elimination Method
            ar = 1.0d0 / a(i, i)
            a(i, i) = 1.0d0
            a(i, i+1:n) = ar * a(i, i+1:n)
            b(i) = ar * b(i)
            do j = i+1, n
                a(j, i+1:n) = a(j, i+1:n) - a(j, i) * a(i, i+1:n)
                b(j) = b(j) - a(j, i) * b(i)
                a(j,     i) = 0.0d0
            end do
        end do
        x(n) = b(n)
        do i = n-1, 1, -1
            x(i) = b(i)
            do j = i+1, n
                x(i) = x(i) - a(i, j) * x(j)  
            end do 
        end do
    end function Gaussian_Elimination
end module solver