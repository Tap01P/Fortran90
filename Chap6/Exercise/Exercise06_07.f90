module solver
    contains
    function Gauss_Jordan(a, b, n) result(x)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: a(n, n), b(n)
        integer i, j, m, s, k
        real(8) adiagmax, t, ar, adummy(n, n), x(n), w(n), det, adiag(n)
        adummy(:, :) = a(:, :)
        x(:) = b(:)
        s = 1
        det = 1.0d0
        do i = 1, n
            m = i
            adiagmax = abs(a(i, i))
            do j = i+1, n
                if (abs(a( j, i)) > adiagmax) then
                    m = j
                    adiagmax = a(j, i)
                end if
            end do
            if ( adiagmax < 1.0e-15) then
                det = 0.0d0
                write(*, *) 'det(A) =', det
                stop 'A is singular'
                ! stopping when A is a singular matrix
            end if
            if (i /= m) then
            ! swapping k the row and m the row
                w(   i:n) = adummy(i, i:n)
                adummy(i, i:n) = adummy(m, i:n)
                adummy(m, i:n) = w(   i:n)
                t    = x(i)
                x(i) = x(m)
                x(m) = t  
                s = s * (-1)
                write(*, *) 'swapping', i, m
            end if
            ! Gauuss-Jordan Method
            ar = 1.0d0 / adummy(i, i)
            adiag(i) = adummy(i, i)
            adummy(i, i) = 1.0d0
            adummy(i, i+1:n) = ar * adummy(i, i+1:n)
            x(i) = ar * x(i)
            do j = 1, n
                if (j /= i) then
                    adummy(j, i+1:n) = adummy(j, i+1:n) - adummy(j, i) * adummy(i, i+1:n)
                    x(j) = x(j) - adummy(j, i) * x(i)
                    adummy(j, i) = 0.0d0
                end if
            end do
        end do
                
        do k = 1, n
            det = det * adiag(k)
        end do
        det = dble(s) * det
        write(*, *) 'det(A) =', det
                
    end function Gauss_Jordan
end module solver
    
    program main
        use solver
        implicit none
        integer, parameter :: n = 3
        real(8) a(n, n), b(n), x(n)
    
        a(1, 1) = 1.0d0
        a(1, 2) = 1.0d0
        a(1, 3) = 1.0d0
    
        a(2, 1) = 1.0d0
        a(2, 2) = 0.0d0
        a(2, 3) = 1.0d0
    
        a(3, 1) = 2.0d0
        a(3, 2) = 5.0d0
        a(3, 3) = 1.0d0
    
        b(1) = 0.2d0
        b(2) = 1.3d0
        b(3) = 1.1d0
    
        x(:) = Gauss_Jordan(a, b, n)
        write(*, *) x
        write(*, *) matmul(a,x)
    end program main