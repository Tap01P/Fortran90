! 演習6.7で生成したプロットデータに最小2乗法を用いて、
! 3次の回帰多項式の係数を求める
module solver
implicit none
contains
    subroutine readarray(x, y, m)
        integer, intent(out) :: m
        real(8), allocatable, intent(out) :: x(:), y(:)
        integer i
        write(*, '(a)', advance = 'no') 'data number :'
        read(*, *) m
        allocate(x(m), y(m))
        open(10, file = 'plotdata.dat')
        do i = 1, m
            read(10, *) x(i), y(i)
        end do
        close(10)
        ! do i = 1, m
        !     write(*, '(3e12.4)') x(i), y(i)
        ! end do
    end subroutine readarray

    subroutine creatematrix(x, y, m, n, c, b)
    integer, intent(in) :: m, n
    real(8), intent(in) :: x(m), y(m)
    real(8), intent(out) :: c(n+1, n+1), b(n+1)
    integer i, j, k
    c(:, :) = 0.0d0
    b(:) = 0.0d0
    do i = 1, n+1
        do j = 1, n+1
            do k = 1, m
                c(i, j) = c(i, j) + x(k) ** ((i - 1) + (j - 1))
            end do
        end do
    end do
    do i = 1, n+1
        do k = 1, m
            b(i) = b(i) + y(k) * ( x(k) ** (i - 1) )
        end do
    end do
    end subroutine creatematrix

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

    subroutine plotPolynomial(a, x, m, n)
        integer, intent(in) :: m, n
        real(8), intent(in) :: a(n+1), x(m)
        integer i
        real(8) yplot(m)
        do i = 1, m
            yplot(i) = a(4) * x(i) * x(i) * x(i) &
                      + a(3) * x(i) * x(i) &
                      + a(2) * x(i) &
                      + a(1)
        end do
        open(11, file = 'regression.dat')
        do i = 1, m
            write(11, *) x(i), yplot(i)
        end do
        write(11, *) ''
        close(11)
    end subroutine plotPolynomial

end module solver

program main
    use solver
    implicit none
    integer m, n, i
    real(8), allocatable :: x(:), y(:), c(:, :), b(:), a(:)
    call readarray(x, y, m)
    n = 3
    allocate(c(n+1, n+1), b(n+1), a(n+1))
    call creatematrix(x, y, m, n, c, b)
    a = Gaussian_Elimination(c, b, n+1)
    do i = 1, n+1
        write(*, *) a(i)
    end do
    call plotPolynomial(a, x, m, n)
    ! ------- Gnuplot -------
    open(12, file = 'plotEx10.plt')
        write(12, *) "reset"
        write(12, *) "plot 'plotdata.dat'"
        write(12, *) "replot 'regression.dat' smooth csplines"
        write(12, *) "pause -1"
    close(12)
    call system("gnuplot plotEx10.plt")
end program main