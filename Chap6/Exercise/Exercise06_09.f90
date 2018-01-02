! 関数f(x) = 0.1 * x^3 + 0.2 * x^2 + 0.5 * x + 1.0を定義し、
! この関数の周りに分布するm個の点(x_i, y_i)を乱数を用いて生成する。
! ただし、x_iには -10 <= x <= 10 の区間を m - 1 等分した値を入れる。
! 生成した点はファイルに出力する。
program main
    implicit none
    integer m, i
    real(8), allocatable :: x(:), y(:)
    write(*, '(a)', advance = 'no') 'input data number m:'
    read(*, *) m
    allocate(x(m), y(m))
    do i = 1, m
        x(i) = -10.0 + 20.0d0 * dble(i - 1) / dble(m - 1) 
    end do
    do i = 1, m
        y(i) = sample(x(i))
    end do
    open(10, file = 'plotdata.dat')
    do i = 1, m
        write(10, *) x(i), y(i)
    end do
    write(10,*) ''
    close(10)
contains
    function sample(x) result(y)
        real(8) :: x, y, r0, r, s, t, u, v
        call RANDOM_SEED
        call RANDOM_NUMBER(r0)
        r = r0 * 2.0d0 - 1.0d0
        s = 0.1d0 * x * x * x
        t = 0.2d0 * x * x
        u = 0.5d0 * x
        v = 1.0d0 + r
        y = s + t + u + v 
    end function sample
end program main