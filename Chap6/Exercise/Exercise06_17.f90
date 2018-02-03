! 加速パラメーター omegaを動かし、omegaと反復回数のグラフを出力する
module SetMethods
    implicit none
    contains
    subroutine SetLatticeandPram(n1, n2, x, c, d)
        integer i
        real(8) dx1, dx2
        integer, intent(in) :: n1, n2
        real(8), intent(out) :: x(2, n1, n2), c, d
        dx1 = 1.0d0 / dble(n1 - 1)
        dx2 = 1.0d0 / dble(n2 - 1)
        do i = 1, n1
            x(1, i, 1:n2) = dble(i - 1) * dx1
        end do
        do i = 1, n2
            x(2, 1:n1, i) = dble(i - 1) * dx2
        end do
        c = - (dx2 * dx2) / (2.0d0 * (dx1 * dx1 + dx2 * dx2))
        d = - (dx1 * dx1) / (2.0d0 * (dx1 * dx1 + dx2 * dx2))
    end subroutine SetLatticeandPram
    subroutine SetDBC(n1, n2, x, phi)
        integer i
        real(8) :: pi = 2.0d0 * acos(0.0d0)
        integer, intent(in) :: n1, n2
        real(8), intent(in) :: x(2, n1, n2)
        real(8), intent(out) :: phi(n1, n2)
        phi(1, 1:n2) = 0.0d0
        phi(n1, 1:n2) = 0.0d0
        phi(1:n1, n2) = 0.0d0
        do i = 1, n1
            phi(i, 1) = sin(pi * x(1, i, 1))
        end do
    end subroutine SetDBC
end module SetMethods

program main
    use SetMethods
    implicit none
    integer i, j, itr, k
    integer, parameter :: n1 = 501, n2 = 501, itrmax = 1000
    real(8) x(2, n1, n2), phi(n1, n2), dx1, dx2, c, d, er, rhs, omega0, omega
    real(8) :: pi = 2.0 * acos(0.0d0)
    real(8) :: er0 = 10E-6
    omega0 = 2.0d0 / (1.0d0 + sin(pi / dble(n1 - 1)))
    call SetLatticeandPram(n1, n2, x, c, d)
    open(10, file = 'sol06_17.dat')
    do k = 1, 100
        phi = 0.0
        call SetDBC(n1, n2, x, phi)
        omega = 2.0d0 - dble(k - 1) * 0.01d0 
        ! SOR法の反復計算
        do itr = 1, itrmax
            er = 0.0d0
            do j = 2, n2 - 1
                do i = 2, n1 - 1
                    rhs = - c * (phi(i-1,   j) + phi(i+1,  j)) &
                          - d * (phi(i  , j-1) + phi(i  ,j+1))
                    er = er + (rhs - phi(i, j)) * (rhs - phi(i,j))
                    phi(i, j) = phi(i, j) + omega * (rhs - phi(i, j))
                end do
            end do
            if (er < er0) then
            write(10, *) omega, itr 
            exit
            end if
            if (itr == itrmax) write(10, *) omega, itr
        end do
    end do
    write(10, *) ' '
    close(10)
    write(*, *) omega0
    open(11, file = 'plot06_17.plt')
        write(11, *) 'reset'
        write(11, *) 'plot "sol06_17.dat"'
        write(11, *) 'pause -1'
    close(11)
    call system('gnuplot plot06_17.plt')
end program main