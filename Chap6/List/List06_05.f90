! SOR法を用いた5点差分式の計算プログラム
! 各所のサブルーチン化は演習6.16で行う
! このリストではひとまず全ての処理をメインプログラム内に記述する
program main
    implicit none
    integer i, j, k, itr
    integer, parameter :: n1 = 31, n2 = 31, itrmax = 100
    real(8) x(2, n1, n2), phi(n1, n2)
    real(8) dx1, dx2, rhs, er, c, d, omega
    real(8) :: pi = 2.0d0 * acos(0.0d0) ,er0 = 10e-8
    omega = 2.0d0 / (1.0d0 + sin(pi/dble(n1 - 1)) )
    dx1 = 1.0d0 / dble(n1 - 1)
    dx2 = 1.0d0 / dble(n2 - 1)
    ! 係数c d の設定
    c = - ( dx2 * dx2 ) /( 2.0d0 * ( dx1 * dx1 + dx2 * dx2 ) )
    d = - ( dx1 * dx1 ) /( 2.0d0 * ( dx1 * dx1 + dx2 * dx2 ) )
    ! 格子点の設定
    do i = 1, n1
        x(1, i, 1:n2) = dble(i - 1) * dx1
    end do
    do i = 1, n2
        x(2, 1:n1, i) = dble(i - 1) * dx2
    end do
    ! ディレクレ境界条件の設定
    phi = 1.0d0
    phi(1, 1:n2) = 0.0d0
    phi(n1, 1:n2) = 0.0d0
    phi(1:n1, n2) = 0.0d0
    do i = 1, n1
        phi(i, 1) = sin(pi * i)
    end do
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
        write(*, *) 'itr, er =', itr, er
        if (er < er0) exit
    end do
end program main