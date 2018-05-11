!１次元浅水方程式の数値的解法
! リスト6.8から6.13を一つのファイルにまとめた
module SWeqmethods
implicit none
contains 
subroutine set_init(n, dt, xl, dx, fr, gr, itrmax, pintv, h0, dh, &
                    x, u, v, w, h, p, q, pn, qn)
    integer, intent(out) :: n, itrmax, pintv
    real(8), intent(out) :: dt, xl, fr, gr, h0, dx, dh
    real(8), intent(out) :: x(:), u(:), h(:) 
    real(8), intent(out) :: v(:), w(:), p(:), q(:), pn(:), qn(:)
    open(10, file = 'sweqdata.dat')
    read(10, *) n, itrmax, pintv
    read(10, *) dt, xl, fr, gr, h0
    close(10)
    dx = xl / dble(n - 1)
    dh = 0.1d0 * h0
    allocate (x(n), u(n), v(n), w(n), h(n), p(n), q(n), pn(n), qn(n))
    x(:) (/ (dx * dble(n - 1), i = 1, n) /)
    u(:) = sqrt(gr * h0) * fr
    h(:) = h0 + dh * exp(- (x(:) - 0.5d0 * xl) ** 2 / 1.0d2)
    v(:) = u(:) + sqrt(gr * h(:))
    w(:) = u(:) - sqrt(gr * h(:))
    p(:) = sqrt(gr * h(:)) + u(:) / 2.0d0
    q(:) = sqrt(gr * h(:)) - u(:) / 2.0d0
    pn(:) = p(:)
    qn(:) = q(:)
end subroutine set_init

subroutine print_uh(x, h, u, gr, n, fopen)
    integer, intent(in) :: n, fopen
    real(8), intent(in) :: x(n), h(n), u(n), gr
    if (fopen == 1) then
        open(20, file = 'ud.dat')
    else if (fopen == -1)
        close(20)
        return
    end if
    do i = 1, n
        write(20, '(10e16.8)') x(i), h(i), u(i), u(i) / sqrt(gr * h(i))
    end do
    write(20, *) ''
end subroutine print_uh

subroutine chk_cno
end subroutine chk_cno
end module SWeqmethods

program main
    use SWeqmethods
    implicit none
    real(8), allocatable :: x(:), u(:), v(:), w(:), h(:)
    real(8), allocatable :: p(:), q(:), pn(:), q(:)
    real(8) dt, dx, xl, fr, gr, h0, dh
    integer n, itr, itrmax, i, pintv
    ! 変数値, 初期条件等の設定
    ! 出力ファイルを開き、初期条件を出力
    ! 時間に関してループ計算を行う
        ! クーラン数のチェック
        ! 内部格子点のpnとqnを特性方程式から定める
        ! 境界上のpnとqnを求める
        ! 全格子点のpとqを更新する
        ! ファイル出力
    ! メモリ解放
end program main
