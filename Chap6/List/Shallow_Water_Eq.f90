!１次元浅水方程式の数値的解法
! リスト6.8から6.13を一つのファイルにまとめた
module SWeqmethods
implicit none
contains 
subroutine set_init(n, dt, xl, dx, fr, gr, itrmax, pintv, h0, dh, &
                    x, u, v, w, h, p, q, pn, qn)
    integer i
    integer, intent(out) :: n, itrmax, pintv
    real(8), intent(out) :: dt, xl, fr, gr, h0, dx, dh
    real(8), allocatable, intent(out) :: x(:), u(:), h(:) 
    real(8), allocatable, intent(out) :: v(:), w(:), p(:), q(:), pn(:), qn(:)
    open(10, file = 'sweqdata')
    read(10, *) n, itrmax, pintv
    read(10, *) dt, xl, fr, gr, h0
    close(10)
    dx = xl / dble(n - 1)
    dh = 0.1d0 * h0
    allocate (x(n), u(n), v(n), w(n), h(n), p(n), q(n), pn(n), qn(n))
    x(:) = (/ (dx * dble(i - 1), i = 1, n) /)
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
    integer i
    integer, intent(in) :: n, fopen
    real(8), intent(in) :: x(n), h(n), u(n), gr
    if (fopen == 1) then ! fopen=1ならファイルをopenする
        open(20, file = 'ud.dat')
    else if (fopen == -1) then ! fopen=-1ならファイルをcloseしてreturn
        close(20)
        return
    end if
    do i = 1, n
        ! 各格子点上のx, h, u, フルード数の出力
        write(20, '(10e16.8)') x(i), h(i), u(i), u(i) / sqrt(gr * h(i))
    end do
    write(20, *) '' ! gnuplotによる描画のため空行を入れておく
end subroutine print_uh

subroutine chk_cno(n, dx, dt, v, w)
    ! クーラン数をチェックするサブルーチン
    integer, intent(in) :: n
    real(8), intent(in) :: dx, dt, v(n), w(n)
    real(8) cno1, cno2, cno
    cno1 = maxval(abs(v(:)) * dt / dx )
    cno2 = maxval(abs(w(:)) * dt / dx )
    cno = max(cno1, cno2)
    if (cno >= 1.0d0) then
        write(*, *) 'stop, cno >=1, cno=', cno
        stop
    end if
end subroutine chk_cno

subroutine cm1d(pnext, pprev, i, v, dt, dx, n)
    ! 内部格子点のpn, qnを特性方程式から求めるサブルーチン
    integer, intent(in) :: i, n
    real(8), intent(in) :: pprev(n), v(n), dt, dx
    real(8), intent(out) :: pnext(n)
    real(8) cno
    cno = v(i) * dt / dx ! クーラン数の算出
    if (cno >= 0.0d0) then ! 特性曲線の出発点位置に応じて空間内挿
        pnext(i) = pprev(i) - cno * (pprev(i) - pprev(i-1))
    else
        pnext(i) = pprev(i) - cno * (pprev(i+1) - pprev(i))
    end if
end subroutine cm1d

subroutine bc_thru(pn,qn, p, q, v, w, n, dt, dx)
    ! 境界のpn, qnを定めるサブルーチン
    integer, intent(in) :: n
    real(8), intent(in) :: p(n), q(n), v(n), w(n), dt, dx
    real(8), intent(out) :: pn(n), qn(n)
    integer i
    real(8) cno1, cno2
    do i = 1, n, n - 1 ! このループにおいて、iは1とnという値のみを取る
        pn(i) = p(i) ! 直前のステップをデフォルト値として設定する
        qn(i) = q(i) ! 同上
        cno1 = v(i) * dt / dx ! クーラン数の算出
        cno2 = w(i) * dt / dx ! 同上
        if (i == 1) then ! 上流端（特性速度が負ならpn,qnを計算:正ならデフォルト値を採用）
            if (cno1 < 0.0d0) pn(i) = p(i) - cno1 * (p(i+1) - p(i))
            if (cno2 < 0.0d0) qn(i) = q(i) - cno2 * (q(i+1) - q(i))
        else if (i == n) then ! 下流端(特性速度が正ならpn,qnを計算:負ならデフォルト値を採用)
            if (cno1 > 0.0d0) pn(i) = p(i) - cno1 * (p(i) - p(i-1))
            if (cno2 > 0.0d0) qn(i) = q(i) - cno2 * (q(i) - q(i-1))
        end if
    end do
end subroutine

subroutine pq2uhvw(pn, qn, gr, u, h, v, w, n)
    integer, intent(in) :: n
    real(8), intent(in) :: pn(n), qn(n), gr
    real(8), intent(out) :: u(n), h(n), v(n), w(n)
    integer i
    real(8) c
    do i = 1, n
        ! u, c, hをpn, qnから算出する
        u(i) = pn(i) - qn(i)
        c    = 0.5d0 * (pn(i) + qn(i))
        h(i) = c * c / gr
        ! v, wをu, cから算出する
        v(i) = u(i) + c
        w(i) = u(i) - c
    end do
end subroutine pq2uhvw

end module SWeqmethods

program main
    use SWeqmethods
    implicit none
    real(8), allocatable :: x(:), u(:), v(:), w(:), h(:)
    real(8), allocatable :: p(:), q(:), pn(:), qn(:)
    real(8) dt, dx, xl, fr, gr, h0, dh
    integer n, itr, itrmax, i, pintv
    ! 変数値, 初期条件等の設定
    call set_init(n, dt, xl, dx, fr, gr, itrmax, pintv, h0, dh, &
                  x, u, v, w, h, p, q, pn, qn)
    ! 出力ファイルを開き、初期条件を出力
    call print_uh(x, h, u, gr, n, 1)
    ! 時間に関してループ計算を行う
    do itr = 1, itrmax
        ! クーラン数のチェック
        call chk_cno(n, dx, dt, v, w)
        do i = 2, n - 1 ! 内部格子点のpnとqnを特性方程式から定める
            call cm1d(pn, p, i, v, dt, dx, n)
            call cm1d(qn, q, i, w, dt, dx, n)
        end do
        ! 境界上のpnとqnを求める
        call bc_thru(pn, qn, p, q, v, w, n, dt, dx)
        call pq2uhvw(pn, qn, gr, u, h, v, w, n)
        ! 全格子点のpとqを更新する
        p(:) = pn(:)
        q(:) = qn(:)
        ! ファイル出力
        if (mod(itr, pintv) == 0) call print_uh(x, h, u, gr, n, 0)
    end do
    ! 出力ファイルをcloseする
    call print_uh(x, h, u, gr, n, -1)
    ! メモリ解放
    deallocate(x, u, v, w, h, p, q, pn, qn)
    ! pltファイルの作成
    open(30, file='sweqplot.plt')
    write(30, *) 'reset'
    write(30, *) 'plot "ud.dat" with lines'
    write(30, *) 'pause -1'
    close(30)
    call system('gnuplot sweqplot.plt')
end program main
