! モンテカルロ法によって円周率piを計算する
program monte
    implicit none
    real(8) x, y, pi, pi0
    integer :: n, i, im = 2 ** 20
    pi0 = 2.0d0 * acos(0.0d0)
    n = 0
    do i = 1, im
        call random_number(x)
        call random_number(y)
        if (x * x + y * y <= 1.0d0) n = n + 1
    end do
    pi = 4.0d0 * dble(n) / dble(im)
    write(*, *) 'pi, pi0, er=', pi, pi0, pi - pi0
end program monte