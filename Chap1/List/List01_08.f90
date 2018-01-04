! ニュートン法(1変数ver.)
! x^2 == aの正の解 sqrt(a) を求める
program Newthon_Method
implicit none
    integer :: itr, itrmax = 100 
    real(8) :: x1, x2, a, er, er0 = 1.0d-6
    write(*, '(a)', advance = 'no') 'input a:'
    read(*, *) a
    if (a <= 0.0d0) stop 'a <= 0.0d0'
    x1 = a
    do itr = 1, itrmax
        x2 = x1 - 0.5d0 * (x1 * x1 - a) / x1
        er = abs(x2 - x1)
        if (er < er0) exit
        x1 = x2
    end do
    write(*, *) 'sol, itr, er', x2, itr, er
end program Newthon_Method