! 偶数の和と奇数の和を求めるプログラム
program loop_odd_even
    implicit none
    integer i, n, s_even, s_odd
    s_even = 0
    s_odd = 0
    read (*, *) n
    do i = 1, n
        if (mod(i, 2) == 0) then
            s_even = s_even + i
        else if (mod(i, 2) == 1) then
            s_odd = s_odd + i
        else
            stop 'something is wrong !!'
        end if
    end do
    write(*, *) 'Sum of even is', s_even
    write(*, *) 'Sum of odd is', s_odd
end program loop_odd_even