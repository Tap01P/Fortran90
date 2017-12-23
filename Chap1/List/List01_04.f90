! 1から入力値までの和を計算するプログラム(無限ループを使用)
program loop_inf
    implicit none
    integer i, n, sum
    do
        write(*, '(a)', advance = 'no') 'input n (if n <= 0, stop):'
        read(*, *) n
        if (n <= 0) stop 'Good bye...'
        sum = 0
        do i = 1, n
            sum = sum + i
        end do
        write(*, *) 'The sum is ', sum
    end do
end program loop_inf