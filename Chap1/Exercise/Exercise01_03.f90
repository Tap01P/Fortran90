!配列を用いる
program main
implicit none
    integer i, a(10)
    a(1) = 1
    a(2) = 2
    do i = 3 , 10
    a(i) = a(i - 1) + a(i - 2)
    end do
    write(*, *) a(:)
end program main