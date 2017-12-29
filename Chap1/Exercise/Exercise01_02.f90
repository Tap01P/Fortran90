program main
implicit none
    integer i, n, s
    integer :: wa = 0, wasq = 0, wacu = 0
    write(*, *) 'input positive number n :'
    read(*, *) n
    do i = 1, n
        wa = wa + i
        wasq = wasq + i ** 2
        wacu = wacu + i ** 3
    end do
    write(*, *) 'wa =', wa
    s = n * (n + 1) / 2
    write(*, *) 's =', s
    write(*, *) 'wasq =', wasq
    s = n * (n + 1) * (2*n +1) / 6
    write(*, *) 's =', s
    write(*, *) 'wacu =', wacu
    s = n * n * (n + 1) * (n + 1) / 4
    write(*, *) 's =', s
end program main