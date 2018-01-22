program main
    implicit none
    integer i
    do i = 1, 4, 2
        write(*, *) i
    end do
    write(*, *) i
    do i = 4, 1, -2
        write(*, *) i
    end do
    write(*, *) i
    do i = 0, 0, -1
        write(*, *) i
    end do
    write(*, *) i
end program main