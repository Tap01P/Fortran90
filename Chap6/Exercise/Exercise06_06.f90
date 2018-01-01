! 右辺ベクトルがm個ある場合のソルバー
! ひとまずピボット選択無しで作成
module solver
    implicit none
    contains
    function Gauss_Jordan(a0, b0, n) result(x)
    integer, intent(in) :: n
    real(8), intent(in) :: a0(n, n), b0(n, n)
    integer i, j
    real(8) ar, a(n, n), x(n, n)
    a(:, :) = a0(:, :)
    x(:, :) = b0(:, :)
    do i = 1, n
        ar = 1.0d0 / a(i, i)
        a(i, i) = 1.0d0
        a(i, i+1:n) = ar * a(i, i+1:n)
        x(i, :) = ar * x(i, :)
        do j = 1, n
            if(j /= i) then
                a(j, i+1:n) = a(j, i+1:n) - a(j, i) * a(i, i+1:n)
                x(j, :) = x(j, :) - a(j, i) * x(i, :)
                a(j, i) = 0.0d0
            end if 
        end do
    end do
    end function Gauss_Jordan

    subroutine set_random_matrix(a, b, n)
        integer, intent(in) :: n
        real(8), intent(out) :: a(n, n), b(n, n)
        call RANDOM_SEED
        call RANDOM_NUMBER(a)
        call RANDOM_NUMBER(b)
    end subroutine set_random_matrix

end module solver

program main
    use solver
    implicit none
    integer n
    real(8), allocatable :: a(:, :), b(:, :), x(:, :)
    write(*, '(a)', advance = 'no') 'input n:'
    read(*, *) n
    allocate(a(n, n), b(n, n), x(n, n))
    call set_random_matrix(a, b, n)
    x = Gauss_Jordan(a, b, n)
    write(*, *) b
    write(*, *) matmul(a, x)
end program main