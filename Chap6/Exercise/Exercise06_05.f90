! ガウス・ジョルダン法を用いて逆行列を求める
module solver
    implicit none
    contains
    function Gauss_Jordan(a0, I0, n) result(inverse)
    integer, intent(in) :: n
    real(8), intent(in) :: a0(n, n), I0(n, n)
    integer i, j
    real(8) ar, a(n, n), inverse(n, n)
    a(:, :) = a0(:, :)
    inverse(:, :) = I0(:, :)
    do i = 1, n
        ar = 1.0d0 / a(i, i)
        a(i, i) = 1.0d0
        a(i, i+1:n) = ar * a(i, i+1:n)
        inverse(i, :) = ar * inverse(i, :)
        do j = 1, n
            if(j /= i) then
                a(j, i+1:n) = a(j, i+1:n) - a(j, i) * a(i, i+1:n)
                inverse(j, :) = inverse(j, :) - a(j, i) * inverse(i, :)
                a(j, i) = 0.0d0
            end if 
        end do
    end do
    end function Gauss_Jordan

    subroutine set_random_matrix(a, n)
        integer, intent(in) :: n
        real(8), intent(out) :: a(n, n)
        call RANDOM_SEED
        call RANDOM_NUMBER(a)
    end subroutine set_random_matrix

end module solver

program main
    use solver
    implicit none
    integer n, k
    real(8), allocatable :: a(:, :), I(:, :), ainverse(:, :)
    write(*, '(a)', advance = 'no') 'input n:'
    read(*, *) n
    allocate(a(n, n), I(n, n), ainverse(n, n))
    call set_random_matrix(a, n)
    I(:, :) = 0.0d0
    do k = 1, n
        I(k, k) = 1.0d0
    end do
    ainverse = Gauss_Jordan(a, I, n)
    write(*, *) matmul(a, ainverse)
end program main
