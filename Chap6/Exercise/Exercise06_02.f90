! ガウス・ジョルダン法とガウスの消去法の計算時間を比較する
! ピボット有り
! サブルーチンではなく関数でソルバーを記述してみた
module solvers
implicit none
contains
    function Gauss_Jordan(a, b, n) result(x)
    integer, intent(in) :: n
    real(8), intent(in) :: a(n, n), b(n)
    integer i, j, m
    real(8) adiagmax, t, ar, adummy(n, n), x(n), w(n)
    adummy(:, :) = a(:, :)
    x(:) = b(:)
    do i = 1, n
        m = i
        adiagmax = abs(a(i, i))
        do j = i+1, n
            if (abs(a( j, i)) > adiagmax) then
                m = j
                adiagmax = a(j, i)
            end if
        end do
        if ( adiagmax < 1.0e-15) stop 'A is singular'
        ! stopping when A is a singular matrix
        if (i /= m) then
        ! swapping k the row and m the row
            w(   i:n) = adummy(i, i:n)
            adummy(i, i:n) = adummy(m, i:n)
            adummy(m, i:n) = w(   i:n)
            t    = x(i)
            x(i) = x(m)
            x(m) = t  
        end if
        ! Gauuss-Jordan Method
        ar = 1.0d0 / adummy(i, i)
        adummy(i, i) = 1.0d0
        adummy(i, i+1:n) = ar * adummy(i, i+1:n)
        x(i) = ar * x(i)
        do j = 1, n
            if (j /= i) then
                adummy(j, i+1:n) = adummy(j, i+1:n) - adummy(j, i) * adummy(i, i+1:n)
                x(j) = x(j) - adummy(j, i) * x(i)
                adummy(j,     i) = 0.0d0
            end if
        end do
    end do
    end function Gauss_Jordan

    function Gaussian_Elimination(a, b, n) result(x)
        integer, intent(in) :: n
        real(8), intent(in) :: a(n, n), b(n)
        integer i, j, m
    real(8) adiagmax, t, ar, adummy(n, n), bdummy(n), x(n), w(n)
    adummy(:, :) = a(:, :)
    bdummy(:) = b(:)
    do i = 1, n
        m = i
        adiagmax = abs(a(i, i))
        do j = i+1, n
            if (abs(a( j, i)) > adiagmax) then
                m = j
                adiagmax = a(j, i)
            end if
        end do
        if ( adiagmax < 1.0e-15) stop 'A is singular'
        ! stopping when A is a singular matrix
        if (i /= m) then
        ! swapping k the row and m the row
            w(   i:n) = adummy(i, i:n)
            adummy(i, i:n) = adummy(m, i:n)
            adummy(m, i:n) = w(   i:n)
            t    = bdummy(i)
            bdummy(i) = bdummy(m)
            bdummy(m) = t  
        end if
        ! Gauussian Elimination Method
        ar = 1.0d0 / adummy(i, i)
        adummy(i, i) = 1.0d0
        adummy(i, i+1:n) = ar * adummy(i, i+1:n)
        bdummy(i) = ar * bdummy(i)
        do j = i+1, n
            adummy(j, i+1:n) = adummy(j, i+1:n) - adummy(j, i) * adummy(i, i+1:n)
            bdummy(j) = bdummy(j) - adummy(j, i) * bdummy(i)
            adummy(j,     i) = 0.0d0
        end do
    end do
    x(n) = bdummy(n)
    do i = n-1, 1, -1
        x(i) = bdummy(i)
        do j = i+1, n
            x(i) = x(i) - adummy(i, j) * x(j)  
        end do 
    end do
    end function Gaussian_Elimination

    subroutine set_random_matrix(a, b, n)
        integer, intent(in) :: n
        real(8), intent(out) :: a(n, n), b(n)
        call RANDOM_SEED
        call RANDOM_NUMBER(a)
        call RANDOM_NUMBER(b)
    end subroutine set_random_matrix
end module solvers

program main
use solvers
implicit none
    integer n
    real(8) tstart, tgj, tge, tre
    real(8), allocatable :: a(:, :), b(:), x(:)
    write(*, '(a)', advance = 'no') 'input n:'
    read(*, *) n
    allocate(a(n, n), b(n), x(n))
    call set_random_matrix(a, b, n)
    !write(*, *) 'vector : ', b(:)
    call CPU_TIME(tstart)
    x(:) = Gauss_Jordan(a, b, n)
    call CPU_TIME(tgj)
    !write(*, *) 'solution GJ:', x(:)
    !write(*, *) 'comfirmation:', matmul(a, x)
    call CPU_TIME(tre)
    x(:) = Gaussian_Elimination(a, b, n)
    call CPU_TIME(tge)
    !write(*, *) 'solution Ge:', x(:)
    !write(*, *) 'comfirmation', matmul(a, x)
    write(*,*) tgj - tstart
    write(*,*) tge - tre
end program main