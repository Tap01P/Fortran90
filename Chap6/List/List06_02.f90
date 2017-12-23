!------------------------------------------------------------------------------
! Institution, Affiliation
!------------------------------------------------------------------------------
!
! MODULE:  Gauss_Jordan_pivoting
!
!> @author
!> Author Name}
!
! DESCRIPTION: 
!>  Gauss-Jordan method with partial pivoting
!
! REVISION HISTORY:
! 26 Nov 2017 - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
module subprogs
implicit none
contains
    subroutine Gauss_Jordan_pivoting(a0, x, b, n)
        ! Gauss-Jordan Method with partinal povoting
        integer, intent(in) :: n ! size of matrices
        real(8), intent(in) :: a0(n, n), b(n) 
        ! a0: coeficient mat, b: vector in rhs
        real(8), intent(out) :: x(n)
        ! x : variable vector
        integer i, j, k, m ! for loop
        real(8) ar, am, t, a(n, n), w(n)
        a(:, :) = a0(:, :)
        x(:) = b(:)
        do k = 1, n
            m = k
            am = abs(a(k, k))
            do i = k+1, n 
            ! serching m th row in which a(i, k) is the maximum
                if (abs(a(i, k)) > am) then
                    am = abs(a(i,k))
                    m = i
                end if
            end do
            if (am == 0.0d0) stop 'A is singular'
            ! stopping when A is a singular matrix
            if( k /= m) then
            ! swap k th row and m th row
                w(   k:n) = a(k, k:n)
                a(k, k:n) = a(m, k:n)
                a(m, k:n) = w(   k:n)
                t = x(k)
                x(k) = x(m)
                x(m) = t
            end if
            ! calculation with Gauss-Jordam method
            ar = 1.0d0 / a(k, k)
            a(k, k) = 1.0d0
            a(k, k+1:n) = ar * a(k, k+1:n)
            x(k) = ar * x(k)
            do i = 1, n
                if (i /= k) then
                    a(i, k+1:n) = a(i, k+1:n) - a(i, k) * a(k, k+1:n)
                    x(i)        = x(i)        - a(i, k) * x(k)
                    a(i, k    ) = 0.0d0
                end if
            end do
        end do


    end subroutine Gauss_Jordan_pivoting 

    subroutine set_random_ab(a, x, b, n)
        real(8), allocatable, intent(out) :: a(:, :), x(:), b(:)
        integer n
        write(*, '(a)', advance = 'no') 'input n :'
        read(*, *) n
        if( n < 1 .or. n > 1001) stop 'n must be 0 < n < 1001'
        allocate (a(n, n), b(n), x(n))
        call random_seed
        call random_number(a)
        call random_number(b)
   end subroutine set_random_ab
end module subprogs

program main
    use subprogs
    implicit none
    real(8), allocatable :: a(:, :), x(:), b(:), r(:)
    integer i, n
    call set_random_ab(a,x,b,n)
    call Gauss_Jordan_pivoting(a,x,b,n)
    allocate(r(n))
    r(:) = b(:) - matmul(a, x)
    write(*, *) 'Gauss-Jordan error = ', dot_product(r, r)
    deallocate(a, b, x)
end program main
