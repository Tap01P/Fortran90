! FileName : List06_07.f90
!************************************************
!
! Program : Solution plot
!
! Purpose : Solving for Difusion Equations
! 2D system, 1.0 * 1.0 square 
!************************************************

module Setting_IandO
implicit none
contains
    subroutine set_dbc(phi, x, n1, n2)
    ! Setting Dirichlet boudary condition
        integer, intent(in) :: n1, n2
        ! The lattice number in the system (n1 + 1) * (n2 + 1)
        integer i
        real(8), intent(in) :: x(2, n1, n2)
        real(8), intent(out) :: phi(n1, n2)
        real(8) :: pi = 2.0d0 * acos(0.0d0)
        do i = 1, n1
            phi(i, 1) = sin(pi * x(1, i, 1))
        end do
        phi(1, 2:n2) = 0.0d0
        phi(1:n1, n2) = 0.0d0
        phi(n1, 2:n2) = 0.0d0
    end subroutine set_dbc

    subroutine set_nbc(phi, x, n1, n2)
    ! setting Neumann boundary condition
        integer, intent(in) :: n1, n2
        integer i
        real(8), intent(in) :: x(2, n1, n2)
        real(8), intent(out) :: phi(n1, n2)
        real(8) :: pi = 2.0 * acos(0.0d0)
    phi(1, 1:n2) = 0.0d0
    phi(1:n1, 1) = 0.0d0
    phi(n1, 1:n2) = 0.0d0
    phi(:, n2) = phi(:, n2 - 1)
    end subroutine set_nbc

    subroutine output(phi, x, n1, n2)
    ! output results
        integer, intent(in) :: n1, n2
        integer i, j
        real(8), intent(in) :: phi(n1, n2), x(2, n1, n2)
        !write(a, *) fo
        !open(fo, file = 'result.dat')
        do i = 1, n1
            do j = 1, n2
                write(*, *) x(1, i, j), x(2, i, j), phi(i, j)
            end do 
        end do
    end subroutine output

    function check_steady(phi, phi2, n1, n2) result(er)
        integer, intent(in) :: n1, n2
        integer i, j
        real(8), intent(in) :: phi(n1, n2), phi2(n1, n2)
        real(8) r(n1, n2), er
        r(1:n1-1, 1:n2-1) = 0.0d0
        er = 0.0d0
        do j = 2, n2 - 1
            do i = 2, n1 - 1
                r(i, j) = phi(n1, n2) - phi2(n1, n2)
                er = er + r(i, j) * r(i, j)
            end do
        end do
    end function check_steady
end module Setting_IandO


program main
use Setting_IandO
implicit none
    integer n1, n2, nstep, pstep, i, j, istep
    real(8) dt, dn1, dn2, er, er0, d1, d2
    real(8), allocatable :: x(:, :, :), phi(:, :), phi2(:, :)
    write(*, '(a)') 'input n1, n2, nstep, pstep :'
    read(*, *) n1, n2, nstep, pstep
    dn1 = 1.0d0 / dble(n1 - 1)
    dn2 = 1.0d0 / dble(n2 - 1)
    allocate(x(2, n1, n2))
    allocate(phi(n1, n2), phi2(n1, n2))
    do i = 1, n1
        x(1, i, 1:n2) = dble(i - 1) * dn1
    end do
    do j = 1, n2
        x(2, 1:n1, j) = dble(j - 1) * dn2
    end do
    er0 = 1.0e-6
    dt = 5.0e-6
    d1 = dt / (dn1 * dn1)
    d2 = dt / (dn2 * dn2)
    phi(2:n1 - 1, 2:n2 - 1) = 0.0d0
    call set_dbc(phi, x, n1, n2)
    do istep = 1, nstep
        do j = 2, n2 - 1
            do i = 2, n1 - 1
                phi2(i, j) = phi(i, j) &
                            + d1 * (phi(i-1, j) - 2.0d0 * phi(i, j) + phi(i+1, j)) &
                            + d2 * (phi(i, j-1) - 2.0d0 * phi(i, j) * phi(i, j+1))
            end do
        end do
        call set_nbc(phi2, x, n1, n2)
        er = check_steady(phi, phi2, n1, n2)
        phi(:, :) = phi2(:, :)
        if(mod(istep, pstep) == 0) call output(phi, x, n1, n2)
        !if(er < er0) exit
    end do

end program main