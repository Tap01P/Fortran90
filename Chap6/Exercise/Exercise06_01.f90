!ガウスの消去法によって連立一次方程式の解を求める
!行列の各要素は乱数を格納する
!
module subprogs
  implicit none
contains
  subroutine gauss_shokyo(a0, x, b0, n)
    integer, intent(in) :: n
    real(8), intent(in) :: a0(n, n), b0(n)
    real(8), intent(out) :: x(n)
    integer i, k, j
    real(8) ar, a(n, n), b(n)
    a(:, :) = a0(:, :)
    b(:)    = b0(:)
    do k = 1, n
      if ( a(k, k) == 0.0d0) stop 'pivot = 0'
      ar = 1.0d0 / a(k, k)
      a(k, k) = 1.0d0
      a(k, k+1:n) = ar * a(k, k+1:n)
      b(k) = ar * b(k)
      do i = k+1, n
        a(i, k+1:n) = a(i, k+1:n) - a(i, k) * a(k, k+1:n)
        b(i)        = b(i)        - a(i, k) * b(k)
        a(i, k ) = 0.0d0
      enddo
    enddo
    x(n) = b(n)
    do k = n-1, 1, -1
      x(k) = b(k)
      do i = k+1, n
        x(k) = x(k) - a(k, i) * x(i)
      end do
    enddo
    do i = 1, n
      write(*, *) a0(i, 1:n)
    enddo
    do i = 1, n
      write(*, *) a(i, 1:n)
    end do
    do i = 1, n
      write(*, *) b0(i)
    end do
    write(*, *) matmul(a0, x)
  end subroutine gauss_shokyo

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
  call gauss_shokyo(a,x,b,n)
  allocate(r(n))
  r(:) = b(:) - matmul(a, x)
  write(*, *) 'error = ', dot_product(r, r)
  deallocate(a, b, x)
end program main
