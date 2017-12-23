! Bi-CGSTAB Method
module sbprogs
  implicit none
  contains
  subroutine Bi_CGSTAB(a, b, x, n, itrmax, er0)
    !
    integer, intent(in) :: n, itrmax
    real(8), intent(in) :: a(n, n), b(n), er0
    real(8), intent(out) :: x(n)
    integer itr
    real(8) alp, bet, c1, c2, c3, ev, vv, rr
    real(8) r(n), r0(n), p(n), y(n), e(n), v(n)
    ! setting the initial values 
    x(:) = 0.0d0
    r(:) = b - MATMUL(a, x)
    c1 = DOT_PRODUCT(r, r)
    if(c1 < er0) return
    p(:) = r(:)
    r0(:) = r(:)
    do itr = 1, itrmax
      y(:) = matmul(a, p)
      c2 = dot_product(r0, y)
      alp = c1 / c2
      e(:) = r(:) - alp * y(:)
      v(:) = matmul(a, e)
      ev = dot_product(e, v)
      vv = dot_product(v, v)
      c3 = ev / vv
      x(:) = x(:) + alp * p(:) + c3 * e(:)
      r(:) = e(:) - c3 * v(:)
      rr = dot_product(r, r)
      write(*, *) 'itr, er = ', itr, rr
      if(rr < er0) exit
      c1 = dot_product(r0, r)
      bet = c1 / (c2 * c3)
      p(:) = r(:) + bet * (p(:) - c3 * y(:))
    end do
  end subroutine Bi_CGSTAB
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

end module sbprogs

program main
  use sbprogs
  implicit none
  real(8), allocatable :: a(:, :), b(:), x(:), r(:)
  integer :: i, n, itrmax = 100
  real(8) :: er0 = 1.0d-6
  call set_random_ab(a, x, b, n)
  call Bi_CGSTAB(a, b, x, n, itrmax, er0) !ガウス・ザイデル法
  !数値解xの出力及び誤差の表示
  do i = 1, n
  write(*, *) a(i, 1:n)
  enddo
  write(*, *) b(1:n)
  write(*, *) x(1:n)
  allocate(r(n))
  r(1:n) = b(1:n) - matmul(a, x)
  write(*, *) 'error = ', dot_product(r, r)
  deallocate(a, b, x)
end program main
