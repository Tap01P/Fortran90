! List06_03.f90
! Gauss-Seidel Method
!
module subprogs
  implicit none
contains
  subroutine Gauss_Seidel(a, b, x, n, itrmax, er0)
    ! a:係数行列, b:右辺列ベクトル, itrmax:最大反復回数, er0=誤差の閾値
    integer, intent(in) :: n, itrmax
    real(8), intent(in) :: a(n, n), b(n), er0
    real(8), intent(out) :: x(n)
    real(8) s, er, rd(n), r(n)
    integer  i, j, itr
    do i = 1, n
      if (a(i, i) == 0.0d0 ) stop 'a(i, i) == 0.0d0'
      rd(i) = 1.0d0 / a(i, i) !対角要素がノンゼロならばその逆数をrdとする
    enddo
    x(1:n) = 0.0d0 !初期解をゼロとする
    do itr = 1, itrmax !反復計算のループ（最大itrmax回の反復)
      do i = 1, n
        s = 0.0d0
        do j = 1, i-1
          s = s +  a(i, j) * x(j)
        enddo
        do j = i+1, n
          s = s + a(i, j) * x(j)
        enddo
        x(i) = rd(i) * (b(i) - s)
      end do
      r(1:n) = b(1:n) - matmul(a, x) !残差ベクトル
      er = dot_product(r, r)         !残差ベクトルの内積を誤差erとする
      write(*, *) 'itr = ', itr, 'err = ', er !途中経過の出力
      if (er <= er0) then
        write(*, *) '# converged #' !収束したことを表示
        exit !収束した場合は反復計算のループから抜ける
      endif
    enddo
end subroutine Gauss_Seidel

subroutine set_random_ab(a, b, n, x)
  integer n, i, j
  real(8), allocatable,intent(out) :: a(:, :), b(:), x(:)
  write(*, '(a)', advance = 'no') 'input n :'
  read(*, *) n
  if (n < 1 .or. n > 100) stop 'n must be 0 < n < 101'
  !計算時間の関係上、小規模の行列計算のみを行う
  allocate (a(n, n), b(n), x(n))
  !call random_seed
  !call random_number(a)
  !call random_seed
  call random_number(b)
  !call random_seed
  !call random_number(x)
  do i = 1, n
     do j = 1, n
       if (j == i) then
         a(i, i) = 1.0d0
       else
         a(i, j) = 1.0d0 / (dble(n) * dble(n))
       endif
     enddo
  enddo
end subroutine set_random_ab
end module subprogs

program main
  use subprogs
  implicit none
  real(8), allocatable :: a(:, :), b(:), x(:), r(:)
  integer :: i, n, itrmax = 100
  real(8) :: er0 = 1.0d-6
  call set_random_ab(a, b, n, x)!係数行列a, 右辺列ベクトルbの設定
  call Gauss_Seidel(a, b, x, n, itrmax, er0) !ガウス・ザイデル法
  !数値解xの出力及び誤差の表示
  do i = 1, n
  write(*, *) a(i, 1:n)
  enddo
  write(*, *) b(1:n)
  write(*, *) x(1:n)
  allocate(r(n))
  r(1:n) = b(1:n) - matmul(a, x)
  write(*, *) 'Gauss-Jordan error = ', dot_product(r, r)
  deallocate(a, b, x)
end program main
