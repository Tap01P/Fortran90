!List01_02.f90でimplicit none宣言をした場合
program loop_err
  implicit none
  integer i, wa
  wa  = 0
  do i = 1, 100
    wa = va + i
  end do
  write(*, *) 'wa = ', wa
end program loop_err
!************************************************
!明示的に宣言していない変数vaがあるため、
!以下の様なエラーが出る
!ensyu01_01.f90:7:11:
!
!     wa = va + i
!           1
!Error: Symbol 'va' at (1) has no IMPLICIT type
!************************************************
