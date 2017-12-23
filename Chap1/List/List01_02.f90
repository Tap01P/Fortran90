! implicit none 宣言をしないとえらい目に遭うことを体験しよう
!
program loop_err
  !implicit none
  !宣言をしない
  integer i, wa
  wa = 0
  do i = 1, 100
    wa = va + i
    !右辺のwaをvaと打ち間違えたと想定
  enddo

  write(*, *) 'wa =', wa
end program loop_err
!********************************
!wa = 100と出力されてしまう
!暗黙のうちに入っている変数vaが原因
!********************************
