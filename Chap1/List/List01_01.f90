! 1から100までの自然数の総和をdoループで求めるプログラム
program loop !プログラムの名称を指定する
  implicit none !暗黙の型変換の宣言を無効化する
  integer i, wa !整数型変数iとwaを宣言
  wa = 0 !合計を表す変数waの初期値をゼロにする（初期化）
  do i = 1, 100 !doループによる反復演算
      wa = wa + i !waにiを加算する
  enddo !doループの終端

  write(*, *) 'wa = ', wa !結果の出力
end program loop !end文で終了する
