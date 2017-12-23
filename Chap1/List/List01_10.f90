! 入出力ファイルを利用するプログラム
program file_io
    implicit none
    real(8) d, x, y, z
    integer :: n, i, j, fi = 10, fo = 11 
    ! 初期値設定する場合は::を付ける
    open(fi, file = 'input.dat')  ! 入力ファイルinput.dを開く
    open(fo, file = 'output.dat') ! 出力ファイルoutput.dを開く
    read(fi, *) n               ! 入力ファイルから読み込んだ数値をnに格納する
    close(fi)
    if (n < 3) stop 'stop, n < 3'
    d = 10.0d0 / dble(n - 1)
    do j = 1, n
        y = -5.0d0 + dble(j - 1) * d
        do i = 1, n
        x = -5.0d0 + dble(i - 1) * d
        z = sin(x) * sin(y)
        write(fo, '(3e12.4)') x, y, z
        end do
        write(fo, *) '' ! gnuplot描画の際に直線が引かれないように空行を入れる
    end do
    close(fo) ! 出力ファイルを閉じる
    !----------Gnuplotを用いて描画する----------------------------------
    open(13, file = 'gnuptest01.plt')
        write(13, *) 'reset'
        write(13, *) "set hidden3d"
        write(13, *) "set contour"
        write(13, *) "splot 'output.dat' with lines"
        write(13, *) "pause -1"
    close(13)
    call system("gnuplot gnuptest01.plt")
end program file_io