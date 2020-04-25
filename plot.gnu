reset
set style data lines
set grid
#set key left
set nokey
set xlabel 'x (m)'
set ylabel 'profondeur (m)'
plot './resu/double_shock_000.out' u 1:3 lt -1 lw 2
pause -1
set ylabel 'u (m/s)'
plot './resu/double_shock_000.out' u 1:4 lt -1 lw 2
pause -1
set ylabel 'v (m/s)'
plot './resu/double_shock_000.out' u 1:5 lt -1 lw 2
pause -1
set ylabel 'p11 (m^{2}/s^{2})'
plot './resu/double_shock_000.out' u 1:6 lt -1 lw 2
pause -1
set ylabel 'p12 (m^{2}/s^{2})'
plot './resu/double_shock_000.out' u 1:7 lt -1 lw 2
pause -1
set ylabel 'p22 (m^{2}/s^{2})'
plot './resu/double_shock_000.out' u 1:8 lt -1 lw 2
pause -1


