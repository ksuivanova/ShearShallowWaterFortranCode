b = 1000
reset
set terminal postscript enhanced eps  font 'Helvetica, 20'
#set style data linespoints
#set terminal postscript enhanced 
#set border linewidth 1.5
set grid
#set parametric
set nokey
set key left # bottom 
set xrange[0.0:1.01]
#set yrange[-0.01:0.05]  
#v= sprintf("%g",a)
set output './image/comparison1.eps'   
#set multiplot layout 3, 1 columnsfirst title " " 
set xlabel 'x/{/Symbol l}'
set ylabel "h/h_{0}"
plot './resu/1Droll.out' u ($1)/(1.3) + 0.01:($3)/(0.00798) every:::a::a  lt -1 lw 2 with lines title 'Numerical results', 'brockCase1_10052017.out' u 1:2 every:::a::a  lt -1 pt 7 ps 2  title "Brock's experimental results" #, './resu/convSHEAR1000p.out' u 1:5 every:::a::a  lt -1   pt 2 lw 0.3 title '1000 p', './resu/convSHEAR10000p.out' u 1:5 every:::a::a  lt -1   pt 7 ps 0.3 title '10000 p'
set ylabel "P_{12} (m^{2}/s^{2})"
#plot './resu/convSHEAR500p.out' u 1:7 every:::a::a  lt -1   pt 1 ps 0.4 lw 0.3 title '500 p', './resu/convSHEAR1000p.out' u 1:7 every:::a::a  lt -1   pt 2 ps 0.4 lw 0.3 title '1000 p','./resu/convSHEAR10000p.out' u 1:7 every:::a::a  lt -1   pt 7 ps 0.3 title '10000 p'
set ylabel "P_{22} (m^{2}/s^{2})"
#plot './resu/convSHEAR500p.out' u 1:8 every:::a::a lt -1 pt 1 ps 0.4 lw 0.3 title '500 p', './resu/convSHEAR1000p.out' u 1:8 every:::a::a lt -1 pt 2 ps 0.4 lw 0.3 title '1000 p', './resu/convSHEAR10000p.out' u 1:8 every:::a::a lt -1 pt 7 ps 0.3 lw 0.3 title '10000 p'
set ylabel "(P_{11}P_{22}-P_{12}^2)/h^2"
#plot './resu/convSHEAR500p.out'  u 1:6 every:::a::a  lt -1  pt 1 ps 0.4  lw 0.3  title '500 p', './resu/convSHEAR1000p.out' u 1:6 every:::a::a  lt -1  pt 2 ps 0.4  lw 0.3   title '1000 p', './resu/convSHEAR10000p.out' u 1:6 every:::a::a  lt -1   pt 7 ps 0.3 lw 0.3  title '10000 p'
set ylabel "h"
#plot './resu/convSHEAR500p.out'  u 1:3 every:::a::a  lt -1  pt 1 ps 0.4 lw 0.3  title '500 p', './resu/convSHEAR1000p.out' u 1:3 every:::a::a  lt -1  pt 2 ps 0.4 lw 0.3  title '1000 p', './resu/convSHEAR10000p.out' u 1:3 every:::a::a  lt -1   pt 7 ps 0.3 lw 0.3  title '10000 p'
set ylabel "u"
#plot './resu/convSHEAR500p.out'  u 1:4 every:::a::a  lt -1   pt 1 ps 0.4 lw 0.3  title '500 p', './resu/convSHEAR1000p.out' u 1:4 every:::a::a   lt -1 pt 2 ps 0.4 lw 0.3  title '1000 p', './resu/convSHEAR10000p.out' u 1:4 every:::a::a  lt -1   pt 7 ps 0.3  lw 0.3  title '10000 p'
#unset multiplot
pause -1
a=a+1
if(a<b) reread



#plot './resu/shear100p.out' u 1:5 every:::a::a  lt -1  lw 0.4 pt 1 ps 0.4  title '100 p', './resu/shear1000p.out' u 1:5 every:::a::a  lt -1 lw 0.3  pt 2 ps 0.4 title '1000 p' , './resu/shear10000p.out' u 1:5 every:::a::a  lt -1   pt 7 ps 0.3 title '10000 p'
set ylabel "P_{12} (m^{2}/s^{2})"
#plot './resu/shear100p.out' u 1:7 every:::a::a  lt -1   pt 1 ps 0.4 lw 0.3 title '100 p', './resu/shear1000p.out' u 1:7 every:::a::a  lt -1   pt 2 ps 0.4 lw 0.3 title '1000 p','./resu/shear10000p.out' u 1:7 every:::a::a  lt -1   pt 7 ps 0.3 title '10000 p'
set ylabel "P_{22} (m^{2}/s^{2})"
#plot './resu/shear100p.out' u 1:8 every:::a::a lt -1 pt 1 ps 0.4 lw 0.3 title '100 p', './resu/shear1000p.out' u 1:8 every:::a::a lt -1 pt 2 ps 0.4 lw 0.3 title '1000 p', './resu/shear10000p.out' u 1:8 every:::a::a lt -1 pt 7 ps 0.3 lw 0.3 title '10000 p'
set ylabel "(P_{11}P_{22}-P_{12}^2)/h^2"
#plot './resu/shear100p.out' u 1:6 every:::a::a  lt -1   lw 0.4 pt 1 ps 0.4 title '100 p', './resu/shear1000p.out' u 1:6 every:::a::a  lt -1 lw 0.3  pt 2 ps 0.4   title '1000 p', './resu/shear10000p.out' u 1:6 every:::a::a  lt -1   pt 7 ps 0.3 title '10000 p'
set ylabel "h"
#plot './resu/shear100p.out' u 1:3 every:::a::a  lt -1   lw 0.4 pt 1 ps 0.4 title '100 p', './resu/shear1000p.out' u 1:3 every:::a::a  lt -1 lw 0.3  pt 2 ps 0.4  title '1000 p', './resu/shear10000p.out' u 1:3 every:::a::a  lt -1   pt 7 ps 0.3 title '10000 p'
set ylabel "u"
#plot './resu/shear100p.out' u 1:4 every:::a::a  lt -1   lw 0.4 pt 1 ps 0.4 title '100 p', './resu/shear1000p.out' u 1:4 every:::a::a   lt -1 lw 0.3  pt 2 ps 0.4 title '1000 p', './resu/shear10000p.out' u 1:4 every:::a::a  lt -1   pt 7 ps 0.3 title '10000 p'