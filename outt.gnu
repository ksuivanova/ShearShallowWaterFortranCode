b = 1000
reset
#set terminal png
#set terminal postscript enhanced eps  font 'Helvetica,16'
 #set style data lines
set grid
set key bottom right
#set xrange[0.06:1]
#set yrange[0.0:0.5]
#v= sprintf("%g",a+1)
#set output './image/ConvOmega1WG_'.v.'.eps' #'roll1D'.v.'.eps'
#set multiplot layout 2, 1 columnsfirst title " " 
set xlabel 'r (m)'
set ylabel 'p22 '
plot  './resu/test_0000.out' u 1:3 every :::a::a  lt -1 pt 8 ps 0.7 title 'numerical, 100 p', './resu/test_0001.out' u 1:3 every :::a::a  lt -1 lw 3 title 'analytical' with lines
#set ylabel 'U_{/Symbol q} (m/s) '
#plot  './resu/waterNum100p_0000.out' u 1:4 every :::a::a  lt -1 pt 8 ps 0.7 title 'numerical, 100 p', './resu/ANALITYCAL_SOLUTION_WATER_GLASS_0000.out' u 1:5 every :::a::a lt -1 lw 3 title 'analytical' with lines
#unset multiplot
pause -1
a=a+1
if(a<b) reread


#plot  './resu/err_NonStatSol_0000.out' u 1:2 every :::a::a lt -1 
#set ylabel 'u '
#plot  './resu/StSol_0000.out' u 1:3 every :::a::a lt -1 
#set ylabel 'p11 '
#plot  './resu/StSol_0000.out' u 1:4 every :::a::a lt -1 
#unset multiplot
#pause -1
#a=a+1
#if(a<b) reread

