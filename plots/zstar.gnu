set terminal epslatex size 4,2.47 standalone color colortext 10
set output "zstar.tex"
set key at -.8,1.3

set multiplot

set xlabel '$\Lambda/\Lambda_c$'
set ylabel '$z_*$'

f(x) = -1.46576
g(x) = 1.46576
h(x) = -1.09868

set xrange [-3:3]
set yrange [-1.6:1.6]

set xzeroaxis
set yzeroaxis

plot '../data/zstar.dat' using 1:2 w l lw 5 dt 1 lc 1 lt -1 title '$\mathrm{Re} z_*$',\
 '../data/zstar.dat' using 1:3 w l lw 5 dt 4 title '$\mathrm{Im} z_*$',\
 f(x) w l lw 5 dt 3 lc 6 notitle, g(x) w l lw 5  dt 3 notitle, h(x) w l lw 5 dt 3 lc 7 notitle

#set nokey
#unset xlabel
#unset ylabel
#unset title

#set logscale y
#set origin 0.35, 0.3
#set size 0.6, 0.4
#set xrange [0:120]
#set yrange [-1:1]
#set xtics 5
#set ytics 100
#replot

#unset multiplot


# instruction to generate .pdf file from epslatex terminal: 
# gnuplot file.gnu
# pdflatex file.tex
