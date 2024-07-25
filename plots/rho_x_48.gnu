set terminal epslatex size 3.5,2.62 standalone color colortext 10
set output "rho_x_48.tex"
set key right Right

set multiplot

set xlabel '$x$'
set ylabel '$\rho$'

#set xrange [0:10]

plot '../data/rho_x_48.dat' using 1:2 w l lw 3 dt 1 title '$\mathrm{Re}(\rho)$',\
 '../data/rho_x_48.dat' using 1:3 w l lw 3 dt 2 title '$\mathrm{Im}(\rho)$'

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
