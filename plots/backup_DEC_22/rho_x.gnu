set terminal epslatex size 4,2.47 standalone color colortext 10
set output "rho_x_a1_b1.tex"
set key right Right

set multiplot

set xlabel '$x$'
set ylabel '$\rho$'

#set xrange [0:10]

plot '../data/rho_x_exact_a1_b1.dat' using 1:2 w l lw 5 dt 1 title '$\mathrm{Re}(\rho_{\mathrm{ex}})$',\
 '../data/rho_x_exact_a1_b1.dat' using 1:3 w l lw 5 dt 2 title '$\mathrm{Im}(\rho_{\mathrm{ex}})$', \
 '../data/rho_x_effective_a1_b1.dat' using 1:2 w l lw 5 dt 3 lt -1 title '$\mathrm{Re}(\rho_{\mathrm{eff}})$', \
 '../data/rho_x_effective_a1_b1.dat' using 1:3 w l lw 5 dt 4 title '$\mathrm{Im}(\rho_{\mathrm{eff}})$'

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
