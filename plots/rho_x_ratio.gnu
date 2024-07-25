set terminal epslatex size 4,2.47 standalone color colortext 8
set output "rho_x_ratio_b5.tex"
set key right Right

set multiplot

set xlabel '$x$'
set ylabel '$\rho$'

#set xrange [0:10]

plot '../data/rho_x_exact_a1_b5.dat' using 1:2 w l lw 3 dt 1 title '$\mathrm{Re}(\rho_{\mathrm{ex}})$',\
 '../data/rho_x_exact_a1_b5.dat' using 1:3 w l lw 3 dt 2 title '$\mathrm{Im}(\rho_{\mathrm{ex}})$', \
 '../data/rho_x_effective_a1_b5.dat' using 1:2 w l lw 3 dt 3 lt -1 title '$\mathrm{Re}(\rho_{\mathrm{eff}})$', \
 '../data/rho_x_effective_a1_b5.dat' using 1:3 w l lw 3 dt 4 title '$\mathrm{Im}(\rho_{\mathrm{eff}})$'

set nokey

set xlabel '$\beta$'
set ylabel '$\rho_{\mathrm{eff}}(0)/\rho_{\mathrm{ex}}(0)$'

set xrange [2:10]
set yrange [1:1.05]
set origin 0.35, 0.35
set size 0.6, 0.4
set xtics 1
set ytics 0.05

plot '../data/ratio.dat' w l lw 3 dt 1 lt 7

replot 
unset multiplot

# instruction to generate .pdf file from epslatex terminal: 
# gnuplot file.gnu
# pdflatex file.tex
