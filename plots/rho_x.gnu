set terminal epslatex size 4,2.47 standalone color colortext 10
set output "rho_x_a0.tex"
set key right Right

set multiplot

set xlabel '$x$'
set ylabel '$\rho_{_T}$'

set xrange [0:7]

#set title '$\alpha=1, \, \beta = 2$'


plot '../data/rho_x_exact_a0_b1.dat' using 1:2 w l lw 5 dt 1 lc 4 title '$\rho_{\mathrm{ex}}, \, \beta = 1$',\
'../data/rho_x_effective_a0_b1.dat' using 1:2 w l lw 5 dt 2 lc 7 title '$\rho_{\mathrm{eff}}, \, \beta = 1$', \
'../data/rho_x_exact_a0_b2.dat' using 1:2 w l lw 5 dt 1 title  '$\rho_{\mathrm{ex}}, \, \beta = 2$',\
'../data/rho_x_effective_a0_b2.dat' using 1:2 w l lw 5 dt 2 lc 9 title '$\rho_{\mathrm{eff}}, \, \beta =2$',\
'../data/rho_x_exact_a0_b5.dat' using 1:2 w l lw 5 dt 1 title  '$\rho_{\mathrm{ex}}, \, \beta = 5$',\
'../data/rho_x_effective_a0_b5.dat' using 1:2 w l lw 5 dt 2 title '$\rho_{\mathrm{eff}}, \, \beta =5$',\
'../data/rho_x_exact_a0_b10.dat' using 1:2 w l lw 5 dt 1 title  '$\rho_{\mathrm{ex}}, \, \beta = 10$',\
'../data/rho_x_effective_a0_b10.dat' using 1:2 w l lw 5 dt 2 title '$\rho_{\mathrm{eff}}, \, \beta =10$'





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
