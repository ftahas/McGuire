set terminal epslatex size 4,2.47 standalone color colortext 8
set output "rho_x_lambda_a1_b1.tex"
set key right Right

#set multiplot

set xlabel '$x$'
set ylabel '$\rho(x,\Lambda)$'

set xrange [0:20]

#set title '$\alpha=1, \; \beta = 1$'

#plot '../data/rho_x_lambda05_a1_b1_exact.dat' using 1:2 w l lw 3 title '$\mathrm{Re}\rho_{\mathrm{ex}}, \Lambda = 0.5$', \
#'../data/rho_x_lambda05_a1_b1_effective.dat' using 1:2 w p ps .4 pt 2 title '$\mathrm{Re}\rho_{\mathrm{eff}}, \Lambda = 0.5$', \
#'../data/rho_x_lambda08_a1_b1_exact.dat' using 1:2 w l lw 3 title '$\mathrm{Re}\rho_{\mathrm{ex}}, \Lambda = 0.8$', \
#'../data/rho_x_lambda08_a1_b1_effective.dat' using 1:2 w p ps .4 pt 2 title '$\mathrm{Re}\rho_{\mathrm{eff}}, \Lambda = 0.8$', \
#'../data/rho_x_lambda09_a1_b1_exact.dat' using 1:2 w l lw 3 title '$\mathrm{Re}\rho_{\mathrm{ex}}, \Lambda = 0.9$', \
#'../data/rho_x_lambda09_a1_b1_effective.dat' using 1:2 w p ps .4 pt 2 title '$\mathrm{Re}\rho_{\mathrm{eff}}, \Lambda = 0.9$', \
#'../data/rho_x_lambda11_a1_b1_exact.dat' using 1:2 w l lw 3 title '$\mathrm{Re}\rho_{\mathrm{ex}}, \Lambda = 1.1$', \
#'../data/rho_x_lambda11_a1_b1_effective.dat' using 1:2 w p ps .4 pt 2 title '$\mathrm{Re}\rho_{\mathrm{eff}}, \Lambda = 1.1$', \
#'../data/rho_x_lambda12_a1_b1_exact.dat' using 1:2 w l lw 3 title '$\mathrm{Re}\rho_{\mathrm{ex}}, \Lambda = 1.2$', \
#'../data/rho_x_lambda12_a1_b1_effective.dat' using 1:2 w p ps .4 pt 2 title '$\mathrm{Re}\rho_{\mathrm{eff}}, \Lambda = 1.2$'



plot '../data/rho_x_lambda05_a1_b1_exact.dat' using 1:2 w l lw 3 title '$\Lambda = 0.5, \, \mathrm{Re}\rho_{\mathrm{ex}}$', \
'../data/rho_x_lambda05_a1_b1_exact.dat' using 1:3 w l lw 3 dt 1 title '$\mathrm{Im}\rho_{\mathrm{ex}}$', \
'../data/rho_x_lambda05_a1_b1_effective.dat' using 1:2 w p ps .4 pt 2  title '$\mathrm{Re}\rho_{\mathrm{eff}}$', \
'../data/rho_x_lambda05_a1_b1_effective.dat' using 1:3 w p ps .4 pt 2 title '$\mathrm{Im}\rho_{\mathrm{eff}}$', \
'../data/rho_x_lambda2_a1_b1_exact.dat' using 1:2 w l lw 3 lt -1 title '$\Lambda = 2,\, \mathrm{Re}\rho_{\mathrm{ex}}$', \
'../data/rho_x_lambda2_a1_b1_exact.dat' using 1:3 w l lw 3 dt 1 title '$\mathrm{Im}\rho_{\mathrm{ex}}$', \
'../data/rho_x_lambda2_a1_b1_effective.dat' using 1:2 w p ps .4 pt 2 title '$\mathrm{Re}\rho_{\mathrm{eff}}$', \
'../data/rho_x_lambda2_a1_b1_effective.dat' using 1:3 w p ps .4 pt 2 title '$\mathrm{Im}\rho_{\mathrm{eff}}$', \


#set nokey
#unset xlabel
#unset ylabel
#unset title

#set logscale y
#set origin 0.35, 0.3
#set size 0.6, 0.4
#set xrange [0:20]
#set yrange [-1:1]
#set xtics 40
#set ytics 100
#replot

#unset multiplot


# instruction to generate .pdf file from epslatex terminal: 
# gnuplot file.gnu
# pdflatex file.tex
