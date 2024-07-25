set terminal epslatex size 4,2.47 standalone color colortext 10
set output "n_k_a1.tex"
#set key right Right
set key at -.8,.78

set multiplot

set xlabel '$k$'
set ylabel '$n_{\mathrm{imp}}$'

set xrange [-4:4]

#set title '$\alpha = 0$'

plot '../data/n_k_exact_a1_b1.dat' using 1:2 w l lw 5 lc 4 title '$n_{\mathrm{ex}}, \, \beta = 1$',\
'../data/n_k_effective_a1_b1.dat' using 1:2 w l lw 5 dt 2 lc 7  title '$n_{\mathrm{eff}}, \, \beta = 1$', \
'../data/n_k_exact_a1_b2.dat' using 1:2 w l lw 5 title '$n_{\mathrm{ex}}, \, \beta = 2$',\
'../data/n_k_effective_a1_b2.dat' using 1:2 w l lw 5 dt 2 lc 9 title '$n_{\mathrm{eff}}, \, \beta = 2$',

#'../data/n_k_exact_a1_b5.dat' using 1:2 w l lw 5 lc 10 title '$n_{\mathrm{ex}}, \, \beta = 5$',\
#'../data/n_k_effective_a1_b5.dat' using 1:2 w l lw 4 dt 2 lc -1  title '$n_{\mathrm{eff}}, \, \beta = 5$'


set nokey
unset xlabel
unset ylabel
unset title

set format y "10^{%L}"

set logscale y
set origin 0.56, 0.6
set size 0.4, 0.309
set xrange [2:4]
#set yrange [-1:1]
set xtics 1
set ytics 100
#replot


plot '../data/n_k_exact_a0_b1.dat' using 1:2 w l lw 2 lc 4 title '$n_{\mathrm{ex}}, \, \beta = 1$',\
'../data/n_k_effective_a0_b1.dat' using 1:2 w l lw 2 dt 2 lc 7  title '$n_{\mathrm{eff}}, \, \beta = 1$', \
'../data/n_k_exact_a0_b2.dat' using 1:2 w l lw 2 title '$n_{\mathrm{ex}}, \, \beta = 2$',\
'../data/n_k_effective_a0_b2.dat' using 1:2 w l lw 2 dt 2 lc 9 title '$n_{\mathrm{eff}}, \, \beta = 2$',

unset multiplot


# instruction to generate .pdf file from epslatex terminal: 
# gnuplot file.gnu
# pdflatex file.tex
