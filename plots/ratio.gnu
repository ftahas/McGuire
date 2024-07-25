set terminal epslatex size 3.5,2.62 standalone color colortext 10
set output "ratio.tex"
set key right Right

#set multiplot

set nokey

set xlabel '$\beta$'
set ylabel '$\rho_{\mathrm{eff}}(0)/\rho_{\mathrm{ex}}(0)$'

#set xrange [0.2:10]

plot '../data/ratio.dat' w p pt 3 dt 3 lt -1   

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
