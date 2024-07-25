set terminal epslatex size 4,2.47 standalone color colortext 12
set output "fig1_nu.tex"
#set key left Left
set key at 0.2,1
#set multiplot

set xlabel '$q$'
set ylabel '$\nu(q)$'

set xrange [-4:4]
set yrange [-0.2:1.1]

set xzeroaxis # linetype 3 linewidth 2.5
set yzeroaxis

plot '../data/fig1_winding.dat' using 3:4 w l lw 5 dt 1 title '$w = 1, \, \mathrm{Re}\nu$',\
 '../data/fig1_winding.dat' using 3:5 w l lw 5 dt 3 lt -1 title '$\mathrm{Im} \nu $', \
'../data/fig1_nonwinding.dat' using 3:4 w l lw 5 dt 2 lc 2 title '$w=0, \, \mathrm{Re} \nu $',\
 '../data/fig1_nonwinding.dat' using 3:5 w l lw 5 dt 4 title '$\mathrm{Im} \nu $'

#'../data/origin.dat' w p pt 3 dt 3 lt -1 notitle,  

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
