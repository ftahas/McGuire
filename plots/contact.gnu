set terminal epslatex size 4,2.47 standalone color colortext 10
set output "contact.tex"
set key right Right

set multiplot

set xlabel '$\beta$'
set ylabel '$C$'

set logscale y

#set xrange [0:1]

plot '../data/contact_a01.dat' using 1:2 w l lw 5 dt 1 title '$\alpha = 0.1$',\
 '../data/contact_a1.dat' using 1:2 w l lw 5 dt 2 title '$\alpha = 1$', \
 '../data/contact_a10.dat' using 1:2 w l lw 5 dt 3 lt -1 title '$\alpha = 10$'

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
