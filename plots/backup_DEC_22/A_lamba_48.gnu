set terminal epslatex size 3.5,2.62 standalone color colortext 10
set output "A_lambda_48_g1_b1.tex"
set key right Right

set xlabel '$\Lambda$'
set ylabel '$A$'

#set title '$g=1, \; \beta = 1$'

#set logscale x

set xrange [-3:3]

plot '../data/xi_lambda_48_g1_b1.dat' using 1:3 w l lw 3 title '$g=1, \beta = 1$'
#'../data/xi_lambda_48_g1_b01.dat' using 1:2 w l lw 3 dt 2 title '$g = 10 , \beta = 1$',\
#'../data/xi_lambda_48_g1_b01.dat' using 1:2 w l lw 3 dt 2 title '$g = 1 , \beta = 0.1$'
#'../data/xi_lambda_48_g10_b01.dat' using 1:2 w l lw 3 dt 2 title '$g = 10 , \beta = 0.1$'


# instruction to generate .pdf file from epslatex terminal: 
# gnuplot file.gnu
# pdflatex file.tex
