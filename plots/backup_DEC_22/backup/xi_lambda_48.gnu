set terminal epslatex size 3.5,2.62 standalone color colortext 10
set output "xi_lambda_48.tex"
set key right Right
set xlabel '$\Lambda$'
set ylabel '$\xi$'


plot '../data/xi_lambda_48_x1=5_x2=6.dat' using 1:3 w l lw 3 title '', \
'../data/xi_lambda_48_x1=8_x2=10.dat' using 1:3 w l lw 3 dt 2 title '', \
'../data/xi_lambda_48_x1=10_x2=20.dat' using 1:3 w l lw 3 dt 4 title '', \
'../data/xi_lambda_48_x1=15_x2=20.dat' using 1:3 w l lw 3 dt 5 title ''

# instruction to generate .pdf file from epslatex terminal: 
# gnuplot file.gnu
# pdflatex file.tex
