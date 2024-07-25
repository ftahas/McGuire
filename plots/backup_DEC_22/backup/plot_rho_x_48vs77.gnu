set terminal epslatex size 3.5,2.62 standalone color colortext 10
set output "plot_rho_x_48vs77.tex"
set key right Right
set xlabel '$x$'
set ylabel '$\rho(x)$'



plot '../data/rho_x_48-b01g10.dat' w l lw 3 title 'Eq. (48), $\beta = 0.1$', \
'../data/rho_x_48-b1g10.dat' w l lw 3 title 'Eq. (48), $\beta=1.0$', \
'../data/rho_x_77-b01g10.dat' w l lw 3 dt 2 title 'Eq. (77), $\beta=0.1$', \
'../data/rho_x_77-b1g10.dat' w l lw 3 dt 4  title 'Eq. (77), $\beta=1.0$', \




# instruction to generate .pdf file from epslatex terminal: 
# gnuplot file.gnu
# pdflatex file.tex
