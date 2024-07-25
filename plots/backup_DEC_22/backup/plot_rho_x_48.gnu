set terminal epslatex size 3.5,2.62 standalone color colortext 10
set output "plot_rho_x_48_g10_b1.tex"
set key right Right
set xlabel '$x$'
set ylabel '$\rho(x)$'



plot '../data/rho_x_48_g10_b1.dat' using 1:2 w l lw 3 title '$\mathrm{Re}(\rho)$', \
'../data/rho_x_48_g10_b1.dat' using 1:3 w l lw 3 dt 2 title '$\mathrm{Im}(\rho)$', \
'../data/rho_x_48_g10_b1.dat' using 1:4 w l lw 3 dt 4 title '$\vert \rho \vert$'




# instruction to generate .pdf file from epslatex terminal: 
# gnuplot file.gnu
# pdflatex file.tex
