set terminal epslatex size 3.5,2.62 standalone color colortext 10
set output "plot_rho_x_lambda10_48vs77.tex"
set key right Right
set xlabel '$x$'
set ylabel '$\rho(x,\lambda)$'



plot '../data/rho_x_lambda_48.dat' w l lw 3 title '$\lambda = 10$, Eq. (48)', \
'../data/rho_x_lambda_77.dat' w l lw 3 title '$\lambda=10$, Eq. (77)'

# instruction to generate .pdf file from epslatex terminal: 
# gnuplot file.gnu
# pdflatex file.tex
