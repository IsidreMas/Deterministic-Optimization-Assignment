#!/bin/usr/gnuplot  

# The Rosenbrock function [HH Rosenbrock (1960) The Computer Journal 3:175–184]
#    (1-x)**2 + 100 * (y - x**2)**2
# is a well-known difficult test case for general purpose minimizers. 
#
reset
set hidden3d
set isosamples 400
set ticslevel 0.
set xlabel '$x_1$' font ",14" 
set ylabel '$x_2$' font ",14"
Rosenbrock(x,y) = (1-x)**2 + 100*(y - x**2)**2


set key right box opaque width 0.1 height 1 bottom
set pm3d map interpolate 5,5 
set logscale cb
set format cb  "10^{%T}"
set multiplot
splot [-2:2] [-2:2] Rosenbrock(x,y)
#splot [-2:2] [-2:2] "../results/ConjugateGradientMethod.dat" u 1:2:3 w linespoints lc "black" lw 2 pt 2 title "Congugate Gradient Method"
splot [-2:2] [-2:2] "../results/SteepestGradientDescent.dat" u 1:2:3 w linespoints lc "black" lw 2 pt 2 title "Gradient Method"
#splot [-2:2] [-2:2] "../results/LevenbergMarquardt.dat" u 1:2:3 w linespoints lc "black" lw 2 pt 2 title "Levenberg-Marquardt "

unset multiplot