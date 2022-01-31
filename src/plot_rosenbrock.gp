#
# The Rosenbrock function [HH Rosenbrock (1960) The Computer Journal 3:175–184]
#    (1-x)**2 + 100 * (y - x**2)**2
# is a well-known difficult test case for general purpose minimizers. 
#
reset
set hidden3d
set logscale z
set isosamples 400
set ticslevel 0.
set view 56,13   #HBB: ,1,2
set xlabel '$x_1$' font ",14" 
set ylabel '$x_2$' font ",14"
set zlabel '$f(x_1,x_2)$' font ",14" offset -5,0,0
set format z  "10^{%T}"
Rosenbrock(x,y) = (1-x)**2 + 100*(y - x**2)**2

unset key
set pm3d
set logscale cb
set format cb  "10^{%T}"
splot [-2:2] [-2:2] Rosenbrock(x,y) title "Rosenbrock's function"