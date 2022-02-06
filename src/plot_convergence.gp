#
# The Rosenbrock function [HH Rosenbrock (1960) The Computer Journal 3:175–184]
#    (1-x)**2 + 100 * (y - x**2)**2
# is a well-known difficult test case for general purpose minimizers. 
#
reset
set hidden3d
set ticslevel 0.
set view 56,13   #HBB: ,1,2
set xlabel '$x_1$' font ",14" 
set ylabel '$x_2$' font ",14"
unset key
set pm3d
set dgrid3d 200,200,50

splot [-450:450] [-450:450] "../results/Convergence1.dat" u 1:2:3 title "Rosenbrock's function"
#splot [-450:450] [-450:450] "../results/Convergence2.dat" u 1:2:3 title "Rosenbrock's function"
#splot [-400:400] [-450:450] "../results/Convergence3.dat" u 1:2:3 title "Rosenbrock's function"