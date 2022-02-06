# Deterministic Optimization Assignment
 Minimization of Rosenbrock’s function using different deterministic optimization methods.

To compile this single script it is required to execute the Makefile from the ./src folder and then
the resulting program called rosenbrock_opt can be executed. It accepts two argument to change
the optimization starting point at run time. En example execution starting at (𝑥1,𝑥2) =(−5,20)
would look like this:

make
./rosenbrock_opt -5 20