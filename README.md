# space-extension-of-explicit-ODEs

This program solves the initial value problem of explicit ordinary differential equations with multinomial right hand side functions. The ordinary differential equation set is transformed into an ordinary differential equation set by with purely second degree right hand side functions by introducing new unknown functions. Beam search is utilized for this transformation. The new ordinary differential equation set is solved by probabilistic evolution theory. Probabilistic evolution theory produces a finite series expansion of the solution. By plugging in different time values and using different truncation levels, the numerical results are shown.

The arithmetic is performed using exact arithmetic (rational numbers), therefore avoiding error buildup because of multiplication of very small numbers with very large numbers. The sparsity of the coefficient matrix is also taken into account, improving the performance.

The systems that are used in the implementation are van der Pol, classical quartic anharmonic oscillator, Henon-Heiles and Rabinovich-Fabrikant. 

The program is compiled with g++ (Ubuntu) 11.3.0 using libraries libgmp 6.2.1 and C++ bindings. 

The command to compile is 
>  g++ main.cpp -lgmpxx -lgmp

which would create the executable a.out.
