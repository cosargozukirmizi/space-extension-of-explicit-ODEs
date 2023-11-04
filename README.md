# space-extension-of-explicit-ODEs

This program solves the initial value problem of explicit ordinary differential equations with multinomial right hand side functions. The ordinary differential equation set is transformed into an ordinary differential equation set by with purely second degree right hand side functions by introducing new unknown functions. Beam search is utilized for this transformation. The new ordinary differential equation set is solved by probabilistic evolution theory. Probabilistic evolution theory produces a finite series expansion of the solution. By plugging in different time values and using different truncation levels, the numerical results are shown.

The arithmetic is performed using exact arithmetic (rational numbers), therefore avoiding error buildup because of multiplication of very small numbers with very large numbers. The sparsity of the coefficient matrix is also taken into account, improving the performance.

The systems that are used in the implementation are van der Pol, classical quartic anharmonic oscillator, Henon-Heiles and Rabinovich-Fabrikant. 

The program is compiled with g++ (Ubuntu) 11.3.0 using libraries libgmp 6.2.1 and libgmp C++ bindings. The program has been tested with GNU compiler g++-10 and above. 

The command to compile is 
>  g++ -std=c++20 main.cpp -lgmp -lgmpxx -o main.x

which would create the executable main.x.

# Ordinary differential equations under consideration

van der Pol
-----------
The van der Pol ODE is
```math
\begin{eqnarray}
 \dot{x} &=& \mu x - \frac{1}{3}\mu x^{3} - \mu y \\
 \dot{y} &=& \frac{1}{\mu} x 
\end{eqnarray}
```
where there are two sought functions. 
The definition
```math
\begin{equation}
 u^{(\ell_{1},\ell_{2})}(t) \equiv x(t)^{\ell_{1}}y(t)^{\ell_{2}} 
\end{equation}
```
can be used to rewrite in the form
```math
\begin{eqnarray}
  \dot{u}^{(1,0)} &=& \mu \, u^{(0,0)} u^{(1,0)} 
  - \frac{1}{3}\mu \, u^{(1,0)} u^{(2,0)} 
  - \mu \, u^{(0,0)} u^{(0,1)}  \\
  \dot{u}^{(0,1)} &=& \frac{1}{\mu}\, u^{(0,0)} u^{(1,0)} 
\end{eqnarray}
```
using the heuristic H4 used in beam search for pure quadratization. We will take $\mu$ as 1. After pure quadratization, the ODE set is 
```math
\begin{eqnarray}
  \dot{u}^{(1,0)} &=& 1\, u^{(0,0)} u^{(1,0)} 
  - 1/3\, u^{(1,0)} u^{(2,0)} 
  - 1\, u^{(0,0)} u^{(0,1)}
  \\
  \dot{u}^{(0,1)} &=& 1\, u^{(0,0)} u^{(1,0)}\\
  \dot{u}^{(0,0)} &=& 0\\
  \dot{u}^{(2,0)} &=& 2\, u^{(1,0)} u^{(1,0)} 
  - 2/3\, u^{(2,0)} u^{(2,0)} 
  - 2\, u^{(0,1)} u^{(1,0)}
\end{eqnarray}
```

Quartic anharmonic oscillator
-----------------------------
The Quartic anharmonic oscillator ODE is
```math
\begin{eqnarray}
 \dot{q} &=& \frac{1}{\mu} p \\
 \dot{p} &=& - k_{1} q - k_{2} q^{3}\\
\end{eqnarray}
``` 
where there are two sought functions. 
The definition
```math
\begin{equation}
 u^{(\ell_{1},\ell_{2})}(t) \equiv q(t)^{\ell_{1}}p(t)^{\ell_{2}} 
\end{equation}
```
can be used to rewrite in the form
```math
\begin{eqnarray}
  \dot{u}^{(1,0)} &=& \frac{1}{\mu} \, u^{(0,0)} u^{(0,1)}\\
  \dot{u}^{(0,1)} &=& -k_{1}\, u^{(0,0)} u^{(1,0)} \nonumber \\
 &-& k_{2}\, u^{(1,0)} u^{(2,0)}\\
\end{eqnarray}
```
using the heuristic H4 used in beam search for pure quadratization. We will take $k_{1}$, $k_{2}$ and $\mu$ as 1. After pure quadratization, the ODE set is 
```math
\begin{eqnarray}
  \dot{u}^{(1,0)} &=& 1\, u^{(0,0)} u^{(0,1)}\\
  \dot{u}^{(0,1)} &=& -1\, u^{(0,0)} u^{(1,0)} \nonumber \\
 &-& 1\, u^{(1,0)} u^{(2,0)}\\
  \dot{u}^{(0,0)} &=& 0\\
  \dot{u}^{(2,0)} &=& 2\, u^{(0,1)} u^{(1,0)}
\end{eqnarray}
```

Henon-Heiles
------------
The Henon-Heiles ODE is
```math
\begin{eqnarray}
  \dot{x} &=& p_{x}\\
  \dot{p}_{x} &=& -x \nonumber \\
 &-& 2\, \lambda xy\\
  \dot{y} &=& p_{y}\\
  \dot{p}_{y} &=& -y \nonumber \\
 &-& \lambda x^{2} \nonumber \\
 &+& \lambda y^{2} 
\end{eqnarray}
```
where there are four sought functions. 
The definition
```math
\begin{equation}
 u^{(\ell_{1},\ell_{2},\ell_{3},\ell_{4})}(t) 
 \equiv x(t)^{\ell_{1}}p_{x}(t)^{\ell_{2}}y(t)^{\ell_{3}}
p_{y}(t)^{\ell_{4}} 
\end{equation}
```
can be used to rewrite in the form
```math
\begin{eqnarray}
  \dot{u}^{(1,0,0,0)} &=& 1\, u^{(0,0,0,0)} u^{(0,1,0,0)}\\
  \dot{u}^{(0,1,0,0)} &=& -1\, u^{(0,0,0,0)} u^{(1,0,0,0)} \nonumber \\
 &-& 2\lambda \, u^{(0,0,1,0)} u^{(1,0,0,0)}\\
  \dot{u}^{(0,0,1,0)} &=& 1\, u^{(0,0,0,0)} u^{(0,0,0,1)}\\
  \dot{u}^{(0,0,0,1)} &=& -1\, u^{(0,0,0,0)} u^{(0,0,1,0)} \nonumber \\
 &-& \lambda \, u^{(1,0,0,0)} u^{(1,0,0,0)} \nonumber \\
 &+& \lambda \, u^{(0,0,1,0)} u^{(0,0,1,0)}
\end{eqnarray}
``` 
using the heuristic H4 used in beam search for pure quadratization. We will take $\lambda$ as 1. After pure quadratization, the ODE set is 
```math
\begin{eqnarray}
  \dot{u}^{(1,0,0,0)} &=& 1\, u^{(0,0,0,0)} u^{(0,1,0,0)}\\
  \dot{u}^{(0,1,0,0)} &=& -1\, u^{(0,0,0,0)} u^{(1,0,0,0)} \nonumber \\
 &-& 2\, u^{(0,0,1,0)} u^{(1,0,0,0)}\\
  \dot{u}^{(0,0,1,0)} &=& 1\, u^{(0,0,0,0)} u^{(0,0,0,1)}\\
  \dot{u}^{(0,0,0,1)} &=& -1\, u^{(0,0,0,0)} u^{(0,0,1,0)} \nonumber \\
 &-& 1\, u^{(1,0,0,0)} u^{(1,0,0,0)} \nonumber \\
 &+& 1\, u^{(0,0,1,0)} u^{(0,0,1,0)}\\
  \dot{u}^{(0,0,0,0)} &=& 0
\end{eqnarray}
```

Rabinovich-Fabrikant
--------------------
The Rabinovich-Fabrikant ODE is 
```math
\begin{eqnarray}
 \dot{x} &=& yz - y + x^{2}y + \gamma x \\
 \dot{y} &=& 3xz + x - x^{3} + \gamma y \\
 \dot{z} &=& -2\alpha z - 2 xyz 
\end{eqnarray}
```
where there are two sought functions. 
The definition
```math
\begin{equation}
 u^{(\ell_{1},\ell_{2},\ell_{3})}(t) 
 \equiv x(t)^{\ell_{1}}y(t)^{\ell_{2}}z(t)^{\ell_{3}} 
\end{equation}
```
can be used to rewrite in the form
```math
\begin{eqnarray}
\dot{u}^{(1,0,0)} &=& 1\, u^{(0,0,1)} u^{(0,1,0)} 
- 1\, u^{(0,0,0)} u^{(0,1,0)} \nonumber\\
&+& 1\, u^{(1,0,0)} u^{(1,1,0)} 
+ \gamma \, u^{(0,0,0)} u^{(1,0,0)}\\
\dot{u}^{(0,1,0)} &=& 3\, u^{(0,0,1)} u^{(1,0,0)} 
+ 1\, u^{(0,0,0)} u^{(1,0,0)} \nonumber\\
&-& 1\, u^{(1,0,0)} u^{(2,0,0)} 
+ \gamma \, u^{(0,0,0)} u^{(0,1,0)} \\
\dot{u}^{(0,0,1)} &=& -2\alpha \, u^{(0,0,0)} u^{(0,0,1)} 
- 2\, u^{(0,1,0)} u^{(1,0,1)} 
\end{eqnarray}
```
using the heuristic H4 used in beam search for pure quadratization. We will take $\alpha$ and $\gamma$ as 1. After pure quadratization, the ODE set is 
```math
\begin{align}
\dot{u}^{(1,0,0)} &= 1\, u^{(0,0,1)} u^{(0,1,0)}  \\
 &- 1\, u^{(0,0,0)} u^{(0,1,0)}  \\
 &+ 1\, u^{(1,0,0)} u^{(1,1,0)} \\
 &+ 1\, u^{(0,0,0)} u^{(1,0,0)}\\
  \dot{u}^{(0,1,0)} &= 3\, u^{(0,0,1)} u^{(1,0,0)} \\
 &+ 1\, u^{(0,0,0)} u^{(1,0,0)}  \\
 &- 1\, u^{(1,0,0)} u^{(2,0,0)}  \\
 &+ 1\, u^{(0,0,0)} u^{(0,1,0)}\\
  \dot{u}^{(0,0,1)} &= -2\, u^{(0,0,0)} u^{(0,0,1)}  \\
 &- 2\, u^{(0,1,0)} u^{(1,0,1)}\\
  \dot{u}^{(0,0,0)} &= 0
\end{align}
```
```math
\begin{align}
  \dot{u}^{(1,1,0)} &= 1\, u^{(0,1,0)} u^{(0,1,1)}  \\
 &- 1\, u^{(0,1,0)} u^{(0,1,0)}  \\
 &+ 1\, u^{(1,1,0)} u^{(1,1,0)}  \\
 &+ 1\, u^{(0,1,0)} u^{(1,0,0)}  \\
 &+ 3\, u^{(1,0,0)} u^{(1,0,1)}  \\
 &+ 1\, u^{(1,0,0)} u^{(1,0,0)}  \\
 &- 1\, u^{(2,0,0)} u^{(2,0,0)}  \\
 &+ 1\, u^{(0,1,0)} u^{(1,0,0)}\\
  \dot{u}^{(2,0,0)} &= 2\, u^{(0,1,0)} u^{(1,0,1)}  \\
 &- 2\, u^{(0,1,0)} u^{(1,0,0)}  \\
 &+ 2\, u^{(1,1,0)} u^{(2,0,0)}  \\
 &+ 2\, u^{(1,0,0)} u^{(1,0,0)}\\
  \dot{u}^{(1,0,1)} &= 1\, u^{(0,0,1)} u^{(0,1,1)}  \\
 &- 1\, u^{(0,0,1)} u^{(0,1,0)}  \\
 &+ 1\, u^{(1,0,1)} u^{(1,1,0)}  \\
 &+ 1\, u^{(0,0,1)} u^{(1,0,0)}  \\
 &- 2\, u^{(0,0,1)} u^{(1,0,0)}  \\
 &- 2\, u^{(1,0,1)} u^{(1,1,0)}\\
  \dot{u}^{(0,1,1)} &= 3\, u^{(0,0,1)} u^{(1,0,1)}  \\
 &+ 1\, u^{(0,0,1)} u^{(1,0,0)}  \\
 &- 1\, u^{(1,0,1)} u^{(2,0,0)}  \\
 &+ 1\, u^{(0,0,1)} u^{(0,1,0)}  \\
 &- 2\, u^{(0,0,1)} u^{(0,1,0)}  \\
 &- 2\, u^{(0,1,1)} u^{(1,1,0)}
\end{align}
```
