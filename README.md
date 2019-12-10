# HopfNormalForm

HopfNormalForm is a Julia package to compute a center manifold and Hopf normal form of a nonlinear autonomous differential equation. For the current version, only non-semisimple Jacobian at Hopf point is available with the dimension of the center manifold of 2 and 3.

This package has center manifold and Hopf normal form and nonlinear parameter optimization result of the journal "Exploring features of a dynamical system near Hopf bifurcation using control based continuation".

module

## Installation

Within Julia, you can install the HopfNormalForm.jl package with the package manager:
```Julia
using Pkg
Pkg.add("https://github.com/Kyounghyunlee/HopfNormalForm.jl")
```


## Usage
To compute the center manifold and Hopf normal form function use function calculate_normal_form

```Julia
function calculate_normal_form(
    f::Vector{Basic}, # Define equation of motion here
    v::Vector{Basic}, # this is a state-space vector
    norder::Integer,  # Choose the order of power series representation of Centermanifold
    cm::Integer, # Dimension of Centermanifold of the dynamical system
    m1::Integer, # Number of the zero eigenvalue corresponding to the center space of the Jacobian
    m2::Integer, # Number of pair of complex conjugate eigenvalue corresponding to the center space of the Jacobian
    m3::Integer, # Dimension of the stable manifold of the dynmical system
    n::Integer, # Total dimension of the dynamical system
    ω_h::AbstractFloat, # Eigenvalue that crosses unit circle at Hopf bifurcation point ± ω_h
    zero_tol::AbstractFloat # zero tolerance to chop out very small coeffs of symbolic computation
)
```
Input function of
"Exploring features of a dynamical system near Hopf bifurcation using control based continuation"
for system 1 and system 2 is automatically generated from the matlab code (can be found in [System_ID_flutter](https://github.com/Kyounghyunlee/System_ID_flutter))

To see the result of the center manifold and normal form computation of mathematical model provided in "Exploring features of a dynamical system near Hopf bifurcation using control based continuation":

```Julia
using HopfNormalForm
Systems.example1() # returns center manifold and Hopf normal form of system 1
Systems.example2() # returns center manifold and Hopf normal form of system 2
```

To see the result of the nonlinear parameter optimization results:

```Julia
using HopfNormalForm
Flutter_opt.optimize_NLstiff(ka2_0,ka3_0) # returns optimized value of nonlinear stifness of system 1 with initial searching point [ka2_0,ka3_0]
Flutter_opt.optimize_NLstiff2(ka2_0,ka3_0) # returns optimized value of nonlinear stifness of system 2 with initial searching point [ka2_0,ka3_0]
```

## License
[MIT](https://choosealicense.com/licenses/mit/)
