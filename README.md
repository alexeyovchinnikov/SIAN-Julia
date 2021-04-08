# SIAN-Julia
![SIAN Tests for Julia](https://github.com/iliailmer/SIAN-Julia/workflows/SIAN%20Tests%20for%20Julia/badge.svg)

## Installing

The current stable version of SIAN.jl can be installed via the following command:

```zsh
>] add "SIAN"
> using SIAN
```

The installation from source is possible via GitHub and SSH or HTTPS:
```zsh
> ]add https://github.com/alexeyovchinnikov/SIAN-Julia.git
```
or
```zsh
> ]add git@github.com:alexeyovchinnikov/SIAN-Julia.git
```

## Example: getting started

In this example we would like to consider the following simple non-linear ODE system:
```
x1'(t) = r1*x1(t)*(1 - x1(t)/k1 + x2(t)/k2),
x2'(t) = r2*x2(t)*(1 - x1(t)/k1 + x2(t)/k2)
y(t) = x1(t)
```

To this end, we can run:
```julia
using SIAN

ode = @ODEmodel(
    x1'(t) = r1*x1(t)*(1 - x1(t)/k1 + x2(t)/k2),
    x2'(t) = r2*x2(t)*(1 - x1(t)/k1 + x2(t)/k2),
    y(t) = x1(t)
);

output = identifiability_ode(ode, get_parameters(ode));
```

The last command prints the following:
```
Solving the problem
Constructing the maximal system
Truncating
Assessing local identifiability
Locally identifiable parameters: [r1, k1, r2, x1]
Not identifiable parameters:     [k2, x2]
Randomizing
GB computation
Remainder computation

=== Summary ===
Globally identifiable parameters:                 [x1, k1, r1, r2]
Locally but not globally identifiable parameters: []
Not identifiable parameters:                      [k2, x2]
===============
```

## Scoping

Once the module is imported via `using SIAN`, the following functions are available immediately via the `export` of the module: `@ODEmodel, identifiability_ode, get_parameters`.

Other SIAN functions are available via prefix call, such as `SIAN.<any_function_name>`.

### Getting started

We recommend checking the `examples` folder to get started with using SIAN, see [this readme file](./examples/README.md).
## Implementation of SIAN in Julia

The algorithm is based on the following papers:
* [Global Identifiability of Differential Models](https://onlinelibrary.wiley.com/doi/abs/10.1002/cpa.21921) (Communications on Pure and Applied Mathematics, Volume 73, Issue 9, Pages 1831-1879, 2020.) 
* [SIAN: software for structural identifiability analysis of ODE models](https://doi.org/10.1093/bioinformatics/bty1069) (Bioinformatics, Volume 35, Issue 16, Pages 2873–2874, 2019)

The original Maple implementation is located [here](https://github.com/pogudingleb/SIAN).

The main function "identifiability_ode" has two required arguments:

1) an ODE model (created by the `@ODEmodel` macros)
2) array of parameters for which the identifiability analysis is requested

and three optional keys:

1) `p`, which is the probability of correctness, with the default value `p = 0.99`, 
2) `p_mod`, which is a prime number and is the characteristic of the field over which the computation of Groebner basis will occur. If `p = 0` (the default value), the computation will be over the rational numbers. If `p > 0`, then the computation will be over `Z/pZ`. The current limit for Groebner basis modular arithmetic implementation suggests that prime numbers bigger than `2^29 - 3` are not to be used. When using `p_mod>0`, the same probability of correctness is _no longer guaranteed_ and the program will raise a warning message in that case.
3) `nthrds` is the maximal number of threads for the Groebner basis computation. The default value is 1.

The function `get_parameters` has one required argument, an ODE model (created by the `@ODEmodel` macros), and one optional key, `initial_conditions`. If the key is set to `true` (the default value), then the function will return the set of all parameters and state variables. If the key is not set to `true`, then the function will return the set of all parameters.

The function "generate_replica" has two required arguments:

1) an ODE model (created by the `@ODEmodel` macros)
2) an integer `r`

and returns the `r`-fold replica of the ODE model (the state, output, and input variables are replicated and the parameters are not replicated). This function can be used to check the `r`-experiment identifiability of the parameters.

The folder `examples/` contains examples of using this.

The folder `without-macros/` contains an earlier version on this implementation that did not use the macros.

If an ODE model has been entered with parameters for some of which it is desirable to further specify their values, this can be done using the `SetParameterValues` function, which accepts:

1) an ODE model (created by the `@ODEmodel` macros)
2) a dictionary (or ordered dictionary) of values such as (taken from the `NFkB.jl` example)

```julia
OrderedDict(

    a1 => Nemo.QQ(1, 2),   
    
    a2 => Nemo.QQ(1, 5),
    
    a3 => Nemo.QQ(1),
    
    c_1a => Nemo.QQ(5, 10^(7)),
    
    c_2a => Nemo.QQ(0),
    
    c_5a => Nemo.QQ(1, 10^(4)),
    
    c_6a => Nemo.QQ(2, 10^(5)),
    
    c1 => Nemo.QQ(5, 10^(7)),
    
    c2 => Nemo.QQ(0),
    
    c3 => Nemo.QQ(4, 10^(4)),
    
    c4 => Nemo.QQ(1, 2),
    
    kv => Nemo.QQ(5),
    
    e_1a => Nemo.QQ(5, 10^(4)),
    
    c_1c => Nemo.QQ(5, 10^(7)),
    
    c_2c => Nemo.QQ(0),
    
    c_3c => Nemo.QQ(4, 10^(4))

)
```

for instance, to specify that a1 is 1/2, etc.