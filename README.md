# SIAN-Julia
Implementation of SIAN in Julia
The main function "identifiability_ode" has two required arguments:

1) an ODE model (created by the @ODEmodel macros)
2) array of parameters for which the identifiability analysis is requested

and three optional keys:

1) "p", which is the probabilily of correctness, with the default value p = 0.99, 
2) "p_mod", which is a prime number and is the characteristic of the field over which the computation of Groebner basis will occur. If p = 0 (the default value), the computation will be over the rational numbers. If p > 0, then the computation will be over Z/pZ. The current limit for Groebner basis modular arithmetic implementation suggests that prime numbers bigger than 2^29 - 3 are not to be used.
3) "nthrds" is the maximal number of threads for the Groebner basis computation. The default value is 64.

The function "get_parameters" has one required argument, an ODE model (created by the @ODEmodel macros), and one optional key, "initial_conditions". If the key is set to "true" (the default value), then the function will return the set of all parameters and state variables. If the key is not set to "true", then the function will return the set of all parameters.

The function "generate_replica" has two required arguments:

1) an ODE model (created by the @ODEmodel macros)
2) an integer "r"

and returns the r-fold replica of the ODE model (the state, output, and input variables are replicated and the parameters are not replicated). This function can be used to check the r-experiment identifiability of the parameters.

The folder "examples" contains examples of using this.

The folder "without-macros" contains an earlier version on this implementation that did not use the macros.

If an ODE model has been entered with parameters for some of which it is desirable to further specify their values, this can be done using the "SetParameterValues" function, which inputs:

1) an ODE model (created by the @ODEmodel macros)
2) a dictionary (or odered dictionary) of values such as (taken from the NFkB.jl example)

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
