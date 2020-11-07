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
