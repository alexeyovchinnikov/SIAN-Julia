# SIAN Examples

This folder contains a collection of examples of how SIAN can be applied to systems of ordinary differential equations.

These examples can be run in the terminal mode by calling 
```
julia <filename>.jl
```
where `<filename>` is the name of the example file.

Each example can be run in the interactive mode as well by copying and pasting the content of the file in question. 

__Note:__ to supress the unnecessary outputs you can use the semicolon `;` after each code line.

Overall usage recipe is very simple:

 1. Define an ODE system via `@ODEmodel` macro:
    ```julia
        ode = @ODEmodel(x1'(t) = a*x1(t), y(t) = x1(t));
    ```
2. If identifiability of all parameters is desired, use a simple call to `identifiability_ode` with `get_parameters` function like so
    ```julia
        result = identifiability_ode(ode,  get_parameters(ode))
    ```
    More information about input arguments of `identifiability_ode` is available by calling 
    ```julia
        ? identifiability_ode
    ```
3. If a specific parameter identifiability property is of interest, say `a` or `x1`, the call can be modified very simply
    ```julia
       result = identifiability_ode(ode, [a])
       result = identifiability_ode(ode, [x1])  
    ```
