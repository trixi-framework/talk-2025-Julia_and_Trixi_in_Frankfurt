# An Overview of Atmospheric Applications in Trixi.jl

The Pluto notebook and the codes have been used on Julia version 1.10.2. 

Please see [introduction_to_julia/README.md](https://github.com/trixi-framework/talk-2025-Julia_and_Trixi_in_Frankfurt/blob/main/introduction_to_julia/README.md) 
for instructions on running Pluto notebooks. However is it strongly recommended to run the Pluto notebook with multiple threads modifying the Pluto call as follows
```julia
Pluto.run(threads = 8)
```
In order to run the examples, navigate to the `codes` directory, activate the environment within that folder, and install the required dependencies:

```julia
'import Pkg; Pkg.activate(pwd()); Pkg.instantiate();
```


