# Changes

* Formatted code using JuliaFormatter.jl
* Generate a package using PkgTemplates; along with Revise, this makes code development easy
* diff!: Remove since diff does the same thing
* coeff2: Condense the loop (the special case for N[1] is unnecessary)
* coeff2: Removed @vars v
* coeff2: Use eachindex(N) rather than 1:length(N)
* chop_expr: Use zero(Basic) to get the zero of that type rather than Basic(0)
* real!: Changed name to recursive_real (no ! since it is not mutating)
* real!: Changed to use real directly on the Basic type
* chop_expr, real!, imag!: Change to use recursive_apply to reduce duplicated code
* Use #== and ==# for multi-line comments rather than a string
* Put all the global scoped code into a function
* Remove all global variables
* Put input.jl into a module (Systems)
* Fix the mistake with `n = j - l - m` overwriting the variable n.
* Make sure that there aren't any unqualified `[]`, e.g., `ind_list = []` to `ind_list = Vector{Vector{Int}}()`
* Replace the ind_l code with a generator expression
* @code_warntype indicates type instabilities for xs and δ; replace with the generic expressions



# To do

* Verify the examples in a test
