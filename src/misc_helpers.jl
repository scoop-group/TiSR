
""" return one weighted sample. Benefit wrt StatsBase, weights can be a tuple.
"""
function weighted_sample(items, weights)
    @assert length(items) == length(weights)
    rand_value = rand() * sum(weights)
    cumulative_weights = cumsum(weights)
    selected_index = findfirst(cumulative_weights .>= rand_value)
    return items[selected_index]
end

""" Recursive functions using multiple dispatch to extract a value from arbitrarily
    nested dual numbers.
"""
extract_from_dual(v::AbstractFloat) = v
extract_from_dual(v::Number) = extract_from_dual(v.value)

""" return a positive or negative float with an absolute within the specified range.
"""
rand_mult(;minn=0.5, maxx=2.0) = (rand() * (maxx - minn) + minn) * rand((1.0, -1.0))
