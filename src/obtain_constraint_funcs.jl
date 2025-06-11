
function obtain_value_constr_func(constr_data, smaller_larger::String, target::Float64; replace_inf=1e100)
    @assert smaller_larger in ("<", ">")
    return let constr_data = constr_data, smaller_larger = smaller_larger, target = target
        (node_func, params) -> begin
            constr_data_ = convert.(typeof(params), constr_data)
            pred = node_func(params, constr_data_)
            clamp!(pred, -replace_inf, replace_inf)
            if smaller_larger == "<"
                return max.(0.0, pred .- target)
            else
                return max.(0.0, target .- pred)
            end
        end
    end
end

function obtain_monotonicity_constr_func(constr_data, var_ind::Int64, smaller_larger::String, target::Float64; replace_inf=1e100)
    @assert smaller_larger in ("<", ">")
    slope = similar(constr_data[1])
    return let constr_data = constr_data, smaller_larger = smaller_larger, target = target, var_ind = var_ind, slope=slope
        (node_func, params) -> begin
            constr_data_ = convert.(typeof(params), constr_data)
            slope_ = convert(typeof(params), slope)
            ForwardDiff.derivative!(
                slope_,
                x -> node_func(params,
                    promote_and_replace_element(constr_data_, var_ind, x)
                ),
                constr_data_[var_ind]
            )
            clamp!(slope_, -replace_inf, replace_inf)
            if smaller_larger == "<"
                return max.(0.0, slope_ .- target)
            else
                return max.(0.0, target .- slope_)
            end
        end
    end
end

function promote_and_replace_element(arr::Vector{T1}, ind, elem::T2) where {T1, T2}
    T = promote_type(T1, T2)
    narr = T.(arr)
    narr[ind] = elem
    return narr
end
